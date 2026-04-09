// sgp4_plugin.cpp - SGP4/SDP4 Propagator Plugin for OrbPro
// =============================================================================
// Uses Vallado's reference SGP4 implementation (SGP4.cpp / SGP4.h) directly.
// No modifications to the core algorithm — line-for-line Vallado 2020-07-13.
//
// Exports:
//   plugin_init_omm(records, count)           → init from binary OMM records
//   plugin_init_omm_flatbuffer_stream(data,len) → init from OMM FlatBuffer bytes
//   plugin_propagate(jd, entity_index, out)   → propagate single entity
//   plugin_propagate_batch(jd, out, count)    → propagate all entities
//   plugin_destroy()                          → cleanup
//   get_entity_info(entity_index, out)        → get NORAD ID
//
//   Multi-OMM per entity (added for ephemeris update blending):
//   plugin_entity_add_omm(entity_index, record)    → add OMM to entity
//   plugin_entity_add_omm_flatbuffer(i, bytes,len) → add `$OMM` payload
//   plugin_query_entity_indices_by_name(q,out,max) → SQL name query (entity indices)
//   plugin_query_entity_rows_by_name(q,out,max)    → SQL name query (rows + names)
//   plugin_query_visibility_mask_by_name(q,mask,n) → SQL name query (0/1 mask)
//   plugin_get_entity_catalog_row(i,out)           → row metadata by entity index
//   plugin_catalog_upsert_cat_record(...)          → CAT ingest/index upsert
//   plugin_catalog_upsert_cat_flatbuffer_stream(...) → CAT/REC FlatBuffer ingest
//   plugin_catalog_upsert_cat_records_packed(...)  → bulk CAT ingest/index upsert
//   plugin_get_cat_record_json_*                   → CAT payload retrieval
//   plugin_entity_omm_count(entity_index)           → count of OMMs
//   plugin_entity_list_omm(entity_index, out, max)  → list epoch JDs
//   plugin_entity_remove_omm(entity_index, jd)      → remove by epoch
//   plugin_entity_remove_all_omm_except(i, jd)      → keep only one
//   plugin_entity_clear_omm(entity_index)            → clear all OMMs
//   plugin_entity_set_mode(entity_index, mode)       → 0=nearest, 1=interp
//   plugin_entity_get_mode(entity_index)             → get mode
//   plugin_get_omm_record_by_pointer(ptr,out)        → OMM payload by sqlite id
//
// All positions are output in ECEF frame, in meters.
// =============================================================================

#include "orbpro_plugin.h"
#include "orbpro_propagator.h"
#include "SGP4.h"
#include "generated/CatalogQueryRequest_generated.h"
#include "generated/CatalogQueryResult_generated.h"
#include "generated/EntityMetadata_generated.h"
#include "generated/PluginMessage_generated.h"
#include "generated/PropagatorTrajectorySegments_generated.h"
#include "generated/PropagatorState_generated.h"
#include "generated/StateVector_generated.h"
#include "generated/TypedArenaBuffer_generated.h"
#include <flatbuffers/flatbuffers.h>
#include <sqlite3.h>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <map>
#include <cstdint>
#include <string>
#include <cstdio>
#include <algorithm>
#include <array>
#include <limits>
#include <map>
#include <unordered_map>

// =============================================================================
// Binary OMM Record — written directly into linear memory by JavaScript.
// Matches SpaceDataStandards OMM schema fields needed for SGP4 initialization.
// All angular values are in DEGREES, mean motion in REV/DAY.
// JS parses OMM JSON and packs these structs into WASM heap.
// =============================================================================
struct OrbProOMMRecord {
    double epoch_jd;          //  0: Epoch as Julian Date (JS converts ISO 8601 → JD)
    double mean_motion;       //  8: Mean motion (rev/day)
    double eccentricity;      // 16: Eccentricity (unitless)
    double inclination;       // 24: Inclination (degrees)
    double ra_of_asc_node;    // 32: RA of ascending node (degrees)
    double arg_of_pericenter; // 40: Argument of pericenter (degrees)
    double mean_anomaly;      // 48: Mean anomaly (degrees)
    double bstar;             // 56: B* drag term (1/earth radii)
    double mean_motion_dot;   // 64: First deriv of mean motion (rev/day²)
    double mean_motion_ddot;  // 72: Second deriv of mean motion (rev/day³)
    uint32_t norad_cat_id;    // 80: NORAD catalog ID
    uint32_t _pad;            // 84: Padding for 8-byte alignment
};                            // Total: 88 bytes

static constexpr size_t ORBPRO_NAME_BUFFER_BYTES = 96;

struct OrbProCatalogRow {
    uint32_t entity_index;
    uint32_t norad_cat_id;
    char object_name[ORBPRO_NAME_BUFFER_BYTES];
    char object_id[ORBPRO_NAME_BUFFER_BYTES];
    char cat_object_name[ORBPRO_NAME_BUFFER_BYTES];
    char cat_object_id[ORBPRO_NAME_BUFFER_BYTES];
};

struct PendingCatalogMetadata {
    std::string objectName;
    std::string objectId;
};

struct PendingCatCatalogRecord {
    uint32_t noradCatId = 0;
    std::string objectName;
    std::string objectId;
};

// =============================================================================
// Constants
// =============================================================================

static constexpr double PI = 3.14159265358979323846;
static constexpr double TWOPI = 2.0 * PI;
static constexpr double DEG_TO_RAD = PI / 180.0;
static constexpr double OMEGA_EARTH = 7.2921151467e-5; // rad/s
static constexpr int CHEBY_N = 12;
static constexpr int CHEBY_NPTS = CHEBY_N + 1;
static constexpr double CHEBY_SEG_SEC = 600.0;
static constexpr double CHEBY_SEG_DAYS = CHEBY_SEG_SEC / 86400.0;
static constexpr double CHEBY_MIN_SEG_SEC = 60.0;
static constexpr double CHEBY_MIN_SEG_DAYS = CHEBY_MIN_SEG_SEC / 86400.0;
static constexpr double CHEBY_CERT_POSITION_KM = 1.0e-3;
static constexpr double CHEBY_CERT_VELOCITY_KMS = 1.0e-6;

// Propagation modes for multi-OMM entities
static constexpr uint8_t MODE_NEAREST_EPOCH = 0;
static constexpr uint8_t MODE_INTERPOLATED_EPOCH = 1;

// =============================================================================
// Plugin Metadata
// =============================================================================

static OrbProPluginInfo PLUGIN_INFO = {
    "com.orbpro.sgp4",
    "SGP4/SDP4 Propagator",
    "2.0.0",
    ORBPRO_ABI_VERSION
};

// =============================================================================
// Satellite Entity Storage
// =============================================================================

struct SatelliteEntity {
    uint32_t noradId;
    elsetrec satrec;       // Current/active satellite record
    double r[3];           // Last position (km, TEME)
    double v[3];           // Last velocity (km/s, TEME)
    double lastEpochJD;
    bool valid;

    // Multi-OMM storage: map of epoch JD → satrec
    // When empty, single-satrec mode (entity.satrec is the sole source).
    // When populated, entity.satrec is a copy of the currently active satrec.
    std::map<double, elsetrec> satrecs;
    std::map<double, int64_t> ommPointers; // epoch JD -> sqlite row id
    int64_t currentOmmPointer;
    double currentSatrecEpoch;  // JD of the currently active satrec
    uint8_t mode;               // MODE_NEAREST_EPOCH or MODE_INTERPOLATED_EPOCH
};

struct CertifiedTrajectorySegment {
    uint32_t sourceHandle;
    double startJd;
    double endJd;
    std::array<std::array<double, CHEBY_NPTS>, 6> coeffs;
    double maxPositionErrorKm;
    double maxVelocityErrorKmS;
};

struct PreparedTrajectorySegmentSet {
    uint32_t handle;
    uint32_t catalogHandle;
    double startJd;
    double durationDays;
    std::string profile;
    std::vector<uint32_t> sourceHandles;
    std::vector<CertifiedTrajectorySegment> segments;
    double maxPositionErrorKm;
    double maxVelocityErrorKmS;
    bool coverageComplete;
};

static bool g_initialized = false;
static std::vector<SatelliteEntity> g_satellites;
static std::unordered_map<uint32_t, PreparedTrajectorySegmentSet> g_segmentSets;
static uint32_t g_nextSegmentSetHandle = 1;
static sqlite3* g_ommDb = nullptr;
static sqlite3_stmt* g_insertOmmStmt = nullptr;
static sqlite3_stmt* g_upsertEntityStmt = nullptr;
static sqlite3_stmt* g_upsertCatStmt = nullptr;
static sqlite3_stmt* g_queryEntityByNameStmt = nullptr;
static sqlite3_stmt* g_queryEntityAllStmt = nullptr;
static sqlite3_stmt* g_queryEntityRowsByNameStmt = nullptr;
static sqlite3_stmt* g_queryEntityRowsAllStmt = nullptr;
static sqlite3_stmt* g_getEntityCatalogRowStmt = nullptr;
static sqlite3_stmt* g_getCatPayloadStmt = nullptr;
static sqlite3_stmt* g_getOmmPayloadByIdStmt = nullptr;
static std::map<uint32_t, PendingCatalogMetadata> g_pendingEntityMetadataByNorad;
static bool g_pendingEntityNamesFromFlatbuffer = false;

extern "C" int32_t plugin_init_omm(const OrbProOMMRecord* records, uint32_t count);
extern "C" double plugin_entity_add_omm(uint32_t entity_index, const OrbProOMMRecord* record);
extern "C" int32_t plugin_propagate(double julian_date, uint32_t entity_index, OrbProStateVector* out);
extern "C" int32_t plugin_get_entity_catalog_row(uint32_t entity_index, OrbProCatalogRow* out_row);
extern "C" int32_t plugin_get_omm_record_by_entity_index(uint32_t entity_index, OrbProOMMRecord* out_record);
extern "C" int32_t plugin_get_omm_record_by_pointer(double omm_pointer, OrbProOMMRecord* out_record);

// =============================================================================
// Internal Helpers
// =============================================================================

namespace {

bool ensureOmmDatabase() {
    if (
        g_ommDb != nullptr &&
        g_insertOmmStmt != nullptr &&
        g_upsertEntityStmt != nullptr &&
        g_upsertCatStmt != nullptr &&
        g_queryEntityByNameStmt != nullptr &&
        g_queryEntityAllStmt != nullptr &&
        g_queryEntityRowsByNameStmt != nullptr &&
        g_queryEntityRowsAllStmt != nullptr &&
        g_getEntityCatalogRowStmt != nullptr &&
        g_getCatPayloadStmt != nullptr &&
        g_getOmmPayloadByIdStmt != nullptr
    ) {
        return true;
    }

    if (g_ommDb == nullptr) {
        int openRc = sqlite3_open(":memory:", &g_ommDb);
        if (openRc != SQLITE_OK || g_ommDb == nullptr) {
            if (g_ommDb != nullptr) {
                sqlite3_close(g_ommDb);
                g_ommDb = nullptr;
            }
            return false;
        }
    }

    const char* createSql =
        "CREATE TABLE IF NOT EXISTS omm_records ("
        " id INTEGER PRIMARY KEY AUTOINCREMENT,"
        " norad_cat_id INTEGER NOT NULL,"
        " epoch_jd REAL NOT NULL,"
        " payload BLOB NOT NULL"
        ");"
        "CREATE INDEX IF NOT EXISTS idx_omm_records_norad_epoch "
        "ON omm_records(norad_cat_id, epoch_jd);"
        "CREATE TABLE IF NOT EXISTS entity_catalog ("
        " entity_index INTEGER PRIMARY KEY,"
        " norad_cat_id INTEGER NOT NULL,"
        " object_name TEXT NOT NULL,"
        " object_name_lc TEXT NOT NULL,"
        " object_id TEXT NOT NULL,"
        " object_id_lc TEXT NOT NULL"
        ");"
        "CREATE INDEX IF NOT EXISTS idx_entity_catalog_name_lc "
        "ON entity_catalog(object_name_lc);"
        "CREATE INDEX IF NOT EXISTS idx_entity_catalog_object_id_lc "
        "ON entity_catalog(object_id_lc);"
        "CREATE TABLE IF NOT EXISTS cat_catalog ("
        " norad_cat_id INTEGER PRIMARY KEY,"
        " object_name TEXT NOT NULL,"
        " object_name_lc TEXT NOT NULL,"
        " object_id TEXT NOT NULL,"
        " object_id_lc TEXT NOT NULL,"
        " payload_json TEXT NOT NULL"
        ");"
        "CREATE INDEX IF NOT EXISTS idx_cat_catalog_name_lc "
        "ON cat_catalog(object_name_lc);"
        "CREATE INDEX IF NOT EXISTS idx_cat_catalog_object_id_lc "
        "ON cat_catalog(object_id_lc);";

    int createRc = sqlite3_exec(g_ommDb, createSql, nullptr, nullptr, nullptr);
    if (createRc != SQLITE_OK) {
        return false;
    }

    if (g_insertOmmStmt == nullptr) {
        const char* insertSql =
            "INSERT INTO omm_records (norad_cat_id, epoch_jd, payload) "
            "VALUES (?1, ?2, ?3);";
        int prepRc = sqlite3_prepare_v2(g_ommDb, insertSql, -1, &g_insertOmmStmt, nullptr);
        if (prepRc != SQLITE_OK || g_insertOmmStmt == nullptr) {
            if (g_insertOmmStmt != nullptr) {
                sqlite3_finalize(g_insertOmmStmt);
                g_insertOmmStmt = nullptr;
            }
            return false;
        }
    }

    if (g_upsertEntityStmt == nullptr) {
        const char* upsertEntitySql =
            "INSERT INTO entity_catalog "
            "(entity_index, norad_cat_id, object_name, object_name_lc, object_id, object_id_lc) "
            "VALUES (?1, ?2, ?3, ?4, ?5, ?6) "
            "ON CONFLICT(entity_index) DO UPDATE SET "
            " norad_cat_id = excluded.norad_cat_id,"
            " object_name = excluded.object_name,"
            " object_name_lc = excluded.object_name_lc,"
            " object_id = excluded.object_id,"
            " object_id_lc = excluded.object_id_lc;";
        int prepRc = sqlite3_prepare_v2(g_ommDb, upsertEntitySql, -1, &g_upsertEntityStmt, nullptr);
        if (prepRc != SQLITE_OK || g_upsertEntityStmt == nullptr) {
            if (g_upsertEntityStmt != nullptr) {
                sqlite3_finalize(g_upsertEntityStmt);
                g_upsertEntityStmt = nullptr;
            }
            return false;
        }
    }

    if (g_upsertCatStmt == nullptr) {
        const char* upsertCatSql =
            "INSERT INTO cat_catalog "
            "(norad_cat_id, object_name, object_name_lc, object_id, object_id_lc, payload_json) "
            "VALUES (?1, ?2, ?3, ?4, ?5, ?6) "
            "ON CONFLICT(norad_cat_id) DO UPDATE SET "
            " object_name = excluded.object_name,"
            " object_name_lc = excluded.object_name_lc,"
            " object_id = excluded.object_id,"
            " object_id_lc = excluded.object_id_lc,"
            " payload_json = excluded.payload_json;";
        int prepRc = sqlite3_prepare_v2(g_ommDb, upsertCatSql, -1, &g_upsertCatStmt, nullptr);
        if (prepRc != SQLITE_OK || g_upsertCatStmt == nullptr) {
            if (g_upsertCatStmt != nullptr) {
                sqlite3_finalize(g_upsertCatStmt);
                g_upsertCatStmt = nullptr;
            }
            return false;
        }
    }

    if (g_queryEntityByNameStmt == nullptr) {
        const char* queryByNameSql =
            "SELECT ec.entity_index "
            "FROM entity_catalog ec "
            "LEFT JOIN cat_catalog cc ON cc.norad_cat_id = ec.norad_cat_id "
            "WHERE ec.object_name_lc LIKE ?1 "
            "   OR ec.object_id_lc LIKE ?1 "
            "   OR cc.object_name_lc LIKE ?1 "
            "   OR cc.object_id_lc LIKE ?1 "
            "   OR (?3 = 1 AND ec.norad_cat_id = ?2) "
            "ORDER BY ec.entity_index "
            "LIMIT ?4;";
        int prepRc = sqlite3_prepare_v2(g_ommDb, queryByNameSql, -1, &g_queryEntityByNameStmt, nullptr);
        if (prepRc != SQLITE_OK || g_queryEntityByNameStmt == nullptr) {
            if (g_queryEntityByNameStmt != nullptr) {
                sqlite3_finalize(g_queryEntityByNameStmt);
                g_queryEntityByNameStmt = nullptr;
            }
            return false;
        }
    }

    if (g_queryEntityRowsByNameStmt == nullptr) {
        const char* queryRowsByNameSql =
            "SELECT ec.entity_index, ec.norad_cat_id, ec.object_name, ec.object_id, "
            "       COALESCE(cc.object_name, ''), COALESCE(cc.object_id, '') "
            "FROM entity_catalog ec "
            "LEFT JOIN cat_catalog cc ON cc.norad_cat_id = ec.norad_cat_id "
            "WHERE ec.object_name_lc LIKE ?1 "
            "   OR ec.object_id_lc LIKE ?1 "
            "   OR cc.object_name_lc LIKE ?1 "
            "   OR cc.object_id_lc LIKE ?1 "
            "   OR (?3 = 1 AND ec.norad_cat_id = ?2) "
            "ORDER BY ec.entity_index "
            "LIMIT ?4;";
        int prepRc = sqlite3_prepare_v2(
            g_ommDb,
            queryRowsByNameSql,
            -1,
            &g_queryEntityRowsByNameStmt,
            nullptr
        );
        if (prepRc != SQLITE_OK || g_queryEntityRowsByNameStmt == nullptr) {
            if (g_queryEntityRowsByNameStmt != nullptr) {
                sqlite3_finalize(g_queryEntityRowsByNameStmt);
                g_queryEntityRowsByNameStmt = nullptr;
            }
            return false;
        }
    }

    if (g_queryEntityAllStmt == nullptr) {
        const char* queryAllSql =
            "SELECT entity_index FROM entity_catalog "
            "ORDER BY entity_index "
            "LIMIT ?1;";
        int prepRc = sqlite3_prepare_v2(g_ommDb, queryAllSql, -1, &g_queryEntityAllStmt, nullptr);
        if (prepRc != SQLITE_OK || g_queryEntityAllStmt == nullptr) {
            if (g_queryEntityAllStmt != nullptr) {
                sqlite3_finalize(g_queryEntityAllStmt);
                g_queryEntityAllStmt = nullptr;
            }
            return false;
        }
    }

    if (g_queryEntityRowsAllStmt == nullptr) {
        const char* queryRowsAllSql =
            "SELECT ec.entity_index, ec.norad_cat_id, ec.object_name, ec.object_id, "
            "       COALESCE(cc.object_name, ''), COALESCE(cc.object_id, '') "
            "FROM entity_catalog ec "
            "LEFT JOIN cat_catalog cc ON cc.norad_cat_id = ec.norad_cat_id "
            "ORDER BY ec.entity_index "
            "LIMIT ?1;";
        int prepRc = sqlite3_prepare_v2(
            g_ommDb,
            queryRowsAllSql,
            -1,
            &g_queryEntityRowsAllStmt,
            nullptr
        );
        if (prepRc != SQLITE_OK || g_queryEntityRowsAllStmt == nullptr) {
            if (g_queryEntityRowsAllStmt != nullptr) {
                sqlite3_finalize(g_queryEntityRowsAllStmt);
                g_queryEntityRowsAllStmt = nullptr;
            }
            return false;
        }
    }

    if (g_getEntityCatalogRowStmt == nullptr) {
        const char* getCatalogRowSql =
            "SELECT ec.entity_index, ec.norad_cat_id, ec.object_name, ec.object_id, "
            "       COALESCE(cc.object_name, ''), COALESCE(cc.object_id, '') "
            "FROM entity_catalog ec "
            "LEFT JOIN cat_catalog cc ON cc.norad_cat_id = ec.norad_cat_id "
            "WHERE ec.entity_index = ?1 "
            "LIMIT 1;";
        int prepRc = sqlite3_prepare_v2(
            g_ommDb,
            getCatalogRowSql,
            -1,
            &g_getEntityCatalogRowStmt,
            nullptr
        );
        if (prepRc != SQLITE_OK || g_getEntityCatalogRowStmt == nullptr) {
            if (g_getEntityCatalogRowStmt != nullptr) {
                sqlite3_finalize(g_getEntityCatalogRowStmt);
                g_getEntityCatalogRowStmt = nullptr;
            }
            return false;
        }
    }

    if (g_getCatPayloadStmt == nullptr) {
        const char* getCatPayloadSql =
            "SELECT payload_json FROM cat_catalog "
            "WHERE norad_cat_id = ?1 "
            "LIMIT 1;";
        int prepRc = sqlite3_prepare_v2(
            g_ommDb,
            getCatPayloadSql,
            -1,
            &g_getCatPayloadStmt,
            nullptr
        );
        if (prepRc != SQLITE_OK || g_getCatPayloadStmt == nullptr) {
            if (g_getCatPayloadStmt != nullptr) {
                sqlite3_finalize(g_getCatPayloadStmt);
                g_getCatPayloadStmt = nullptr;
            }
            return false;
        }
    }

    if (g_getOmmPayloadByIdStmt == nullptr) {
        const char* getOmmPayloadSql =
            "SELECT payload FROM omm_records "
            "WHERE id = ?1 "
            "LIMIT 1;";
        int prepRc = sqlite3_prepare_v2(
            g_ommDb,
            getOmmPayloadSql,
            -1,
            &g_getOmmPayloadByIdStmt,
            nullptr
        );
        if (prepRc != SQLITE_OK || g_getOmmPayloadByIdStmt == nullptr) {
            if (g_getOmmPayloadByIdStmt != nullptr) {
                sqlite3_finalize(g_getOmmPayloadByIdStmt);
                g_getOmmPayloadByIdStmt = nullptr;
            }
            return false;
        }
    }

    return true;
}

void clearOmmDatabase() {
    if (!ensureOmmDatabase()) {
        return;
    }
    sqlite3_exec(g_ommDb, "DELETE FROM omm_records;", nullptr, nullptr, nullptr);
    sqlite3_exec(g_ommDb, "DELETE FROM entity_catalog;", nullptr, nullptr, nullptr);
    sqlite3_exec(g_ommDb, "DELETE FROM cat_catalog;", nullptr, nullptr, nullptr);
}

void closeOmmDatabase() {
    if (g_insertOmmStmt != nullptr) {
        sqlite3_finalize(g_insertOmmStmt);
        g_insertOmmStmt = nullptr;
    }
    if (g_upsertEntityStmt != nullptr) {
        sqlite3_finalize(g_upsertEntityStmt);
        g_upsertEntityStmt = nullptr;
    }
    if (g_upsertCatStmt != nullptr) {
        sqlite3_finalize(g_upsertCatStmt);
        g_upsertCatStmt = nullptr;
    }
    if (g_queryEntityByNameStmt != nullptr) {
        sqlite3_finalize(g_queryEntityByNameStmt);
        g_queryEntityByNameStmt = nullptr;
    }
    if (g_queryEntityAllStmt != nullptr) {
        sqlite3_finalize(g_queryEntityAllStmt);
        g_queryEntityAllStmt = nullptr;
    }
    if (g_queryEntityRowsByNameStmt != nullptr) {
        sqlite3_finalize(g_queryEntityRowsByNameStmt);
        g_queryEntityRowsByNameStmt = nullptr;
    }
    if (g_queryEntityRowsAllStmt != nullptr) {
        sqlite3_finalize(g_queryEntityRowsAllStmt);
        g_queryEntityRowsAllStmt = nullptr;
    }
    if (g_getEntityCatalogRowStmt != nullptr) {
        sqlite3_finalize(g_getEntityCatalogRowStmt);
        g_getEntityCatalogRowStmt = nullptr;
    }
    if (g_getCatPayloadStmt != nullptr) {
        sqlite3_finalize(g_getCatPayloadStmt);
        g_getCatPayloadStmt = nullptr;
    }
    if (g_getOmmPayloadByIdStmt != nullptr) {
        sqlite3_finalize(g_getOmmPayloadByIdStmt);
        g_getOmmPayloadByIdStmt = nullptr;
    }
    if (g_ommDb != nullptr) {
        sqlite3_close(g_ommDb);
        g_ommDb = nullptr;
    }
    g_pendingEntityMetadataByNorad.clear();
}

int64_t insertOmmRecord(const OrbProOMMRecord& omm, double epochJD) {
    if (!ensureOmmDatabase()) {
        return -1;
    }

    sqlite3_reset(g_insertOmmStmt);
    sqlite3_clear_bindings(g_insertOmmStmt);

    if (sqlite3_bind_int64(g_insertOmmStmt, 1, static_cast<sqlite3_int64>(omm.norad_cat_id)) != SQLITE_OK) {
        sqlite3_reset(g_insertOmmStmt);
        return -1;
    }
    if (sqlite3_bind_double(g_insertOmmStmt, 2, epochJD) != SQLITE_OK) {
        sqlite3_reset(g_insertOmmStmt);
        return -1;
    }
    if (sqlite3_bind_blob(g_insertOmmStmt, 3, &omm, sizeof(OrbProOMMRecord), SQLITE_TRANSIENT) != SQLITE_OK) {
        sqlite3_reset(g_insertOmmStmt);
        return -1;
    }

    int stepRc = sqlite3_step(g_insertOmmStmt);
    if (stepRc != SQLITE_DONE) {
        sqlite3_reset(g_insertOmmStmt);
        return -1;
    }

    int64_t rowId = static_cast<int64_t>(sqlite3_last_insert_rowid(g_ommDb));
    sqlite3_reset(g_insertOmmStmt);
    sqlite3_clear_bindings(g_insertOmmStmt);
    return rowId > 0 ? rowId : -1;
}

std::string toLowerAscii(const std::string& text) {
    std::string lowered;
    lowered.resize(text.size());
    for (size_t i = 0; i < text.size(); i++) {
        const unsigned char ch = static_cast<unsigned char>(text[i]);
        if (ch >= 'A' && ch <= 'Z') {
            lowered[i] = static_cast<char>(ch - 'A' + 'a');
        } else {
            lowered[i] = static_cast<char>(ch);
        }
    }
    return lowered;
}

std::string trimAsciiWhitespace(const char* text) {
    if (text == nullptr) {
        return "";
    }
    std::string value(text);
    size_t start = 0;
    while (start < value.size()) {
        const unsigned char ch = static_cast<unsigned char>(value[start]);
        if (ch != ' ' && ch != '\t' && ch != '\n' && ch != '\r' && ch != '\f' && ch != '\v') {
            break;
        }
        start++;
    }
    size_t end = value.size();
    while (end > start) {
        const unsigned char ch = static_cast<unsigned char>(value[end - 1]);
        if (ch != ' ' && ch != '\t' && ch != '\n' && ch != '\r' && ch != '\f' && ch != '\v') {
            break;
        }
        end--;
    }
    return value.substr(start, end - start);
}

bool tryParseStrictUint32(const std::string& text, uint32_t& out) {
    if (text.empty()) {
        return false;
    }
    uint64_t value = 0;
    for (size_t i = 0; i < text.size(); i++) {
        const unsigned char ch = static_cast<unsigned char>(text[i]);
        if (ch < '0' || ch > '9') {
            return false;
        }
        value = value * 10ULL + static_cast<uint64_t>(ch - '0');
        if (value > static_cast<uint64_t>(std::numeric_limits<uint32_t>::max())) {
            return false;
        }
    }
    out = static_cast<uint32_t>(value);
    return true;
}

std::string sanitizeCatalogText(const std::string& candidate) {
    std::string text = candidate;
    text.erase(std::remove_if(text.begin(), text.end(), [](unsigned char c) {
        return c == '\r' || c == '\n' || c == '\t';
    }), text.end());
    if (!text.empty()) {
        const size_t first = text.find_first_not_of(' ');
        if (first != std::string::npos) {
            const size_t last = text.find_last_not_of(' ');
            text = text.substr(first, last - first + 1);
        } else {
            text.clear();
        }
    }
    return text;
}

std::string sanitizeDisplayName(const std::string& candidate, uint32_t noradCatId, uint32_t entityIndex) {
    const std::string name = sanitizeCatalogText(candidate);
    if (!name.empty()) {
        return name;
    }
    if (noradCatId > 0) {
        return "Satellite " + std::to_string(entityIndex) + " (NORAD " + std::to_string(noradCatId) + ")";
    }
    return "Satellite " + std::to_string(entityIndex);
}

std::string sanitizeObjectId(const std::string& candidate) {
    return sanitizeCatalogText(candidate);
}

PendingCatalogMetadata resolvePendingCatalogMetadata(uint32_t noradCatId, uint32_t entityIndex) {
    PendingCatalogMetadata metadata;
    auto pending = g_pendingEntityMetadataByNorad.find(noradCatId);
    if (pending != g_pendingEntityMetadataByNorad.end()) {
        metadata.objectName = sanitizeDisplayName(pending->second.objectName, noradCatId, entityIndex);
        metadata.objectId = sanitizeObjectId(pending->second.objectId);
        return metadata;
    }
    metadata.objectName = sanitizeDisplayName("", noradCatId, entityIndex);
    metadata.objectId = "";
    return metadata;
}

bool upsertEntityCatalog(
    uint32_t entityIndex,
    uint32_t noradCatId,
    const std::string& objectName,
    const std::string& objectId
) {
    if (!ensureOmmDatabase()) {
        return false;
    }
    if (g_upsertEntityStmt == nullptr) {
        return false;
    }

    const std::string safeName = sanitizeDisplayName(objectName, noradCatId, entityIndex);
    const std::string lowerName = toLowerAscii(safeName);
    const std::string safeObjectId = sanitizeObjectId(objectId);
    const std::string lowerObjectId = toLowerAscii(safeObjectId);

    sqlite3_reset(g_upsertEntityStmt);
    sqlite3_clear_bindings(g_upsertEntityStmt);

    if (sqlite3_bind_int(g_upsertEntityStmt, 1, static_cast<int>(entityIndex)) != SQLITE_OK) {
        sqlite3_reset(g_upsertEntityStmt);
        return false;
    }
    if (sqlite3_bind_int64(g_upsertEntityStmt, 2, static_cast<sqlite3_int64>(noradCatId)) != SQLITE_OK) {
        sqlite3_reset(g_upsertEntityStmt);
        return false;
    }
    if (sqlite3_bind_text(g_upsertEntityStmt, 3, safeName.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK) {
        sqlite3_reset(g_upsertEntityStmt);
        return false;
    }
    if (sqlite3_bind_text(g_upsertEntityStmt, 4, lowerName.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK) {
        sqlite3_reset(g_upsertEntityStmt);
        return false;
    }
    if (sqlite3_bind_text(g_upsertEntityStmt, 5, safeObjectId.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK) {
        sqlite3_reset(g_upsertEntityStmt);
        return false;
    }
    if (sqlite3_bind_text(g_upsertEntityStmt, 6, lowerObjectId.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK) {
        sqlite3_reset(g_upsertEntityStmt);
        return false;
    }

    const int rc = sqlite3_step(g_upsertEntityStmt);
    sqlite3_reset(g_upsertEntityStmt);
    sqlite3_clear_bindings(g_upsertEntityStmt);
    return rc == SQLITE_DONE;
}

std::string sanitizeCatalogName(const std::string& candidate, uint32_t noradCatId) {
    std::string name = sanitizeCatalogText(candidate);
    if (!name.empty()) {
        return name;
    }
    if (noradCatId > 0) {
        return "NORAD " + std::to_string(noradCatId);
    }
    return "";
}

std::string escapeJsonString(const std::string& input) {
    std::string escaped;
    escaped.reserve(input.size() + 8);
    for (size_t i = 0; i < input.size(); i++) {
        const unsigned char c = static_cast<unsigned char>(input[i]);
        switch (c) {
            case '\"': escaped += "\\\""; break;
            case '\\': escaped += "\\\\"; break;
            case '\b': escaped += "\\b"; break;
            case '\f': escaped += "\\f"; break;
            case '\n': escaped += "\\n"; break;
            case '\r': escaped += "\\r"; break;
            case '\t': escaped += "\\t"; break;
            default:
                if (c < 0x20) {
                    char hex[7];
                    std::snprintf(hex, sizeof(hex), "\\u%04x", static_cast<unsigned int>(c));
                    escaped += hex;
                } else {
                    escaped.push_back(static_cast<char>(c));
                }
                break;
        }
    }
    return escaped;
}

bool upsertCatCatalog(
    uint32_t noradCatId,
    const std::string& objectName,
    const std::string& objectId,
    const std::string& payloadJson
) {
    if (!ensureOmmDatabase() || g_upsertCatStmt == nullptr || noradCatId == 0) {
        return false;
    }

    const std::string safeName = sanitizeCatalogName(objectName, noradCatId);
    const std::string safeObjectId = sanitizeObjectId(objectId);
    const std::string safeNameJson = escapeJsonString(safeName);
    const std::string safeObjectIdJson = escapeJsonString(safeObjectId);
    const std::string safePayload =
        payloadJson.empty()
            ? ("{\"NORAD_CAT_ID\":" + std::to_string(noradCatId) +
              ",\"OBJECT_NAME\":\"" + safeNameJson + "\"" +
              ",\"OBJECT_ID\":\"" + safeObjectIdJson + "\"}")
            : payloadJson;
    const std::string lowerName = toLowerAscii(safeName);
    const std::string lowerObjectId = toLowerAscii(safeObjectId);

    sqlite3_reset(g_upsertCatStmt);
    sqlite3_clear_bindings(g_upsertCatStmt);

    if (sqlite3_bind_int64(g_upsertCatStmt, 1, static_cast<sqlite3_int64>(noradCatId)) != SQLITE_OK) {
        sqlite3_reset(g_upsertCatStmt);
        return false;
    }
    if (sqlite3_bind_text(g_upsertCatStmt, 2, safeName.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK) {
        sqlite3_reset(g_upsertCatStmt);
        return false;
    }
    if (sqlite3_bind_text(g_upsertCatStmt, 3, lowerName.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK) {
        sqlite3_reset(g_upsertCatStmt);
        return false;
    }
    if (sqlite3_bind_text(g_upsertCatStmt, 4, safeObjectId.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK) {
        sqlite3_reset(g_upsertCatStmt);
        return false;
    }
    if (sqlite3_bind_text(g_upsertCatStmt, 5, lowerObjectId.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK) {
        sqlite3_reset(g_upsertCatStmt);
        return false;
    }
    if (sqlite3_bind_text(g_upsertCatStmt, 6, safePayload.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK) {
        sqlite3_reset(g_upsertCatStmt);
        return false;
    }

    const int rc = sqlite3_step(g_upsertCatStmt);
    sqlite3_reset(g_upsertCatStmt);
    sqlite3_clear_bindings(g_upsertCatStmt);
    return rc == SQLITE_DONE;
}

void writeFixedString(char* destination, size_t capacity, const char* text) {
    if (destination == nullptr || capacity == 0) {
        return;
    }
    destination[0] = '\0';
    if (text == nullptr) {
        return;
    }
    size_t i = 0;
    for (; i + 1 < capacity && text[i] != '\0'; i++) {
        destination[i] = text[i];
    }
    destination[i] = '\0';
}

bool bindNameQueryStatement(
    sqlite3_stmt* stmt,
    const std::string& queryText,
    uint32_t maxCount
) {
    if (stmt == nullptr) {
        return false;
    }
    const int sqliteLimit = maxCount > 0x7fffffffU ? 0x7fffffff : static_cast<int>(maxCount);
    const std::string loweredQuery = toLowerAscii(queryText);
    const std::string pattern = "%" + loweredQuery + "%";
    uint32_t parsedNorad = 0;
    const bool hasNumericNorad = tryParseStrictUint32(queryText, parsedNorad);

    sqlite3_reset(stmt);
    sqlite3_clear_bindings(stmt);
    if (sqlite3_bind_text(stmt, 1, pattern.c_str(), -1, SQLITE_TRANSIENT) != SQLITE_OK) {
        sqlite3_reset(stmt);
        sqlite3_clear_bindings(stmt);
        return false;
    }
    if (sqlite3_bind_int64(stmt, 2, hasNumericNorad ? static_cast<sqlite3_int64>(parsedNorad) : 0) != SQLITE_OK) {
        sqlite3_reset(stmt);
        sqlite3_clear_bindings(stmt);
        return false;
    }
    if (sqlite3_bind_int(stmt, 3, hasNumericNorad ? 1 : 0) != SQLITE_OK) {
        sqlite3_reset(stmt);
        sqlite3_clear_bindings(stmt);
        return false;
    }
    if (sqlite3_bind_int(stmt, 4, sqliteLimit) != SQLITE_OK) {
        sqlite3_reset(stmt);
        sqlite3_clear_bindings(stmt);
        return false;
    }

    return true;
}

void fillCatalogRowFromSqlite(sqlite3_stmt* stmt, OrbProCatalogRow& row) {
    std::memset(&row, 0, sizeof(OrbProCatalogRow));
    row.entity_index = static_cast<uint32_t>(sqlite3_column_int(stmt, 0));
    row.norad_cat_id = static_cast<uint32_t>(sqlite3_column_int64(stmt, 1));
    writeFixedString(
        row.object_name,
        ORBPRO_NAME_BUFFER_BYTES,
        reinterpret_cast<const char*>(sqlite3_column_text(stmt, 2))
    );
    writeFixedString(
        row.object_id,
        ORBPRO_NAME_BUFFER_BYTES,
        reinterpret_cast<const char*>(sqlite3_column_text(stmt, 3))
    );
    writeFixedString(
        row.cat_object_name,
        ORBPRO_NAME_BUFFER_BYTES,
        reinterpret_cast<const char*>(sqlite3_column_text(stmt, 4))
    );
    writeFixedString(
        row.cat_object_id,
        ORBPRO_NAME_BUFFER_BYTES,
        reinterpret_cast<const char*>(sqlite3_column_text(stmt, 5))
    );
}

bool readCatPayloadJson(uint32_t noradCatId, std::string& jsonOut) {
    jsonOut.clear();
    if (!ensureOmmDatabase() || g_getCatPayloadStmt == nullptr || noradCatId == 0) {
        return false;
    }
    sqlite3_stmt* stmt = g_getCatPayloadStmt;
    sqlite3_reset(stmt);
    sqlite3_clear_bindings(stmt);
    if (sqlite3_bind_int64(stmt, 1, static_cast<sqlite3_int64>(noradCatId)) != SQLITE_OK) {
        sqlite3_reset(stmt);
        sqlite3_clear_bindings(stmt);
        return false;
    }

    const int rc = sqlite3_step(stmt);
    if (rc == SQLITE_ROW) {
        const unsigned char* payload = sqlite3_column_text(stmt, 0);
        if (payload != nullptr) {
            jsonOut.assign(reinterpret_cast<const char*>(payload));
        }
    }

    sqlite3_reset(stmt);
    sqlite3_clear_bindings(stmt);
    return rc == SQLITE_ROW;
}

bool executeCatalogSql(const char* sql) {
    if (sql == nullptr || !ensureOmmDatabase()) {
        return false;
    }
    return sqlite3_exec(g_ommDb, sql, nullptr, nullptr, nullptr) == SQLITE_OK;
}

// -----------------------------------------------------------------------------
// FlatBuffer OMM parsing helpers (direct binary ingest into WASM/sqlite)
// -----------------------------------------------------------------------------

static constexpr uint16_t OMM_FB_OBJECT_NAME_VT_OFFSET = 10;
static constexpr uint16_t OMM_FB_OBJECT_ID_VT_OFFSET = 12;
static constexpr uint16_t OMM_FB_EPOCH_VT_OFFSET = 26;
static constexpr uint16_t OMM_FB_MEAN_MOTION_VT_OFFSET = 30;
static constexpr uint16_t OMM_FB_ECCENTRICITY_VT_OFFSET = 32;
static constexpr uint16_t OMM_FB_INCLINATION_VT_OFFSET = 34;
static constexpr uint16_t OMM_FB_RA_ASC_NODE_VT_OFFSET = 36;
static constexpr uint16_t OMM_FB_ARG_PERICENTER_VT_OFFSET = 38;
static constexpr uint16_t OMM_FB_MEAN_ANOMALY_VT_OFFSET = 40;
static constexpr uint16_t OMM_FB_NORAD_CAT_ID_VT_OFFSET = 58;
static constexpr uint16_t OMM_FB_BSTAR_VT_OFFSET = 64;
static constexpr uint16_t OMM_FB_MEAN_MOTION_DOT_VT_OFFSET = 66;
static constexpr uint16_t OMM_FB_MEAN_MOTION_DDOT_VT_OFFSET = 68;
static constexpr uint16_t CAT_FB_OBJECT_NAME_VT_OFFSET = 4;
static constexpr uint16_t CAT_FB_OBJECT_ID_VT_OFFSET = 6;
static constexpr uint16_t CAT_FB_NORAD_CAT_ID_VT_OFFSET = 8;
static constexpr uint16_t REC_FB_RECORDS_VT_OFFSET = 6;
static constexpr uint16_t REC_RECORD_FB_VALUE_VT_OFFSET = 6;
static constexpr uint16_t REC_RECORD_FB_STANDARD_VT_OFFSET = 8;
static constexpr double J2000_JULIAN_DATE = 2451545.0;
static constexpr double MILLISECONDS_PER_DAY = 86400000.0;

inline bool readU16LE(const uint8_t* data, size_t len, size_t offset, uint16_t& out) {
    if (offset + 2 > len) return false;
    out = static_cast<uint16_t>(
        static_cast<uint16_t>(data[offset]) |
        (static_cast<uint16_t>(data[offset + 1]) << 8)
    );
    return true;
}

inline bool readU32LE(const uint8_t* data, size_t len, size_t offset, uint32_t& out) {
    if (offset + 4 > len) return false;
    out =
        static_cast<uint32_t>(data[offset]) |
        (static_cast<uint32_t>(data[offset + 1]) << 8) |
        (static_cast<uint32_t>(data[offset + 2]) << 16) |
        (static_cast<uint32_t>(data[offset + 3]) << 24);
    return true;
}

inline bool readI32LE(const uint8_t* data, size_t len, size_t offset, int32_t& out) {
    uint32_t u = 0;
    if (!readU32LE(data, len, offset, u)) return false;
    out = static_cast<int32_t>(u);
    return true;
}

inline bool readF64LE(const uint8_t* data, size_t len, size_t offset, double& out) {
    if (offset + 8 > len) return false;
    std::memcpy(&out, data + offset, sizeof(double));
    return true;
}

bool getFlatBufferTableBounds(
    const uint8_t* payload,
    size_t payloadLen,
    size_t& tablePos,
    size_t& vtablePos,
    uint16_t& vtableLen
) {
    int32_t rootOffset = 0;
    if (!readI32LE(payload, payloadLen, 0, rootOffset)) return false;
    if (rootOffset <= 0) return false;
    tablePos = static_cast<size_t>(rootOffset);
    if (tablePos + 4 > payloadLen) return false;

    int32_t vtableDistance = 0;
    if (!readI32LE(payload, payloadLen, tablePos, vtableDistance)) return false;
    if (vtableDistance <= 0 || static_cast<size_t>(vtableDistance) > tablePos) return false;
    vtablePos = tablePos - static_cast<size_t>(vtableDistance);

    if (!readU16LE(payload, payloadLen, vtablePos, vtableLen)) return false;
    return vtableLen >= 4 && vtablePos + static_cast<size_t>(vtableLen) <= payloadLen;
}

bool getFlatBufferTableMetadataAtPosition(
    const uint8_t* payload,
    size_t payloadLen,
    size_t tablePos,
    size_t& vtablePos,
    uint16_t& vtableLen
) {
    if (tablePos + 4 > payloadLen) return false;
    int32_t vtableDistance = 0;
    if (!readI32LE(payload, payloadLen, tablePos, vtableDistance)) return false;
    if (vtableDistance <= 0 || static_cast<size_t>(vtableDistance) > tablePos) return false;
    vtablePos = tablePos - static_cast<size_t>(vtableDistance);
    if (!readU16LE(payload, payloadLen, vtablePos, vtableLen)) return false;
    return vtableLen >= 4 && vtablePos + static_cast<size_t>(vtableLen) <= payloadLen;
}

bool getFlatBufferFieldAbsolutePos(
    const uint8_t* payload,
    size_t payloadLen,
    size_t tablePos,
    size_t vtablePos,
    uint16_t vtableLen,
    uint16_t fieldVtableOffset,
    size_t& fieldAbsPos
) {
    if (fieldVtableOffset + 2 > vtableLen) return false;
    uint16_t fieldRelOffset = 0;
    if (!readU16LE(payload, payloadLen, vtablePos + fieldVtableOffset, fieldRelOffset)) return false;
    if (fieldRelOffset == 0) return false;
    fieldAbsPos = tablePos + static_cast<size_t>(fieldRelOffset);
    return fieldAbsPos < payloadLen;
}

bool readFlatBufferF64Field(
    const uint8_t* payload,
    size_t payloadLen,
    size_t tablePos,
    size_t vtablePos,
    uint16_t vtableLen,
    uint16_t fieldVtableOffset,
    double defaultValue,
    double& out
) {
    size_t fieldAbsPos = 0;
    if (!getFlatBufferFieldAbsolutePos(
            payload, payloadLen, tablePos, vtablePos, vtableLen, fieldVtableOffset, fieldAbsPos)) {
        out = defaultValue;
        return true;
    }
    return readF64LE(payload, payloadLen, fieldAbsPos, out);
}

bool readFlatBufferU32Field(
    const uint8_t* payload,
    size_t payloadLen,
    size_t tablePos,
    size_t vtablePos,
    uint16_t vtableLen,
    uint16_t fieldVtableOffset,
    uint32_t defaultValue,
    uint32_t& out
) {
    size_t fieldAbsPos = 0;
    if (!getFlatBufferFieldAbsolutePos(
            payload, payloadLen, tablePos, vtablePos, vtableLen, fieldVtableOffset, fieldAbsPos)) {
        out = defaultValue;
        return true;
    }
    return readU32LE(payload, payloadLen, fieldAbsPos, out);
}

bool readFlatBufferStringField(
    const uint8_t* payload,
    size_t payloadLen,
    size_t tablePos,
    size_t vtablePos,
    uint16_t vtableLen,
    uint16_t fieldVtableOffset,
    std::string& out
) {
    size_t fieldAbsPos = 0;
    if (!getFlatBufferFieldAbsolutePos(
            payload, payloadLen, tablePos, vtablePos, vtableLen, fieldVtableOffset, fieldAbsPos)) {
        return false;
    }
    uint32_t stringRel = 0;
    if (!readU32LE(payload, payloadLen, fieldAbsPos, stringRel)) return false;
    size_t stringPos = fieldAbsPos + static_cast<size_t>(stringRel);
    uint32_t stringLen = 0;
    if (!readU32LE(payload, payloadLen, stringPos, stringLen)) return false;
    size_t dataPos = stringPos + 4;
    if (dataPos + static_cast<size_t>(stringLen) > payloadLen) return false;
    out.assign(reinterpret_cast<const char*>(payload + dataPos), stringLen);
    return true;
}

bool getFlatBufferVectorInfo(
    const uint8_t* payload,
    size_t payloadLen,
    size_t tablePos,
    size_t vtablePos,
    uint16_t vtableLen,
    uint16_t fieldVtableOffset,
    size_t& vectorPos,
    size_t& vectorDataPos,
    uint32_t& vectorLength
) {
    size_t fieldAbsPos = 0;
    if (!getFlatBufferFieldAbsolutePos(
            payload, payloadLen, tablePos, vtablePos, vtableLen, fieldVtableOffset, fieldAbsPos)) {
        return false;
    }

    uint32_t vectorRel = 0;
    if (!readU32LE(payload, payloadLen, fieldAbsPos, vectorRel)) return false;
    vectorPos = fieldAbsPos + static_cast<size_t>(vectorRel);
    if (!readU32LE(payload, payloadLen, vectorPos, vectorLength)) return false;
    vectorDataPos = vectorPos + 4;
    return vectorDataPos + static_cast<size_t>(vectorLength) * 4 <= payloadLen;
}

bool getFlatBufferIndirectTableField(
    const uint8_t* payload,
    size_t payloadLen,
    size_t tablePos,
    size_t vtablePos,
    uint16_t vtableLen,
    uint16_t fieldVtableOffset,
    size_t& nestedTablePos,
    size_t& nestedVtablePos,
    uint16_t& nestedVtableLen
) {
    size_t fieldAbsPos = 0;
    if (!getFlatBufferFieldAbsolutePos(
            payload, payloadLen, tablePos, vtablePos, vtableLen, fieldVtableOffset, fieldAbsPos)) {
        return false;
    }

    uint32_t nestedRel = 0;
    if (!readU32LE(payload, payloadLen, fieldAbsPos, nestedRel)) return false;
    nestedTablePos = fieldAbsPos + static_cast<size_t>(nestedRel);
    return getFlatBufferTableMetadataAtPosition(
        payload,
        payloadLen,
        nestedTablePos,
        nestedVtablePos,
        nestedVtableLen
    );
}

bool isOmmStandardName(const std::string& value) {
    return value == "OMM" || value == "$OMM" || value == "orbpro.sds.omm";
}

bool isCatStandardName(const std::string& value) {
    return value == "CAT" || value == "$CAT" || value == "orbpro.sds.cat";
}

bool parseEpochIso8601ToJulianDate(const std::string& epochText, double& julianDate) {
    int year = 0, mon = 0, day = 0, hr = 0, minute = 0;
    double sec = 0.0;

    int matched = std::sscanf(
        epochText.c_str(),
        "%d-%d-%dT%d:%d:%lf",
        &year,
        &mon,
        &day,
        &hr,
        &minute,
        &sec
    );
    if (matched < 6) {
        matched = std::sscanf(
            epochText.c_str(),
            "%d-%d-%d %d:%d:%lf",
            &year,
            &mon,
            &day,
            &hr,
            &minute,
            &sec
        );
    }
    if (matched < 6) return false;

    double jd = 0.0;
    double jdFrac = 0.0;
    SGP4Funcs::jday_SGP4(year, mon, day, hr, minute, sec, jd, jdFrac);
    julianDate = jd + jdFrac;
    return std::isfinite(julianDate);
}

bool parseOmmFlatBufferTable(
    const uint8_t* payload,
    size_t payloadLen,
    size_t tablePos,
    size_t vtablePos,
    uint16_t vtableLen,
    OrbProOMMRecord& out,
    PendingCatalogMetadata* catalogMetadataOut = nullptr
) {
    if (!payload || payloadLen < 8) {
        return false;
    }

    std::memset(&out, 0, sizeof(OrbProOMMRecord));

    std::string epochText;
    if (!readFlatBufferStringField(
            payload,
            payloadLen,
            tablePos,
            vtablePos,
            vtableLen,
            OMM_FB_EPOCH_VT_OFFSET,
            epochText)) {
        return false;
    }
    if (!parseEpochIso8601ToJulianDate(epochText, out.epoch_jd)) {
        return false;
    }

    if (catalogMetadataOut != nullptr) {
        std::string objectName;
        if (readFlatBufferStringField(
                payload,
                payloadLen,
                tablePos,
                vtablePos,
                vtableLen,
                OMM_FB_OBJECT_NAME_VT_OFFSET,
                objectName)) {
            catalogMetadataOut->objectName = objectName;
        } else {
            catalogMetadataOut->objectName.clear();
        }

        std::string objectId;
        if (readFlatBufferStringField(
                payload,
                payloadLen,
                tablePos,
                vtablePos,
                vtableLen,
                OMM_FB_OBJECT_ID_VT_OFFSET,
                objectId)) {
            catalogMetadataOut->objectId = objectId;
        } else {
            catalogMetadataOut->objectId.clear();
        }
    }

    if (!readFlatBufferF64Field(payload, payloadLen, tablePos, vtablePos, vtableLen, OMM_FB_MEAN_MOTION_VT_OFFSET, 0.0, out.mean_motion)) return false;
    if (!readFlatBufferF64Field(payload, payloadLen, tablePos, vtablePos, vtableLen, OMM_FB_ECCENTRICITY_VT_OFFSET, 0.0, out.eccentricity)) return false;
    if (!readFlatBufferF64Field(payload, payloadLen, tablePos, vtablePos, vtableLen, OMM_FB_INCLINATION_VT_OFFSET, 0.0, out.inclination)) return false;
    if (!readFlatBufferF64Field(payload, payloadLen, tablePos, vtablePos, vtableLen, OMM_FB_RA_ASC_NODE_VT_OFFSET, 0.0, out.ra_of_asc_node)) return false;
    if (!readFlatBufferF64Field(payload, payloadLen, tablePos, vtablePos, vtableLen, OMM_FB_ARG_PERICENTER_VT_OFFSET, 0.0, out.arg_of_pericenter)) return false;
    if (!readFlatBufferF64Field(payload, payloadLen, tablePos, vtablePos, vtableLen, OMM_FB_MEAN_ANOMALY_VT_OFFSET, 0.0, out.mean_anomaly)) return false;
    if (!readFlatBufferF64Field(payload, payloadLen, tablePos, vtablePos, vtableLen, OMM_FB_BSTAR_VT_OFFSET, 0.0, out.bstar)) return false;
    if (!readFlatBufferF64Field(payload, payloadLen, tablePos, vtablePos, vtableLen, OMM_FB_MEAN_MOTION_DOT_VT_OFFSET, 0.0, out.mean_motion_dot)) return false;
    if (!readFlatBufferF64Field(payload, payloadLen, tablePos, vtablePos, vtableLen, OMM_FB_MEAN_MOTION_DDOT_VT_OFFSET, 0.0, out.mean_motion_ddot)) return false;
    if (!readFlatBufferU32Field(payload, payloadLen, tablePos, vtablePos, vtableLen, OMM_FB_NORAD_CAT_ID_VT_OFFSET, 0, out.norad_cat_id)) return false;
    out._pad = 0;

    return true;
}

bool parseOmmFlatBufferPayload(
    const uint8_t* payload,
    size_t payloadLen,
    OrbProOMMRecord& out,
    PendingCatalogMetadata* catalogMetadataOut = nullptr
) {
    if (!payload || payloadLen < 8) return false;
    if (payload[4] != '$' || payload[5] != 'O' || payload[6] != 'M' || payload[7] != 'M') {
        return false;
    }

    size_t tablePos = 0;
    size_t vtablePos = 0;
    uint16_t vtableLen = 0;
    if (!getFlatBufferTableBounds(payload, payloadLen, tablePos, vtablePos, vtableLen)) {
        return false;
    }

    return parseOmmFlatBufferTable(
        payload,
        payloadLen,
        tablePos,
        vtablePos,
        vtableLen,
        out,
        catalogMetadataOut
    );
}

bool parseCatFlatBufferTable(
    const uint8_t* payload,
    size_t payloadLen,
    size_t tablePos,
    size_t vtablePos,
    uint16_t vtableLen,
    PendingCatCatalogRecord& out
) {
    out = PendingCatCatalogRecord{};
    if (!readFlatBufferU32Field(
            payload,
            payloadLen,
            tablePos,
            vtablePos,
            vtableLen,
            CAT_FB_NORAD_CAT_ID_VT_OFFSET,
            0,
            out.noradCatId)) {
        return false;
    }
    if (out.noradCatId == 0) {
        return false;
    }

    if (!readFlatBufferStringField(
            payload,
            payloadLen,
            tablePos,
            vtablePos,
            vtableLen,
            CAT_FB_OBJECT_NAME_VT_OFFSET,
            out.objectName)) {
        out.objectName.clear();
    }
    if (!readFlatBufferStringField(
            payload,
            payloadLen,
            tablePos,
            vtablePos,
            vtableLen,
            CAT_FB_OBJECT_ID_VT_OFFSET,
            out.objectId)) {
        out.objectId.clear();
    }
    return true;
}

bool parseCatFlatBufferPayload(
    const uint8_t* payload,
    size_t payloadLen,
    PendingCatCatalogRecord& out
) {
    if (!payload || payloadLen < 8) return false;
    if (payload[4] != '$' || payload[5] != 'C' || payload[6] != 'A' || payload[7] != 'T') {
        return false;
    }

    size_t tablePos = 0;
    size_t vtablePos = 0;
    uint16_t vtableLen = 0;
    if (!getFlatBufferTableBounds(payload, payloadLen, tablePos, vtablePos, vtableLen)) {
        return false;
    }

    return parseCatFlatBufferTable(
        payload,
        payloadLen,
        tablePos,
        vtablePos,
        vtableLen,
        out
    );
}

bool appendParsedOmmRecord(
    const OrbProOMMRecord& record,
    const PendingCatalogMetadata& metadata,
    std::vector<OrbProOMMRecord>& recordsOut,
    std::map<uint32_t, PendingCatalogMetadata>& metadataOut
) {
    recordsOut.push_back(record);
    if (
        record.norad_cat_id > 0 &&
        (!metadata.objectName.empty() || !metadata.objectId.empty())
    ) {
        metadataOut[record.norad_cat_id] = metadata;
    }
    return true;
}

bool parseRecFlatBufferPayload(
    const uint8_t* payload,
    size_t payloadLen,
    std::vector<OrbProOMMRecord>& recordsOut,
    std::map<uint32_t, PendingCatalogMetadata>& metadataOut
) {
    if (!payload || payloadLen < 8) return false;
    if (payload[4] != '$' || payload[5] != 'R' || payload[6] != 'E' || payload[7] != 'C') {
        return false;
    }

    size_t recTablePos = 0;
    size_t recVtablePos = 0;
    uint16_t recVtableLen = 0;
    if (!getFlatBufferTableBounds(payload, payloadLen, recTablePos, recVtablePos, recVtableLen)) {
        return false;
    }

    size_t recordsVectorPos = 0;
    size_t recordsVectorDataPos = 0;
    uint32_t recordsLength = 0;
    if (!getFlatBufferVectorInfo(
            payload,
            payloadLen,
            recTablePos,
            recVtablePos,
            recVtableLen,
            REC_FB_RECORDS_VT_OFFSET,
            recordsVectorPos,
            recordsVectorDataPos,
            recordsLength)) {
        return false;
    }

    bool parsedAny = false;
    for (uint32_t i = 0; i < recordsLength; i++) {
        const size_t entryPos = recordsVectorDataPos + static_cast<size_t>(i) * 4;
        uint32_t recordRel = 0;
        if (!readU32LE(payload, payloadLen, entryPos, recordRel)) {
            return false;
        }

        const size_t recordTablePos = entryPos + static_cast<size_t>(recordRel);
        size_t recordVtablePos = 0;
        uint16_t recordVtableLen = 0;
        if (!getFlatBufferTableMetadataAtPosition(
                payload,
                payloadLen,
                recordTablePos,
                recordVtablePos,
                recordVtableLen)) {
            return false;
        }

        std::string standard;
        if (!readFlatBufferStringField(
                payload,
                payloadLen,
                recordTablePos,
                recordVtablePos,
                recordVtableLen,
                REC_RECORD_FB_STANDARD_VT_OFFSET,
                standard) ||
            !isOmmStandardName(standard)) {
            continue;
        }

        size_t ommTablePos = 0;
        size_t ommVtablePos = 0;
        uint16_t ommVtableLen = 0;
        if (!getFlatBufferIndirectTableField(
                payload,
                payloadLen,
                recordTablePos,
                recordVtablePos,
                recordVtableLen,
                REC_RECORD_FB_VALUE_VT_OFFSET,
                ommTablePos,
                ommVtablePos,
                ommVtableLen)) {
            return false;
        }

        OrbProOMMRecord record;
        PendingCatalogMetadata metadata;
        if (!parseOmmFlatBufferTable(
                payload,
                payloadLen,
                ommTablePos,
                ommVtablePos,
                ommVtableLen,
                record,
                &metadata)) {
            return false;
        }

        appendParsedOmmRecord(record, metadata, recordsOut, metadataOut);
        parsedAny = true;
    }

    return parsedAny;
}

bool parseRecFlatBufferCatPayload(
    const uint8_t* payload,
    size_t payloadLen,
    std::vector<PendingCatCatalogRecord>& recordsOut
) {
    if (!payload || payloadLen < 8) return false;
    if (payload[4] != '$' || payload[5] != 'R' || payload[6] != 'E' || payload[7] != 'C') {
        return false;
    }

    size_t recTablePos = 0;
    size_t recVtablePos = 0;
    uint16_t recVtableLen = 0;
    if (!getFlatBufferTableBounds(payload, payloadLen, recTablePos, recVtablePos, recVtableLen)) {
        return false;
    }

    size_t recordsVectorPos = 0;
    size_t recordsVectorDataPos = 0;
    uint32_t recordsLength = 0;
    if (!getFlatBufferVectorInfo(
            payload,
            payloadLen,
            recTablePos,
            recVtablePos,
            recVtableLen,
            REC_FB_RECORDS_VT_OFFSET,
            recordsVectorPos,
            recordsVectorDataPos,
            recordsLength)) {
        return false;
    }

    bool parsedAny = false;
    for (uint32_t i = 0; i < recordsLength; i++) {
        const size_t entryPos = recordsVectorDataPos + static_cast<size_t>(i) * 4;
        uint32_t recordRel = 0;
        if (!readU32LE(payload, payloadLen, entryPos, recordRel)) {
            return false;
        }

        const size_t recordTablePos = entryPos + static_cast<size_t>(recordRel);
        size_t recordVtablePos = 0;
        uint16_t recordVtableLen = 0;
        if (!getFlatBufferTableMetadataAtPosition(
                payload,
                payloadLen,
                recordTablePos,
                recordVtablePos,
                recordVtableLen)) {
            return false;
        }

        std::string standard;
        if (!readFlatBufferStringField(
                payload,
                payloadLen,
                recordTablePos,
                recordVtablePos,
                recordVtableLen,
                REC_RECORD_FB_STANDARD_VT_OFFSET,
                standard) ||
            !isCatStandardName(standard)) {
            continue;
        }

        size_t catTablePos = 0;
        size_t catVtablePos = 0;
        uint16_t catVtableLen = 0;
        if (!getFlatBufferIndirectTableField(
                payload,
                payloadLen,
                recordTablePos,
                recordVtablePos,
                recordVtableLen,
                REC_RECORD_FB_VALUE_VT_OFFSET,
                catTablePos,
                catVtablePos,
                catVtableLen)) {
            return false;
        }

        PendingCatCatalogRecord record;
        if (!parseCatFlatBufferTable(
                payload,
                payloadLen,
                catTablePos,
                catVtablePos,
                catVtableLen,
                record)) {
            return false;
        }
        recordsOut.push_back(record);
        parsedAny = true;
    }

    return parsedAny;
}

bool appendOmmRecordsFromFlatBufferPayload(
    const uint8_t* payload,
    size_t payloadLen,
    std::vector<OrbProOMMRecord>& recordsOut,
    std::map<uint32_t, PendingCatalogMetadata>& metadataOut
) {
    if (!payload || payloadLen < 8) return false;

    if (payload[4] == '$' && payload[5] == 'O' && payload[6] == 'M' && payload[7] == 'M') {
        OrbProOMMRecord record;
        PendingCatalogMetadata metadata;
        if (!parseOmmFlatBufferPayload(payload, payloadLen, record, &metadata)) {
            return false;
        }
        return appendParsedOmmRecord(record, metadata, recordsOut, metadataOut);
    }

    if (payload[4] == '$' && payload[5] == 'R' && payload[6] == 'E' && payload[7] == 'C') {
        return parseRecFlatBufferPayload(payload, payloadLen, recordsOut, metadataOut);
    }

    return false;
}

bool appendCatRecordsFromFlatBufferPayload(
    const uint8_t* payload,
    size_t payloadLen,
    std::vector<PendingCatCatalogRecord>& recordsOut
) {
    if (!payload || payloadLen < 8) return false;

    if (payload[4] == '$' && payload[5] == 'C' && payload[6] == 'A' && payload[7] == 'T') {
        PendingCatCatalogRecord record;
        if (!parseCatFlatBufferPayload(payload, payloadLen, record)) {
            return false;
        }
        recordsOut.push_back(record);
        return true;
    }

    if (payload[4] == '$' && payload[5] == 'R' && payload[6] == 'E' && payload[7] == 'C') {
        return parseRecFlatBufferCatPayload(payload, payloadLen, recordsOut);
    }

    return false;
}

bool parseCatFlatBufferStream(
    const uint8_t* data,
    size_t len,
    std::vector<PendingCatCatalogRecord>& recordsOut
) {
    recordsOut.clear();
    if (data == nullptr || len == 0) {
        return true;
    }

    if (len >= 8 && data[4] == '$') {
        std::vector<PendingCatCatalogRecord> directRecords;
        if (appendCatRecordsFromFlatBufferPayload(data, len, directRecords)) {
            recordsOut = std::move(directRecords);
            return true;
        }
    }

    size_t offset = 0;
    while (offset + 4 <= len) {
        uint32_t messageSize = 0;
        if (!readU32LE(data, len, offset, messageSize)) {
            return false;
        }
        if (messageSize == 0 || offset + 4 + static_cast<size_t>(messageSize) > len) {
            return false;
        }

        const uint8_t* payload = data + offset + 4;
        if (!appendCatRecordsFromFlatBufferPayload(
                payload,
                static_cast<size_t>(messageSize),
                recordsOut)) {
            return false;
        }

        offset += 4 + static_cast<size_t>(messageSize);
    }

    return offset == len;
}

bool parseOmmFlatBufferStream(
    const uint8_t* data,
    size_t len,
    std::vector<OrbProOMMRecord>& recordsOut
) {
    recordsOut.clear();
    g_pendingEntityMetadataByNorad.clear();
    g_pendingEntityNamesFromFlatbuffer = false;
    if (data == nullptr || len == 0) {
        g_pendingEntityNamesFromFlatbuffer = true;
        return true;
    }

    // Single `$OMM` or `$REC` payload.
    if (len >= 8 && data[4] == '$') {
        std::vector<OrbProOMMRecord> directRecords;
        std::map<uint32_t, PendingCatalogMetadata> directMetadata;
        if (appendOmmRecordsFromFlatBufferPayload(
                data,
                len,
                directRecords,
                directMetadata)) {
            recordsOut = std::move(directRecords);
            g_pendingEntityMetadataByNorad = std::move(directMetadata);
            g_pendingEntityNamesFromFlatbuffer = true;
            return true;
        }
    }

    // Size-prefixed stream of `$OMM` / `$REC` payloads.
    size_t offset = 0;
    while (offset + 4 <= len) {
        uint32_t messageSize = 0;
        if (!readU32LE(data, len, offset, messageSize)) {
            return false;
        }
        if (messageSize == 0 || offset + 4 + static_cast<size_t>(messageSize) > len) {
            return false;
        }

        const uint8_t* payload = data + offset + 4;
        if (!appendOmmRecordsFromFlatBufferPayload(
                payload,
                static_cast<size_t>(messageSize),
                recordsOut,
                g_pendingEntityMetadataByNorad)) {
            g_pendingEntityMetadataByNorad.clear();
            g_pendingEntityNamesFromFlatbuffer = false;
            return false;
        }

        offset += 4 + static_cast<size_t>(messageSize);
    }

    if (offset != len) {
        g_pendingEntityMetadataByNorad.clear();
        g_pendingEntityNamesFromFlatbuffer = false;
        return false;
    }

    if (recordsOut.empty()) {
        g_pendingEntityMetadataByNorad.clear();
        g_pendingEntityNamesFromFlatbuffer = false;
        return false;
    }

    g_pendingEntityNamesFromFlatbuffer = true;
    return true;
}

// Convert TEME position/velocity to ECEF via GMST rotation and output in meters.
void temeToEcef(const double posTeme[3], const double velTeme[3],
                double gmst,
                double posEcef[3], double velEcef[3]) {
    double cosG = cos(gmst);
    double sinG = sin(gmst);

    // Rotate position: R_z(-gmst) * r_teme, km → meters
    posEcef[0] = ( cosG * posTeme[0] + sinG * posTeme[1]) * 1000.0;
    posEcef[1] = (-sinG * posTeme[0] + cosG * posTeme[1]) * 1000.0;
    posEcef[2] = posTeme[2] * 1000.0;

    // Rotate velocity, account for Earth rotation
    double vRotX = ( cosG * velTeme[0] + sinG * velTeme[1]) * 1000.0;
    double vRotY = (-sinG * velTeme[0] + cosG * velTeme[1]) * 1000.0;
    double vRotZ = velTeme[2] * 1000.0;

    velEcef[0] = vRotX + OMEGA_EARTH * posEcef[1];
    velEcef[1] = vRotY - OMEGA_EARTH * posEcef[0];
    velEcef[2] = vRotZ;
}

// Write state vector to OrbProStateVector buffer (ECEF, meters)
void writeStateVector(const SatelliteEntity& entity, double epochJulian,
                      OrbProStateVector* out) {
    double gmst = SGP4Funcs::gstime_SGP4(epochJulian);

    double posEcef[3], velEcef[3];
    temeToEcef(entity.r, entity.v, gmst, posEcef, velEcef);

    out->epoch = epochJulian;
    out->position[0] = posEcef[0];
    out->position[1] = posEcef[1];
    out->position[2] = posEcef[2];
    out->velocity[0] = velEcef[0];
    out->velocity[1] = velEcef[1];
    out->velocity[2] = velEcef[2];
    out->reference_frame = ORBPRO_FRAME_ECEF;
    out->flags = entity.valid ? ORBPRO_STATE_VALID : 0;
}

// -------------------------------------------------------------------------
// Multi-OMM helpers — ported from WasmSource/main_old
// -------------------------------------------------------------------------

// Find the closest satrec epoch JD to a given julian_date.
// Uses std::map::lower_bound for O(log n) binary search.
// Returns the epoch JD key, or -1.0 if the map is empty.
double findSatrecIndex(const std::map<double, elsetrec>& satrecs, double jd) {
    if (satrecs.empty()) return -1.0;

    auto it = satrecs.lower_bound(jd);

    if (it == satrecs.end()) {
        return std::prev(it)->first;
    }
    if (it == satrecs.begin()) {
        return it->first;
    }

    auto prevIt = std::prev(it);
    return (jd - prevIt->first < it->first - jd) ? prevIt->first : it->first;
}

// Find two bracketing satrec epoch JDs around a given julian_date.
// Returns {before, after}. Either may be -1.0 if the query is outside range.
std::pair<double, double> findSatrecIndices(const std::map<double, elsetrec>& satrecs, double jd) {
    if (satrecs.empty()) return {-1.0, -1.0};

    auto it = satrecs.lower_bound(jd);

    if (it == satrecs.end()) {
        return {std::prev(it)->first, -1.0};
    }
    if (it == satrecs.begin()) {
        return {-1.0, it->first};
    }

    return {std::prev(it)->first, it->first};
}

// Initialize a satrec from a binary OrbProOMMRecord.
// Returns true on success. Sets jdsatepoch/jdsatepochF on the satrec.
bool initSatrecFromOMM(const OrbProOMMRecord& omm, elsetrec& satrec) {
    int year, mon, day, hr, minute;
    double sec;
    SGP4Funcs::invjday_SGP4(omm.epoch_jd, 0.0, year, mon, day, hr, minute, sec);

    double jd, jdFrac;
    SGP4Funcs::jday_SGP4(year, mon, day, hr, minute, sec, jd, jdFrac);

    double epoch = jd + jdFrac - 2433281.5;

    double xno     = omm.mean_motion     * TWOPI / 1440.0;
    double xndot   = omm.mean_motion_dot * TWOPI / (1440.0 * 1440.0);
    double xnddot  = omm.mean_motion_ddot * TWOPI / (1440.0 * 1440.0 * 1440.0);
    double xinclo  = omm.inclination      * DEG_TO_RAD;
    double xnodeo  = omm.ra_of_asc_node   * DEG_TO_RAD;
    double xargpo  = omm.arg_of_pericenter * DEG_TO_RAD;
    double xmo     = omm.mean_anomaly     * DEG_TO_RAD;

    char satn[6];
    snprintf(satn, sizeof(satn), "%05u", omm.norad_cat_id);

    bool initOk = SGP4Funcs::sgp4init(
        wgs72, 'i', satn, epoch,
        omm.bstar, xndot, xnddot,
        omm.eccentricity, xargpo, xinclo, xmo, xno, xnodeo,
        satrec
    );

    if (!initOk || satrec.error != 0) return false;

    satrec.jdsatepoch  = jd;
    satrec.jdsatepochF = jdFrac;

    // Validate with test propagation at epoch
    double r[3], v[3];
    bool propOk = SGP4Funcs::sgp4(satrec, 0.0, r, v);
    return propOk && satrec.error == 0;
}

// -------------------------------------------------------------------------
// Core propagation helper — handles single and multi-OMM entities.
//
// Writes TEME position (km) into r[3] and velocity (km/s) into v[3].
// Updates entity.satrec to the currently active satrec for metadata queries.
// Returns true on success.
// -------------------------------------------------------------------------
bool propagateEntityTEME(SatelliteEntity& entity, double julian_date,
                         double r[3], double v[3]) {
    // Fast path: 0 or 1 OMM — use entity.satrec directly (no map lookups,
    // no struct copies). This covers the common case of single-OMM entities
    // and keeps ~29K-satellite constellations at full speed.
    if (entity.satrecs.size() <= 1) {
        double tleEpoch = entity.satrec.jdsatepoch + entity.satrec.jdsatepochF;
        double tsince = (julian_date - tleEpoch) * 1440.0;
        bool ok = SGP4Funcs::sgp4(entity.satrec, tsince, r, v);
        return ok && entity.satrec.error == 0;
    }

    // Interpolated epoch mode: blend between two bracketing OMMs
    if (entity.mode == MODE_INTERPOLATED_EPOCH) {
        auto indices = findSatrecIndices(entity.satrecs, julian_date);

        if (indices.first < 0 || indices.second < 0) {
            // Outside bracket range — fall back to nearest
            double nearest = findSatrecIndex(entity.satrecs, julian_date);
            auto& sat = entity.satrecs[nearest];
            double tleEpoch = sat.jdsatepoch + sat.jdsatepochF;
            double tsince = (julian_date - tleEpoch) * 1440.0;
            bool ok = SGP4Funcs::sgp4(sat, tsince, r, v);
            if (ok && sat.error == 0) {
                entity.satrec = sat;
                entity.currentSatrecEpoch = nearest;
            }
            return ok && sat.error == 0;
        }

        // Propagate both bracketing satrecs to the query time
        auto& sat1 = entity.satrecs[indices.first];
        auto& sat2 = entity.satrecs[indices.second];

        double epoch1 = sat1.jdsatepoch + sat1.jdsatepochF;
        double epoch2 = sat2.jdsatepoch + sat2.jdsatepochF;

        double r1[3], v1[3], r2[3], v2[3];
        bool ok1 = SGP4Funcs::sgp4(sat1, (julian_date - epoch1) * 1440.0, r1, v1);
        bool ok2 = SGP4Funcs::sgp4(sat2, (julian_date - epoch2) * 1440.0, r2, v2);

        if (!ok1 || sat1.error != 0 || !ok2 || sat2.error != 0) {
            // Error — fall back to nearest
            double nearest = findSatrecIndex(entity.satrecs, julian_date);
            auto& sat = entity.satrecs[nearest];
            double tleEpoch = sat.jdsatepoch + sat.jdsatepochF;
            bool ok = SGP4Funcs::sgp4(sat, (julian_date - tleEpoch) * 1440.0, r, v);
            return ok && sat.error == 0;
        }

        // Linear interpolation between the two TEME states.
        // With 2-point natural cubic spline (as in old WasmSource), this is
        // mathematically identical. Each satrec is propagated TO the query
        // time, so r1/r2 are both at julian_date but from different TLE
        // epochs. The blend weight controls smooth handoff between TLEs.
        double t = (julian_date - indices.first) / (indices.second - indices.first);
        r[0] = r1[0] + t * (r2[0] - r1[0]);
        r[1] = r1[1] + t * (r2[1] - r1[1]);
        r[2] = r1[2] + t * (r2[2] - r1[2]);
        v[0] = v1[0] + t * (v2[0] - v1[0]);
        v[1] = v1[1] + t * (v2[1] - v1[1]);
        v[2] = v1[2] + t * (v2[2] - v1[2]);

        // Update entity.satrec to the nearest bracket for metadata queries
        if (julian_date - indices.first < indices.second - julian_date) {
            entity.satrec = sat1;
            entity.currentSatrecEpoch = indices.first;
        } else {
            entity.satrec = sat2;
            entity.currentSatrecEpoch = indices.second;
        }
        return true;
    }

    // Nearest epoch mode (default): find closest satrec and propagate from it
    double nearest = findSatrecIndex(entity.satrecs, julian_date);
    auto& sat = entity.satrecs[nearest];
    double tleEpoch = sat.jdsatepoch + sat.jdsatepochF;
    double tsince = (julian_date - tleEpoch) * 1440.0;
    bool ok = SGP4Funcs::sgp4(sat, tsince, r, v);
    if (ok && sat.error == 0) {
        entity.satrec = sat;
        entity.currentSatrecEpoch = nearest;
    }
    return ok && sat.error == 0;
}

int32_t findEntityIndexByNorad(uint32_t noradCatId) {
    for (uint32_t i = 0; i < static_cast<uint32_t>(g_satellites.size()); i++) {
        if (g_satellites[i].noradId == noradCatId) {
            return static_cast<int32_t>(i);
        }
    }
    return -1;
}

bool addNewEntityFromOmmRecord(
    const OrbProOMMRecord& omm,
    const PendingCatalogMetadata* metadata = nullptr
) {
    elsetrec satrec;
    if (!initSatrecFromOMM(omm, satrec)) {
        return false;
    }

    const double epochJD = satrec.jdsatepoch + satrec.jdsatepochF;
    const int64_t ommPointer = insertOmmRecord(omm, epochJD);
    if (ommPointer <= 0) {
        return false;
    }

    double r[3];
    double v[3];
    const bool propOk = SGP4Funcs::sgp4(satrec, 0.0, r, v);
    if (!propOk || satrec.error != 0) {
        return false;
    }

    SatelliteEntity entity;
    entity.noradId = omm.norad_cat_id;
    entity.valid = true;
    entity.lastEpochJD = epochJD;
    entity.r[0] = r[0]; entity.r[1] = r[1]; entity.r[2] = r[2];
    entity.v[0] = v[0]; entity.v[1] = v[1]; entity.v[2] = v[2];
    entity.currentSatrecEpoch = epochJD;
    entity.currentOmmPointer = ommPointer;
    entity.mode = MODE_NEAREST_EPOCH;
    entity.satrec = satrec;

    const uint32_t entityIndex = static_cast<uint32_t>(g_satellites.size());
    g_satellites.push_back(entity);

    const std::string objectName = metadata ? metadata->objectName : std::string();
    const std::string objectId = metadata ? metadata->objectId : std::string();
    return upsertEntityCatalog(entityIndex, omm.norad_cat_id, objectName, objectId);
}

bool ingestOmmRecords(
    const std::vector<OrbProOMMRecord>& records,
    const std::map<uint32_t, PendingCatalogMetadata>& metadataByNorad
) {
    if (records.empty()) {
        return true;
    }

    if (!g_initialized || g_satellites.empty()) {
        std::map<uint32_t, PendingCatalogMetadata> previousMetadata = g_pendingEntityMetadataByNorad;
        const bool previousPendingNames = g_pendingEntityNamesFromFlatbuffer;
        g_pendingEntityMetadataByNorad = metadataByNorad;
        g_pendingEntityNamesFromFlatbuffer = true;
        const int32_t initResult = plugin_init_omm(records.data(), static_cast<uint32_t>(records.size()));
        if (initResult <= 0) {
            g_pendingEntityMetadataByNorad = previousMetadata;
            g_pendingEntityNamesFromFlatbuffer = previousPendingNames;
            return false;
        }
        return true;
    }

    uint32_t appliedCount = 0;
    for (const auto& record : records) {
        const auto metadataIt = metadataByNorad.find(record.norad_cat_id);
        const PendingCatalogMetadata* metadata =
            metadataIt != metadataByNorad.end() ? &metadataIt->second : nullptr;

        const int32_t entityIndex = findEntityIndexByNorad(record.norad_cat_id);
        if (entityIndex >= 0) {
            const double epochJD = plugin_entity_add_omm(static_cast<uint32_t>(entityIndex), &record);
            if (epochJD < 0.0) {
                continue;
            }

            if (metadata != nullptr) {
                if (!upsertEntityCatalog(
                        static_cast<uint32_t>(entityIndex),
                        record.norad_cat_id,
                        metadata->objectName,
                        metadata->objectId)) {
                    return false;
                }
            }
            appliedCount += 1;
            continue;
        }

        if (!addNewEntityFromOmmRecord(record, metadata)) {
            continue;
        }
        appliedCount += 1;
    }

    g_initialized = !g_satellites.empty();
    return appliedCount > 0;
}

int64_t julianDateToJ2000Milliseconds(double julianDate) {
    const double milliseconds = (julianDate - J2000_JULIAN_DATE) * MILLISECONDS_PER_DAY;
    return static_cast<int64_t>(std::llround(milliseconds));
}

uint8_t* copyFlatBufferToHeap(const flatbuffers::FlatBufferBuilder& builder, uint32_t* sizeOut) {
    if (sizeOut != nullptr) {
        *sizeOut = static_cast<uint32_t>(builder.GetSize());
    }
    if (builder.GetSize() == 0) {
        return nullptr;
    }
    uint8_t* bytes = static_cast<uint8_t*>(orbpro_malloc(builder.GetSize()));
    if (bytes == nullptr) {
        if (sizeOut != nullptr) {
            *sizeOut = 0;
        }
        return nullptr;
    }
    std::memcpy(bytes, builder.GetBufferPointer(), builder.GetSize());
    return bytes;
}

uint8_t* buildErrorStreamInvokeResponse(
    int32_t errorCode,
    const char* errorMessage,
    uint32_t* responseSizeOut
) {
    flatbuffers::FlatBufferBuilder builder(256);
    const auto response = orbpro::plugin::CreateStreamInvokeResponseDirect(
        builder,
        nullptr,
        0,
        false,
        errorCode,
        errorMessage
    );
    builder.Finish(response);
    return copyFlatBufferToHeap(builder, responseSizeOut);
}

uint8_t* buildEmptyStreamInvokeResponse(uint32_t* responseSizeOut) {
    flatbuffers::FlatBufferBuilder builder(128);
    const auto response = orbpro::plugin::CreateStreamInvokeResponseDirect(
        builder,
        nullptr,
        0,
        false,
        0,
        nullptr
    );
    builder.Finish(response);
    return copyFlatBufferToHeap(builder, responseSizeOut);
}

bool decodePropagatorBatchRequest(
    const uint8_t* payload,
    size_t payloadSize,
    const orbpro::propagator::PropagatorBatchRequest*& requestOut
) {
    requestOut = nullptr;
    if (payload == nullptr || payloadSize == 0) {
        return false;
    }

    flatbuffers::Verifier verifier(payload, payloadSize);
    if (!verifier.VerifyBuffer<orbpro::propagator::PropagatorBatchRequest>(nullptr)) {
        return false;
    }

    requestOut = flatbuffers::GetRoot<orbpro::propagator::PropagatorBatchRequest>(payload);
    return requestOut != nullptr;
}

bool decodeDescribeSourcesBatchRequest(
    const uint8_t* payload,
    size_t payloadSize,
    const orbpro::propagator::PropagatorDescribeSourcesBatchRequest*& requestOut
) {
    requestOut = nullptr;
    if (payload == nullptr || payloadSize == 0) {
        return false;
    }

    flatbuffers::Verifier verifier(payload, payloadSize);
    if (!verifier.VerifyBuffer<orbpro::propagator::PropagatorDescribeSourcesBatchRequest>(nullptr)) {
        return false;
    }

    requestOut =
        flatbuffers::GetRoot<orbpro::propagator::PropagatorDescribeSourcesBatchRequest>(payload);
    return requestOut != nullptr;
}

bool decodePrepareTrajectorySegmentsRequest(
    const uint8_t* payload,
    size_t payloadSize,
    const orbpro::propagator::PropagatorPrepareTrajectorySegmentsRequest*& requestOut
) {
    requestOut = nullptr;
    if (payload == nullptr || payloadSize == 0) {
        return false;
    }

    flatbuffers::Verifier verifier(payload, payloadSize);
    if (!verifier.VerifyBuffer<orbpro::propagator::PropagatorPrepareTrajectorySegmentsRequest>(nullptr)) {
        return false;
    }

    requestOut = flatbuffers::GetRoot<orbpro::propagator::PropagatorPrepareTrajectorySegmentsRequest>(payload);
    return requestOut != nullptr;
}

bool decodeSampleTrajectoryStatesRequest(
    const uint8_t* payload,
    size_t payloadSize,
    const orbpro::propagator::PropagatorSampleTrajectoryStatesRequest*& requestOut
) {
    requestOut = nullptr;
    if (payload == nullptr || payloadSize == 0) {
        return false;
    }

    flatbuffers::Verifier verifier(payload, payloadSize);
    if (!verifier.VerifyBuffer<orbpro::propagator::PropagatorSampleTrajectoryStatesRequest>(nullptr)) {
        return false;
    }

    requestOut = flatbuffers::GetRoot<orbpro::propagator::PropagatorSampleTrajectoryStatesRequest>(payload);
    return requestOut != nullptr;
}

bool decodeDescribeTrajectorySegmentsRequest(
    const uint8_t* payload,
    size_t payloadSize,
    const orbpro::propagator::PropagatorDescribeTrajectorySegmentsRequest*& requestOut
) {
    requestOut = nullptr;
    if (payload == nullptr || payloadSize == 0) {
        return false;
    }

    flatbuffers::Verifier verifier(payload, payloadSize);
    if (!verifier.VerifyBuffer<orbpro::propagator::PropagatorDescribeTrajectorySegmentsRequest>(nullptr)) {
        return false;
    }

    requestOut = flatbuffers::GetRoot<orbpro::propagator::PropagatorDescribeTrajectorySegmentsRequest>(payload);
    return requestOut != nullptr;
}

bool decodeCatalogQueryRequest(
    const uint8_t* payload,
    size_t payloadSize,
    const orbpro::query::CatalogQueryRequest*& requestOut
) {
    requestOut = nullptr;
    if (payload == nullptr || payloadSize == 0) {
        return false;
    }

    flatbuffers::Verifier verifier(payload, payloadSize);
    if (!orbpro::query::VerifyCatalogQueryRequestBuffer(verifier)) {
        return false;
    }

    requestOut = orbpro::query::GetCatalogQueryRequest(payload);
    return requestOut != nullptr;
}

std::string readFixedCatalogText(const char* value, size_t capacity) {
    if (value == nullptr || capacity == 0) {
        return std::string();
    }

    size_t length = 0;
    while (length < capacity && value[length] != '\0') {
        length++;
    }
    return std::string(value, length);
}

std::string buildCatalogEntityId(const OrbProCatalogRow& row) {
    if (row.norad_cat_id != 0) {
        return std::string("sat-") + std::to_string(row.norad_cat_id);
    }
    return std::string("sgp4-entity-") + std::to_string(row.entity_index);
}

std::string buildCatalogSearchText(
    const std::string& entityId,
    const std::string& name,
    const std::string& objectName,
    const std::string& objectId,
    const std::string& catObjectName,
    const std::string& catObjectId,
    uint32_t noradCatId
) {
    std::string searchText;
    const auto appendToken = [&searchText](const std::string& token) {
        if (token.empty()) {
            return;
        }
        if (!searchText.empty()) {
            searchText.push_back(' ');
        }
        searchText += toLowerAscii(token);
    };

    appendToken(entityId);
    appendToken(name);
    appendToken(objectName);
    appendToken(objectId);
    appendToken(catObjectName);
    appendToken(catObjectId);
    if (noradCatId != 0) {
        appendToken(std::to_string(noradCatId));
    }
    return searchText;
}

flatbuffers::Offset<orbpro::entity::EntityMetadata> createCatalogEntityMetadata(
    flatbuffers::FlatBufferBuilder& builder,
    const OrbProCatalogRow& row
) {
    const std::string objectName = readFixedCatalogText(
        row.object_name,
        sizeof(row.object_name)
    );
    const std::string objectId = readFixedCatalogText(
        row.object_id,
        sizeof(row.object_id)
    );
    const std::string catObjectName = readFixedCatalogText(
        row.cat_object_name,
        sizeof(row.cat_object_name)
    );
    const std::string catObjectId = readFixedCatalogText(
        row.cat_object_id,
        sizeof(row.cat_object_id)
    );
    const std::string entityId = buildCatalogEntityId(row);
    const std::string name = !objectName.empty()
        ? objectName
        : (!catObjectName.empty() ? catObjectName : entityId);
    const std::string searchText = buildCatalogSearchText(
        entityId,
        name,
        objectName,
        objectId,
        catObjectName,
        catObjectId,
        row.norad_cat_id
    );

    return orbpro::entity::CreateEntityMetadataDirect(
        builder,
        entityId.c_str(),
        name.c_str(),
        orbpro::entity::EntityKind_SPACE,
        "SGP4Entity",
        nullptr,
        nullptr,
        0.0,
        row.entity_index,
        0,
        0,
        0,
        0,
        0,
        0,
        row.norad_cat_id,
        objectName.empty() ? nullptr : objectName.c_str(),
        objectId.empty() ? nullptr : objectId.c_str(),
        catObjectName.empty() ? nullptr : catObjectName.c_str(),
        catObjectId.empty() ? nullptr : catObjectId.c_str(),
        nullptr,
        searchText.empty() ? nullptr : searchText.c_str()
    );
}

uint8_t* encodeCatalogQueryResultPayload(
    orbpro::query::CatalogQueryKind queryKind,
    const std::vector<OrbProCatalogRow>* rows,
    const std::vector<uint32_t>* entityIndices,
    const std::vector<uint8_t>* mask,
    uint32_t visibleCount,
    uint32_t entityIndex,
    const OrbProCatalogRow* row,
    uint32_t* payloadSizeOut
) {
    flatbuffers::FlatBufferBuilder builder(1024);

    flatbuffers::Offset<flatbuffers::Vector<flatbuffers::Offset<orbpro::entity::EntityMetadata>>> rowsOffset = 0;
    if (rows != nullptr && !rows->empty()) {
        std::vector<flatbuffers::Offset<orbpro::entity::EntityMetadata>> rowOffsets;
        rowOffsets.reserve(rows->size());
        for (const OrbProCatalogRow& catalogRow : *rows) {
            rowOffsets.push_back(createCatalogEntityMetadata(builder, catalogRow));
        }
        rowsOffset = builder.CreateVector(rowOffsets);
    }

    const auto entityIndicesOffset =
        entityIndices != nullptr && !entityIndices->empty()
            ? builder.CreateVector(*entityIndices)
            : flatbuffers::Offset<flatbuffers::Vector<uint32_t>>(0);
    const auto maskOffset =
        mask != nullptr && !mask->empty()
            ? builder.CreateVector(*mask)
            : flatbuffers::Offset<flatbuffers::Vector<uint8_t>>(0);
    const auto rowOffset =
        row != nullptr
            ? createCatalogEntityMetadata(builder, *row)
            : flatbuffers::Offset<orbpro::entity::EntityMetadata>(0);

    const auto result = orbpro::query::CreateCatalogQueryResult(
        builder,
        queryKind,
        rowsOffset,
        entityIndicesOffset,
        maskOffset,
        visibleCount,
        entityIndex,
        rowOffset
    );
    orbpro::query::FinishCatalogQueryResultBuffer(builder, result);
    return copyFlatBufferToHeap(builder, payloadSizeOut);
}

uint8_t* encodePropagatorStatePayload(
    const OrbProStateVector& state,
    uint32_t entityIndex,
    uint32_t catalogNumber,
    bool valid,
    uint32_t* payloadSizeOut
) {
    flatbuffers::FlatBufferBuilder builder(256);
    const std::vector<double> position = {
        state.position[0],
        state.position[1],
        state.position[2],
    };
    const std::vector<double> velocity = {
        state.velocity[0],
        state.velocity[1],
        state.velocity[2],
    };

    const auto positionOffset = builder.CreateVector(position);
    const auto velocityOffset = builder.CreateVector(velocity);
    const auto stateOffset = orbpro::plugins::CreatePropagatorState(
        builder,
        positionOffset,
        velocityOffset,
        julianDateToJ2000Milliseconds(state.epoch),
        orbpro::plugins::ReferenceFrame_ECEF,
        0,
        0.0,
        0.0,
        catalogNumber,
        entityIndex,
        valid
    );
    orbpro::plugins::FinishPropagatorStateBuffer(builder, stateOffset);
    return copyFlatBufferToHeap(builder, payloadSizeOut);
}

bool resolveRequestedSatelliteHandles(
    const ::flatbuffers::Vector<uint32_t>* requestedHandles,
    std::vector<uint32_t>& handlesOut
) {
    handlesOut.clear();
    const uint32_t satelliteCount = static_cast<uint32_t>(g_satellites.size());
    if (requestedHandles != nullptr && requestedHandles->size() > 0) {
        handlesOut.reserve(requestedHandles->size());
        for (uint32_t handle : *requestedHandles) {
            if (handle >= satelliteCount) {
                return false;
            }
            handlesOut.push_back(handle);
        }
        return true;
    }

    handlesOut.reserve(satelliteCount);
    for (uint32_t handle = 0; handle < satelliteCount; handle++) {
        handlesOut.push_back(handle);
    }
    return true;
}

bool resolveRequestedSegmentHandles(
    const ::flatbuffers::Vector<uint32_t>* requestedHandles,
    const PreparedTrajectorySegmentSet& segmentSet,
    std::vector<uint32_t>& handlesOut
) {
    handlesOut.clear();
    if (requestedHandles != nullptr && requestedHandles->size() > 0) {
        handlesOut.reserve(requestedHandles->size());
        for (uint32_t handle : *requestedHandles) {
            if (std::find(segmentSet.sourceHandles.begin(), segmentSet.sourceHandles.end(), handle) ==
                segmentSet.sourceHandles.end()) {
                continue;
            }
            handlesOut.push_back(handle);
        }
        return true;
    }

    handlesOut = segmentSet.sourceHandles;
    return true;
}

double computePerigeeKm(const SatelliteEntity& satellite) {
    return satellite.satrec.altp * satellite.satrec.radiusearthkm;
}

double computeApogeeKm(const SatelliteEntity& satellite) {
    return satellite.satrec.alta * satellite.satrec.radiusearthkm;
}

OrbProOMMRecord buildFallbackDescribeOmmRecord(const SatelliteEntity& satellite) {
    OrbProOMMRecord record = {};
    record.epoch_jd = satellite.satrec.jdsatepoch + satellite.satrec.jdsatepochF;
    record.mean_motion = satellite.satrec.no_unkozai * (24.0 * 60.0) / TWOPI;
    record.eccentricity = satellite.satrec.ecco;
    record.inclination = satellite.satrec.inclo / DEG_TO_RAD;
    record.ra_of_asc_node = satellite.satrec.nodeo / DEG_TO_RAD;
    record.arg_of_pericenter = satellite.satrec.argpo / DEG_TO_RAD;
    record.mean_anomaly = satellite.satrec.mo / DEG_TO_RAD;
    record.bstar = satellite.satrec.bstar;
    record.mean_motion_dot = 0.0;
    record.mean_motion_ddot = 0.0;
    record.norad_cat_id = satellite.noradId;
    return record;
}

OrbProOMMRecord getDescribeOmmRecord(uint32_t handle, const SatelliteEntity& satellite) {
    OrbProOMMRecord record = {};
    if (
        satellite.currentOmmPointer > 0 &&
        plugin_get_omm_record_by_pointer(
            static_cast<double>(satellite.currentOmmPointer),
            &record
        ) == 0
    ) {
        return record;
    }

    if (!satellite.ommPointers.empty()) {
        auto pointerIt = satellite.ommPointers.find(satellite.currentSatrecEpoch);
        if (
            pointerIt != satellite.ommPointers.end() &&
            pointerIt->second > 0 &&
            plugin_get_omm_record_by_pointer(
                static_cast<double>(pointerIt->second),
                &record
            ) == 0
        ) {
            return record;
        }
    }

    return buildFallbackDescribeOmmRecord(satellite);
}

PendingCatalogMetadata getDescribeCatalogMetadata(uint32_t handle, const SatelliteEntity& satellite) {
    OrbProCatalogRow row;
    if (plugin_get_entity_catalog_row(handle, &row) == 0) {
        return PendingCatalogMetadata{
            sanitizeDisplayName(row.object_name, satellite.noradId, handle),
            sanitizeObjectId(row.object_id),
        };
    }

    return resolvePendingCatalogMetadata(satellite.noradId, handle);
}

uint8_t* encodeDescribeSourcesBatchResultPayload(
    uint32_t catalogHandle,
    const std::vector<uint32_t>& handles,
    uint32_t* payloadSizeOut
) {
    flatbuffers::FlatBufferBuilder builder(1024);
    std::vector<flatbuffers::Offset<orbpro::propagator::PropagatorSourceDescription>> descriptions;
    descriptions.reserve(handles.size());

    for (uint32_t handle : handles) {
        const SatelliteEntity& satellite = g_satellites[handle];
        const OrbProOMMRecord omm = getDescribeOmmRecord(handle, satellite);
        const PendingCatalogMetadata metadata =
            getDescribeCatalogMetadata(handle, satellite);
        const char classificationType =
            satellite.satrec.classification != '\0'
                ? satellite.satrec.classification
                : 'U';
        const std::string classificationValue(1, classificationType);
        descriptions.push_back(
            orbpro::propagator::CreatePropagatorSourceDescriptionDirect(
                builder,
                handle,
                orbpro::propagator::PropagatorSourceKind_SGP4,
                metadata.objectName.c_str(),
                metadata.objectId.c_str(),
                satellite.noradId,
                omm.epoch_jd,
                omm.mean_motion,
                omm.eccentricity,
                omm.inclination,
                omm.ra_of_asc_node,
                omm.arg_of_pericenter,
                omm.mean_anomaly,
                satellite.satrec.ephtype,
                classificationValue.c_str(),
                static_cast<uint32_t>(std::max<long>(0, satellite.satrec.elnum)),
                static_cast<uint32_t>(std::max<long>(0, satellite.satrec.revnum)),
                omm.bstar,
                omm.mean_motion_dot,
                omm.mean_motion_ddot,
                computePerigeeKm(satellite),
                computeApogeeKm(satellite)
            )
        );
    }

    const auto descriptionsOffset = builder.CreateVector(descriptions);
    const auto result = orbpro::propagator::CreatePropagatorDescribeSourcesBatchResult(
        builder,
        catalogHandle,
        descriptionsOffset
    );
    builder.Finish(result);
    return copyFlatBufferToHeap(builder, payloadSizeOut);
}

struct TemeStateVector {
    double position[3];
    double velocity[3];
};

double chebyshevEval(const std::array<double, CHEBY_NPTS>& coeffs, double tau) {
    double bKPlusOne = 0.0;
    double bKPlusTwo = 0.0;
    for (int i = CHEBY_N; i >= 1; i--) {
        const double bK = 2.0 * tau * bKPlusOne - bKPlusTwo + coeffs[static_cast<size_t>(i)];
        bKPlusTwo = bKPlusOne;
        bKPlusOne = bK;
    }
    return tau * bKPlusOne - bKPlusTwo + coeffs[0];
}

void chebyshevFit(
    const std::array<double, CHEBY_NPTS>& values,
    std::array<double, CHEBY_NPTS>& coeffs
) {
    for (int j = 0; j <= CHEBY_N; j++) {
        double sum = 0.0;
        for (int k = 0; k <= CHEBY_N; k++) {
            const double weight = (k == 0 || k == CHEBY_N) ? 0.5 : 1.0;
            sum += weight * values[static_cast<size_t>(k)] * std::cos(j * k * PI / CHEBY_N);
        }
        coeffs[static_cast<size_t>(j)] = 2.0 * sum / CHEBY_N;
    }
    coeffs[0] *= 0.5;
    coeffs[CHEBY_N] *= 0.5;
}

void evaluateSegmentState(
    const CertifiedTrajectorySegment& segment,
    double julianDate,
    TemeStateVector& stateOut
) {
    const double halfSpan = (segment.endJd - segment.startJd) * 0.5;
    const double mid = (segment.startJd + segment.endJd) * 0.5;
    const double tau = halfSpan > 0.0 ? (julianDate - mid) / halfSpan : 0.0;
    for (size_t axis = 0; axis < 3; axis++) {
        stateOut.position[axis] = chebyshevEval(segment.coeffs[axis], tau);
        stateOut.velocity[axis] = chebyshevEval(segment.coeffs[axis + 3], tau);
    }
}

std::array<double, CHEBY_NPTS> buildChebyshevNodeJulianDates(
    double startJd,
    double endJd
) {
    std::array<double, CHEBY_NPTS> sampleJds = {};
    const double mid = (startJd + endJd) * 0.5;
    const double halfSpan = (endJd - startJd) * 0.5;
    for (int k = 0; k <= CHEBY_N; k++) {
        sampleJds[static_cast<size_t>(k)] =
            mid + std::cos(k * PI / CHEBY_N) * halfSpan;
    }
    return sampleJds;
}

bool propagateStateForSegmentSource(
    uint32_t handle,
    double julianDate,
    TemeStateVector& stateOut
) {
    if (handle >= g_satellites.size()) {
        return false;
    }

    SatelliteEntity entity = g_satellites[handle];
    return propagateEntityTEME(entity, julianDate, stateOut.position, stateOut.velocity);
}

bool buildCertifiedSegment(
    uint32_t handle,
    double segmentStartJd,
    double segmentEndJd,
    CertifiedTrajectorySegment& segmentOut
) {
    if (!(segmentEndJd > segmentStartJd)) {
        return false;
    }

    std::array<double, CHEBY_NPTS> componentSamples[6];
    const double mid = (segmentStartJd + segmentEndJd) * 0.5;
    const double halfSpan = (segmentEndJd - segmentStartJd) * 0.5;

    for (int k = 0; k <= CHEBY_N; k++) {
        const double nodeJd = mid + std::cos(k * PI / CHEBY_N) * halfSpan;
        TemeStateVector exactState = {};
        if (!propagateStateForSegmentSource(handle, nodeJd, exactState)) {
            return false;
        }
        componentSamples[0][static_cast<size_t>(k)] = exactState.position[0];
        componentSamples[1][static_cast<size_t>(k)] = exactState.position[1];
        componentSamples[2][static_cast<size_t>(k)] = exactState.position[2];
        componentSamples[3][static_cast<size_t>(k)] = exactState.velocity[0];
        componentSamples[4][static_cast<size_t>(k)] = exactState.velocity[1];
        componentSamples[5][static_cast<size_t>(k)] = exactState.velocity[2];
    }

    segmentOut.sourceHandle = handle;
    segmentOut.startJd = segmentStartJd;
    segmentOut.endJd = segmentEndJd;
    for (size_t component = 0; component < 6; component++) {
        chebyshevFit(componentSamples[component], segmentOut.coeffs[component]);
    }

    double maxPositionErrorKm = 0.0;
    double maxVelocityErrorKmS = 0.0;
    for (int sampleIndex = 1; sampleIndex <= 5; sampleIndex++) {
        const double fraction = static_cast<double>(sampleIndex) / 6.0;
        const double sampleJd = segmentStartJd + (segmentEndJd - segmentStartJd) * fraction;
        TemeStateVector exactState = {};
        TemeStateVector approxState = {};
        if (!propagateStateForSegmentSource(handle, sampleJd, exactState)) {
            return false;
        }
        evaluateSegmentState(segmentOut, sampleJd, approxState);

        const double dx = exactState.position[0] - approxState.position[0];
        const double dy = exactState.position[1] - approxState.position[1];
        const double dz = exactState.position[2] - approxState.position[2];
        const double dvx = exactState.velocity[0] - approxState.velocity[0];
        const double dvy = exactState.velocity[1] - approxState.velocity[1];
        const double dvz = exactState.velocity[2] - approxState.velocity[2];

        maxPositionErrorKm = std::max(maxPositionErrorKm, std::sqrt(dx * dx + dy * dy + dz * dz));
        maxVelocityErrorKmS = std::max(maxVelocityErrorKmS, std::sqrt(dvx * dvx + dvy * dvy + dvz * dvz));
    }

    segmentOut.maxPositionErrorKm = maxPositionErrorKm;
    segmentOut.maxVelocityErrorKmS = maxVelocityErrorKmS;
    return true;
}

bool appendCertifiedSegmentsForRange(
    uint32_t handle,
    double segmentStartJd,
    double segmentEndJd,
    PreparedTrajectorySegmentSet& setOut
) {
    CertifiedTrajectorySegment segment = {};
    if (!buildCertifiedSegment(handle, segmentStartJd, segmentEndJd, segment)) {
        return false;
    }

    if (segment.maxPositionErrorKm <= CHEBY_CERT_POSITION_KM &&
        segment.maxVelocityErrorKmS <= CHEBY_CERT_VELOCITY_KMS) {
        setOut.maxPositionErrorKm = std::max(setOut.maxPositionErrorKm, segment.maxPositionErrorKm);
        setOut.maxVelocityErrorKmS = std::max(setOut.maxVelocityErrorKmS, segment.maxVelocityErrorKmS);
        setOut.segments.push_back(std::move(segment));
        return true;
    }

    if ((segmentEndJd - segmentStartJd) <= CHEBY_MIN_SEG_DAYS) {
        setOut.maxPositionErrorKm = std::max(setOut.maxPositionErrorKm, segment.maxPositionErrorKm);
        setOut.maxVelocityErrorKmS = std::max(setOut.maxVelocityErrorKmS, segment.maxVelocityErrorKmS);
        setOut.segments.push_back(std::move(segment));
        return true;
    }

    const double midJd = (segmentStartJd + segmentEndJd) * 0.5;
    return appendCertifiedSegmentsForRange(handle, segmentStartJd, midJd, setOut) &&
           appendCertifiedSegmentsForRange(handle, midJd, segmentEndJd, setOut);
}

bool buildPreparedTrajectorySegmentSet(
    uint32_t catalogHandle,
    const std::vector<uint32_t>& handles,
    double startJd,
    double durationDays,
    const std::string& profile,
    PreparedTrajectorySegmentSet& setOut
) {
    setOut.handle = g_nextSegmentSetHandle++;
    setOut.catalogHandle = catalogHandle;
    setOut.startJd = startJd;
    setOut.durationDays = durationDays;
    setOut.profile = profile;
    setOut.sourceHandles.clear();
    setOut.sourceHandles.reserve(handles.size());
    setOut.segments.clear();
    setOut.maxPositionErrorKm = 0.0;
    setOut.maxVelocityErrorKmS = 0.0;
    setOut.coverageComplete = false;

    const double endJd = startJd + std::max(0.0, durationDays);
    for (uint32_t handle : handles) {
        PreparedTrajectorySegmentSet sourceSet = {};
        sourceSet.maxPositionErrorKm = 0.0;
        sourceSet.maxVelocityErrorKmS = 0.0;
        bool certified = true;
        for (double segmentStartJd = startJd; segmentStartJd < endJd; segmentStartJd += CHEBY_SEG_DAYS) {
            const double segmentEndJd = std::min(segmentStartJd + CHEBY_SEG_DAYS, endJd);
            if (!appendCertifiedSegmentsForRange(handle, segmentStartJd, segmentEndJd, sourceSet)) {
                certified = false;
                break;
            }
        }
        if (!certified || sourceSet.segments.empty()) {
            continue;
        }

        setOut.sourceHandles.push_back(handle);
        setOut.maxPositionErrorKm =
            std::max(setOut.maxPositionErrorKm, sourceSet.maxPositionErrorKm);
        setOut.maxVelocityErrorKmS =
            std::max(setOut.maxVelocityErrorKmS, sourceSet.maxVelocityErrorKmS);
        setOut.segments.insert(
            setOut.segments.end(),
            std::make_move_iterator(sourceSet.segments.begin()),
            std::make_move_iterator(sourceSet.segments.end()));
    }

    setOut.coverageComplete = !setOut.sourceHandles.empty();
    return setOut.coverageComplete;
}

uint8_t* encodePrepareTrajectorySegmentsResultPayload(
    const PreparedTrajectorySegmentSet& segmentSet,
    uint32_t* payloadSizeOut
) {
    flatbuffers::FlatBufferBuilder builder(256);
    const auto result = orbpro::propagator::CreatePropagatorPrepareTrajectorySegmentsResult(
        builder,
        segmentSet.handle,
        segmentSet.coverageComplete,
        segmentSet.maxPositionErrorKm,
        segmentSet.maxVelocityErrorKmS
    );
    builder.Finish(result);
    return copyFlatBufferToHeap(builder, payloadSizeOut);
}

uint8_t* encodeSampleTrajectoryStatesResultPayload(
    uint32_t catalogHandle,
    const std::vector<uint32_t>& handles,
    double startJd,
    double durationDays,
    uint32_t* payloadSizeOut
) {
    flatbuffers::FlatBufferBuilder builder(4096);
    const double endJd = startJd + std::max(0.0, durationDays);
    std::vector<double> sampleJds;
    if (endJd <= startJd) {
        const auto segmentSampleJds = buildChebyshevNodeJulianDates(startJd, startJd);
        sampleJds.assign(segmentSampleJds.begin(), segmentSampleJds.end());
    } else {
        const size_t segmentCount = static_cast<size_t>(
            std::ceil(std::max(0.0, endJd - startJd) / CHEBY_SEG_DAYS));
        sampleJds.reserve(std::max<size_t>(1, segmentCount) * CHEBY_NPTS);
        for (double segmentStartJd = startJd;
             segmentStartJd < endJd;
             segmentStartJd += CHEBY_SEG_DAYS) {
            const double segmentEndJd = std::min(segmentStartJd + CHEBY_SEG_DAYS, endJd);
            const auto segmentSampleJds =
                buildChebyshevNodeJulianDates(segmentStartJd, segmentEndJd);
            sampleJds.insert(
                sampleJds.end(),
                segmentSampleJds.begin(),
                segmentSampleJds.end());
        }
    }

    std::vector<uint32_t> sampledHandles;
    sampledHandles.reserve(handles.size());
    std::vector<orbpro::propagator::StateVector> sampledStates;
    sampledStates.reserve(handles.size() * sampleJds.size());

    for (uint32_t handle : handles) {
        std::vector<TemeStateVector> sourceStates(sampleJds.size());
        bool valid = true;
        for (size_t sampleIndex = 0; sampleIndex < sampleJds.size(); sampleIndex++) {
            if (!propagateStateForSegmentSource(handle, sampleJds[sampleIndex], sourceStates[sampleIndex])) {
                valid = false;
                break;
            }
        }
        if (!valid) {
            continue;
        }

        sampledHandles.push_back(handle);
        for (size_t sampleIndex = 0; sampleIndex < sampleJds.size(); sampleIndex++) {
            const auto& sample = sourceStates[sampleIndex];
            sampledStates.emplace_back(
                sampleJds[sampleIndex],
                orbpro::Vec3(sample.position[0], sample.position[1], sample.position[2]),
                orbpro::Vec3(sample.velocity[0], sample.velocity[1], sample.velocity[2]),
                orbpro::propagator::ReferenceFrame_TEME,
                orbpro::propagator::StateFlags_VALID
            );
        }
    }

    const auto sampleJdsOffset = builder.CreateVector(sampleJds.data(), sampleJds.size());
    const auto sourceHandlesOffset = builder.CreateVector(sampledHandles);
    const auto statesOffset = builder.CreateVectorOfStructs(sampledStates);
    const auto result = orbpro::propagator::CreatePropagatorSampleTrajectoryStatesResult(
        builder,
        catalogHandle,
        startJd,
        durationDays,
        sampleJdsOffset,
        sourceHandlesOffset,
        orbpro::propagator::ReferenceFrame_TEME,
        statesOffset
    );
    builder.Finish(result, "PSSR");
    return copyFlatBufferToHeap(builder, payloadSizeOut);
}

uint8_t* encodeDescribeTrajectorySegmentsResultPayload(
    const PreparedTrajectorySegmentSet& segmentSet,
    const std::vector<uint32_t>& handles,
    uint32_t* payloadSizeOut
) {
    flatbuffers::FlatBufferBuilder builder(4096);
    std::vector<flatbuffers::Offset<orbpro::propagator::PropagatorTrajectorySegment>> segments;
    segments.reserve(segmentSet.segments.size());

    for (const CertifiedTrajectorySegment& segment : segmentSet.segments) {
        if (std::find(handles.begin(), handles.end(), segment.sourceHandle) == handles.end()) {
            continue;
        }

        const auto xCoefficients = builder.CreateVector(
            segment.coeffs[0].data(),
            segment.coeffs[0].size()
        );
        const auto yCoefficients = builder.CreateVector(
            segment.coeffs[1].data(),
            segment.coeffs[1].size()
        );
        const auto zCoefficients = builder.CreateVector(
            segment.coeffs[2].data(),
            segment.coeffs[2].size()
        );
        const auto vxCoefficients = builder.CreateVector(
            segment.coeffs[3].data(),
            segment.coeffs[3].size()
        );
        const auto vyCoefficients = builder.CreateVector(
            segment.coeffs[4].data(),
            segment.coeffs[4].size()
        );
        const auto vzCoefficients = builder.CreateVector(
            segment.coeffs[5].data(),
            segment.coeffs[5].size()
        );

        segments.push_back(orbpro::propagator::CreatePropagatorTrajectorySegment(
            builder,
            segment.sourceHandle,
            segment.startJd,
            segment.endJd,
            static_cast<uint32_t>(CHEBY_N),
            orbpro::propagator::ReferenceFrame_TEME,
            xCoefficients,
            yCoefficients,
            zCoefficients,
            vxCoefficients,
            vyCoefficients,
            vzCoefficients,
            segment.maxPositionErrorKm,
            segment.maxVelocityErrorKmS
        ));
    }

    const auto segmentsOffset = builder.CreateVector(segments);
    const auto result = orbpro::propagator::CreatePropagatorDescribeTrajectorySegmentsResult(
        builder,
        segmentSet.handle,
        segmentsOffset
    );
    builder.Finish(result);
    return copyFlatBufferToHeap(builder, payloadSizeOut);
}

}  // anonymous namespace

// =============================================================================
// Plugin Exports
// =============================================================================

extern "C" {

ORBPRO_EXPORT int32_t plugin_catalog_begin_transaction();
ORBPRO_EXPORT int32_t plugin_catalog_commit_transaction();
ORBPRO_EXPORT int32_t plugin_catalog_rollback_transaction();
ORBPRO_EXPORT int32_t plugin_query_entity_indices_by_name(
    const char* query,
    uint32_t* out_indices,
    uint32_t max_count
);
ORBPRO_EXPORT int32_t plugin_query_visibility_mask_by_name(
    const char* query,
    uint8_t* mask_out,
    uint32_t mask_count
);
ORBPRO_EXPORT int32_t plugin_query_entity_rows_by_name(
    const char* query,
    OrbProCatalogRow* out_rows,
    uint32_t max_count
);
ORBPRO_EXPORT int32_t plugin_get_entity_catalog_row(
    uint32_t entity_index,
    OrbProCatalogRow* out_row
);

// -----------------------------------------------------------------------------
// plugin_init_omm — Initialize from binary OMM records in linear memory.
// Records are grouped by NORAD_CAT_ID: each unique NORAD ID creates one entity
// with all its OMMs stored in the satrecs map (keyed by epoch JD).
// entity.satrec is set to the most recent (highest epoch) OMM for each entity.
// Returns: number of entities initialized (>0), or negative error code
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_init_omm(const OrbProOMMRecord* records, uint32_t count) {
    g_satellites.clear();
    g_segmentSets.clear();
    g_nextSegmentSetHandle = 1;
    g_initialized = false;
    if (!g_pendingEntityNamesFromFlatbuffer) {
        g_pendingEntityMetadataByNorad.clear();
    }

    if (!ensureOmmDatabase()) {
        g_pendingEntityMetadataByNorad.clear();
        g_pendingEntityNamesFromFlatbuffer = false;
        return -static_cast<int32_t>(ORBPRO_ERROR_INIT_FAILED);
    }
    clearOmmDatabase();

    if (records == nullptr || count == 0) {
        g_initialized = true;
        g_pendingEntityMetadataByNorad.clear();
        g_pendingEntityNamesFromFlatbuffer = false;
        return 0;
    }

    // Group records by NORAD_CAT_ID, preserving insertion order
    // Use a map from noradId -> entity index for lookup
    std::map<uint32_t, uint32_t> noradToIndex;

    for (uint32_t i = 0; i < count; i++) {
        const OrbProOMMRecord& omm = records[i];

        elsetrec satrec;
        if (!initSatrecFromOMM(omm, satrec)) {
            continue;
        }

        double epochJD = satrec.jdsatepoch + satrec.jdsatepochF;
        int64_t ommPointer = insertOmmRecord(omm, epochJD);
        if (ommPointer <= 0) {
            continue;
        }

        // Test propagation at epoch (t=0)
        double r[3], v[3];
        bool propOk = SGP4Funcs::sgp4(satrec, 0.0, r, v);
        if (!propOk || satrec.error != 0) {
            continue;
        }

        auto it = noradToIndex.find(omm.norad_cat_id);
        if (it == noradToIndex.end()) {
            // First OMM for this NORAD ID — create new entity.
            // Leave satrecs map EMPTY for the common single-OMM case.
            // The fast path in propagateEntityTEME (satrecs.size() <= 1)
            // uses entity.satrec directly, avoiding map overhead entirely.
            SatelliteEntity entity;
            entity.noradId = omm.norad_cat_id;
            entity.valid = true;
            entity.lastEpochJD = epochJD;
            entity.r[0] = r[0]; entity.r[1] = r[1]; entity.r[2] = r[2];
            entity.v[0] = v[0]; entity.v[1] = v[1]; entity.v[2] = v[2];
            entity.currentSatrecEpoch = epochJD;
            entity.currentOmmPointer = ommPointer;
            entity.mode = MODE_NEAREST_EPOCH;
            entity.satrec = satrec;

            uint32_t idx = static_cast<uint32_t>(g_satellites.size());
            noradToIndex[omm.norad_cat_id] = idx;
            g_satellites.push_back(entity);
        } else {
            // Additional OMM for existing entity — add to satrecs map.
            // If the map is empty (single-OMM fast path), seed it with the
            // existing satrec before adding the new one.
            SatelliteEntity& entity = g_satellites[it->second];
            if (entity.satrecs.empty()) {
                entity.satrecs[entity.currentSatrecEpoch] = entity.satrec;
                if (entity.currentOmmPointer > 0) {
                    entity.ommPointers[entity.currentSatrecEpoch] = entity.currentOmmPointer;
                }
            }
            entity.satrecs[epochJD] = satrec;
            entity.ommPointers[epochJD] = ommPointer;

            // Keep entity.satrec pointing to the most recent epoch
            if (epochJD > entity.currentSatrecEpoch) {
                entity.satrec = satrec;
                entity.currentSatrecEpoch = epochJD;
                entity.currentOmmPointer = ommPointer;
                entity.lastEpochJD = epochJD;
                entity.r[0] = r[0]; entity.r[1] = r[1]; entity.r[2] = r[2];
                entity.v[0] = v[0]; entity.v[1] = v[1]; entity.v[2] = v[2];
            }
        }
    }

    g_initialized = true;

    if (g_satellites.empty()) {
        g_pendingEntityMetadataByNorad.clear();
        g_pendingEntityNamesFromFlatbuffer = false;
        return -static_cast<int32_t>(ORBPRO_ERROR_INIT_FAILED);
    }

    bool catalogOk = true;
    for (uint32_t i = 0; i < static_cast<uint32_t>(g_satellites.size()); i++) {
        const auto& entity = g_satellites[i];
        const PendingCatalogMetadata metadata = resolvePendingCatalogMetadata(
            entity.noradId,
            i
        );
        if (
            !upsertEntityCatalog(
                i,
                entity.noradId,
                metadata.objectName,
                metadata.objectId
            )
        ) {
            catalogOk = false;
            break;
        }
    }
    if (!catalogOk) {
        g_satellites.clear();
        g_initialized = false;
        clearOmmDatabase();
        g_pendingEntityMetadataByNorad.clear();
        g_pendingEntityNamesFromFlatbuffer = false;
        return -static_cast<int32_t>(ORBPRO_ERROR_INIT_FAILED);
    }

    g_pendingEntityMetadataByNorad.clear();
    g_pendingEntityNamesFromFlatbuffer = false;
    return static_cast<int32_t>(g_satellites.size());
}

// -----------------------------------------------------------------------------
// plugin_init — Generic init (delegates to plugin_init_omm)
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_init(const uint8_t* data, size_t len) {
    if (len < sizeof(OrbProOMMRecord)) {
        g_initialized = true;
        return 0;
    }
    uint32_t count = static_cast<uint32_t>(len / sizeof(OrbProOMMRecord));
    return plugin_init_omm(reinterpret_cast<const OrbProOMMRecord*>(data), count);
}

// -----------------------------------------------------------------------------
// plugin_init_omm_flatbuffer_stream — Initialize from OMM FlatBuffer bytes.
// Accepts either:
//   1) a single `$OMM` payload
//   2) a size-prefixed stream of `$OMM` payloads
// Returns: number of entities initialized (>0), or negative error code
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_init_omm_flatbuffer_stream(const uint8_t* data, size_t len) {
    std::vector<OrbProOMMRecord> records;
    if (!parseOmmFlatBufferStream(data, len, records)) {
        g_satellites.clear();
        g_initialized = false;
        g_pendingEntityMetadataByNorad.clear();
        g_pendingEntityNamesFromFlatbuffer = false;
        return -static_cast<int32_t>(ORBPRO_ERROR_INIT_FAILED);
    }

    if (records.empty()) {
        return plugin_init_omm(nullptr, 0);
    }

    return plugin_init_omm(records.data(), static_cast<uint32_t>(records.size()));
}

// -----------------------------------------------------------------------------
// plugin_stream_invoke — Canonical stream-based method entrypoint.
// Supported methods:
//   - ingest_omm: accepts direct `$OMM`, size-prefixed `$OMM` streams, or
//                 `$REC` payloads containing OMM records
//   - upsert_cat: accepts direct `$CAT`, size-prefixed `$CAT` streams, or
//                 `$REC` payloads containing CAT records
//   - propagate_state: accepts PropagatorBatchRequest frames and emits
//                      PropagatorState frames
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
uint8_t* plugin_stream_invoke(
    const uint8_t* request_data,
    size_t request_size,
    uint32_t* response_size_out
) {
    if (response_size_out != nullptr) {
        *response_size_out = 0;
    }
    if (request_data == nullptr || request_size == 0) {
        return buildErrorStreamInvokeResponse(-1, "Missing stream invoke request bytes.", response_size_out);
    }

    flatbuffers::Verifier verifier(request_data, request_size);
    if (!verifier.VerifyBuffer<orbpro::plugin::StreamInvokeRequest>(nullptr)) {
        return buildErrorStreamInvokeResponse(-1, "Invalid StreamInvokeRequest payload.", response_size_out);
    }

    const auto* request = flatbuffers::GetRoot<orbpro::plugin::StreamInvokeRequest>(request_data);
    if (request == nullptr || request->method_id() == nullptr) {
        return buildErrorStreamInvokeResponse(-1, "StreamInvokeRequest is missing method_id.", response_size_out);
    }

    const std::string methodId = request->method_id()->str();
    const auto* inputs = request->inputs();

    if (methodId == "ingest_omm") {
        std::vector<OrbProOMMRecord> records;
        std::map<uint32_t, PendingCatalogMetadata> metadataByNorad;
        if (inputs != nullptr) {
            for (flatbuffers::uoffset_t inputIndex = 0; inputIndex < inputs->size(); inputIndex++) {
                const auto* input = inputs->Get(inputIndex);
                if (input == nullptr || input->size() == 0 || input->offset() == 0) {
                    continue;
                }

                const auto* payload = reinterpret_cast<const uint8_t*>(
                    static_cast<uintptr_t>(input->offset())
                );
                std::vector<OrbProOMMRecord> inputRecords;
                if (!parseOmmFlatBufferStream(payload, input->size(), inputRecords)) {
                    g_pendingEntityMetadataByNorad.clear();
                    g_pendingEntityNamesFromFlatbuffer = false;
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "ingest_omm expects direct $OMM, size-prefixed OMM stream, or $REC payloads.",
                        response_size_out
                    );
                }

                records.insert(records.end(), inputRecords.begin(), inputRecords.end());
                for (const auto& entry : g_pendingEntityMetadataByNorad) {
                    metadataByNorad[entry.first] = entry.second;
                }
                g_pendingEntityMetadataByNorad.clear();
                g_pendingEntityNamesFromFlatbuffer = false;
            }
        }

        if (!ingestOmmRecords(records, metadataByNorad)) {
            return buildErrorStreamInvokeResponse(-1, "ingest_omm failed to apply any OMM records.", response_size_out);
        }
        return buildEmptyStreamInvokeResponse(response_size_out);
    }

    if (methodId == "upsert_cat") {
        std::vector<PendingCatCatalogRecord> records;
        if (inputs != nullptr) {
            for (flatbuffers::uoffset_t inputIndex = 0; inputIndex < inputs->size(); inputIndex++) {
                const auto* input = inputs->Get(inputIndex);
                if (input == nullptr || input->size() == 0 || input->offset() == 0) {
                    continue;
                }

                const auto* payload = reinterpret_cast<const uint8_t*>(
                    static_cast<uintptr_t>(input->offset())
                );
                std::vector<PendingCatCatalogRecord> inputRecords;
                if (!parseCatFlatBufferStream(payload, input->size(), inputRecords)) {
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "upsert_cat expects direct $CAT, size-prefixed CAT stream, or $REC payloads.",
                        response_size_out
                    );
                }
                records.insert(records.end(), inputRecords.begin(), inputRecords.end());
            }
        }

        if (records.empty()) {
            return buildEmptyStreamInvokeResponse(response_size_out);
        }

        bool inTransaction = false;
        if (plugin_catalog_begin_transaction() == 0) {
            inTransaction = true;
        }

        for (const PendingCatCatalogRecord& record : records) {
            if (record.noradCatId == 0) {
                continue;
            }
            if (!upsertCatCatalog(
                    record.noradCatId,
                    record.objectName,
                    record.objectId,
                    std::string())) {
                if (inTransaction) {
                    plugin_catalog_rollback_transaction();
                }
                return buildErrorStreamInvokeResponse(
                    -1,
                    "upsert_cat failed to apply one or more CAT records.",
                    response_size_out
                );
            }
        }

        if (inTransaction && plugin_catalog_commit_transaction() != 0) {
            plugin_catalog_rollback_transaction();
            return buildErrorStreamInvokeResponse(
                -1,
                "upsert_cat failed to commit CAT transaction.",
                response_size_out
            );
        }

        return buildEmptyStreamInvokeResponse(response_size_out);
    }

    if (methodId == "describe_sources_batch") {
        struct OutputFrame {
            uint8_t* payload;
            uint32_t payloadSize;
            uint64_t traceId;
            uint32_t streamId;
            uint64_t sequence;
        };

        std::vector<OutputFrame> outputFrames;
        const uint32_t outputCap = request->output_stream_cap();
        const auto releaseOutputFrames = [&outputFrames]() {
            for (const auto& output : outputFrames) {
                orbpro_free(output.payload);
            }
        };

        if (inputs != nullptr) {
            for (flatbuffers::uoffset_t inputIndex = 0; inputIndex < inputs->size(); inputIndex++) {
                const auto* input = inputs->Get(inputIndex);
                if (input == nullptr || input->size() == 0 || input->offset() == 0) {
                    continue;
                }

                const auto* payload = reinterpret_cast<const uint8_t*>(
                    static_cast<uintptr_t>(input->offset())
                );
                const orbpro::propagator::PropagatorDescribeSourcesBatchRequest* describeRequest = nullptr;
                if (!decodeDescribeSourcesBatchRequest(payload, input->size(), describeRequest) ||
                    describeRequest == nullptr) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "describe_sources_batch expects PropagatorDescribeSourcesBatchRequest input frames.",
                        response_size_out
                    );
                }

                if (outputCap > 0 &&
                    outputFrames.size() + 1 > static_cast<size_t>(outputCap)) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "describe_sources_batch output_stream_cap is smaller than the requested result count.",
                        response_size_out
                    );
                }

                std::vector<uint32_t> handles;
                if (!resolveRequestedSatelliteHandles(describeRequest->sourceHandles(), handles)) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "describe_sources_batch received one or more invalid source handles.",
                        response_size_out
                    );
                }

                uint32_t payloadSize = 0;
                uint8_t* resultPayload = encodeDescribeSourcesBatchResultPayload(
                    describeRequest->catalogHandle(),
                    handles,
                    &payloadSize
                );
                if (resultPayload == nullptr) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "Failed to encode describe_sources_batch output.",
                        response_size_out
                    );
                }

                outputFrames.push_back({
                    resultPayload,
                    payloadSize,
                    input->trace_id(),
                    input->stream_id(),
                    input->sequence(),
                });
            }
        }

        flatbuffers::FlatBufferBuilder builder(1024);
        std::vector<flatbuffers::Offset<orbpro::stream::TypedArenaBuffer>> outputs;
        outputs.reserve(outputFrames.size());
        const auto typeRef = orbpro::stream::CreateFlatBufferTypeRefDirect(
            builder,
            "orbpro.propagator.PropagatorDescribeSourcesBatchResult",
            "PSRS",
            nullptr,
            false,
            orbpro::stream::PayloadWireFormat_AlignedBinary,
            "PropagatorDescribeSourcesBatchResult",
            0,
            0,
            8
        );
        for (const auto& output : outputFrames) {
            outputs.push_back(orbpro::stream::CreateTypedArenaBufferDirect(
                builder,
                typeRef,
                "result",
                8,
                static_cast<uint32_t>(reinterpret_cast<uintptr_t>(output.payload)),
                output.payloadSize,
                orbpro::stream::BufferOwnership_BORROWED,
                0,
                orbpro::stream::BufferMutability_IMMUTABLE,
                output.traceId,
                output.streamId,
                output.sequence,
                false
            ));
        }

        const auto outputsVector = builder.CreateVector(outputs);
        const auto response = orbpro::plugin::CreateStreamInvokeResponse(
            builder,
            outputsVector,
            0,
            false,
            0,
            0
        );
        builder.Finish(response);
        uint8_t* responseBytes = copyFlatBufferToHeap(builder, response_size_out);
        if (responseBytes == nullptr) {
            releaseOutputFrames();
            return nullptr;
        }
        return responseBytes;
    }

    if (methodId == "sample_trajectory_states") {
        struct OutputFrame {
            uint8_t* payload;
            uint32_t payloadSize;
            uint64_t traceId;
            uint32_t streamId;
            uint64_t sequence;
        };

        std::vector<OutputFrame> outputFrames;
        const uint32_t outputCap = request->output_stream_cap();
        const auto releaseOutputFrames = [&outputFrames]() {
            for (const auto& output : outputFrames) {
                orbpro_free(output.payload);
            }
        };

        if (inputs != nullptr) {
            for (flatbuffers::uoffset_t inputIndex = 0; inputIndex < inputs->size(); inputIndex++) {
                const auto* input = inputs->Get(inputIndex);
                if (input == nullptr || input->size() == 0 || input->offset() == 0) {
                    continue;
                }

                const auto* payload = reinterpret_cast<const uint8_t*>(
                    static_cast<uintptr_t>(input->offset())
                );
                const orbpro::propagator::PropagatorSampleTrajectoryStatesRequest* sampleRequest = nullptr;
                if (!decodeSampleTrajectoryStatesRequest(payload, input->size(), sampleRequest) ||
                    sampleRequest == nullptr) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "sample_trajectory_states expects PropagatorSampleTrajectoryStatesRequest input frames.",
                        response_size_out
                    );
                }

                if (outputCap > 0 &&
                    outputFrames.size() + 1 > static_cast<size_t>(outputCap)) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "sample_trajectory_states output_stream_cap is smaller than the requested result count.",
                        response_size_out
                    );
                }

                std::vector<uint32_t> handles;
                if (!resolveRequestedSatelliteHandles(sampleRequest->sourceHandles(), handles)) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "sample_trajectory_states received one or more invalid source handles.",
                        response_size_out
                    );
                }

                uint32_t payloadSize = 0;
                uint8_t* resultPayload = encodeSampleTrajectoryStatesResultPayload(
                    sampleRequest->catalogHandle(),
                    handles,
                    sampleRequest->startJd(),
                    sampleRequest->durationDays(),
                    &payloadSize
                );
                if (resultPayload == nullptr) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "Failed to encode sample_trajectory_states output.",
                        response_size_out
                    );
                }

                outputFrames.push_back({
                    resultPayload,
                    payloadSize,
                    input->trace_id(),
                    input->stream_id(),
                    input->sequence(),
                });
            }
        }

        flatbuffers::FlatBufferBuilder builder(512);
        std::vector<flatbuffers::Offset<orbpro::stream::TypedArenaBuffer>> outputs;
        outputs.reserve(outputFrames.size());
        const auto typeRef = orbpro::stream::CreateFlatBufferTypeRefDirect(
            builder,
            "orbpro.propagator.PropagatorSampleTrajectoryStatesResult",
            "PSSR",
            nullptr,
            false,
            orbpro::stream::PayloadWireFormat_AlignedBinary,
            "PropagatorSampleTrajectoryStatesResult",
            0,
            0,
            8
        );
        for (const auto& output : outputFrames) {
            outputs.push_back(orbpro::stream::CreateTypedArenaBufferDirect(
                builder,
                typeRef,
                "result",
                8,
                static_cast<uint32_t>(reinterpret_cast<uintptr_t>(output.payload)),
                output.payloadSize,
                orbpro::stream::BufferOwnership_BORROWED,
                0,
                orbpro::stream::BufferMutability_IMMUTABLE,
                output.traceId,
                output.streamId,
                output.sequence,
                false
            ));
        }

        const auto outputsVector = builder.CreateVector(outputs);
        const auto response = orbpro::plugin::CreateStreamInvokeResponse(
            builder,
            outputsVector,
            0,
            false,
            0,
            0
        );
        builder.Finish(response);
        uint8_t* responseBytes = copyFlatBufferToHeap(builder, response_size_out);
        if (responseBytes == nullptr) {
            releaseOutputFrames();
            return nullptr;
        }
        return responseBytes;
    }

    if (methodId == "prepare_trajectory_segments") {
        struct OutputFrame {
            uint8_t* payload;
            uint32_t payloadSize;
            uint64_t traceId;
            uint32_t streamId;
            uint64_t sequence;
        };

        std::vector<OutputFrame> outputFrames;
        const uint32_t outputCap = request->output_stream_cap();
        const auto releaseOutputFrames = [&outputFrames]() {
            for (const auto& output : outputFrames) {
                orbpro_free(output.payload);
            }
        };

        if (inputs != nullptr) {
            for (flatbuffers::uoffset_t inputIndex = 0; inputIndex < inputs->size(); inputIndex++) {
                const auto* input = inputs->Get(inputIndex);
                if (input == nullptr || input->size() == 0 || input->offset() == 0) {
                    continue;
                }

                const auto* payload = reinterpret_cast<const uint8_t*>(
                    static_cast<uintptr_t>(input->offset())
                );
                const orbpro::propagator::PropagatorPrepareTrajectorySegmentsRequest* segmentRequest = nullptr;
                if (!decodePrepareTrajectorySegmentsRequest(payload, input->size(), segmentRequest) ||
                    segmentRequest == nullptr) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "prepare_trajectory_segments expects PropagatorPrepareTrajectorySegmentsRequest input frames.",
                        response_size_out
                    );
                }

                if (outputCap > 0 &&
                    outputFrames.size() + 1 > static_cast<size_t>(outputCap)) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "prepare_trajectory_segments output_stream_cap is smaller than the requested result count.",
                        response_size_out
                    );
                }

                std::vector<uint32_t> handles;
                if (!resolveRequestedSatelliteHandles(segmentRequest->sourceHandles(), handles)) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "prepare_trajectory_segments received one or more invalid source handles.",
                        response_size_out
                    );
                }

                PreparedTrajectorySegmentSet segmentSet = {};
                if (!buildPreparedTrajectorySegmentSet(
                        segmentRequest->catalogHandle(),
                        handles,
                        segmentRequest->startJd(),
                        segmentRequest->durationDays(),
                        segmentRequest->profile() != nullptr
                            ? segmentRequest->profile()->str()
                            : std::string("conjunction-screening"),
                        segmentSet)) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "prepare_trajectory_segments failed to certify full-window coverage.",
                        response_size_out
                    );
                }

                const uint32_t segmentSetHandle = segmentSet.handle;
                g_segmentSets[segmentSetHandle] = segmentSet;

                uint32_t payloadSize = 0;
                uint8_t* resultPayload = encodePrepareTrajectorySegmentsResultPayload(
                    g_segmentSets[segmentSetHandle],
                    &payloadSize
                );
                if (resultPayload == nullptr) {
                    g_segmentSets.erase(segmentSetHandle);
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "Failed to encode prepare_trajectory_segments output.",
                        response_size_out
                    );
                }

                outputFrames.push_back({
                    resultPayload,
                    payloadSize,
                    input->trace_id(),
                    input->stream_id(),
                    input->sequence(),
                });
            }
        }

        flatbuffers::FlatBufferBuilder builder(512);
        std::vector<flatbuffers::Offset<orbpro::stream::TypedArenaBuffer>> outputs;
        outputs.reserve(outputFrames.size());
        const auto typeRef = orbpro::stream::CreateFlatBufferTypeRefDirect(
            builder,
            "orbpro.propagator.PropagatorPrepareTrajectorySegmentsResult",
            "PTSS",
            nullptr,
            false,
            orbpro::stream::PayloadWireFormat_AlignedBinary,
            "PropagatorPrepareTrajectorySegmentsResult",
            0,
            0,
            8
        );
        for (const auto& output : outputFrames) {
            outputs.push_back(orbpro::stream::CreateTypedArenaBufferDirect(
                builder,
                typeRef,
                "result",
                8,
                static_cast<uint32_t>(reinterpret_cast<uintptr_t>(output.payload)),
                output.payloadSize,
                orbpro::stream::BufferOwnership_BORROWED,
                0,
                orbpro::stream::BufferMutability_IMMUTABLE,
                output.traceId,
                output.streamId,
                output.sequence,
                false
            ));
        }

        const auto outputsVector = builder.CreateVector(outputs);
        const auto response = orbpro::plugin::CreateStreamInvokeResponse(
            builder,
            outputsVector,
            0,
            false,
            0,
            0
        );
        builder.Finish(response);
        uint8_t* responseBytes = copyFlatBufferToHeap(builder, response_size_out);
        if (responseBytes == nullptr) {
            releaseOutputFrames();
            return nullptr;
        }
        return responseBytes;
    }

    if (methodId == "propagate_state") {
        struct OutputFrame {
            uint8_t* payload;
            uint32_t payloadSize;
            uint64_t traceId;
            uint32_t streamId;
            uint64_t sequence;
        };

        std::vector<OutputFrame> outputFrames;
        const uint32_t outputCap = request->output_stream_cap();

        if (inputs != nullptr) {
            for (flatbuffers::uoffset_t inputIndex = 0; inputIndex < inputs->size(); inputIndex++) {
                const auto* input = inputs->Get(inputIndex);
                if (input == nullptr || input->size() == 0 || input->offset() == 0) {
                    continue;
                }

                const auto* payload = reinterpret_cast<const uint8_t*>(
                    static_cast<uintptr_t>(input->offset())
                );
                const orbpro::propagator::PropagatorBatchRequest* batchRequest = nullptr;
                if (!decodePropagatorBatchRequest(payload, input->size(), batchRequest) || batchRequest == nullptr) {
                    for (const auto& output : outputFrames) {
                        orbpro_free(output.payload);
                    }
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "propagate_state expects PropagatorBatchRequest input frames.",
                        response_size_out
                    );
                }

                std::vector<uint32_t> handles;
                if (batchRequest->entity_handles() != nullptr && batchRequest->entity_handles()->size() > 0) {
                    handles.reserve(batchRequest->entity_handles()->size());
                    for (uint32_t handle : *batchRequest->entity_handles()) {
                        handles.push_back(handle);
                    }
                } else {
                    const uint32_t entityCount = static_cast<uint32_t>(g_satellites.size());
                    const uint32_t maxCount =
                        batchRequest->max_count() > 0
                            ? std::min(batchRequest->max_count(), entityCount)
                            : entityCount;
                    handles.reserve(maxCount);
                    for (uint32_t handle = 0; handle < maxCount; handle++) {
                        handles.push_back(handle);
                    }
                }

                if (outputCap > 0 &&
                    outputFrames.size() + handles.size() > static_cast<size_t>(outputCap)) {
                    for (const auto& output : outputFrames) {
                        orbpro_free(output.payload);
                    }
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "propagate_state output_stream_cap is smaller than the requested state count.",
                        response_size_out
                    );
                }

                for (size_t handleIndex = 0; handleIndex < handles.size(); handleIndex++) {
                    const uint32_t entityIndex = handles[handleIndex];
                    OrbProStateVector state;
                    orbpro_state_init(&state);
                    state.epoch = batchRequest->epoch();
                    state.reference_frame = ORBPRO_FRAME_ECEF;

                    const int32_t propagateResult = plugin_propagate(batchRequest->epoch(), entityIndex, &state);
                    const bool valid = propagateResult == 0 &&
                        entityIndex < static_cast<uint32_t>(g_satellites.size()) &&
                        (state.flags & ORBPRO_STATE_VALID) != 0;
                    const uint32_t catalogNumber =
                        entityIndex < static_cast<uint32_t>(g_satellites.size())
                            ? g_satellites[entityIndex].noradId
                            : 0;

                    uint32_t payloadSize = 0;
                    uint8_t* statePayload = encodePropagatorStatePayload(
                        state,
                        entityIndex,
                        catalogNumber,
                        valid,
                        &payloadSize
                    );
                    if (statePayload == nullptr) {
                        for (const auto& output : outputFrames) {
                            orbpro_free(output.payload);
                        }
                        return buildErrorStreamInvokeResponse(-1, "Failed to encode PropagatorState output.", response_size_out);
                    }

                    outputFrames.push_back({
                        statePayload,
                        payloadSize,
                        input->trace_id(),
                        input->stream_id(),
                        input->sequence() + static_cast<uint64_t>(handleIndex),
                    });
                }
            }
        }

        flatbuffers::FlatBufferBuilder builder(1024);
        std::vector<flatbuffers::Offset<orbpro::stream::TypedArenaBuffer>> outputs;
        outputs.reserve(outputFrames.size());
        const auto typeRef = orbpro::stream::CreateFlatBufferTypeRefDirect(
            builder,
            "orbpro.plugins.PropagatorState",
            orbpro::plugins::PropagatorStateIdentifier(),
            nullptr,
            false,
            orbpro::stream::PayloadWireFormat_AlignedBinary,
            "PropagatorState",
            0,
            0,
            8
        );
        for (const auto& output : outputFrames) {
            outputs.push_back(orbpro::stream::CreateTypedArenaBufferDirect(
                builder,
                typeRef,
                "state",
                8,
                static_cast<uint32_t>(reinterpret_cast<uintptr_t>(output.payload)),
                output.payloadSize,
                orbpro::stream::BufferOwnership_BORROWED,
                0,
                orbpro::stream::BufferMutability_IMMUTABLE,
                output.traceId,
                output.streamId,
                output.sequence,
                false
            ));
        }

        const auto outputsVector = builder.CreateVector(outputs);
        const auto response = orbpro::plugin::CreateStreamInvokeResponse(
            builder,
            outputsVector,
            0,
            false,
            0,
            0
        );
        builder.Finish(response);
        uint8_t* responseBytes = copyFlatBufferToHeap(builder, response_size_out);
        if (responseBytes == nullptr) {
            for (const auto& output : outputFrames) {
                orbpro_free(output.payload);
            }
            return nullptr;
        }
        return responseBytes;
    }

    if (methodId == "describe_trajectory_segments") {
        struct OutputFrame {
            uint8_t* payload;
            uint32_t payloadSize;
            uint64_t traceId;
            uint32_t streamId;
            uint64_t sequence;
        };

        std::vector<OutputFrame> outputFrames;
        const uint32_t outputCap = request->output_stream_cap();
        const auto releaseOutputFrames = [&outputFrames]() {
            for (const auto& output : outputFrames) {
                orbpro_free(output.payload);
            }
        };

        if (inputs != nullptr) {
            for (flatbuffers::uoffset_t inputIndex = 0; inputIndex < inputs->size(); inputIndex++) {
                const auto* input = inputs->Get(inputIndex);
                if (input == nullptr || input->size() == 0 || input->offset() == 0) {
                    continue;
                }

                const auto* payload = reinterpret_cast<const uint8_t*>(
                    static_cast<uintptr_t>(input->offset())
                );
                const orbpro::propagator::PropagatorDescribeTrajectorySegmentsRequest* describeRequest = nullptr;
                if (!decodeDescribeTrajectorySegmentsRequest(payload, input->size(), describeRequest) ||
                    describeRequest == nullptr) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "describe_trajectory_segments expects PropagatorDescribeTrajectorySegmentsRequest input frames.",
                        response_size_out
                    );
                }

                if (outputCap > 0 &&
                    outputFrames.size() + 1 > static_cast<size_t>(outputCap)) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "describe_trajectory_segments output_stream_cap is smaller than the requested result count.",
                        response_size_out
                    );
                }

                const auto segmentSetIt = g_segmentSets.find(describeRequest->segmentSetHandle());
                if (segmentSetIt == g_segmentSets.end()) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "describe_trajectory_segments received an unknown segment set handle.",
                        response_size_out
                    );
                }

                std::vector<uint32_t> handles;
                if (!resolveRequestedSegmentHandles(
                        describeRequest->sourceHandles(),
                        segmentSetIt->second,
                        handles)) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "describe_trajectory_segments received one or more invalid source handles.",
                        response_size_out
                    );
                }

                uint32_t payloadSize = 0;
                uint8_t* resultPayload = encodeDescribeTrajectorySegmentsResultPayload(
                    segmentSetIt->second,
                    handles,
                    &payloadSize
                );
                if (resultPayload == nullptr) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "Failed to encode describe_trajectory_segments output.",
                        response_size_out
                    );
                }

                outputFrames.push_back({
                    resultPayload,
                    payloadSize,
                    input->trace_id(),
                    input->stream_id(),
                    input->sequence(),
                });
            }
        }

        flatbuffers::FlatBufferBuilder builder(1024);
        std::vector<flatbuffers::Offset<orbpro::stream::TypedArenaBuffer>> outputs;
        outputs.reserve(outputFrames.size());
        const auto typeRef = orbpro::stream::CreateFlatBufferTypeRefDirect(
            builder,
            "orbpro.propagator.PropagatorDescribeTrajectorySegmentsResult",
            "PTDS",
            nullptr,
            false,
            orbpro::stream::PayloadWireFormat_AlignedBinary,
            "PropagatorDescribeTrajectorySegmentsResult",
            0,
            0,
            8
        );
        for (const auto& output : outputFrames) {
            outputs.push_back(orbpro::stream::CreateTypedArenaBufferDirect(
                builder,
                typeRef,
                "result",
                8,
                static_cast<uint32_t>(reinterpret_cast<uintptr_t>(output.payload)),
                output.payloadSize,
                orbpro::stream::BufferOwnership_BORROWED,
                0,
                orbpro::stream::BufferMutability_IMMUTABLE,
                output.traceId,
                output.streamId,
                output.sequence,
                false
            ));
        }

        const auto outputsVector = builder.CreateVector(outputs);
        const auto response = orbpro::plugin::CreateStreamInvokeResponse(
            builder,
            outputsVector,
            0,
            false,
            0,
            0
        );
        builder.Finish(response);
        uint8_t* responseBytes = copyFlatBufferToHeap(builder, response_size_out);
        if (responseBytes == nullptr) {
            releaseOutputFrames();
            return nullptr;
        }
        return responseBytes;
    }

    if (methodId == "catalog_query") {
        struct OutputFrame {
            uint8_t* payload;
            uint32_t payloadSize;
            uint64_t traceId;
            uint32_t streamId;
            uint64_t sequence;
        };

        std::vector<OutputFrame> outputFrames;
        const uint32_t outputCap = request->output_stream_cap();
        const auto releaseOutputFrames = [&outputFrames]() {
            for (const auto& output : outputFrames) {
                orbpro_free(output.payload);
            }
        };

        if (inputs != nullptr) {
            for (flatbuffers::uoffset_t inputIndex = 0; inputIndex < inputs->size(); inputIndex++) {
                const auto* input = inputs->Get(inputIndex);
                if (input == nullptr || input->size() == 0 || input->offset() == 0) {
                    continue;
                }
                if (outputCap > 0 &&
                    outputFrames.size() >= static_cast<size_t>(outputCap)) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "catalog_query output_stream_cap is smaller than the requested query count.",
                        response_size_out
                    );
                }

                const auto* payload = reinterpret_cast<const uint8_t*>(
                    static_cast<uintptr_t>(input->offset())
                );
                const orbpro::query::CatalogQueryRequest* queryRequest = nullptr;
                if (!decodeCatalogQueryRequest(payload, input->size(), queryRequest) ||
                    queryRequest == nullptr) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "catalog_query expects CatalogQueryRequest input frames.",
                        response_size_out
                    );
                }

                const std::string queryText =
                    queryRequest->query() != nullptr
                        ? queryRequest->query()->str()
                        : std::string();
                const uint32_t entityCount =
                    static_cast<uint32_t>(g_satellites.size());
                uint32_t payloadSize = 0;
                uint8_t* resultPayload = nullptr;

                switch (queryRequest->query_kind()) {
                    case orbpro::query::CatalogQueryKind_ROWS: {
                        const uint32_t maxCount =
                            queryRequest->max_count() > 0
                                ? queryRequest->max_count()
                                : entityCount;
                        std::vector<OrbProCatalogRow> rows(maxCount);
                        int32_t written = 0;
                        if (maxCount > 0) {
                            written = plugin_query_entity_rows_by_name(
                                queryText.c_str(),
                                rows.data(),
                                maxCount
                            );
                            if (written < 0) {
                                releaseOutputFrames();
                                return buildErrorStreamInvokeResponse(
                                    -1,
                                    "catalog_query rows lookup failed.",
                                    response_size_out
                                );
                            }
                            rows.resize(static_cast<size_t>(written));
                        }
                        resultPayload = encodeCatalogQueryResultPayload(
                            orbpro::query::CatalogQueryKind_ROWS,
                            &rows,
                            nullptr,
                            nullptr,
                            0,
                            0,
                            nullptr,
                            &payloadSize
                        );
                        break;
                    }
                    case orbpro::query::CatalogQueryKind_ENTITY_INDICES: {
                        const uint32_t maxCount =
                            queryRequest->max_count() > 0
                                ? queryRequest->max_count()
                                : entityCount;
                        std::vector<uint32_t> entityIndices(maxCount);
                        int32_t written = 0;
                        if (maxCount > 0) {
                            written = plugin_query_entity_indices_by_name(
                                queryText.c_str(),
                                entityIndices.data(),
                                maxCount
                            );
                            if (written < 0) {
                                releaseOutputFrames();
                                return buildErrorStreamInvokeResponse(
                                    -1,
                                    "catalog_query entity-index lookup failed.",
                                    response_size_out
                                );
                            }
                            entityIndices.resize(static_cast<size_t>(written));
                        }
                        resultPayload = encodeCatalogQueryResultPayload(
                            orbpro::query::CatalogQueryKind_ENTITY_INDICES,
                            nullptr,
                            &entityIndices,
                            nullptr,
                            0,
                            0,
                            nullptr,
                            &payloadSize
                        );
                        break;
                    }
                    case orbpro::query::CatalogQueryKind_VISIBILITY_MASK: {
                        const uint32_t maskCount =
                            queryRequest->entity_count() > 0
                                ? queryRequest->entity_count()
                                : entityCount;
                        std::vector<uint8_t> mask(maskCount, 0);
                        int32_t visibleCount = 0;
                        if (maskCount > 0) {
                            visibleCount = plugin_query_visibility_mask_by_name(
                                queryText.c_str(),
                                mask.data(),
                                maskCount
                            );
                            if (visibleCount < 0) {
                                releaseOutputFrames();
                                return buildErrorStreamInvokeResponse(
                                    -1,
                                    "catalog_query visibility-mask lookup failed.",
                                    response_size_out
                                );
                            }
                        }
                        resultPayload = encodeCatalogQueryResultPayload(
                            orbpro::query::CatalogQueryKind_VISIBILITY_MASK,
                            nullptr,
                            nullptr,
                            &mask,
                            static_cast<uint32_t>(visibleCount),
                            0,
                            nullptr,
                            &payloadSize
                        );
                        break;
                    }
                    case orbpro::query::CatalogQueryKind_CATALOG_ROW: {
                        OrbProCatalogRow row;
                        std::memset(&row, 0, sizeof(row));
                        const bool found = plugin_get_entity_catalog_row(
                            queryRequest->entity_index(),
                            &row
                        ) == 0;
                        resultPayload = encodeCatalogQueryResultPayload(
                            orbpro::query::CatalogQueryKind_CATALOG_ROW,
                            nullptr,
                            nullptr,
                            nullptr,
                            0,
                            queryRequest->entity_index(),
                            found ? &row : nullptr,
                            &payloadSize
                        );
                        break;
                    }
                    default:
                        releaseOutputFrames();
                        return buildErrorStreamInvokeResponse(
                            -1,
                            "catalog_query received an unsupported query kind.",
                            response_size_out
                        );
                }

                if (resultPayload == nullptr) {
                    releaseOutputFrames();
                    return buildErrorStreamInvokeResponse(
                        -1,
                        "Failed to encode CatalogQueryResult output.",
                        response_size_out
                    );
                }

                outputFrames.push_back({
                    resultPayload,
                    payloadSize,
                    input->trace_id(),
                    input->stream_id(),
                    input->sequence(),
                });
            }
        }

        flatbuffers::FlatBufferBuilder builder(1024);
        std::vector<flatbuffers::Offset<orbpro::stream::TypedArenaBuffer>> outputs;
        outputs.reserve(outputFrames.size());
        const auto typeRef = orbpro::stream::CreateFlatBufferTypeRefDirect(
            builder,
            "orbpro.query.CatalogQueryResult",
            orbpro::query::CatalogQueryResultIdentifier(),
            nullptr,
            false,
            orbpro::stream::PayloadWireFormat_AlignedBinary,
            "CatalogQueryResult",
            0,
            0,
            8
        );
        for (const auto& output : outputFrames) {
            outputs.push_back(orbpro::stream::CreateTypedArenaBufferDirect(
                builder,
                typeRef,
                "results",
                8,
                static_cast<uint32_t>(reinterpret_cast<uintptr_t>(output.payload)),
                output.payloadSize,
                orbpro::stream::BufferOwnership_BORROWED,
                0,
                orbpro::stream::BufferMutability_IMMUTABLE,
                output.traceId,
                output.streamId,
                output.sequence,
                false
            ));
        }

        const auto outputsVector = builder.CreateVector(outputs);
        const auto response = orbpro::plugin::CreateStreamInvokeResponse(
            builder,
            outputsVector,
            0,
            false,
            0,
            0
        );
        builder.Finish(response);
        uint8_t* responseBytes = copyFlatBufferToHeap(builder, response_size_out);
        if (responseBytes == nullptr) {
            releaseOutputFrames();
            return nullptr;
        }
        return responseBytes;
    }

    return buildErrorStreamInvokeResponse(-1, "Unsupported stream method for SGP4 plugin.", response_size_out);
}

// -----------------------------------------------------------------------------
// plugin_propagate — Propagate a single entity to a Julian date.
// Output is ECEF meters.
// Returns: 0 on success, negative error code on failure
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_propagate(double julian_date, uint32_t entity_index,
                          OrbProStateVector* out) {
    if (!g_initialized || entity_index >= g_satellites.size()) {
        return -1;
    }

    SatelliteEntity& entity = g_satellites[entity_index];
    if (!entity.valid) {
        orbpro_state_init(out);
        return -1;
    }

    bool ok = propagateEntityTEME(entity, julian_date, entity.r, entity.v);

    if (!ok) {
        orbpro_state_init(out);
        out->flags = ORBPRO_STATE_DECAYED;
        // Do NOT set entity.valid = false — SGP4 errors are often transient
        // (e.g. clock far from TLE epoch). The entity may propagate fine at
        // other times. Permanent invalidation was causing dots to disappear
        // and camera tracking to break.
        return -1;
    }

    entity.lastEpochJD = julian_date;
    writeStateVector(entity, julian_date, out);
    return 0;
}

// -----------------------------------------------------------------------------
// plugin_propagate_batch — Propagate all entities to a Julian date
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_propagate_batch(double julian_date, OrbProStateVector* out,
                                uint32_t count) {
    if (!g_initialized) return -1;

    uint32_t n = count;
    if (n > g_satellites.size()) n = static_cast<uint32_t>(g_satellites.size());

    for (uint32_t i = 0; i < n; i++) {
        plugin_propagate(julian_date, i, &out[i]);
    }

    return 0;
}

// -----------------------------------------------------------------------------
// plugin_propagate_path — Propagate one entity at many uniformly-spaced times.
// Output is positions-only (3 × float64 per sample, ECEF meters).
// Single WASM call replaces N individual plugin_propagate calls.
// Returns: number of samples written, or negative error code
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_propagate_path(uint32_t entity_index, double start_jd,
                               double step_days, uint32_t count,
                               double* positions_out) {
    if (!g_initialized || entity_index >= g_satellites.size()) return -1;

    SatelliteEntity& entity = g_satellites[entity_index];
    if (!entity.valid) return -1;

    uint32_t written = 0;

    for (uint32_t i = 0; i < count; i++) {
        double jd = start_jd + i * step_days;
        double r[3], v[3];
        bool ok = propagateEntityTEME(entity, jd, r, v);

        if (!ok) {
            positions_out[written * 3 + 0] = 0.0 / 0.0;
            positions_out[written * 3 + 1] = 0.0 / 0.0;
            positions_out[written * 3 + 2] = 0.0 / 0.0;
            written++;
            continue;
        }

        double gmst = SGP4Funcs::gstime_SGP4(jd);
        double cosG = cos(gmst);
        double sinG = sin(gmst);

        positions_out[written * 3 + 0] = ( cosG * r[0] + sinG * r[1]) * 1000.0;
        positions_out[written * 3 + 1] = (-sinG * r[0] + cosG * r[1]) * 1000.0;
        positions_out[written * 3 + 2] = r[2] * 1000.0;
        written++;
    }

    entity.lastEpochJD = start_jd + (count - 1) * step_days;
    return static_cast<int32_t>(written);
}

// -----------------------------------------------------------------------------
// get_entity_info — Get entity metadata (NORAD ID at offset 0)
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t get_entity_info(uint32_t entity_index, uint8_t* out) {
    if (!g_initialized || entity_index >= g_satellites.size() || !out) return -1;

    uint32_t noradId = g_satellites[entity_index].noradId;
    std::memcpy(out, &noradId, sizeof(uint32_t));
    return 0;
}

// -----------------------------------------------------------------------------
// plugin_query_entity_indices_by_name — Query entity indices by object name.
// Empty query returns all entity indices, sorted by entity_index.
// Returns: number of indices written, or -1 on error.
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_query_entity_indices_by_name(
    const char* query,
    uint32_t* out_indices,
    uint32_t max_count
) {
    if (!g_initialized || out_indices == nullptr || max_count == 0) {
        return 0;
    }
    if (!ensureOmmDatabase()) {
        return -1;
    }

    const std::string queryText = trimAsciiWhitespace(query);
    sqlite3_stmt* stmt = nullptr;

    if (queryText.empty()) {
        stmt = g_queryEntityAllStmt;
        sqlite3_reset(stmt);
        sqlite3_clear_bindings(stmt);
        const int sqliteLimit =
            max_count > 0x7fffffffU ? 0x7fffffff : static_cast<int>(max_count);
        if (sqlite3_bind_int(stmt, 1, sqliteLimit) != SQLITE_OK) {
            sqlite3_reset(stmt);
            sqlite3_clear_bindings(stmt);
            return -1;
        }
    } else {
        stmt = g_queryEntityByNameStmt;
        if (!bindNameQueryStatement(stmt, queryText, max_count)) {
            return -1;
        }
    }

    uint32_t written = 0;
    while (written < max_count) {
        const int rc = sqlite3_step(stmt);
        if (rc == SQLITE_DONE) {
            break;
        }
        if (rc != SQLITE_ROW) {
            sqlite3_reset(stmt);
            sqlite3_clear_bindings(stmt);
            return -1;
        }
        const int entityIndex = sqlite3_column_int(stmt, 0);
        if (entityIndex >= 0) {
            out_indices[written++] = static_cast<uint32_t>(entityIndex);
        }
    }

    sqlite3_reset(stmt);
    sqlite3_clear_bindings(stmt);
    return static_cast<int32_t>(written);
}

// -----------------------------------------------------------------------------
// plugin_query_visibility_mask_by_name — Query a show/hide mask by object name.
// Empty query sets all entries to visible.
// Returns: number of entries set to visible, or -1 on error.
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_query_visibility_mask_by_name(
    const char* query,
    uint8_t* mask_out,
    uint32_t mask_count
) {
    if (mask_out == nullptr || mask_count == 0) {
        return 0;
    }

    std::memset(mask_out, 0, static_cast<size_t>(mask_count));

    if (!g_initialized) {
        return 0;
    }
    if (!ensureOmmDatabase()) {
        return -1;
    }

    const std::string queryText = trimAsciiWhitespace(query);
    if (queryText.empty()) {
        const uint32_t visibleCount = std::min(mask_count, static_cast<uint32_t>(g_satellites.size()));
        std::memset(mask_out, 1, static_cast<size_t>(visibleCount));
        return static_cast<int32_t>(visibleCount);
    }

    sqlite3_stmt* stmt = g_queryEntityByNameStmt;
    if (!bindNameQueryStatement(stmt, queryText, mask_count)) {
        return -1;
    }

    uint32_t visibleCount = 0;
    while (true) {
        const int rc = sqlite3_step(stmt);
        if (rc == SQLITE_DONE) {
            break;
        }
        if (rc != SQLITE_ROW) {
            sqlite3_reset(stmt);
            sqlite3_clear_bindings(stmt);
            return -1;
        }
        const int entityIndex = sqlite3_column_int(stmt, 0);
        if (
            entityIndex >= 0 &&
            static_cast<uint32_t>(entityIndex) < mask_count &&
            mask_out[entityIndex] == 0
        ) {
            mask_out[entityIndex] = 1;
            visibleCount++;
        }
    }

    sqlite3_reset(stmt);
    sqlite3_clear_bindings(stmt);
    return static_cast<int32_t>(visibleCount);
}

// -----------------------------------------------------------------------------
// plugin_query_entity_rows_by_name — Query searchable rows with names.
// Empty query returns all rows, sorted by entity_index.
// Returns: number of rows written, or -1 on error.
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_query_entity_rows_by_name(
    const char* query,
    OrbProCatalogRow* out_rows,
    uint32_t max_count
) {
    if (!g_initialized || out_rows == nullptr || max_count == 0) {
        return 0;
    }
    if (!ensureOmmDatabase()) {
        return -1;
    }

    const std::string queryText = trimAsciiWhitespace(query);
    sqlite3_stmt* stmt = nullptr;

    if (queryText.empty()) {
        stmt = g_queryEntityRowsAllStmt;
        sqlite3_reset(stmt);
        sqlite3_clear_bindings(stmt);
        const int sqliteLimit =
            max_count > 0x7fffffffU ? 0x7fffffff : static_cast<int>(max_count);
        if (sqlite3_bind_int(stmt, 1, sqliteLimit) != SQLITE_OK) {
            sqlite3_reset(stmt);
            sqlite3_clear_bindings(stmt);
            return -1;
        }
    } else {
        stmt = g_queryEntityRowsByNameStmt;
        if (!bindNameQueryStatement(stmt, queryText, max_count)) {
            return -1;
        }
    }

    uint32_t written = 0;
    while (written < max_count) {
        const int rc = sqlite3_step(stmt);
        if (rc == SQLITE_DONE) {
            break;
        }
        if (rc != SQLITE_ROW) {
            sqlite3_reset(stmt);
            sqlite3_clear_bindings(stmt);
            return -1;
        }
        fillCatalogRowFromSqlite(stmt, out_rows[written++]);
    }

    sqlite3_reset(stmt);
    sqlite3_clear_bindings(stmt);
    return static_cast<int32_t>(written);
}

// -----------------------------------------------------------------------------
// plugin_get_entity_catalog_row — Get one catalog row by entity index.
// Returns: 0 on success, -1 on failure/not found.
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_get_entity_catalog_row(uint32_t entity_index, OrbProCatalogRow* out_row) {
    if (!g_initialized || out_row == nullptr) {
        return -1;
    }
    if (!ensureOmmDatabase() || g_getEntityCatalogRowStmt == nullptr) {
        return -1;
    }

    sqlite3_stmt* stmt = g_getEntityCatalogRowStmt;
    sqlite3_reset(stmt);
    sqlite3_clear_bindings(stmt);
    if (sqlite3_bind_int(stmt, 1, static_cast<int>(entity_index)) != SQLITE_OK) {
        sqlite3_reset(stmt);
        sqlite3_clear_bindings(stmt);
        return -1;
    }

    const int rc = sqlite3_step(stmt);
    if (rc != SQLITE_ROW) {
        sqlite3_reset(stmt);
        sqlite3_clear_bindings(stmt);
        return -1;
    }

    fillCatalogRowFromSqlite(stmt, *out_row);
    sqlite3_reset(stmt);
    sqlite3_clear_bindings(stmt);
    return 0;
}

// -----------------------------------------------------------------------------
// plugin_catalog_upsert_cat_record — Upsert CAT metadata row keyed by NORAD.
// Returns: 0 on success, -1 on failure.
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_catalog_begin_transaction() {
    return executeCatalogSql("BEGIN IMMEDIATE TRANSACTION;") ? 0 : -1;
}

ORBPRO_EXPORT
int32_t plugin_catalog_commit_transaction() {
    return executeCatalogSql("COMMIT;") ? 0 : -1;
}

ORBPRO_EXPORT
int32_t plugin_catalog_rollback_transaction() {
    return executeCatalogSql("ROLLBACK;") ? 0 : -1;
}

// -----------------------------------------------------------------------------
// plugin_catalog_upsert_cat_record — Upsert CAT metadata row keyed by NORAD.
// Returns: 0 on success, -1 on failure.
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_catalog_upsert_cat_record(
    uint32_t norad_cat_id,
    const char* object_name,
    const char* object_id,
    const char* payload_json
) {
    if (norad_cat_id == 0) {
        return -1;
    }
    const std::string name = object_name != nullptr ? std::string(object_name) : std::string();
    const std::string objectId = object_id != nullptr ? std::string(object_id) : std::string();
    const std::string payload = payload_json != nullptr ? std::string(payload_json) : std::string();
    return upsertCatCatalog(norad_cat_id, name, objectId, payload) ? 0 : -1;
}

// -----------------------------------------------------------------------------
// plugin_catalog_upsert_cat_records_packed — Upsert many CAT rows in one call.
//
// Packed byte format (little-endian, repeated until end of buffer):
//   u32 norad_cat_id
//   u32 object_name_len
//   u32 object_id_len
//   u32 payload_json_len
//   u8[object_name_len] object_name UTF-8 bytes (no NUL)
//   u8[object_id_len] object_id UTF-8 bytes (no NUL)
//   u8[payload_json_len] payload JSON UTF-8 bytes (no NUL)
//
// Returns: number of records upserted, or -1 on malformed input/fatal error.
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_catalog_upsert_cat_records_packed(
    const uint8_t* packed_bytes,
    uint32_t packed_len
) {
    if (packed_bytes == nullptr || packed_len == 0) {
        return 0;
    }
    if (!ensureOmmDatabase() || g_upsertCatStmt == nullptr) {
        return -1;
    }

    auto readU32Le = [](const uint8_t* ptr) -> uint32_t {
        return static_cast<uint32_t>(ptr[0]) |
               (static_cast<uint32_t>(ptr[1]) << 8) |
               (static_cast<uint32_t>(ptr[2]) << 16) |
               (static_cast<uint32_t>(ptr[3]) << 24);
    };

    bool inTransaction = false;
    if (plugin_catalog_begin_transaction() == 0) {
        inTransaction = true;
    }

    size_t offset = 0;
    int32_t upserted = 0;
    const size_t totalLen = static_cast<size_t>(packed_len);

    while (offset + 16 <= totalLen) {
        const uint8_t* header = packed_bytes + offset;
        const uint32_t noradCatId = readU32Le(header + 0);
        const uint32_t objectNameLen = readU32Le(header + 4);
        const uint32_t objectIdLen = readU32Le(header + 8);
        const uint32_t payloadLen = readU32Le(header + 12);
        offset += 16;

        const size_t fieldBytes =
            static_cast<size_t>(objectNameLen) +
            static_cast<size_t>(objectIdLen) +
            static_cast<size_t>(payloadLen);
        if (offset + fieldBytes > totalLen) {
            if (inTransaction) {
                plugin_catalog_rollback_transaction();
            }
            return -1;
        }

        const std::string objectName(
            reinterpret_cast<const char*>(packed_bytes + offset),
            static_cast<size_t>(objectNameLen)
        );
        offset += static_cast<size_t>(objectNameLen);

        const std::string objectId(
            reinterpret_cast<const char*>(packed_bytes + offset),
            static_cast<size_t>(objectIdLen)
        );
        offset += static_cast<size_t>(objectIdLen);

        const std::string payload(
            reinterpret_cast<const char*>(packed_bytes + offset),
            static_cast<size_t>(payloadLen)
        );
        offset += static_cast<size_t>(payloadLen);

        if (noradCatId == 0) {
            continue;
        }
        if (!upsertCatCatalog(noradCatId, objectName, objectId, payload)) {
            if (inTransaction) {
                plugin_catalog_rollback_transaction();
            }
            return upserted;
        }
        upserted++;
    }

    if (offset != totalLen) {
        if (inTransaction) {
            plugin_catalog_rollback_transaction();
        }
        return -1;
    }

    if (inTransaction && plugin_catalog_commit_transaction() != 0) {
        plugin_catalog_rollback_transaction();
        return upserted;
    }

    return upserted;
}

// -----------------------------------------------------------------------------
// plugin_catalog_upsert_cat_flatbuffer_stream — Upsert CAT rows directly from
// `$CAT`, `$REC`, or size-prefixed `$CAT`/`$REC` streams.
// Returns: number of records upserted, or -1 on malformed input/fatal error.
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_catalog_upsert_cat_flatbuffer_stream(
    const uint8_t* data,
    uint32_t len
) {
    if (data == nullptr || len == 0) {
        return 0;
    }
    if (!ensureOmmDatabase() || g_upsertCatStmt == nullptr) {
        return -1;
    }

    std::vector<PendingCatCatalogRecord> records;
    if (!parseCatFlatBufferStream(data, static_cast<size_t>(len), records)) {
        return -1;
    }
    if (records.empty()) {
        return 0;
    }

    bool inTransaction = false;
    if (plugin_catalog_begin_transaction() == 0) {
        inTransaction = true;
    }

    int32_t upserted = 0;
    for (const PendingCatCatalogRecord& record : records) {
        if (record.noradCatId == 0) {
            continue;
        }
        if (!upsertCatCatalog(
                record.noradCatId,
                record.objectName,
                record.objectId,
                std::string())) {
            if (inTransaction) {
                plugin_catalog_rollback_transaction();
            }
            return upserted;
        }
        upserted++;
    }

    if (inTransaction && plugin_catalog_commit_transaction() != 0) {
        plugin_catalog_rollback_transaction();
        return upserted;
    }

    return upserted;
}

// -----------------------------------------------------------------------------
// plugin_get_cat_record_json_size — Get CAT JSON payload bytes (+NUL).
// Returns: required bytes including trailing NUL, or 0 if not found.
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_get_cat_record_json_size(uint32_t norad_cat_id) {
    std::string json;
    if (!readCatPayloadJson(norad_cat_id, json)) {
        return 0;
    }
    return static_cast<int32_t>(json.size() + 1);
}

// -----------------------------------------------------------------------------
// plugin_get_cat_record_json — Copy CAT JSON payload into caller buffer.
// Returns: bytes written (excluding trailing NUL), 0 if not found, or -1 on error.
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_get_cat_record_json(
    uint32_t norad_cat_id,
    char* out_json,
    uint32_t out_len
) {
    if (out_json == nullptr || out_len == 0) {
        return -1;
    }
    std::string json;
    if (!readCatPayloadJson(norad_cat_id, json)) {
        out_json[0] = '\0';
        return 0;
    }
    if (json.size() + 1 > static_cast<size_t>(out_len)) {
        return -1;
    }
    std::memcpy(out_json, json.c_str(), json.size() + 1);
    return static_cast<int32_t>(json.size());
}

// -----------------------------------------------------------------------------
// plugin_get_omm_record_by_entity_index — Resolve the active OMM payload for an
// entity index.
// Returns: 0 on success, -1 on failure.
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_get_omm_record_by_entity_index(uint32_t entity_index, OrbProOMMRecord* out_record) {
    if (!g_initialized || out_record == nullptr || entity_index >= g_satellites.size()) {
        return -1;
    }

    const auto& entity = g_satellites[entity_index];
    if (entity.currentOmmPointer <= 0) {
        return -1;
    }

    return plugin_get_omm_record_by_pointer(
        static_cast<double>(entity.currentOmmPointer),
        out_record);
}

// -----------------------------------------------------------------------------
// plugin_get_omm_record_by_pointer — Resolve sqlite OMM pointer to payload.
// Returns: 0 on success, -1 on failure.
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_get_omm_record_by_pointer(double omm_pointer, OrbProOMMRecord* out_record) {
    if (!g_initialized || out_record == nullptr || !std::isfinite(omm_pointer) || omm_pointer <= 0.0) {
        return -1;
    }
    if (!ensureOmmDatabase() || g_getOmmPayloadByIdStmt == nullptr) {
        return -1;
    }

    const int64_t pointer = static_cast<int64_t>(omm_pointer);
    if (pointer <= 0) {
        return -1;
    }

    sqlite3_stmt* stmt = g_getOmmPayloadByIdStmt;
    sqlite3_reset(stmt);
    sqlite3_clear_bindings(stmt);
    if (sqlite3_bind_int64(stmt, 1, static_cast<sqlite3_int64>(pointer)) != SQLITE_OK) {
        sqlite3_reset(stmt);
        sqlite3_clear_bindings(stmt);
        return -1;
    }

    const int rc = sqlite3_step(stmt);
    if (rc != SQLITE_ROW) {
        sqlite3_reset(stmt);
        sqlite3_clear_bindings(stmt);
        return -1;
    }

    const void* blob = sqlite3_column_blob(stmt, 0);
    const int blobBytes = sqlite3_column_bytes(stmt, 0);
    if (blob == nullptr || blobBytes < static_cast<int>(sizeof(OrbProOMMRecord))) {
        sqlite3_reset(stmt);
        sqlite3_clear_bindings(stmt);
        return -1;
    }

    std::memcpy(out_record, blob, sizeof(OrbProOMMRecord));
    sqlite3_reset(stmt);
    sqlite3_clear_bindings(stmt);
    return 0;
}

ORBPRO_EXPORT
void plugin_destroy(void) {
    g_satellites.clear();
    g_segmentSets.clear();
    g_nextSegmentSetHandle = 1;
    closeOmmDatabase();
    g_initialized = false;
}

ORBPRO_EXPORT
uint32_t get_satellite_count(void) {
    return static_cast<uint32_t>(g_satellites.size());
}

ORBPRO_EXPORT
double get_orbital_period(uint32_t entity_index) {
    if (entity_index >= g_satellites.size()) return 0.0;
    const auto& e = g_satellites[entity_index];
    if (!e.valid) return 0.0;
    return TWOPI / e.satrec.no_unkozai;
}

ORBPRO_EXPORT
double get_tle_epoch(uint32_t entity_index) {
    if (entity_index >= g_satellites.size()) return 0.0;
    const auto& e = g_satellites[entity_index];
    return e.satrec.jdsatepoch + e.satrec.jdsatepochF;
}

// -----------------------------------------------------------------------------
// plugin_propagate_path_sv — Propagate one entity at many uniformly-spaced
// times, outputting full state vectors (pos+vel) in TEME frame, in km and km/s.
// Layout: 6 × float64 per sample [px, py, pz, vx, vy, vz]
// Returns: number of samples written, or negative error code
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_propagate_path_sv(uint32_t entity_index, double start_jd,
                                  double step_days, uint32_t count,
                                  double* sv_out) {
    if (!g_initialized || entity_index >= g_satellites.size()) return -1;

    SatelliteEntity& entity = g_satellites[entity_index];
    if (!entity.valid) return -1;

    for (uint32_t i = 0; i < count; i++) {
        double jd = start_jd + i * step_days;

        double r[3], v[3];
        bool ok = propagateEntityTEME(entity, jd, r, v);
        if (!ok) {
            double nan = 0.0 / 0.0;
            for (int j = 0; j < 6; j++) sv_out[i * 6 + j] = nan;
            continue;
        }

        sv_out[i * 6 + 0] = r[0];
        sv_out[i * 6 + 1] = r[1];
        sv_out[i * 6 + 2] = r[2];
        sv_out[i * 6 + 3] = v[0];
        sv_out[i * 6 + 4] = v[1];
        sv_out[i * 6 + 5] = v[2];
    }

    return static_cast<int32_t>(count);
}

// -----------------------------------------------------------------------------
// plugin_propagate_path_adaptive — True anomaly-based adaptive path sampling.
//
// Performs all orbital mechanics in WASM: propagates, computes eccentricity
// vector and true anomaly, steps by angular increments, injects epoch
// hypersamples, and outputs ECEF positions (meters) ready for rendering.
//
// Parameters:
//   entity_index    — satellite index
//   start_jd        — path start Julian Date
//   stop_jd         — path stop Julian Date
//   update_jd       — current clock time (for epoch hypersampling)
//   samples_per_period — angular samples per orbit (e.g. 360)
//   extra_epoch_samples — dense samples near update_jd (e.g. 30)
//   resolution_at_epoch — time step (seconds) for epoch hypersamples
//   max_samples     — output buffer capacity
//   positions_out   — output: 3 × float64 per sample [x, y, z] ECEF meters
//
// Returns: number of positions written
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_propagate_path_adaptive(
    uint32_t entity_index,
    double start_jd, double stop_jd, double update_jd,
    uint32_t samples_per_period,
    uint32_t extra_epoch_samples,
    double resolution_at_epoch,
    uint32_t max_samples,
    double* positions_out)
{
    if (!g_initialized || entity_index >= g_satellites.size()) return -1;

    SatelliteEntity& entity = g_satellites[entity_index];
    if (!entity.valid) return -1;

    // For multi-OMM (>1), ensure entity.satrec is set to the nearest for start_jd
    // so period/mus/no_unkozai are correct for stepping parameters.
    // Skip for single-OMM entities — entity.satrec is already correct and
    // the map lookup + struct copy would be wasteful every frame.
    if (entity.satrecs.size() > 1) {
        double nearest = findSatrecIndex(entity.satrecs, start_jd);
        if (nearest >= 0) {
            entity.satrec = entity.satrecs[nearest];
            entity.currentSatrecEpoch = nearest;
        }
    }

    double periodSec = TWOPI / entity.satrec.no_unkozai;  // minutes → need *60
    periodSec *= 60.0;
    double periodDays = periodSec / 86400.0;

    double nominalStep = TWOPI / static_cast<double>(samples_per_period);

    // Epoch hypersampling window (in days)
    double epochResolutionDays = resolution_at_epoch / 86400.0;
    double epochWindowHalf = (extra_epoch_samples / 2.0) * epochResolutionDays;
    double epochBegin = update_jd - epochWindowHalf;
    bool doneWithEpoch = (update_jd <= start_jd || update_jd >= stop_jd);

    // Gravity parameter for rv2coe (km³/s²)
    double mus = entity.satrec.mus;

    uint32_t written = 0;
    double currentJD = start_jd;
    double lastTrueAnom = -999.0;

    auto writeEcefPosition = [&](double jd) -> bool {
        if (written >= max_samples) return false;
        double r[3], v[3];
        bool ok = propagateEntityTEME(entity, jd, r, v);
        if (!ok) return true; // skip but continue

        double gmst = SGP4Funcs::gstime_SGP4(jd);
        double cosG = cos(gmst);
        double sinG = sin(gmst);

        positions_out[written * 3 + 0] = ( cosG * r[0] + sinG * r[1]) * 1000.0;
        positions_out[written * 3 + 1] = (-sinG * r[0] + cosG * r[1]) * 1000.0;
        positions_out[written * 3 + 2] = r[2] * 1000.0;
        written++;
        return true;
    };

    while (currentJD < stop_jd && written < max_samples) {
        // Inject epoch hypersamples before we pass the epoch window
        if (!doneWithEpoch && currentJD >= epochBegin) {
            for (uint32_t i = 0; i < extra_epoch_samples && written < max_samples; i++) {
                double epochJD = epochBegin + i * epochResolutionDays;
                if (epochJD >= start_jd && epochJD <= stop_jd) {
                    writeEcefPosition(epochJD);
                }
            }
            doneWithEpoch = true;
            double epochEnd = epochBegin + extra_epoch_samples * epochResolutionDays;
            if (epochEnd > currentJD) {
                currentJD = epochEnd;
            }
            continue;
        }

        // Propagate at current time
        double r[3], v[3];
        bool ok = propagateEntityTEME(entity, currentJD, r, v);
        if (!ok) {
            currentJD += periodDays / samples_per_period;
            continue;
        }

        // Write ECEF position
        {
            double gmst = SGP4Funcs::gstime_SGP4(currentJD);
            double cosG = cos(gmst);
            double sinG = sin(gmst);
            if (written < max_samples) {
                positions_out[written * 3 + 0] = ( cosG * r[0] + sinG * r[1]) * 1000.0;
                positions_out[written * 3 + 1] = (-sinG * r[0] + cosG * r[1]) * 1000.0;
                positions_out[written * 3 + 2] = r[2] * 1000.0;
                written++;
            }
        }

        // Compute true anomaly and eccentricity from state vector (TEME)
        double p, a, ecc, incl, omega, argp, nu, m, arglat, truelon, lonper;
        SGP4Funcs::rv2coe_SGP4(r, v, mus, p, a, ecc, incl, omega, argp,
                                nu, m, arglat, truelon, lonper);

        if (nu != nu || ecc != ecc || ecc > 0.99) {
            currentJD += periodDays / samples_per_period;
            lastTrueAnom = -999.0;
            continue;
        }

        // Adaptive step correction for true anomaly drift
        double step1 = nominalStep;
        if (lastTrueAnom > -900.0) {
            double actualStep = fabs(nu - lastTrueAnom);
            if (actualStep > PI) actualStep = TWOPI - actualStep;
            step1 = 2.0 * nominalStep - actualStep;
            if (step1 <= 0.0 || !isfinite(step1)) step1 = nominalStep;
        }

        double nextNu = nu + step1;
        if (nextNu > TWOPI) nextNu -= TWOPI;

        lastTrueAnom = nu;

        // Convert true anomaly step to mean anomaly step, then to time
        double e0_old, m_old, e0_new, m_new;
        SGP4Funcs::newtonnu_SGP4(ecc, nu, e0_old, m_old);
        SGP4Funcs::newtonnu_SGP4(ecc, nextNu, e0_new, m_new);

        double dm = fabs(m_new - m_old);
        if (dm > PI) dm = TWOPI - dm;

        double timeEstimateSec = periodSec * (dm / TWOPI);

        double maxStepSec = periodSec / 2.0;
        double minStepSec = periodSec / (samples_per_period * 4.0);
        if (!isfinite(timeEstimateSec) || timeEstimateSec <= 0.0) {
            timeEstimateSec = periodSec / samples_per_period;
        } else if (timeEstimateSec > maxStepSec) {
            timeEstimateSec = maxStepSec;
        } else if (timeEstimateSec < minStepSec) {
            timeEstimateSec = minStepSec;
        }

        currentJD += timeEstimateSec / 86400.0;
    }

    // Final sample at stop
    if (written < max_samples && currentJD >= stop_jd) {
        writeEcefPosition(stop_jd);
    }

    return static_cast<int32_t>(written);
}

// -----------------------------------------------------------------------------
// plugin_propagate_all_positions — Propagate ALL entities to one Julian date.
// Output is positions-only (3 × float64 per entity, ECEF meters), indexed by
// entity_index. Designed for bulk scatter-write into the wasm-engine position
// buffer. Single WASM call for all satellites.
// Returns: number of entities written
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_propagate_all_positions(double julian_date,
                                        double* positions_out,
                                        uint32_t max_count) {
    if (!g_initialized) return 0;

    uint32_t n = static_cast<uint32_t>(g_satellites.size());
    if (n > max_count) n = max_count;

    double gmst = SGP4Funcs::gstime_SGP4(julian_date);
    double cosG = cos(gmst);
    double sinG = sin(gmst);

    for (uint32_t i = 0; i < n; i++) {
        SatelliteEntity& entity = g_satellites[i];
        double* out = positions_out + i * 3;

        if (!entity.valid) {
            out[0] = out[1] = out[2] = 0.0;
            continue;
        }

        double r[3], v[3];
        bool ok = propagateEntityTEME(entity, julian_date, r, v);
        if (!ok) {
            out[0] = out[1] = out[2] = 0.0;
            continue;
        }

        out[0] = ( cosG * r[0] + sinG * r[1]) * 1000.0;
        out[1] = (-sinG * r[0] + cosG * r[1]) * 1000.0;
        out[2] = r[2] * 1000.0;

        entity.lastEpochJD = julian_date;
    }

    return static_cast<int32_t>(n);
}

// -----------------------------------------------------------------------------
// plugin_propagate_to_registry — Scatter-write ALL entity positions and
// finite-differenced velocities directly into the wasm-engine's EntityRegistry
// buffers using a handle map.
// Each satellite_index maps to a registry handle via handle_map[i].
// Positions: positions_out[handle_map[i] * 3 + {0,1,2}]
// Velocities: velocities_out[handle_map[i] * 3 + {0,1,2}] (ECEF m/s)
// velocities_out may be null (backwards compatible — positions only).
// GMST is computed once and shared across all entities.
// Returns: number of entities written
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_propagate_to_registry(
    double julian_date,
    double* positions_out,
    const uint32_t* handle_map,
    uint32_t count,
    double* velocities_out
) {
    if (!g_initialized || !positions_out || !handle_map) return 0;

    uint32_t n = count;
    if (n > static_cast<uint32_t>(g_satellites.size()))
        n = static_cast<uint32_t>(g_satellites.size());

    double gmst = SGP4Funcs::gstime_SGP4(julian_date);

    uint32_t written = 0;
    for (uint32_t i = 0; i < n; i++) {
        SatelliteEntity& entity = g_satellites[i];
        uint32_t handle = handle_map[i];
        double* posOut = positions_out + handle * 3;

        if (!entity.valid) {
            posOut[0] = posOut[1] = posOut[2] = 0.0;
            if (velocities_out) {
                double* velOut = velocities_out + handle * 3;
                velOut[0] = velOut[1] = velOut[2] = 0.0;
            }
            continue;
        }

        double r[3], v[3];
        bool ok = propagateEntityTEME(entity, julian_date, r, v);
        if (!ok) {
            posOut[0] = posOut[1] = posOut[2] = 0.0;
            if (velocities_out) {
                double* velOut = velocities_out + handle * 3;
                velOut[0] = velOut[1] = velOut[2] = 0.0;
            }
            continue;
        }

        double posEcef[3];
        double velEcef[3];
        temeToEcef(r, v, gmst, posEcef, velEcef);

        posOut[0] = posEcef[0];
        posOut[1] = posEcef[1];
        posOut[2] = posEcef[2];

        if (velocities_out) {
            double* velOut = velocities_out + handle * 3;
            velOut[0] = velEcef[0];
            velOut[1] = velEcef[1];
            velOut[2] = velEcef[2];
        }

        entity.lastEpochJD = julian_date;
        written++;
    }

    return static_cast<int32_t>(written);
}

// -----------------------------------------------------------------------------
// get_eccentricity — Return eccentricity for a satellite
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
double get_eccentricity(uint32_t entity_index) {
    if (entity_index >= g_satellites.size()) return -1.0;
    return g_satellites[entity_index].satrec.ecco;
}

// =============================================================================
// Multi-OMM Management — per-entity OMM add/remove/list/mode
// =============================================================================

// -----------------------------------------------------------------------------
// plugin_entity_add_omm — Add an additional OMM to an existing entity.
// Returns: epoch JD of the added OMM on success, -1.0 on error
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
double plugin_entity_add_omm(uint32_t entity_index, const OrbProOMMRecord* record) {
    if (!g_initialized || entity_index >= g_satellites.size() || !record) return -1.0;

    SatelliteEntity& entity = g_satellites[entity_index];

    elsetrec newSatrec;
    if (!initSatrecFromOMM(*record, newSatrec)) {
        return -1.0;
    }

    double epochJD = newSatrec.jdsatepoch + newSatrec.jdsatepochF;
    int64_t ommPointer = insertOmmRecord(*record, epochJD);
    if (ommPointer <= 0) {
        return -1.0;
    }

    // Seed the map with the existing single satrec before adding the new one
    if (entity.satrecs.empty()) {
        entity.satrecs[entity.currentSatrecEpoch] = entity.satrec;
        if (entity.currentOmmPointer > 0) {
            entity.ommPointers[entity.currentSatrecEpoch] = entity.currentOmmPointer;
        }
    }
    entity.satrecs[epochJD] = newSatrec;
    entity.ommPointers[epochJD] = ommPointer;

    // Only update entity.satrec if this is the most recent epoch
    if (epochJD > entity.currentSatrecEpoch) {
        entity.satrec = newSatrec;
        entity.currentSatrecEpoch = epochJD;
        entity.currentOmmPointer = ommPointer;
    }
    entity.valid = true;

    return epochJD;
}

// -----------------------------------------------------------------------------
// plugin_entity_add_omm_flatbuffer — Add a single `$OMM` payload to an entity.
// Returns: epoch JD on success, -1.0 on error
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
double plugin_entity_add_omm_flatbuffer(
    uint32_t entity_index,
    const uint8_t* payload,
    size_t payload_len
) {
    PendingCatalogMetadata metadata;
    OrbProOMMRecord record;
    if (!parseOmmFlatBufferPayload(payload, payload_len, record, &metadata)) {
        return -1.0;
    }
    const double epochJD = plugin_entity_add_omm(entity_index, &record);
    if (
        epochJD >= 0.0 &&
        (!metadata.objectName.empty() || !metadata.objectId.empty()) &&
        entity_index < static_cast<uint32_t>(g_satellites.size())
    ) {
        const auto& entity = g_satellites[entity_index];
        upsertEntityCatalog(
            entity_index,
            entity.noradId,
            metadata.objectName,
            metadata.objectId
        );
    }
    return epochJD;
}

// -----------------------------------------------------------------------------
// plugin_entity_omm_count — Get the number of OMMs stored for an entity
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
uint32_t plugin_entity_omm_count(uint32_t entity_index) {
    if (!g_initialized || entity_index >= g_satellites.size()) return 0;
    const auto& entity = g_satellites[entity_index];
    // satrecs map is empty for single-OMM entities (fast path optimization)
    return entity.satrecs.empty() ? (entity.valid ? 1 : 0)
                                  : static_cast<uint32_t>(entity.satrecs.size());
}

// -----------------------------------------------------------------------------
// plugin_entity_list_omm — List all OMM epoch JDs for an entity.
// Writes epoch JDs to epochs_out (sorted ascending). Returns count written.
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_entity_list_omm(uint32_t entity_index, double* epochs_out,
                                 uint32_t max_count) {
    if (!g_initialized || entity_index >= g_satellites.size() || !epochs_out) return 0;

    const auto& entity = g_satellites[entity_index];
    if (entity.satrecs.empty()) {
        // Single-OMM fast path: map is empty, report the sole epoch
        if (max_count > 0 && entity.valid) {
            epochs_out[0] = entity.currentSatrecEpoch;
            return 1;
        }
        return 0;
    }
    uint32_t written = 0;
    for (const auto& pair : entity.satrecs) {
        if (written >= max_count) break;
        epochs_out[written++] = pair.first;
    }
    return static_cast<int32_t>(written);
}

// -----------------------------------------------------------------------------
// plugin_entity_remove_omm — Remove a single OMM by its epoch JD
// Returns: 0 on success, -1 if not found
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_entity_remove_omm(uint32_t entity_index, double epoch_jd) {
    if (!g_initialized || entity_index >= g_satellites.size()) return -1;

    auto& entity = g_satellites[entity_index];
    auto it = entity.satrecs.find(epoch_jd);
    if (it == entity.satrecs.end()) return -1;

    entity.satrecs.erase(it);
    entity.ommPointers.erase(epoch_jd);

    // If removed the current satrec, switch to first remaining
    if (!entity.satrecs.empty() && entity.currentSatrecEpoch == epoch_jd) {
        auto first = entity.satrecs.begin();
        entity.satrec = first->second;
        entity.currentSatrecEpoch = first->first;
        auto ptrIt = entity.ommPointers.find(first->first);
        entity.currentOmmPointer = ptrIt != entity.ommPointers.end() ? ptrIt->second : -1;
    } else if (entity.satrecs.empty()) {
        entity.currentOmmPointer = -1;
    }

    return 0;
}

// -----------------------------------------------------------------------------
// plugin_entity_remove_all_omm_except — Keep only the OMM at epoch_jd
// Returns: 0 on success, -1 if not found
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_entity_remove_all_omm_except(uint32_t entity_index, double epoch_jd) {
    if (!g_initialized || entity_index >= g_satellites.size()) return -1;

    auto& entity = g_satellites[entity_index];
    auto it = entity.satrecs.find(epoch_jd);
    if (it == entity.satrecs.end()) return -1;
    auto ptrIt = entity.ommPointers.find(epoch_jd);
    int64_t keeperPointer = ptrIt != entity.ommPointers.end() ? ptrIt->second : -1;

    elsetrec keeper = it->second;
    entity.satrecs.clear();
    entity.ommPointers.clear();
    entity.satrecs[epoch_jd] = keeper;
    if (keeperPointer > 0) {
        entity.ommPointers[epoch_jd] = keeperPointer;
    }
    entity.satrec = keeper;
    entity.currentSatrecEpoch = epoch_jd;
    entity.currentOmmPointer = keeperPointer;

    return 0;
}

// -----------------------------------------------------------------------------
// plugin_entity_clear_omm — Clear all OMMs for an entity.
// Reverts to single-satrec mode using the current satrec.
// Returns: 0 on success
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_entity_clear_omm(uint32_t entity_index) {
    if (!g_initialized || entity_index >= g_satellites.size()) return -1;

    auto& entity = g_satellites[entity_index];
    entity.satrecs.clear();
    entity.ommPointers.clear();
    entity.currentSatrecEpoch = 0.0;
    entity.currentOmmPointer = -1;
    entity.mode = MODE_NEAREST_EPOCH;

    return 0;
}

// -----------------------------------------------------------------------------
// plugin_entity_set_mode — Set propagation mode for a multi-OMM entity.
//   0 = nearest_epoch (default): use closest OMM for propagation
//   1 = interpolated_epoch: blend between two bracketing OMMs
// Returns: 0 on success, -1 on error
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
int32_t plugin_entity_set_mode(uint32_t entity_index, uint32_t mode) {
    if (!g_initialized || entity_index >= g_satellites.size()) return -1;
    if (mode > MODE_INTERPOLATED_EPOCH) return -1;

    g_satellites[entity_index].mode = static_cast<uint8_t>(mode);
    return 0;
}

// -----------------------------------------------------------------------------
// plugin_entity_get_mode — Get propagation mode for an entity.
// Returns: mode (0 or 1), or 0xFF on error
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
uint32_t plugin_entity_get_mode(uint32_t entity_index) {
    if (!g_initialized || entity_index >= g_satellites.size()) return 0xFF;
    return g_satellites[entity_index].mode;
}

// -----------------------------------------------------------------------------
// plugin_entity_get_omm_pointer — Get active OMM sqlite row pointer for entity.
// Returns: sqlite row id (>0), or -1 when unavailable.
// -----------------------------------------------------------------------------
ORBPRO_EXPORT
double plugin_entity_get_omm_pointer(uint32_t entity_index) {
    if (!g_initialized || entity_index >= g_satellites.size()) return -1.0;
    const auto& entity = g_satellites[entity_index];
    return entity.currentOmmPointer > 0
        ? static_cast<double>(entity.currentOmmPointer)
        : -1.0;
}

}  // extern "C"
