/**
 * @orbpro/plugin-sgp4
 * OrbPro SGP4/SDP4 Orbital Propagator Plugin
 *
 * Accepts pre-decrypted WASM bytes delivered via the SDN plugin-delivery system
 * (ecies-x25519-hkdf-sha256-aes-256-gcm), or loads raw WASM from dist/ directly.
 */

import { fileURLToPath } from "node:url";
import { readFileSync } from "node:fs";
import path from "node:path";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

/**
 * Load the raw WASM bytes from the package dist/ directory.
 * @param {string} filename
 * @returns {Uint8Array}
 */
function loadPackageWasmBytes(filename) {
  const wasmPath = path.join(__dirname, "dist", filename);
  return readFileSync(wasmPath);
}

/**
 * Create an SGP4/SDP4 propagator instance.
 *
 * @param {object} options
 * @param {Uint8Array} [options.wasmBytes]  Pre-decrypted WASM bytes (from ecies delivery).
 *                                           If omitted, loads from package dist/.
 * @param {Function}  [options.decryptFn]   Legacy AES-256-GCM decrypt function (protection-runtime).
 * @param {boolean}   [options.lowMemory]   Use reduced memory configuration.
 * @returns {Promise<object>} Plugin instance
 */
export async function createSGP4Propagator(options = {}) {
  // Resolve WASM bytes: ecies-delivered > decryptFn > raw from dist
  let wasmBytes = options.wasmBytes;

  if (!wasmBytes && typeof options.decryptFn === "function") {
    // Legacy path: use protection-runtime decryptFn
    const { encryptedData } = await import("./dist/sgp4-encrypted.js").catch(() => ({}));
    if (encryptedData) {
      const result = await options.decryptFn(encryptedData);
      wasmBytes = result?.data ?? result;
    }
  }

  if (!wasmBytes) {
    // Development/direct path: load raw WASM from dist
    wasmBytes = loadPackageWasmBytes("sgp4.wasm");
  }

  // Delegate to OrbPro's internal loader (when installed in OrbPro monorepo)
  // This import is resolved via the OrbPro workspace during development
  const loader = await import("./loader.js").catch(() => null);
  if (loader?.createSGP4Propagator) {
    return loader.createSGP4Propagator({ ...options, wasmBytes });
  }

  // Standalone mode: return minimal plugin interface
  throw new Error(
    "createSGP4Propagator requires the OrbPro engine for full functionality. Install via the OrbPro workspace."
  );
}
