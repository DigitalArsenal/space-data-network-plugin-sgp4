/**
 * @orbpro/plugin-sgp4
 * OrbPro SGP4/SDP4 Orbital Propagator Plugin
 *
 * Delivered via SDN ecies-x25519-hkdf-sha256-aes-256-gcm artifact envelope,
 * or loads raw WASM from dist/ for local development.
 *
 * The Emscripten WASM exports the full SGP4/SDP4 propagation engine and
 * implements the OrbPro plugin_stream_invoke ABI.
 */

import { fileURLToPath } from "node:url";
import { readFileSync } from "node:fs";
import path from "node:path";

const __dirname = path.dirname(fileURLToPath(import.meta.url));

function loadPackageWasmBytes(filename) {
  const wasmPath = path.join(__dirname, "dist", filename);
  return readFileSync(wasmPath);
}

/**
 * Create an SGP4/SDP4 propagator instance.
 *
 * @param {object}    options
 * @param {Uint8Array} [options.wasmBytes]  Decrypted WASM bytes from SDN ecies delivery.
 *                                           If omitted, loads from package dist/.
 * @param {boolean}   [options.lowMemory]   Use reduced memory configuration.
 * @returns {Promise<object>} Propagator instance
 */
export async function createSGP4Propagator(options = {}) {
  let wasmBytes = options.wasmBytes;

  if (!wasmBytes) {
    wasmBytes = loadPackageWasmBytes("sgp4.wasm");
  }

  // Delegate to OrbPro's internal loader when installed in the OrbPro workspace
  const loader = await import("./loader.js").catch(() => null);
  if (loader?.createSGP4Propagator) {
    return loader.createSGP4Propagator({ ...options, wasmBytes });
  }

  throw new Error(
    "createSGP4Propagator requires the OrbPro engine. Install via the OrbPro workspace.",
  );
}
