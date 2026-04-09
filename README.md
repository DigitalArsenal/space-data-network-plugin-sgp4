# @orbpro/plugin-sgp4

OrbPro SGP4/SDP4 Orbital Propagator Plugin.

Implements the SGP4 and SDP4 algorithms for propagating satellite state vectors from Two-Line Element (TLE) sets. Compiled to WebAssembly for high-performance, cross-platform use.

## Installation

```bash
npm install @orbpro/plugin-sgp4
```

This package is intended to be used within an OrbPro workspace or alongside the OrbPro engine. Standalone use requires the OrbPro plugin-sdk.

## Building

Build the WASM artifact using the OrbPro plugin-sdk:

```bash
# From within the OrbPro plugin-sdk workspace:
npm run build:sgp4
# Output: dist/sgp4.wasm
```

## Usage

### Via SDN Plugin Delivery (ecies-decrypted bytes)

```javascript
import { createSGP4Propagator } from "@orbpro/plugin-sgp4";

// wasmBytes are delivered pre-decrypted by the SDN plugin-delivery system
// (ecies-x25519-hkdf-sha256-aes-256-gcm)
const propagator = await createSGP4Propagator({ wasmBytes });
```

### Direct / Development

```javascript
import { createSGP4Propagator } from "@orbpro/plugin-sgp4";

// Without wasmBytes, loads raw WASM from dist/sgp4.wasm
const propagator = await createSGP4Propagator();
```

### Options

| Option | Type | Description |
|--------|------|-------------|
| `wasmBytes` | `Uint8Array` | Pre-decrypted WASM bytes from the SDN delivery system. |
| `decryptFn` | `Function` | Legacy AES-256-GCM decrypt function (protection-runtime). |
| `lowMemory` | `boolean` | Use reduced memory configuration. |

## License

UNLICENSED — Proprietary. All rights reserved by DigitalArsenal.io, Inc.
