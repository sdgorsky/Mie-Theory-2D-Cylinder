import { defineConfig } from 'vite'
import type { Plugin } from 'vite'
import react from '@vitejs/plugin-react'
import wasm from 'vite-plugin-wasm'
import topLevelAwait from 'vite-plugin-top-level-await'
import path from 'path'
import { readFileSync } from 'fs'

const pkg = JSON.parse(readFileSync('./package.json', 'utf-8'))

// Plugin to provide "env" module for WASM imports (math functions from libm)
function wasmEnvPlugin(): Plugin {
  const envModuleId = 'env'
  const resolvedEnvModuleId = '\0' + envModuleId

  return {
    name: 'wasm-env',
    resolveId(id) {
      if (id === envModuleId) {
        return resolvedEnvModuleId
      }
    },
    load(id) {
      if (id === resolvedEnvModuleId) {
        // Provide math functions that WASM modules compiled from C/Fortran may need
        return `
          export const sin = Math.sin;
          export const cos = Math.cos;
          export const tan = Math.tan;
          export const asin = Math.asin;
          export const acos = Math.acos;
          export const atan = Math.atan;
          export const atan2 = Math.atan2;
          export const sinh = Math.sinh;
          export const cosh = Math.cosh;
          export const tanh = Math.tanh;
          export const exp = Math.exp;
          export const log = Math.log;
          export const log10 = Math.log10;
          export const pow = Math.pow;
          export const sqrt = Math.sqrt;
          export const floor = Math.floor;
          export const ceil = Math.ceil;
          export const fabs = Math.abs;
          export const abs = Math.abs;
          export const fmod = (x, y) => x % y;
          export const round = Math.round;
        `
      }
    },
  }
}

// https://vite.dev/config/
export default defineConfig({
  base: process.env.CI ? '/Mie-Theory-2D-Cylinder/' : '/',
  define: {
    __APP_VERSION__: JSON.stringify(pkg.version),
  },
  plugins: [
    react(),
    wasm(),
    topLevelAwait(),
    wasmEnvPlugin(),
  ],
  resolve: {
    alias: {
      'scattering-core': path.resolve(__dirname, 'scattering-core/pkg'),
    },
  },
  optimizeDeps: {
    exclude: ['scattering-core'],
  },
})
