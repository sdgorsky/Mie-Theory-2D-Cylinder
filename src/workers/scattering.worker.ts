import init, {
  compute_all,
  get_field_grid_size,
  get_parameter_bounds,
} from "scattering-core";

// ── Types ────────────────────────────────────────────────────────────

export type VisualizationMode =
  | "magnitude"
  | "log_magnitude"
  | "real"
  | "imag"
  | "phase";

export interface ComputeParams {
  wavelength: number;
  permittivityRe: number;
  permittivityIm: number;
  permeabilityRe: number;
  permeabilityIm: number;
  sourceType: number; // 0=PlaneWave TM, 1=PlaneWave TE, 2=DipoleEz, 3=DipoleExy
  dipoleXs: number;
  dipoleYs: number;
  dipoleAlpha: number;
  maxOrder: number;
  viewSize: number;
}

export interface ImageStats {
  min: number;
  max: number;
}

export interface ParameterBounds {
  wavelength_min: number;
  wavelength_max: number;
  permittivity_re_min: number;
  permittivity_re_max: number;
  permittivity_im_min: number;
  permittivity_im_max: number;
  permeability_re_min: number;
  permeability_re_max: number;
  permeability_im_min: number;
  permeability_im_max: number;
  max_order_min: number;
  max_order_max: number;
}

export interface ScatteringMeta {
  numOrders: number;
  minOrder: number;
  maxOrder: number;
}

export type WorkerRequest =
  | { type: "init" }
  | {
      type: "compute";
      params: ComputeParams;
      mode: VisualizationMode;
      recycledBuffer?: ArrayBuffer;
    }
  | {
      type: "recolor";
      mode: VisualizationMode;
      recycledBuffer?: ArrayBuffer;
    };

export type WorkerResponse =
  | { type: "ready"; bounds: ParameterBounds }
  | {
      type: "result";
      scattering: ScatteringMeta;
      imageData: Uint8ClampedArray;
      gridSize: number;
      viewSize: number;
      stats: ImageStats;
    }
  | { type: "recolored"; imageData: Uint8ClampedArray; stats: ImageStats }
  | { type: "error"; message: string };

// ── Pre-allocated buffers (allocated once at init, never GC'd) ──────

// Field data buffers — passed to WASM, also used for recolor
let fieldRealBuf: Float64Array | null = null;
let fieldImagBuf: Float64Array | null = null;
let cachedGridSize = 0;

// Reused every frame — never transferred, never GC'd
let reusableValues: Float64Array | null = null;

// RGBA buffer recycled from main thread (ping-pong)
let recycledRgbaBuffer: ArrayBuffer | null = null;

// Cached compute params — needed by recolor for dipole exclusion zone
let cachedComputeParams: ComputeParams | null = null;

function recycleBuffer(buf: ArrayBuffer | undefined) {
  if (buf && buf.byteLength > 0) {
    recycledRgbaBuffer = buf;
  }
}

// ── Colormap functions ───────────────────────────────────────────────
// Write directly into RGBA buffer to avoid per-pixel tuple allocation.
// Uint8ClampedArray handles clamping and rounding automatically.

function writeViridis(rgba: Uint8ClampedArray, idx: number, t: number): void {
  if (t < 0) t = 0;
  else if (t > 1) t = 1;
  rgba[idx] =
    255 *
    (0.267004 + t * (0.329415 + t * (-0.814464 + t * (2.28653 - t * 1.06868))));
  rgba[idx + 1] =
    255 *
    (0.004874 +
      t * (0.873449 + t * (0.107514 + t * (-0.631923 + t * 0.645732))));
  rgba[idx + 2] =
    255 *
    (0.329415 + t * (1.01541 + t * (-1.67917 + t * (1.59456 - t * 0.594634))));
  rgba[idx + 3] = 255;
}

function writePhase(rgba: Uint8ClampedArray, idx: number, phase: number): void {
  const h = ((phase + Math.PI) / (2 * Math.PI)) * 360;
  const c = 0.8; // s=0.8, l=0.5 → c = (1 - |2*0.5 - 1|) * 0.8
  const x = c * (1 - Math.abs(((h / 60) % 2) - 1));
  const m = 0.1; // l - c/2 = 0.5 - 0.4
  let r = 0,
    g = 0,
    b = 0;
  if (h < 60) {
    r = c;
    g = x;
  } else if (h < 120) {
    r = x;
    g = c;
  } else if (h < 180) {
    g = c;
    b = x;
  } else if (h < 240) {
    g = x;
    b = c;
  } else if (h < 300) {
    r = x;
    b = c;
  } else {
    r = c;
    b = x;
  }
  rgba[idx] = (r + m) * 255;
  rgba[idx + 1] = (g + m) * 255;
  rgba[idx + 2] = (b + m) * 255;
  rgba[idx + 3] = 255;
}

function writeDiverging(rgba: Uint8ClampedArray, idx: number, t: number): void {
  if (t < 0) t = 0;
  else if (t > 1) t = 1;
  if (t < 0.5) {
    const s = t * 2;
    rgba[idx] = 59 + s * 196;
    rgba[idx + 1] = 76 + s * 179;
    rgba[idx + 2] = 192 + s * 63;
  } else {
    const s = (t - 0.5) * 2;
    rgba[idx] = 255 - s * 75;
    rgba[idx + 1] = 255 - s * 196;
    rgba[idx + 2] = 255 - s * 196;
  }
  rgba[idx + 3] = 255;
}

// ── Render field to RGBA ─────────────────────────────────────────────

/** Physical exclusion radius around dipole (in world units).
 *  DipoleExy diverges as H_1(kR)/R (stronger than DipoleEz's log singularity),
 *  so it needs a larger exclusion zone to avoid washing out the color range. */
const DIPOLE_EXCLUSION_RADIUS: Record<number, number> = {
  2: 0.15, // DipoleEz: log singularity
  3: 0.4, // DipoleExy: 1/R singularity
};

function renderFieldToRGBA(
  fieldReal: Float64Array,
  fieldImag: Float64Array,
  gridSize: number,
  mode: VisualizationMode,
  params: ComputeParams | null,
): { imageData: Uint8ClampedArray; stats: ImageStats } {
  const n = gridSize * gridSize;

  // Reuse values buffer — never freed, never GC'd
  if (!reusableValues || reusableValues.length !== n) {
    reusableValues = new Float64Array(n);
  }
  const values = reusableValues;

  // Dipole exclusion zone: skip pixels near source when computing min/max
  const isDipole = params != null && params.sourceType >= 2;
  let exclIx = -1;
  let exclIy = -1;
  let exclR2 = 0; // exclusion radius squared, in pixel units
  if (isDipole && params != null) {
    const vs = params.viewSize;
    const pxPerUnit = gridSize / vs;
    exclIx = Math.round((params.dipoleXs / vs + 0.5) * gridSize);
    exclIy = Math.round((0.5 - params.dipoleYs / vs) * gridSize);
    const rPx =
      (DIPOLE_EXCLUSION_RADIUS[params.sourceType] ?? 0.15) * pxPerUnit;
    exclR2 = rPx * rPx;
  }

  let minVal = Infinity;
  let maxVal = -Infinity;

  // Pass 1: extract values and find range
  for (let i = 0; i < n; i++) {
    const re = fieldReal[i];
    const im = fieldImag[i];
    let val: number;
    switch (mode) {
      case "magnitude":
        val = Math.sqrt(re * re + im * im);
        break;
      case "log_magnitude": {
        const mag = Math.sqrt(re * re + im * im);
        val = mag > 0 ? Math.log10(mag) : -20; // floor at -20 for zero
        break;
      }
      case "real":
        val = re;
        break;
      case "imag":
        val = im;
        break;
      case "phase":
        val = Math.atan2(im, re);
        break;
    }
    values[i] = val;
    if (mode !== "phase") {
      // Skip pixels near dipole source for range calculation
      if (isDipole) {
        const px = i % gridSize;
        const py = (i / gridSize) | 0;
        const ddx = px - exclIx;
        const ddy = py - exclIy;
        if (ddx * ddx + ddy * ddy < exclR2) continue;
      }
      if (val < minVal) minVal = val;
      if (val > maxVal) maxVal = val;
    }
  }
  if (mode === "phase") {
    minVal = -Math.PI;
    maxVal = Math.PI;
  }

  // Pass 2: colormap → RGBA (reuse recycled buffer from main thread)
  let rgba: Uint8ClampedArray;
  const requiredBytes = n * 4;
  if (recycledRgbaBuffer && recycledRgbaBuffer.byteLength === requiredBytes) {
    rgba = new Uint8ClampedArray(recycledRgbaBuffer);
    recycledRgbaBuffer = null;
  } else {
    rgba = new Uint8ClampedArray(requiredBytes);
  }

  if (mode === "phase") {
    for (let i = 0; i < n; i++) {
      writePhase(rgba, i * 4, values[i]);
    }
  } else if (mode === "real" || mode === "imag") {
    const absMax = Math.max(Math.abs(minVal), Math.abs(maxVal));
    const invTwoAbsMax = absMax > 0 ? 1 / (2 * absMax) : 0;
    for (let i = 0; i < n; i++) {
      const t = absMax > 0 ? (values[i] + absMax) * invTwoAbsMax : 0.5;
      writeDiverging(rgba, i * 4, t);
    }
  } else {
    const range = maxVal - minVal;
    const invRange = range > 0 ? 1 / range : 0;
    for (let i = 0; i < n; i++) {
      const t = range > 0 ? (values[i] - minVal) * invRange : 0;
      writeViridis(rgba, i * 4, t);
    }
  }

  return { imageData: rgba, stats: { min: minVal, max: maxVal } };
}

// ── Message handler ──────────────────────────────────────────────────

self.onmessage = async (e: MessageEvent<WorkerRequest>) => {
  const msg = e.data;

  if (msg.type === "init") {
    try {
      await init();

      // Pre-allocate field buffers once — never GC'd
      cachedGridSize = get_field_grid_size();
      const n = cachedGridSize * cachedGridSize;
      fieldRealBuf = new Float64Array(n);
      fieldImagBuf = new Float64Array(n);

      const bounds = get_parameter_bounds();
      (self as unknown as Worker).postMessage({
        type: "ready",
        bounds,
      } satisfies WorkerResponse);
    } catch (err: unknown) {
      const message = err instanceof Error ? err.message : String(err);
      (self as unknown as Worker).postMessage({
        type: "error",
        message: `Failed to load WASM: ${message}`,
      } satisfies WorkerResponse);
    }
    return;
  }

  if (msg.type === "compute") {
    recycleBuffer(msg.recycledBuffer);
    try {
      const p = msg.params;
      cachedComputeParams = p;

      // Single WASM call — void return, no JS objects created
      compute_all(
        p.wavelength,
        p.permittivityRe,
        p.permittivityIm,
        p.permeabilityRe,
        p.permeabilityIm,
        p.sourceType,
        p.dipoleXs,
        p.dipoleYs,
        p.dipoleAlpha,
        p.maxOrder,
        p.viewSize,
        fieldRealBuf!,
        fieldImagBuf!,
      );

      // Render to RGBA
      const { imageData, stats } = renderFieldToRGBA(
        fieldRealBuf!,
        fieldImagBuf!,
        cachedGridSize,
        msg.mode,
        p,
      );

      // Order metadata derived from input — no WASM return needed
      const response: WorkerResponse = {
        type: "result",
        scattering: {
          numOrders: 2 * p.maxOrder + 1,
          minOrder: -p.maxOrder,
          maxOrder: p.maxOrder,
        },
        imageData,
        gridSize: cachedGridSize,
        viewSize: p.viewSize,
        stats,
      };

      (self as unknown as Worker).postMessage(response, [imageData.buffer]);
    } catch (err: unknown) {
      const message = err instanceof Error ? err.message : String(err);
      (self as unknown as Worker).postMessage({
        type: "error",
        message: `Computation error: ${message}`,
      } satisfies WorkerResponse);
    }
    return;
  }

  if (msg.type === "recolor") {
    recycleBuffer(msg.recycledBuffer);
    try {
      if (!fieldRealBuf || !fieldImagBuf) {
        return;
      }

      const { imageData, stats } = renderFieldToRGBA(
        fieldRealBuf,
        fieldImagBuf,
        cachedGridSize,
        msg.mode,
        cachedComputeParams,
      );

      const response: WorkerResponse = {
        type: "recolored",
        imageData,
        stats,
      };

      (self as unknown as Worker).postMessage(response, [imageData.buffer]);
    } catch (err: unknown) {
      const message = err instanceof Error ? err.message : String(err);
      (self as unknown as Worker).postMessage({
        type: "error",
        message: `Recolor error: ${message}`,
      } satisfies WorkerResponse);
    }
  }
};
