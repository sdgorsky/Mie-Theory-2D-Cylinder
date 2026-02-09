/**
 * Complex number as object (used for frontend input).
 */
export interface ComplexObj {
  re: number;
  im: number;
}

/**
 * Complex number as tuple [re, im] (WASM serialization format from num_complex::Complex64).
 */
export type ComplexTuple = [number, number];

/** Helper to extract real part from ComplexTuple */
export const re = (c: ComplexTuple): number => c[0];

/** Helper to extract imaginary part from ComplexTuple */
export const im = (c: ComplexTuple): number => c[1];

/**
 * Material properties of the cylinder.
 * Both permittivity and permeability are complex to account for losses.
 */
export interface MaterialProperties {
  /** Complex relative electric permittivity (ε_r = ε' + iε'') */
  permittivity: ComplexObj;
  /** Complex relative magnetic permeability (μ_r = μ' + iμ'') */
  permeability: ComplexObj;
}

/**
 * Polarization of the incident electromagnetic wave.
 */
export type Polarization = "TM" | "TE";

/**
 * Complete specification for 2D electromagnetic scattering simulation.
 * The cylinder has a fixed diameter of 1 (unitless) and is centered at the origin.
 */
export interface ScatteringParams {
  /** Wavelength (unitless, relative to cylinder diameter = 1) */
  wavelength: number;
  /** Material properties of the cylinder */
  material: MaterialProperties;
  /** Polarization of incident wave (TM or TE) */
  polarization: Polarization;
  /** Maximum Bessel order for the expansion (computes -N to +N) */
  maxOrder: number;
}

/**
 * Result of the scattering calculation from WASM.
 * Complex numbers are serialized as [re, im] tuples by serde.
 */
export interface ScatteringResult {
  /** Incident wave coefficients a_n */
  incident_coefficients: ComplexTuple[];
  /** Scattering coefficients b_n (exterior field) */
  scattering_coefficients: ComplexTuple[];
  /** Internal coefficients c_n (interior field) */
  internal_coefficients: ComplexTuple[];
  /** Orders included: from -max_order to +max_order */
  orders: number[];
}

/**
 * Result of field computation from WASM.
 * Contains the complex electric field values on a 1024x1024 grid.
 */
export interface FieldResult {
  /** Real part of complex field values (row-major, 1024x1024) */
  field_real: Float64Array | number[];
  /** Imaginary part of complex field values */
  field_imag: Float64Array | number[];
  /** Grid size (1024) */
  grid_size: number;
  /** Physical extent of the grid (5.0 means -2.5 to +2.5) */
  view_size: number;
  /** Minimum x coordinate */
  x_min: number;
  /** Maximum x coordinate */
  x_max: number;
  /** Minimum y coordinate */
  y_min: number;
  /** Maximum y coordinate */
  y_max: number;
}

/** Grid resolution for field computation */
export const FIELD_GRID_SIZE = 1024;

/** View size in units of cylinder diameter */
export const FIELD_VIEW_SIZE = 5.0;

/**
 * Creates default scattering parameters.
 */
export function createDefaultParams(): ScatteringParams {
  return {
    wavelength: 1.0,
    material: {
      permittivity: { re: 4.0, im: 0.0 },
      permeability: { re: 1.0, im: 0.0 },
    },
    polarization: "TM",
    maxOrder: 10,
  };
}

/**
 * Calculate the size parameter x = π * d / λ
 */
export function calculateSizeParameter(wavelength: number): number {
  const diameter = 1.0;
  return (Math.PI * diameter) / wavelength;
}

/**
 * Calculate the complex refractive index n = sqrt(εr * μr)
 */
export function calculateRefractiveIndex(
  material: MaterialProperties,
): ComplexObj {
  const eps = material.permittivity;
  const mu = material.permeability;

  // Complex multiplication: (a + bi)(c + di) = (ac - bd) + (ad + bc)i
  const prodReal = eps.re * mu.re - eps.im * mu.im;
  const prodImag = eps.re * mu.im + eps.im * mu.re;

  // Complex square root
  const magnitude = Math.sqrt(prodReal * prodReal + prodImag * prodImag);
  const phase = Math.atan2(prodImag, prodReal);

  return {
    re: Math.sqrt(magnitude) * Math.cos(phase / 2),
    im: Math.sqrt(magnitude) * Math.sin(phase / 2),
  };
}

/**
 * Format a complex number for display.
 */
export function formatComplex(c: ComplexObj): string {
  const realStr = c.re.toFixed(3);
  const imagAbs = Math.abs(c.im).toFixed(3);
  const sign = c.im >= 0 ? "+" : "-";
  return `${realStr} ${sign} ${imagAbs}i`;
}
