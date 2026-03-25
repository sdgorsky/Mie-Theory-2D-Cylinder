/**
 * Complex number as object (used for frontend input).
 */
export interface ComplexObj {
  re: number;
  im: number;
}

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
 * Polarization of the electromagnetic field.
 */
export type Polarization = "TM" | "TE";

/**
 * Combined source + polarization selection.
 * Each option implies a specific polarization:
 *   planewave_tm / dipole_ez → TM (Ez is the scalar field)
 *   planewave_te / dipole_exy → TE (Hz is the scalar field)
 */
export type SourceType =
  | "planewave_tm"
  | "planewave_te"
  | "dipole_ez"
  | "dipole_exy";

/** Dipole source position and orientation. */
export interface DipoleParams {
  xs: number;
  ys: number;
  /** Orientation angle in radians (only used for dipole_exy). */
  alpha: number;
}

/** Derive polarization from source type. */
export function getPolarization(sourceType: SourceType): Polarization {
  switch (sourceType) {
    case "planewave_tm":
    case "dipole_ez":
      return "TM";
    case "planewave_te":
    case "dipole_exy":
      return "TE";
  }
}

/** Whether the source type is a dipole. */
export function isDipoleSource(sourceType: SourceType): boolean {
  return sourceType === "dipole_ez" || sourceType === "dipole_exy";
}

/**
 * Complete specification for 2D electromagnetic scattering simulation.
 * The cylinder has a fixed diameter of 1 (unitless) and is centered at the origin.
 */
export interface ScatteringParams {
  /** Wavelength (unitless, relative to cylinder diameter = 1) */
  wavelength: number;
  /** Material properties of the cylinder */
  material: MaterialProperties;
  /** Source type (determines both source and polarization) */
  sourceType: SourceType;
  /** Dipole position and orientation (ignored for plane wave sources) */
  dipole: DipoleParams;
  /** Maximum Bessel order for the expansion (computes -N to +N) */
  maxOrder: number;
}

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
    sourceType: "planewave_tm",
    dipole: { xs: 1.5, ys: 0.0, alpha: 0.0 },
    maxOrder: 21,
  };
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
export function formatComplex(c: ComplexObj, decimals = 3): string {
  const realStr = c.re.toFixed(decimals);
  const imagAbs = Math.abs(c.im).toFixed(decimals);
  const sign = c.im >= 0 ? "+" : "-";
  return `${realStr} ${sign} ${imagAbs}i`;
}
