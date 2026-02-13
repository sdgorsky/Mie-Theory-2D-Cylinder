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
