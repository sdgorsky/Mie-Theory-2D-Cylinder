//! Electromagnetic scattering coefficients for 2D cylinder.
//!
//! Computes the scattering (bn) and internal (cn) coefficients for
//! plane wave scattering from an infinite dielectric/magnetic cylinder.

use crate::bessel::{bessel_j, bessel_j_derivative, hankel1, hankel1_derivative};
use num_complex::Complex64;
use std::f64::consts::PI;

/// Polarization of the incident electromagnetic wave.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Polarization {
    TM, // Transverse-Magnetic: Electric field parallel to cylinder axis (Ez)
    TE, // Transverse-Electric: Magnetic field parallel to cylinder axis (Hz)
}

/// Complex material properties.
#[derive(Debug, Clone, Copy)]
pub struct Material {
    pub permittivity_real: f64,
    pub permittivity_imag: f64,
    pub permeability_real: f64,
    pub permeability_imag: f64,
}

impl Material {
    pub fn permittivity(&self) -> Complex64 {
        Complex64::new(self.permittivity_real, self.permittivity_imag)
    }

    pub fn permeability(&self) -> Complex64 {
        Complex64::new(self.permeability_real, self.permeability_imag)
    }

    /// Complex refractive index: n = sqrt(εr * μr)
    pub fn refractive_index(&self) -> Complex64 {
        (self.permittivity() * self.permeability()).sqrt()
    }
}

/// Input parameters for scattering calculation.
#[derive(Debug, Clone)]
pub struct ScatteringParams {
    /// Wavelength (unitless, diameter is fixed at 1)
    pub wavelength: f64,
    /// Cylinder material properties
    pub material: Material,
    /// Polarization (TM or TE)
    pub polarization: Polarization,
    /// Maximum Bessel order to include in the expansion
    pub max_order: i32,
}

/// Result of the scattering calculation.
#[derive(Debug, Clone)]
pub struct ScatteringResult {
    /// Incident wave coefficients a_n (n = -max_order to +max_order)
    pub incident_coefficients: Vec<Complex64>,
    /// Scattering coefficients b_n (exterior field)
    pub scattering_coefficients: Vec<Complex64>,
    /// Internal coefficients c_n (interior field)
    pub internal_coefficients: Vec<Complex64>,
    /// Orders included: from -max_order to +max_order
    pub orders: Vec<i32>,
}

/// Cylinder radius (diameter = 1, so radius = 0.5)
pub const RADIUS: f64 = 0.5;

// Parameter bounds — single source of truth for UI sliders and tests
pub const WAVELENGTH_MIN: f64 = 0.3;
pub const WAVELENGTH_MAX: f64 = 5.0;
pub const PERMITTIVITY_RE_MIN: f64 = -5.0;
pub const PERMITTIVITY_RE_MAX: f64 = 5.0;
pub const PERMITTIVITY_IM_MIN: f64 = -5.0;
pub const PERMITTIVITY_IM_MAX: f64 = 5.0;
pub const PERMEABILITY_RE_MIN: f64 = -5.0;
pub const PERMEABILITY_RE_MAX: f64 = 5.0;
pub const PERMEABILITY_IM_MIN: f64 = -5.0;
pub const PERMEABILITY_IM_MAX: f64 = 5.0;
pub const MAX_ORDER_MIN: i32 = 1;
pub const MAX_ORDER_MAX: i32 = 50;

/// Calculate scattering coefficients for a 2D cylinder.
///
/// The cylinder has diameter = 1 (unitless), centered at the origin.
/// Returns coefficients for orders n = -max_order to +max_order.
pub fn calculate_scattering(params: &ScatteringParams) -> ScatteringResult {
    // Size parameters: k*a where a is the cylinder radius
    let k0 = 2.0 * PI / params.wavelength; // Free-space wavenumber
    let kor = Complex64::new(k0 * RADIUS, 0.0); // k0 * a
    let knr = kor * params.material.refractive_index(); // k1 * a = k0 * m * a

    let mut an_vec = Vec::new();
    let mut bn_vec = Vec::new();
    let mut cn_vec = Vec::new();
    let mut l_vec = Vec::new();

    // Calculate coefficients for l = -max_order to +max_order
    for l in -params.max_order..=params.max_order {
        let an = Complex64::i().powi(l);
        let (bn, cn) = calculate_coefficients_for_order(params, l, kor, knr, params.polarization);

        l_vec.push(l);
        an_vec.push(an);
        bn_vec.push(bn);
        cn_vec.push(cn);
    }

    ScatteringResult {
        incident_coefficients: an_vec,
        scattering_coefficients: bn_vec,
        internal_coefficients: cn_vec,
        orders: l_vec,
    }
}

/// Calculate b_n and c_n for a single order n.
/// Note: a_n = i^n is the incident wave coefficient for a plane wave.
fn calculate_coefficients_for_order(
    params: &ScatteringParams,
    l: i32,         // angular order
    kor: Complex64, // free-space size parameterf
    knr: Complex64, // internal size parameter
    polarization: Polarization,
) -> (Complex64, Complex64) {
    let an = Complex64::i().powi(l);

    let jlo = bessel_j(l, kor);
    let jlo_prime = bessel_j_derivative(l, kor);
    let hlo = hankel1(l, kor);
    let hlo_prime = hankel1_derivative(l, kor);

    let jln = bessel_j(l, knr);
    let jln_prime = bessel_j_derivative(l, knr);

    let psi: Complex64 = match polarization {
        Polarization::TM => 1.0 / params.material.permittivity(),
        Polarization::TE => 1.0 / params.material.permeability(),
    };

    let gamma = psi * (knr / kor) * (jln_prime / jln);

    let numerator = jlo_prime - gamma * jlo;
    let denominator = hlo_prime - gamma * hlo;
    let bn = -an * numerator / denominator;
    let cn = (an * jlo + bn * hlo) / jln;
    (bn, cn)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: Complex64, b: Complex64, tol: f64) -> bool {
        (a - b).norm() < tol
    }

    #[test]
    fn test_dielectric_cylinder() {
        let params = ScatteringParams {
            wavelength: 1.0,
            material: Material {
                permittivity_real: 4.0,
                permittivity_imag: 0.0,
                permeability_real: 1.0,
                permeability_imag: 0.0,
            },
            polarization: Polarization::TM,
            max_order: 5,
        };

        let result = calculate_scattering(&params);

        // Should have 11 coefficients: -5 to +5
        assert_eq!(result.scattering_coefficients.len(), 11);
        assert_eq!(result.orders.len(), 11);
        assert_eq!(result.orders[0], -5);
        assert_eq!(result.orders[10], 5);

        // Expected values computed with complex-bessel-rs crate
        // Size parameter: x = k*a = (2π/λ) * 0.5 = π for λ=1, radius=0.5
        let expected = vec![
            Complex64::new(0.03729800378834639, 0.001393081763394997), // b_{-5}
            Complex64::new(-0.21924538868704604, 0.4137352392853578),  // b_{-4}
            Complex64::new(0.4959048849072327, -0.43613807765855905),  // b_{-3}
            Complex64::new(0.07495853048067841, 0.2633244181401634),   // b_{-2}
            Complex64::new(-0.4110146488234517, 0.2152774008749309),   // b_{-1}
            Complex64::new(-0.06642015681023666, 0.24901509909951297), // b_0
            Complex64::new(0.4110146488234517, -0.2152774008749309),   // b_1
            Complex64::new(0.07495853048067841, 0.2633244181401634),   // b_2
            Complex64::new(-0.4959048849072327, 0.43613807765855905),  // b_3
            Complex64::new(-0.21924538868704604, 0.4137352392853578),  // b_4
            Complex64::new(-0.03729800378834639, -0.001393081763394997), // b_5
        ];

        for (i, (computed, exp)) in result
            .scattering_coefficients
            .iter()
            .zip(expected.iter())
            .enumerate()
        {
            assert!(
                approx_eq(*computed, *exp, 1e-10),
                "Coefficient {} mismatch: computed {} vs expected {}",
                result.orders[i],
                computed,
                exp
            );
        }

        // Verify symmetry: b_{-n} = (-1)^n * conj(b_n) for real materials
        // For TM mode with real ε, μ: b_{-n} should relate to b_n by symmetry
        let n_coeffs = result.scattering_coefficients.len();
        let mid = n_coeffs / 2; // index of b_0
        for i in 1..=5 {
            let b_neg = result.scattering_coefficients[mid - i]; // b_{-i}
            let b_pos = result.scattering_coefficients[mid + i]; // b_i
                                                                 // For integer orders: b_{-n} = (-1)^n * b_n
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            assert!(
                approx_eq(b_neg, sign * b_pos, 1e-10),
                "Symmetry check failed for order {}: b_{} = {} vs (-1)^{} * b_{} = {}",
                i,
                -(i as i32),
                b_neg,
                i,
                i,
                sign * b_pos
            );
        }
    }

    #[test]
    fn test_vacuum_cylinder() {
        // Cylinder with same properties as vacuum (εr=1, μr=1)
        // should have zero scattering
        let params = ScatteringParams {
            wavelength: 1.0,
            material: Material {
                permittivity_real: 1.0,
                permittivity_imag: 0.0,
                permeability_real: 1.0,
                permeability_imag: 0.0,
            },
            polarization: Polarization::TM,
            max_order: 3,
        };

        let result = calculate_scattering(&params);

        // Scattering coefficients should be nearly zero
        for coeff in &result.scattering_coefficients {
            let magnitude = coeff.norm();
            assert!(
                magnitude < 1e-10,
                "Expected zero scattering for vacuum, got {}",
                magnitude
            );
        }
    }
}
