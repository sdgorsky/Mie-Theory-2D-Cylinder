//! Electromagnetic scattering coefficients for 2D cylinder.
//!
//! Computes the scattering (bn) and internal (cn) coefficients for
//! plane wave scattering from an infinite dielectric/magnetic cylinder.

use crate::bessel::{bessel_j, bessel_j_derivative, hankel1, hankel1_derivative};
use crate::sources::{compute_source_coefficient_pair, Polarization, Source};
use num_complex::Complex64;
use std::f64::consts::PI;

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
    /// Maximum Bessel order to include in the expansion
    pub max_order: i32,
    /// Incident field source (carries polarization and domain)
    pub source: Source,
}

/// Result of the scattering calculation.
#[derive(Debug, Clone)]
pub struct ScatteringResult {
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

    let mut bn_vec = Vec::new();
    let mut cn_vec = Vec::new();
    let mut l_vec = Vec::new();

    let orders: Vec<i32> = (-params.max_order..=params.max_order).collect();
    let k0c = Complex64::new(k0, 0.0);
    let k1 = k0c * params.material.refractive_index();
    let (a_n_ext, a_n_int) = compute_source_coefficient_pair(params.source, &orders, k0c, k1);

    let polarization = params.source.polarization();

    // Calculate coefficients for l = -max_order to +max_order
    for (idx, &l) in orders.iter().enumerate() {
        let (bn, cn) = calculate_coefficients_for_order(
            params,
            l,
            kor,
            knr,
            polarization,
            a_n_ext[idx],
            a_n_int[idx],
        );

        l_vec.push(l);
        bn_vec.push(bn);
        cn_vec.push(cn);
    }

    ScatteringResult {
        scattering_coefficients: bn_vec,
        internal_coefficients: cn_vec,
        orders: l_vec,
    }
}

/// Calculate b_n and c_n for a single order n.
fn calculate_coefficients_for_order(
    params: &ScatteringParams,
    l: i32,         // angular order
    kor: Complex64, // free-space size parameter
    knr: Complex64, // internal size parameter
    polarization: Polarization,
    an_exterior: Complex64, // exterior source coefficients
    an_interior: Complex64, // interior source coefficients
) -> (Complex64, Complex64) {
    let jlo = bessel_j(l, kor);
    let jlo_prime = bessel_j_derivative(l, kor);
    let jln = bessel_j(l, knr);
    let jln_prime = bessel_j_derivative(l, knr);

    let hlo = hankel1(l, kor);
    let hlo_prime = hankel1_derivative(l, kor);
    let hln = hankel1(l, knr);
    let hln_prime = hankel1_derivative(l, knr);

    let psi: Complex64 = match polarization {
        Polarization::TM => 1.0 / params.material.permittivity(),
        Polarization::TE => 1.0 / params.material.permeability(),
    };

    let factor = psi * knr / kor;
    let gamma = factor * jln_prime / jln;

    let sn = -(jlo_prime - gamma * jlo) / (hlo_prime - gamma * hlo);
    let tn = (factor * hln_prime - gamma * hln) / (hlo_prime - gamma * hlo);

    let bn = an_exterior * sn + an_interior * tn;
    let cn = (an_exterior * jlo - an_interior * hln + bn * hlo) / jln;
    (bn, cn)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sources::{Source, SourceKind};

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
            max_order: 5,
            source: Source::new(SourceKind::PlaneWaveTM, RADIUS),
        };

        let result = calculate_scattering(&params);

        // Should have 11 coefficients: -5 to +5
        assert_eq!(result.scattering_coefficients.len(), 11);
        assert_eq!(result.orders.len(), 11);
        assert_eq!(result.orders[0], -5);
        assert_eq!(result.orders[10], 5);

        // Expected values computed with complex-bessel-rs crate
        // Size parameter: x = k*a = (2pi/lambda) * 0.5 = pi for lambda=1, radius=0.5
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
        let n_coeffs = result.scattering_coefficients.len();
        let mid = n_coeffs / 2; // index of b_0
        for i in 1..=5 {
            let b_neg = result.scattering_coefficients[mid - i]; // b_{-i}
            let b_pos = result.scattering_coefficients[mid + i]; // b_i
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
            max_order: 3,
            source: Source::new(SourceKind::PlaneWaveTM, RADIUS),
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

    /// Shared material and wavelength for dipole coefficient tests,
    /// matching the boundary condition test in field.rs.
    // Real dielectric+magnetic material — avoids complex Bessel arguments
    // where our hankel1 implementation has a known bug (NaN for |Im(z)| > ~0.01).
    fn dipole_test_material() -> Material {
        Material {
            permittivity_real: 2.3,
            permittivity_imag: 0.0,
            permeability_real: 1.4,
            permeability_imag: 0.0,
        }
    }
    const DIPOLE_TEST_WAVELENGTH: f64 = 1.2;

    /// Helper: run calculate_scattering for a dipole config and check b_n
    /// against hardcoded reference values for orders -3..=3.
    fn check_dipole_coefficients(
        label: &str,
        kind: SourceKind,
        max_order: i32,
        expected_b: &[Complex64],
    ) {
        let source = Source::new(kind, RADIUS);
        let params = ScatteringParams {
            wavelength: DIPOLE_TEST_WAVELENGTH,
            material: dipole_test_material(),
            max_order,
            source,
        };
        let result = calculate_scattering(&params);
        let expected_n = (2 * max_order + 1) as usize;

        assert_eq!(
            result.scattering_coefficients.len(),
            expected_n,
            "{label}: count"
        );
        assert_eq!(result.orders[0], -max_order, "{label}: min order");

        // Compare the central 7 coefficients (orders -3..=3)
        let mid = result.orders.iter().position(|&o| o == 0).unwrap();
        for (j, &n) in (-3_i32..=3).collect::<Vec<_>>().iter().enumerate() {
            let idx = (mid as i32 + n) as usize;
            let computed = result.scattering_coefficients[idx];
            let diff = (computed - expected_b[j]).norm();
            let scale = expected_b[j].norm().max(1e-15);
            assert!(
                diff / scale < 1e-8,
                "{label} order {n}: computed {computed} vs expected {}, rel err = {:.2e}",
                expected_b[j],
                diff / scale
            );
        }
    }

    #[test]
    fn test_dipole_ez_exterior_coefficients() {
        #[rustfmt::skip]
        let expected_b = vec![
            Complex64::new(-1.92688951246365210e-2, -3.25985959742882761e-2),
            Complex64::new(-4.42897979981103432e-2, -4.45036788702432573e-2),
            Complex64::new(-5.56114045090237291e-2,  6.26494057845542948e-3),
            Complex64::new(-1.06392947571174370e-2,  5.84459410048261893e-2),
            Complex64::new( 3.64439025807683711e-2,  4.24702219984933321e-2),
            Complex64::new( 5.76779109585631791e-2, -2.48077049957309177e-2),
            Complex64::new(-2.59238904099264658e-2, -2.76027296627208117e-2),
        ];
        check_dipole_coefficients(
            "DipoleEz exterior",
            SourceKind::DipoleEz { xs: 1.7, ys: -0.9 },
            5,
            &expected_b,
        );
    }

    #[test]
    fn test_dipole_ez_interior_coefficients() {
        #[rustfmt::skip]
        let expected_b = vec![
            Complex64::new(-1.23082570259623082e-2, -5.74852023304760380e-3),
            Complex64::new(-2.67352016511678697e-2,  4.47690720080264962e-2),
            Complex64::new( 8.67203082943883430e-2, -9.73677090469979301e-3),
            Complex64::new(-6.38803503069468731e-2, -2.19496278598972083e-2),
            Complex64::new(-4.23417532560260712e-2, -7.63046034622431729e-2),
            Complex64::new(-1.29633115057645049e-2, -5.05073595686112362e-2),
            Complex64::new(-9.23419585565310452e-3, -9.96334797075497207e-3),
        ];
        check_dipole_coefficients(
            "DipoleEz interior",
            SourceKind::DipoleEz { xs: 0.15, ys: -0.1 },
            5,
            &expected_b,
        );
    }

    #[test]
    fn test_dipole_exy_exterior_coefficients() {
        #[rustfmt::skip]
        let expected_b = vec![
            Complex64::new( 1.61646549031968234e-1,  1.30787698709159017e-1),
            Complex64::new(-2.79497861090538313e-1, -6.60633996553582209e-2),
            Complex64::new( 1.23353760465490792e-1, -2.91383366061061577e-1),
            Complex64::new( 1.28452640464221751e-1,  2.54211804300940647e-1),
            Complex64::new(-3.37956554335932080e-1,  1.09528660580023848e-1),
            Complex64::new( 1.32928408144887072e-1, -3.39486060855426042e-1),
            Complex64::new( 3.01751684266101439e-1,  3.03869495675694601e-2),
        ];
        check_dipole_coefficients(
            "DipoleExy exterior",
            SourceKind::DipoleExy {
                xs: -1.3,
                ys: 0.6,
                alpha: 0.73,
            },
            5,
            &expected_b,
        );
    }

    /// Canary for the hankel1 complex-argument NaN bug.
    /// Uses the exact reproducer from project_hankel1_complex_bug.md:
    ///   eps_r = 2.3 - 0.7i, mu_r = 1.4 + 0.3i, lambda = 1.2
    /// which gives knr ≈ 4.853 - 0.205i, triggering NaN in hankel1.
    #[test]
    fn test_lossy_material_no_nan() {
        let params = ScatteringParams {
            wavelength: 1.2,
            material: Material {
                permittivity_real: 2.3,
                permittivity_imag: -0.7,
                permeability_real: 1.4,
                permeability_imag: 0.3,
            },
            max_order: 5,
            source: Source::new(SourceKind::PlaneWaveTM, RADIUS),
        };

        let result = calculate_scattering(&params);

        // Every coefficient must be finite (not NaN or Inf)
        for (i, bn) in result.scattering_coefficients.iter().enumerate() {
            assert!(
                bn.re.is_finite() && bn.im.is_finite(),
                "b_n[{}] (order {}) is not finite: {}",
                i,
                result.orders[i],
                bn
            );
        }
        for (i, cn) in result.internal_coefficients.iter().enumerate() {
            assert!(
                cn.re.is_finite() && cn.im.is_finite(),
                "c_n[{}] (order {}) is not finite: {}",
                i,
                result.orders[i],
                cn
            );
        }

        // Scattering must actually occur — not all zeros
        let max_bn = result
            .scattering_coefficients
            .iter()
            .map(|c| c.norm())
            .fold(0.0_f64, f64::max);
        assert!(
            max_bn > 1e-10,
            "All scattering coefficients are ~zero; scattering didn't happen (max |b_n| = {:.2e})",
            max_bn
        );
    }

    #[test]
    fn test_dipole_exy_interior_coefficients() {
        #[rustfmt::skip]
        let expected_b = vec![
            Complex64::new(-4.06447172257134615e-1, -4.67865479192744738e-1),
            Complex64::new( 6.16161393188024986e-1, -1.19998513938925327e-1),
            Complex64::new( 3.16483538781634599e-1,  3.68734517640899540e-1),
            Complex64::new(-6.36782914696539526e-2, -4.59490211118910863e-2),
            Complex64::new(-4.64851760952698201e-1,  1.41562054092150003e-1),
            Complex64::new( 4.07405336744003199e-1,  4.77572400160921839e-1),
            Complex64::new(-4.02377582913878207e-1, -4.71370015194499448e-1),
        ];
        check_dipole_coefficients(
            "DipoleExy interior",
            SourceKind::DipoleExy {
                xs: -0.12,
                ys: 0.2,
                alpha: -1.1,
            },
            5,
            &expected_b,
        );
    }
}
