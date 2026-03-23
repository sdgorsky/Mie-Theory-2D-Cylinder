//! Source definitions for 2D cylinder scattering.
//!
//! Each source type provides:
//! - Coefficients a_n: how the incident field decomposes into cylindrical harmonics
//! - Incident field: the spatial field value at any point (x, y)

use num_complex::Complex64;

/// Supported incident field source types.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Source {
    PlaneWave,
}

/// Compute source coefficients a_n for all orders in the given list.
///
/// Routes to the appropriate source-specific function based on `source`.
pub fn compute_source_coefficients(source: Source, orders: &[i32]) -> Vec<Complex64> {
    match source {
        Source::PlaneWave => compute_planewave_coefficients(orders),
    }
}

/// Plane-wave source: a_n = i^n.
///
/// This is the cylindrical harmonic expansion of exp(i k x),
/// i.e. the Jacobi–Anger identity coefficients.
fn compute_planewave_coefficients(orders: &[i32]) -> Vec<Complex64> {
    orders.iter().map(|&n| Complex64::i().powi(n)).collect()
}

/// Compute the incident field at a single point (x, y).
///
/// Returns the physical field value in real space, not a cylindrical
/// harmonic coefficient.
pub fn compute_incident_field(source: Source, k0: f64, x: f64, y: f64) -> Complex64 {
    match source {
        Source::PlaneWave => compute_planewave_field(k0, x, y),
    }
}

/// Plane-wave incident field: exp(i k0 x).
///
/// The plane wave propagates in the +x direction and is independent of y.
fn compute_planewave_field(k0: f64, x: f64, _y: f64) -> Complex64 {
    let phase = k0 * x;
    Complex64::new(phase.cos(), phase.sin())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_planewave_coefficients() {
        let orders: Vec<i32> = (-3..=3).collect();
        let coeffs = compute_source_coefficients(Source::PlaneWave, &orders);

        assert_eq!(coeffs.len(), 7);

        let tol = 1e-15;
        for (i, &n) in orders.iter().enumerate() {
            let expected = Complex64::i().powi(n);
            let diff = (coeffs[i] - expected).norm();
            assert!(
                diff < tol,
                "Order {n}: got {} expected {}",
                coeffs[i],
                expected
            );
        }
    }
}
