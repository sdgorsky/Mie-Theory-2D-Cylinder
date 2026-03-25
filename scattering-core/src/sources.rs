//! Source definitions for 2D cylinder scattering.
//!
//! Each source type provides:
//! - Coefficients a_n: how the incident field decomposes into cylindrical harmonics
//! - Incident field: the spatial field value at any point (x, y)
//!
//! The wavenumber k passed to these functions is the wavenumber of the medium
//! containing the source. For sources in free space this is k0 (real); for
//! sources inside the cylinder this is k1 = k0*m (generally complex).
//! The caller is responsible for passing the appropriate k.

use crate::bessel::{bessel_j, hankel1};
use num_complex::Complex64;

// ============================================================================
// Core types
// ============================================================================

/// Polarization of the electromagnetic field.
/// Transverse-Magnetic (TM): Electric field parallel to cylinder axis (Ez)
/// Transverse-Electric (TE): Magnetic field parallel to cylinder axis (Hz)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Polarization {
    TM,
    TE,
}

/// Which domain the source resides in relative to the cylinder boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Domain {
    Exterior,
    Interior,
}

/// The physical source type and its source-specific parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SourceKind {
    PlaneWaveTM,
    PlaneWaveTE,
    /// Z-directed electric line source (TM). Isotropic radiation.
    DipoleEz {
        xs: f64,
        ys: f64,
    },
    /// In-plane electric dipole (TE). Orientation angle alpha in radians.
    DipoleExy {
        xs: f64,
        ys: f64,
        alpha: f64,
    },
}

/// Complete source specification: kind + domain.
///
/// Constructed via `Source::new(kind, cylinder_radius)` which automatically
/// determines the domain from the source position.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Source {
    pub kind: SourceKind,
    pub domain: Domain,
}

impl Source {
    /// Construct a source, computing domain from position relative to cylinder.
    pub fn new(kind: SourceKind, cylinder_radius: f64) -> Self {
        let domain = match &kind {
            SourceKind::PlaneWaveTM | SourceKind::PlaneWaveTE => Domain::Exterior,
            SourceKind::DipoleEz { xs, ys } | SourceKind::DipoleExy { xs, ys, .. } => {
                if (xs * xs + ys * ys).sqrt() < cylinder_radius {
                    Domain::Interior
                } else {
                    Domain::Exterior
                }
            }
        };
        Source { kind, domain }
    }

    /// The polarization implied by this source type.
    pub fn polarization(&self) -> Polarization {
        match self.kind {
            SourceKind::PlaneWaveTM | SourceKind::DipoleEz { .. } => Polarization::TM,
            SourceKind::PlaneWaveTE | SourceKind::DipoleExy { .. } => Polarization::TE,
        }
    }

    /// Whether the source is inside the cylinder.
    pub fn is_interior(&self) -> bool {
        self.domain == Domain::Interior
    }
}

// ============================================================================
// Source coefficients a_n (cylindrical harmonic decomposition)
// ============================================================================

/// Compute paired source coefficients (exterior, interior) for all orders.
///
/// For an exterior source, the exterior vector has the coefficients and the
/// interior vector is all zeros (and vice versa for an interior source).
/// `k0` is the free-space wavenumber; `k1` is the interior wavenumber.
pub fn compute_source_coefficient_pair(
    source: Source,
    orders: &[i32],
    k0: Complex64,
    k1: Complex64,
) -> (Vec<Complex64>, Vec<Complex64>) {
    let zeros = vec![Complex64::new(0.0, 0.0); orders.len()];
    match source.domain {
        Domain::Interior => {
            let interior = compute_source_coefficients(source, orders, k0, k1);
            (zeros, interior)
        }
        Domain::Exterior => {
            let exterior = compute_source_coefficients(source, orders, k0, k1);
            (exterior, zeros)
        }
    }
}

/// Compute source coefficients a_n for all orders in the given list.
///
/// `k` is the wavenumber of the medium containing the source.
pub fn compute_source_coefficients(
    source: Source,
    orders: &[i32],
    k0: Complex64,
    kn: Complex64,
) -> Vec<Complex64> {
    match source.kind {
        SourceKind::PlaneWaveTM | SourceKind::PlaneWaveTE => compute_planewave_coefficients(orders),
        SourceKind::DipoleEz { xs, ys } => {
            compute_dipole_ez_coefficients(orders, k0, kn, xs, ys, source.is_interior())
        }
        SourceKind::DipoleExy { xs, ys, alpha } => {
            compute_dipole_exy_coefficients(orders, k0, kn, xs, ys, alpha, source.is_interior())
        }
    }
}

/// Plane-wave source: a_n = i^n.
///
/// Cylindrical harmonic expansion of exp(i k x) (Jacobi-Anger identity).
fn compute_planewave_coefficients(orders: &[i32]) -> Vec<Complex64> {
    orders.iter().map(|&n| Complex64::i().powi(n)).collect()
}

/// DipoleEz source: a_n = (i/4) H_n^(1)(k r0) exp(-i n theta0).
///
/// From Graf's addition theorem applied to the 2D Green's function
/// (i/4) H_0^(1)(k |r - r0|), valid for |r| < r0.
fn compute_dipole_ez_coefficients(
    orders: &[i32],
    k0: Complex64,
    kn: Complex64,
    xs: f64,
    ys: f64,
    is_interior: bool,
) -> Vec<Complex64> {
    let r0 = (xs * xs + ys * ys).sqrt();
    let theta0 = ys.atan2(xs);
    let prefactor = Complex64::i() / 4.0;

    if is_interior {
        orders
            .iter()
            .map(|&n| {
                let hn = bessel_j(n, kn * r0);
                let phase = Complex64::new(0.0, -(n as f64) * theta0).exp();
                prefactor * hn * phase
            })
            .collect()
    } else {
        orders
            .iter()
            .map(|&n| {
                let hn = hankel1(n, k0 * r0);
                let phase = Complex64::new(0.0, -(n as f64) * theta0).exp();
                prefactor * hn * phase
            })
            .collect()
    }
}

/// DipoleExy source:
///   a_n = (k/8) [e^{i alpha}  H_{n+1}^(1)(k r0) e^{-i(n+1) theta0}
///              + e^{-i alpha} H_{n-1}^(1)(k r0) e^{-i(n-1) theta0}]
///
/// Derived by applying the curl of the in-plane vector potential (at
/// orientation alpha) to the addition theorem for the 2D Green's function.
/// Hz = sin(alpha) dG/dx - cos(alpha) dG/dy, using the observation-point
/// gradient identities for J_n(kr) e^{in theta}.
fn compute_dipole_exy_coefficients(
    orders: &[i32],
    k0: Complex64,
    kn: Complex64,
    xs: f64,
    ys: f64,
    alpha: f64,
    is_interior: bool,
) -> Vec<Complex64> {
    let r0 = (xs * xs + ys * ys).sqrt();
    let theta0 = ys.atan2(xs);
    let exp_neg_ia = Complex64::new(alpha.cos(), -alpha.sin()); // e^{-i alpha}
    let exp_pos_ia = Complex64::new(alpha.cos(), alpha.sin()); // e^{i alpha}

    if is_interior {
        orders
            .iter()
            .map(|&n| {
                let jn_plus = bessel_j(n + 1, kn * r0);
                let phase_plus = Complex64::new(0.0, -((n + 1) as f64) * theta0).exp();

                let jn_minus = bessel_j(n - 1, kn * r0);
                let phase_minus = Complex64::new(0.0, -((n - 1) as f64) * theta0).exp();

                kn / 8.0 * (exp_pos_ia * jn_plus * phase_plus + exp_neg_ia * jn_minus * phase_minus)
            })
            .collect()
    } else {
        orders
            .iter()
            .map(|&n| {
                let hn_plus = hankel1(n + 1, k0 * r0);
                let phase_plus = Complex64::new(0.0, -((n + 1) as f64) * theta0).exp();

                let hn_minus = hankel1(n - 1, k0 * r0);
                let phase_minus = Complex64::new(0.0, -((n - 1) as f64) * theta0).exp();

                k0 / 8.0 * (exp_pos_ia * hn_plus * phase_plus + exp_neg_ia * hn_minus * phase_minus)
            })
            .collect()
    }
}

// ============================================================================
// Incident field in real space
// ============================================================================

/// Compute the incident field at a single point (x, y).
///
/// `k` is the wavenumber of the medium containing the source.
/// The caller should only evaluate this in the spatial domain where the
/// source resides (exterior for outside sources, interior for inside sources).
pub fn compute_incident_field(source: Source, k: Complex64, x: f64, y: f64) -> Complex64 {
    match source.kind {
        SourceKind::PlaneWaveTM | SourceKind::PlaneWaveTE => compute_planewave_field(k, x),
        SourceKind::DipoleEz { xs, ys } => compute_dipole_ez_field(k, x, y, xs, ys),
        SourceKind::DipoleExy { xs, ys, alpha } => compute_dipole_exy_field(k, x, y, xs, ys, alpha),
    }
}

/// Plane-wave incident field: exp(i k x).
fn compute_planewave_field(k: Complex64, x: f64) -> Complex64 {
    (Complex64::i() * k * x).exp()
}

/// DipoleEz incident field: (i/4) H_0^(1)(k R).
fn compute_dipole_ez_field(k: Complex64, x: f64, y: f64, xs: f64, ys: f64) -> Complex64 {
    let dx = x - xs;
    let dy = y - ys;
    let dist = (dx * dx + dy * dy).sqrt();
    let prefactor = Complex64::i() / 4.0;
    prefactor * hankel1(0, k * dist)
}

/// DipoleExy incident field: (ik/4) H_1^(1)(k R) sin(phi_R - alpha).
///
/// phi_R = atan2(y - ys, x - xs) is the angle from source to observation point.
fn compute_dipole_exy_field(
    k: Complex64,
    x: f64,
    y: f64,
    xs: f64,
    ys: f64,
    alpha: f64,
) -> Complex64 {
    let dx = x - xs;
    let dy = y - ys;
    let dist = (dx * dx + dy * dy).sqrt();
    let phi_r = dy.atan2(dx);
    let prefactor = Complex64::i() * k / 4.0;
    prefactor * hankel1(1, k * dist) * (phi_r - alpha).sin()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bessel::bessel_j;
    use std::f64::consts::PI;

    /// Helper: build a Source for tests (exterior by default for large r0).
    fn exterior_source(kind: SourceKind) -> Source {
        Source {
            kind,
            domain: Domain::Exterior,
        }
    }

    #[test]
    fn test_planewave_coefficients() {
        let orders: Vec<i32> = (-3..=3).collect();
        let k = Complex64::new(1.0, 0.0);
        let source = exterior_source(SourceKind::PlaneWaveTM);
        let coeffs = compute_source_coefficients(source, &orders, k, k);

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

    /// Verify DipoleEz coefficients reproduce the field via the addition theorem:
    ///   Sigma_n a_n J_n(k r) exp(i n theta)  ~  (i/4) H_0^(1)(k |r - r0|)
    /// at a test point inside the expansion region (r < r0).
    #[test]
    fn test_dipole_ez_addition_theorem() {
        let xs = 2.0;
        let ys = 1.0;
        let source = exterior_source(SourceKind::DipoleEz { xs, ys });
        let k = Complex64::new(2.0 * PI, 0.0);
        let max_order = 25;
        let orders: Vec<i32> = (-max_order..=max_order).collect();

        let coeffs = compute_source_coefficients(source, &orders, k, k);

        // Test point inside expansion region (r < r0)
        let x: f64 = 0.3;
        let y: f64 = -0.2;
        let r = (x * x + y * y).sqrt();
        let theta = y.atan2(x);

        // Reconstruct field from coefficients
        let mut field_sum = Complex64::new(0.0, 0.0);
        for (i, &n) in orders.iter().enumerate() {
            let jn = bessel_j(n, k * r);
            let phase = Complex64::new(0.0, n as f64 * theta).exp();
            field_sum += coeffs[i] * jn * phase;
        }

        // Direct field evaluation
        let field_direct = compute_incident_field(source, k, x, y);

        let err = (field_sum - field_direct).norm();
        assert!(
            err < 1e-10,
            "DipoleEz addition theorem mismatch: err={err:.2e}, sum={field_sum}, direct={field_direct}"
        );
    }

    /// Verify DipoleExy coefficients reproduce the field via the addition theorem.
    #[test]
    fn test_dipole_exy_addition_theorem() {
        let xs = 1.5;
        let ys = -0.8;
        let alpha = PI / 3.0;
        let source = exterior_source(SourceKind::DipoleExy { xs, ys, alpha });
        let k = Complex64::new(2.0 * PI, 0.0);
        let max_order = 25;
        let orders: Vec<i32> = (-max_order..=max_order).collect();

        let coeffs = compute_source_coefficients(source, &orders, k, k);

        // Test point inside expansion region (r < r0)
        let x: f64 = 0.2;
        let y: f64 = 0.1;
        let r = (x * x + y * y).sqrt();
        let theta = y.atan2(x);

        // Reconstruct field from coefficients
        let mut field_sum = Complex64::new(0.0, 0.0);
        for (i, &n) in orders.iter().enumerate() {
            let jn = bessel_j(n, k * r);
            let phase = Complex64::new(0.0, n as f64 * theta).exp();
            field_sum += coeffs[i] * jn * phase;
        }

        // Direct field evaluation
        let field_direct = compute_incident_field(source, k, x, y);

        let err = (field_sum - field_direct).norm();
        assert!(
            err < 1e-10,
            "DipoleExy addition theorem mismatch: err={err:.2e}, sum={field_sum}, direct={field_direct}"
        );
    }

    /// Verify DipoleEz with complex k (interior source scenario).
    #[test]
    fn test_dipole_ez_complex_k() {
        let xs = 0.2;
        let ys = 0.1;
        let source = Source {
            kind: SourceKind::DipoleEz { xs, ys },
            domain: Domain::Interior,
        };
        let k0 = Complex64::new(2.0 * PI, 0.0); // free-space (unused for interior)
        let kn = Complex64::new(4.0 * PI, 0.5); // lossy medium interior wavenumber
        let max_order = 25;
        let orders: Vec<i32> = (-max_order..=max_order).collect();

        let coeffs = compute_source_coefficients(source, &orders, k0, kn);

        // Test point farther from origin than source (r > r0)
        // For interior source expansion: a_n ~ J_n(kn*r0), field ~ J_n(kn*r)
        let x: f64 = 0.4;
        let y: f64 = -0.2;
        let r = (x * x + y * y).sqrt();
        let theta = y.atan2(x);

        // For r > r0, addition theorem: H_0(k|r-r0|) = Σ J_n(kr0) H_n(kr) e^{in(θ-θ0)}
        // So reconstruction uses H_n at observation point
        let mut field_sum = Complex64::new(0.0, 0.0);
        for (i, &n) in orders.iter().enumerate() {
            let hn = hankel1(n, kn * r);
            let phase = Complex64::new(0.0, n as f64 * theta).exp();
            field_sum += coeffs[i] * hn * phase;
        }

        let field_direct = compute_incident_field(source, kn, x, y);

        let err = (field_sum - field_direct).norm();
        assert!(
            err < 1e-8,
            "DipoleEz complex-k mismatch: err={err:.2e}, sum={field_sum}, direct={field_direct}"
        );
    }

    #[test]
    fn test_source_new_computes_domain() {
        let radius = 0.5;

        // Plane waves are always exterior
        let pw = Source::new(SourceKind::PlaneWaveTM, radius);
        assert_eq!(pw.domain, Domain::Exterior);

        // Dipole outside
        let d_out = Source::new(SourceKind::DipoleEz { xs: 1.0, ys: 0.0 }, radius);
        assert_eq!(d_out.domain, Domain::Exterior);

        // Dipole inside
        let d_in = Source::new(SourceKind::DipoleEz { xs: 0.2, ys: 0.1 }, radius);
        assert_eq!(d_in.domain, Domain::Interior);
    }

    #[test]
    fn test_source_polarization() {
        let r = 0.5;
        assert_eq!(
            Source::new(SourceKind::PlaneWaveTM, r).polarization(),
            Polarization::TM
        );
        assert_eq!(
            Source::new(SourceKind::PlaneWaveTE, r).polarization(),
            Polarization::TE
        );
        assert_eq!(
            Source::new(SourceKind::DipoleEz { xs: 1.0, ys: 0.0 }, r).polarization(),
            Polarization::TM
        );
        assert_eq!(
            Source::new(
                SourceKind::DipoleExy {
                    xs: 1.0,
                    ys: 0.0,
                    alpha: 0.0
                },
                r
            )
            .polarization(),
            Polarization::TE
        );
    }
}
