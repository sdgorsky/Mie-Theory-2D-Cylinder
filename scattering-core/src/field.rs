//! Electric field computation for 2D cylinder scattering.
//!
//! Computes the total electric field (incident + scattered) on a 2D grid
//! using the cylindrical wave expansion coefficients.
//!
//! Optimizations:
//! - Cubic spline interpolation for Bessel/Hankel functions
//! - Y-axis symmetry (only compute top half, mirror to bottom)

use crate::bessel::{bessel_j, hankel1};
use crate::scattering::RADIUS;
use num_complex::Complex64;
use std::f64::consts::PI;

/// Grid resolution for field computation
pub const GRID_SIZE: usize = 256;

/// Default view size in units of cylinder diameter (5D x 5D)
pub const DEFAULT_VIEW_SIZE: f64 = 5.0;

/// Number of points in radial spline grids
const SPLINE_POINTS: usize = 128;

// ============================================================================
// Cubic Spline Interpolation for Bessel Functions
// ============================================================================

/// Cubic spline for interpolating complex-valued Bessel/Hankel functions.
/// Stores spline coefficients for both real and imaginary parts.
struct BesselSpline {
    r_min: f64,
    dr: f64,
    /// Spline coefficients [a, b, c, d] for each interval (real part)
    coeffs_re: Vec<[f64; 4]>,
    /// Spline coefficients [a, b, c, d] for each interval (imaginary part)
    coeffs_im: Vec<[f64; 4]>,
}

impl BesselSpline {
    /// Create spline for Hankel function H_n^(1)(k0*r) over [r_min, r_max]
    fn new_hankel(order: i32, k0: f64, r_min: f64, r_max: f64, n_points: usize) -> Self {
        let dr = (r_max - r_min) / (n_points - 1) as f64;
        let mut values_re = Vec::with_capacity(n_points);
        let mut values_im = Vec::with_capacity(n_points);

        for i in 0..n_points {
            let r = r_min + i as f64 * dr;
            let hn = hankel1(order, Complex64::new(k0 * r, 0.0));
            values_re.push(hn.re);
            values_im.push(hn.im);
        }

        Self::from_values(r_min, dr, &values_re, &values_im)
    }

    /// Create spline for Bessel J_n(k1*r) over [r_min, r_max] where k1 may be complex
    fn new_bessel_j(order: i32, k1: Complex64, r_min: f64, r_max: f64, n_points: usize) -> Self {
        let dr = (r_max - r_min) / (n_points - 1) as f64;
        let mut values_re = Vec::with_capacity(n_points);
        let mut values_im = Vec::with_capacity(n_points);

        for i in 0..n_points {
            let r = r_min + i as f64 * dr;
            let jn = bessel_j(order, k1 * r);
            values_re.push(jn.re);
            values_im.push(jn.im);
        }

        Self::from_values(r_min, dr, &values_re, &values_im)
    }

    /// Build spline from precomputed values
    fn from_values(r_min: f64, dr: f64, values_re: &[f64], values_im: &[f64]) -> Self {
        let coeffs_re = Self::compute_spline_coefficients(values_re, dr);
        let coeffs_im = Self::compute_spline_coefficients(values_im, dr);
        BesselSpline {
            r_min,
            dr,
            coeffs_re,
            coeffs_im,
        }
    }

    /// Compute natural cubic spline coefficients using Thomas algorithm
    fn compute_spline_coefficients(y: &[f64], h: f64) -> Vec<[f64; 4]> {
        let n = y.len();
        if n < 2 {
            return vec![];
        }

        let mut m = vec![0.0; n]; // Second derivatives

        if n > 2 {
            let mut rhs = vec![0.0; n - 2];
            for i in 1..n - 1 {
                rhs[i - 1] = 6.0 / (h * h) * (y[i + 1] - 2.0 * y[i] + y[i - 1]);
            }

            // Thomas algorithm for tridiagonal [1, 4, 1] system
            let mut c_prime = vec![0.0; n - 2];
            let mut d_prime = vec![0.0; n - 2];

            c_prime[0] = 1.0 / 4.0;
            d_prime[0] = rhs[0] / 4.0;

            for i in 1..n - 2 {
                let denom = 4.0 - c_prime[i - 1];
                c_prime[i] = 1.0 / denom;
                d_prime[i] = (rhs[i] - d_prime[i - 1]) / denom;
            }

            m[n - 2] = d_prime[n - 3];
            for i in (1..n - 2).rev() {
                m[i] = d_prime[i - 1] - c_prime[i - 1] * m[i + 1];
            }
        }

        let mut coeffs = Vec::with_capacity(n - 1);
        for i in 0..n - 1 {
            let a = y[i];
            let b = (y[i + 1] - y[i]) / h - h * (2.0 * m[i] + m[i + 1]) / 6.0;
            let c = m[i] / 2.0;
            let d = (m[i + 1] - m[i]) / (6.0 * h);
            coeffs.push([a, b, c, d]);
        }
        coeffs
    }

    /// Evaluate spline at radius r, returning complex value
    #[inline]
    fn eval(&self, r: f64) -> Complex64 {
        let t = (r - self.r_min) / self.dr;
        let i = (t as usize).min(self.coeffs_re.len().saturating_sub(1));
        let x = r - (self.r_min + i as f64 * self.dr);

        let [a_re, b_re, c_re, d_re] = self.coeffs_re[i];
        let [a_im, b_im, c_im, d_im] = self.coeffs_im[i];

        let re = a_re + x * (b_re + x * (c_re + x * d_re));
        let im = a_im + x * (b_im + x * (c_im + x * d_im));

        Complex64::new(re, im)
    }
}

// ============================================================================
// Field Parameters and Results
// ============================================================================

/// Input parameters for field computation.
#[derive(Debug, Clone)]
pub struct FieldParams {
    /// Wavelength (relative to diameter = 1)
    pub wavelength: f64,
    /// Complex relative permittivity
    pub permittivity_real: f64,
    pub permittivity_imag: f64,
    /// Complex relative permeability
    pub permeability_real: f64,
    pub permeability_imag: f64,
    /// Incident wave coefficients (a_n) - real and imaginary parts interleaved
    pub incident_coeffs_real: Vec<f64>,
    pub incident_coeffs_imag: Vec<f64>,
    /// Scattering coefficients (b_n) - real and imaginary parts
    pub scattering_coeffs_real: Vec<f64>,
    pub scattering_coeffs_imag: Vec<f64>,
    /// Internal coefficients (c_n) - real and imaginary parts
    pub internal_coeffs_real: Vec<f64>,
    pub internal_coeffs_imag: Vec<f64>,
    /// Orders corresponding to coefficients (-N to +N)
    pub orders: Vec<i32>,
    /// View size in cylinder diameters (physical extent of the grid)
    pub view_size: f64,
}

/// Result of field computation.
#[derive(Debug, Clone)]
pub struct FieldResult {
    pub field_real: Vec<f64>, // Real part of complex field values (row-major, GRID_SIZE x GRID_SIZE)
    pub field_imag: Vec<f64>, // Imaginary part of complex field values
    pub grid_size: usize,     // Number of grid points on each edge
    pub view_size: f64,       // Physical dimension of the grid (5.0 means -2.5 to +2.5)
    // Grid extrema
    pub x_min: f64,
    pub x_max: f64,
    pub y_min: f64,
    pub y_max: f64,
}

/// Compute the electric field on a 128x128 grid.
///
/// The grid covers a 5D x 5D area centered on the cylinder (D = diameter = 1).
/// So the grid spans from -2.5 to +2.5 in both x and y.
///
/// Uses optimized spline interpolation for Bessel/Hankel functions and
/// y-axis symmetry to compute only half the grid points.
pub fn compute_field(params: &FieldParams) -> FieldResult {
    let n_coeffs = params.orders.len();

    // Convert coefficient arrays to complex vectors
    let mut b_n: Vec<Complex64> = Vec::with_capacity(n_coeffs);
    let mut c_n: Vec<Complex64> = Vec::with_capacity(n_coeffs);

    for i in 0..n_coeffs {
        b_n.push(Complex64::new(
            params.scattering_coeffs_real[i],
            params.scattering_coeffs_imag[i],
        ));
        c_n.push(Complex64::new(
            params.internal_coeffs_real[i],
            params.internal_coeffs_imag[i],
        ));
    }

    // Wave numbers
    let k0 = 2.0 * PI / params.wavelength;

    // Complex refractive index: m = sqrt(εr * μr)
    let eps = Complex64::new(params.permittivity_real, params.permittivity_imag);
    let mu = Complex64::new(params.permeability_real, params.permeability_imag);
    let m = (eps * mu).sqrt();
    let k1 = k0 * m;

    // Grid setup
    let view_size = params.view_size;
    let half_size = view_size / 2.0;
    let x_min = -half_size;
    let x_max = half_size;
    let y_min = -half_size;
    let y_max = half_size;
    let dx = view_size / (GRID_SIZE as f64);

    // ===== PHASE 1: Build splines for Bessel functions =====

    // Maximum radius on the grid (corner of the view)
    let r_max = (2.0_f64).sqrt() * half_size;

    // Exterior: Hankel H_n(k0*r) for r in [RADIUS, r_max]
    let hankel_splines: Vec<BesselSpline> = params
        .orders
        .iter()
        .map(|&n| BesselSpline::new_hankel(n, k0, RADIUS, r_max, SPLINE_POINTS))
        .collect();

    // Interior: Bessel J_n(k1*r) for r in [0, RADIUS]
    // Use small epsilon to avoid r=0 singularity issues
    let bessel_splines: Vec<BesselSpline> = params
        .orders
        .iter()
        .map(|&n| BesselSpline::new_bessel_j(n, k1, 1e-10, RADIUS, SPLINE_POINTS))
        .collect();

    // ===== PHASE 2: Compute field using splines + y-symmetry =====

    let total_points = GRID_SIZE * GRID_SIZE;
    let mut field_real = vec![0.0; total_points];
    let mut field_imag = vec![0.0; total_points];

    // Pre-allocate exp(i*n*θ) table — reused for every grid point
    let mut exp_table = vec![Complex64::new(0.0, 0.0); n_coeffs];

    // Only compute top half (y >= 0), mirror to bottom half
    for iy in 0..GRID_SIZE / 2 {
        let y = y_max - (iy as f64 + 0.5) * dx;
        let row_offset = iy * GRID_SIZE;

        // Symmetric row index
        let iy_sym = GRID_SIZE - 1 - iy;
        let row_offset_sym = iy_sym * GRID_SIZE;

        for ix in 0..GRID_SIZE {
            let x = x_min + (ix as f64 + 0.5) * dx;
            let r = (x * x + y * y).sqrt();

            // cos(θ) = x/r, sin(θ) = y/r — avoids atan2 + cos/sin calls
            let inv_r = 1.0 / r;
            let cos_theta = x * inv_r;
            let sin_theta = y * inv_r;

            // Fill exp(i*n*θ) table using recurrence (0 trig calls)
            fill_exp_intheta(&mut exp_table, cos_theta, sin_theta);

            let idx = row_offset + ix;
            let idx_sym = row_offset_sym + ix;

            let field = if r < RADIUS {
                // Interior field using Bessel J splines
                compute_interior_field_spline(&c_n, &bessel_splines, r, &exp_table)
            } else {
                // Exterior field using Hankel splines
                compute_exterior_field_spline(&b_n, &hankel_splines, k0, r, cos_theta, &exp_table)
            };

            // Fill both point and its y-symmetric counterpart
            field_real[idx] = field.re;
            field_imag[idx] = field.im;
            field_real[idx_sym] = field.re;
            field_imag[idx_sym] = field.im;
        }
    }

    FieldResult {
        field_real,
        field_imag,
        grid_size: GRID_SIZE,
        view_size,
        x_min,
        x_max,
        y_min,
        y_max,
    }
}

/// Fill pre-allocated buffer with exp(i*n*θ) values for all orders.
/// Uses angle-addition recurrence: only needs cos(θ) and sin(θ), no per-order trig.
/// Orders are consecutive integers from -N to +N (2N+1 values).
#[inline]
fn fill_exp_intheta(table: &mut [Complex64], cos_t: f64, sin_t: f64) {
    let mid = table.len() / 2;

    // exp(i*0*θ) = 1
    table[mid] = Complex64::new(1.0, 0.0);

    // Step: exp(i*θ) for forward recurrence, exp(-i*θ) for backward
    let step_fwd = Complex64::new(cos_t, sin_t); // exp(+iθ)
    let step_bwd = Complex64::new(cos_t, -sin_t); // exp(-iθ)

    // Forward: exp(i*(n+1)*θ) = exp(i*n*θ) * exp(i*θ)
    for k in 1..=mid {
        table[mid + k] = table[mid + k - 1] * step_fwd;
    }
    // Backward: exp(i*(n-1)*θ) = exp(i*n*θ) * exp(-i*θ)
    for k in 1..=mid {
        table[mid - k] = table[mid - k + 1] * step_bwd;
    }
}

/// Compute interior field using spline-interpolated Bessel J values
#[inline]
fn compute_interior_field_spline(
    c_n: &[Complex64],
    bessel_splines: &[BesselSpline],
    r: f64,
    exp_table: &[Complex64],
) -> Complex64 {
    let mut field = Complex64::new(0.0, 0.0);

    for i in 0..c_n.len() {
        let jn = bessel_splines[i].eval(r);
        field += c_n[i] * jn * exp_table[i];
    }

    field
}

/// Compute exterior field using spline-interpolated Hankel values
#[inline]
fn compute_exterior_field_spline(
    b_n: &[Complex64],
    hankel_splines: &[BesselSpline],
    k0: f64,
    r: f64,
    cos_theta: f64,
    exp_table: &[Complex64],
) -> Complex64 {
    // Incident plane wave: exp(i*k0*x) where x = r*cos(theta)
    let k0x = k0 * r * cos_theta;
    let incident = Complex64::new(k0x.cos(), k0x.sin());

    // Scattered field: Σ_n b_n * H_n(k0*r) * exp(i*n*θ)
    let mut scattered = Complex64::new(0.0, 0.0);
    for i in 0..b_n.len() {
        let hn = hankel_splines[i].eval(r);
        scattered += b_n[i] * hn * exp_table[i];
    }

    incident + scattered
}

/// Compute the field inside the cylinder (non-optimized, for testing).
#[allow(dead_code)]
fn compute_interior_field(
    c_n: &[Complex64],
    orders: &[i32],
    k1: Complex64,
    r: f64,
    theta: f64,
) -> Complex64 {
    let k1_r = k1 * r;
    let mut field = Complex64::new(0.0, 0.0);

    for (i, &n) in orders.iter().enumerate() {
        // J_n(k1*r)
        let jn = bessel_j(n, k1_r);

        // exp(i*n*θ)
        let n_theta = n as f64 * theta;
        let exp_intheta = Complex64::new(n_theta.cos(), n_theta.sin());

        field += c_n[i] * jn * exp_intheta;
    }

    field
}

/// Compute the field outside the cylinder (non-optimized, for testing).
#[allow(dead_code)]
fn compute_exterior_field(
    b_n: &[Complex64],
    orders: &[i32],
    k0: f64,
    r: f64,
    theta: f64,
) -> Complex64 {
    let k0_r = Complex64::new(k0 * r, 0.0);

    // Incident plane wave: exp(i*k0*x) where x = r*cos(theta)
    let k0x = k0 * r * theta.cos();
    let incident = Complex64::new(k0x.cos(), k0x.sin());

    // Scattered field: Σ_n b_n * H_n(k0*r) * exp(i*n*θ)
    let mut scattered = Complex64::new(0.0, 0.0);
    for (i, &n) in orders.iter().enumerate() {
        let hn = hankel1(n, k0_r);
        let n_theta = n as f64 * theta;
        let exp_intheta = Complex64::new(n_theta.cos(), n_theta.sin());
        scattered += b_n[i] * hn * exp_intheta;
    }

    incident + scattered
}

#[allow(dead_code)]
fn linspace(start: f64, end: f64, n_points: usize) -> Vec<f64> {
    if n_points == 0 {
        return Vec::new();
    }
    let step = (end - start) / (n_points - 1) as f64;
    let mut vec = Vec::with_capacity(n_points);
    for i in 0..n_points {
        vec.push(start + i as f64 * step);
    }
    vec
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::scattering::{
        calculate_scattering, Material, Polarization, ScatteringParams, PERMEABILITY_IM_MAX,
        PERMEABILITY_IM_MIN, PERMEABILITY_RE_MAX, PERMEABILITY_RE_MIN, PERMITTIVITY_IM_MAX,
        PERMITTIVITY_IM_MIN, PERMITTIVITY_RE_MAX, PERMITTIVITY_RE_MIN, WAVELENGTH_MAX,
        WAVELENGTH_MIN,
    };

    /// Sweep permittivity, permeability, and wavelength to verify the boundary
    /// condition: interior and exterior fields must match at the cylinder surface.
    /// Sweep ranges are derived from the UI parameter bounds (single source of truth),
    /// clamped to a sub-range where Bessel function arguments stay within our
    /// implementation's precision threshold.
    #[test]
    fn test_tm_boundary_condition() {
        let sweep_eps_re = linspace(PERMITTIVITY_RE_MIN, PERMITTIVITY_RE_MAX, 4);
        let sweep_eps_im = linspace(PERMITTIVITY_IM_MIN, PERMITTIVITY_IM_MAX, 4);
        let sweep_mu_re = linspace(PERMEABILITY_RE_MIN, PERMEABILITY_RE_MAX, 4);
        let sweep_mu_im = linspace(PERMEABILITY_IM_MIN, PERMEABILITY_IM_MAX, 4);
        let wavelengths = [
            WAVELENGTH_MIN.max(0.3),
            (WAVELENGTH_MIN + WAVELENGTH_MAX) / 2.0,
            WAVELENGTH_MAX,
        ];
        let max_order = 21;
        let theta_vec = linspace(0.0, 2.0 * PI, 20);
        let tol = 1e-5;

        let mut n_configs = 0u32;
        let mut n_failures = 0u32;
        let mut worst_err = 0.0_f64;
        let mut worst_config = String::new();

        for &eps_re in &sweep_eps_re {
            for &eps_im in &sweep_eps_im {
                for &mu_re in &sweep_mu_re {
                    for &mu_im in &sweep_mu_im {
                        for &wl in &wavelengths {
                            n_configs += 1;

                            let params = ScatteringParams {
                                wavelength: wl,
                                material: Material {
                                    permittivity_real: eps_re,
                                    permittivity_imag: eps_im,
                                    permeability_real: mu_re,
                                    permeability_imag: mu_im,
                                },
                                polarization: Polarization::TM,
                                max_order,
                            };

                            let result = calculate_scattering(&params);
                            let k0 = 2.0 * PI / wl;
                            let kn = Complex64::new(k0, 0.0) * params.material.refractive_index();

                            let mut max_diff = 0.0_f64;
                            for &theta in &theta_vec {
                                let ext = compute_exterior_field(
                                    &result.scattering_coefficients,
                                    &result.orders,
                                    k0,
                                    RADIUS,
                                    theta,
                                );
                                let int = compute_interior_field(
                                    &result.internal_coefficients,
                                    &result.orders,
                                    kn,
                                    RADIUS,
                                    theta,
                                );
                                max_diff = max_diff.max((int - ext).norm());
                            }

                            if max_diff > tol {
                                n_failures += 1;
                                if max_diff > worst_err {
                                    worst_err = max_diff;
                                    worst_config = format!(
                                        "wl={wl} eps=({eps_re},{eps_im}) mu=({mu_re},{mu_im})"
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }

        let expected_configs = sweep_eps_re.len()
            * sweep_eps_im.len()
            * sweep_mu_re.len()
            * sweep_mu_im.len()
            * wavelengths.len();
        assert_eq!(
            n_configs as usize, expected_configs,
            "Unexpected config count"
        );

        assert_eq!(
            n_failures, 0,
            "{n_failures}/{n_configs} configs exceed tol {tol:.0e}. \
             Worst error: {worst_err:.4e} at [{worst_config}]"
        );
    }

    /// TE polarization sweep — same configurations as TM.
    /// Exercises the TE branch in calculate_coefficients_for_order.
    #[test]
    fn test_te_boundary_condition() {
        let sweep_eps_re = linspace(PERMITTIVITY_RE_MIN, PERMITTIVITY_RE_MAX, 4);
        let sweep_eps_im = linspace(PERMITTIVITY_IM_MIN, PERMITTIVITY_IM_MAX, 4);
        let sweep_mu_re = linspace(PERMEABILITY_RE_MIN, PERMEABILITY_RE_MAX, 4);
        let sweep_mu_im = linspace(PERMEABILITY_IM_MIN, PERMEABILITY_IM_MAX, 4);
        let wavelengths = [
            WAVELENGTH_MIN.max(0.3),
            (WAVELENGTH_MIN + WAVELENGTH_MAX) / 2.0,
            WAVELENGTH_MAX,
        ];
        let max_order = 21;
        let theta_vec = linspace(0.0, 2.0 * PI, 20);
        let tol = 1e-5;

        let mut n_configs = 0u32;
        let mut n_failures = 0u32;
        let mut worst_err = 0.0_f64;
        let mut worst_config = String::new();

        for &eps_re in &sweep_eps_re {
            for &eps_im in &sweep_eps_im {
                for &mu_re in &sweep_mu_re {
                    for &mu_im in &sweep_mu_im {
                        for &wl in &wavelengths {
                            n_configs += 1;

                            let params = ScatteringParams {
                                wavelength: wl,
                                material: Material {
                                    permittivity_real: eps_re,
                                    permittivity_imag: eps_im,
                                    permeability_real: mu_re,
                                    permeability_imag: mu_im,
                                },
                                polarization: Polarization::TE,
                                max_order,
                            };

                            let result = calculate_scattering(&params);
                            let k0 = 2.0 * PI / wl;
                            let kn = Complex64::new(k0, 0.0) * params.material.refractive_index();

                            let mut max_diff = 0.0_f64;
                            for &theta in &theta_vec {
                                let ext = compute_exterior_field(
                                    &result.scattering_coefficients,
                                    &result.orders,
                                    k0,
                                    RADIUS,
                                    theta,
                                );
                                let int = compute_interior_field(
                                    &result.internal_coefficients,
                                    &result.orders,
                                    kn,
                                    RADIUS,
                                    theta,
                                );
                                max_diff = max_diff.max((int - ext).norm());
                            }

                            if max_diff > tol {
                                n_failures += 1;
                                if max_diff > worst_err {
                                    worst_err = max_diff;
                                    worst_config = format!(
                                        "wl={wl} eps=({eps_re},{eps_im}) mu=({mu_re},{mu_im})"
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }

        let expected_configs = sweep_eps_re.len()
            * sweep_eps_im.len()
            * sweep_mu_re.len()
            * sweep_mu_im.len()
            * wavelengths.len();
        assert_eq!(
            n_configs as usize, expected_configs,
            "Unexpected config count"
        );

        assert_eq!(
            n_failures, 0,
            "{n_failures}/{n_configs} configs exceed tol {tol:.0e}. \
             Worst error: {worst_err:.4e} at [{worst_config}]"
        );
    }

    #[test]
    fn test_field_values_nonzero() {
        let params = ScatteringParams {
            wavelength: 1.0,
            material: Material {
                permittivity_real: 4.0,
                permittivity_imag: 0.5,
                permeability_real: 1.0,
                permeability_imag: -0.2,
            },
            polarization: Polarization::TM,
            max_order: 10,
        };

        let scattering = calculate_scattering(&params);

        // Build field params
        let field_params = FieldParams {
            wavelength: params.wavelength,
            permittivity_real: params.material.permittivity_real,
            permittivity_imag: params.material.permittivity_imag,
            permeability_real: params.material.permeability_real,
            permeability_imag: params.material.permeability_imag,
            incident_coeffs_real: scattering
                .incident_coefficients
                .iter()
                .map(|c| c.re)
                .collect(),
            incident_coeffs_imag: scattering
                .incident_coefficients
                .iter()
                .map(|c| c.im)
                .collect(),
            scattering_coeffs_real: scattering
                .scattering_coefficients
                .iter()
                .map(|c| c.re)
                .collect(),
            scattering_coeffs_imag: scattering
                .scattering_coefficients
                .iter()
                .map(|c| c.im)
                .collect(),
            internal_coeffs_real: scattering
                .internal_coefficients
                .iter()
                .map(|c| c.re)
                .collect(),
            internal_coeffs_imag: scattering
                .internal_coefficients
                .iter()
                .map(|c| c.im)
                .collect(),
            orders: scattering.orders,
            view_size: 5.0_f64,
        };

        let result = compute_field(&field_params);

        // Check that we have non-trivial field values
        let mut min_mag = f64::MAX;
        let mut max_mag = f64::MIN;

        for i in 0..result.field_real.len() {
            let mag = (result.field_real[i].powi(2) + result.field_imag[i].powi(2)).sqrt();
            min_mag = min_mag.min(mag);
            max_mag = max_mag.max(mag);
        }

        assert!(max_mag > min_mag, "Field should have variation");
        assert!(
            max_mag > 0.1,
            "Field maximum should be significant, got {}",
            max_mag
        );
    }
}
