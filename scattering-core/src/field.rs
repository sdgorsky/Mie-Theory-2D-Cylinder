//! Electric field computation for 2D cylinder scattering.
//!
//! Computes the total electric field (incident + scattered) on a 2D grid
//! using the cylindrical wave expansion coefficients.
//!
//! Optimizations:
//! - Cubic spline interpolation for Bessel/Hankel functions

use crate::bessel::{bessel_j, hankel1};
use crate::scattering::RADIUS;
use crate::sources::{Source, SourceKind};
use num_complex::Complex64;
use std::f64::consts::PI;

/// Grid resolution for field computation
pub const GRID_SIZE: usize = 512;

/// Default view size in units of cylinder diameter (5D x 5D)
pub const DEFAULT_VIEW_SIZE: f64 = 5.0;

/// Number of points in radial spline grids
const SPLINE_POINTS: usize = 128;
/// Fewer knots for interior J_n(k1*r) splines (smaller range, complex k1 is expensive)
const SPLINE_POINTS_INTERIOR: usize = 64;

// ============================================================================
// Cubic Spline Interpolation for Bessel Functions
// ============================================================================

/// Cubic spline for interpolating complex-valued Bessel/Hankel functions.
/// Stores spline coefficients for both real and imaginary parts.
struct BesselSpline {
    r_min: f64,
    inv_dr: f64,
    dr: f64,
    /// Spline coefficients [a, b, c, d] for each interval (real part)
    coeffs_re: Vec<[f64; 4]>,
    /// Spline coefficients [a, b, c, d] for each interval (imaginary part)
    coeffs_im: Vec<[f64; 4]>,
}

impl BesselSpline {
    /// Build a spline by sampling `f(r)` at `n_points` uniformly spaced points
    /// over `[r_min, r_max]`. All specialized constructors delegate here.
    fn from_fn(r_min: f64, r_max: f64, n_points: usize, f: impl Fn(f64) -> Complex64) -> Self {
        let dr = (r_max - r_min) / (n_points - 1) as f64;
        let mut values_re = Vec::with_capacity(n_points);
        let mut values_im = Vec::with_capacity(n_points);

        for i in 0..n_points {
            let r = r_min + i as f64 * dr;
            let v = f(r);
            values_re.push(v.re);
            values_im.push(v.im);
        }

        Self::from_values(r_min, dr, &values_re, &values_im)
    }

    /// Build spline from precomputed values
    fn from_values(r_min: f64, dr: f64, values_re: &[f64], values_im: &[f64]) -> Self {
        let coeffs_re = Self::compute_spline_coefficients(values_re, dr);
        let coeffs_im = Self::compute_spline_coefficients(values_im, dr);
        BesselSpline {
            r_min,
            inv_dr: 1.0 / dr,
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
        let t = (r - self.r_min) * self.inv_dr;
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
// Modal Spline Cache
// ============================================================================

/// Key identifying the parameters that determine modal spline values.
/// If the key matches, the cached splines can be reused.
#[derive(Clone, PartialEq)]
struct ModalSplineKey {
    k0: u64,    // f64 bits for exact comparison
    k1_re: u64, // complex k1 real part bits
    k1_im: u64, // complex k1 imag part bits
    max_order: usize,
    view_size: u64, // affects r_max for exterior splines
}

impl ModalSplineKey {
    fn new(k0: f64, k1: Complex64, max_order: usize, view_size: f64) -> Self {
        Self {
            k0: k0.to_bits(),
            k1_re: k1.re.to_bits(),
            k1_im: k1.im.to_bits(),
            max_order,
            view_size: view_size.to_bits(),
        }
    }
}

/// Minimum distance for incident field splines. Points closer than this to the
/// dipole source get NaN (handled gracefully in the FE colormap).
/// Must not be too small: the Miller backward recurrence for K_n(z) needs
/// ~O(1/|z|) iterations, so very small R causes expensive spline construction.
/// 0.02 gives ~3-4 pixel exclusion zone at 512×512 on a 5×5 grid.
const INCIDENT_R_MIN: f64 = 0.02;

/// Cached modal splines — reused across frames when material/wavelength/order don't change.
pub struct ModalSplineCache {
    key: ModalSplineKey,
    k0: f64,
    k1: Complex64,
    /// Max radius for exterior incident splines (diagonal of 2× the view).
    ext_r_max: f64,
    /// Max radius for interior incident splines (2× cylinder radius).
    int_r_max: f64,
    /// Exterior modal: H_n(k0*r) for n=0..N, r in [RADIUS, r_max] (real k0)
    hankel_splines: Vec<BesselSpline>,
    /// Interior modal: J_n(k1*r) for n=0..N, r in [0, RADIUS] (complex k1)
    bessel_splines: Vec<BesselSpline>,
    /// Incident field splines — built lazily on first dipole use.
    /// Exterior dipole splines (real k0, cheap to build):
    incident_h0_k0: Option<BesselSpline>,
    incident_rh1_k0: Option<BesselSpline>,
    /// Interior dipole splines (complex k1, expensive to build):
    incident_h0_k1: Option<BesselSpline>,
    incident_rh1_k1: Option<BesselSpline>,
}

/// Knot counts for incident field splines.
const INCIDENT_SPLINE_POINTS_EXT: usize = 256;
const INCIDENT_SPLINE_POINTS_INT: usize = 64;

impl ModalSplineCache {
    fn build(k0: f64, k1: Complex64, max_order: usize, view_size: f64) -> Self {
        let half_size = view_size / 2.0;
        let r_max = (2.0_f64).sqrt() * half_size;
        let ext_r_max = 2.0 * r_max;
        let int_r_max = 2.0 * RADIUS;

        // Always needed: modal splines for the scattering expansion
        let hankel_splines: Vec<BesselSpline> = (0..=max_order as i32)
            .map(|n| {
                BesselSpline::from_fn(RADIUS, r_max, SPLINE_POINTS, |r| {
                    hankel1(n, Complex64::new(k0 * r, 0.0))
                })
            })
            .collect();

        let bessel_splines: Vec<BesselSpline> = (0..=max_order as i32)
            .map(|n| {
                BesselSpline::from_fn(1e-10, RADIUS, SPLINE_POINTS_INTERIOR, |r| {
                    bessel_j(n, k1 * r)
                })
            })
            .collect();

        // Incident field splines built lazily — only when a dipole source is used
        Self {
            key: ModalSplineKey::new(k0, k1, max_order, view_size),
            k0,
            k1,
            ext_r_max,
            int_r_max,
            hankel_splines,
            bessel_splines,
            incident_h0_k0: None,
            incident_rh1_k0: None,
            incident_h0_k1: None,
            incident_rh1_k1: None,
        }
    }

    /// Get or build the exterior H_0(k0*R) incident spline.
    fn get_incident_h0_k0(&mut self) -> &BesselSpline {
        let k0 = self.k0;
        let ext_r_max = self.ext_r_max;
        self.incident_h0_k0.get_or_insert_with(|| {
            BesselSpline::from_fn(INCIDENT_R_MIN, ext_r_max, INCIDENT_SPLINE_POINTS_EXT, |r| {
                hankel1(0, Complex64::new(k0 * r, 0.0))
            })
        })
    }

    /// Get or build the exterior R*H_1(k0*R) regularized incident spline.
    fn get_incident_rh1_k0(&mut self) -> &BesselSpline {
        let k0 = self.k0;
        let ext_r_max = self.ext_r_max;
        self.incident_rh1_k0.get_or_insert_with(|| {
            BesselSpline::from_fn(INCIDENT_R_MIN, ext_r_max, INCIDENT_SPLINE_POINTS_EXT, |r| {
                hankel1(1, Complex64::new(k0 * r, 0.0)) * r
            })
        })
    }

    /// Get or build the interior H_0(k1*R) incident spline.
    fn get_incident_h0_k1(&mut self) -> &BesselSpline {
        let k1 = self.k1;
        let int_r_max = self.int_r_max;
        self.incident_h0_k1.get_or_insert_with(|| {
            BesselSpline::from_fn(INCIDENT_R_MIN, int_r_max, INCIDENT_SPLINE_POINTS_INT, |r| {
                hankel1(0, k1 * r)
            })
        })
    }

    /// Get or build the interior R*H_1(k1*R) regularized incident spline.
    fn get_incident_rh1_k1(&mut self) -> &BesselSpline {
        let k1 = self.k1;
        let int_r_max = self.int_r_max;
        self.incident_rh1_k1.get_or_insert_with(|| {
            BesselSpline::from_fn(INCIDENT_R_MIN, int_r_max, INCIDENT_SPLINE_POINTS_INT, |r| {
                hankel1(1, k1 * r) * r
            })
        })
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
    /// Grid resolution (number of points per edge). Defaults to GRID_SIZE.
    pub grid_size: usize,
    /// Incident field source type
    pub source: Source,
}

/// Result of field computation.
#[derive(Debug, Clone)]
pub struct FieldResult {
    pub field_real: Vec<f64>, // Real part of complex field values (row-major, grid_size x grid_size)
    pub field_imag: Vec<f64>, // Imaginary part of complex field values
    pub grid_size: usize,     // Number of grid points on each edge
    pub view_size: f64,       // Physical dimension of the grid (5.0 means -2.5 to +2.5)
    // Grid extrema
    pub x_min: f64,
    pub x_max: f64,
    pub y_min: f64,
    pub y_max: f64,
}

/// Compute the electric field on a grid, reusing cached modal splines when possible.
///
/// If `cache` is `Some` and its key matches the current parameters, the cached
/// splines are reused (skipping expensive Bessel function evaluations).
/// The cache is updated in-place when a rebuild is needed.
///
/// Uses optimized spline interpolation for Bessel/Hankel functions.
pub fn compute_field(params: &FieldParams, cache: &mut Option<ModalSplineCache>) -> FieldResult {
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
    let grid_size = params.grid_size;
    let view_size = params.view_size;
    let half_size = view_size / 2.0;
    let x_min = -half_size;
    let x_max = half_size;
    let y_min = -half_size;
    let y_max = half_size;
    let dx = view_size / (grid_size as f64);

    // ===== PHASE 1: Modal splines (cached) =====
    // Exploit symmetry: f_{-n}(r) = (-1)^n f_n(r) for both J_n and H_n^(1).
    // Build splines only for non-negative orders 0..N (N+1 splines instead of 2N+1).
    // Reuse cached splines when material/wavelength/order haven't changed.

    let max_order = n_coeffs / 2; // N, where orders span -N..+N

    let key = ModalSplineKey::new(k0, k1, max_order, view_size);
    if cache.as_ref().is_none_or(|c| c.key != key) {
        *cache = Some(ModalSplineCache::build(k0, k1, max_order, view_size));
    }

    // ===== PHASE 2: Build dipole incident field spline (if applicable) =====

    let source = params.source;
    let source_interior = source.is_interior();

    let incident_mode = build_incident_mode(source, source_interior, k0, k1);

    // Lazily build the incident splines needed for this source type.
    {
        let modal = cache.as_mut().unwrap();
        match source.kind {
            SourceKind::DipoleEz { .. } => {
                if source_interior {
                    modal.get_incident_h0_k1();
                } else {
                    modal.get_incident_h0_k0();
                }
            }
            SourceKind::DipoleExy { .. } => {
                if source_interior {
                    modal.get_incident_rh1_k1();
                } else {
                    modal.get_incident_rh1_k0();
                }
            }
            _ => {}
        }
    }

    // Immutable borrow for the grid loop
    let modal = cache.as_ref().unwrap();
    let hankel_splines = &modal.hankel_splines;
    let bessel_splines = &modal.bessel_splines;
    // Dummy fallback for plane wave (never actually evaluated)
    let dummy = &modal.hankel_splines[0];
    let h0_spline = if source_interior {
        modal.incident_h0_k1.as_ref().unwrap_or(dummy)
    } else {
        modal.incident_h0_k0.as_ref().unwrap_or(dummy)
    };
    let rh1_spline = if source_interior {
        modal.incident_rh1_k1.as_ref().unwrap_or(dummy)
    } else {
        modal.incident_rh1_k0.as_ref().unwrap_or(dummy)
    };

    // ===== PHASE 3: Compute field on full grid =====

    let total_points = grid_size * grid_size;
    let mut field_real = vec![0.0; total_points];
    let mut field_imag = vec![0.0; total_points];

    // Pre-allocate exp(i*k*θ) table for k = 0..N — reused for every grid point
    let mut exp_table = vec![Complex64::new(0.0, 0.0); max_order + 1];

    // Precompute x values (same for every row)
    let x_vals: Vec<f64> = (0..grid_size)
        .map(|ix| x_min + (ix as f64 + 0.5) * dx)
        .collect();

    // mid = N, the index of order 0 in the coefficient arrays
    let mid = max_order;

    for iy in 0..grid_size {
        let y = y_max - (iy as f64 + 0.5) * dx;
        let y_sq = y * y;
        let row_offset = iy * grid_size;

        for (ix, &x) in x_vals.iter().enumerate() {
            let r = (x * x + y_sq).sqrt();

            // cos(θ) = x/r, sin(θ) = y/r — avoids atan2 + cos/sin calls
            let inv_r = 1.0 / r;
            let cos_theta = x * inv_r;
            let sin_theta = y * inv_r;

            // Fill exp(i*k*θ) for k = 0..N using forward recurrence only
            fill_exp_intheta_half(&mut exp_table, cos_theta, sin_theta);

            let idx = row_offset + ix;

            let field = if r < RADIUS {
                let modal_c = compute_modal_sum_symmetric(&c_n, bessel_splines, r, &exp_table, mid);
                if source_interior {
                    modal_c + eval_incident(&incident_mode, h0_spline, rh1_spline, x, y)
                } else {
                    modal_c
                }
            } else {
                let scattered =
                    compute_modal_sum_symmetric(&b_n, hankel_splines, r, &exp_table, mid);
                if source_interior {
                    scattered
                } else {
                    scattered + eval_incident(&incident_mode, h0_spline, rh1_spline, x, y)
                }
            };

            field_real[idx] = field.re;
            field_imag[idx] = field.im;
        }
    }

    FieldResult {
        field_real,
        field_imag,
        grid_size,
        view_size,
        x_min,
        x_max,
        y_min,
        y_max,
    }
}

/// Fill pre-allocated buffer with exp(i*k*θ) for k = 0..N (forward only).
/// Negative-order angular factors are obtained via conjugation in the modal sum.
#[inline]
fn fill_exp_intheta_half(table: &mut [Complex64], cos_t: f64, sin_t: f64) {
    table[0] = Complex64::new(1.0, 0.0);
    let step = Complex64::new(cos_t, sin_t); // exp(+iθ)
    for k in 1..table.len() {
        table[k] = table[k - 1] * step;
    }
}

/// Compute modal sum exploiting Bessel order symmetry:
///   f_{-n}(r) = (-1)^n f_n(r)
///
/// The coefficient array `coeffs` has 2N+1 entries for orders -N..+N,
/// with index `mid` corresponding to order 0.
/// The `splines` array has N+1 entries for orders 0..N.
/// The `exp_table` has N+1 entries: exp(i*k*θ) for k = 0..N.
///
/// Sum = coeffs[mid] * f_0(r)
///     + Σ_{k=1}^{N} f_k(r) * [coeffs[mid+k]*exp(ikθ) + (-1)^k * coeffs[mid-k]*exp(-ikθ)]
#[inline]
fn compute_modal_sum_symmetric(
    coeffs: &[Complex64],
    splines: &[BesselSpline],
    r: f64,
    exp_table: &[Complex64],
    mid: usize,
) -> Complex64 {
    // k = 0 term
    let mut sum = coeffs[mid] * splines[0].eval(r);

    // Paired terms k = 1..N
    let n_max = exp_table.len() - 1; // = N
    for k in 1..=n_max {
        let fk = splines[k].eval(r);
        let exp_pos = exp_table[k]; // exp(+ikθ)
                                    // exp(-ikθ) = conj(exp(+ikθ))
        let exp_neg = Complex64::new(exp_pos.re, -exp_pos.im);

        let combined = if k % 2 == 0 {
            // (-1)^k = +1
            coeffs[mid + k] * exp_pos + coeffs[mid - k] * exp_neg
        } else {
            // (-1)^k = -1
            coeffs[mid + k] * exp_pos - coeffs[mid - k] * exp_neg
        };

        sum += fk * combined;
    }

    sum
}

// ============================================================================
// Incident Field — spline-accelerated for dipole sources
// ============================================================================

/// Precomputed incident field evaluation strategy.
/// Dipole sources use cached spline-interpolated Hankel functions over distance R
/// from the dipole. The splines are built in `ModalSplineCache` and reused across
/// frames — only the dipole position (xs, ys) changes per frame.
enum IncidentMode {
    /// exp(ikx) — cheap, no spline needed.
    PlaneWave { k: Complex64 },
    /// (i/4) H_0(k R) with cached spline over R = distance from (xs, ys).
    DipoleEz {
        xs: f64,
        ys: f64,
        prefactor: Complex64,
    },
    /// (ik/4) H_1(k R) sin(φ_R - α) with cached regularized spline.
    /// sin(φ_R - α) = (dy cos α − dx sin α) / R  (avoids atan2 + sin).
    DipoleExy {
        xs: f64,
        ys: f64,
        prefactor: Complex64,
        cos_alpha: f64,
        sin_alpha: f64,
    },
}

fn build_incident_mode(
    source: Source,
    source_interior: bool,
    k0: f64,
    k1: Complex64,
) -> IncidentMode {
    match source.kind {
        SourceKind::PlaneWaveTM | SourceKind::PlaneWaveTE => {
            let k = if source_interior {
                k1
            } else {
                Complex64::new(k0, 0.0)
            };
            IncidentMode::PlaneWave { k }
        }
        SourceKind::DipoleEz { xs, ys } => IncidentMode::DipoleEz {
            xs,
            ys,
            prefactor: Complex64::i() / 4.0,
        },
        SourceKind::DipoleExy { xs, ys, alpha } => {
            let k = if source_interior {
                k1
            } else {
                Complex64::new(k0, 0.0)
            };
            IncidentMode::DipoleExy {
                xs,
                ys,
                prefactor: Complex64::i() * k / 4.0,
                cos_alpha: alpha.cos(),
                sin_alpha: alpha.sin(),
            }
        }
    }
}

/// Evaluate incident field at (x, y) using the precomputed mode and cached splines.
/// Returns NaN for points too close to the dipole source (R < INCIDENT_R_MIN).
#[inline]
fn eval_incident(
    mode: &IncidentMode,
    h0_spline: &BesselSpline,
    rh1_spline: &BesselSpline,
    x: f64,
    y: f64,
) -> Complex64 {
    match mode {
        IncidentMode::PlaneWave { k } => (Complex64::i() * k * x).exp(),
        IncidentMode::DipoleEz { xs, ys, prefactor } => {
            let ddx = x - xs;
            let ddy = y - ys;
            let dist = (ddx * ddx + ddy * ddy).sqrt();
            if dist < INCIDENT_R_MIN {
                return Complex64::new(f64::NAN, f64::NAN);
            }
            prefactor * h0_spline.eval(dist)
        }
        IncidentMode::DipoleExy {
            xs,
            ys,
            prefactor,
            cos_alpha,
            sin_alpha,
        } => {
            let ddx = x - xs;
            let ddy = y - ys;
            let r_sq = ddx * ddx + ddy * ddy;
            if r_sq < INCIDENT_R_MIN * INCIDENT_R_MIN {
                return Complex64::new(f64::NAN, f64::NAN);
            }
            // Spline stores r * H_1(k*r) (regularized — smooth near r=0).
            // sin(φ_R - α) = (dy cos α − dx sin α) / R.
            // Total: prefactor * [r * H_1] * (dy cos α − dx sin α) / r²
            let dist = r_sq.sqrt();
            let rh1 = rh1_spline.eval(dist);
            let angular = ddy * cos_alpha - ddx * sin_alpha;
            prefactor * rh1 * angular / r_sq
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bessel::{bessel_j, hankel1};

    /// Compute the field inside the cylinder (non-optimized, for testing).
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

    use crate::scattering::{
        calculate_scattering, Material, ScatteringParams, PERMEABILITY_IM_MAX, PERMEABILITY_IM_MIN,
        PERMEABILITY_RE_MAX, PERMEABILITY_RE_MIN, PERMITTIVITY_IM_MAX, PERMITTIVITY_IM_MIN,
        PERMITTIVITY_RE_MAX, PERMITTIVITY_RE_MIN, WAVELENGTH_MAX, WAVELENGTH_MIN,
    };
    use crate::sources::{Source, SourceKind};

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
                                max_order,
                                source: Source::new(SourceKind::PlaneWaveTM, RADIUS),
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
                                max_order,
                                source: Source::new(SourceKind::PlaneWaveTE, RADIUS),
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
            max_order: 10,
            source: Source::new(SourceKind::PlaneWaveTM, RADIUS),
        };

        let scattering = calculate_scattering(&params);

        // Build field params
        let field_params = FieldParams {
            wavelength: params.wavelength,
            permittivity_real: params.material.permittivity_real,
            permittivity_imag: params.material.permittivity_imag,
            permeability_real: params.material.permeability_real,
            permeability_imag: params.material.permeability_imag,
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
            grid_size: GRID_SIZE,
            source: Source::new(SourceKind::PlaneWaveTM, RADIUS),
        };

        let result = compute_field(&field_params, &mut None);

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

    use crate::sources::{compute_incident_field, Domain};

    /// Generalized exterior total field at (r, theta) for any source.
    /// Total = incident (if exterior source) + scattered modal sum.
    fn compute_total_exterior_field(
        b_n: &[Complex64],
        orders: &[i32],
        k0: f64,
        r: f64,
        theta: f64,
        source: Source,
    ) -> Complex64 {
        let k0_r = Complex64::new(k0 * r, 0.0);

        let mut scattered = Complex64::new(0.0, 0.0);
        for (i, &n) in orders.iter().enumerate() {
            let hn = hankel1(n, k0_r);
            let n_theta = n as f64 * theta;
            let exp_intheta = Complex64::new(n_theta.cos(), n_theta.sin());
            scattered += b_n[i] * hn * exp_intheta;
        }

        let incident = if source.domain == Domain::Exterior {
            let x = r * theta.cos();
            let y = r * theta.sin();
            compute_incident_field(source, Complex64::new(k0, 0.0), x, y)
        } else {
            Complex64::new(0.0, 0.0)
        };

        incident + scattered
    }

    /// Generalized interior total field at (r, theta) for any source.
    /// Total = incident (if interior source) + internal modal sum.
    fn compute_total_interior_field(
        c_n: &[Complex64],
        orders: &[i32],
        k1: Complex64,
        r: f64,
        theta: f64,
        source: Source,
    ) -> Complex64 {
        let k1_r = k1 * r;
        let mut modal = Complex64::new(0.0, 0.0);

        for (i, &n) in orders.iter().enumerate() {
            let jn = bessel_j(n, k1_r);
            let n_theta = n as f64 * theta;
            let exp_intheta = Complex64::new(n_theta.cos(), n_theta.sin());
            modal += c_n[i] * jn * exp_intheta;
        }

        let incident = if source.domain == Domain::Interior {
            let x = r * theta.cos();
            let y = r * theta.sin();
            compute_incident_field(source, k1, x, y)
        } else {
            Complex64::new(0.0, 0.0)
        };

        incident + modal
    }

    /// Verify boundary condition (interior == exterior at r = RADIUS) for
    /// four dipole configurations with a single non-trivial material.
    #[test]
    fn test_dipole_boundary_conditions() {
        // Real dielectric+magnetic material — avoids complex Bessel arguments
        // where our hankel1 implementation has a known bug (NaN for |Im(z)| > ~0.01).
        let material = Material {
            permittivity_real: 2.3,
            permittivity_imag: 0.0,
            permeability_real: 1.4,
            permeability_imag: 0.0,
        };
        let wavelength = 1.2;
        let max_order = 25;
        let theta_vec = linspace(0.0, 2.0 * PI, 40);
        let tol = 1e-5;

        let k0 = 2.0 * PI / wavelength;
        let kn = Complex64::new(k0, 0.0) * material.refractive_index();

        // Four dipole configurations: (label, SourceKind)
        let configs: Vec<(&str, SourceKind)> = vec![
            (
                "DipoleEz exterior",
                SourceKind::DipoleEz { xs: 1.7, ys: -0.9 },
            ),
            (
                "DipoleEz interior",
                SourceKind::DipoleEz { xs: 0.15, ys: -0.1 },
            ),
            (
                "DipoleExy exterior",
                SourceKind::DipoleExy {
                    xs: -1.3,
                    ys: 0.6,
                    alpha: 0.73,
                },
            ),
            (
                "DipoleExy interior",
                SourceKind::DipoleExy {
                    xs: -0.12,
                    ys: 0.2,
                    alpha: -1.1,
                },
            ),
        ];

        for (label, kind) in &configs {
            let source = Source::new(*kind, RADIUS);

            let params = ScatteringParams {
                wavelength,
                material,
                max_order,
                source,
            };

            let result = calculate_scattering(&params);

            let mut max_diff = 0.0_f64;
            let mut worst_theta = 0.0_f64;
            for &theta in &theta_vec {
                let ext = compute_total_exterior_field(
                    &result.scattering_coefficients,
                    &result.orders,
                    k0,
                    RADIUS,
                    theta,
                    source,
                );
                let int = compute_total_interior_field(
                    &result.internal_coefficients,
                    &result.orders,
                    kn,
                    RADIUS,
                    theta,
                    source,
                );
                let diff = (int - ext).norm();
                if diff > max_diff {
                    max_diff = diff;
                    worst_theta = theta;
                }
            }

            assert!(
                max_diff < tol,
                "{label}: boundary condition violated, max error = {max_diff:.4e} at theta = {worst_theta:.3}"
            );
        }
    }
}
