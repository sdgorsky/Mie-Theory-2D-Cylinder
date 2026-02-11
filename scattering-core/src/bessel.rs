//! Complex Bessel functions for electromagnetic scattering calculations.
//!
//! Pure Rust implementation following the architecture of the AMOS Fortran library
//! (D.E. Amos, Sandia National Laboratories), which underlies SciPy's Bessel routines.
//!
//! # Architecture
//!
//! The central design decision is that H^(1)_n(z) for complex z is computed via the
//! modified Bessel function K_n, rather than J_n + iY_n (which suffers from catastrophic
//! cancellation in the upper half-plane). The identity used is:
//!
//!   H^(1)_n(z) = (2/pi*i) * exp(-i*n*pi/2) * K_n(-iz)
//!
//! with two sub-paths for K_n depending on where w = -iz falls:
//!
//! - **Re(w) > 0 (ZBKNU path):** K_n computed directly via Miller backward recurrence
//!   for |z| > 2, or power series for |z| <= 2.
//! - **Re(w) <= 0 (ZACON path):** Analytic continuation into the left half-plane using
//!   K_n(w) = (-1)^n * K_n(-w) + pi*i * I_n(-w), where I_n is obtained through the
//!   connection formula I_n(z) = i^{-n} * J_n(iz).
//!
//! For real arguments (Im(z) = 0), H^(1)_n uses the simpler J_n + iY_n path, which is
//! exact for real inputs and avoids unnecessary trips through K.
//!
//! # References
//!
//! - D.E. Amos, "A Portable Package for Bessel Functions of a Complex Argument and
//!   Nonnegative Order", ACM Trans. Math. Software, Vol. 12, No. 3, 1986, pp. 265-273.
//!   <https://dl.acm.org/doi/10.1145/7921.214331>
//!
//! - D.E. Amos, "Algorithm 644: A Portable Package for Bessel Functions of a Complex
//!   Argument and Nonnegative Order", ACM Trans. Math. Software, Vol. 12, No. 3, 1986,
//!   pp. 265-273. (Accompanying source code for ZBESH, ZBKNU, ZACON, et al.)
//!
//! - DLMF (NIST Digital Library of Mathematical Functions), Chapter 10: Bessel Functions.
//!   <https://dlmf.nist.gov/10>
//!   Key identities:
//!     - 10.27.6 (I-J connection)
//!     - 10.34.2 (K analytic continuation)
//!     - 10.4.4 (J reflection)
//!     - 10.2.2 (J power series),
//!     - 10.17.5 (Hankel asymptotic)
//!
//! - Gonum (Go Numerical Library) translation of AMOS, used as a readable reference for
//!   the ZBKNU Miller backward recurrence algorithm.
//!   <https://github.com/gonum/gonum>

use num_complex::Complex64;
use std::f64::consts::PI;

/// Euler-Mascheroni constant
const EULER_GAMMA: f64 = 0.5772156649015329;

/// Threshold for switching from power series to asymptotic expansion.
/// For Hankel H1 (single series, no cancellation), this threshold works directly.
/// For Bessel J/Y (computed as H1±H2, cancellation for complex z), the effective
/// threshold is increased by the order to ensure accuracy.
const ASYMPTOTIC_THRESHOLD: f64 = 25.0;

/// Factorial function returning f64 (avoids integer overflow for n > 20)
fn factorial(n: u32) -> f64 {
    match n {
        0 | 1 => 1.0,
        _ => (2..=n as u64).fold(1.0_f64, |acc, k| acc * k as f64),
    }
}

/// Digamma function ψ(n) for positive integers
/// ψ(n) = -γ + Σ_{k=1}^{n-1} 1/k
fn digamma(n: u32) -> f64 {
    if n == 0 {
        return f64::NEG_INFINITY;
    }
    let mut result = -EULER_GAMMA;
    for k in 1..n {
        result += 1.0 / k as f64;
    }
    result
}

/// Asymptotic expansion for Hankel function of the first kind H^(1)_n(z).
///
/// DLMF 10.17.5: H^(1)_n(z) ~ sqrt(2/(πz)) * exp(i*ω) * Σ_k a_k(n) * (i/(8z))^k / k!
/// where ω = z - nπ/2 - π/4.
///
/// Recurrence: term_k = term_{k-1} * (4n² - (2k-1)²) * i / (8kz)
///
/// Uses optimal truncation: stops when terms start growing.
fn hankel1_asymptotic(n: i32, z: Complex64) -> Complex64 {
    let nf = n as f64;
    let four_n_sq = 4.0 * nf * nf;

    // ω = z - nπ/2 - π/4
    let omega = z - Complex64::new(nf * PI / 2.0 + PI / 4.0, 0.0);

    // sqrt(2/(πz))
    let prefactor = (Complex64::new(2.0 / PI, 0.0) / z).sqrt();

    // exp(iω)
    let phase = (Complex64::i() * omega).exp();

    // Recurrence factor: multiply by (4n²-(2k-1)²) * i / (8kz)
    // Note: i/z, NOT 1/(iz) — the sign matters since 1/(iz) = -i/z
    let i_over_z = Complex64::i() / z;
    let mut term = Complex64::new(1.0, 0.0);
    let mut sum = term;
    let mut prev_norm = term.norm();
    let mut decreasing = false;

    for k in 1..=60u32 {
        let twok_minus_1 = (2 * k - 1) as f64;
        let numerator = four_n_sq - twok_minus_1 * twok_minus_1;
        term = term * Complex64::new(numerator / (8.0 * k as f64), 0.0) * i_over_z;

        let term_norm = term.norm();
        if term_norm < prev_norm {
            decreasing = true;
        } else if decreasing {
            // Terms have started growing; stop before adding this term
            break;
        }

        sum += term;
        prev_norm = term_norm;

        if term_norm < 1e-15 * sum.norm() {
            break;
        }
    }

    prefactor * phase * sum
}

/// Asymptotic expansion for Hankel function of the second kind H^(2)_n(z).
///
/// H^(2)_n(z) ~ sqrt(2/(πz)) * exp(-i*ω) * Σ_k a_k(n) * (-i/(8z))^k / k!
/// where ω = z - nπ/2 - π/4.
///
/// Recurrence: term_k = term_{k-1} * (4n² - (2k-1)²) * (-i) / (8kz)
fn hankel2_asymptotic(n: i32, z: Complex64) -> Complex64 {
    let nf = n as f64;
    let four_n_sq = 4.0 * nf * nf;

    let omega = z - Complex64::new(nf * PI / 2.0 + PI / 4.0, 0.0);
    let prefactor = (Complex64::new(2.0 / PI, 0.0) / z).sqrt();
    let phase = (-Complex64::i() * omega).exp();

    // Recurrence factor: multiply by (4n²-(2k-1)²) * (-i) / (8kz)
    let neg_i_over_z = -Complex64::i() / z;
    let mut term = Complex64::new(1.0, 0.0);
    let mut sum = term;
    let mut prev_norm = term.norm();
    let mut decreasing = false;

    for k in 1..=60u32 {
        let twok_minus_1 = (2 * k - 1) as f64;
        let numerator = four_n_sq - twok_minus_1 * twok_minus_1;
        term = term * Complex64::new(numerator / (8.0 * k as f64), 0.0) * neg_i_over_z;

        let term_norm = term.norm();
        if term_norm < prev_norm {
            decreasing = true;
        } else if decreasing {
            break;
        }

        sum += term;
        prev_norm = term_norm;

        if term_norm < 1e-15 * sum.norm() {
            break;
        }
    }

    prefactor * phase * sum
}

/// J_n(z) via asymptotic Hankel expansions: J = (H1 + H2) / 2
///
/// For complex z with significant Im(z), one Hankel function dominates:
///   H^(1) ~ exp(iz) → decays as exp(-Im(z)) for Im(z) > 0
///   H^(2) ~ exp(-iz) → grows as exp(Im(z)) for Im(z) > 0
/// So when Im(z) > 0: J ≈ H^(2)/2, when Im(z) < 0: J ≈ H^(1)/2.
fn bessel_j_asymptotic(n: i32, z: Complex64) -> Complex64 {
    if z.im > 15.0 {
        // H^(2) dominates (grows as exp(Im(z))), H^(1) is negligible
        return hankel2_asymptotic(n, z) / 2.0;
    }
    if z.im < -15.0 {
        // H^(1) dominates (grows as exp(-Im(z)) = exp(|Im(z)|)), H^(2) is negligible
        return hankel1_asymptotic(n, z) / 2.0;
    }
    let h1 = hankel1_asymptotic(n, z);
    let h2 = hankel2_asymptotic(n, z);
    (h1 + h2) / 2.0
}

/// Y_n(z) via asymptotic Hankel expansions: Y = (H1 - H2) / (2i)
///
/// Dominant Hankel: same logic as bessel_j_asymptotic.
/// Y = (H1 - H2)/(2i). When H2 dominates: Y ≈ -H2/(2i).
/// When H1 dominates: Y ≈ H1/(2i).
fn bessel_y_asymptotic(n: i32, z: Complex64) -> Complex64 {
    if z.im > 15.0 {
        return -hankel2_asymptotic(n, z) / (2.0 * Complex64::i());
    }
    if z.im < -15.0 {
        return hankel1_asymptotic(n, z) / (2.0 * Complex64::i());
    }
    let h1 = hankel1_asymptotic(n, z);
    let h2 = hankel2_asymptotic(n, z);
    (h1 - h2) / (2.0 * Complex64::i())
}

/// Miller's backward recurrence for J_n(z).
///
/// Computes J_n(z) by starting from J_M ≈ 0 at a large order M and using the
/// recurrence backwards. Normalizes using J_0(z) + 2*Σ J_{2k}(z) = 1.
///
/// Numerically stable for real and moderate-complex z. Includes rescaling
/// to prevent overflow for complex z.
fn bessel_j_miller(n: i32, z: Complex64) -> Complex64 {
    let n_abs = n.unsigned_abs() as usize;
    // Starting order: must be well above both n and |z| for convergence
    let m = n_abs + 2 * z.norm() as usize + 40;

    let two_over_z = Complex64::new(2.0, 0.0) / z;

    // We iterate k from m down to 0, maintaining:
    //   j_prev = J_{k+1} (unnormalized)
    //   j_curr = J_{k+2} (unnormalized)
    // At each step we compute J_k = (2*(k+1)/z)*J_{k+1} - J_{k+2}
    let mut j_prev = Complex64::new(1.0, 0.0); // J_{m} (arbitrary start)
    let mut j_curr = Complex64::new(0.0, 0.0); // J_{m+1} = 0

    let mut j_n_unnorm = Complex64::new(0.0, 0.0);
    let mut norm_sum = Complex64::new(0.0, 0.0);

    for k in (0..=m).rev() {
        // At this point j_prev = J_{k} (for k=m, this is our initial J_m = 1)
        // Capture J_n
        if k == n_abs {
            j_n_unnorm = j_prev;
        }

        // Accumulate normalization: J_0 + 2*(J_2 + J_4 + ...)
        if k % 2 == 0 {
            if k == 0 {
                norm_sum += j_prev;
            } else {
                norm_sum += j_prev * 2.0;
            }
        }

        if k == 0 {
            break;
        }

        // Compute J_{k-1} = (2k/z)*J_k - J_{k+1}
        let j_new = two_over_z * k as f64 * j_prev - j_curr;

        // Rescale to prevent overflow
        let mag = j_new.norm();
        if mag > 1e200 {
            let inv_mag = 1.0 / mag;
            j_curr = j_prev * inv_mag;
            j_prev = j_new * inv_mag;
            j_n_unnorm *= inv_mag;
            norm_sum *= inv_mag;
        } else {
            j_curr = j_prev;
            j_prev = j_new;
        }
    }

    j_n_unnorm / norm_sum
}

/// Power series for J_n(z), accurate when |z| is small enough to avoid cancellation.
fn bessel_j_series(n: i32, z: Complex64) -> Complex64 {
    let z_half = z / 2.0;
    let neg_z_half_sq = -z_half * z_half;

    // (z/2)^n
    let mut prefix = Complex64::new(1.0, 0.0);
    for _ in 0..n {
        prefix *= z_half;
    }

    // k=0 term: 1/n!
    let mut term = Complex64::new(1.0 / factorial(n as u32), 0.0);
    let mut sum = term;

    // Sum terms until convergence using ratio recurrence
    for k in 1i32..500 {
        let denom = k as f64 * (n + k) as f64;
        term = term * neg_z_half_sq / Complex64::new(denom, 0.0);
        sum += term;

        if term.norm() < 1e-15 * sum.norm() {
            break;
        }
    }

    prefix * sum
}

/// Power series for modified Bessel function I_n(z).
///
/// I_n(z) = (z/2)^n * Σ_k (z²/4)^k / (k! * (n+k)!)
///
/// Same structure as J_n series but with +z²/4 instead of -z²/4.
/// No alternating signs means no cancellation — converges for all z.
fn bessel_i_series(n: i32, z: Complex64) -> Complex64 {
    let z_half = z / 2.0;
    let z_half_sq = z_half * z_half;

    // (z/2)^n
    let mut prefix = Complex64::new(1.0, 0.0);
    for _ in 0..n {
        prefix *= z_half;
    }

    // k=0 term: 1/n!
    let mut term = Complex64::new(1.0 / factorial(n as u32), 0.0);
    let mut sum = term;

    for k in 1i32..500 {
        let denom = k as f64 * (n + k) as f64;
        term = term * z_half_sq / Complex64::new(denom, 0.0);
        sum += term;

        if term.norm() < 1e-15 * sum.norm() {
            break;
        }
    }

    prefix * sum
}

/// I_n(z) via forward recurrence from I_0, I_1.
///
/// I_{k+1}(z) = I_{k-1}(z) - (2k/z)*I_k(z)
///
/// Forward recurrence is stable for I_n (the dominant solution for Re(z)>0).
/// K_0(z) and K_1(z) via Miller's backward recurrence (AMOS ZBKNU algorithm).
///
/// This is the algorithm AMOS uses for |z| > 2. It computes:
///   K_0(z) = coef * P_0 / CS
///   K_1(z) = K_0(z) * ((0.5 - P_1/P_0) / z + 1)
/// where coef = sqrt(π/(2z)) * exp(-z), and P_k/CS come from backward recurrence.
///
/// Converges to machine precision for all arg(z) in the RHP — no Stokes errors.
fn bessel_k01_miller(z: Complex64) -> (Complex64, Complex64) {
    let caz = z.norm();
    let tol = f64::EPSILON; // ~2.22e-16

    // Prefactor: sqrt(pi/(2z)) * exp(-z)
    let coef = (Complex64::new(PI / 2.0, 0.0) / z).sqrt() * (-z).exp();

    // --- Determine starting index FK ---
    let fhs_init = 0.25_f64; // |0.25 - dnu²| for dnu=0
    let mut fk: f64;

    // Precision parameters (for f64: 53-bit mantissa)
    let t1_prec = 52.8_f64; // 52 * log10(2) * 3.322 ≈ 52.8
    let t2 = (2.0 / 3.0) * t1_prec - 6.0; // ≈ 29.2

    if caz >= t2 {
        // ETEST forward loop (AMOS labels 140-160)
        let etest = 1.0 / (PI * caz * tol);
        fk = 1.0;
        if etest >= 1.0 {
            let mut fks = 2.0_f64;
            let mut ckr = 2.0 * caz + 2.0;
            let mut p1r = 0.0_f64;
            let mut p2r = 1.0_f64;
            let mut fhs = fhs_init;

            for _ in 1..=30 {
                let ak = fhs / fks;
                let cbr = ckr / (fk + 1.0);
                let ptr = p2r;
                p2r = cbr * p2r - p1r * ak;
                p1r = ptr;
                ckr += 2.0;
                fks += 2.0 * fk + 2.0;
                fhs += 2.0 * fk;
                fk += 1.0;
                if etest < p2r.abs() * fk {
                    break;
                }
            }

            // Safety margin based on angle of z
            let t1_angle = if z.re != 0.0 {
                (z.im / z.re).abs().atan()
            } else {
                PI / 2.0
            };
            // SPI = sqrt(6/pi) ≈ 1.9099
            fk += 1.90986 * t1_angle * (t2 / caz).sqrt();
        }
    } else {
        // Heuristic FK formula (AMOS label 170)
        let a2 = caz.sqrt();
        // FPI ≈ pi * 2/sqrt(3) ≈ 1.8977
        let mut ak = 1.89770 / (tol * a2.sqrt());
        let aa = 3.0 * t1_prec / (1.0 + caz);
        let bb = 14.7 * t1_prec / (28.0 + caz);
        ak = (ak.ln() + caz * aa.cos() / (1.0 + 0.008 * caz)) / bb.cos();
        fk = 0.12125 * ak * ak / caz + 1.5;
    }

    // --- Backward recurrence ---
    let k = fk as usize;
    let mut fk = k as f64;
    let mut fks = fk * fk;
    let fhs = 0.25_f64; // reset for backward recurrence (dnu=0)

    let mut p1 = Complex64::new(0.0, 0.0); // P_{k+1}
    let mut p2 = Complex64::new(tol, 0.0); // P_k
    let mut cs = p2; // normalization sum

    for _ in 0..k {
        let a1 = fks - fk; // fk*(fk-1)
        let ak = (fks + fk) / (a1 + fhs); // (fk²+fk)/(fk²-fk+0.25)
        let rak = 2.0 / (fk + 1.0);
        let cb = (Complex64::new(fk, 0.0) + z) * rak; // (fk + z) * 2/(fk+1)

        let p2_new = (cb * p2 - p1) * ak;
        p1 = p2;
        p2 = p2_new;
        cs += p2;

        fks = a1 - fk + 1.0; // (fk-1)²
        fk -= 1.0;
    }

    // Normalize: K_0 = coef * P_0 / CS
    let k0 = coef * p2 / cs;

    // K_1 from backward recurrence ratio: K_1 = K_0 * ((0.5 - P_1/P_0)/z + 1)
    let ratio = p1 / p2;
    let k1 = k0 * ((Complex64::new(0.5, 0.0) - ratio) / z + Complex64::new(1.0, 0.0));

    (k0, k1)
}

/// K_0(z) via power series (DLMF 10.31.2 for n=0).
///
/// K_0(z) = -(ln(z/2) + γ) * I_0(z) + Σ_{k=1}^∞ H_k * (z/2)^{2k} / (k!)²
///
/// where H_k = 1 + 1/2 + ... + 1/k (harmonic numbers).
fn bessel_k0_series(z: Complex64) -> Complex64 {
    let i0 = bessel_i_series(0, z);
    let ln_z_half = (z / 2.0).ln();

    let z_half_sq = (z / 2.0) * (z / 2.0);

    // Sum: Σ_{k=1}^∞ H_k * (z²/4)^k / (k!)²
    // base_k = (z²/4)^k / (k!)², ratio: base_k = base_{k-1} * z_half_sq / k²
    let mut base = Complex64::new(1.0, 0.0);
    let mut harmonic = 0.0_f64;
    let mut sum = Complex64::new(0.0, 0.0);

    for k in 1..500 {
        let kf = k as f64;
        base = base * z_half_sq / Complex64::new(kf * kf, 0.0);
        harmonic += 1.0 / kf;
        let term = base * harmonic;
        sum += term;

        if k > 5 && term.norm() < 1e-15 * sum.norm() {
            break;
        }
    }

    -(ln_z_half + Complex64::new(EULER_GAMMA, 0.0)) * i0 + sum
}

/// I_n(z) for complex z via connection to J_n.
///
/// Uses the identity I_n(z) = i^{-n} * J_n(iz) (DLMF 10.27.6).
/// This leverages the well-tested J_n code (series, Miller, asymptotic) and
/// avoids I_n-specific issues: series cancellation for complex z near the
/// imaginary axis, asymptotic Stokes errors, and forward recurrence instability.
pub fn bessel_i(n: i32, z: Complex64) -> Complex64 {
    let iz = Complex64::i() * z;
    // i^{-n} = (-i)^n: cycle of length 4
    let phase = match n.rem_euclid(4) {
        0 => Complex64::new(1.0, 0.0),  // i^0 = 1
        1 => Complex64::new(0.0, -1.0), // i^{-1} = -i
        2 => Complex64::new(-1.0, 0.0), // i^{-2} = -1
        3 => Complex64::new(0.0, 1.0),  // i^{-3} = i
        _ => unreachable!(),
    };
    phase * bessel_j(n, iz)
}

/// K_n(z) for Re(z) > 0.
///
/// Always computes K_0/K_1 first, then forward-recurs to K_n.
/// For |z| > 2: Miller backward recurrence (AMOS ZBKNU) — machine precision.
/// For |z| <= 2: K_0 from series, K_1 from Wronskian I_0*K_1 + I_1*K_0 = 1/z.
/// Forward recurrence K_{k+1} = K_{k-1} + (2k/z)*K_k is stable (K is dominant).
pub fn bessel_k(n: i32, z: Complex64) -> Complex64 {
    let n_abs = n.unsigned_abs() as i32; // K_{-n} = K_n
    let znorm = z.norm();

    // Compute K_0 and K_1, then forward-recur
    let (k0, k1) = if znorm > 2.0 {
        // Miller backward recurrence (AMOS ZBKNU): machine precision for all arg(z)
        bessel_k01_miller(z)
    } else {
        // K_0 from series, K_1 from Wronskian: I_0*K_1 + I_1*K_0 = 1/z
        let k0 = bessel_k0_series(z);
        let i0 = bessel_i_series(0, z);
        let i1 = bessel_i_series(1, z);
        let k1 = (Complex64::new(1.0, 0.0) / z - i1 * k0) / i0;
        (k0, k1)
    };

    if n_abs == 0 {
        return k0;
    }
    if n_abs == 1 {
        return k1;
    }

    // Forward recurrence: K_{k+1} = K_{k-1} + (2k/z)*K_k
    // K is dominant, so forward recurrence is stable.
    let two_over_z = Complex64::new(2.0, 0.0) / z;
    let mut k_prev = k0;
    let mut k_curr = k1;

    for k in 1..n_abs {
        let k_next = k_prev + two_over_z * k as f64 * k_curr;
        k_prev = k_curr;
        k_curr = k_next;
    }

    k_curr
}

/// Bessel function of the first kind J_n(z) for complex argument.
///
/// Uses three algorithms depending on the regime:
/// 1. Hankel asymptotic expansion for large |z| with small Im(z)
/// 2. Miller's backward recurrence for large real |z| (transition region)
/// 3. Power series for small |z| and complex z with large Im(z)
pub fn bessel_j(n: i32, z: Complex64) -> Complex64 {
    // Handle negative orders: J_{-n}(z) = (-1)^n * J_n(z)
    if n < 0 {
        let sign = if (-n) % 2 == 0 { 1.0 } else { -1.0 };
        return sign * bessel_j(-n, z);
    }

    // Reflect to positive real half-plane: J_n(-z) = (-1)^n * J_n(z) (integer n)
    // Avoids asymptotic expansion issues near the negative real axis (sqrt branch cut).
    if z.re < 0.0 {
        let sign = if n % 2 == 0 { 1.0 } else { -1.0 };
        return sign * bessel_j(n, -z);
    }

    // For very small z, use leading term
    if z.norm() < 1e-15 {
        return if n == 0 {
            Complex64::new(1.0, 0.0)
        } else {
            Complex64::new(0.0, 0.0)
        };
    }

    let znorm = z.norm();
    let nf = n as f64;

    // Asymptotic with dominant Hankel for complex z with large |Im(z)|.
    // One Hankel dominates by exp(2|Im(z)|), avoiding cancellation.
    // Requires both |z| > 2n+5 (enough asymptotic terms) and |z| > 10
    // (expansion needs |z| absolutely large enough for the series to converge).
    if z.im.abs() >= 5.0 && znorm >= (2.0 * nf + 5.0).max(10.0) {
        return bessel_j_asymptotic(n, z);
    }

    // Asymptotic J = (H1+H2)/2 for near-real z with |z| >> n²
    if znorm > nf * nf / 2.0 + nf + ASYMPTOTIC_THRESHOLD {
        return bessel_j_asymptotic(n, z);
    }

    // Power series: works when cancellation is bounded.
    // Cancellation depends on Re(z)²-Im(z)²: when Im(z) dominates, terms don't
    // alternate and the series converges without digit loss. When Re(z) dominates,
    // alternating sign cancellation loses ~sqrt(Re²-Im²)/2.3 digits.
    // With 15 f64 digits and needing 1e-6 accuracy (9 digits), allow up to 6 digits loss.
    let cancel_param = z.re * z.re - z.im * z.im;
    if cancel_param < 196.0 {
        // Effective real modulus < 14 → loses < 6 digits → keeps 9+ digits
        return bessel_j_series(n, z);
    }

    // Transition region: |Re(z)| large (series cancels), |z| not large enough
    // relative to n for asymptotic. Use Miller's backward recurrence.
    // Miller works well when |Im(z)| is moderate (normalization sum cancellation
    // loses ~|Im(z)|/2.3 digits, acceptable for |Im(z)| < ~20).
    bessel_j_miller(n, z)
}

/// Bessel function of the second kind Y_n(z) for complex argument.
pub fn bessel_y(n: i32, z: Complex64) -> Complex64 {
    // Handle negative orders: Y_{-n}(z) = (-1)^n * Y_n(z)
    if n < 0 {
        let sign = if (-n) % 2 == 0 { 1.0 } else { -1.0 };
        return sign * bessel_y(-n, z);
    }

    if z.norm() > ASYMPTOTIC_THRESHOLD {
        return bessel_y_asymptotic(n, z);
    }

    // Y_n(z) = (J_n(z) * cos(n*pi) - J_{-n}(z)) / sin(n*pi)
    // For integer n, use limit form:
    // Y_n(z) = (2/pi) * J_n(z) * ln(z/2) - (1/pi) * sum of correction terms

    let j_n = bessel_j(n, z);
    let ln_z_half = (z / 2.0).ln();

    // First part: (2/π) * J_n(z) * ln(z/2)
    let part1 = (2.0 / PI) * j_n * ln_z_half;

    // Second part: negative power series
    let z_half = z / 2.0;

    // (z/2)^n term
    let mut prefix = Complex64::new(1.0, 0.0);
    for _ in 0..n {
        prefix *= z_half;
    }

    // Sum for ψ(k+1) + ψ(n+k+1) terms using ratio recurrence to avoid overflow
    let neg_z_half_sq = -z_half * z_half;
    let mut sum1 = Complex64::new(0.0, 0.0);

    // ratio_term tracks (-z²/4)^k / (k! * (n+k)!) via incremental multiplication
    let mut ratio_term = Complex64::new(1.0 / factorial(n as u32), 0.0);

    for k in 0i32..300 {
        let psi_sum = digamma(k as u32 + 1) + digamma((n + k) as u32 + 1);
        sum1 += ratio_term * psi_sum;

        // Update ratio_term for next iteration: multiply by (-z²/4) / ((k+1)*(n+k+1))
        let denom = (k + 1) as f64 * (n + k + 1) as f64;
        ratio_term = ratio_term * neg_z_half_sq / Complex64::new(denom, 0.0);

        if k > 0 && (ratio_term * psi_sum).norm() < 1e-15 * sum1.norm() {
            break;
        }
    }

    let part2 = -(1.0 / PI) * prefix * sum1;

    // Third part: (z/2)^{-n} * sum for k < n
    let mut part3 = Complex64::new(0.0, 0.0);
    if n > 0 {
        let mut z_half_neg_n = Complex64::new(1.0, 0.0);
        for _ in 0..n {
            z_half_neg_n /= z_half;
        }

        let z_half_sq = z_half * z_half;
        let mut sum2 = Complex64::new(0.0, 0.0);
        let mut term2 = Complex64::new(factorial((n - 1) as u32), 0.0); // k=0 term
        for k in 0..n {
            if k > 0 {
                // term2 *= z_half_sq * (n-k) / k  [from ratio of factorials]
                term2 = term2 * z_half_sq / Complex64::new((k * (n - k)) as f64, 0.0);
            }
            sum2 += term2;
        }
        part3 = (1.0 / PI) * z_half_neg_n * sum2;
    }

    part1 + part2 - part3
}

/// H^(1)_n(z) via modified Bessel K_n (full AMOS architecture).
///
/// Uses the identity: H^(1)_n(z) = (2/(πi)) * exp(-inπ/2) * K_n(-iz)
///
/// Two sub-paths based on w = -iz:
/// - Re(w) >= 0 (i.e. Im(z) >= 0): K_n(w) directly (ZBKNU path)
/// - Re(w) < 0  (i.e. Im(z) < 0):  analytic continuation (ZACON path)
///   K_n(w) = (-1)^n * K_n(u) + πi * I_n(u)  where u = -w, Re(u) > 0
fn hankel1_via_k(n: i32, z: Complex64) -> Complex64 {
    let nf = n as f64;
    let w = -Complex64::i() * z; // w = -iz, Re(w) = Im(z)

    let k_n_w = if w.re >= 0.0 {
        // RHP (ZBKNU path): K_n directly
        bessel_k(n, w)
    } else {
        // LHP (ZACON path): analytic continuation of K_n for the H^(1) identity.
        // H^(1)(z) = (2/(πi)) * e^{-inπ/2} * K_n(-iz), where K is analytically
        // continued clockwise (m=-1) from the RHP into the LHP:
        //   K_n(w) = (-1)^n K_n(u) + πi I_n(u),  u = -w (in RHP)
        let u = -w;
        let k_n_u = bessel_k(n, u);
        let i_n_u = bessel_i(n, u);
        let sign_n = if n % 2 == 0 { 1.0 } else { -1.0 };
        Complex64::new(sign_n, 0.0) * k_n_u + Complex64::new(0.0, PI) * i_n_u
    };

    let exp_factor = Complex64::new(0.0, -nf * PI / 2.0).exp();
    let prefactor = Complex64::new(2.0 / PI, 0.0) / Complex64::i();
    prefactor * exp_factor * k_n_w
}

/// Hankel function of the first kind H^(1)_n(z) = J_n(z) + i*Y_n(z)
///
/// For real z (Im=0): J+iY is simpler and exact (no cancellation for real args).
/// For complex z: routes through K_n via the full AMOS architecture (ZBKNU/ZACON).
pub fn hankel1(n: i32, z: Complex64) -> Complex64 {
    if z.im == 0.0 {
        return bessel_j(n, z) + Complex64::i() * bessel_y(n, z);
    }

    hankel1_via_k(n, z)
}

/// Derivative of Bessel function J'_n(z) using recurrence relation:
/// J'_n(z) = (J_{n-1}(z) - J_{n+1}(z)) / 2
pub fn bessel_j_derivative(n: i32, z: Complex64) -> Complex64 {
    (bessel_j(n - 1, z) - bessel_j(n + 1, z)) / 2.0
}

/// Derivative of Hankel function H^(1)'_n(z) using recurrence relation:
/// H'_n(z) = (H_{n-1}(z) - H_{n+1}(z)) / 2
pub fn hankel1_derivative(n: i32, z: Complex64) -> Complex64 {
    (hankel1(n - 1, z) - hankel1(n + 1, z)) / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx_eq(a: Complex64, b: Complex64, tol: f64) -> bool {
        (a - b).norm() < tol
    }

    #[test]
    fn test_j0_real() {
        // J_0(0) ≈ 1
        let j0_0 = bessel_j(0, Complex64::new(1e-15, 0.0));
        assert!(approx_eq(j0_0, Complex64::new(1.0, 0.0), 1e-6));

        // J_0(1) ≈ 0.7651976866
        let j0_1 = bessel_j(0, Complex64::new(1.0, 0.0));
        assert!(approx_eq(j0_1, Complex64::new(0.7651976866, 0.0), 1e-6));

        // J_0(5) ≈ -0.1775967713
        let j0_5 = bessel_j(0, Complex64::new(5.0, 0.0));
        assert!(approx_eq(j0_5, Complex64::new(-0.1775967713, 0.0), 1e-6));
    }

    #[test]
    fn test_j1_real() {
        // J_1(1) ≈ 0.4400505857
        let j1_1 = bessel_j(1, Complex64::new(1.0, 0.0));
        assert!(approx_eq(j1_1, Complex64::new(0.4400505857, 0.0), 1e-6));
    }

    #[test]
    fn test_negative_order() {
        // J_{-n}(z) = (-1)^n J_n(z) for integer n
        let z = Complex64::new(2.5, 0.0);

        // Even order: J_{-2} = J_2
        let j2 = bessel_j(2, z);
        let j_neg2 = bessel_j(-2, z);
        assert!(
            approx_eq(j2, j_neg2, 1e-10),
            "J_2 != J_-2: {} vs {}",
            j2,
            j_neg2
        );

        // Odd order: J_{-3} = -J_3
        let j3 = bessel_j(3, z);
        let j_neg3 = bessel_j(-3, z);
        assert!(
            approx_eq(j3, -j_neg3, 1e-10),
            "J_3 != -J_-3: {} vs {}",
            j3,
            j_neg3
        );
    }
}
