//! Complex Bessel functions for electromagnetic scattering calculations.
//!
//! Pure Rust implementation using series expansions and asymptotic formulas.

use num_complex::Complex64;
use std::f64::consts::PI;

/// Euler-Mascheroni constant
const EULER_GAMMA: f64 = 0.5772156649015329;

/// Threshold for switching from power series to asymptotic expansion.
/// For Hankel H1 (single series, no cancellation), this threshold works directly.
/// For Bessel J/Y (computed as H1±H2, cancellation for complex z), the effective
/// threshold is increased by the order to ensure accuracy.
const ASYMPTOTIC_THRESHOLD: f64 = 25.0;

/// Factorial function for small integers
fn factorial(n: u32) -> u64 {
    match n {
        0 | 1 => 1,
        _ => (2..=n as u64).product(),
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
    let mut term = Complex64::new(1.0 / factorial(n as u32) as f64, 0.0);
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
    // Requires |z| > 2n+5 so the asymptotic series has enough converging terms.
    if z.im.abs() >= 5.0 && znorm >= 2.0 * nf + 5.0 {
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

    // Use asymptotic expansion for large |z|
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
    let mut ratio_term = Complex64::new(1.0 / factorial(n as u32) as f64, 0.0);

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
        let mut term2 = Complex64::new(factorial((n - 1) as u32) as f64, 0.0); // k=0 term
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

/// Hankel function of the first kind H^(1)_n(z) = J_n(z) + i*Y_n(z)
pub fn hankel1(n: i32, z: Complex64) -> Complex64 {
    // H1 asymptotic computes H1 directly (single series, no cancellation).
    //
    // The J+iY fallback suffers catastrophic cancellation for complex z when
    // |Im(z)| is large (J and Y each ~ exp(|Im(z)|) while H1 ~ exp(-|Im(z)|)).
    // Cancellation loses ~2|Im(z)|/ln(10) digits. With 15 digits of f64 precision,
    // this is catastrophic when |Im(z)| > ~17.
    //
    // Use asymptotic when either:
    // (a) |z| is large enough for good asymptotic accuracy, OR
    // (b) |Im(z)| is large enough that J+iY cancellation would be worse
    let znorm = z.norm();
    let abs_n = n.unsigned_abs() as f64;
    let asymptotic_threshold = (abs_n * abs_n / 2.0 + abs_n).max(ASYMPTOTIC_THRESHOLD);
    let jiy_cancellation_bad = z.im.abs() >= 5.0 && znorm >= abs_n + 3.0;
    if znorm > asymptotic_threshold || jiy_cancellation_bad {
        return hankel1_asymptotic(n, z);
    }
    bessel_j(n, z) + Complex64::i() * bessel_y(n, z)
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
    fn diagnose_asymptotic_bug() {
        // Compare asymptotic at z with positive vs negative real part
        let cases: Vec<(&str, i32, Complex64, Complex64)> = vec![
            // z = +30 + 5i (first quadrant - far from branch cut)
            (
                "Q1 +30+5i",
                0,
                Complex64::new(30.0, 5.0),
                bessel_j_miller(0, Complex64::new(30.0, 5.0)),
            ),
            // z = -30 + 5i (second quadrant - near negative real axis)
            (
                "Q2 -30+5i",
                0,
                Complex64::new(-30.0, 5.0),
                bessel_j_miller(0, Complex64::new(-30.0, 5.0)),
            ),
            // z = 30 (pure real positive)
            (
                "Real +30",
                0,
                Complex64::new(30.0, 0.0),
                bessel_j_miller(0, Complex64::new(30.0, 0.0)),
            ),
            // z = -30 (pure real negative)
            (
                "Real -30",
                0,
                Complex64::new(-30.0, 0.0),
                bessel_j_miller(0, Complex64::new(-30.0, 0.0)),
            ),
        ];

        println!("\n=== Asymptotic vs Miller for n=0 at various z ===");
        for (label, n, z, ref_val) in &cases {
            let j_asym = bessel_j_asymptotic(*n, *z);
            let h1 = hankel1_asymptotic(*n, *z);
            let h2 = hankel2_asymptotic(*n, *z);
            let err = if ref_val.norm() < 1e-15 {
                (j_asym - ref_val).norm()
            } else {
                (j_asym - ref_val).norm() / ref_val.norm()
            };
            println!("{:14}: z=({:>6.1},{:>5.1}) asym_err={:.2e}  |H1|={:.4e} |H2|={:.4e} |(H1+H2)/2|={:.4e} |ref|={:.4e}",
                label, z.re, z.im, err, h1.norm(), h2.norm(), ((h1+h2)/2.0).norm(), ref_val.norm());
        }
    }

    #[test]
    fn diagnose_algo_comparison() {
        // n=0 at the problematic angle
        let r = 30.39195382313201_f64;
        let theta = 2.9762456718219092_f64;
        let z = Complex64::new(r * theta.cos(), r * theta.sin());
        let expected = Complex64::new(-5.8537694126084805, -9.042308585991176);

        let j_asym = bessel_j_asymptotic(0, z);
        let j_mill = bessel_j_miller(0, z);
        let j_ser = bessel_j_series(0, z);
        let h1 = hankel1_asymptotic(0, z);
        let h2 = hankel2_asymptotic(0, z);

        println!("\nn=0, z=({:.4},{:.4}), |z|={:.4}", z.re, z.im, z.norm());
        println!("expected:  ({:.15e}, {:.15e})", expected.re, expected.im);
        println!(
            "asym(h1h2):({:.15e}, {:.15e}) err={:.2e}",
            j_asym.re,
            j_asym.im,
            (j_asym - expected).norm() / expected.norm()
        );
        println!(
            "miller:    ({:.15e}, {:.15e}) err={:.2e}",
            j_mill.re,
            j_mill.im,
            (j_mill - expected).norm() / expected.norm()
        );
        println!(
            "series:    ({:.15e}, {:.15e}) err={:.2e}",
            j_ser.re,
            j_ser.im,
            (j_ser - expected).norm() / expected.norm()
        );
        println!(
            "|H1|={:.6e}, |H2|={:.6e}, |H2/H1|={:.1}",
            h1.norm(),
            h2.norm(),
            h2.norm() / h1.norm()
        );

        // Also test n=0 at r=6.21
        let r2 = 6.2101694189156165_f64;
        let theta2 = 0.3306939635357677_f64;
        let z2 = Complex64::new(r2 * theta2.cos(), r2 * theta2.sin());
        let exp2 = Complex64::new(0.6060728635027623, 1.0272445013386244);
        let j2_asym = bessel_j_asymptotic(0, z2);
        let j2_mill = bessel_j_miller(0, z2);
        let j2_ser = bessel_j_series(0, z2);
        println!("\nn=0, z=({:.4},{:.4}), |z|={:.4}", z2.re, z2.im, z2.norm());
        println!("expected:  ({:.15e}, {:.15e})", exp2.re, exp2.im);
        println!(
            "asym(h1h2):({:.15e}, {:.15e}) err={:.2e}",
            j2_asym.re,
            j2_asym.im,
            (j2_asym - exp2).norm() / exp2.norm()
        );
        println!(
            "miller:    ({:.15e}, {:.15e}) err={:.2e}",
            j2_mill.re,
            j2_mill.im,
            (j2_mill - exp2).norm() / exp2.norm()
        );
        println!(
            "series:    ({:.15e}, {:.15e}) err={:.2e}",
            j2_ser.re,
            j2_ser.im,
            (j2_ser - exp2).norm() / exp2.norm()
        );
    }

    /// Diagnose which angles cause errors at a given radius and order.
    #[test]
    fn diagnose_per_angle_errors() {
        let r = 30.39195382313201_f64;
        let angles_and_expected: Vec<(f64, (f64, f64))> = vec![
            (0.0, (-0.03477788848272825, 3.469446951953614e-18)),
            (
                0.3306939635357677,
                (-1382.6275021454458, -213.3234763703907),
            ),
            (0.6613879270715354, (-338160.5788035898, 9275993.117730552)),
            (0.992081890607303, (-6613236194.86076, 4753157673.184689)),
            (1.3227758541430708, (224585675204.88672, -393842137797.5656)),
            (1.6534698176788385, (-809681933409.5487, 646143275259.1244)),
            (1.984163781214606, (74950974939.76328, -47666655870.90963)),
            (2.3148577447503738, (78239842.2365667, 364702230.0513121)),
            (2.6455517082861415, (68448.89656846572, 120760.0286668394)),
            (
                2.9762456718219092,
                (-5.8537694126084805, -9.042308585991176),
            ),
            (3.306939635357677, (-5.853769412608403, 9.042308585991124)),
            (3.6376335988934447, (68448.89656846496, -120760.02866683899)),
            (3.968327562429212, (78239842.23656052, -364702230.051308)),
            (4.29902152596498, (74950974939.76364, 47666655870.90786)),
            (4.6297154895007475, (-809681933409.5579, -646143275259.1133)),
            (4.960409453036515, (224585675204.8931, 393842137797.564)),
            (5.291103416572283, (-6613236194.860839, -4753157673.184674)),
            (5.621797380108051, (-338160.5788036227, -9275993.117730618)),
            (5.9524913436438185, (-1382.627502145458, 213.32347637039246)),
            (
                6.283185307179586,
                (-0.03477788848272825, -1.052977149917922e-15),
            ),
        ];

        println!("\n=== n=0, r={:.4} per-angle diagnosis ===", r);
        for (theta, (re_exp, im_exp)) in &angles_and_expected {
            let z = Complex64::new(r * theta.cos(), r * theta.sin());
            let expected = Complex64::new(*re_exp, *im_exp);
            let computed = bessel_j(0, z);
            let rel_err = if expected.norm() < 1e-12 {
                (computed - expected).norm()
            } else {
                (computed - expected).norm() / expected.norm()
            };

            // Determine which algorithm is used
            let znorm = z.norm();
            let cancel_param = z.re * z.re - z.im * z.im;
            let algo = if z.im.abs() >= 5.0 && znorm >= 5.0 {
                "asym(dom)"
            } else if znorm > 25.0 {
                "asym(h1h2)"
            } else if cancel_param < 196.0 {
                "series"
            } else {
                "miller"
            };

            if rel_err > 1e-10 {
                println!(
                    "  θ={:.3} z=({:>8.2},{:>8.2}) |Im|={:.2} algo={:<12} rel_err={:.2e}",
                    theta,
                    z.re,
                    z.im,
                    z.im.abs(),
                    algo,
                    rel_err
                );
            }
        }
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
