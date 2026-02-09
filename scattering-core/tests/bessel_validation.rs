//! Validation tests for Bessel function implementations against SciPy reference data.

use num_complex::Complex64;
use scattering_core::bessel::{bessel_j, hankel1};
use std::collections::BTreeMap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/// A single reference point from SciPy data.
#[derive(Debug, Clone)]
struct ReferencePoint {
    order: i32,
    radius: f64,
    z: Complex64,
    expected: Complex64,
}

/// Statistics for a group of test points.
#[derive(Debug, Clone, Default)]
struct GroupStats {
    count: usize,
    max_absolute_error: f64,
    max_relative_error: f64,
    sum_absolute_error: f64,
    sum_relative_error: f64,
}

impl GroupStats {
    fn mean_relative_error(&self) -> f64 {
        if self.count == 0 {
            0.0
        } else {
            self.sum_relative_error / self.count as f64
        }
    }

    fn update(&mut self, abs_err: f64, rel_err: f64) {
        self.count += 1;
        self.max_absolute_error = self.max_absolute_error.max(abs_err);
        self.max_relative_error = self.max_relative_error.max(rel_err);
        self.sum_absolute_error += abs_err;
        self.sum_relative_error += rel_err;
    }
}

/// Details about a failing or high-error case.
#[derive(Debug, Clone)]
struct FailureDetail {
    order: i32,
    z: Complex64,
    expected: Complex64,
    computed: Complex64,
    absolute_error: f64,
    relative_error: f64,
}

/// Full accuracy report for a Bessel function.
#[derive(Debug)]
struct AccuracyReport {
    function_name: String,
    total_points: usize,
    points_passing: usize,
    max_absolute_error: f64,
    max_relative_error: f64,
    mean_absolute_error: f64,
    mean_relative_error: f64,
    max_rel_error_point: Option<FailureDetail>,
    /// Key: (order, radius_bits) — radius_bits is f64::to_bits() for exact grouping
    errors_by_order_radius: BTreeMap<(i32, u64), GroupStats>,
    worst_cases: Vec<FailureDetail>,
}

/// Parse a line in format: (order, radius, angle, re_result, im_result)
fn parse_line(line: &str) -> Option<ReferencePoint> {
    let line = line.trim();
    if line.is_empty() || line.starts_with('#') {
        return None;
    }

    let inner = line.trim_start_matches('(').trim_end_matches(')');
    let parts: Vec<&str> = inner.split(',').map(|s| s.trim()).collect();

    if parts.len() != 5 {
        return None;
    }

    let order: i32 = parts[0].parse().ok()?;
    let z_radius: f64 = parts[1].parse().ok()?;
    let z_angle: f64 = parts[2].parse().ok()?;
    let re_result: f64 = parts[3].parse().ok()?;
    let im_result: f64 = parts[4].parse().ok()?;

    Some(ReferencePoint {
        order,
        radius: z_radius,
        z: Complex64::new(z_radius * z_angle.cos(), z_radius * z_angle.sin()),
        expected: Complex64::new(re_result, im_result),
    })
}

/// Compute absolute and relative error between computed and expected values.
fn compare(computed: Complex64, expected: Complex64) -> (f64, f64) {
    let abs_err = (computed - expected).norm();
    let expected_norm = expected.norm();

    let rel_err = if expected_norm < 1e-12 {
        if abs_err < 1e-10 {
            0.0
        } else {
            abs_err
        }
    } else {
        abs_err / expected_norm
    };

    (abs_err, rel_err)
}

/// Run validation against a data file using the provided function.
fn run_validation<F>(
    data_path: &str,
    function_name: &str,
    compute_fn: F,
    sample_size: Option<usize>,
) -> AccuracyReport
where
    F: Fn(i32, Complex64) -> Complex64,
{
    let file = File::open(data_path).unwrap_or_else(|_| panic!("Failed to open {}", data_path));
    let reader = BufReader::new(file);

    let rel_err_threshold = 1e-6;

    let mut total_points = 0usize;
    let mut points_passing = 0usize;
    let mut max_absolute_error = 0.0f64;
    let mut max_relative_error = 0.0f64;
    let mut sum_absolute_error = 0.0f64;
    let mut sum_relative_error = 0.0f64;
    let mut max_rel_error_point: Option<FailureDetail> = None;
    let mut errors_by_order_radius: BTreeMap<(i32, u64), GroupStats> = BTreeMap::new();
    let mut worst_cases: Vec<FailureDetail> = Vec::new();

    let all_points: Vec<ReferencePoint> = reader
        .lines()
        .filter_map(|line| line.ok())
        .filter_map(|line| parse_line(&line))
        .collect();

    let points_to_test: Vec<&ReferencePoint> = match sample_size {
        Some(n) if n < all_points.len() => {
            let step = all_points.len() / n;
            all_points.iter().step_by(step).take(n).collect()
        }
        _ => all_points.iter().collect(),
    };

    for point in points_to_test {
        if !point.expected.re.is_finite() || !point.expected.im.is_finite() {
            continue;
        }

        let computed = compute_fn(point.order, point.z);

        if !computed.re.is_finite() || !computed.im.is_finite() {
            total_points += 1;
            let detail = FailureDetail {
                order: point.order,
                z: point.z,
                expected: point.expected,
                computed,
                absolute_error: f64::INFINITY,
                relative_error: f64::INFINITY,
            };
            max_absolute_error = f64::INFINITY;
            max_relative_error = f64::INFINITY;
            max_rel_error_point = Some(detail.clone());
            worst_cases.push(detail);
            continue;
        }

        let (abs_err, rel_err) = compare(computed, point.expected);

        total_points += 1;
        sum_absolute_error += abs_err;
        sum_relative_error += rel_err;

        if rel_err < rel_err_threshold {
            points_passing += 1;
        }

        if rel_err > max_relative_error {
            max_relative_error = rel_err;
            max_rel_error_point = Some(FailureDetail {
                order: point.order,
                z: point.z,
                expected: point.expected,
                computed,
                absolute_error: abs_err,
                relative_error: rel_err,
            });
        }
        if abs_err > max_absolute_error {
            max_absolute_error = abs_err;
        }

        // Track by (order, radius)
        let key = (point.order, point.radius.to_bits());
        errors_by_order_radius
            .entry(key)
            .or_default()
            .update(abs_err, rel_err);

        // Collect candidates for worst-cases list
        if rel_err > rel_err_threshold * 0.1 {
            worst_cases.push(FailureDetail {
                order: point.order,
                z: point.z,
                expected: point.expected,
                computed,
                absolute_error: abs_err,
                relative_error: rel_err,
            });
        }
    }

    worst_cases.sort_by(|a, b| {
        b.relative_error
            .partial_cmp(&a.relative_error)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    worst_cases.truncate(10);

    AccuracyReport {
        function_name: function_name.to_string(),
        total_points,
        points_passing,
        max_absolute_error,
        max_relative_error,
        mean_absolute_error: if total_points > 0 {
            sum_absolute_error / total_points as f64
        } else {
            0.0
        },
        mean_relative_error: if total_points > 0 {
            sum_relative_error / total_points as f64
        } else {
            0.0
        },
        max_rel_error_point,
        errors_by_order_radius,
        worst_cases,
    }
}

/// Collect sorted unique radii from the (order, radius_bits) keys.
fn collect_radii(report: &AccuracyReport) -> Vec<f64> {
    let mut radii_bits: Vec<u64> = report
        .errors_by_order_radius
        .keys()
        .map(|&(_, rb)| rb)
        .collect();
    radii_bits.sort();
    radii_bits.dedup();
    radii_bits.iter().map(|&b| f64::from_bits(b)).collect()
}

/// Collect sorted unique orders from the keys.
fn collect_orders(report: &AccuracyReport) -> Vec<i32> {
    let mut orders: Vec<i32> = report
        .errors_by_order_radius
        .keys()
        .map(|&(o, _)| o)
        .collect();
    orders.sort();
    orders.dedup();
    orders
}

/// Format and print the accuracy report.
fn print_report(report: &AccuracyReport) {
    println!("\n{:=<60}", "");
    println!("  {} Validation Report", report.function_name);
    println!("{:=<60}", "");
    println!("Total points tested: {}", report.total_points);
    println!(
        "Points passing (rel_err < 1e-6): {} ({:.2}%)",
        report.points_passing,
        100.0 * report.points_passing as f64 / report.total_points as f64
    );
    println!("Max relative error: {:.2e}", report.max_relative_error);
    if let Some(ref point) = report.max_rel_error_point {
        println!(
            "  at order={}, z=({:.4}, {:.4}), |z|={:.2}",
            point.order,
            point.z.re,
            point.z.im,
            point.z.norm()
        );
    }
    println!("Mean relative error: {:.2e}", report.mean_relative_error);

    // === Errors by (order, radius) ===
    let radii = collect_radii(report);
    let orders = collect_orders(report);

    println!("\n=== Errors by (Order, Radius) — only showing max_rel > 1e-8 ===");

    for &r in &radii {
        let rb = r.to_bits();

        // Check if any order at this radius has notable error
        let has_notable = orders.iter().any(|&o| {
            report
                .errors_by_order_radius
                .get(&(o, rb))
                .map_or(false, |s| s.max_relative_error > 1e-8)
        });

        if !has_notable {
            // Summarize the whole ring in one line
            let ring_max: f64 = orders
                .iter()
                .filter_map(|&o| report.errors_by_order_radius.get(&(o, rb)))
                .map(|s| s.max_relative_error)
                .fold(0.0f64, f64::max);
            println!(
                "r={:<10.4}  all orders OK (ring max_rel={:.2e})",
                r, ring_max
            );
            continue;
        }

        println!("r={:<10.4}", r);
        for &o in &orders {
            if let Some(stats) = report.errors_by_order_radius.get(&(o, rb)) {
                if stats.max_relative_error > 1e-8 {
                    let flag = if stats.max_relative_error > 1e-6 {
                        "  ** FAIL"
                    } else {
                        ""
                    };
                    println!(
                        "  n={:>2}: max_rel={:.2e}  mean_rel={:.2e}  (n={} pts){}",
                        o,
                        stats.max_relative_error,
                        stats.mean_relative_error(),
                        stats.count,
                        flag
                    );
                }
            }
        }
    }

    // Worst cases
    if !report.worst_cases.is_empty() {
        println!("\n=== Worst {} Cases ===", report.worst_cases.len());
        for (i, case) in report.worst_cases.iter().enumerate() {
            println!(
                "{}. n={}, z=({:.4}, {:.4}), |z|={:.2}: rel_err={:.2e}, abs_err={:.2e}",
                i + 1,
                case.order,
                case.z.re,
                case.z.im,
                case.z.norm(),
                case.relative_error,
                case.absolute_error
            );
            println!(
                "   expected=({:.8e}, {:.8e})",
                case.expected.re, case.expected.im
            );
            println!(
                "   computed=({:.8e}, {:.8e})",
                case.computed.re, case.computed.im
            );
        }
    }

    println!();
}

/// Get the sample size from environment variable, if set.
fn get_sample_size() -> Option<usize> {
    env::var("BESSEL_SAMPLE_SIZE")
        .ok()
        .and_then(|s| s.parse().ok())
}

/// Get the path to the artifacts directory.
fn artifacts_path(filename: &str) -> String {
    let manifest_dir = env!("CARGO_MANIFEST_DIR");
    format!("{}/artifacts/{}", manifest_dir, filename)
}

#[test]
#[ignore]
fn validate_bessel_j_against_scipy() {
    let data_path = artifacts_path("bessel_j_evaluations.txt");
    let sample_size = get_sample_size();

    if sample_size.is_some() {
        println!("Running with sample size: {:?}", sample_size);
    }

    let report = run_validation(&data_path, "Bessel J", bessel_j, sample_size);
    print_report(&report);

    assert!(
        report.max_relative_error < 1e-6,
        "Max relative error {:.2e} exceeds threshold 1e-6",
        report.max_relative_error
    );
}

#[test]
#[ignore]
fn validate_hankel1_against_scipy() {
    let data_path = artifacts_path("hankel1_evaluations.txt");
    let sample_size = get_sample_size();

    if sample_size.is_some() {
        println!("Running with sample size: {:?}", sample_size);
    }

    let report = run_validation(&data_path, "Hankel H1", hankel1, sample_size);
    print_report(&report);

    assert!(
        report.max_relative_error < 1e-6,
        "Max relative error {:.2e} exceeds threshold 1e-6",
        report.max_relative_error
    );
}

#[test]
#[ignore] // Run with: cargo test validate_bessel_detailed_report -- --ignored --nocapture
fn validate_bessel_detailed_report() {
    let bessel_j_path = artifacts_path("bessel_j_evaluations.txt");
    let hankel1_path = artifacts_path("hankel1_evaluations.txt");

    let bessel_j_report = run_validation(&bessel_j_path, "Bessel J", bessel_j, None);
    let hankel1_report = run_validation(&hankel1_path, "Hankel H1", hankel1, None);

    print_report(&bessel_j_report);
    print_report(&hankel1_report);

    // Write summary to file
    let report_path = artifacts_path("validation_report.txt");
    let mut file = File::create(&report_path).expect("Failed to create report file");

    writeln!(file, "Bessel Function Validation Report").unwrap();
    writeln!(file, "==================================\n").unwrap();

    for report in [&bessel_j_report, &hankel1_report] {
        writeln!(file, "Function: {}", report.function_name).unwrap();
        writeln!(file, "Total points: {}", report.total_points).unwrap();
        writeln!(
            file,
            "Points passing: {} ({:.3}%)",
            report.points_passing,
            100.0 * report.points_passing as f64 / report.total_points as f64
        )
        .unwrap();
        writeln!(
            file,
            "Max relative error: {:.6e}",
            report.max_relative_error
        )
        .unwrap();
        writeln!(
            file,
            "Mean relative error: {:.6e}",
            report.mean_relative_error
        )
        .unwrap();
        writeln!(file).unwrap();
    }

    println!("Detailed report written to: {}", report_path);
}
