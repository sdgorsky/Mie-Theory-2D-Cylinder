//! End-to-end benchmark: scattering coefficients + field computation.
//!
//! Sweeps grid edge sizes [2, 4, 8, ..., 1024] × max orders [5, 10, 20]
//! and prints a formatted report of timing results.

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use scattering_core::field::{compute_field, FieldParams};
use scattering_core::scattering::{calculate_scattering, Material, ScatteringParams, RADIUS};
use scattering_core::sources::{Source, SourceKind};
use std::fmt::Write as FmtWrite;
use std::fs;
use std::time::Instant;

/// Build FieldParams for a given grid_size and max_order, running the
/// scattering coefficient solve as part of the setup.
fn build_field_params(grid_size: usize, max_order: i32) -> FieldParams {
    let wavelength = 1.0;
    let material = Material {
        permittivity_real: 4.0,
        permittivity_imag: 0.0,
        permeability_real: 1.0,
        permeability_imag: 0.0,
    };
    let source = Source::new(SourceKind::PlaneWaveTM, RADIUS);

    let scat_params = ScatteringParams {
        wavelength,
        material,
        max_order,
        source,
    };
    let scattering = calculate_scattering(&scat_params);

    FieldParams {
        wavelength,
        permittivity_real: material.permittivity_real,
        permittivity_imag: material.permittivity_imag,
        permeability_real: material.permeability_real,
        permeability_imag: material.permeability_imag,
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
        view_size: 5.0,
        grid_size,
        source,
    }
}

fn bench_field(c: &mut Criterion) {
    let edge_sizes: Vec<usize> = (1..=10).map(|e| 1 << e).collect(); // 2, 4, 8, ..., 1024
    let max_orders = [5i32, 10, 20];

    for &order in &max_orders {
        let mut group = c.benchmark_group(format!("max_order_{order}"));
        for &edge in &edge_sizes {
            let params = build_field_params(edge, order);
            group.bench_with_input(BenchmarkId::new("compute_field", edge), &params, |b, p| {
                b.iter(|| compute_field(p))
            });
        }
        group.finish();
    }
}

/// Run a single-pass timing sweep and print a formatted report.
/// This is NOT a criterion benchmark — it runs once per combo for a quick overview.
fn print_report() {
    let edge_sizes: Vec<usize> = (1..=10).map(|e| 1 << e).collect();
    let max_orders = [5i32, 10, 20];
    // Warm-up iterations for small grids so timing is stable
    let warmup_iters = 3;

    let mut report = String::new();

    writeln!(report, "{:=<70}", "").unwrap();
    writeln!(report, "  End-to-end field benchmark report").unwrap();
    writeln!(report, "{:=<70}", "").unwrap();
    writeln!(
        report,
        "{:>8} {:>8} {:>14} {:>16}",
        "order", "edge", "total (ms)", "per point (µs)"
    )
    .unwrap();
    writeln!(report, "{:-<70}", "").unwrap();

    for &order in &max_orders {
        for &edge in &edge_sizes {
            let params = build_field_params(edge, order);
            let n_points = edge * edge;

            // Warm up
            for _ in 0..warmup_iters {
                let _ = compute_field(&params);
            }

            // Timed run (average of several iterations for small grids)
            let iters = if n_points < 10_000 { 10 } else { 3 };
            let start = Instant::now();
            for _ in 0..iters {
                let _ = compute_field(&params);
            }
            let elapsed = start.elapsed().as_secs_f64() / iters as f64;

            let total_ms = elapsed * 1e3;
            let per_point_us = elapsed * 1e6 / n_points as f64;

            writeln!(
                report,
                "{:>8} {:>8} {:>14.4} {:>16.4}",
                order, edge, total_ms, per_point_us
            )
            .unwrap();
        }
        writeln!(report, "{:-<70}", "").unwrap();
    }

    // Print to stdout
    print!("{report}");

    // Write to file
    let out_path = "benches/field_benchmark_report.txt";
    fs::write(out_path, &report).expect("failed to write benchmark report");
    println!("Report written to {out_path}");
}

fn bench_with_report(c: &mut Criterion) {
    bench_field(c);
    print_report();
}

criterion_group!(benches, bench_with_report);
criterion_main!(benches);
