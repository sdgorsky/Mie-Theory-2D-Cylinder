//! End-to-end benchmark: scattering coefficients + field computation.
//!
//! Sweeps grid edge sizes [2, 4, 8, ..., 1024] × max orders [5, 10, 20]
//! and prints a formatted report of timing results.

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use scattering_core::field::{compute_field, FieldParams, ModalSplineCache, GRID_SIZE};
use scattering_core::scattering::{calculate_scattering, Material, ScatteringParams, RADIUS};
use scattering_core::sources::{Source, SourceKind};
use scattering_core::{compute_all_native, compute_all_native_cached};
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
                let mut cache = None;
                b.iter(|| compute_field(p, &mut cache))
            });
        }
        group.finish();
    }
}

fn build_scattering_params(max_order: i32) -> ScatteringParams {
    ScatteringParams {
        wavelength: 1.0,
        material: Material {
            permittivity_real: 4.0,
            permittivity_imag: 0.0,
            permeability_real: 1.0,
            permeability_imag: 0.0,
        },
        max_order,
        source: Source::new(SourceKind::PlaneWaveTM, RADIUS),
    }
}

/// Run a single-pass timing sweep and print a formatted report.
/// This is NOT a criterion benchmark — it runs once per combo for a quick overview.
fn print_report() {
    let edge_sizes: Vec<usize> = (1..=10).map(|e| 1 << e).collect();
    let max_orders = [5i32, 10, 20];
    let warmup_iters = 3;

    let mut report = String::new();

    // --- Section 1: calculate_scattering cost ---
    writeln!(report, "{:=<70}", "").unwrap();
    writeln!(report, "  calculate_scattering timing").unwrap();
    writeln!(report, "{:=<70}", "").unwrap();
    writeln!(report, "{:>8} {:>14}", "order", "time (ms)").unwrap();
    writeln!(report, "{:-<70}", "").unwrap();

    for &order in &max_orders {
        let scat_params = build_scattering_params(order);

        // Warm up
        for _ in 0..warmup_iters {
            let _ = calculate_scattering(&scat_params);
        }

        let iters = 10;
        let start = Instant::now();
        for _ in 0..iters {
            let _ = calculate_scattering(&scat_params);
        }
        let elapsed_ms = start.elapsed().as_secs_f64() * 1e3 / iters as f64;

        writeln!(report, "{:>8} {:>14.4}", order, elapsed_ms).unwrap();
    }
    writeln!(report).unwrap();

    // --- Section 2: compute_all_native by source type at production grid size ---
    let source_configs: Vec<(&str, Source)> = vec![
        ("PlaneWaveTM", Source::new(SourceKind::PlaneWaveTM, RADIUS)),
        (
            "DipoleEz ext",
            Source::new(SourceKind::DipoleEz { xs: 1.7, ys: -0.9 }, RADIUS),
        ),
        (
            "DipoleEz int",
            Source::new(SourceKind::DipoleEz { xs: 0.15, ys: -0.1 }, RADIUS),
        ),
        (
            "DipoleExy ext",
            Source::new(
                SourceKind::DipoleExy {
                    xs: -1.3,
                    ys: 0.6,
                    alpha: 0.73,
                },
                RADIUS,
            ),
        ),
        (
            "DipoleExy int",
            Source::new(
                SourceKind::DipoleExy {
                    xs: -0.12,
                    ys: 0.2,
                    alpha: -1.1,
                },
                RADIUS,
            ),
        ),
    ];

    writeln!(report, "{:=<70}", "").unwrap();
    writeln!(
        report,
        "  compute_all_native (scattering + field) at {GRID_SIZE}x{GRID_SIZE}"
    )
    .unwrap();
    writeln!(report, "{:=<70}", "").unwrap();
    writeln!(
        report,
        "{:>8} {:>14} {:>8} {:>14}",
        "material", "source", "order", "time (ms)"
    )
    .unwrap();
    writeln!(report, "{:-<70}", "").unwrap();

    // Material configs: (label, eps_re, eps_im, mu_re, mu_im)
    let material_configs: Vec<(&str, f64, f64, f64, f64)> = vec![
        ("real", 4.0, 0.0, 1.0, 0.0),
        ("complex", 4.0, 0.5, 1.0, -0.2),
    ];

    for (mat_label, eps_re, eps_im, mu_re, mu_im) in &material_configs {
        for (src_label, source) in &source_configs {
            for &order in &max_orders {
                // Use a persistent cache — warm it up, then time with cache hits
                let mut cache = None;
                for _ in 0..warmup_iters {
                    let _ = compute_all_native_cached(
                        1.0, *eps_re, *eps_im, *mu_re, *mu_im, *source, order, 5.0, GRID_SIZE,
                        &mut cache,
                    );
                }

                let iters = 3;
                let start = Instant::now();
                for _ in 0..iters {
                    let _ = compute_all_native_cached(
                        1.0, *eps_re, *eps_im, *mu_re, *mu_im, *source, order, 5.0, GRID_SIZE,
                        &mut cache,
                    );
                }
                let elapsed_ms = start.elapsed().as_secs_f64() * 1e3 / iters as f64;

                writeln!(
                    report,
                    "{:>8} {:>14} {:>8} {:>14.4}",
                    mat_label, src_label, order, elapsed_ms
                )
                .unwrap();
            }
        }
        writeln!(report, "{:-<70}", "").unwrap();
    }
    writeln!(report).unwrap();

    // --- Section 3: compute_field cost ---
    writeln!(report, "{:=<70}", "").unwrap();
    writeln!(report, "  compute_field timing").unwrap();
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

            let mut cache = None;
            for _ in 0..warmup_iters {
                let _ = compute_field(&params, &mut cache);
            }

            let iters = if n_points < 10_000 { 10 } else { 3 };
            cache = None; // fresh cache for timing
            let start = Instant::now();
            for _ in 0..iters {
                let _ = compute_field(&params, &mut cache);
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
