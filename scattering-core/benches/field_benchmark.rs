use criterion::{criterion_group, criterion_main, Criterion};
use scattering_core::field::{compute_field, FieldParams};
use scattering_core::scattering::{calculate_scattering, Material, Polarization, ScatteringParams};

fn calculate_field() {
    let params = ScatteringParams {
        wavelength: 1.0,
        material: Material {
            permittivity_real: 4.0,
            permittivity_imag: 0.0,
            permeability_real: 1.0,
            permeability_imag: 0.0,
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
    };

    let _result = compute_field(&field_params);
}

fn criterion_benchmark(c: &mut Criterion) {
    let mut group = c.benchmark_group("Field Benchmark Group");

    group.sample_size(10);
    group.bench_function("calculate_field", |b| b.iter(|| calculate_field()));

    group.finish();
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
