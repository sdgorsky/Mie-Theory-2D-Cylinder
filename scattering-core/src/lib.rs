pub mod bessel;
pub mod field;
pub mod scattering;
pub mod sources;
use field::{compute_field, FieldParams, ModalSplineCache, DEFAULT_VIEW_SIZE, GRID_SIZE};
use scattering::{
    calculate_scattering, Material, ScatteringParams, MAX_ORDER_MAX, MAX_ORDER_MIN,
    PERMEABILITY_IM_MAX, PERMEABILITY_IM_MIN, PERMEABILITY_RE_MAX, PERMEABILITY_RE_MIN,
    PERMITTIVITY_IM_MAX, PERMITTIVITY_IM_MIN, PERMITTIVITY_RE_MAX, PERMITTIVITY_RE_MIN, RADIUS,
    WAVELENGTH_MAX, WAVELENGTH_MIN,
};
use sources::{Source, SourceKind};
use wasm_bindgen::prelude::*;

/// Initialize panic hook for better error messages in browser console.
#[wasm_bindgen(start)]
pub fn init() {
    #[cfg(feature = "console_error_panic_hook")]
    console_error_panic_hook::set_once();
}

/// Get the grid size used for field computation.
#[wasm_bindgen]
pub fn get_field_grid_size() -> usize {
    GRID_SIZE
}

/// Get the default view size (in cylinder diameters) for field computation.
#[wasm_bindgen]
pub fn get_field_view_size() -> f64 {
    DEFAULT_VIEW_SIZE
}

/// Get parameter bounds (single source of truth for UI sliders).
#[wasm_bindgen]
pub fn get_parameter_bounds() -> JsValue {
    let obj = js_sys::Object::new();
    let _ = js_sys::Reflect::set(
        &obj,
        &"wavelength_min".into(),
        &JsValue::from(WAVELENGTH_MIN),
    );
    let _ = js_sys::Reflect::set(
        &obj,
        &"wavelength_max".into(),
        &JsValue::from(WAVELENGTH_MAX),
    );
    let _ = js_sys::Reflect::set(
        &obj,
        &"permittivity_re_min".into(),
        &JsValue::from(PERMITTIVITY_RE_MIN),
    );
    let _ = js_sys::Reflect::set(
        &obj,
        &"permittivity_re_max".into(),
        &JsValue::from(PERMITTIVITY_RE_MAX),
    );
    let _ = js_sys::Reflect::set(
        &obj,
        &"permittivity_im_min".into(),
        &JsValue::from(PERMITTIVITY_IM_MIN),
    );
    let _ = js_sys::Reflect::set(
        &obj,
        &"permittivity_im_max".into(),
        &JsValue::from(PERMITTIVITY_IM_MAX),
    );
    let _ = js_sys::Reflect::set(
        &obj,
        &"permeability_re_min".into(),
        &JsValue::from(PERMEABILITY_RE_MIN),
    );
    let _ = js_sys::Reflect::set(
        &obj,
        &"permeability_re_max".into(),
        &JsValue::from(PERMEABILITY_RE_MAX),
    );
    let _ = js_sys::Reflect::set(
        &obj,
        &"permeability_im_min".into(),
        &JsValue::from(PERMEABILITY_IM_MIN),
    );
    let _ = js_sys::Reflect::set(
        &obj,
        &"permeability_im_max".into(),
        &JsValue::from(PERMEABILITY_IM_MAX),
    );
    let _ = js_sys::Reflect::set(&obj, &"max_order_min".into(), &JsValue::from(MAX_ORDER_MIN));
    let _ = js_sys::Reflect::set(&obj, &"max_order_max".into(), &JsValue::from(MAX_ORDER_MAX));
    obj.into()
}

/// Full pipeline: scattering coefficients + field computation.
/// Callable from native code (benchmarks, tests) without JS dependencies.
#[allow(clippy::too_many_arguments)]
pub fn compute_all_native(
    wavelength: f64,
    permittivity_real: f64,
    permittivity_imag: f64,
    permeability_real: f64,
    permeability_imag: f64,
    source: Source,
    max_order: i32,
    view_size: f64,
    grid_size: usize,
) -> field::FieldResult {
    compute_all_native_cached(
        wavelength,
        permittivity_real,
        permittivity_imag,
        permeability_real,
        permeability_imag,
        source,
        max_order,
        view_size,
        grid_size,
        &mut None,
    )
}

/// Full pipeline with spline cache for reuse across frames.
#[allow(clippy::too_many_arguments)]
pub fn compute_all_native_cached(
    wavelength: f64,
    permittivity_real: f64,
    permittivity_imag: f64,
    permeability_real: f64,
    permeability_imag: f64,
    source: Source,
    max_order: i32,
    view_size: f64,
    grid_size: usize,
    spline_cache: &mut Option<ModalSplineCache>,
) -> field::FieldResult {
    let scat_params = ScatteringParams {
        wavelength,
        material: Material {
            permittivity_real,
            permittivity_imag,
            permeability_real,
            permeability_imag,
        },
        max_order,
        source,
    };
    let scattering = calculate_scattering(&scat_params);

    let n_coeffs = scattering.orders.len();
    let mut scat_re = Vec::with_capacity(n_coeffs);
    let mut scat_im = Vec::with_capacity(n_coeffs);
    let mut int_re = Vec::with_capacity(n_coeffs);
    let mut int_im = Vec::with_capacity(n_coeffs);

    for i in 0..n_coeffs {
        scat_re.push(scattering.scattering_coefficients[i].re);
        scat_im.push(scattering.scattering_coefficients[i].im);
        int_re.push(scattering.internal_coefficients[i].re);
        int_im.push(scattering.internal_coefficients[i].im);
    }

    let field_params = FieldParams {
        wavelength,
        permittivity_real,
        permittivity_imag,
        permeability_real,
        permeability_imag,
        scattering_coeffs_real: scat_re,
        scattering_coeffs_imag: scat_im,
        internal_coeffs_real: int_re,
        internal_coeffs_imag: int_im,
        orders: scattering.orders,
        view_size,
        grid_size,
        source,
    };

    compute_field(&field_params, spline_cache)
}

/// WASM wrapper around `compute_all_native` that writes into pre-allocated
/// JS Float64Arrays. Zero per-frame JS heap allocations.
///
/// # Arguments
/// * `wavelength` - Wavelength in units where cylinder diameter = 1
/// * `permittivity_real/imag` - Complex relative permittivity
/// * `permeability_real/imag` - Complex relative permeability
/// * `source_type` - 0 = PlaneWave TM, 1 = PlaneWave TE, 2 = DipoleEz, 3 = DipoleExy
/// * `dipole_xs/ys` - Dipole position (ignored for plane wave)
/// * `dipole_alpha` - Dipole orientation in radians (ignored except for DipoleExy)
/// * `max_order` - Maximum Bessel order N (computes -N to +N)
/// * `out_field_real` - Pre-allocated Float64Array for real part of field
/// * `out_field_imag` - Pre-allocated Float64Array for imaginary part of field
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn compute_all(
    wavelength: f64,
    permittivity_real: f64,
    permittivity_imag: f64,
    permeability_real: f64,
    permeability_imag: f64,
    source_type: u32,
    dipole_xs: f64,
    dipole_ys: f64,
    dipole_alpha: f64,
    max_order: i32,
    view_size: f64,
    out_field_real: &js_sys::Float64Array,
    out_field_imag: &js_sys::Float64Array,
) {
    let source_kind = match source_type {
        0 => SourceKind::PlaneWaveTM,
        1 => SourceKind::PlaneWaveTE,
        2 => SourceKind::DipoleEz {
            xs: dipole_xs,
            ys: dipole_ys,
        },
        3 => SourceKind::DipoleExy {
            xs: dipole_xs,
            ys: dipole_ys,
            alpha: dipole_alpha,
        },
        _ => SourceKind::PlaneWaveTM,
    };
    let source = Source::new(source_kind, RADIUS);

    // Persistent spline cache — survives across WASM calls.
    // Modal splines are only rebuilt when material/wavelength/order change,
    // not on every dipole position update.
    thread_local! {
        static SPLINE_CACHE: std::cell::RefCell<Option<ModalSplineCache>> =
            const { std::cell::RefCell::new(None) };
    }

    let field = SPLINE_CACHE.with(|cell| {
        let mut cache = cell.borrow_mut();
        compute_all_native_cached(
            wavelength,
            permittivity_real,
            permittivity_imag,
            permeability_real,
            permeability_imag,
            source,
            max_order,
            view_size,
            GRID_SIZE,
            &mut cache,
        )
    });

    out_field_real.copy_from(&field.field_real);
    out_field_imag.copy_from(&field.field_imag);
}
