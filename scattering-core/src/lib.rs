pub mod bessel;
pub mod field;
pub mod scattering;
pub mod sources;
use field::{compute_field, FieldParams, DEFAULT_VIEW_SIZE, GRID_SIZE};
use scattering::{
    calculate_scattering, Material, Polarization, ScatteringParams, MAX_ORDER_MAX, MAX_ORDER_MIN,
    PERMEABILITY_IM_MAX, PERMEABILITY_IM_MIN, PERMEABILITY_RE_MAX, PERMEABILITY_RE_MIN,
    PERMITTIVITY_IM_MAX, PERMITTIVITY_IM_MIN, PERMITTIVITY_RE_MAX, PERMITTIVITY_RE_MIN,
    WAVELENGTH_MAX, WAVELENGTH_MIN,
};
use sources::Source;
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

/// Combined scattering + field computation that writes directly into
/// pre-allocated JS Float64Arrays. Zero per-frame JS heap allocations:
/// no return object, no string params, just two Float64Array.set() calls.
///
/// # Arguments
/// * `wavelength` - Wavelength in units where cylinder diameter = 1
/// * `permittivity_real/imag` - Complex relative permittivity
/// * `permeability_real/imag` - Complex relative permeability
/// * `polarization` - 0 = TM, 1 = TE (integer to avoid string encoding)
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
    polarization: u32,
    max_order: i32,
    view_size: f64,
    out_field_real: &js_sys::Float64Array,
    out_field_imag: &js_sys::Float64Array,
) {
    let pol = if polarization == 0 {
        Polarization::TM
    } else {
        Polarization::TE
    };

    let scat_params = ScatteringParams {
        wavelength,
        material: Material {
            permittivity_real,
            permittivity_imag,
            permeability_real,
            permeability_imag,
        },
        polarization: pol,
        max_order,
        source: Source::PlaneWave,
    };
    let scattering = calculate_scattering(&scat_params);

    // Decompose Complex64 → real/imag for FieldParams (stays on WASM heap)
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
        source: Source::PlaneWave,
    };

    let field = compute_field(&field_params);

    // Copy into pre-allocated JS arrays — only JS-side cost is 2 Float64Array.set() calls
    out_field_real.copy_from(&field.field_real);
    out_field_imag.copy_from(&field.field_imag);
}
