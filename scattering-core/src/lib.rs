pub mod bessel;
pub mod field;
pub mod scattering;
mod utils;

use field::{compute_field, FieldParams, GRID_SIZE, VIEW_SIZE};
use scattering::{calculate_scattering, Material, Polarization, ScatteringParams};
use wasm_bindgen::prelude::*;

/// Initialize panic hook for better error messages in browser console.
#[wasm_bindgen(start)]
pub fn init() {
    utils::set_panic_hook();
}

/// Compute scattering coefficients for a 2D cylinder.
///
/// # Arguments
/// * `wavelength` - Wavelength in units where cylinder diameter = 1
/// * `permittivity_real` - Real part of relative permittivity
/// * `permittivity_imag` - Imaginary part of relative permittivity
/// * `permeability_real` - Real part of relative permeability
/// * `permeability_imag` - Imaginary part of relative permeability
/// * `polarization` - "TM" or "TE"
/// * `max_order` - Maximum Bessel order N (computes -N to +N)
///
/// # Returns
/// A JsValue containing the ScatteringResult object
#[wasm_bindgen]
pub fn compute_scattering(
    wavelength: f64,
    permittivity_real: f64,
    permittivity_imag: f64,
    permeability_real: f64,
    permeability_imag: f64,
    polarization: &str,
    max_order: i32,
) -> Result<JsValue, JsValue> {
    let pol = match polarization.to_uppercase().as_str() {
        "TM" => Polarization::TM,
        "TE" => Polarization::TE,
        _ => return Err(JsValue::from_str("Invalid polarization. Use 'TM' or 'TE'.")),
    };

    let params = ScatteringParams {
        wavelength,
        material: Material {
            permittivity_real,
            permittivity_imag,
            permeability_real,
            permeability_imag,
        },
        polarization: pol,
        max_order,
    };

    let result = calculate_scattering(&params);

    serde_wasm_bindgen::to_value(&result).map_err(|e| JsValue::from_str(&e.to_string()))
}

/// Compute the electric field on a grid (currently 128x128).
///
/// # Arguments
/// * `wavelength` - Wavelength in units where cylinder diameter = 1
/// * `permittivity_real` - Real part of relative permittivity
/// * `permittivity_imag` - Imaginary part of relative permittivity
/// * `permeability_real` - Real part of relative permeability
/// * `permeability_imag` - Imaginary part of relative permeability
/// * `incident_coeffs_real` - Real parts of incident coefficients (a_n)
/// * `incident_coeffs_imag` - Imaginary parts of incident coefficients
/// * `scattering_coeffs_real` - Real parts of scattering coefficients (b_n)
/// * `scattering_coeffs_imag` - Imaginary parts of scattering coefficients
/// * `internal_coeffs_real` - Real parts of internal coefficients (c_n)
/// * `internal_coeffs_imag` - Imaginary parts of internal coefficients
/// * `orders` - Array of orders (-N to +N)
///
/// # Returns
/// A JsValue containing the FieldResult object with field_real and field_imag arrays
#[wasm_bindgen]
#[allow(clippy::too_many_arguments)]
pub fn compute_electric_field(
    wavelength: f64,
    permittivity_real: f64,
    permittivity_imag: f64,
    permeability_real: f64,
    permeability_imag: f64,
    incident_coeffs_real: Vec<f64>,
    incident_coeffs_imag: Vec<f64>,
    scattering_coeffs_real: Vec<f64>,
    scattering_coeffs_imag: Vec<f64>,
    internal_coeffs_real: Vec<f64>,
    internal_coeffs_imag: Vec<f64>,
    orders: Vec<i32>,
) -> Result<JsValue, JsValue> {
    let params = FieldParams {
        wavelength,
        permittivity_real,
        permittivity_imag,
        permeability_real,
        permeability_imag,
        incident_coeffs_real,
        incident_coeffs_imag,
        scattering_coeffs_real,
        scattering_coeffs_imag,
        internal_coeffs_real,
        internal_coeffs_imag,
        orders,
    };

    let result = compute_field(&params);

    // Build JS object manually with Float64Array bulk copies instead of serde
    let obj = js_sys::Object::new();
    let field_real = js_sys::Float64Array::from(result.field_real.as_slice());
    let field_imag = js_sys::Float64Array::from(result.field_imag.as_slice());
    js_sys::Reflect::set(&obj, &"field_real".into(), &field_real)?;
    js_sys::Reflect::set(&obj, &"field_imag".into(), &field_imag)?;
    js_sys::Reflect::set(
        &obj,
        &"grid_size".into(),
        &JsValue::from(result.grid_size as u32),
    )?;
    js_sys::Reflect::set(&obj, &"view_size".into(), &JsValue::from(result.view_size))?;
    js_sys::Reflect::set(&obj, &"x_min".into(), &JsValue::from(result.x_min))?;
    js_sys::Reflect::set(&obj, &"x_max".into(), &JsValue::from(result.x_max))?;
    js_sys::Reflect::set(&obj, &"y_min".into(), &JsValue::from(result.y_min))?;
    js_sys::Reflect::set(&obj, &"y_max".into(), &JsValue::from(result.y_max))?;

    Ok(obj.into())
}

/// Get the grid size used for field computation.
#[wasm_bindgen]
pub fn get_field_grid_size() -> usize {
    GRID_SIZE
}

/// Get the view size (in cylinder diameters) used for field computation.
#[wasm_bindgen]
pub fn get_field_view_size() -> f64 {
    VIEW_SIZE
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
    };
    let scattering = calculate_scattering(&scat_params);

    // Decompose Complex64 → real/imag for FieldParams (stays on WASM heap)
    let n_coeffs = scattering.orders.len();
    let mut inc_re = Vec::with_capacity(n_coeffs);
    let mut inc_im = Vec::with_capacity(n_coeffs);
    let mut scat_re = Vec::with_capacity(n_coeffs);
    let mut scat_im = Vec::with_capacity(n_coeffs);
    let mut int_re = Vec::with_capacity(n_coeffs);
    let mut int_im = Vec::with_capacity(n_coeffs);

    for i in 0..n_coeffs {
        inc_re.push(scattering.incident_coefficients[i].re);
        inc_im.push(scattering.incident_coefficients[i].im);
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
        incident_coeffs_real: inc_re,
        incident_coeffs_imag: inc_im,
        scattering_coeffs_real: scat_re,
        scattering_coeffs_imag: scat_im,
        internal_coeffs_real: int_re,
        internal_coeffs_imag: int_im,
        orders: scattering.orders,
    };

    let field = compute_field(&field_params);

    // Copy into pre-allocated JS arrays — only JS-side cost is 2 Float64Array.set() calls
    out_field_real.copy_from(&field.field_real);
    out_field_imag.copy_from(&field.field_imag);
}

/// Get information about the scattering computation.
#[wasm_bindgen]
pub fn get_info() -> String {
    String::from(
        "2D Electromagnetic Scattering Calculator\n\
         Computes scattering from an infinite cylinder.\n\
         Cylinder diameter = 1 (unitless)\n\
         Field grid: 128x128 over 5D x 5D area\n\
         Supports TM and TE polarizations.",
    )
}

// Re-export types for testing
pub use field::FieldResult as FieldResultType;
pub use scattering::ScatteringResult as ScatteringResultType;
