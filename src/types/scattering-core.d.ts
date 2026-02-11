declare module "scattering-core" {
  export function compute_scattering(
    wavelength: number,
    permittivity_re: number,
    permittivity_im: number,
    permeability_re: number,
    permeability_im: number,
    polarization: string,
    max_order: number,
  ): {
    incident_coefficients: Array<[number, number]>;
    scattering_coefficients: Array<[number, number]>;
    internal_coefficients: Array<[number, number]>;
    orders: number[];
  };

  export function compute_electric_field(
    wavelength: number,
    permittivity_re: number,
    permittivity_im: number,
    permeability_re: number,
    permeability_im: number,
    incident_coeffs_re: Float64Array,
    incident_coeffs_im: Float64Array,
    scattering_coeffs_re: Float64Array,
    scattering_coeffs_im: Float64Array,
    internal_coeffs_re: Float64Array,
    internal_coeffs_im: Float64Array,
    orders: Int32Array,
    view_size: number,
  ): {
    field_real: number[];
    field_imag: number[];
    grid_size: number;
    view_size: number;
    x_min: number;
    x_max: number;
    y_min: number;
    y_max: number;
  };

  export function get_field_grid_size(): number;
  export function get_field_view_size(): number;
  export function get_info(): string;
  export function init(): void;

  export default function initWasm(): Promise<void>;
}
