declare module "scattering-core" {
  export function compute_all(
    wavelength: number,
    permittivity_real: number,
    permittivity_imag: number,
    permeability_real: number,
    permeability_imag: number,
    polarization: number,
    max_order: number,
    out_field_real: Float64Array,
    out_field_imag: Float64Array,
  ): void;

  export function compute_scattering(
    wavelength: number,
    permittivity_real: number,
    permittivity_imag: number,
    permeability_real: number,
    permeability_imag: number,
    polarization: string,
    max_order: number,
  ): unknown;

  export function get_field_grid_size(): number;
  export function get_field_view_size(): number;
  export function get_info(): string;
  export function init(): void;

  export default function __wbg_init(): Promise<unknown>;
}
