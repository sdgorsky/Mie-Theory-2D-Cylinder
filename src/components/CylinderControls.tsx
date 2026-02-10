import { useState, memo } from "react";
import type {
  ScatteringParams,
  ComplexObj,
  Polarization,
} from "../types/cylinder";
import { createDefaultParams } from "../types/cylinder";
import type { ParameterBounds } from "../hooks/useScattering";
import { ComplexPlane } from "./ComplexPlane";
import "./CylinderControls.css";

interface CylinderControlsProps {
  onChange: (params: ScatteringParams) => void;
  onShowOverlayChange: (show: boolean) => void;
  bounds: ParameterBounds | null;
}

export const CylinderControls = memo(function CylinderControls({
  onChange,
  onShowOverlayChange,
  bounds,
}: CylinderControlsProps) {
  const [showOverlay, setShowOverlay] = useState(true);
  // CylinderControls owns its state — no params prop means React.memo blocks
  // all parent re-renders (onChange is a stable ref-writing callback).
  const [local, setLocal] = useState(createDefaultParams);

  const update = (next: ScatteringParams) => {
    setLocal(next);
    onChange(next); // just a ref write in App — instant, no React work
  };

  const updateWavelength = (wavelength: number) => {
    update({ ...local, wavelength });
  };

  const updatePolarization = (polarization: Polarization) => {
    update({ ...local, polarization });
  };

  const updateMaxOrder = (maxOrder: number) => {
    update({ ...local, maxOrder });
  };

  const updatePermittivity = (permittivity: ComplexObj) => {
    update({
      ...local,
      material: { ...local.material, permittivity },
    });
  };

  const updatePermeability = (permeability: ComplexObj) => {
    update({
      ...local,
      material: { ...local.material, permeability },
    });
  };

  // Wait for WASM bounds before rendering controls
  if (!bounds) {
    return (
      <div className="cylinder-controls">
        <h2>Simulation Parameters</h2>
        <p>Loading...</p>
      </div>
    );
  }

  return (
    <div className="cylinder-controls">
      <h2>Simulation Parameters</h2>

      <div className="control-section">
        <label className="section-label">Polarization</label>
        <div className="polarization-buttons">
          <button
            className={local.polarization === "TM" ? "active" : ""}
            onClick={() => updatePolarization("TM")}
          >
            TM (E∥z)
          </button>
          <button
            className={local.polarization === "TE" ? "active" : ""}
            onClick={() => updatePolarization("TE")}
          >
            TE (H∥z)
          </button>
        </div>
      </div>

      <div className="control-section">
        <label className="section-label">Wavelength (λ / d)</label>
        <div className="wavelength-control">
          <input
            type="range"
            min={bounds.wavelength_min}
            max={bounds.wavelength_max}
            step={0.01}
            value={local.wavelength}
            onChange={(e) => updateWavelength(parseFloat(e.target.value))}
          />
          <input
            type="number"
            min={bounds.wavelength_min}
            max={bounds.wavelength_max}
            step={0.01}
            value={local.wavelength}
            onChange={(e) =>
              updateWavelength(
                parseFloat(e.target.value) || bounds.wavelength_min,
              )
            }
          />
        </div>
      </div>

      <div className="control-section">
        <label className="section-label">Material Properties</label>
        <ComplexPlane
          permittivity={local.material.permittivity}
          permeability={local.material.permeability}
          onPermittivityChange={updatePermittivity}
          onPermeabilityChange={updatePermeability}
          bounds={bounds}
        />

        <div className="complex-inline-row eps-row">
          <span className="complex-inline-label">εᵣ</span>
          <span className="complex-inline-eq">=</span>
          <span className="complex-inline-sign">
            {local.material.permittivity.re < 0 ? "−" : "\u2007"}
          </span>
          <input
            type="number"
            min={bounds.permittivity_re_min}
            max={bounds.permittivity_re_max}
            step={0.1}
            value={Math.abs(local.material.permittivity.re).toFixed(1)}
            onChange={(e) => {
              const abs = parseFloat(e.target.value) || 0;
              const sign = local.material.permittivity.re < 0 ? -1 : 1;
              updatePermittivity({
                ...local.material.permittivity,
                re: sign * abs,
              });
            }}
          />
          <span className="complex-inline-op">
            {local.material.permittivity.im < 0 ? "−" : "+"}
          </span>
          <span className="complex-inline-i">i</span>
          <input
            type="number"
            min={bounds.permittivity_im_min}
            max={bounds.permittivity_im_max}
            step={0.1}
            value={Math.abs(local.material.permittivity.im).toFixed(1)}
            onChange={(e) => {
              const abs = parseFloat(e.target.value) || 0;
              const sign = local.material.permittivity.im < 0 ? -1 : 1;
              updatePermittivity({
                ...local.material.permittivity,
                im: sign * abs,
              });
            }}
          />
        </div>

        <div className="complex-inline-row mu-row">
          <span className="complex-inline-label">μᵣ</span>
          <span className="complex-inline-eq">=</span>
          <span className="complex-inline-sign">
            {local.material.permeability.re < 0 ? "−" : "\u2007"}
          </span>
          <input
            type="number"
            min={bounds.permeability_re_min}
            max={bounds.permeability_re_max}
            step={0.1}
            value={Math.abs(local.material.permeability.re).toFixed(1)}
            onChange={(e) => {
              const abs = parseFloat(e.target.value) || 0;
              const sign = local.material.permeability.re < 0 ? -1 : 1;
              updatePermeability({
                ...local.material.permeability,
                re: sign * abs,
              });
            }}
          />
          <span className="complex-inline-op">
            {local.material.permeability.im < 0 ? "−" : "+"}
          </span>
          <span className="complex-inline-i">i</span>
          <input
            type="number"
            min={bounds.permeability_im_min}
            max={bounds.permeability_im_max}
            step={0.1}
            value={Math.abs(local.material.permeability.im).toFixed(1)}
            onChange={(e) => {
              const abs = parseFloat(e.target.value) || 0;
              const sign = local.material.permeability.im < 0 ? -1 : 1;
              updatePermeability({
                ...local.material.permeability,
                im: sign * abs,
              });
            }}
          />
        </div>
      </div>

      <details className="advanced-section">
        <summary className="section-label">Advanced</summary>
        <div className="control-section">
          <label className="section-label">Max Bessel Order (N)</label>
          <div className="order-control">
            <input
              type="range"
              min={bounds.max_order_min}
              max={bounds.max_order_max}
              step={1}
              value={local.maxOrder}
              onChange={(e) => updateMaxOrder(parseInt(e.target.value))}
            />
            <input
              type="number"
              min={bounds.max_order_min}
              max={bounds.max_order_max}
              step={1}
              value={local.maxOrder}
              onChange={(e) =>
                updateMaxOrder(parseInt(e.target.value) || bounds.max_order_min)
              }
            />
          </div>
        </div>
        <label className="checkbox-row">
          <input
            type="checkbox"
            checked={!showOverlay}
            onChange={(e) => {
              const hide = e.target.checked;
              setShowOverlay(!hide);
              onShowOverlayChange(!hide);
            }}
          />
          Hide cylinder overlay
        </label>
      </details>
    </div>
  );
});
