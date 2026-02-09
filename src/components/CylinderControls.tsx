import { useState, memo } from "react";
import type {
  ScatteringParams,
  ComplexObj,
  Polarization,
} from "../types/cylinder";
import { createDefaultParams, calculateSizeParameter } from "../types/cylinder";
import "./CylinderControls.css";

interface CylinderControlsProps {
  onChange: (params: ScatteringParams) => void;
}

interface ComplexInputProps {
  label: string;
  value: ComplexObj;
  onChange: (value: ComplexObj) => void;
  realMin?: number;
  realMax?: number;
  imagMin?: number;
  imagMax?: number;
  realStep?: number;
  imagStep?: number;
}

function ComplexInput({
  label,
  value,
  onChange,
  realMin = -10,
  realMax = 10,
  imagMin = -10,
  imagMax = 10,
  realStep = 0.1,
  imagStep = 0.1,
}: ComplexInputProps) {
  return (
    <div className="complex-input">
      <label className="complex-label">{label}</label>
      <div className="complex-fields">
        <div className="field-group">
          <label className="field-label">Real</label>
          <input
            type="range"
            min={realMin}
            max={realMax}
            step={realStep}
            value={value.re}
            onChange={(e) =>
              onChange({ ...value, re: parseFloat(e.target.value) })
            }
          />
          <input
            type="number"
            min={realMin}
            max={realMax}
            step={realStep}
            value={value.re}
            onChange={(e) =>
              onChange({ ...value, re: parseFloat(e.target.value) || 0 })
            }
          />
        </div>
        <div className="field-group">
          <label className="field-label">Imag</label>
          <input
            type="range"
            min={imagMin}
            max={imagMax}
            step={imagStep}
            value={value.im}
            onChange={(e) =>
              onChange({ ...value, im: parseFloat(e.target.value) })
            }
          />
          <input
            type="number"
            min={imagMin}
            max={imagMax}
            step={imagStep}
            value={value.im}
            onChange={(e) =>
              onChange({ ...value, im: parseFloat(e.target.value) || 0 })
            }
          />
        </div>
      </div>
    </div>
  );
}

export const CylinderControls = memo(function CylinderControls({
  onChange,
}: CylinderControlsProps) {
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

  const sizeParameter = calculateSizeParameter(local.wavelength);

  return (
    <div className="cylinder-controls">
      <h2>Simulation Parameters</h2>

      <div className="control-section">
        <label className="section-label">Wavelength (λ / d)</label>
        <p className="section-hint">Relative to cylinder diameter = 1</p>
        <div className="wavelength-control">
          <input
            type="range"
            min={0.1}
            max={5}
            step={0.01}
            value={local.wavelength}
            onChange={(e) => updateWavelength(parseFloat(e.target.value))}
          />
          <input
            type="number"
            min={0.01}
            max={10}
            step={0.01}
            value={local.wavelength}
            onChange={(e) =>
              updateWavelength(parseFloat(e.target.value) || 0.1)
            }
          />
        </div>
        <div className="derived-value">
          Size parameter: x = πd/λ = {sizeParameter.toFixed(3)}
        </div>
      </div>

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
        <label className="section-label">Max Bessel Order (N)</label>
        <p className="section-hint">Computes orders -N to +N</p>
        <div className="order-control">
          <input
            type="range"
            min={1}
            max={30}
            step={1}
            value={local.maxOrder}
            onChange={(e) => updateMaxOrder(parseInt(e.target.value))}
          />
          <input
            type="number"
            min={1}
            max={50}
            step={1}
            value={local.maxOrder}
            onChange={(e) => updateMaxOrder(parseInt(e.target.value) || 1)}
          />
        </div>
        <div className="derived-value">Terms: {2 * local.maxOrder + 1}</div>
      </div>

      <div className="control-section">
        <ComplexInput
          label="Relative Permittivity (εᵣ)"
          value={local.material.permittivity}
          onChange={updatePermittivity}
          realMin={-10}
          realMax={20}
          imagMin={-5}
          imagMax={5}
        />
      </div>

      <div className="control-section">
        <ComplexInput
          label="Relative Permeability (μᵣ)"
          value={local.material.permeability}
          onChange={updatePermeability}
          realMin={-10}
          realMax={10}
          imagMin={-5}
          imagMax={5}
        />
      </div>
    </div>
  );
});
