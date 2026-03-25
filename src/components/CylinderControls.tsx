import { useState, useRef, useEffect, memo } from "react";
import type { JSX } from "react";
import type {
  ScatteringParams,
  ComplexObj,
  SourceType,
  DipoleParams,
} from "../types/cylinder";
import { createDefaultParams, isDipoleSource } from "../types/cylinder";
import type { ParameterBounds } from "../hooks/useScattering";
import { ComplexPlane } from "./ComplexPlane";
import "./CylinderControls.css";

const SOURCE_OPTIONS: { value: SourceType; label: JSX.Element }[] = [
  {
    value: "planewave_tm",
    label: (
      <span>
        Plane Wave TM (E<sub>z</sub>)
      </span>
    ),
  },
  {
    value: "planewave_te",
    label: (
      <span>
        Plane Wave TE (H<sub>z</sub>)
      </span>
    ),
  },
  {
    value: "dipole_ez",
    label: (
      <span>
        Dipole E<sub>z</sub>
      </span>
    ),
  },
  {
    value: "dipole_exy",
    label: (
      <span>
        Dipole E<sub>xy</sub>
      </span>
    ),
  },
];

export type DipoleUpdateFn = (dipole: DipoleParams) => void;

interface CylinderControlsProps {
  onChange: (params: ScatteringParams) => void;
  onShowOverlayChange: (show: boolean) => void;
  bounds: ParameterBounds | null;
  /** Ref that CylinderControls populates with its updateDipole function,
   *  allowing FieldVisualization to route drag updates through CylinderControls' state. */
  dipoleUpdateRef: React.RefObject<DipoleUpdateFn | null>;
}

function SourceDropdown({
  value,
  onChange,
}: {
  value: SourceType;
  onChange: (v: SourceType) => void;
}) {
  const [open, setOpen] = useState(false);
  const ref = useRef<HTMLDivElement>(null);

  // Close on outside click
  useEffect(() => {
    if (!open) return;
    const handler = (e: MouseEvent) => {
      if (ref.current && !ref.current.contains(e.target as Node)) {
        setOpen(false);
      }
    };
    document.addEventListener("mousedown", handler);
    return () => document.removeEventListener("mousedown", handler);
  }, [open]);

  const selected = SOURCE_OPTIONS.find((o) => o.value === value)!;

  return (
    <div className="control-section">
      <label className="section-label">Source</label>
      <div className="source-dropdown" ref={ref}>
        <button
          className="source-dropdown-trigger"
          onClick={() => setOpen(!open)}
        >
          {selected.label}
          <span className="source-dropdown-arrow">{open ? "▴" : "▾"}</span>
        </button>
        {open && (
          <div className="source-dropdown-menu">
            {SOURCE_OPTIONS.map((opt) => (
              <button
                key={opt.value}
                className={
                  "source-dropdown-item" +
                  (opt.value === value ? " active" : "")
                }
                onClick={() => {
                  onChange(opt.value);
                  setOpen(false);
                }}
              >
                {opt.label}
              </button>
            ))}
          </div>
        )}
      </div>
    </div>
  );
}

export const CylinderControls = memo(function CylinderControls({
  onChange,
  onShowOverlayChange,
  bounds,
  dipoleUpdateRef,
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

  const updateSourceType = (sourceType: SourceType) => {
    update({ ...local, sourceType });
  };

  const updateDipole = (dipole: DipoleParams) => {
    update({ ...local, dipole });
  };

  // Expose updateDipole so FieldVisualization drag routes through our state
  useEffect(() => {
    dipoleUpdateRef.current = updateDipole;
  });

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

  const showDipole = isDipoleSource(local.sourceType);

  return (
    <div className="cylinder-controls">
      <div className="controls-header">
        <h2>Simulation Parameters</h2>
        <button
          className="reset-button"
          title="Reset simulation to default settings"
          onClick={() =>
            update({ ...createDefaultParams(), dipole: local.dipole })
          }
        >
          ↻
        </button>
      </div>

      <SourceDropdown value={local.sourceType} onChange={updateSourceType} />

      {showDipole && (
        <div className="control-section dipole-section">
          <label className="section-label">Dipole Position</label>
          <div className="dipole-coords">
            <label>
              x₀
              <input
                type="number"
                step={0.1}
                value={local.dipole.xs.toFixed(2)}
                onChange={(e) => {
                  const val = parseFloat(e.target.value);
                  if (!isNaN(val)) updateDipole({ ...local.dipole, xs: val });
                }}
              />
            </label>
            <label>
              y₀
              <input
                type="number"
                step={0.1}
                value={local.dipole.ys.toFixed(2)}
                onChange={(e) => {
                  const val = parseFloat(e.target.value);
                  if (!isNaN(val)) updateDipole({ ...local.dipole, ys: val });
                }}
              />
            </label>
          </div>
          {local.sourceType === "dipole_exy" && (
            <div className="dipole-orientation">
              <label>
                α
                <input
                  type="range"
                  min={-Math.PI}
                  max={Math.PI}
                  step={0.01}
                  value={local.dipole.alpha}
                  onChange={(e) =>
                    updateDipole({
                      ...local.dipole,
                      alpha: parseFloat(e.target.value),
                    })
                  }
                />
                <span className="alpha-display">
                  {((local.dipole.alpha * 180) / Math.PI).toFixed(0)}°
                </span>
              </label>
            </div>
          )}
        </div>
      )}

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
          <div className="stepper-wrap">
            <span className="stepper-display">
              {Math.abs(local.material.permittivity.re).toFixed(1)}
            </span>
            <input
              type="number"
              min={bounds.permittivity_re_min}
              max={bounds.permittivity_re_max}
              step={0.1}
              value={local.material.permittivity.re.toFixed(1)}
              onChange={(e) => {
                const val = parseFloat(e.target.value);
                if (!isNaN(val))
                  updatePermittivity({
                    ...local.material.permittivity,
                    re: val,
                  });
              }}
            />
          </div>
          <span className="complex-inline-op">
            {local.material.permittivity.im < 0 ? "−" : "+"}
          </span>
          <span className="complex-inline-i">i</span>
          <div className="stepper-wrap">
            <span className="stepper-display">
              {Math.abs(local.material.permittivity.im).toFixed(1)}
            </span>
            <input
              type="number"
              min={bounds.permittivity_im_min}
              max={bounds.permittivity_im_max}
              step={0.1}
              value={local.material.permittivity.im.toFixed(1)}
              onChange={(e) => {
                const val = parseFloat(e.target.value);
                if (!isNaN(val))
                  updatePermittivity({
                    ...local.material.permittivity,
                    im: val,
                  });
              }}
            />
          </div>
        </div>

        <div className="complex-inline-row mu-row">
          <span className="complex-inline-label">μᵣ</span>
          <span className="complex-inline-eq">=</span>
          <span className="complex-inline-sign">
            {local.material.permeability.re < 0 ? "−" : "\u2007"}
          </span>
          <div className="stepper-wrap">
            <span className="stepper-display">
              {Math.abs(local.material.permeability.re).toFixed(1)}
            </span>
            <input
              type="number"
              min={bounds.permeability_re_min}
              max={bounds.permeability_re_max}
              step={0.1}
              value={local.material.permeability.re.toFixed(1)}
              onChange={(e) => {
                const val = parseFloat(e.target.value);
                if (!isNaN(val))
                  updatePermeability({
                    ...local.material.permeability,
                    re: val,
                  });
              }}
            />
          </div>
          <span className="complex-inline-op">
            {local.material.permeability.im < 0 ? "−" : "+"}
          </span>
          <span className="complex-inline-i">i</span>
          <div className="stepper-wrap">
            <span className="stepper-display">
              {Math.abs(local.material.permeability.im).toFixed(1)}
            </span>
            <input
              type="number"
              min={bounds.permeability_im_min}
              max={bounds.permeability_im_max}
              step={0.1}
              value={local.material.permeability.im.toFixed(1)}
              onChange={(e) => {
                const val = parseFloat(e.target.value);
                if (!isNaN(val))
                  updatePermeability({
                    ...local.material.permeability,
                    im: val,
                  });
              }}
            />
          </div>
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
