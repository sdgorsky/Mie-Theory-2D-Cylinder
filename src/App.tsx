import { useState, useEffect, useCallback } from "react";
import { CylinderControls } from "./components/CylinderControls";
import { FieldVisualization } from "./components/FieldVisualization";
import { useScattering } from "./hooks/useScattering";
import {
  createDefaultParams,
  calculateRefractiveIndex,
  formatComplex,
} from "./types/cylinder";
import type { ScatteringParams } from "./types/cylinder";
import "./App.css";

function App() {
  const {
    isLoading,
    isReady,
    error,
    parameterBounds,
    scatteringResult,
    computedParams,
    imageStats,
    hasImage,
    paintRef,
    computeAll,
    recolor,
  } = useScattering();

  const [showOverlay, setShowOverlay] = useState(true);

  // Called on every slider event — non-blocking, posts to worker.
  // No setState here: CylinderControls owns its own display state,
  // and App's info section updates when worker results arrive.
  const handleParamsChange = useCallback(
    (p: ScatteringParams) => {
      computeAll(p);
    },
    [computeAll],
  );

  // Initial computation when WASM worker is ready
  useEffect(() => {
    if (isReady) {
      computeAll(createDefaultParams());
    }
  }, [isReady, computeAll]);

  // Use computedParams (set when worker returns) for the info section.
  // Falls back to defaults before first result arrives.
  const displayParams = computedParams ?? createDefaultParams();

  return (
    <div className="app">
      <header className="app-header">
        <h1>2D Electromagnetic Scattering</h1>
        <p>
          Interactive visualization of EM wave scattering from an infinite
          cylinder
        </p>
      </header>

      <main className="app-main">
        <div className="visualization-panel">
          <FieldVisualization
            paintRef={paintRef}
            imageStats={imageStats}
            hasImage={hasImage}
            polarization={displayParams.polarization}
            showOverlay={showOverlay}
            onModeChange={recolor}
            width={512}
            height={512}
          />

          <div className="status-section">
            {isLoading && (
              <span className="status">Loading WASM module...</span>
            )}
            {error && <span className="status error">{error}</span>}
          </div>

          <div className="scattering-info">
            <h3>Scattering Info</h3>
            <p>
              Permittivity: εᵣ ={" "}
              {formatComplex(displayParams.material.permittivity, 1)}
            </p>
            <p>
              Permeability: μᵣ ={" "}
              {formatComplex(displayParams.material.permeability, 1)}
            </p>
            <p>
              Refractive index: n ={" "}
              {formatComplex(
                calculateRefractiveIndex(displayParams.material),
                2,
              )}
            </p>
            {scatteringResult && (
              <p>
                Orders: {scatteringResult.minOrder} to{" "}
                {scatteringResult.maxOrder} ({scatteringResult.numOrders} terms)
              </p>
            )}
          </div>
        </div>

        <div className="controls-panel">
          <CylinderControls
            onChange={handleParamsChange}
            onShowOverlayChange={setShowOverlay}
            bounds={parameterBounds}
          />
        </div>
      </main>
    </div>
  );
}

export default App;
