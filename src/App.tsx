import { useState, useEffect, useCallback, useRef } from "react";
import { CylinderControls } from "./components/CylinderControls";
import { FieldVisualization } from "./components/FieldVisualization";
import { useScattering } from "./hooks/useScattering";
import {
  createDefaultParams,
  calculateRefractiveIndex,
  formatComplex,
  getPolarization,
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
    setViewSize,
  } = useScattering();

  const [showOverlay, setShowOverlay] = useState(true);
  const dipoleUpdateRef = useRef<
    ((d: import("./types/cylinder").DipoleParams) => void) | null
  >(null);

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
        <h1>
          2D Electromagnetic Scattering
          <a
            href="https://github.com/sdgorsky/webapp-single-cylinder-scattering-2D"
            target="_blank"
            rel="noopener noreferrer"
            className="github-link"
            title="View on GitHub"
          >
            <svg
              viewBox="0 0 16 16"
              width="24"
              height="24"
              fill="currentColor"
              aria-hidden="true"
            >
              <path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z" />
            </svg>
          </a>
          <span className="version-tag">v{__APP_VERSION__}</span>
        </h1>
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
            polarization={getPolarization(displayParams.sourceType)}
            sourceType={displayParams.sourceType}
            dipole={displayParams.dipole}
            showOverlay={showOverlay}
            onModeChange={recolor}
            onZoomChange={setViewSize}
            dipoleUpdateRef={dipoleUpdateRef}
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
            dipoleUpdateRef={dipoleUpdateRef}
          />
        </div>
      </main>
    </div>
  );
}

export default App;
