import type { JSX, MutableRefObject } from "react";
import { useRef, useEffect, useState, memo } from "react";
import type { Polarization } from "../types/cylinder";
import type {
  VisualizationMode,
  ImageStats,
  PaintCallback,
} from "../hooks/useScattering";
import "./FieldVisualization.css";

export type { VisualizationMode };

interface FieldVisualizationProps {
  paintRef: MutableRefObject<PaintCallback | null>;
  imageStats: ImageStats | null;
  hasImage: boolean;
  polarization: Polarization;
  showOverlay: boolean;
  onModeChange: (mode: VisualizationMode) => void;
  width?: number;
  height?: number;
}

export const FieldVisualization = memo(function FieldVisualization({
  paintRef,
  imageStats,
  hasImage,
  polarization,
  showOverlay,
  onModeChange,
  width = 512,
  height = 512,
}: FieldVisualizationProps) {
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const tempCanvasRef = useRef<HTMLCanvasElement>(
    document.createElement("canvas"),
  );
  const showOverlayRef = useRef(showOverlay);
  const lastViewSizeRef = useRef(0);
  const [mode, setMode] = useState<VisualizationMode>("magnitude");

  const handleModeClick = (m: VisualizationMode) => {
    setMode(m);
    onModeChange(m);
  };

  // Register imperative paint callback — called from rAF, bypasses React render
  useEffect(() => {
    paintRef.current = (
      imageData: Uint8ClampedArray,
      gridSize: number,
      viewSize: number,
    ) => {
      const canvas = canvasRef.current;
      if (!canvas || gridSize === 0) return;

      const ctx = canvas.getContext("2d");
      if (!ctx) return;

      // Blit RGBA to temp canvas at full resolution
      const tempCanvas = tempCanvasRef.current;
      if (tempCanvas.width !== gridSize) tempCanvas.width = gridSize;
      if (tempCanvas.height !== gridSize) tempCanvas.height = gridSize;
      const tempCtx = tempCanvas.getContext("2d")!;
      tempCtx.putImageData(
        new ImageData(
          new Uint8ClampedArray(imageData.buffer as ArrayBuffer),
          gridSize,
          gridSize,
        ),
        0,
        0,
      );

      // Scale to display canvas
      ctx.imageSmoothingEnabled = false;
      ctx.drawImage(tempCanvas, 0, 0, width, height);
      lastViewSizeRef.current = viewSize;

      if (showOverlayRef.current) {
        drawOverlay(ctx, width, height, viewSize);
      }
    };

    return () => {
      paintRef.current = null;
    };
  }, [paintRef, width, height]);

  // Repaint from temp canvas when overlay toggle changes
  useEffect(() => {
    showOverlayRef.current = showOverlay;
    const canvas = canvasRef.current;
    const tempCanvas = tempCanvasRef.current;
    if (!canvas || tempCanvas.width === 0) return;
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(tempCanvas, 0, 0, width, height);
    if (showOverlay) {
      drawOverlay(ctx, width, height, lastViewSizeRef.current);
    }
  }, [showOverlay, width, height]);

  const fieldSymbol = polarization === "TM" ? "E" : "H";
  const modeLabels: Record<VisualizationMode, JSX.Element> = {
    magnitude: (
      <span>
        |{fieldSymbol}
        <sub>z</sub>|
      </span>
    ),
    real: (
      <span>
        Re({fieldSymbol}
        <sub>z</sub>)
      </span>
    ),
    imag: (
      <span>
        Im({fieldSymbol}
        <sub>z</sub>)
      </span>
    ),
    phase: (
      <span>
        Phase({fieldSymbol}
        <sub>z</sub>)
      </span>
    ),
  };

  return (
    <div className="field-visualization">
      <div className="visualization-controls">
        <label>Display:</label>
        <div className="mode-buttons">
          {(Object.keys(modeLabels) as VisualizationMode[]).map((m) => (
            <button
              key={m}
              className={mode === m ? "active" : ""}
              onClick={() => handleModeClick(m)}
            >
              {modeLabels[m]}
            </button>
          ))}
        </div>
      </div>

      <div className="canvas-container">
        <canvas ref={canvasRef} width={width} height={height} />
      </div>

      {imageStats && (
        <div className="field-stats">
          <span>
            Range: [{imageStats.min.toExponential(2)},{" "}
            {imageStats.max.toExponential(2)}]
          </span>
        </div>
      )}

      {!hasImage && (
        <div className="no-data">
          <p>Loading...</p>
        </div>
      )}
    </div>
  );
});

function drawOverlay(
  ctx: CanvasRenderingContext2D,
  width: number,
  height: number,
  viewSize: number,
) {
  const centerX = width / 2;
  const centerY = height / 2;
  const cylinderRadius = (0.5 / viewSize) * width;

  // Cylinder outline
  ctx.strokeStyle = "rgba(255, 255, 255, 0.8)";
  ctx.lineWidth = 2;
  ctx.beginPath();
  ctx.arc(centerX, centerY, cylinderRadius, 0, 2 * Math.PI);
  ctx.stroke();

  // Axis markers
  ctx.strokeStyle = "rgba(255, 255, 255, 0.3)";
  ctx.lineWidth = 1;
  ctx.setLineDash([5, 5]);

  ctx.beginPath();
  ctx.moveTo(0, centerY);
  ctx.lineTo(width, centerY);
  ctx.stroke();

  ctx.beginPath();
  ctx.moveTo(centerX, 0);
  ctx.lineTo(centerX, height);
  ctx.stroke();

  ctx.setLineDash([]);
}
