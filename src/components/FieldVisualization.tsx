import type { JSX, MutableRefObject } from "react";
import { useRef, useEffect, useState, memo, useCallback } from "react";
import type { Polarization, SourceType, DipoleParams } from "../types/cylinder";
import type { DipoleUpdateFn } from "./CylinderControls";
import { isDipoleSource } from "../types/cylinder";
import type {
  VisualizationMode,
  ImageStats,
  PaintCallback,
} from "../hooks/useScattering";
import "./FieldVisualization.css";

export type { VisualizationMode };

const ZOOM_OPTIONS = [
  { label: "0.5x", viewSize: 10 },
  { label: "1.0x", viewSize: 5 },
  { label: "2.0x", viewSize: 2.5 },
];

/** Length of the orientation arm in pixels. */
const ARM_LENGTH_PX = 24;

interface FieldVisualizationProps {
  paintRef: MutableRefObject<PaintCallback | null>;
  imageStats: ImageStats | null;
  hasImage: boolean;
  polarization: Polarization;
  sourceType: SourceType;
  dipole: DipoleParams;
  showOverlay: boolean;
  onModeChange: (mode: VisualizationMode) => void;
  onZoomChange: (viewSize: number) => void;
  dipoleUpdateRef: React.RefObject<DipoleUpdateFn | null>;
  width?: number;
  height?: number;
}

type DragTarget = "position" | "orientation" | null;

export const FieldVisualization = memo(function FieldVisualization({
  paintRef,
  imageStats,
  hasImage,
  polarization,
  sourceType,
  dipole,
  showOverlay,
  onModeChange,
  onZoomChange,
  dipoleUpdateRef,
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
  const [zoomIndex, setZoomIndex] = useState(1);

  // Refs for dipole interaction (avoid stale closures)
  const dipoleRef = useRef(dipole);
  const sourceTypeRef = useRef(sourceType);
  const dragTargetRef = useRef<DragTarget>(null);
  const hoveringCanvasRef = useRef(false);

  // Sync props into refs for use in imperative callbacks (rAF, mouse handlers)
  useEffect(() => {
    dipoleRef.current = dipole;
    sourceTypeRef.current = sourceType;
  });

  const handleModeClick = (m: VisualizationMode) => {
    setMode(m);
    onModeChange(m);
  };

  // Convert physical coordinates to canvas pixels
  const physToCanvas = useCallback(
    (x: number, y: number, viewSize: number) => {
      const scale = width / viewSize;
      return {
        cx: width / 2 + x * scale,
        cy: height / 2 - y * scale, // y-axis inverted
      };
    },
    [width, height],
  );

  // Convert canvas pixels to physical coordinates
  const canvasToPhys = useCallback(
    (cx: number, cy: number, viewSize: number) => {
      const scale = width / viewSize;
      return {
        x: (cx - width / 2) / scale,
        y: -(cy - height / 2) / scale, // y-axis inverted
      };
    },
    [width, height],
  );

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
        if (
          isDipoleSource(sourceTypeRef.current) &&
          hoveringCanvasRef.current
        ) {
          drawDipoleMarker(
            ctx,
            dipoleRef.current,
            sourceTypeRef.current,
            viewSize,
            width,
            height,
          );
        }
      }
    };

    return () => {
      paintRef.current = null;
    };
  }, [paintRef, width, height]);

  // Repaint from temp canvas when overlay toggle or dipole params change
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
      if (isDipoleSource(sourceType) && hoveringCanvasRef.current) {
        drawDipoleMarker(
          ctx,
          dipole,
          sourceType,
          lastViewSizeRef.current,
          width,
          height,
        );
      }
    }
    // dipole intentionally excluded: marker position is handled by
    // dipoleRef + repaintOverlay (immediate) and paintRef (on worker result).
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [showOverlay, width, height, sourceType]);

  // Repaint overlay from temp canvas (used for dipole marker updates and hover)
  const repaintOverlay = useCallback(() => {
    const canvas = canvasRef.current;
    const tempCanvas = tempCanvasRef.current;
    if (!canvas || tempCanvas.width === 0) return;
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    ctx.imageSmoothingEnabled = false;
    ctx.drawImage(tempCanvas, 0, 0, width, height);
    if (showOverlayRef.current) {
      drawOverlay(ctx, width, height, lastViewSizeRef.current);
      if (isDipoleSource(sourceTypeRef.current) && hoveringCanvasRef.current) {
        drawDipoleMarker(
          ctx,
          dipoleRef.current,
          sourceTypeRef.current,
          lastViewSizeRef.current,
          width,
          height,
        );
      }
    }
  }, [width, height]);

  // Mouse interaction for dipole dragging
  const handleMouseDown = useCallback(
    (e: React.MouseEvent<HTMLCanvasElement>) => {
      if (!isDipoleSource(sourceTypeRef.current)) return;
      const viewSize = lastViewSizeRef.current;
      if (!viewSize) return;

      const rect = canvasRef.current!.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const my = e.clientY - rect.top;

      const d = dipoleRef.current;
      const { cx, cy } = physToCanvas(d.xs, d.ys, viewSize);

      // Check orientation arm endpoint first (smaller target)
      if (sourceTypeRef.current === "dipole_exy") {
        const armX = cx + ARM_LENGTH_PX * Math.cos(d.alpha);
        const armY = cy - ARM_LENGTH_PX * Math.sin(d.alpha); // y inverted
        const armDist = Math.hypot(mx - armX, my - armY);
        if (armDist < 12) {
          dragTargetRef.current = "orientation";
          return;
        }
      }

      // Check dipole position
      const posDist = Math.hypot(mx - cx, my - cy);
      if (posDist < 14) {
        dragTargetRef.current = "position";
      }
    },
    [physToCanvas],
  );

  const handleMouseMove = useCallback(
    (e: React.MouseEvent<HTMLCanvasElement>) => {
      const target = dragTargetRef.current;
      if (!target) return;

      const viewSize = lastViewSizeRef.current;
      if (!viewSize) return;

      const rect = canvasRef.current!.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const my = e.clientY - rect.top;

      const d = dipoleRef.current;

      let newDipole: typeof d;
      if (target === "position") {
        const { x, y } = canvasToPhys(mx, my, viewSize);
        newDipole = {
          ...d,
          xs: Math.round(x * 20) / 20,
          ys: Math.round(y * 20) / 20,
        };
      } else {
        const { cx, cy } = physToCanvas(d.xs, d.ys, viewSize);
        const alpha = Math.atan2(-(my - cy), mx - cx); // y inverted
        newDipole = { ...d, alpha: Math.round(alpha * 100) / 100 };
      }

      // Update ref immediately so the marker tracks the mouse without
      // waiting for the worker round-trip, then repaint overlay.
      dipoleRef.current = newDipole;
      repaintOverlay();

      // Dispatch to controls → worker (async)
      dipoleUpdateRef.current?.(newDipole);
    },
    [canvasToPhys, physToCanvas, dipoleUpdateRef, repaintOverlay],
  );

  const handleMouseEnter = useCallback(() => {
    hoveringCanvasRef.current = true;
    repaintOverlay();
  }, [repaintOverlay]);

  const handleMouseLeave = useCallback(() => {
    hoveringCanvasRef.current = false;
    dragTargetRef.current = null;
    repaintOverlay();
  }, [repaintOverlay]);

  const handleMouseUp = useCallback(() => {
    dragTargetRef.current = null;
  }, []);

  // Set cursor based on hover target
  const handleMouseHover = useCallback(
    (e: React.MouseEvent<HTMLCanvasElement>) => {
      if (!isDipoleSource(sourceTypeRef.current) || dragTargetRef.current)
        return;
      const viewSize = lastViewSizeRef.current;
      if (!viewSize) return;

      const rect = canvasRef.current!.getBoundingClientRect();
      const mx = e.clientX - rect.left;
      const my = e.clientY - rect.top;
      const d = dipoleRef.current;
      const { cx, cy } = physToCanvas(d.xs, d.ys, viewSize);

      let cursor = "default";
      if (sourceTypeRef.current === "dipole_exy") {
        const armX = cx + ARM_LENGTH_PX * Math.cos(d.alpha);
        const armY = cy - ARM_LENGTH_PX * Math.sin(d.alpha);
        if (Math.hypot(mx - armX, my - armY) < 12) cursor = "grab";
      }
      if (cursor === "default" && Math.hypot(mx - cx, my - cy) < 14) {
        cursor = "move";
      }
      canvasRef.current!.style.cursor = cursor;
    },
    [physToCanvas],
  );

  const fieldSymbol = polarization === "TM" ? "E" : "H";
  const modeLabels: Record<VisualizationMode, JSX.Element> = {
    magnitude: (
      <span>
        |{fieldSymbol}
        <sub>z</sub>|
      </span>
    ),
    log_magnitude: (
      <span>
        log|{fieldSymbol}
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
        <div className="mode-pill">
          <div
            className="pill-highlight pill-highlight-5"
            style={{
              transform: `translateX(${(Object.keys(modeLabels) as VisualizationMode[]).indexOf(mode) * 100}%)`,
            }}
          />
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

      <div className="zoom-control">
        <label>Zoom:</label>
        <div className="mode-pill">
          <div
            className="pill-highlight pill-highlight-3"
            style={{ transform: `translateX(${zoomIndex * 100}%)` }}
          />
          {ZOOM_OPTIONS.map((opt, i) => (
            <button
              key={i}
              className={i === zoomIndex ? "active" : ""}
              onClick={() => {
                setZoomIndex(i);
                onZoomChange(opt.viewSize);
              }}
            >
              {opt.label}
            </button>
          ))}
        </div>
      </div>

      <div className="canvas-container">
        <canvas
          ref={canvasRef}
          width={width}
          height={height}
          onMouseDown={handleMouseDown}
          onMouseMove={(e) => {
            handleMouseMove(e);
            handleMouseHover(e);
          }}
          onMouseUp={handleMouseUp}
          onMouseEnter={handleMouseEnter}
          onMouseLeave={handleMouseLeave}
        />
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

function drawDipoleMarker(
  ctx: CanvasRenderingContext2D,
  dipole: DipoleParams,
  sourceType: SourceType,
  viewSize: number,
  width: number,
  height: number,
) {
  if (!viewSize) return;

  const scale = width / viewSize;
  const cx = width / 2 + dipole.xs * scale;
  const cy = height / 2 - dipole.ys * scale;

  // Dipole position dot
  ctx.fillStyle = "#ff4444";
  ctx.beginPath();
  ctx.arc(cx, cy, 5, 0, 2 * Math.PI);
  ctx.fill();

  ctx.strokeStyle = "rgba(255, 255, 255, 0.9)";
  ctx.lineWidth = 1.5;
  ctx.beginPath();
  ctx.arc(cx, cy, 5, 0, 2 * Math.PI);
  ctx.stroke();

  // Orientation arm (only for DipoleExy)
  if (sourceType === "dipole_exy") {
    const armX = cx + ARM_LENGTH_PX * Math.cos(dipole.alpha);
    const armY = cy - ARM_LENGTH_PX * Math.sin(dipole.alpha); // y inverted

    ctx.strokeStyle = "#ff4444";
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(cx, cy);
    ctx.lineTo(armX, armY);
    ctx.stroke();

    // Arrowhead
    const headLen = 7;
    const angle = Math.atan2(-(armY - cy), armX - cx);
    ctx.fillStyle = "#ff4444";
    ctx.beginPath();
    ctx.moveTo(armX, armY);
    ctx.lineTo(
      armX - headLen * Math.cos(angle - 0.4),
      armY + headLen * Math.sin(angle - 0.4),
    );
    ctx.lineTo(
      armX - headLen * Math.cos(angle + 0.4),
      armY + headLen * Math.sin(angle + 0.4),
    );
    ctx.closePath();
    ctx.fill();
  }
}
