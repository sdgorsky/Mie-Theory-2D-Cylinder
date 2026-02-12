import type { JSX } from "react";
import { useRef, useCallback, useState } from "react";
import type { ComplexObj } from "../types/cylinder";
import type { ParameterBounds } from "../hooks/useScattering";
import "./ComplexPlane.css";

interface ComplexPlaneProps {
  permittivity: ComplexObj;
  permeability: ComplexObj;
  onPermittivityChange: (v: ComplexObj) => void;
  onPermeabilityChange: (v: ComplexObj) => void;
  bounds: ParameterBounds;
}

type DragTarget = "permittivity" | "permeability" | null;

const PADDING = 32; // px around the plot area for labels/ticks
const POINT_RADIUS = 8;

function snap(value: number): number {
  return Math.round(value * 10) / 10;
}

function clamp(value: number, min: number, max: number): number {
  return Math.max(min, Math.min(max, value));
}

export function ComplexPlane({
  permittivity,
  permeability,
  onPermittivityChange,
  onPermeabilityChange,
  bounds,
}: ComplexPlaneProps) {
  const svgRef = useRef<SVGSVGElement>(null);
  const [dragging, setDragging] = useState<DragTarget>(null);

  // Axis ranges: union of both material bounds
  const xMin = Math.min(bounds.permittivity_re_min, bounds.permeability_re_min);
  const xMax = Math.max(bounds.permittivity_re_max, bounds.permeability_re_max);
  const yMin = Math.min(bounds.permittivity_im_min, bounds.permeability_im_min);
  const yMax = Math.max(bounds.permittivity_im_max, bounds.permeability_im_max);

  // SVG viewBox dimensions (square)
  const size = 300;
  const plotSize = size - 2 * PADDING;

  // Data → SVG coordinate transforms
  const dataToSvgX = useCallback(
    (re: number) => PADDING + ((re - xMin) / (xMax - xMin)) * plotSize,
    [xMin, xMax, plotSize],
  );
  const dataToSvgY = useCallback(
    (im: number) => PADDING + ((yMax - im) / (yMax - yMin)) * plotSize,
    [yMin, yMax, plotSize],
  );

  // SVG → Data coordinate transforms
  const svgToDataX = useCallback(
    (sx: number) => xMin + ((sx - PADDING) / plotSize) * (xMax - xMin),
    [xMin, xMax, plotSize],
  );
  const svgToDataY = useCallback(
    (sy: number) => yMax - ((sy - PADDING) / plotSize) * (yMax - yMin),
    [yMin, yMax, plotSize],
  );

  const getSvgPoint = useCallback(
    (e: React.PointerEvent) => {
      const svg = svgRef.current;
      if (!svg) return null;
      const rect = svg.getBoundingClientRect();
      const scaleX = size / rect.width;
      const scaleY = size / rect.height;
      return {
        x: (e.clientX - rect.left) * scaleX,
        y: (e.clientY - rect.top) * scaleY,
      };
    },
    [size],
  );

  const handlePointerDown = useCallback(
    (target: DragTarget, e: React.PointerEvent) => {
      e.preventDefault();
      setDragging(target);
      (e.target as Element).setPointerCapture(e.pointerId);
    },
    [],
  );

  const handlePointerMove = useCallback(
    (e: React.PointerEvent) => {
      if (!dragging) return;
      const pt = getSvgPoint(e);
      if (!pt) return;

      const rawRe = svgToDataX(pt.x);
      const rawIm = svgToDataY(pt.y);

      if (dragging === "permittivity") {
        onPermittivityChange({
          re: snap(
            clamp(
              rawRe,
              bounds.permittivity_re_min,
              bounds.permittivity_re_max,
            ),
          ),
          im: snap(
            clamp(
              rawIm,
              bounds.permittivity_im_min,
              bounds.permittivity_im_max,
            ),
          ),
        });
      } else {
        onPermeabilityChange({
          re: snap(
            clamp(
              rawRe,
              bounds.permeability_re_min,
              bounds.permeability_re_max,
            ),
          ),
          im: snap(
            clamp(
              rawIm,
              bounds.permeability_im_min,
              bounds.permeability_im_max,
            ),
          ),
        });
      }
    },
    [
      dragging,
      getSvgPoint,
      svgToDataX,
      svgToDataY,
      bounds,
      onPermittivityChange,
      onPermeabilityChange,
    ],
  );

  const handlePointerUp = useCallback(() => {
    setDragging(null);
  }, []);

  // Generate integer grid lines
  const gridLines: JSX.Element[] = [];
  for (let x = Math.ceil(xMin); x <= Math.floor(xMax); x++) {
    const sx = dataToSvgX(x);
    const isAxis = x === 0;
    gridLines.push(
      <line
        key={`vg-${x}`}
        x1={sx}
        y1={PADDING}
        x2={sx}
        y2={PADDING + plotSize}
        stroke={isAxis ? "#666" : "#333"}
        strokeWidth={isAxis ? 1.5 : 0.5}
      />,
    );
  }
  for (let y = Math.ceil(yMin); y <= Math.floor(yMax); y++) {
    const sy = dataToSvgY(y);
    const isAxis = y === 0;
    gridLines.push(
      <line
        key={`hg-${y}`}
        x1={PADDING}
        y1={sy}
        x2={PADDING + plotSize}
        y2={sy}
        stroke={isAxis ? "#666" : "#333"}
        strokeWidth={isAxis ? 1.5 : 0.5}
      />,
    );
  }

  // Tick labels
  const tickLabels: JSX.Element[] = [];
  for (let x = Math.ceil(xMin); x <= Math.floor(xMax); x++) {
    tickLabels.push(
      <text
        key={`xt-${x}`}
        x={dataToSvgX(x)}
        y={PADDING + plotSize + 14}
        textAnchor="middle"
        fill="#888"
        fontSize={9}
      >
        {x}
      </text>,
    );
  }
  for (let y = Math.ceil(yMin); y <= Math.floor(yMax); y++) {
    tickLabels.push(
      <text
        key={`yt-${y}`}
        x={PADDING - 6}
        y={dataToSvgY(y) + 3}
        textAnchor="end"
        fill="#888"
        fontSize={9}
      >
        {y}
      </text>,
    );
  }

  const epsX = dataToSvgX(permittivity.re);
  const epsY = dataToSvgY(permittivity.im);
  const muX = dataToSvgX(permeability.re);
  const muY = dataToSvgY(permeability.im);

  return (
    <div className="complex-plane">
      <svg
        ref={svgRef}
        viewBox={`0 0 ${size} ${size}`}
        className={`complex-plane-svg${dragging ? " dragging" : ""}`}
        onPointerMove={handlePointerMove}
        onPointerUp={handlePointerUp}
        onPointerLeave={handlePointerUp}
      >
        {/* Plot background */}
        <rect
          x={PADDING}
          y={PADDING}
          width={plotSize}
          height={plotSize}
          fill="#1a1a2a"
        />

        {/* Grid lines and axes */}
        {gridLines}

        {/* Axis labels */}
        <text
          x={PADDING + plotSize / 2}
          y={size - 2}
          textAnchor="middle"
          fill="#888"
          fontSize={11}
        >
          Re
        </text>
        <text
          x={6}
          y={PADDING + plotSize / 2}
          textAnchor="middle"
          fill="#888"
          fontSize={11}
          transform={`rotate(-90, 6, ${PADDING + plotSize / 2})`}
        >
          Im
        </text>

        {/* Tick labels */}
        {tickLabels}

        {/* Permittivity point (εᵣ) — cyan */}
        <circle
          cx={epsX}
          cy={epsY}
          r={POINT_RADIUS}
          fill="#00d4ff"
          stroke="#fff"
          strokeWidth={1.5}
          className="complex-plane-point"
          onPointerDown={(e) => handlePointerDown("permittivity", e)}
        />

        {/* Permeability point (μᵣ) — orange */}
        <circle
          cx={muX}
          cy={muY}
          r={POINT_RADIUS}
          fill="#ffaa00"
          stroke="#fff"
          strokeWidth={1.5}
          className="complex-plane-point"
          onPointerDown={(e) => handlePointerDown("permeability", e)}
        />

        {/* Legend */}
        <circle cx={PADDING + 8} cy={PADDING - 16} r={5} fill="#00d4ff" />
        <text x={PADDING + 18} y={PADDING - 12} fill="#ccc" fontSize={10}>
          εᵣ
        </text>
        <circle cx={PADDING + 42} cy={PADDING - 16} r={5} fill="#ffaa00" />
        <text x={PADDING + 52} y={PADDING - 12} fill="#ccc" fontSize={10}>
          μᵣ
        </text>
      </svg>
    </div>
  );
}
