import { useState, useEffect, useCallback, useRef } from "react";
import type { MutableRefObject } from "react";
import type { ScatteringParams } from "../types/cylinder";
import type {
  ComputeParams,
  ImageStats,
  ScatteringMeta,
  VisualizationMode,
  WorkerRequest,
  WorkerResponse,
} from "../workers/scattering.worker";

export type {
  VisualizationMode,
  ImageStats,
  ScatteringMeta,
} from "../workers/scattering.worker";

// Callback for imperative canvas painting — bypasses React render cycle
export type PaintCallback = (
  imageData: Uint8ClampedArray,
  gridSize: number,
  viewSize: number,
) => void;

// Metadata-only state — no imageData, so React renders are lightweight (text only)
interface MetaState {
  scatteringResult: ScatteringMeta | null;
  computedParams: ScatteringParams | null;
  imageStats: ImageStats | null;
  hasImage: boolean;
}

const INITIAL_META: MetaState = {
  scatteringResult: null,
  computedParams: null,
  imageStats: null,
  hasImage: false,
};

export interface UseScatteringResult {
  isLoading: boolean;
  isReady: boolean;
  error: string | null;
  scatteringResult: ScatteringMeta | null;
  computedParams: ScatteringParams | null;
  imageStats: ImageStats | null;
  hasImage: boolean;
  paintRef: MutableRefObject<PaintCallback | null>;
  computeAll: (params: ScatteringParams) => void;
  recolor: (mode: VisualizationMode) => void;
}

function toComputeParams(p: ScatteringParams): ComputeParams {
  return {
    wavelength: p.wavelength,
    permittivityRe: p.material.permittivity.re,
    permittivityIm: p.material.permittivity.im,
    permeabilityRe: p.material.permeability.re,
    permeabilityIm: p.material.permeability.im,
    polarization: p.polarization === "TM" ? 0 : 1,
    maxOrder: p.maxOrder,
  };
}

export function useScattering(): UseScatteringResult {
  const [isLoading, setIsLoading] = useState(true);
  const [isReady, setIsReady] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [meta, setMeta] = useState<MetaState>(INITIAL_META);

  const workerRef = useRef<Worker | null>(null);
  const busyRef = useRef(false);
  const pendingParamsRef = useRef<ScatteringParams | null>(null);
  const pendingModeRef = useRef<VisualizationMode | null>(null);
  const modeRef = useRef<VisualizationMode>("magnitude");
  const lastSentParamsRef = useRef<ScatteringParams | null>(null);

  // Imperative paint callback — set by FieldVisualization
  const paintRef = useRef<PaintCallback | null>(null);

  // Cached constants from worker (set once, never change)
  const cachedGridSizeRef = useRef(0);
  const cachedViewSizeRef = useRef(0);
  const cachedScatteringRef = useRef<ScatteringMeta | null>(null);

  // Buffer recycling
  const recycleBufferRef = useRef<ArrayBuffer | null>(null);
  const currentImageDataRef = useRef<Uint8ClampedArray | null>(null);

  // rAF coalescing — latest result wins, one rAF at a time
  const pendingPaintRef = useRef<{
    imageData: Uint8ClampedArray;
    gridSize: number;
    viewSize: number;
    scattering: ScatteringMeta | null;
    stats: ImageStats;
  } | null>(null);
  const rafIdRef = useRef(0);

  // Helper: postMessage with recycled buffer attached
  const postToWorker = useCallback((worker: Worker, msg: WorkerRequest) => {
    if (msg.type === "init") {
      worker.postMessage(msg);
      return;
    }
    const transfers: Transferable[] = [];
    const buf = recycleBufferRef.current;
    if (buf && buf.byteLength > 0) {
      (msg as { recycledBuffer?: ArrayBuffer }).recycledBuffer = buf;
      transfers.push(buf);
      recycleBufferRef.current = null;
    }
    worker.postMessage(msg, transfers);
  }, []);

  // Schedule canvas paint + metadata update for next animation frame.
  // Coalesces: if multiple results arrive between frames, only the latest paints.
  const schedulePaint = useCallback(
    (data: {
      imageData: Uint8ClampedArray;
      gridSize: number;
      viewSize: number;
      scattering: ScatteringMeta | null;
      stats: ImageStats;
    }) => {
      pendingPaintRef.current = data;
      if (!rafIdRef.current) {
        rafIdRef.current = requestAnimationFrame(() => {
          rafIdRef.current = 0;
          const result = pendingPaintRef.current;
          if (!result) return;
          pendingPaintRef.current = null;

          // 1. Paint canvas imperatively (bypasses React entirely)
          paintRef.current?.(
            result.imageData,
            result.gridSize,
            result.viewSize,
          );

          // 2. Recycle old buffer (after paint has read it)
          const old = currentImageDataRef.current;
          if (old && old.buffer.byteLength > 0) {
            recycleBufferRef.current = old.buffer as ArrayBuffer;
          }
          currentImageDataRef.current = result.imageData;

          // 3. Update lightweight metadata state (triggers text-only React render)
          setMeta({
            scatteringResult: result.scattering,
            computedParams: lastSentParamsRef.current,
            imageStats: result.stats,
            hasImage: true,
          });
        });
      }
    },
    [],
  );

  const computeAll = useCallback(
    (params: ScatteringParams) => {
      if (!workerRef.current) return;
      if (busyRef.current) {
        pendingParamsRef.current = params;
      } else {
        pendingParamsRef.current = null;
        pendingModeRef.current = null;
        busyRef.current = true;
        lastSentParamsRef.current = params;
        postToWorker(workerRef.current, {
          type: "compute",
          params: toComputeParams(params),
          mode: modeRef.current,
        });
      }
    },
    [postToWorker],
  );

  const recolor = useCallback(
    (mode: VisualizationMode) => {
      modeRef.current = mode;
      if (!workerRef.current) return;
      if (busyRef.current) {
        if (!pendingParamsRef.current) {
          pendingModeRef.current = mode;
        }
      } else {
        pendingModeRef.current = null;
        busyRef.current = true;
        postToWorker(workerRef.current, { type: "recolor", mode });
      }
    },
    [postToWorker],
  );

  const dispatchPending = useCallback(
    (worker: Worker) => {
      const pendingParams = pendingParamsRef.current;
      if (pendingParams) {
        pendingParamsRef.current = null;
        pendingModeRef.current = null;
        busyRef.current = true;
        lastSentParamsRef.current = pendingParams;
        postToWorker(worker, {
          type: "compute",
          params: toComputeParams(pendingParams),
          mode: modeRef.current,
        });
      } else {
        const pendingMode = pendingModeRef.current;
        if (pendingMode) {
          pendingModeRef.current = null;
          busyRef.current = true;
          postToWorker(worker, { type: "recolor", mode: pendingMode });
        } else {
          busyRef.current = false;
        }
      }
    },
    [postToWorker],
  );

  useEffect(() => {
    const worker = new Worker(
      new URL("../workers/scattering.worker.ts", import.meta.url),
      { type: "module" },
    );
    workerRef.current = worker;

    worker.onmessage = (e: MessageEvent<WorkerResponse>) => {
      const msg = e.data;

      if (msg.type === "ready") {
        setIsReady(true);
        setIsLoading(false);
      } else if (msg.type === "result") {
        // Cache constants (these never change after first result)
        cachedGridSizeRef.current = msg.gridSize;
        cachedViewSizeRef.current = msg.viewSize;
        cachedScatteringRef.current = msg.scattering;

        // Schedule paint for next animation frame (non-blocking)
        schedulePaint({
          imageData: msg.imageData,
          gridSize: msg.gridSize,
          viewSize: msg.viewSize,
          scattering: msg.scattering,
          stats: msg.stats,
        });

        // Dispatch next pending work immediately (don't wait for paint)
        dispatchPending(worker);
      } else if (msg.type === "recolored") {
        // Use cached gridSize/viewSize/scattering (refs, not stale state)
        schedulePaint({
          imageData: msg.imageData,
          gridSize: cachedGridSizeRef.current,
          viewSize: cachedViewSizeRef.current,
          scattering: cachedScatteringRef.current,
          stats: msg.stats,
        });

        dispatchPending(worker);
      } else if (msg.type === "error") {
        setError(msg.message);
        busyRef.current = false;
      }
    };

    worker.onerror = (e) => {
      setError(`Worker error: ${e.message}`);
      setIsLoading(false);
    };

    worker.postMessage({ type: "init" });

    return () => {
      worker.terminate();
      workerRef.current = null;
      if (rafIdRef.current) {
        cancelAnimationFrame(rafIdRef.current);
        rafIdRef.current = 0;
      }
    };
  }, [dispatchPending, schedulePaint]);

  return {
    isLoading,
    isReady,
    error,
    scatteringResult: meta.scatteringResult,
    computedParams: meta.computedParams,
    imageStats: meta.imageStats,
    hasImage: meta.hasImage,
    paintRef,
    computeAll,
    recolor,
  };
}
