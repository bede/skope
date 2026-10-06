// Messages between the page and worker.js
export const MSG = Object.freeze({
  READY: "ready",
  RUN: "run",
  PROGRESS: "progress",
  TARGETS: "targets",
  BACKGROUND: "background",
  SAMPLE: "sample",
  RESULT: "result",
  ERROR: "error",
});

// Phases of a run, in order
export const PHASE = Object.freeze({
  TARGETS: "targets",
  BACKGROUND: "background",
  SAMPLES: "samples",
});
