// Runs a query off the main thread. FileReaderSync lets wasm pull file slices as it parses
import init, { Query, matesName, sampleName } from "./pkg/skope_wasm.js?v=dev";
import { MSG, PHASE } from "./protocol.js?v=dev";

const PROGRESS_MS = 100;
const reader = new FileReaderSync();
const bytes = (blob) => new Uint8Array(reader.readAsArrayBuffer(blob));
const read = (file) => (offset, length) => bytes(file.slice(offset, offset + length));
const post = (type, data = {}) => self.postMessage({ type, ...data });
const isIndex = (file) => new TextDecoder().decode(bytes(file.slice(0, 4))) === "SKPE";
// Summaries lead with their stream, e.g. "Sample: ", which the page already shows
const unlabel = (summary) => summary.replace(/^(Targets|Background|Sample): /gm, "");

// Throttled progress of item `item` in `phase`
function progress(phase, item) {
  let last = 0;
  post(MSG.PROGRESS, { phase, item, bytes: 0, bases: 0 });
  return (records, bases, read) => {
    const now = performance.now();
    if (now - last < PROGRESS_MS) return;
    last = now;
    post(MSG.PROGRESS, { phase, item, records, bases, bytes: read });
  };
}

function run({ targets, background, samples, interleaved, options }) {
  const started = performance.now();
  const query = new Query(options);
  try {
    let info;
    if (targets.length === 1 && isIndex(targets[0])) {
      query.loadIndex(bytes(targets[0]));
      info = `Query index, k=${query.kmer}, ${query.smer ? `s=${query.smer}` : "all k-mers"}`;
    } else {
      targets.forEach((file, i) => {
        query.addTargets(sampleName(file.name), read(file), progress(PHASE.TARGETS, i));
      });
    }
    background.forEach((file, i) => {
      const summary = query.addBackground(read(file), progress(PHASE.BACKGROUND, i));
      post(MSG.BACKGROUND, { item: i, summary: unlabel(summary) });
    });
    const summary = unlabel(query.prepare());
    post(MSG.TARGETS, { summary: info ? `${info}\n${summary}` : summary });
    samples.forEach((files, i) => {
      const [r1, r2] = files;
      const name = r2 ? matesName(r1.name) : sampleName(r1.name);
      const report = progress(PHASE.SAMPLES, i);
      const summary = query.addSample(name, read(r1), r2 && read(r2), !r2 && interleaved, report);
      post(MSG.SAMPLE, { item: i, name, summary: unlabel(summary) });
    });
    post(MSG.RESULT, { tsv: query.tsv(), seconds: (performance.now() - started) / 1000 });
  } finally {
    query.free();
  }
}

self.onmessage = ({ data }) => {
  if (data.type !== MSG.RUN) return;
  try {
    run(data);
  } catch (err) {
    post(MSG.ERROR, { message: err?.message ?? String(err) });
  }
};

// Name the wasm too, as the glue would resolve it without a cache token
const wasm = new URL("./pkg/skope_wasm_bg.wasm?v=dev", import.meta.url);
init({ module_or_path: wasm }).then(() => post(MSG.READY), (err) => post(MSG.ERROR, { message: `Failed to load wasm: ${err.message}` }));
