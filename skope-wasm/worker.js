// Runs a query off the main thread. FileReaderSync lets wasm pull file slices as it parses
import init, { Query, matesName, sampleName } from "./pkg/skope_wasm.js";
import { MSG } from "./protocol.js";

const PROGRESS_MS = 100;
const files = new FileReaderSync();
const bytes = (blob) => new Uint8Array(files.readAsArrayBuffer(blob));
const read = (file) => (offset, length) => bytes(file.slice(offset, offset + length));
const post = (type, data = {}) => self.postMessage({ type, ...data });
const log = (line) => post(MSG.LOG, { line });

const isIndex = (file) => new TextDecoder().decode(bytes(file.slice(0, 4))) === "SKPE";

// Progress of one stage over `size` input bytes, throttled
function progress(stage, size) {
  let last = 0;
  post(MSG.PROGRESS, { stage, records: 0, bases: 0, bytes: 0, size });
  return (records, bases, read) => {
    const now = performance.now();
    if (now - last < PROGRESS_MS) return;
    last = now;
    post(MSG.PROGRESS, { stage, records, bases, bytes: read, size });
  };
}

// Samples as single files, consecutive R1/R2 files sorted by name, or interleaved files
function samples(list, layout) {
  if (layout !== "paired") {
    return list.map((file) => ({ name: sampleName(file.name), files: [file] }));
  }
  const sorted = [...list].sort((a, b) => a.name.localeCompare(b.name));
  if (sorted.length % 2) throw new Error("Paired samples need an R1 and R2 file each");
  return sorted.flatMap((file, i) => (i % 2 ? [] : [{ name: matesName(file.name), files: [file, sorted[i + 1]] }]));
}

function run({ targets, background, sampleFiles, layout, options, saveIndex }) {
  const started = performance.now();
  const list = samples(sampleFiles, layout);
  const query = new Query(options);
  try {
    if (targets.length === 1 && isIndex(targets[0])) {
      query.loadIndex(bytes(targets[0]));
      log(`Index: ${targets[0].name} (k=${query.kmer}, ${query.smer ? `s=${query.smer}` : "all k-mers"})`);
    } else {
      for (const file of targets) {
        query.addTargets(sampleName(file.name), read(file), progress(`Targets: ${file.name}`, file.size));
      }
    }
    for (const file of background) {
      log(query.addBackground(read(file), progress(`Background: ${file.name}`, file.size)));
    }
    if (saveIndex) {
      const index = query.index();
      self.postMessage({ type: MSG.INDEX, index }, [index.buffer]);
    }
    for (const line of query.prepare().split("\n")) log(line);
    for (const sample of list) {
      const [r1, r2] = sample.files;
      const size = sample.files.reduce((sum, file) => sum + file.size, 0);
      const report = progress(`Sample: ${sample.name}`, size);
      log(query.addSample(sample.name, read(r1), r2 && read(r2), layout === "interleaved", report));
    }
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

init().then(() => post(MSG.READY), (err) => post(MSG.ERROR, { message: `Failed to load wasm: ${err.message}` }));
