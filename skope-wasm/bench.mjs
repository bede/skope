// Single-threaded wasm throughput: node skope-wasm/bench.mjs <targets> <sample> [options-json]
import { closeSync, fstatSync, openSync, readSync } from "node:fs";
import { basename } from "node:path";
import { Query, sampleName } from "./pkg-node/skope_wasm.js";

const [targets, sample, json = "{}"] = process.argv.slice(2);
if (!sample) throw new Error("usage: node bench.mjs <targets> <sample> [options-json]");

function timed(label, path, use) {
  const fd = openSync(path);
  const size = fstatSync(fd).size;
  const read = (offset, length) => {
    const chunk = Buffer.allocUnsafe(Math.max(0, Math.min(length, size - offset)));
    readSync(fd, chunk, 0, chunk.length, offset);
    return chunk;
  };
  let bases = 0;
  const started = performance.now();
  const summary = use(read, (_records, bp) => (bases = bp));
  const seconds = (performance.now() - started) / 1000;
  closeSync(fd);
  console.log(`${label}: ${(size / 1e6 / seconds).toFixed(0)} MB/s, ${(bases / 1e6 / seconds).toFixed(0)} Mbp/s (${seconds.toFixed(2)} s)${summary ? `\n  ${summary}` : ""}`);
}

const query = new Query(JSON.parse(json));
timed("Targets", targets, (read, progress) => query.addTargets(sampleName(basename(targets)), read, progress));
query.prepare();
timed("Sample", sample, (read, progress) => query.addSample(sampleName(basename(sample)), read, undefined, false, progress));
