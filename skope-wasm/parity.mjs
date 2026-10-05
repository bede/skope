// Compiled wasm (pkg-node) against the release CLI: node skope-wasm/parity.mjs
import { strict as assert } from "node:assert";
import { spawnSync } from "node:child_process";
import { closeSync, fstatSync, mkdtempSync, openSync, readdirSync, readFileSync, readSync, statSync, writeFileSync } from "node:fs";
import { tmpdir } from "node:os";
import { basename, dirname, join } from "node:path";
import { fileURLToPath } from "node:url";
import { gzipSync, zstdDecompressSync } from "node:zlib";
import { Query, matesName, sampleName } from "./pkg-node/skope_wasm.js";

const root = join(dirname(fileURLToPath(import.meta.url)), "..");
const skope = join(root, "target/release/skope");
const tmp = mkdtempSync(join(tmpdir(), "skope-wasm-parity-"));
const targets = join(root, "data/zmrp21.viruses.fa");
const targetDir = join(root, "data/target-dir");
const dirFiles = readdirSync(targetDir).filter((f) => f.endsWith(".fa")).sort().map((f) => join(targetDir, f));

// Plain, gzipped and mate-split copies of the zstd sample
const fastq = zstdDecompressSync(readFileSync(join(root, "data/rsviruses17900.1k.fastq.zst"))).toString();
const records = fastq.trimEnd().split("\n").reduce((acc, line, i) => (i % 4 ? (acc[acc.length - 1] += `\n${line}`) : acc.push(line), acc), []);
const write = (name, data) => (writeFileSync(join(tmp, name), data), join(tmp, name));
const sample = write("sample.fq", fastq);
const sampleGz = write("sample.fq.gz", gzipSync(fastq));
const targetsGz = write("targets.fa.gz", gzipSync(readFileSync(targets)));
const r1 = write("mates_R1.fq", records.filter((_, i) => i % 2 === 0).map((r) => `${r}\n`).join(""));
const r2 = write("mates_R2.fq", records.filter((_, i) => i % 2 === 1).map((r) => `${r}\n`).join(""));
const interleaved = write("interleaved.fq", records.slice(0, records.length & ~1).map((r) => `${r}\n`).join(""));

function reader(path) {
  const fd = openSync(path);
  const size = fstatSync(fd).size;
  const read = (offset, length) => {
    const chunk = Buffer.allocUnsafe(Math.max(0, Math.min(length, size - offset)));
    readSync(fd, chunk, 0, chunk.length, offset);
    return chunk;
  };
  return { read, close: () => closeSync(fd) };
}

function stream(path, use) {
  const { read, close } = reader(path);
  try {
    return use(read);
  } finally {
    close();
  }
}

const isDir = (path) => statSync(path).isDirectory();
const sampleNameOf = (s) => s.files.length === 2 ? matesName(basename(s.files[0])) : sampleName(basename(s.files[0]));

function wasm({ targets: t, index, background = [], samples, options, saveIndex = false }) {
  const query = new Query(options);
  try {
    if (index) query.loadIndex(readFileSync(index));
    for (const path of index ? [] : isDir(t) ? dirFiles : [t]) {
      stream(path, (read) => query.addTargets(sampleName(basename(path)), read));
    }
    for (const path of background) stream(path, (read) => query.addBackground(read));
    const bytes = saveIndex ? query.index() : undefined;
    query.prepare();
    for (const s of samples) {
      const [a, b] = s.files.map(reader);
      try {
        query.addSample(sampleNameOf(s), a.read, b?.read, !!s.interleaved, () => {});
      } finally {
        a.close();
        b?.close();
      }
    }
    return { tsv: query.tsv(), index: bytes };
  } finally {
    query.free();
  }
}

function flags(options = {}) {
  const args = [];
  if (options.kmer) args.push("-k", options.kmer);
  if (options.smer) args.push("-s", options.smer);
  if (options.allKmers) args.push("--all-kmers");
  if (options.individual) args.push("-i");
  if (options.fraction) args.push("-f", options.fraction);
  if (options.complexity) args.push("--complexity", options.complexity);
  return args;
}

function cli({ targets: t, index, background = [], samples, options = {} }) {
  const out = join(tmp, "out.tsv");
  const args = ["query", "-q", "-t", "1", "-o", out, ...flags(options)];
  if (options.discriminatory) args.push("-d");
  if (options.confidence) args.push("-c");
  if (options.abundanceThresholds) args.push("-a", options.abundanceThresholds.join(","));
  if (options.sort) args.push("--sort", options.sort);
  if (options.noTotal) args.push("--no-total");
  if (options.limit) args.push("-l", options.limit);
  if (samples.some((s) => s.interleaved)) args.push("--interleaved");
  for (const path of background) args.push("-b", path);
  args.push(index ?? t, ...samples.map((s) => s.files.join(",")));
  run(args);
  return readFileSync(out, "utf8");
}

function cliIndex({ targets: t, background = [], options = {} }) {
  const out = join(tmp, "cli.sk");
  const args = ["index", "build-query", "-q", "-t", "1", "-o", out, ...flags(options)];
  if (options.confidence) args.push("-p");
  for (const path of background) args.push("-b", path);
  run([...args, t]);
  return out;
}

function run(args) {
  const result = spawnSync(skope, args.map(String), { encoding: "utf8" });
  assert.equal(result.status, 0, `skope ${args.join(" ")}\n${result.stderr}`);
}

const single = [{ files: [sample] }];
const scenarios = {
  "defaults": { targets, samples: single },
  "individual discriminatory confidence": {
    targets, samples: single,
    options: { individual: true, discriminatory: true, confidence: true, abundanceThresholds: [2, 10], sort: "target" },
  },
  "all k-mers, u128, fraction, complexity": {
    targets, samples: single,
    options: { individual: true, kmer: 35, allKmers: true, fraction: 0.5, complexity: 0.3, noTotal: true, sort: "input" },
  },
  "gzip targets and sample": { targets: targetsGz, samples: [{ files: [sampleGz] }], options: { individual: true } },
  "directory of targets": { targets: targetDir, samples: single, options: { confidence: true } },
  "background masking": { targets, background: [dirFiles[0]], samples: single, options: { individual: true } },
  "paired, interleaved and single samples": {
    targets, options: { individual: true },
    samples: [{ files: [r1, r2] }, { files: [interleaved], interleaved: true }],
  },
  "base limit": { targets, samples: single, options: { individual: true, limit: 40000 } },
};

for (const [name, spec] of Object.entries(scenarios)) {
  const { tsv } = wasm(spec);
  assert.equal(tsv, cli(spec), name);
  assert.ok(tsv.split("\n").length > 3, name);
  console.log(`ok ${name}`);
}

// Indexes: wasm writes the CLI's bytes, and each reads the other's
const indexSpec = { targets, background: [dirFiles[0]], samples: single, options: { individual: true, confidence: true, fraction: 0.5 } };
const cliSk = cliIndex(indexSpec);
const { index, tsv } = wasm({ ...indexSpec, saveIndex: true });
assert.deepEqual(Buffer.from(index), readFileSync(cliSk), "index bytes");
const wasmSk = write("wasm.sk", index);
const fromIndex = { samples: single, options: { confidence: true } };
assert.equal(wasm({ ...fromIndex, index: cliSk }).tsv, cli({ ...fromIndex, index: wasmSk }), "index query");
assert.equal(wasm({ ...fromIndex, index: wasmSk }).tsv, tsv, "index matches targets");
console.log("ok index round trip");

assert.throws(() => new Query({ kmer: 30 }), /k=30/);
assert.throws(() => new Query({ allKmers: true, smer: 9 }), /conflicts/);
assert.throws(() => new Query({ sort: "size" }), /sort order/);
assert.throws(() => new Query({}).addTargets("x", () => 7), /must return a Uint8Array/);
assert.throws(() => new Query({}).addTargets("x", () => { throw new Error("gone"); }), /Read failed: gone/);
console.log("compiled wasm matches the CLI");
