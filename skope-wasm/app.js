// Page logic for the skope-wasm demo. Queries run in worker.js
import { MSG, PHASE } from "./protocol.js?v=dev";
import { SEQ_RE, describeGroups, pairSequenceFiles } from "./pairing.js?v=dev";

// Source commit, stamped by deploy.sh like every ?v= cache token
const BUILD = "dev";

const MAX_ROWS = 2000;
const INDEX_RE = /\.sk$/i;
const UNSUPPORTED_RE = /\.(zst|xz|bz2)$/i;
const PHASE_LABELS = {
  [PHASE.TARGETS]: "Collecting target k-mers",
  [PHASE.BACKGROUND]: "Masking background",
  [PHASE.SAMPLES]: "Counting sample",
};

const $ = (id) => document.getElementById(id);
const runBtn = $("run-btn");
const fileListEl = $("file-list");

// --- State -------------------------------------------------------------------
let worker = null;
let ready = false;
let running = false;
let targets = []; // target files, or a lone query index
let background = [];
let samples = []; // work groups: single file or paired R1/R2
let rows = []; // { li, st } per sample group
let run = null; // byte offsets and the active item of the current run
let outputs = {}; // downloadable Blobs

// --- Helpers -----------------------------------------------------------------
function setStatus(msg, kind) {
  $("status").textContent = msg;
  $("status").className = kind || "";
}

function humanBytes(n) {
  const units = ["B", "KB", "MB", "GB", "TB"];
  let i = 0;
  while (n >= 1024 && i < units.length - 1) {
    n /= 1024;
    i++;
  }
  return `${n.toFixed(i === 0 ? 0 : 1)}${units[i]}`;
}

const humanSeconds = (s) => (s < 1 ? `${Math.round(s * 1000)} ms` : `${s.toFixed(1)} s`);
const plural = (n, word) => `${n} ${word}${n === 1 ? "" : "s"}`;
const relativePath = (file) => file.webkitRelativePath || file.name;
const totalSize = (files) => files.reduce((sum, file) => sum + (file.size || 0), 0);

function setButtonDisabled(btn, disabled) {
  btn.disabled = disabled;
  btn.classList.toggle("pulse-glow", !disabled);
}

// What the user still has to do
function readyStatus() {
  if (!ready) return "Loading wasm…";
  if (!targets.length) return "Select targets to begin";
  if (!samples.length) return "Select samples to query against the targets";
  return "Press 'Run' to begin";
}

function update() {
  setButtonDisabled(runBtn, !ready || running || !targets.length || !samples.length);
  $("smer").disabled = $("all-kmers").checked;
}

// An error naming files the browser build cannot decompress, if any
function unsupported(files) {
  const bad = files.filter((file) => UNSUPPORTED_RE.test(file.name));
  if (!bad.length) return "";
  const more = bad.length > 1 ? ` and ${bad.length - 1} more` : "";
  return `${bad[0].name}${more}: zstd, xz and bzip2 are unsupported in the browser. Use plain or gzipped fasta/fastq`;
}

function setZone(kind, files, summary, error) {
  $(`${kind}-zone`).classList.toggle("loaded", files.length > 0);
  $(`${kind}-summary`).textContent = summary;
  $(`${kind}-info`)?.replaceChildren();
  setStatus(error || readyStatus(), error ? "error" : "");
  update();
}

// --- Selection ---------------------------------------------------------------
function setTargets(files) {
  const indexes = files.filter((file) => INDEX_RE.test(file.name));
  const seqs = files.filter((file) => SEQ_RE.test(file.name));
  let error = unsupported(files);
  if (!error && indexes.length && files.length > 1) error = "Drop a query index (.sk) on its own, as it replaces target files";
  targets = error ? [] : indexes.length ? indexes : seqs.sort((a, b) => (relativePath(a) < relativePath(b) ? -1 : 1));

  let summary = "";
  if (indexes.length && targets.length) summary = `${targets[0].name} · query index · ${humanBytes(targets[0].size)}`;
  else if (targets.length === 1) summary = `${targets[0].name} · ${humanBytes(targets[0].size)}`;
  else if (targets.length) summary = `${targets.length} target files, one target each · ${humanBytes(totalSize(targets))}`;
  else if (files.length && !error) summary = "No fasta/fastq files or query index found";
  setZone("targets", targets, summary, error);
}

function setSamples(files) {
  let error = unsupported(files);
  try {
    samples = error ? [] : pairSequenceFiles(files);
  } catch (err) {
    samples = [];
    error = err.message;
  }
  const size = samples.reduce((sum, group) => sum + group.size, 0);
  let summary = "";
  if (samples.length) summary = `${plural(samples.length, "sample")} (${describeGroups(samples)}) · ${humanBytes(size)}`;
  else if (files.length && !error) summary = "No fasta/fastq files found in selection";
  setZone("samples", samples, summary, error);
  renderRows();
}

function setBackground(files) {
  const error = unsupported(files);
  background = error ? [] : files.filter((file) => SEQ_RE.test(file.name));
  let summary = "";
  if (background.length) summary = `${plural(background.length, "background file")} · ${humanBytes(totalSize(background))}`;
  else if (files.length && !error) summary = "No fasta/fastq files found";
  setZone("background", background, summary, error);
}

function renderRows() {
  fileListEl.replaceChildren();
  rows = samples.map((group) => {
    const li = document.createElement("li");
    const name = Object.assign(document.createElement("span"), { className: "name", textContent: group.label });
    const st = Object.assign(document.createElement("span"), { className: "st" });
    li.append(name, st);
    fileListEl.append(li);
    return { li, st };
  });
  fileListEl.hidden = !rows.length;
}

// Flatten a dropped FileSystemEntry tree to File[], stamping webkitRelativePath like the picker
async function collectEntries(entries, prefix = "") {
  const out = [];
  for (const entry of entries) {
    if (entry.isFile) {
      const file = await new Promise((res, rej) => entry.file(res, rej));
      try {
        Object.defineProperty(file, "webkitRelativePath", { value: prefix + entry.name, configurable: true });
      } catch (_) {
        // Read-only in some browsers, so labels fall back to the basename
      }
      out.push(file);
    } else if (entry.isDirectory) {
      const children = await readAllEntries(entry.createReader());
      out.push(...(await collectEntries(children, prefix + entry.name + "/")));
    }
  }
  return out;
}

// readEntries returns ~100 entries per call, so loop until exhausted
function readAllEntries(reader) {
  return new Promise((resolve, reject) => {
    const all = [];
    const next = () =>
      reader.readEntries((batch) => {
        if (!batch.length) resolve(all);
        else {
          all.push(...batch);
          next();
        }
      }, reject);
    next();
  });
}

function setupZone(kind, onFiles) {
  const zone = $(`${kind}-zone`);
  const input = $(`${kind}-input`);
  zone.addEventListener("click", () => {
    if (!running) input.click();
  });
  input.addEventListener("change", () => onFiles([...input.files]));
  zone.addEventListener("dragover", (e) => {
    e.preventDefault();
    if (!running) zone.classList.add("dragover");
  });
  zone.addEventListener("dragleave", () => zone.classList.remove("dragover"));
  zone.addEventListener("drop", async (e) => {
    e.preventDefault();
    zone.classList.remove("dragover");
    if (running) return;
    const entries = [...e.dataTransfer.items].map((item) => item.webkitGetAsEntry?.()).filter(Boolean);
    onFiles(entries.length ? await collectEntries(entries) : [...e.dataTransfer.files]);
  });
  return () => {
    input.value = "";
    onFiles([]);
  };
}

// --- Options -----------------------------------------------------------------
// Bases with an optional K, M or G suffix, as `skope query --limit`
function parseBases(text) {
  const match = text.trim().match(/^(\d+(?:\.\d+)?)\s*([kmg]?)$/i);
  if (!match) throw new Error(`Invalid base limit ${text}: expected bases such as 50M`);
  return Math.round(Number(match[1]) * { "": 1, k: 1e3, m: 1e6, g: 1e9 }[match[2].toLowerCase()]);
}

function readOptions() {
  const allKmers = $("all-kmers").checked;
  const thresholds = $("thresholds").value.split(",").map((t) => t.trim()).filter(Boolean).map(Number);
  if (thresholds.some((t) => !Number.isInteger(t) || t < 0)) throw new Error("Abundance thresholds must be whole numbers, e.g. 2,10");
  return {
    kmer: Number($("kmer").value),
    smer: allKmers ? undefined : Number($("smer").value),
    allKmers,
    individual: $("individual").checked,
    discriminatory: $("discriminatory").checked,
    confidence: $("confidence").checked,
    abundanceThresholds: thresholds,
    sort: $("sort").value,
    noTotal: !$("total").checked,
    fraction: Number($("fraction").value),
    complexity: Number($("complexity").value),
    limit: $("limit").value.trim() ? parseBases($("limit").value) : undefined,
  };
}

// --- Results -----------------------------------------------------------------
function renderResults(tsv) {
  const [header, ...lines] = tsv.trimEnd().split("\n").map((line) => line.split("\t"));
  const table = document.createElement("table");
  const head = table.createTHead().insertRow();
  for (const name of header) head.append(Object.assign(document.createElement("th"), { textContent: name }));
  const body = table.createTBody();
  for (const line of lines.slice(0, MAX_ROWS)) {
    const tr = body.insertRow();
    if (line[0] === "TOTAL") tr.className = "total";
    line.forEach((value, i) => {
      const td = tr.insertCell();
      td.textContent = value;
      if (header[i] === "containment1") {
        td.className = "bar";
        td.style.setProperty("--fill", `${Number(value) * 100}%`);
      }
    });
  }
  $("table-wrap").replaceChildren(table);
  $("table-note").hidden = lines.length <= MAX_ROWS;
  $("table-note").textContent = `Showing ${MAX_ROWS} of ${lines.length} rows. Download the TSV for all`;
  $("results").hidden = false;
}

function downloadBlob(blob, filename) {
  const url = URL.createObjectURL(blob);
  const a = Object.assign(document.createElement("a"), { href: url, download: filename });
  document.body.append(a);
  a.click();
  a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

// --- Running -----------------------------------------------------------------
// Input bytes before each item, phase by phase, for overall progress
function byteOffsets() {
  const sizes = {
    [PHASE.TARGETS]: INDEX_RE.test(targets[0].name) ? [] : targets.map((file) => file.size),
    [PHASE.BACKGROUND]: background.map((file) => file.size),
    [PHASE.SAMPLES]: samples.map((group) => group.size),
  };
  let total = 0;
  const offsets = {};
  for (const [phase, list] of Object.entries(sizes)) {
    offsets[phase] = list.map((size) => (total += size) - size);
  }
  return { sizes, offsets, total };
}

function onProgress({ phase, item, bytes, bases }) {
  run.active = { phase, item };
  const done = run.offsets[phase][item] + bytes;
  const pct = run.total ? Math.min(100, (done / run.total) * 100) : 0;
  const count = run.sizes[phase].length;
  $("progress").value = pct;
  $("progress-label").textContent =
    `${PHASE_LABELS[phase]} ${count > 1 ? `${item + 1} of ${count} ` : ""}· ` +
    `${humanBytes(done)} / ${humanBytes(run.total)} (${pct.toFixed(1)}%)`;
  if (phase === PHASE.SAMPLES) {
    const size = run.sizes[phase][item];
    rows[item].li.className = "active";
    rows[item].st.textContent = `counting ${size ? Math.min(100, (bytes / size) * 100).toFixed(0) : 100}%`;
  }
}

function onMessage({ data }) {
  switch (data.type) {
    case MSG.READY:
      ready = true;
      setStatus(readyStatus());
      break;
    case MSG.PROGRESS:
      onProgress(data);
      break;
    case MSG.TARGETS:
      $("targets-info").textContent = data.summary;
      break;
    case MSG.BACKGROUND:
      $("background-info").textContent += `${data.item ? "\n" : ""}${data.summary}`;
      break;
    case MSG.SAMPLE: {
      const found = /found (\d+) of (\d+)/.exec(data.summary);
      rows[data.item].li.className = "done";
      rows[data.item].st.textContent = found
        ? `done · ${Number(found[1]).toLocaleString()} of ${Number(found[2]).toLocaleString()} k-mers found`
        : "done";
      break;
    }
    case MSG.RESULT:
      running = false;
      outputs.tsv = new Blob([data.tsv], { type: "text/tab-separated-values" });
      $("progress").value = 100;
      $("progress-label").textContent = `Processed ${humanBytes(run.total)} across ${describeGroups(samples)}`;
      setStatus(`Done: ${plural(samples.length, "sample")} queried in ${humanSeconds(data.seconds)}`, "success");
      renderResults(data.tsv);
      break;
    case MSG.ERROR:
      running = false;
      if (run?.active?.phase === PHASE.SAMPLES) {
        rows[run.active.item].li.className = "failed";
        rows[run.active.item].st.textContent = "failed";
      }
      setStatus(`Error: ${data.message}`, "error");
      break;
  }
  update();
}

function startWorker() {
  worker?.terminate();
  ready = running = false;
  worker = new Worker(`./worker.js?v=${BUILD}`, { type: "module" });
  worker.onmessage = onMessage;
  worker.onerror = (e) => {
    e.preventDefault();
    setStatus(`Error: ${e.message || "the worker failed to load"}. Build pkg/ with ./build.sh`, "error");
  };
  update();
}

runBtn.addEventListener("click", () => {
  let options;
  try {
    options = readOptions();
  } catch (err) {
    setStatus(err.message, "error");
    return;
  }
  running = true;
  run = { ...byteOffsets(), active: null };
  outputs = {};
  $("results").hidden = true;
  $("targets-info").replaceChildren();
  $("background-info").replaceChildren();
  $("progress-wrap").hidden = false;
  $("progress").value = 0;
  $("progress-label").textContent = "";
  rows.forEach(({ li, st }) => {
    li.className = "";
    st.textContent = "queued";
  });
  setStatus("Running…");
  worker.postMessage({
    type: MSG.RUN,
    targets,
    background,
    samples: samples.map((group) => (group.kind === "paired" ? [group.file1, group.file2] : [group.file])),
    interleaved: $("interleaved").checked,
    options,
  });
  update();
});

$("download-tsv").addEventListener("click", () => downloadBlob(outputs.tsv, "skope.tsv"));
$("all-kmers").addEventListener("change", update);

const clears = [setupZone("targets", setTargets), setupZone("samples", setSamples), setupZone("background", setBackground)];

// Terminating the worker also cancels a run
$("reset-btn").addEventListener("click", () => {
  startWorker();
  clears.forEach((clear) => clear());
  run = null;
  outputs = {};
  $("results").hidden = true;
  $("progress-wrap").hidden = true;
});

console.log(`skope-wasm ${BUILD}`);
startWorker();
