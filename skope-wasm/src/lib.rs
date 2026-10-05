//! WebAssembly bindings for `skope query`. Inputs are read synchronously through JS callbacks,
//! as `FileReaderSync` allows in a worker, so fastx parsing and decompression match the CLI's
use js_sys::{Function, Uint8Array};
use serde::Deserialize;
use skope::{
    Input, ProcessingStats, QueryOptions, QuerySession, derive_mates_name, derive_sample_name,
};
use std::cell::{Cell, RefCell};
use std::io::{self, BufReader, Read};
use std::path::Path;
use std::rc::Rc;
use wasm_bindgen::prelude::*;

/// Bytes requested per JS read, amortising call overhead
const READ_BYTES: usize = 1 << 20;

/// Query options, named after `skope query` flags. Omitted options take the CLI's defaults
#[derive(Default, Deserialize)]
#[serde(default, rename_all = "camelCase")]
struct Options {
    kmer: Option<u8>,
    smer: Option<u8>,
    all_kmers: bool,
    individual: bool,
    discriminatory: bool,
    confidence: bool,
    abundance_thresholds: Option<Vec<usize>>,
    sort: Option<String>,
    no_total: bool,
    fraction: Option<f64>,
    complexity: Option<f32>,
    limit: Option<u64>,
}

impl Options {
    fn into_query(self) -> Result<QueryOptions, JsError> {
        if self.all_kmers && self.smer.is_some() {
            return Err(JsError::new("allKmers conflicts with smer"));
        }
        let defaults = QueryOptions::default();
        Ok(QueryOptions {
            kmer_length: self.kmer.unwrap_or(defaults.kmer_length),
            smer_length: match self.all_kmers {
                true => 0,
                false => self.smer.unwrap_or(defaults.smer_length),
            },
            individual: self.individual,
            discriminatory: self.discriminatory,
            confidence: self.confidence,
            abundance_thresholds: self
                .abundance_thresholds
                .unwrap_or(defaults.abundance_thresholds),
            sort_order: match self.sort {
                Some(sort) => sort.parse().map_err(js_error)?,
                None => defaults.sort_order,
            },
            no_total: self.no_total,
            fraction: self.fraction.unwrap_or(defaults.fraction),
            complexity: self.complexity.unwrap_or(defaults.complexity),
            limit_bp: self.limit,
        })
    }
}

fn js_error(err: impl std::fmt::Display) -> JsError {
    JsError::new(&format!("{err:#}"))
}

/// Bytes read and the first read failure, shared by one call's sources
#[derive(Default)]
struct Reads {
    bytes: Cell<f64>,
    failure: RefCell<Option<String>>,
}

impl Reads {
    /// A read failure in place of the error it caused, which can hide it
    fn error(&self, err: impl std::fmt::Display) -> JsError {
        match self.failure.take() {
            Some(failure) => JsError::new(&failure),
            None => js_error(err),
        }
    }
}

/// Sequential reads through a JS `read(offset, length) -> Uint8Array` callback, which returns
/// fewer bytes only at the end
struct Source {
    read: Function,
    offset: f64,
    reads: Rc<Reads>,
}

impl Source {
    fn chunk(&self, length: usize) -> Result<Uint8Array, String> {
        let (offset, length) = (self.offset.into(), (length as f64).into());
        let chunk = self.read.call2(&JsValue::NULL, &offset, &length);
        let chunk = chunk.map_err(|err| format!("Read failed: {}", js_message(&err)))?;
        chunk
            .dyn_into()
            .map_err(|_| "Read must return a Uint8Array".to_string())
    }
}

impl Read for Source {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        let chunk = self.chunk(buf.len()).map_err(|failure| {
            let mut first = self.reads.failure.borrow_mut();
            io::Error::other(first.get_or_insert(failure).clone())
        })?;
        let n = (chunk.length() as usize).min(buf.len());
        chunk.subarray(0, n as u32).copy_to(&mut buf[..n]);
        self.offset += n as f64;
        self.reads.bytes.set(self.reads.bytes.get() + n as f64);
        Ok(n)
    }
}

fn js_message(value: &JsValue) -> String {
    match value.dyn_ref::<js_sys::Error>() {
        Some(err) => err.message().into(),
        None => value.as_string().unwrap_or_else(|| format!("{value:?}")),
    }
}

fn source(read: Function, reads: &Rc<Reads>) -> BufReader<Source> {
    let reads = Rc::clone(reads);
    BufReader::with_capacity(
        READ_BYTES,
        Source {
            read,
            offset: 0.0,
            reads,
        },
    )
}

/// Per-batch `progress(records, bases, bytes)` calls, ignoring their errors
fn reporter(progress: Option<Function>, reads: &Rc<Reads>) -> impl FnMut(&ProcessingStats) {
    let reads = Rc::clone(reads);
    move |stats| {
        if let Some(progress) = &progress {
            let (seqs, bp) = (stats.total_seqs as f64, stats.total_bp as f64);
            let bytes = reads.bytes.get();
            let _ = progress.call3(&JsValue::NULL, &seqs.into(), &bp.into(), &bytes.into());
        }
    }
}

/// A query session. Add target files or load an index, optionally add background files, then
/// add samples and take the TSV. Progress callbacks receive running records, bases and bytes
#[wasm_bindgen]
pub struct Query {
    session: QuerySession,
}

#[wasm_bindgen]
impl Query {
    #[wasm_bindgen(constructor)]
    pub fn new(options: JsValue) -> Result<Query, JsError> {
        console_error_panic_hook::set_once();
        let options: Options = match options.is_undefined() || options.is_null() {
            true => Options::default(),
            false => serde_wasm_bindgen::from_value(options)?,
        };
        let session = QuerySession::new(options.into_query()?).map_err(js_error)?;
        Ok(Query { session })
    }

    /// K-mer length, an index's once loaded
    #[wasm_bindgen(getter)]
    pub fn kmer(&self) -> u8 {
        self.session.kmer_length()
    }

    /// S-mer length, an index's once loaded, 0 for all k-mers
    #[wasm_bindgen(getter)]
    pub fn smer(&self) -> u8 {
        self.session.smer_length()
    }

    /// Add a target file. One acts as a file, several as a directory of files
    #[wasm_bindgen(js_name = addTargets)]
    pub fn add_targets(
        &mut self,
        name: &str,
        read: Function,
        progress: Option<Function>,
    ) -> Result<(), JsError> {
        let reads = Rc::default();
        let (stream, report) = (source(read, &reads), reporter(progress, &reads));
        self.session
            .add_targets(name, stream, report)
            .map_err(|err| reads.error(err))
    }

    /// Load a prebuilt query index (.sk) in place of targets
    #[wasm_bindgen(js_name = loadIndex)]
    pub fn load_index(&mut self, bytes: &[u8]) -> Result<(), JsError> {
        self.session.load_index(bytes).map_err(js_error)
    }

    /// Drop target k-mers found in a background file. Returns its summary
    #[wasm_bindgen(js_name = addBackground)]
    pub fn add_background(
        &mut self,
        read: Function,
        progress: Option<Function>,
    ) -> Result<String, JsError> {
        let reads = Rc::default();
        let (stream, report) = (source(read, &reads), reporter(progress, &reads));
        self.session
            .add_background(stream, report)
            .map_err(|err| reads.error(err))
    }

    /// Query index (.sk) of the targets, as `skope index build-query` writes. Precedes samples
    pub fn index(&mut self) -> Result<Vec<u8>, JsError> {
        self.session.index().map_err(js_error)
    }

    /// Ready targets for samples. Returns the targets summary
    pub fn prepare(&mut self) -> Result<String, JsError> {
        self.session.prepare().map_err(js_error)
    }

    /// Add a sample from one file, interleaved mates, or two mate files. Returns its summary
    #[wasm_bindgen(js_name = addSample)]
    pub fn add_sample(
        &mut self,
        name: &str,
        read: Function,
        read2: Option<Function>,
        interleaved: bool,
        progress: Option<Function>,
    ) -> Result<String, JsError> {
        let reads = Rc::default();
        let input = match (read2, interleaved) {
            (Some(read2), _) => Input::Paired(source(read, &reads), source(read2, &reads)),
            (None, true) => Input::Interleaved(source(read, &reads)),
            (None, false) => Input::File(source(read, &reads)),
        };
        let report = reporter(progress, &reads);
        self.session
            .add_sample(name, input, report)
            .map_err(|err| reads.error(err))
    }

    /// Results of every sample so far, as `skope query` writes
    pub fn tsv(&self) -> Result<String, JsError> {
        self.session.tsv().map_err(js_error)
    }
}

/// Target or sample name of a file, as the CLI derives it
#[wasm_bindgen(js_name = sampleName)]
pub fn sample_name(file_name: &str) -> String {
    derive_sample_name(Path::new(file_name), false)
}

/// Sample name of mates, from the R1 file name less its mate tag
#[wasm_bindgen(js_name = matesName)]
pub fn mates_name(r1_file_name: &str) -> String {
    derive_mates_name(Path::new(r1_file_name))
}
