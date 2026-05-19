//! Helicase-backed FASTA/FASTQ parsing with a small producer/consumer parallel layer.
//!
//! `process_parallel` runs a single helicase parser on a dedicated producer thread, batches
//! owned records into a bounded channel, and fans them out to `threads` consumers that each
//! hold their own cloned [`RecordProcessor`]. When `threads <= 1` the producer and processor
//! run inline in the caller's thread.

use anyhow::Result;
use crossbeam_channel::bounded;
use helicase::input::{FromFile, FromStdin};
use helicase::{Config, FastxParser, HelicaseParser, ParserOptions};
use parking_lot::Mutex;
use std::path::Path;
use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

const HELICASE_CONFIG: Config = ParserOptions::default().config();

/// Batch size target in bases — matches what paraseq used via `update_batch_size_in_bp`.
const BATCH_TARGET_BP: usize = 256 * 1024;

pub struct OwnedRecord {
    pub header: Vec<u8>,
    pub seq: Vec<u8>,
}

pub trait RecordProcessor: Clone + Send {
    fn process_record(&mut self, header: &[u8], seq: &[u8]) -> Result<()>;
    fn on_batch_complete(&mut self) -> Result<()>;
}

const SAMPLE_LIMIT_SENTINEL: &str = "__skope_sample_limit_reached__";

pub fn sample_limit_reached_error() -> anyhow::Error {
    anyhow::anyhow!(SAMPLE_LIMIT_SENTINEL)
}

pub fn is_sample_limit_error(err: &anyhow::Error) -> bool {
    err.chain().any(|e| e.to_string() == SAMPLE_LIMIT_SENTINEL)
}

pub fn process_parallel<P: RecordProcessor>(
    in_path: Option<&Path>,
    processor: &mut P,
    threads: usize,
) -> Result<()> {
    if threads <= 1 {
        run_inline(in_path, processor)
    } else {
        run_parallel(in_path, processor, threads)
    }
}

fn run_inline<P: RecordProcessor>(in_path: Option<&Path>, processor: &mut P) -> Result<()> {
    let mut batch_bp: usize = 0;
    let mut pending = false;

    macro_rules! drive {
        ($parser_expr:expr) => {{
            let mut parser = $parser_expr;
            while parser.next().is_some() {
                let header = parser.get_header();
                let seq = parser.get_dna_string();
                processor.process_record(header, seq)?;
                batch_bp += seq.len();
                pending = true;
                if batch_bp >= BATCH_TARGET_BP {
                    processor.on_batch_complete()?;
                    batch_bp = 0;
                    pending = false;
                }
            }
        }};
    }

    match in_path {
        Some(path) => drive!(FastxParser::<HELICASE_CONFIG>::from_file(path)?),
        None => drive!(FastxParser::<HELICASE_CONFIG>::from_stdin()?),
    }

    if pending {
        processor.on_batch_complete()?;
    }
    Ok(())
}

fn run_parallel<P: RecordProcessor>(
    in_path: Option<&Path>,
    processor: &mut P,
    threads: usize,
) -> Result<()> {
    let (tx, rx) = bounded::<Vec<OwnedRecord>>(threads * 2);
    let cancel = Arc::new(AtomicBool::new(false));
    let error: Arc<Mutex<Option<anyhow::Error>>> = Arc::new(Mutex::new(None));

    std::thread::scope(|scope| {
        let cancel_p = Arc::clone(&cancel);
        let error_p = Arc::clone(&error);
        let producer = scope.spawn(move || {
            let result = produce(in_path, &tx, &cancel_p);
            if let Err(e) = result {
                let mut slot = error_p.lock();
                if slot.is_none() {
                    *slot = Some(e);
                }
                cancel_p.store(true, Ordering::Relaxed);
            }
            // `tx` is dropped here, signalling EOF to consumers.
        });

        let mut handles = Vec::with_capacity(threads);
        for _ in 0..threads {
            let mut local = processor.clone();
            let rx = rx.clone();
            let cancel_c = Arc::clone(&cancel);
            let error_c = Arc::clone(&error);
            handles.push(scope.spawn(move || {
                while let Ok(batch) = rx.recv() {
                    if cancel_c.load(Ordering::Relaxed) {
                        break;
                    }
                    let outcome: Result<()> = (|| {
                        for rec in &batch {
                            local.process_record(&rec.header, &rec.seq)?;
                        }
                        local.on_batch_complete()?;
                        Ok(())
                    })();
                    if let Err(e) = outcome {
                        let mut slot = error_c.lock();
                        if slot.is_none() {
                            *slot = Some(e);
                        }
                        cancel_c.store(true, Ordering::Relaxed);
                        break;
                    }
                }
            }));
        }

        drop(rx);

        producer.join().expect("fastx producer panicked");
        for h in handles {
            h.join().expect("fastx consumer panicked");
        }
    });

    if let Some(e) = error.lock().take() {
        Err(e)
    } else {
        Ok(())
    }
}

fn produce(
    in_path: Option<&Path>,
    tx: &crossbeam_channel::Sender<Vec<OwnedRecord>>,
    cancel: &AtomicBool,
) -> Result<()> {
    let mut batch: Vec<OwnedRecord> = Vec::new();
    let mut batch_bp: usize = 0;

    macro_rules! drive {
        ($parser_expr:expr) => {{
            let mut parser = $parser_expr;
            while parser.next().is_some() {
                if cancel.load(Ordering::Relaxed) {
                    return Ok(());
                }
                let header = parser.get_header();
                let seq = parser.get_dna_string();
                batch.push(OwnedRecord {
                    header: header.to_vec(),
                    seq: seq.to_vec(),
                });
                batch_bp += seq.len();
                if batch_bp >= BATCH_TARGET_BP {
                    let to_send = std::mem::take(&mut batch);
                    batch_bp = 0;
                    if tx.send(to_send).is_err() {
                        return Ok(());
                    }
                }
            }
        }};
    }

    match in_path {
        Some(path) => drive!(FastxParser::<HELICASE_CONFIG>::from_file(path)?),
        None => drive!(FastxParser::<HELICASE_CONFIG>::from_stdin()?),
    }

    if !batch.is_empty() {
        let _ = tx.send(batch);
    }
    Ok(())
}
