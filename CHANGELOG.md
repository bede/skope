# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `skope index build-query`: build a reusable query index (`.sk`) from target fastx file(s), so target syncmer extraction (and any background masking) is done once up front. `skope query` accepts such an index in place of a fastx `<TARGETS>` argument (auto-detected). `-p/--positions` bakes in syncmer positions so the index also supports `--confidence`/`--dump-syncmers`.
- `skope query` and `skope index build-query`: `-b/--background` masks syncmers shared with off-target/background sequences out of the targets, reducing false positives from sequences common across the tree of life. Background is streamed in a memory-bounded manner (peak memory stays bounded by the targets, not the background); directories are searched recursively.
- `skope index info <index.sk>`: print metadata for a query or classification index (k, s, target/group counts, syncmer totals, and — for query indexes — whether positions are stored and whether background masking was applied). Reads only the header and metadata, not the syncmer entries.

### Changed

- On-disk indexes now share the `.sk` extension for both classification and query indexes. Passing the wrong index type throws a clear error.
- `skope index build` is now `skope index build-classify` (the old `build` name remains as an alias).
- `--smer-length` is now long-only in `skope classify`
- `skope query`: removed the `--disjoint`/`--spacing` target syncmer downsampling option. Replaced with FracMinHash downsampling of selected syncmers `--fraction` (`-f`).
- `skope query`: a multi-fastx file passed as `<TARGETS>` is now merged into a single target named after the file unless `--individual` (`-i`) is passed. Directory inputs are unaffected (still one target per top-level file or subdirectory).
- `skope query`: skip the shared-syncmer counting pass when there is only one target.
- `skope query`: renamed `--dump-positions` to `--dump-syncmers` and added a third TSV column with the (canonical) k-mer sequence. Output is now `target\tposition\tkmer`. `--dump-syncmers` respects `--discriminatory`.

## [0.3.0] - 2026-05-15

### Added

- `skope classify` and `skope index build`: accept subdirectories as groups, with one group per top-level fastx file or subdirectory of fastx files. Matches `skope query` target-directory behaviour.
- Suggested `.skcl` file extension for classification indexes (documentation only; no on-disk format change).
- `skope lenhist`: per-group length histograms. The first positional argument is now a `.skcl` index or a directory of groups (same input as `classify`); each read is binned into its single matching group, `ambiguous`, or `unclassified`. Output gains a `group` column and `group_seqs`/`group_bases` totals; `seqs_with_hits`/`seqs_without_hits` are removed (derivable from the non-`unclassified` rows). New flags `-m/--min-hits`, `-r/--min-fraction`, `-d/--discriminatory` mirror `classify`. Pass `-` to disable group filtering, in which case all reads go to a single `all` bucket.

### Changed

- Target/class directory discovery is now unified across `skope query`, `skope classify`, and `skope index build`. Nested sub-subdirectories inside a group directory are now an error rather than silently ignored. Duplicate derived names (e.g. `foo.fa` alongside `foo/`) are rejected in both query and classify.
- Classification index magic bytes changed from `SKPE` to `SKCL`. Indexes built with previous versions must be rebuilt.
- Group bitmask widened from `u64` to `u128` in `skope classify` and `skope lenhist`, raising the maximum number of groups from 64 to 128. Index format version bumped to 2; previous indexes must be rebuilt.

## [0.2.0] - 2026-05-13

Default open syncmer parameters: `k=31`, `s=9`.

### Added

- `skope query`: `sample_seqs` and `sample_bases` output columns, respecting `--limit`.
- `skope query`: accept directories as targets, treating each contained fastx as a sample.
- `lenhist.py`: `--mode kde` restores KDE plotting alongside the default histogram.
- Support for process substitution as input.
- Published on crates.io.
- MIT license.

### Changed

- Disjoint syncmers are now opt-in via `--disjoint` (previously default).
- TSV is now the default output format.
- Switched serialization from `bincode` to `wincode`.
- Plot y-axis is divided by bin width.
- Improved axis-limit handling in plots; `lenhist.py` min/max length now filter data as well as axes.
- Startup message is always shown.
- Replaced remaining references to "minimizers" with "syncmers".

### Removed

- Overlapping syncmers.

[Unreleased]: https://github.com/bede/skope/compare/0.3.0...HEAD
[0.3.0]: https://github.com/bede/skope/compare/0.2.0...0.3.0
[0.2.0]: https://github.com/bede/skope/compare/0.1.0...0.2.0
