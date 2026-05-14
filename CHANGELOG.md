# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/bede/skope/compare/0.2.0...HEAD
[0.2.0]: https://github.com/bede/skope/compare/0.1.0...0.2.0
