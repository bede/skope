# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- `skope classify` and `skope index build`: accept subdirectories as groups, with one group per top-level fastx file or subdirectory of fastx files. Matches `skope query` target-directory behaviour.

### Changed

- Target/class directory discovery is now unified across `skope query`, `skope classify`, and `skope index build`. Nested sub-subdirectories inside a group directory are now an error rather than silently ignored. Duplicate derived names (e.g. `foo.fa` alongside `foo/`) are rejected in both query and classify.

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
