# bwTools 0.7.0

## Large BigWig benchmark framework

- Added public `benchmark_bwg()` for reproducible profiling of large local
  BigWig files before performance optimization.
- Added package-metadata-cold and cached lazy-read benchmarks without claiming
  to flush the operating-system file cache.
- Added indexed retrieval benchmarks across user-defined or deterministic
  genomic windows.
- Added zoom-enabled and exact full-resolution statistics benchmarks with a
  configurable `full_stats_max_bases` safety bound.
- Added optional in-memory merge profiling using four non-overlapping pieces of
  the largest selected benchmark region.
- Added optional native BigWig writing benchmarks for the largest selected
  region with zoom generation disabled and enabled.
- Added explicit opt-in full-memory loading benchmarks; full-memory reads are
  excluded from the default suite to reduce out-of-memory risk.
- Added structured per-iteration metrics for elapsed/CPU time, approximate R
  heap maxima, returned object size and rows, genomic throughput, file
  throughput, status, and captured diagnostic messages.
- Added the `bwToolsBenchmark` return contract with system, file, config,
  regions, results, and summary components.
- Added regression coverage for the benchmark contract, bounded exact-stat
  skips, explicit regions, merge/write profiling, full-memory profiling, and
  non-BigWig rejection.

# bwTools 0.6.3

## BigWig dispatch fix

- Fixed a regression introduced during the 0.6.1 API standardization where
  `vapply()` attached file-path names to automatically detected format labels.
- Fixed `read_bwg()` dispatch because `identical(fmt, "bigwig")` returns
  `FALSE` for a named scalar even when its value is `"bigwig"`.
- Prevented valid BigWig inputs from falling through to the WIG text reader,
  which caused `trimws()` to report `input string 1 is invalid UTF-8`.
- Automatic format detection now explicitly returns an unnamed format vector,
  and per-file dispatch uses scalar extraction without retaining attributes.
- Added regression coverage for automatic BigWig dispatch with explicit sample
  IDs.

# bwTools 0.6.2

## Binary-safe format detection

- Hardened format detection so binary BigWig content is checked before the
  text-format probe, reducing invalid UTF-8 failures during format detection.
- Changed `detect_bwg_format()` to inspect the BigWig four-byte magic signature
  before any filename or text decoding.
- Added a binary guard so unrecognized binary inputs fail with a format error
  instead of being passed to `readLines()`.
- Made display-only path handling UTF-8 safe without changing the path used for
  file access.
- Avoided deriving filename-based sample IDs when `sample_ids` is supplied
  explicitly.
- Added regression coverage for BigWig files without a standard filename
  extension and for unrecognized binary input.

# bwTools 0.6.1

## Standardized public API

- Standardized object-taking public functions on the primary argument name `x`.
- Standardized sample identifiers and sample selection on `sample_ids`.
- Renamed the public constructor from `bw_track()` to `bwg_track()`.
- Renamed format detection from `bw_detect_format()` to
  `detect_bwg_format()`.
- Added `validate_bwg()` for strict downstream-package contract validation.
- Added `samples_bwg()` as the public sample-metadata accessor.
- Added `metadata_bwg()` as the public schema, operation, and provenance
  metadata accessor.
- Restricted `seqinfo_bwg()`, `zoominfo_bwg()`, and `stats_bwg()` to validated
  BwgTrack input so file I/O remains centralized in `read_bwg()`.
- Changed `summary_bwg()` to return a one-row data.table with stable fields.

## BwgTrack schema v2

- Standardized `samples` core columns as `sample_id`, `file`, `format`,
  `strand`, and `has_strand`.
- Standardized memory signal columns as `sample_id`, `chrom`, `start`, `end`,
  `value`, and `strand`.
- Standardized chromosome metadata as `sample_id`, `chrom`, and `length`.
- Standardized required metadata fields as `coordinate`, `schema_version`,
  `backend`, and `mode`.
- Added strict validation for sample IDs, interval coordinates, finite signal
  values, strand values, chromosome lengths, lazy BigWig backing files, and
  metadata consistency.
- Fixed schema-v2 core column types so downstream packages receive deterministic
  character, logical, integer, and numeric fields after `bwg_track()` normalization.
- `bwg_track()` now normalizes missing standard sample columns and expands
  sample-independent seqinfo across samples.

## Standardized loading behavior

- Added `mode = "auto"` as the default for `read_bwg()`.
- All-BigWig input defaults to lazy indexed access.
- WIG and bedGraph input defaults to memory mode.
- Explicit `mode = "memory"` and `mode = "lazy"` remain available.

## Integration contract

- Added `inst/API_CONTRACT.md` describing the supported downstream-package
  integration boundary.
- Reworked `README.qmd` around the standardized public workflow and return
  contracts.
- Added API contract regression tests covering schema normalization, strict
  validation, stable return columns, sample and metadata accessors, auto loading
  mode, malformed merge inputs, and object-only analysis interfaces.

## Breaking changes

- Removed public `bw_track()`; use `bwg_track()`.
- Removed public `bw_detect_format()`; use `detect_bwg_format()`.
- Replaced `sample_names` with `sample_ids` in `read_bwg()`.
- Replaced sample selectors named `samples` with `sample_ids`.
- `seqinfo_bwg()`, `zoominfo_bwg()`, and `stats_bwg()` no longer accept file
  paths directly. Use `read_bwg()` first.
- `summary_bwg()` now returns a one-row data.table instead of a named list.

# bwTools 0.5.1

- Unified direct merging and arithmetic aggregation under the single public `merge_bwg()` interface.
- Added conservative `method = "auto"` as the default so users do not need to classify sample identity or regional overlap before merging.
- Auto mode preserves different samples, appends disjoint same-sample regions, automatically consolidates exact duplicates and equal-valued overlaps, and rejects only conflicting same-sample overlaps.
- Moved `mean` and `sum` aggregation into `merge_bwg()` and removed `collapse_bwg()` from the public API.
- Changed the default mean rule to `missing = "ignore"`; `missing = "zero"` remains available when uncovered members should explicitly contribute zero.
- Preserved grouped replicate/treatment aggregation, strand-aware arithmetic, chromosome-length validation, sparse zero handling, provenance metadata, and writer separation.
- Updated README validation and regression tests around the simplified automatic merge workflow.

# bwTools 0.5.0

- Added `merge_bwg()` for direct union of same- or different-sample signal tracks without changing values.
- Added `collapse_bwg()` for `sum` and `mean` aggregation across partially or fully overlapping intervals.
- Added atomic interval segmentation so arithmetic aggregation does not require identical input boundaries.
- Added explicit mean semantics through `missing = "zero"` and `missing = "ignore"`.
- Added sample grouping through `groups` for replicate- or treatment-level aggregation.
- Added strand-aware aggregation, chromosome-length conflict checks, explicit duplicate handling, compact adjacent-segment merging, and provenance metadata.
- Kept file persistence outside merge/collapse; merged objects are written only through `write_bwg()`.
- Added regression coverage for the four core merge scenarios, partial overlaps, grouping, duplicates, strand handling, and sequence-metadata conflicts.

# bwTools 0.4.1

- Consolidated genomic extraction into `retrieve_bwg()` and removed the overlapping `subset_bwg()` public API.
- Added multi-region retrieval through the `regions` argument while preserving the existing single-region `chrom`/`start`/`end` interface.
- Added `result = "track"` for write-ready in-memory `BwgTrack` results while keeping `result = "data"` as the backward-compatible default.
- Removed all file-output, format, compression, and zoom controls from genomic extraction; persistence remains exclusively in `write_bwg()`.
- Kept indexed lazy BigWig retrieval, region normalization, boundary clipping, and original genomic coordinates in the unified retrieval path.
- Added regression tests for single-region compatibility, multi-region unions, lazy indexed retrieval, write-ready track results, and query validation.
- Reworked `README.qmd` and the development roadmap around the separated read/retrieve/stats/write responsibilities.

# bwTools 0.4.0

- Added indexed multi-region genomic extraction from memory and lazy BigWig tracks.
- Added region normalization that merges overlapping or directly adjacent intervals before extraction.
- Added boundary clipping while preserving original genomic coordinates and chromosome lengths.
- This functionality was consolidated into `retrieve_bwg()` in 0.4.1 before API stabilization.

# bwTools 0.3.4

- Fixed WIG chromosome grouping where data.table non-standard evaluation caused every chromosome block to receive all signal intervals.
- Fixed sample-specific chromosome-size selection and in-memory chromosome statistics affected by the same column-name scoping risk.
- Made the documentation fence regression test compatible with older supported testthat releases.
- Added regression tests for WIG chromosome isolation, sample-specific seqinfo selection, and chromosome-specific in-memory statistics.

# bwTools 0.3.3

- Reworked README signal round-trip validation to compare canonicalized signal tables rather than raw data-frame representations.
- Added diagnostic output for WIG round-trip mismatches so coordinate, value, ordering, and type differences can be identified directly.
- Added a bundled-example WIG round-trip regression test covering all 296 intervals and heterogeneous interval spans.

# bwTools 0.3.2

- Standardized all language-tagged Markdown and QMD code fences to braced Quarto/knitr form such as `{r}`, `{bash}`, and `{text}`.
- Added a regression test that rejects unbraced language code fences in package Markdown and QMD files.

# bwTools 0.3.1

- Added public `zoominfo_bwg()` for stable inspection of BigWig zoom levels without accessing internal metadata structures.
- Updated zoom regression tests to use the public zoom metadata interface where appropriate.
- Replaced `README.md` as the documentation source with `README.qmd`; `README.md` is retained as the rendered GitHub-facing output.
- Reworked the README into an executable manual validation workflow covering reader, writer, zoom, statistics, WIG, and bedGraph round trips.
- Removed user-facing examples that accessed `bwTools:::bw_metadata()`.

# bwTools 0.3.0

- Added native-R BigWig zoom-level construction and writing.
- Added BigWig zoom-header parsing and zoom R-tree querying.
- Added `stats_bwg()` for coverage-aware `mean`, `stdev`, `max`, `min`, `coverage`, and `sum` statistics across genomic bins.
- Added automatic use of an appropriate zoom resolution for coarse lazy-BigWig statistics, with `use_zoom = FALSE` for exact full-resolution calculations.
- Added `zoom` and `max_zoom_levels` controls to `write_bwg()`.
- Preserved `zoom = FALSE` as a supported writer mode.
- Added total-summary parsing to native BigWig metadata.
- Added regression tests for zoom headers, zoom records, zoom-disabled output, exact statistics, missing bins, and zoom/full summary agreement.
- Updated all Markdown R code fences to Quarto/knitr form using `{r}` chunk headers.

# bwTools 0.2.0

- Added `write_bwg()` with the GeneTrackR-compatible `outdir`, `format`, `samples`, `chrom_sizes`, `overwrite`, and `compress` interface.
- Added a pure R BigWig v4 writer with fixed-header and total-summary output, chromosome B+ tree construction, zlib-compressed bedGraph-like data blocks, and full-data R-tree construction.
- Removed the need for UCSC `bedGraphToBigWig`; BigWig output does not invoke external executables or compiled BigWig libraries.
- Added direct WIG and bedGraph export from in-memory `BwgTrack` objects.
- Added chromosome-size, coordinate, chromosome-boundary, finite-value, and overlapping-interval validation before BigWig writing.
- Allow `chrom_sizes` to be omitted when complete chromosome lengths are available from `BwgTrack$seqinfo`.
- Preserve lazy same-format file copying without decoding the source track.
- Added native binary packing tests and BigWig writer round-trip regression tests against the bundled example dataset.
- BigWig zoom levels are intentionally not written in 0.2.x and remain planned for 0.3.x.

# bwTools 0.1.3

- Added a bundled three-chromosome BigWig example dataset under `inst/extdata/` for documentation and regression testing.
- Added the matching bedGraph reference and chromosome-size file for coordinate/value validation and future writer round-trip tests.
- Expanded `README.md` with package design, object contract, example-data description, lazy and memory BigWig workflows, bedGraph coordinate conversion, and the development roadmap.

# bwTools 0.1.2

- Fixed package-internal `data.table` awareness so `[.data.table` no longer falls back to data.frame semantics during `devtools::test()` or package checks.
- Replaced package-internal `:=` update expressions with `data.table::set()` for explicit namespace-safe by-reference assignment.
- Copy user-supplied `data.table` inputs in `bwg_track()` before normalization to avoid modifying caller-owned objects by reference.
- Added a regression test for the `.datatable.aware` package marker.

# bwTools 0.1.1

- Converted `NAMESPACE` and manual pages to a roxygen2-managed documentation workflow.
- Added package-level roxygen documentation and explicit documentation for the `print.bwToolsTrack()` S3 method.
- Added a namespace regression test for the public API and S3 registration.

# bwTools 0.1.0

## New features

- Add a standalone R package for genomic signal-track I/O.
- Add a pure R local bigWig reader using header parsing, chromosome B+ tree parsing,
  full-data R-tree interval lookup, zlib block decompression, and decoding of
  bedGraph, variableStep, and fixedStep bigWig blocks.
- Add bedGraph and WIG readers with conversion to a canonical 1-based closed
  in-memory coordinate contract.
- Add a GeneTrackR-compatible `BwgTrack` object contract with memory and lazy
  bigWig modes.
- Add chromosome metadata lookup and interval retrieval APIs.
