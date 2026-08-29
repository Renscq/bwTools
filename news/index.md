# Changelog

## bwTools 0.8.9

### R CMD check portability fixes

- Reworked package tests to validate installed package assets through
  [`system.file()`](https://rdrr.io/r/base/system.file.html) and
  [`packageDescription()`](https://rdrr.io/r/utils/packageDescription.html)
  instead of assuming the test working directory is the source-package
  root.
- Moved source-tree-only README, roxygen, pkgdown, and generated-Rd
  contracts into `tools/check-repository.003.R`, which runs before
  GitHub CI jobs.
- Kept installed compatibility checks pointed at `compatibility/` and
  `COMPATIBILITY.md`, matching how `inst/` assets are installed by R.
- Qualified
  [`utils::object.size()`](https://rdrr.io/r/utils/object.size.html) in
  benchmark helpers and registered the data.table
  non-standard-evaluation symbols reported by `R CMD check`.
- Preserved the 15-function public API, BwgTrack schema v2, binary I/O
  algorithms, and all numerical results.

## bwTools 0.8.8

### pkgdown reference-index fix

- Marked the `print.bwToolsBenchmark` and `print.bwToolsTrack` S3 help
  topics with the `internal` keyword so pkgdown no longer requires them
  in the user-facing reference index.
- Preserved S3 print registration, behavior, and Rd help while keeping
  the 15-function public API index unchanged.
- Added repository-level regression checks for internal print topics and
  advanced the repository contract checker to `check-repository.002.R`.
- Updated GitHub workflows to use the new repository checker.
- No package function body, public signature, BwgTrack schema, or binary
  I/O implementation changed.

## bwTools 0.8.7

### Release-quality checks and GitHub pkgdown deployment

- Added GitHub Actions workflows for cross-platform `R CMD check` and
  automatic pkgdown site deployment from `main`/`master` to the
  `gh-pages` branch.
- Based the workflows on the current `r-lib/actions` v2 examples,
  including current R setup, dependency installation, and pkgdown site
  build commands.
- Added repository and website metadata to DESCRIPTION and
  `_pkgdown.yml` for `Renscq/bwTools` and the GitHub Pages site.
- Added a repository-level `tools/check-repository.001.R` contract
  checker for GitHub workflows, pkgdown configuration, README
  synchronization, code fences, and release-only `.Rbuildignore` rules.
- Moved repository-only documentation assertions out of package testthat
  logic: installed tests now validate `README.md` and shipped vignettes
  only, avoiding failures caused by intentionally Rbuildignored
  pkgdown/GitHub source files.
- Excluded `DEVELOPMENT.md`, `README.qmd`, and `tools/` from
  source-package tarballs while retaining them in the Git repository.
- Updated README with GitHub installation, CI/pkgdown badges, published
  article links, and the package website URL.
- Kept all 15 public signatures, BwgTrack schema v2, and production R
  function bodies unchanged.

## bwTools 0.8.6

### Documentation and user workflow consolidation

- Reduced README to a concise package entry point covering installation,
  the standard data flow, coordinates, multi-sample use, the stable
  public API, errors, diagnostics, and links to focused articles.
- Added six package vignettes for getting started, reading/writing,
  regional retrieval/statistics, merge/multi-sample workflows, file
  formats/coordinates, and performance/compatibility validation.
- Added pkgdown navigation that groups the 15-function public API by
  file I/O, BwgTrack contract, signal analysis, and diagnostics.
- Standardized roxygen details, runnable examples, and cross-references
  for all 15 exported functions without changing function bodies or
  public signatures.
- Added documentation contract tests requiring stable article assets and
  `@return`/`@examples` coverage for every exported function.
- Extended the code-fence regression to R Markdown vignettes.
- Added `knitr` and `rmarkdown` as documentation-only suggested
  dependencies; runtime dependencies remain unchanged.
- Simplified DESCRIPTION to describe current package capabilities rather
  than embedding release-history text.

## bwTools 0.8.5

### Multi-sample and mixed-format workflow hardening

- `retrieve_bwg(strand = ...)` now uses row-level strand values for
  in-memory tracks, so mixed-strand samples produced by merges remain
  retrievable by exact `+`, `-`, or `*` strand.
- Lazy BigWig retrieval continues to use file/sample-level strand
  metadata, preserving indexed access without materializing signal.
- `result = "track"` normalizes returned sample strand metadata to the
  exact requested strand for in-memory strand-filtered retrieval.
- [`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md)
  now preflights all selected in-memory samples and raises
  `bwTools_error` before writing any files if one or more selected
  samples have no signal records, preventing incomplete manifests and
  partial output.
- Added workflow regression coverage for mixed BigWig/bedGraph reads,
  BigWig-only multi-sample lazy reads, duplicate sample IDs,
  mixed-strand merge retrieval, multi-sample writer manifests, and empty
  selected samples.
- Removed the `Development status` section from README as requested.
- Kept all 15 public function signatures, BwgTrack schema v2, native
  BigWig binary I/O, statistics, zoom, and merge algorithms unchanged.

## bwTools 0.8.4

### Cross-implementation BigWig compatibility validation

- Added an installed, optional compatibility harness for
  pyBigWig/libBigWig, Bioconductor rtracklayer, and UCSC BigWig
  command-line utilities.
- The harness validates both `bwTools -> external` and
  `external -> bwTools` directions when a reference implementation is
  available.
- Both zoom-disabled and zoom-enabled bwTools BigWig outputs are
  supplied to external readers.
- Compatibility comparisons use chromosome lengths, missing-signal
  positions, and per-base signal values rather than byte-identical files
  or identical internal block segmentation.
- External implementations remain development-time references only and
  are not added to `Depends`, `Imports`, or normal package runtime code.
- Unavailable implementations are reported as `SKIP`; detected
  disagreements are reported as `FAIL` with a non-zero runner exit
  status.
- Added a machine-readable `compatibility-report.tsv` output and strict
  mode for prepared release-validation environments.
- Added regression tests that protect the compatibility assets and
  runtime dependency policy.
- No public API signatures, BwgTrack schema, or native BigWig
  reader/writer algorithms changed.

## bwTools 0.8.3

### Scalar retrieval validation consistency

- `retrieve_bwg(chrom, start, end)` now reuses the same shared genomic
  interval validator as
  [`stats_bwg()`](https://renscq.github.io/bwTools/reference/stats_bwg.md).
- Invalid scalar intervals such as `end < start` now return the same
  precise `bwTools_error` message across retrieval and statistics.
- Batch `regions=` validation remains unchanged and continues to report
  invalid rows with the table-level interval error.
- No public function signatures, BwgTrack schema, BigWig binary I/O, or
  performance paths changed.

## bwTools 0.8.2

### Boundary and edge-case hardening

- Unified public genomic query bounds across memory and lazy BigWig
  tracks when chromosome lengths are known in `seqinfo`.
- [`retrieve_bwg()`](https://renscq.github.io/bwTools/reference/retrieve_bwg.md)
  now clips query ends to known chromosome lengths, returns a stable
  empty signal table when a query starts beyond the chromosome, and
  records the effective bounded regions in `result = "track"` metadata.
- [`stats_bwg()`](https://renscq.github.io/bwTools/reference/stats_bwg.md)
  now applies the same chromosome-end clipping before bin construction
  and returns a typed empty statistics table for fully out-of-range
  queries.
- Unknown chromosomes now raise `bwTools_error` when selected samples
  provide complete known-length sequence metadata; samples with any
  unknown sequence length remain permissive because their metadata may
  be incomplete.
- Strengthened `BwgTrack` validation so in-memory signal
  sample/chromosome pairs must exist in supplied `seqinfo`, and known
  chromosome lengths cannot be exceeded by signal intervals.
- Added regression coverage for point queries,
  invalid/fractional/out-of-range coordinates, empty signal tracks,
  signed finite values, float32 writer limits, exact chromosome-name
  preservation, region-union normalization, and the R integer coordinate
  boundary.
- Kept BwgTrack schema v2, all 15 public function signatures, native
  BigWig binary reader/writer layout, and optimized 0.7.x performance
  paths unchanged.

## bwTools 0.8.1

### Tests

- Corrected the bounded full-stat benchmark regression so the 1000 bp
  case runs below a 1500 bp limit while the 2000 bp case is skipped.
- Corrected the writer size-limit regression to use a 1000 bp
  bundled-example region with known signal coverage; the previous 100 bp
  centered region was legitimately empty and therefore skipped by
  [`benchmark_bwg()`](https://renscq.github.io/bwTools/reference/benchmark_bwg.md).
- No production R code or public API behavior changed in this patch
  release.

## bwTools 0.8.0

### Public API stability and cross-format hardening

- Established 0.8.0 as the downstream-integration stability baseline
  while retaining BwgTrack schema v2 and the optimized 0.7.5
  reader/writer paths.
- Added a common `bwTools_error` condition class for errors
  intentionally generated by bwTools validation and I/O checks.
- Replaced public [`match.arg()`](https://rdrr.io/r/base/match.arg.html)
  dispatch with package-owned exact choice validation so invalid public
  choices use the same bwTools condition class and partial argument
  matching is no longer accepted.
- Added shared scalar, logical, positive-integer, and genomic-interval
  validators and applied them to public statistics, merge, writer,
  reader, retrieval, and benchmark entry points where appropriate.
- Tightened exact-statistics coordinates so fractional genomic
  coordinates are rejected rather than silently truncated by integer
  coercion.
- Added regression tests that lock public function parameter names and
  order for all 15 exported functions.
- Added cross-format regression coverage for BigWig zoom/no-zoom, WIG,
  WIG.gz, bedGraph, and bedGraph.gz round trips through the standardized
  signal schema.
- Added stable writer-manifest checks for `sample_id`, `file`, and
  `format`.
- Reorganized `README.qmd` around the normal user workflow and moved
  benchmarking into a clearly marked advanced diagnostics section.
- Expanded `inst/API_CONTRACT.md` with stable argument vocabulary,
  return schemas, error-condition behavior, coordinate rules, merge
  semantics, and cross-format compatibility guarantees.

## bwTools 0.7.5

### Hierarchical native BigWig zoom writer

- Kept the validated 0.7.2 reader/retrieval baseline and 0.7.1 streaming
  exact statistics unchanged.
- Replaced the original per-interval first-level zoom builder with a
  vectorized same-window aggregation path and a limited splitter for
  boundary-crossing intervals.
- Added hierarchical zoom aggregation so level 2 and higher are built
  exactly from the preceding zoom summaries rather than rescanning the
  complete signal for every reduction level.
- Preserved `valid_count`, minimum, maximum, sum, and squared-sum
  semantics across hierarchical aggregation.
- Updated the writer to retain only the preceding zoom level while
  emitting the next level, avoiding a full in-memory zoom pyramid during
  file writing.
- Kept the public
  [`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md)
  API, BigWig binary layout, BwgTrack schema v2, and benchmark interface
  unchanged.
- Added regression tests comparing hierarchical summaries with direct
  full-signal summaries across multiple levels and boundary-crossing
  intervals.

## bwTools 0.7.4

### Bounded and observable writer benchmarking

- Fixed benchmark warmup semantics so untimed warmups run once per
  benchmark case instead of once before every timed iteration.
- Added independent `write_iterations` and `write_warmup` controls with
  safe writer defaults of one timed write and no warmup write.
- Added `write_max_intervals` to bound writer profiling by actual signal
  density rather than genomic span alone.
- Added `write_zoom_max_intervals` as a stricter safety guard for the
  current pure-R zoom builder, which scans dense signals across multiple
  zoom levels.
- Added `verbose = TRUE` progress messages for writer materialization,
  skips, writer variants, and timed completion.
- Kept native BigWig writer output, reader behavior, BwgTrack schema v2,
  and all normal analysis APIs unchanged.
- Added regression tests for writer-specific repetition controls and
  interval-density safety guards.

## bwTools 0.7.3

### Native BigWig writer benchmark and profiling

- Kept the 0.7.2 indexed reader, block-streamed exact statistics, zoom
  logic, BwgTrack schema v2, and public analysis API unchanged.
- Expanded the optional `write` benchmark to profile every selected
  benchmark region rather than only the largest region.
- Added `write_max_bases` as an explicit safety bound for writer stress
  tests.
- Benchmarked native BigWig writing with `zoom = FALSE` and
  `zoom = TRUE` for each eligible region while excluding indexed
  retrieval from timed writer work.
- Added writer-input diagnostics for materialized signal size, interval
  count, covered bases, estimated uncompressed type-1 data payload, and
  data-block count.
- Added writer throughput metrics for input signal MB/s, estimated
  payload MB/s, intervals/s, output MB/s, and genomic Mb/s.
- Added output-size diagnostics comparing final BigWig size with the
  in-memory signal object and estimated uncompressed primary data
  payload.
- Added summary-level zoom elapsed-time and output-size overhead
  relative to the corresponding `zoom_off` writer benchmark.
- Added regression tests for multi-region writer profiling, new writer
  metrics, zoom-overhead summaries, and `write_max_bases` skipping
  behavior.

## bwTools 0.7.2

### Indexed retrieval and benchmark memory refinement

- Retained the 0.7.1 block-streamed exact-statistics implementation
  after real large-file benchmarking confirmed a substantial
  elapsed-time improvement.
- Replaced growable full-query vectors in indexed BigWig retrieval with
  per-block primitive-vector chunks followed by one-shot flattening.
- Avoided repeated R vector resizing and copying while continuing to
  avoid per-block data.table construction.
- Added `live_heap_after_mb`, `peak_heap_delta_mb`, and
  `live_heap_delta_mb` to benchmark iteration results.
- Added corresponding median live/peak heap metrics to benchmark
  summaries.
- Retained `heap_delta_mb` as a compatibility alias of
  `peak_heap_delta_mb`.
- Clarified that benchmark heap values are R-heap diagnostics rather
  than operating-system RSS measurements.
- Kept the public API, BwgTrack schema v2, streaming exact-statistics
  behavior, zoom behavior, and coordinate conventions unchanged.
- Added regression coverage for the expanded benchmark memory contract
  and retained indexed-retrieval signal equivalence tests.

## bwTools 0.7.1

### Benchmark-driven BigWig optimization

- Reworked exact lazy-BigWig statistics so `stats_bwg(use_zoom = FALSE)`
  streams overlapping binary blocks directly into bin accumulators
  instead of first materializing a complete full-resolution signal
  table.
- Replaced interval-by-interval full-stat binning with vectorized
  same-bin accumulation and limited loops only for intervals crossing
  bin boundaries.
- Split full-data block decoding into a reusable vector decoder so
  statistics can consume coordinates and values without per-block
  data.table allocation.
- Reduced indexed-retrieval allocations by filling growable typed
  vectors and constructing the returned data.table only once.
- Added a single-result retrieval fast path to avoid unnecessary
  `rbindlist()` and re-sorting for the common one-sample, one-region
  query.
- Removed an unnecessary full copy of memory-mode signal data before
  regional filtering.
- Kept the BwgTrack schema v2, public functions, argument names, zoom
  semantics, and benchmark result contract unchanged.
- Added regression coverage requiring streaming exact BigWig statistics
  to match memory-mode exact statistics for mean, stdev, max, min,
  coverage, and sum, plus indexed-retrieval equivalence against the
  bundled bedGraph signal.

## bwTools 0.7.0

### Large BigWig benchmark framework

- Added public
  [`benchmark_bwg()`](https://renscq.github.io/bwTools/reference/benchmark_bwg.md)
  for reproducible profiling of large local BigWig files before
  performance optimization.
- Added package-metadata-cold and cached lazy-read benchmarks without
  claiming to flush the operating-system file cache.
- Added indexed retrieval benchmarks across user-defined or
  deterministic genomic windows.
- Added zoom-enabled and exact full-resolution statistics benchmarks
  with a configurable `full_stats_max_bases` safety bound.
- Added optional in-memory merge profiling using four non-overlapping
  pieces of the largest selected benchmark region.
- Added optional native BigWig writing benchmarks for the largest
  selected region with zoom generation disabled and enabled.
- Added explicit opt-in full-memory loading benchmarks; full-memory
  reads are excluded from the default suite to reduce out-of-memory
  risk.
- Added structured per-iteration metrics for elapsed/CPU time,
  approximate R heap maxima, returned object size and rows, genomic
  throughput, file throughput, status, and captured diagnostic messages.
- Added the `bwToolsBenchmark` return contract with system, file,
  config, regions, results, and summary components.
- Added regression coverage for the benchmark contract, bounded
  exact-stat skips, explicit regions, merge/write profiling, full-memory
  profiling, and non-BigWig rejection.

## bwTools 0.6.3

### BigWig dispatch fix

- Fixed a regression introduced during the 0.6.1 API standardization
  where [`vapply()`](https://rdrr.io/r/base/lapply.html) attached
  file-path names to automatically detected format labels.
- Fixed
  [`read_bwg()`](https://renscq.github.io/bwTools/reference/read_bwg.md)
  dispatch because `identical(fmt, "bigwig")` returns `FALSE` for a
  named scalar even when its value is `"bigwig"`.
- Prevented valid BigWig inputs from falling through to the WIG text
  reader, which caused [`trimws()`](https://rdrr.io/r/base/trimws.html)
  to report `input string 1 is invalid UTF-8`.
- Automatic format detection now explicitly returns an unnamed format
  vector, and per-file dispatch uses scalar extraction without retaining
  attributes.
- Added regression coverage for automatic BigWig dispatch with explicit
  sample IDs.

## bwTools 0.6.2

### Binary-safe format detection

- Hardened format detection so binary BigWig content is checked before
  the text-format probe, reducing invalid UTF-8 failures during format
  detection.
- Changed
  [`detect_bwg_format()`](https://renscq.github.io/bwTools/reference/detect_bwg_format.md)
  to inspect the BigWig four-byte magic signature before any filename or
  text decoding.
- Added a binary guard so unrecognized binary inputs fail with a format
  error instead of being passed to
  [`readLines()`](https://rdrr.io/r/base/readLines.html).
- Made display-only path handling UTF-8 safe without changing the path
  used for file access.
- Avoided deriving filename-based sample IDs when `sample_ids` is
  supplied explicitly.
- Added regression coverage for BigWig files without a standard filename
  extension and for unrecognized binary input.

## bwTools 0.6.1

### Standardized public API

- Standardized object-taking public functions on the primary argument
  name `x`.
- Standardized sample identifiers and sample selection on `sample_ids`.
- Renamed the public constructor from `bw_track()` to
  [`bwg_track()`](https://renscq.github.io/bwTools/reference/bwg_track.md).
- Renamed format detection from `bw_detect_format()` to
  [`detect_bwg_format()`](https://renscq.github.io/bwTools/reference/detect_bwg_format.md).
- Added
  [`validate_bwg()`](https://renscq.github.io/bwTools/reference/validate_bwg.md)
  for strict downstream-package contract validation.
- Added
  [`samples_bwg()`](https://renscq.github.io/bwTools/reference/samples_bwg.md)
  as the public sample-metadata accessor.
- Added
  [`metadata_bwg()`](https://renscq.github.io/bwTools/reference/metadata_bwg.md)
  as the public schema, operation, and provenance metadata accessor.
- Restricted
  [`seqinfo_bwg()`](https://renscq.github.io/bwTools/reference/seqinfo_bwg.md),
  [`zoominfo_bwg()`](https://renscq.github.io/bwTools/reference/zoominfo_bwg.md),
  and
  [`stats_bwg()`](https://renscq.github.io/bwTools/reference/stats_bwg.md)
  to validated BwgTrack input so file I/O remains centralized in
  [`read_bwg()`](https://renscq.github.io/bwTools/reference/read_bwg.md).
- Changed
  [`summary_bwg()`](https://renscq.github.io/bwTools/reference/summary_bwg.md)
  to return a one-row data.table with stable fields.

### BwgTrack schema v2

- Standardized `samples` core columns as `sample_id`, `file`, `format`,
  `strand`, and `has_strand`.
- Standardized memory signal columns as `sample_id`, `chrom`, `start`,
  `end`, `value`, and `strand`.
- Standardized chromosome metadata as `sample_id`, `chrom`, and
  `length`.
- Standardized required metadata fields as `coordinate`,
  `schema_version`, `backend`, and `mode`.
- Added strict validation for sample IDs, interval coordinates, finite
  signal values, strand values, chromosome lengths, lazy BigWig backing
  files, and metadata consistency.
- Fixed schema-v2 core column types so downstream packages receive
  deterministic character, logical, integer, and numeric fields after
  [`bwg_track()`](https://renscq.github.io/bwTools/reference/bwg_track.md)
  normalization.
- [`bwg_track()`](https://renscq.github.io/bwTools/reference/bwg_track.md)
  now normalizes missing standard sample columns and expands
  sample-independent seqinfo across samples.

### Standardized loading behavior

- Added `mode = "auto"` as the default for
  [`read_bwg()`](https://renscq.github.io/bwTools/reference/read_bwg.md).
- All-BigWig input defaults to lazy indexed access.
- WIG and bedGraph input defaults to memory mode.
- Explicit `mode = "memory"` and `mode = "lazy"` remain available.

### Integration contract

- Added `inst/API_CONTRACT.md` describing the supported
  downstream-package integration boundary.
- Reworked `README.qmd` around the standardized public workflow and
  return contracts.
- Added API contract regression tests covering schema normalization,
  strict validation, stable return columns, sample and metadata
  accessors, auto loading mode, malformed merge inputs, and object-only
  analysis interfaces.

### Breaking changes

- Removed public `bw_track()`; use
  [`bwg_track()`](https://renscq.github.io/bwTools/reference/bwg_track.md).
- Removed public `bw_detect_format()`; use
  [`detect_bwg_format()`](https://renscq.github.io/bwTools/reference/detect_bwg_format.md).
- Replaced `sample_names` with `sample_ids` in
  [`read_bwg()`](https://renscq.github.io/bwTools/reference/read_bwg.md).
- Replaced sample selectors named `samples` with `sample_ids`.
- [`seqinfo_bwg()`](https://renscq.github.io/bwTools/reference/seqinfo_bwg.md),
  [`zoominfo_bwg()`](https://renscq.github.io/bwTools/reference/zoominfo_bwg.md),
  and
  [`stats_bwg()`](https://renscq.github.io/bwTools/reference/stats_bwg.md)
  no longer accept file paths directly. Use
  [`read_bwg()`](https://renscq.github.io/bwTools/reference/read_bwg.md)
  first.
- [`summary_bwg()`](https://renscq.github.io/bwTools/reference/summary_bwg.md)
  now returns a one-row data.table instead of a named list.

## bwTools 0.5.1

- Unified direct merging and arithmetic aggregation under the single
  public
  [`merge_bwg()`](https://renscq.github.io/bwTools/reference/merge_bwg.md)
  interface.
- Added conservative `method = "auto"` as the default so users do not
  need to classify sample identity or regional overlap before merging.
- Auto mode preserves different samples, appends disjoint same-sample
  regions, automatically consolidates exact duplicates and equal-valued
  overlaps, and rejects only conflicting same-sample overlaps.
- Moved `mean` and `sum` aggregation into
  [`merge_bwg()`](https://renscq.github.io/bwTools/reference/merge_bwg.md)
  and removed `collapse_bwg()` from the public API.
- Changed the default mean rule to `missing = "ignore"`;
  `missing = "zero"` remains available when uncovered members should
  explicitly contribute zero.
- Preserved grouped replicate/treatment aggregation, strand-aware
  arithmetic, chromosome-length validation, sparse zero handling,
  provenance metadata, and writer separation.
- Updated README validation and regression tests around the simplified
  automatic merge workflow.

## bwTools 0.5.0

- Added
  [`merge_bwg()`](https://renscq.github.io/bwTools/reference/merge_bwg.md)
  for direct union of same- or different-sample signal tracks without
  changing values.
- Added `collapse_bwg()` for `sum` and `mean` aggregation across
  partially or fully overlapping intervals.
- Added atomic interval segmentation so arithmetic aggregation does not
  require identical input boundaries.
- Added explicit mean semantics through `missing = "zero"` and
  `missing = "ignore"`.
- Added sample grouping through `groups` for replicate- or
  treatment-level aggregation.
- Added strand-aware aggregation, chromosome-length conflict checks,
  explicit duplicate handling, compact adjacent-segment merging, and
  provenance metadata.
- Kept file persistence outside merge/collapse; merged objects are
  written only through
  [`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md).
- Added regression coverage for the four core merge scenarios, partial
  overlaps, grouping, duplicates, strand handling, and sequence-metadata
  conflicts.

## bwTools 0.4.1

- Consolidated genomic extraction into
  [`retrieve_bwg()`](https://renscq.github.io/bwTools/reference/retrieve_bwg.md)
  and removed the overlapping `subset_bwg()` public API.
- Added multi-region retrieval through the `regions` argument while
  preserving the existing single-region `chrom`/`start`/`end` interface.
- Added `result = "track"` for write-ready in-memory `BwgTrack` results
  while keeping `result = "data"` as the backward-compatible default.
- Removed all file-output, format, compression, and zoom controls from
  genomic extraction; persistence remains exclusively in
  [`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md).
- Kept indexed lazy BigWig retrieval, region normalization, boundary
  clipping, and original genomic coordinates in the unified retrieval
  path.
- Added regression tests for single-region compatibility, multi-region
  unions, lazy indexed retrieval, write-ready track results, and query
  validation.
- Reworked `README.qmd` and the development roadmap around the separated
  read/retrieve/stats/write responsibilities.

## bwTools 0.4.0

- Added indexed multi-region genomic extraction from memory and lazy
  BigWig tracks.
- Added region normalization that merges overlapping or directly
  adjacent intervals before extraction.
- Added boundary clipping while preserving original genomic coordinates
  and chromosome lengths.
- This functionality was consolidated into
  [`retrieve_bwg()`](https://renscq.github.io/bwTools/reference/retrieve_bwg.md)
  in 0.4.1 before API stabilization.

## bwTools 0.3.4

- Fixed WIG chromosome grouping where data.table non-standard evaluation
  caused every chromosome block to receive all signal intervals.
- Fixed sample-specific chromosome-size selection and in-memory
  chromosome statistics affected by the same column-name scoping risk.
- Made the documentation fence regression test compatible with older
  supported testthat releases.
- Added regression tests for WIG chromosome isolation, sample-specific
  seqinfo selection, and chromosome-specific in-memory statistics.

## bwTools 0.3.3

- Reworked README signal round-trip validation to compare canonicalized
  signal tables rather than raw data-frame representations.
- Added diagnostic output for WIG round-trip mismatches so coordinate,
  value, ordering, and type differences can be identified directly.
- Added a bundled-example WIG round-trip regression test covering all
  296 intervals and heterogeneous interval spans.

## bwTools 0.3.2

- Standardized all language-tagged Markdown and QMD code fences to
  braced Quarto/knitr form such as `{r}`, `{bash}`, and
  [text](https://r-text.org/).
- Added a regression test that rejects unbraced language code fences in
  package Markdown and QMD files.

## bwTools 0.3.1

- Added public
  [`zoominfo_bwg()`](https://renscq.github.io/bwTools/reference/zoominfo_bwg.md)
  for stable inspection of BigWig zoom levels without accessing internal
  metadata structures.
- Updated zoom regression tests to use the public zoom metadata
  interface where appropriate.
- Replaced `README.md` as the documentation source with `README.qmd`;
  `README.md` is retained as the rendered GitHub-facing output.
- Reworked the README into an executable manual validation workflow
  covering reader, writer, zoom, statistics, WIG, and bedGraph round
  trips.
- Removed user-facing examples that accessed `bwTools:::bw_metadata()`.

## bwTools 0.3.0

- Added native-R BigWig zoom-level construction and writing.
- Added BigWig zoom-header parsing and zoom R-tree querying.
- Added
  [`stats_bwg()`](https://renscq.github.io/bwTools/reference/stats_bwg.md)
  for coverage-aware `mean`, `stdev`, `max`, `min`, `coverage`, and
  `sum` statistics across genomic bins.
- Added automatic use of an appropriate zoom resolution for coarse
  lazy-BigWig statistics, with `use_zoom = FALSE` for exact
  full-resolution calculations.
- Added `zoom` and `max_zoom_levels` controls to
  [`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md).
- Preserved `zoom = FALSE` as a supported writer mode.
- Added total-summary parsing to native BigWig metadata.
- Added regression tests for zoom headers, zoom records, zoom-disabled
  output, exact statistics, missing bins, and zoom/full summary
  agreement.
- Updated all Markdown R code fences to Quarto/knitr form using `{r}`
  chunk headers.

## bwTools 0.2.0

- Added
  [`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md)
  with the GeneTrackR-compatible `outdir`, `format`, `samples`,
  `chrom_sizes`, `overwrite`, and `compress` interface.
- Added a pure R BigWig v4 writer with fixed-header and total-summary
  output, chromosome B+ tree construction, zlib-compressed bedGraph-like
  data blocks, and full-data R-tree construction.
- Removed the need for UCSC `bedGraphToBigWig`; BigWig output does not
  invoke external executables or compiled BigWig libraries.
- Added direct WIG and bedGraph export from in-memory `BwgTrack`
  objects.
- Added chromosome-size, coordinate, chromosome-boundary, finite-value,
  and overlapping-interval validation before BigWig writing.
- Allow `chrom_sizes` to be omitted when complete chromosome lengths are
  available from `BwgTrack$seqinfo`.
- Preserve lazy same-format file copying without decoding the source
  track.
- Added native binary packing tests and BigWig writer round-trip
  regression tests against the bundled example dataset.
- BigWig zoom levels are intentionally not written in 0.2.x and remain
  planned for 0.3.x.

## bwTools 0.1.3

- Added a bundled three-chromosome BigWig example dataset under
  `inst/extdata/` for documentation and regression testing.
- Added the matching bedGraph reference and chromosome-size file for
  coordinate/value validation and future writer round-trip tests.
- Expanded `README.md` with package design, object contract,
  example-data description, lazy and memory BigWig workflows, bedGraph
  coordinate conversion, and the development roadmap.

## bwTools 0.1.2

- Fixed package-internal `data.table` awareness so `[.data.table` no
  longer falls back to data.frame semantics during `devtools::test()` or
  package checks.
- Replaced package-internal `:=` update expressions with
  [`data.table::set()`](https://rdrr.io/pkg/data.table/man/assign.html)
  for explicit namespace-safe by-reference assignment.
- Copy user-supplied `data.table` inputs in
  [`bwg_track()`](https://renscq.github.io/bwTools/reference/bwg_track.md)
  before normalization to avoid modifying caller-owned objects by
  reference.
- Added a regression test for the `.datatable.aware` package marker.

## bwTools 0.1.1

- Converted `NAMESPACE` and manual pages to a roxygen2-managed
  documentation workflow.
- Added package-level roxygen documentation and explicit documentation
  for the
  [`print.bwToolsTrack()`](https://renscq.github.io/bwTools/reference/print.bwToolsTrack.md)
  S3 method.
- Added a namespace regression test for the public API and S3
  registration.

## bwTools 0.1.0

### New features

- Add a standalone R package for genomic signal-track I/O.
- Add a pure R local bigWig reader using header parsing, chromosome B+
  tree parsing, full-data R-tree interval lookup, zlib block
  decompression, and decoding of bedGraph, variableStep, and fixedStep
  bigWig blocks.
- Add bedGraph and WIG readers with conversion to a canonical 1-based
  closed in-memory coordinate contract.
- Add a GeneTrackR-compatible `BwgTrack` object contract with memory and
  lazy bigWig modes.
- Add chromosome metadata lookup and interval retrieval APIs.
