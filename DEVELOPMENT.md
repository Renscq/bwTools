# bwTools development status

Current development release: **0.8.8**

0.8.8 fixes pkgdown reference indexing after the 0.8.7 GitHub deployment
consolidation. GitHub Actions now run cross-platform R CMD check and
build/deploy the pkgdown site, while repository-only checks are
separated from package tests. The 15-function public API, BwgTrack
schema v2, and all production R function bodies remain unchanged.

------------------------------------------------------------------------

# bwTools development roadmap

## 0.8.8 pkgdown reference-index fix

- Marked `print.bwToolsBenchmark` and `print.bwToolsTrack` as internal
  documentation topics so pkgdown does not require them in the
  user-facing reference index.
- Preserved both S3 print registrations and their Rd help pages; the
  change is documentation visibility only and does not change print
  behavior.
- Added repository-contract coverage for internal S3 print topics and
  updated the GitHub workflows to run `tools/check-repository.002.R`.
- Kept all 15 public function signatures, BwgTrack schema v2, and
  production function bodies unchanged apart from roxygen metadata.

## 0.8.7 release-quality checks and GitHub pkgdown deployment

Status: **implemented; GitHub Actions and Pages deployment require
repository push.**

Scope:

- Add cross-platform `R CMD check` CI for release, devel, and oldrel R.
- Build and deploy pkgdown automatically on GitHub pushes and releases.
- Keep GitHub/pkgdown repository files outside the source-package
  tarball.
- Separate repository-only contract checks from installed package
  testthat tests.
- Add repository, issue, and website metadata without changing package
  APIs.
- Keep README source synchronized and retain the no-Development-status
  rule.

Acceptance criteria:

- `tools/check-repository.001.R` validates repository configuration
  before CI work.
- Package testthat logic does not require Rbuildignored `_pkgdown.yml`
  or README.qmd.
- GitHub pkgdown workflow builds `docs/` and deploys it to `gh-pages`.
- R-CMD-check workflow covers macOS release, Windows release, Ubuntu
  devel, Ubuntu release, and Ubuntu oldrel-1.
- `_pkgdown.yml` resolves the GitHub source, issue tracker, and
  published site.
- No public signatures, schema fields, or production R function bodies
  change.

### GitHub pkgdown deployment

After pushing the repository to `https://github.com/Renscq/bwTools`,
confirm that the `pkgdown.yaml` workflow succeeds. On first deployment,
open GitHub **Settings \> Pages** and select **Deploy from a branch**,
branch `gh-pages`, folder `/ (root)`. The published site is configured
as `https://renscq.github.io/bwTools/`.

## 0.8.x - API and boundary hardening

Status: **0.8.5 implemented; local R and external compatibility
verification required.**

Scope:

- Keep the 0.8.0 public API and BwgTrack schema v2 stable.
- Apply known chromosome bounds consistently to memory and lazy queries.
- Reject unknown chromosomes only when selected samples provide complete
  known-length sequence metadata.
- Treat any sample with unknown sequence lengths as incomplete rather
  than as a full genome dictionary.
- Keep empty retrieval and statistics results structurally stable.
- Validate `data` against supplied sample/chromosome `seqinfo` and known
  lengths.
- Preserve chromosome identifiers exactly; never add or remove `chr`
  prefixes.
- Keep the public in-memory coordinate range within R’s signed integer
  range.

Acceptance criteria:

- Memory and lazy BigWig statistics use the same effective interval at
  known chromosome ends.
- Queries beginning beyond a known chromosome return typed empty results
  rather than mode-specific structures.
- Invalid coordinates, unknown complete-reference chromosomes,
  non-finite signal values, and data/seqinfo mismatches raise
  `bwTools_error`.
- Zero and negative finite signal values round-trip through native
  BigWig.
- Chromosome names are preserved byte-for-byte at the R character level
  across BigWig, WIG, and bedGraph round trips.
- Existing public function signatures and schema-v2 core field types do
  not change.
- External BigWig compatibility is reproducibly validated in both
  directions when pyBigWig, rtracklayer, or UCSC utilities are
  available.
- External reference implementations remain outside package runtime
  dependencies and normal test requirements.

## 0.8.6 documentation and user workflow consolidation

Status: **implemented; local vignette build and R CMD check required.**

Scope:

- Keep README concise and focused on first-use workflows.
- Move detailed guidance into focused package vignettes.
- Ensure all 15 exported functions have explicit return and example
  docs.
- Add pkgdown reference and article navigation without changing runtime
  APIs.
- Keep documentation-only dependencies in `Suggests`.
- Preserve all binary I/O, statistics, retrieval, merge, and schema
  behavior.

Acceptance criteria:

- README contains no development-status section or release-history
  narrative.
- Six focused vignettes cover the normal user workflow and advanced
  diagnostics.
- Every public export has `@return` and `@examples` roxygen
  documentation.
- R Markdown code fences follow the same braced language-tag convention.
- `DESCRIPTION` declares `VignetteBuilder: knitr` and documentation-only
  suggested dependencies.
- Public function signatures and BwgTrack schema v2 remain unchanged.

## 0.8.5 multi-sample and mixed-format workflow hardening

Status: **implemented; local R regression testing required.**

Scope:

- Preserve the 15-function public API, BwgTrack schema v2, and native
  BigWig binary algorithms.
- Verify mixed BigWig/bedGraph ingestion keeps source format, sample
  identity, strand metadata, and automatic memory-mode behavior.
- Verify multi-BigWig ingestion remains lazy and sample-selectable.
- Make exact strand retrieval operate on row-level signal for
  memory-mode mixed-strand samples while retaining sample-level
  filtering for lazy BigWig.
- Require one writer-manifest row for every selected in-memory sample
  with signal and fail before file creation when any selected sample is
  empty.
- Remove the README `Development status` section.

Acceptance criteria:

- Mixed BigWig/text reads produce a valid memory-mode BwgTrack with all
  samples.
- BigWig-only multi-sample reads remain lazy.
- Exact `+`, `-`, and `*` retrieval works after a same-sample
  mixed-strand merge.
- Multi-sample writing produces exactly one manifest row per selected
  sample.
- Empty selected in-memory samples cause `bwTools_error` before partial
  files are created.
- Public signatures, schema v2, native BigWig reader/writer, statistics,
  zoom, and merge algorithms do not change.
- README contains no `Development status` heading.

## 0.8.4 external compatibility validation

Status: **harness implemented; reference implementations must be run
locally.**

The installed runner is:

`{text} inst/compatibility/run-compatibility.001.R`

Supported reference implementations:

`{text} pyBigWig / libBigWig rtracklayer UCSC bigWigInfo + bigWigToBedGraph + bedGraphToBigWig`

The runner writes a machine-readable `compatibility-report.tsv`. Missing
tools are `SKIP` by default; detected incompatibilities are `FAIL`.
`--strict` turns missing reference implementations into a non-zero
release-validation result.

## Design contract

- Package code for BigWig parsing and writing is implemented in R only.
- No bundled C/C++ source, Rcpp, libBigWig linkage, or UCSC executable
  is used at runtime.
- `data.table` is retained as the only runtime package dependency
  because GeneTrackR already uses it and its signal APIs consume
  `data.table` objects.
- Public in-memory coordinates are 1-based closed.
- BigWig and bedGraph file coordinates are converted from 0-based
  half-open at the I/O boundary.
- Large BigWig files should use lazy indexed access rather than
  full-memory loading.
- Returned objects inherit from `BwgTrack` for GeneTrackR compatibility.

## 0.1.x - Reader foundation

Scope:

- Package skeleton and stable object contract.
- bedGraph reader.
- fixedStep/variableStep WIG reader.
- Native R BigWig header parsing.
- Chromosome B+ tree parsing.
- Full-data R-tree interval lookup.
- zlib data-block decompression.
- bedGraph, variableStep, and fixedStep BigWig data-block decoding.
- Memory and lazy BigWig read modes.
- Region retrieval and sequence metadata.
- Bundled BigWig/bedGraph/chrom.sizes example dataset for integration
  and writer round-trip tests.

Acceptance criteria:

- Coordinate conversion is deterministic and covered by tests.
- BigWig fixture metadata and region queries match reference values.
- No runtime calls to Rcpp,
  [`.Call()`](https://rdrr.io/r/base/CallExternal.html), libBigWig,
  rtracklayer, Rsamtools, or UCSC executables.

## 0.2.x - Native BigWig writer

Status: **implemented and locally validated.**

Scope:

- Canonical signal validation and sorting.
- Chromosome-size validation.
- BigWig fixed header and total summary writer.
- Chromosome B+ tree writer, including multi-level trees.
- bedGraph-like full-data block encoder.
- zlib block compression using base R.
- Full-data R-tree construction, including multi-level trees.
- [`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md)
  for BigWig, WIG, and bedGraph.
- Lazy same-format file copy without decoding.

The first writer release will intentionally omit zoom levels. A valid
BigWig without zoom levels is sufficient for exact interval access and
isolates writer correctness from summary-index complexity.

Acceptance criteria:

- bwTools can read its own BigWig output.
- Reference readers can read the output when available. A libBigWig
  reference-reader format check is part of development validation, while
  the package test suite remains R-only.
- bedGraph -\> BigWig -\> bwTools region queries preserve coordinates
  and values.
- No temporary UCSC conversion executable is invoked.

## 0.3.x - Zoom levels and statistics

Status: **implemented in 0.3.x; public zoom metadata API added in
0.3.1.**

Scope:

- Zoom-level construction.
- Zoom summary blocks and zoom R-tree indexes.
- Mean, min, max, coverage, sum, and standard-deviation interval
  statistics.
- Optional zoom generation during BigWig writing.

Acceptance criteria:

- Zoom and full-data summaries agree within floating-point tolerance.
- Zoom omission remains a supported writer option.
- Exact statistics remain available through `use_zoom = FALSE`.
- Zoom-accelerated results identify their source and reduction level
  explicitly.
- User-facing zoom inspection uses
  [`zoominfo_bwg()`](https://renscq.github.io/bwTools/reference/zoominfo_bwg.md)
  rather than internal metadata structures.

Implementation notes:

- Zoom summary records store covered bases, minimum, maximum, sum, and
  squared sum for each reduction window.
- Coarse lazy-BigWig statistics select the largest stored reduction
  level not exceeding half of the requested bin width, matching the
  selection strategy used by common BigWig readers.
- Partial overlap between a requested bin and a zoom window is
  necessarily approximate; exact scientific calculations should use
  full-resolution data.

### 0.3.4 regression note

WIG chromosome grouping and other internal data.table filters must avoid
bare external variable names that collide with column names. Package
tests now cover chromosome isolation, sample-specific sequence metadata,
and chromosome-specific in-memory statistics.

## 0.4.x - Unified genomic retrieval

Status: **implemented in 0.4.1.**

Scope:

- [`retrieve_bwg()`](https://renscq.github.io/bwTools/reference/retrieve_bwg.md)
  is the only public genomic extraction interface.
- Single-region queries continue to use `chrom`, `start`, and `end`.
- Multi-region queries use a `regions` table with `chrom`, `start`, and
  `end`.
- `result = "data"` preserves the existing data.table return contract.
- `result = "track"` returns a memory-mode `BwgTrack` suitable for
  downstream
  [`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md)
  calls.
- Lazy BigWig extraction remains R-tree indexed and does not load the
  complete source track.
- File output, compression, format conversion, and BigWig zoom
  construction are intentionally excluded from retrieval and remain
  responsibilities of
  [`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md).

Acceptance criteria:

- Existing single-region
  [`retrieve_bwg()`](https://renscq.github.io/bwTools/reference/retrieve_bwg.md)
  calls remain valid.
- Multi-region coordinates are normalized as a genomic union and clipped
  correctly without coordinate rebasing.
- Retrieval never accepts writer-only controls such as output format,
  zoom, or compression.
- A retrieved `BwgTrack` can be passed explicitly to
  [`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md)
  without an intermediate conversion step.
- Large-file retrieval memory use scales with queried BigWig blocks
  rather than complete file size.

Implementation notes:

- Overlapping and directly adjacent query regions are merged before
  retrieval to prevent duplicate signal records.
- Sample selection validates unknown sample IDs rather than silently
  dropping them.
- `result = "track"` preserves selected sample metadata and chromosome
  lengths for the retrieved chromosomes.
- The former `subset_bwg()` API was removed before the package reached a
  stable 1.0 public API because it duplicated
  [`retrieve_bwg()`](https://renscq.github.io/bwTools/reference/retrieve_bwg.md)
  and mixed extraction with persistence.

## 0.5.x - Unified merge semantics

Status: **implemented in 0.5.1.**

Scope:

- [`merge_bwg()`](https://renscq.github.io/bwTools/reference/merge_bwg.md)
  is the only public track-combination interface.
- `method = "auto"` is the default and performs conservative
  non-arithmetic merging.
- Different sample IDs remain independent in auto mode even when
  coordinates overlap.
- Same-sample disjoint regions are appended without requiring a
  user-selected mode.
- Exact duplicates and equal-valued same-sample overlaps are
  consolidated automatically.
- Conflicting same-sample overlaps stop with a clear instruction to
  select `mean` or `sum`.
- `method = "mean"` and `method = "sum"` perform explicit arithmetic
  aggregation.
- Partial overlaps are segmented at all input boundaries before
  arithmetic.
- `missing = "ignore"` is the default mean denominator rule;
  `missing = "zero"` remains explicit.
- `groups` supports replicate- and treatment-level arithmetic
  aggregation.
- Strand-aware aggregation is separate by default and can be explicitly
  ignored.
- Chromosome-size conflicts are errors for all merge methods.
- Merge always returns a memory-mode `BwgTrack`; persistence remains
  exclusively in
  [`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md).

Acceptance criteria:

- A normal call to `merge_bwg(a, b, ...)` does not require users to
  classify inputs as same/different samples or same/different regions.
- Auto mode never averages or sums signal values.
- Different samples are never collapsed solely because their coordinates
  overlap.
- Same-sample exact duplicates and equal-valued overlaps are resolved
  without extra user parameters.
- Conflicting same-sample values cannot be silently resolved.
- Arithmetic methods handle both identical and partially overlapping
  interval boundaries.
- Mean aggregation has explicit uncovered-signal semantics.
- Grouped aggregation preserves group provenance and deterministic
  output sample IDs.
- Conflicting chromosome lengths fail before records are combined.
- Input and output interval counts are recorded in operation metadata.
- `collapse_bwg()` is not exported.

Implementation notes:

- Auto merging canonicalizes signal independently by sample, chromosome,
  and strand.
- Arithmetic aggregation uses a sweep-line event representation over
  1-based closed intervals.
- Atomic segments are emitted only while at least one input member is
  active.
- `drop_zero = TRUE` affects arithmetic methods only; auto mode
  preserves input zero-valued signal.
- `merge_adjacent = TRUE` joins adjacent equal-valued segments after
  auto or arithmetic merging.
- Weighted aggregation and additional arithmetic methods can be added
  later without adding another public merge function.

## 0.6.1 - Public API and BwgTrack contract standardization

Status: **implemented; local R regression testing required.**

Goals:

- Freeze a single downstream-package object contract before GeneTrackR
  adopts bwTools as a backend.
- Standardize primary object arguments on `x`.
- Standardize sample identifiers and selectors on `sample_ids`.
- Use
  [`bwg_track()`](https://renscq.github.io/bwTools/reference/bwg_track.md)
  as the single public constructor.
- Use
  [`detect_bwg_format()`](https://renscq.github.io/bwTools/reference/detect_bwg_format.md)
  as the public format detector.
- Add
  [`validate_bwg()`](https://renscq.github.io/bwTools/reference/validate_bwg.md),
  [`samples_bwg()`](https://renscq.github.io/bwTools/reference/samples_bwg.md),
  and
  [`metadata_bwg()`](https://renscq.github.io/bwTools/reference/metadata_bwg.md)
  for external package integration.
- Require all analysis/accessor functions to consume BwgTrack rather
  than mix object and file-path semantics.
- Keep file input in
  [`read_bwg()`](https://renscq.github.io/bwTools/reference/read_bwg.md)
  and persistence in
  [`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md).
- Standardize schema v2 for samples, signal data, seqinfo, and metadata.
- Make lazy BigWig loading the safe default through `mode = "auto"`.

Acceptance criteria:

- Public API names and argument conventions are internally consistent.
- All public BwgTrack-producing functions return schema v2 objects.
- [`validate_bwg()`](https://renscq.github.io/bwTools/reference/validate_bwg.md)
  detects malformed downstream objects with actionable errors.
- Public accessors do not require direct list-field assumptions.
- README examples use only the standardized API.
- Old `bw_track`, `bw_detect_format`, `sample_names`, and selector
  `samples` names are absent from the public API.
- GeneTrackR can later integrate by depending only on documented bwTools
  exports and `inst/API_CONTRACT.md`.

## 0.6.2 - Binary-safe reader hardening

- Detect BigWig by raw four-byte magic before any text decoding.
- Reject unknown binary inputs before
  [`readLines()`](https://rdrr.io/r/base/readLines.html) /
  [`trimws()`](https://rdrr.io/r/base/trimws.html).
- Keep path values unchanged for file access while normalizing
  display-only path text safely.
- Skip unnecessary filename-stem derivation when explicit `sample_ids`
  are supplied.
- Lock these behaviors with format-detection regression tests.

## 0.7.0 - Large-file performance profiling

Status: **implemented; real large-file benchmark runs required.**

Scope:

- Add
  [`benchmark_bwg()`](https://renscq.github.io/bwTools/reference/benchmark_bwg.md)
  as a file-I/O diagnostic utility without changing the BwgTrack schema
  or normal analysis contract.
- Benchmark package-metadata-cold and cached lazy BigWig loading
  separately.
- Benchmark indexed retrieval over deterministic or user-supplied
  regions.
- Compare zoom-enabled statistics with bounded exact full-resolution
  statistics.
- Optionally benchmark in-memory merge materialization over four
  non-overlapping pieces of the largest selected region.
- Optionally benchmark native BigWig writing for the largest selected
  region with zoom generation disabled and enabled.
- Keep full-file memory loading explicit and opt-in.
- Return raw per-iteration data and aggregate summaries rather than
  writing benchmark reports automatically.

Acceptance criteria:

- The default benchmark suite can run on a large BigWig without
  full-file materialization.
- `stats_full` cannot accidentally run above `full_stats_max_bases`;
  oversized cases are recorded as skipped.
- Benchmark failures are captured in `status` and `message` without
  discarding successful measurements from the same run.
- Benchmark output records wall time, CPU time, approximate R heap
  maximum, returned object size and rows, region throughput, and
  relevant file sizes.
- Cold-cache labels explicitly refer only to the internal bwTools
  metadata cache and do not claim to clear the operating-system page
  cache.
- Writer profiling is bounded by the selected benchmark region rather
  than silently rewriting the complete input file.
- The benchmark utility does not add persistence controls outside
  existing signal-writing responsibilities.

Implementation notes:

- Default centered windows are 10 kb, 1 Mb, and 10 Mb on the largest
  known chromosome; representative user-supplied regions are preferred
  for real profiling.
- Peak memory uses `gc(reset = TRUE)` and therefore approximates R heap
  usage, not total process RSS.
- `read_memory` is available only as an explicit operation because
  compressed BigWig size is not a safe predictor of in-memory
  interval-table size.
- Benchmark output is a `bwToolsBenchmark` list with `system`, `file`,
  `config`, `regions`, `results`, and `summary` components.

## 0.7.1 - Benchmark-driven optimization

Status: **implemented; repeat real large-file benchmark required.**

Measured bottlenecks from the first real dense BigWig benchmark:

- Lazy metadata loading was already fast and is unchanged.
- Zoom statistics were already fast and are unchanged.
- Full-resolution statistics spent substantially more time and heap than
  indexed retrieval because the old path materialized all intervals and
  then looped across them in R.
- Large indexed retrieval showed acceptable elapsed time but excessive
  temporary heap relative to the final returned signal table.

Implemented optimizations:

- Stream full-resolution BigWig blocks directly into statistical bin
  accumulators without constructing an intermediate signal table.
- Vectorize accumulation for intervals contained within one output bin
  and loop only across the much smaller set of intervals spanning bin
  boundaries.
- Decode full-data blocks into reusable coordinate/value vectors.
- Materialize indexed retrieval into growable typed vectors and build
  one final data.table instead of one data.table per block followed by
  `rbindlist()`.
- Skip redundant list binding and sorting for single-result retrieval.
- Avoid copying the complete memory-mode signal table before region
  filtering.

Acceptance criteria:

- Exact streaming statistics match memory-mode exact statistics for all
  six public statistics.
- Indexed retrieval matches the bundled bedGraph reference after
  optimization.
- Public API, BwgTrack schema v2, coordinate conventions, zoom behavior,
  and benchmark output contract remain unchanged.
- Repeat the same real-file 0.7.0 benchmark configuration and verify
  lower `stats_full` elapsed time/heap and lower retrieval heap without
  regressions in `read_lazy` or `stats_zoom`.
- Defer streaming-writer changes until write benchmarks demonstrate a
  concrete bottleneck.

## 0.7.2 - Indexed retrieval and memory-diagnostic refinement

Status: **implemented and validated on a dense real BigWig benchmark.**

Measured 0.7.1 outcome:

- Block-streamed exact statistics produced the intended large speedup
  and are retained unchanged.
- Growable-vector indexed retrieval did not consistently improve elapsed
  time and could increase temporary R-heap use because repeated vector
  extension may trigger reallocations and copies.
- Existing benchmark peak-heap metrics did not distinguish retained live
  memory from historical maximum R-heap use.

Implemented refinements:

- Store decoded BigWig query blocks as primitive `start`, `end`, and
  `value` vector chunks.
- Flatten each vector family once after decoding rather than repeatedly
  growing large R vectors.
- Construct only one final retrieval data.table and retain the existing
  single-result fast path.
- Record live R heap after garbage collection separately from the
  maximum heap observed since `gc(reset = TRUE)`.
- Report `peak_heap_delta_mb` and `live_heap_delta_mb` separately while
  keeping `heap_delta_mb` as a compatibility alias for the peak delta.

Acceptance criteria:

- Indexed retrieval remains signal-identical to the bundled bedGraph
  reference.
- Streaming exact statistics remain unchanged and continue to match
  memory-mode exact statistics.
- Dense real-file benchmarking confirmed retrieval around 0.28 s for 1
  Mb and 1.85 s for an approximately 9.8 Mb chromosome-wide query on the
  test system.
- The 0.7.1 full-statistics improvement remained stable at approximately
  0.30 s for a 1 Mb / 1000-bin exact mean benchmark.
- Live-heap deltas after GC closely tracked final retrieval object size,
  indicating transient allocation pressure rather than persistent
  retention.
- Reader/retrieval optimization is therefore frozen for the 0.7.3 writer
  profiling stage unless later benchmarks demonstrate a new regression.

## 0.7.4 - Writer benchmark safety and diagnostics

Status: **implemented; real-file rerun required.**

Goals:

- Make writer profiling finish predictably on dense BigWig signal
  tracks.
- Run general benchmark warmups once per benchmark case rather than once
  per timed iteration.
- Separate expensive writer repetition settings from general
  read/statistics benchmark repetition settings.
- Bound writer work by both genomic span and actual interval count.
- Apply a stricter density bound to `zoom_on`, whose current pure-R
  builder rescans signal intervals for each zoom level.
- Print concise progress messages for expensive writer cases.

Acceptance criteria:

- `iterations = 3, warmup = 1` no longer multiplies writer work when the
  caller leaves writer-specific controls at their safe defaults.
- Default writer profiling performs one timed write and zero writer
  warmups per eligible variant/region case.
- A dense region may run `zoom_off` while independently skipping
  `zoom_on` when it exceeds `write_zoom_max_intervals`.
- Skipped writer cases remain visible in benchmark results with
  actionable messages.
- Writer timing still excludes region retrieval/materialization.
- No BigWig binary-format or public BwgTrack contract changes are
  introduced.

## 0.7.3 - Native BigWig writer benchmark and profiling

Status: **implemented; real large-file writer benchmark required.**

Measured 0.7.2 outcome:

- Indexed retrieval returned to the 0.7.0 performance baseline while
  retaining the 0.7.1 block-streamed exact-statistics speedup.
- Live-heap diagnostics showed that large retrievals retain
  approximately the final result size after GC, while the remaining
  concern is transient writer and retrieval allocation rather than
  persistent memory retention.
- Reader/retrieval optimization can therefore pause until a later
  benchmark demonstrates a concrete regression or bottleneck.

Implemented writer profiling:

- Profile every selected benchmark region when `operations` includes
  `write`.
- Compare native BigWig writer behavior with zoom disabled and enabled.
- Exclude indexed retrieval/materialization from the timed writer
  expression so the result reflects
  [`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md)
  rather than read performance.
- Add an explicit `write_max_bases` safety limit; oversized writer cases
  are recorded as skipped.
- Record materialized input size, interval count, covered bases,
  estimated uncompressed type-1 data payload, and expected primary
  data-block count.
- Record output file size, interval throughput, input/payload
  throughput, output throughput, and output-size ratios.
- Derive zoom time overhead and zoom output-size overhead from matching
  `zoom_off`/`zoom_on` summary rows.

Acceptance criteria:

- Writer benchmark output remains deterministic and valid across all
  selected regions within `write_max_bases`.
- `zoom_off` and `zoom_on` output files are valid BigWig files and the
  benchmark captures both variants without timing indexed retrieval.
- Real-file profiling should use at least two region sizes, preferably
  around 1 Mb and 5 Mb or 10 Mb, on a dense signal track.
- Use writer timing, intervals/s, peak/live heap, output size, and zoom
  overhead to decide whether 0.7.5 needs streaming/chunked writer
  changes.
- Do not optimize writer internals until the 0.7.3 real-file benchmark
  identifies a concrete bottleneck.

## 0.7.5 - Hierarchical native zoom-writer optimization

Status: **implemented and validated on dense real-file 1 Mb and 5 Mb
writer benchmarks.**

Measured 0.7.4 outcome:

- Primary `zoom_off` writing scaled acceptably on the dense real signal
  track.
- The 1 Mb writer case completed quickly without zoom summaries, while
  enabling zoom added orders of magnitude more elapsed time.
- The additional file size from zoom summaries was modest relative to
  the elapsed-time overhead, identifying zoom construction as the
  dominant writer bottleneck.
- Primary data-block encoding, compression, chromosome-tree writing, and
  the full-data R-tree therefore remain unchanged in 0.7.5.

Implemented optimizations:

- Build the first zoom level with vectorized same-window aggregation.
- Split only intervals that actually cross a first-level zoom boundary.
- Build every higher zoom level from the preceding level’s summary
  records instead of rescanning all full-resolution intervals.
- Aggregate `valid_count`, `sum_data`, and `sum_squared` by summation
  and preserve `min_value`/`max_value` by reduction over child
  summaries.
- Retain only the preceding zoom level while writing the next level so
  the writer does not keep a complete zoom pyramid resident at once.
- Keep the public writer API, selected zoom levels, binary record
  layout, BwgTrack schema v2, and benchmark interface unchanged.

Acceptance criteria:

- Hierarchical summaries match direct full-signal summaries for every
  stored statistic within floating-point tolerance.
- Boundary-crossing signal intervals are split correctly at first-level
  zoom window boundaries and chromosome ends.
- BigWig zoom metadata and zoom queries remain readable through existing
  bwTools reader paths.
- `zoom = FALSE` writer performance does not regress because the primary
  writer path is untouched.
- Rerun the bounded real-file writer benchmark and compare `zoom_on`
  elapsed time, peak/live heap, output size, and zoom overhead against
  the 0.7.4 baseline before increasing default zoom benchmark density
  limits.

## 0.8.0 - Public API stability and cross-format hardening

Status: **implemented; local R package regression testing required.**

Scope:

- Keep the optimized 0.7.5 BigWig reader, indexed retrieval, streamed
  exact statistics, and hierarchical zoom writer unchanged.
- Treat the 15 exported function signatures as the downstream
  integration boundary and regression-test their parameter names and
  order.
- Standardize package-generated user errors under the `bwTools_error`
  condition class so downstream packages can catch bwTools failures
  without matching full message text.
- Centralize common scalar, flag, positive-integer, choice, and
  genomic-interval validation helpers.
- Reject fractional statistics coordinates rather than silently
  truncating them.
- Verify full-resolution signal equivalence across BigWig zoom/no-zoom,
  WIG, WIG.gz, bedGraph, and bedGraph.gz round trips.
- Keep BwgTrack schema v2 unchanged.
- Reorganize README documentation around the normal read -\> retrieve
  -\> merge/stats -\> write workflow, with benchmarking clearly
  separated as an advanced diagnostic feature.

Acceptance criteria:

- All exported function names and formals match the 0.8.0 API contract.
- BwgTrack schema v2 core fields and data types remain unchanged.
- Package-generated input-validation errors inherit from
  `bwTools_error`.
- Cross-format round trips preserve standardized coordinates, values,
  and sample identity within floating-point tolerance.
- BigWig zoom/no-zoom output preserves identical full-resolution signal
  while differing only in zoom metadata as requested.
- Writer manifests retain `sample_id`, `file`, and `format` core
  columns.
- No new runtime dependency, compiled code, or external executable is
  added.

## 0.8.x - Stability follow-ups

Planned only when testing identifies a concrete need:

- Additional malformed-file or extreme-coordinate regression coverage.
- Further error subclasses if downstream integration requires
  finer-grained condition handling.
- Documentation refinements that do not change the public API contract.

## 0.9.x - Release hardening

Planned:

- Freeze the public API.
- Normalize error classes/messages and documentation.
- Run Linux, macOS, and Windows compatibility checks.
- Complete large-file and malformed-file regression coverage.

## 1.0.0 - Stable release

Release only after reader, writer, retrieval, merge, statistics, zoom,
object contract, large-file behavior, and cross-platform installation
are validated.

## 0.6.3 - Automatic format dispatch regression fix

Status: implemented; local R regression testing required.

- Keep the 0.6.1 public API and BwgTrack schema v2 unchanged.
- Remove names from the vector returned by automatic per-file format
  detection.
- Extract each detected format as an attribute-free scalar before
  dispatch.
- Add a regression test that exercises
  [`read_bwg()`](https://renscq.github.io/bwTools/reference/read_bwg.md)
  with automatic BigWig detection and explicit `sample_ids`.
