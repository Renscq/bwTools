# bwTools development roadmap

## Design contract

- Package code for BigWig parsing and writing is implemented in R only.
- No bundled C/C++ source, Rcpp, libBigWig linkage, or UCSC executable is used at runtime.
- `data.table` is retained as the only runtime package dependency because GeneTrackR already uses it and its signal APIs consume `data.table` objects.
- Public in-memory coordinates are 1-based closed.
- BigWig and bedGraph file coordinates are converted from 0-based half-open at the I/O boundary.
- Large BigWig files should use lazy indexed access rather than full-memory loading.
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
- Bundled BigWig/bedGraph/chrom.sizes example dataset for integration and writer round-trip tests.

Acceptance criteria:

- Coordinate conversion is deterministic and covered by tests.
- BigWig fixture metadata and region queries match reference values.
- No runtime calls to Rcpp, `.Call()`, libBigWig, rtracklayer, Rsamtools, or UCSC executables.

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
- `write_bwg()` for BigWig, WIG, and bedGraph.
- Lazy same-format file copy without decoding.

The first writer release will intentionally omit zoom levels. A valid BigWig without zoom levels is sufficient for exact interval access and isolates writer correctness from summary-index complexity.

Acceptance criteria:

- bwTools can read its own BigWig output.
- Reference readers can read the output when available. A libBigWig reference-reader format check is part of development validation, while the package test suite remains R-only.
- bedGraph -> BigWig -> bwTools region queries preserve coordinates and values.
- No temporary UCSC conversion executable is invoked.

## 0.3.x - Zoom levels and statistics

Status: **implemented in 0.3.x; public zoom metadata API added in 0.3.1.**

Scope:

- Zoom-level construction.
- Zoom summary blocks and zoom R-tree indexes.
- Mean, min, max, coverage, sum, and standard-deviation interval statistics.
- Optional zoom generation during BigWig writing.

Acceptance criteria:

- Zoom and full-data summaries agree within floating-point tolerance.
- Zoom omission remains a supported writer option.
- Exact statistics remain available through `use_zoom = FALSE`.
- Zoom-accelerated results identify their source and reduction level explicitly.
- User-facing zoom inspection uses `zoominfo_bwg()` rather than internal metadata structures.

Implementation notes:

- Zoom summary records store covered bases, minimum, maximum, sum, and squared sum for each reduction window.
- Coarse lazy-BigWig statistics select the largest stored reduction level not exceeding half of the requested bin width, matching the selection strategy used by common BigWig readers.
- Partial overlap between a requested bin and a zoom window is necessarily approximate; exact scientific calculations should use full-resolution data.

### 0.3.4 regression note

WIG chromosome grouping and other internal data.table filters must avoid bare
external variable names that collide with column names. Package tests now cover
chromosome isolation, sample-specific sequence metadata, and chromosome-specific
in-memory statistics.

## 0.4.x - Unified genomic retrieval

Status: **implemented in 0.4.1.**

Scope:

- `retrieve_bwg()` is the only public genomic extraction interface.
- Single-region queries continue to use `chrom`, `start`, and `end`.
- Multi-region queries use a `regions` table with `chrom`, `start`, and `end`.
- `result = "data"` preserves the existing data.table return contract.
- `result = "track"` returns a memory-mode `BwgTrack` suitable for downstream
  `write_bwg()` calls.
- Lazy BigWig extraction remains R-tree indexed and does not load the complete
  source track.
- File output, compression, format conversion, and BigWig zoom construction are
  intentionally excluded from retrieval and remain responsibilities of
  `write_bwg()`.

Acceptance criteria:

- Existing single-region `retrieve_bwg()` calls remain valid.
- Multi-region coordinates are normalized as a genomic union and clipped
  correctly without coordinate rebasing.
- Retrieval never accepts writer-only controls such as output format, zoom, or
  compression.
- A retrieved `BwgTrack` can be passed explicitly to `write_bwg()` without an
  intermediate conversion step.
- Large-file retrieval memory use scales with queried BigWig blocks rather than
  complete file size.

Implementation notes:

- Overlapping and directly adjacent query regions are merged before retrieval to
  prevent duplicate signal records.
- Sample selection validates unknown sample IDs rather than silently dropping
  them.
- `result = "track"` preserves selected sample metadata and chromosome lengths
  for the retrieved chromosomes.
- The former `subset_bwg()` API was removed before the package reached a stable
  1.0 public API because it duplicated `retrieve_bwg()` and mixed extraction
  with persistence.

## 0.5.x - Unified merge semantics

Status: **implemented in 0.5.1.**

Scope:

- `merge_bwg()` is the only public track-combination interface.
- `method = "auto"` is the default and performs conservative non-arithmetic merging.
- Different sample IDs remain independent in auto mode even when coordinates overlap.
- Same-sample disjoint regions are appended without requiring a user-selected mode.
- Exact duplicates and equal-valued same-sample overlaps are consolidated automatically.
- Conflicting same-sample overlaps stop with a clear instruction to select `mean` or `sum`.
- `method = "mean"` and `method = "sum"` perform explicit arithmetic aggregation.
- Partial overlaps are segmented at all input boundaries before arithmetic.
- `missing = "ignore"` is the default mean denominator rule; `missing = "zero"` remains explicit.
- `groups` supports replicate- and treatment-level arithmetic aggregation.
- Strand-aware aggregation is separate by default and can be explicitly ignored.
- Chromosome-size conflicts are errors for all merge methods.
- Merge always returns a memory-mode `BwgTrack`; persistence remains exclusively in `write_bwg()`.

Acceptance criteria:

- A normal call to `merge_bwg(a, b, ...)` does not require users to classify inputs as same/different samples or same/different regions.
- Auto mode never averages or sums signal values.
- Different samples are never collapsed solely because their coordinates overlap.
- Same-sample exact duplicates and equal-valued overlaps are resolved without extra user parameters.
- Conflicting same-sample values cannot be silently resolved.
- Arithmetic methods handle both identical and partially overlapping interval boundaries.
- Mean aggregation has explicit uncovered-signal semantics.
- Grouped aggregation preserves group provenance and deterministic output sample IDs.
- Conflicting chromosome lengths fail before records are combined.
- Input and output interval counts are recorded in operation metadata.
- `collapse_bwg()` is not exported.

Implementation notes:

- Auto merging canonicalizes signal independently by sample, chromosome, and strand.
- Arithmetic aggregation uses a sweep-line event representation over 1-based closed intervals.
- Atomic segments are emitted only while at least one input member is active.
- `drop_zero = TRUE` affects arithmetic methods only; auto mode preserves input zero-valued signal.
- `merge_adjacent = TRUE` joins adjacent equal-valued segments after auto or arithmetic merging.
- Weighted aggregation and additional arithmetic methods can be added later without adding another public merge function.

## 0.6.1 - Public API and BwgTrack contract standardization

Status: **implemented; local R regression testing required.**

Goals:

- Freeze a single downstream-package object contract before GeneTrackR adopts
  bwTools as a backend.
- Standardize primary object arguments on `x`.
- Standardize sample identifiers and selectors on `sample_ids`.
- Use `bwg_track()` as the single public constructor.
- Use `detect_bwg_format()` as the public format detector.
- Add `validate_bwg()`, `samples_bwg()`, and `metadata_bwg()` for external package integration.
- Require all analysis/accessor functions to consume BwgTrack rather than mix
  object and file-path semantics.
- Keep file input in `read_bwg()` and persistence in `write_bwg()`.
- Standardize schema v2 for samples, signal data, seqinfo, and metadata.
- Make lazy BigWig loading the safe default through `mode = "auto"`.

Acceptance criteria:

- Public API names and argument conventions are internally consistent.
- All public BwgTrack-producing functions return schema v2 objects.
- `validate_bwg()` detects malformed downstream objects with actionable errors.
- Public accessors do not require direct list-field assumptions.
- README examples use only the standardized API.
- Old `bw_track`, `bw_detect_format`, `sample_names`, and selector `samples`
  names are absent from the public API.
- GeneTrackR can later integrate by depending only on documented bwTools
  exports and `inst/API_CONTRACT.md`.

## 0.6.2 - Binary-safe reader hardening

- Detect BigWig by raw four-byte magic before any text decoding.
- Reject unknown binary inputs before `readLines()` / `trimws()`.
- Keep path values unchanged for file access while normalizing display-only
  path text safely.
- Skip unnecessary filename-stem derivation when explicit `sample_ids` are
  supplied.
- Lock these behaviors with format-detection regression tests.

## 0.7.0 - Large-file performance profiling

Status: **implemented; real large-file benchmark runs required.**

Scope:

- Add `benchmark_bwg()` as a file-I/O diagnostic utility without changing the
  BwgTrack schema or normal analysis contract.
- Benchmark package-metadata-cold and cached lazy BigWig loading separately.
- Benchmark indexed retrieval over deterministic or user-supplied regions.
- Compare zoom-enabled statistics with bounded exact full-resolution
  statistics.
- Optionally benchmark in-memory merge materialization over four non-overlapping
  pieces of the largest selected region.
- Optionally benchmark native BigWig writing for the largest selected region
  with zoom generation disabled and enabled.
- Keep full-file memory loading explicit and opt-in.
- Return raw per-iteration data and aggregate summaries rather than writing
  benchmark reports automatically.

Acceptance criteria:

- The default benchmark suite can run on a large BigWig without full-file
  materialization.
- `stats_full` cannot accidentally run above `full_stats_max_bases`; oversized
  cases are recorded as skipped.
- Benchmark failures are captured in `status` and `message` without discarding
  successful measurements from the same run.
- Benchmark output records wall time, CPU time, approximate R heap maximum,
  returned object size and rows, region throughput, and relevant file sizes.
- Cold-cache labels explicitly refer only to the internal bwTools metadata cache
  and do not claim to clear the operating-system page cache.
- Writer profiling is bounded by the selected benchmark region rather than
  silently rewriting the complete input file.
- The benchmark utility does not add persistence controls outside existing
  signal-writing responsibilities.

Implementation notes:

- Default centered windows are 10 kb, 1 Mb, and 10 Mb on the largest known
  chromosome; representative user-supplied regions are preferred for real
  profiling.
- Peak memory uses `gc(reset = TRUE)` and therefore approximates R heap usage,
  not total process RSS.
- `read_memory` is available only as an explicit operation because compressed
  BigWig size is not a safe predictor of in-memory interval-table size.
- Benchmark output is a `bwToolsBenchmark` list with `system`, `file`, `config`,
  `regions`, `results`, and `summary` components.

## 0.7.1 - Benchmark-driven optimization

Status: **implemented; repeat real large-file benchmark required.**

Measured bottlenecks from the first real dense BigWig benchmark:

- Lazy metadata loading was already fast and is unchanged.
- Zoom statistics were already fast and are unchanged.
- Full-resolution statistics spent substantially more time and heap than
  indexed retrieval because the old path materialized all intervals and then
  looped across them in R.
- Large indexed retrieval showed acceptable elapsed time but excessive temporary
  heap relative to the final returned signal table.

Implemented optimizations:

- Stream full-resolution BigWig blocks directly into statistical bin
  accumulators without constructing an intermediate signal table.
- Vectorize accumulation for intervals contained within one output bin and loop
  only across the much smaller set of intervals spanning bin boundaries.
- Decode full-data blocks into reusable coordinate/value vectors.
- Materialize indexed retrieval into growable typed vectors and build one final
  data.table instead of one data.table per block followed by `rbindlist()`.
- Skip redundant list binding and sorting for single-result retrieval.
- Avoid copying the complete memory-mode signal table before region filtering.

Acceptance criteria:

- Exact streaming statistics match memory-mode exact statistics for all six
  public statistics.
- Indexed retrieval matches the bundled bedGraph reference after optimization.
- Public API, BwgTrack schema v2, coordinate conventions, zoom behavior, and
  benchmark output contract remain unchanged.
- Repeat the same real-file 0.7.0 benchmark configuration and verify lower
  `stats_full` elapsed time/heap and lower retrieval heap without regressions in
  `read_lazy` or `stats_zoom`.
- Defer streaming-writer changes until write benchmarks demonstrate a concrete
  bottleneck.

## 0.7.2 - Indexed retrieval and memory-diagnostic refinement

Status: **implemented and validated on a dense real BigWig benchmark.**

Measured 0.7.1 outcome:

- Block-streamed exact statistics produced the intended large speedup and are
  retained unchanged.
- Growable-vector indexed retrieval did not consistently improve elapsed time
  and could increase temporary R-heap use because repeated vector extension may
  trigger reallocations and copies.
- Existing benchmark peak-heap metrics did not distinguish retained live memory
  from historical maximum R-heap use.

Implemented refinements:

- Store decoded BigWig query blocks as primitive `start`, `end`, and `value`
  vector chunks.
- Flatten each vector family once after decoding rather than repeatedly growing
  large R vectors.
- Construct only one final retrieval data.table and retain the existing
  single-result fast path.
- Record live R heap after garbage collection separately from the maximum heap
  observed since `gc(reset = TRUE)`.
- Report `peak_heap_delta_mb` and `live_heap_delta_mb` separately while keeping
  `heap_delta_mb` as a compatibility alias for the peak delta.

Acceptance criteria:

- Indexed retrieval remains signal-identical to the bundled bedGraph reference.
- Streaming exact statistics remain unchanged and continue to match memory-mode
  exact statistics.
- Dense real-file benchmarking confirmed retrieval around 0.28 s for 1 Mb and
  1.85 s for an approximately 9.8 Mb chromosome-wide query on the test system.
- The 0.7.1 full-statistics improvement remained stable at approximately 0.30 s
  for a 1 Mb / 1000-bin exact mean benchmark.
- Live-heap deltas after GC closely tracked final retrieval object size,
  indicating transient allocation pressure rather than persistent retention.
- Reader/retrieval optimization is therefore frozen for the 0.7.3 writer
  profiling stage unless later benchmarks demonstrate a new regression.

## 0.7.4 - Writer benchmark safety and diagnostics

Status: **implemented; real-file rerun required.**

Goals:

- Make writer profiling finish predictably on dense BigWig signal tracks.
- Run general benchmark warmups once per benchmark case rather than once per
  timed iteration.
- Separate expensive writer repetition settings from general read/statistics
  benchmark repetition settings.
- Bound writer work by both genomic span and actual interval count.
- Apply a stricter density bound to `zoom_on`, whose current pure-R builder
  rescans signal intervals for each zoom level.
- Print concise progress messages for expensive writer cases.

Acceptance criteria:

- `iterations = 3, warmup = 1` no longer multiplies writer work when the caller
  leaves writer-specific controls at their safe defaults.
- Default writer profiling performs one timed write and zero writer warmups per
  eligible variant/region case.
- A dense region may run `zoom_off` while independently skipping `zoom_on` when
  it exceeds `write_zoom_max_intervals`.
- Skipped writer cases remain visible in benchmark results with actionable
  messages.
- Writer timing still excludes region retrieval/materialization.
- No BigWig binary-format or public BwgTrack contract changes are introduced.

## 0.7.3 - Native BigWig writer benchmark and profiling

Status: **implemented; real large-file writer benchmark required.**

Measured 0.7.2 outcome:

- Indexed retrieval returned to the 0.7.0 performance baseline while retaining
  the 0.7.1 block-streamed exact-statistics speedup.
- Live-heap diagnostics showed that large retrievals retain approximately the
  final result size after GC, while the remaining concern is transient writer and
  retrieval allocation rather than persistent memory retention.
- Reader/retrieval optimization can therefore pause until a later benchmark
  demonstrates a concrete regression or bottleneck.

Implemented writer profiling:

- Profile every selected benchmark region when `operations` includes `write`.
- Compare native BigWig writer behavior with zoom disabled and enabled.
- Exclude indexed retrieval/materialization from the timed writer expression so
  the result reflects `write_bwg()` rather than read performance.
- Add an explicit `write_max_bases` safety limit; oversized writer cases are
  recorded as skipped.
- Record materialized input size, interval count, covered bases, estimated
  uncompressed type-1 data payload, and expected primary data-block count.
- Record output file size, interval throughput, input/payload throughput, output
  throughput, and output-size ratios.
- Derive zoom time overhead and zoom output-size overhead from matching
  `zoom_off`/`zoom_on` summary rows.

Acceptance criteria:

- Writer benchmark output remains deterministic and valid across all selected
  regions within `write_max_bases`.
- `zoom_off` and `zoom_on` output files are valid BigWig files and the benchmark
  captures both variants without timing indexed retrieval.
- Real-file profiling should use at least two region sizes, preferably around
  1 Mb and 5 Mb or 10 Mb, on a dense signal track.
- Use writer timing, intervals/s, peak/live heap, output size, and zoom overhead
  to decide whether 0.7.5 needs streaming/chunked writer changes.
- Do not optimize writer internals until the 0.7.3 real-file benchmark identifies
  a concrete bottleneck.

## 0.7.5 - Hierarchical native zoom-writer optimization

Status: **implemented; real-file writer benchmark rerun required.**

Measured 0.7.4 outcome:

- Primary `zoom_off` writing scaled acceptably on the dense real signal track.
- The 1 Mb writer case completed quickly without zoom summaries, while enabling
  zoom added orders of magnitude more elapsed time.
- The additional file size from zoom summaries was modest relative to the
  elapsed-time overhead, identifying zoom construction as the dominant writer
  bottleneck.
- Primary data-block encoding, compression, chromosome-tree writing, and the
  full-data R-tree therefore remain unchanged in 0.7.5.

Implemented optimizations:

- Build the first zoom level with vectorized same-window aggregation.
- Split only intervals that actually cross a first-level zoom boundary.
- Build every higher zoom level from the preceding level's summary records
  instead of rescanning all full-resolution intervals.
- Aggregate `valid_count`, `sum_data`, and `sum_squared` by summation and
  preserve `min_value`/`max_value` by reduction over child summaries.
- Retain only the preceding zoom level while writing the next level so the
  writer does not keep a complete zoom pyramid resident at once.
- Keep the public writer API, selected zoom levels, binary record layout,
  BwgTrack schema v2, and benchmark interface unchanged.

Acceptance criteria:

- Hierarchical summaries match direct full-signal summaries for every stored
  statistic within floating-point tolerance.
- Boundary-crossing signal intervals are split correctly at first-level zoom
  window boundaries and chromosome ends.
- BigWig zoom metadata and zoom queries remain readable through existing
  bwTools reader paths.
- `zoom = FALSE` writer performance does not regress because the primary writer
  path is untouched.
- Rerun the bounded real-file writer benchmark and compare `zoom_on` elapsed
  time, peak/live heap, output size, and zoom overhead against the 0.7.4
  baseline before increasing default zoom benchmark density limits.

## 0.8.x - Statistics and zoom refinements

Planned:

- Refine statistics or zoom behavior only where real workflows justify it.
- Preserve exact/full-resolution and approximate/zoom semantics explicitly.

## 0.9.x - Release hardening

Planned:

- Freeze the public API.
- Normalize error classes/messages and documentation.
- Run Linux, macOS, and Windows compatibility checks.
- Complete large-file and malformed-file regression coverage.

## 1.0.0 - Stable release

Release only after reader, writer, retrieval, merge, statistics, zoom, object
contract, large-file behavior, and cross-platform installation are validated.

## 0.6.3 - Automatic format dispatch regression fix

Status: implemented; local R regression testing required.

- Keep the 0.6.1 public API and BwgTrack schema v2 unchanged.
- Remove names from the vector returned by automatic per-file format detection.
- Extract each detected format as an attribute-free scalar before dispatch.
- Add a regression test that exercises `read_bwg()` with automatic BigWig
  detection and explicit `sample_ids`.

