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

## 0.7.x - Large-file performance profiling

Planned:

- Benchmark lazy retrieval, statistics, merge materialization, native writing,
  and zoom generation on large RNA-seq and Ribo-seq BigWig files.
- Optimize only measured bottlenecks.
- Introduce streaming writer changes if peak memory requires them.

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

