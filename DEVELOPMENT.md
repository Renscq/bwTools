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

## 0.4.x - Subset and file transformation

Status: **implemented in 0.4.0.**

Scope:

- `retrieve_bwg()` for single-region retrieval.
- `subset_bwg()` for multi-region in-memory and file-to-file extraction.
- Indexed lazy BigWig extraction without loading the complete source track.
- Cross-format subset output to BigWig, WIG, or bedGraph.
- Explicit empty-sample behavior for file output.

Acceptance criteria:

- Subset coordinates are clipped correctly.
- Empty regions produce valid empty results or explicit output behavior.
- Large-file subset memory use scales with the selected intervals and queried
  BigWig blocks, not the complete source file size.

Implementation notes:

- Regions use 1-based closed coordinates and are treated as a genomic union.
- Overlapping and directly adjacent regions are merged before retrieval.
- Signal intervals crossing a subset boundary are clipped without coordinate
  rebasing.
- Lazy BigWig sources are queried through the native R-tree index for each
  normalized region.
- The 0.4.0 writer path materializes only the selected subset before calling the
  existing writer; a whole source track is never loaded solely for subsetting.
- Empty in-memory subsets are valid. File output errors by default when a
  selected sample is empty, or records it explicitly when `empty = "skip"`.

## 0.5.x - Merge semantics

Scope:

- `merge_bwg()` to combine independent samples into one multi-sample BwgTrack without changing values.
- `collapse_bwg()` to mathematically combine signals into one output track with an explicit overlap rule such as sum, mean, min, or max.
- Chromosome-size and coordinate-contract conflict checks.
- Streaming merge where feasible.

Acceptance criteria:

- Sample merge never silently averages or drops signal.
- Arithmetic collapse always requires an explicit aggregation rule when overlaps exist.
- Input counts and resulting interval counts are reported by validation helpers.

## 0.6.x - GeneTrackR integration

Scope:

- GeneTrackR imports bwTools for signal I/O.
- Existing GeneTrackR public `read_bwg()`, `retrieve_bwg()`, `merge_bwg()`, and `write_bwg()` interfaces are retained as compatibility wrappers where practical.
- Remove GeneTrackR's libBigWig/Rcpp signal backend and UCSC BigWig writer dependency after regression validation.

Acceptance criteria:

- Existing `plot_signal_gene()`, `plot_signal_transcript()`, `plot_signal_region()`, and `plot_tracks()` consume bwTools objects without conversion.
- Existing GeneTrackR signal examples run with the same object contract.

## 1.0.0 - Stable release

Scope:

- Windows and Linux regression coverage.
- Corrupt/truncated-file diagnostics.
- Large-file benchmarks.
- Public API stabilization and documentation.
- Optional remote-range access only if it can be implemented without compromising the no-external-binary design.

## 0.2.0 development validation note

The BigWig v4 binary layout used by the native writer was independently checked
against the supplied libBigWig C reader during development. Both a standard
three-chromosome file and a synthetic 300-chromosome file exercising multi-level
chromosome B+ trees and R-trees were readable by the reference C implementation.
This validates the binary layout independently of bwTools' own reader. The R
package tests remain the required runtime validation for the actual R writer.

# Documentation conventions

All language-tagged code fences in Markdown and QMD files must use braced
Quarto/knitr syntax, for example `{r}`, `{bash}`, and `{text}`. Unbraced
forms such as ` ```r`, ` ```bash`, and ` ```text` are not permitted.

