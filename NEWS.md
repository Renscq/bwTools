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
- Copy user-supplied `data.table` inputs in `bw_track()` before normalization to avoid modifying caller-owned objects by reference.
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
