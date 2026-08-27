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
