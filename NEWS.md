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

## Current limitations

- BigWig writing is intentionally deferred to the next development stage.
- Remote HTTP/HTTPS/FTP bigWig access is not implemented in 0.1.0.
- Lazy access is currently supported for bigWig only.
