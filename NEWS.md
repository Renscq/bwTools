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
