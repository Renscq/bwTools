# Package index

## File I/O

Read, detect, and write genomic signal files.

- [`detect_bwg_format()`](https://renscq.github.io/bwTools/reference/detect_bwg_format.md)
  : Detect genomic signal file format
- [`read_bwg()`](https://renscq.github.io/bwTools/reference/read_bwg.md)
  : Read genomic signal tracks
- [`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md)
  : Write genomic signal tracks

## BwgTrack contract

Construct, validate, inspect, and summarize standardized tracks.

- [`bwg_track()`](https://renscq.github.io/bwTools/reference/bwg_track.md)
  : Construct a standardized BwgTrack object
- [`is_bwg_track()`](https://renscq.github.io/bwTools/reference/is_bwg_track.md)
  : Test whether an object follows the BwgTrack contract
- [`validate_bwg()`](https://renscq.github.io/bwTools/reference/validate_bwg.md)
  : Validate a BwgTrack object
- [`samples_bwg()`](https://renscq.github.io/bwTools/reference/samples_bwg.md)
  : Return sample metadata
- [`seqinfo_bwg()`](https://renscq.github.io/bwTools/reference/seqinfo_bwg.md)
  : Return sample-specific chromosome metadata
- [`metadata_bwg()`](https://renscq.github.io/bwTools/reference/metadata_bwg.md)
  : Return BwgTrack metadata
- [`summary_bwg()`](https://renscq.github.io/bwTools/reference/summary_bwg.md)
  : Summarize a BwgTrack object
- [`zoominfo_bwg()`](https://renscq.github.io/bwTools/reference/zoominfo_bwg.md)
  : Return BigWig zoom-level metadata

## Signal analysis

Retrieve, summarize, merge, and aggregate genomic signal.

- [`retrieve_bwg()`](https://renscq.github.io/bwTools/reference/retrieve_bwg.md)
  : Retrieve signal intervals
- [`stats_bwg()`](https://renscq.github.io/bwTools/reference/stats_bwg.md)
  : Calculate interval statistics
- [`merge_bwg()`](https://renscq.github.io/bwTools/reference/merge_bwg.md)
  : Merge genomic signal tracks

## Diagnostics

Profile native BigWig workflows.

- [`benchmark_bwg()`](https://renscq.github.io/bwTools/reference/benchmark_bwg.md)
  : Benchmark a large BigWig workflow
