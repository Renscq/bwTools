# Getting started with bwTools

## Overview

`bwTools` uses one stable data flow for BigWig, WIG, and bedGraph signal
tracks:

    file -> read_bwg() -> BwgTrack -> analysis -> write_bwg()

Public coordinates are always **1-based closed**.

## Load the bundled example

``` r

library(bwTools)

bw_file <- system.file(
  "extdata",
  "bwtools_example.bigwig",
  package = "bwTools",
  mustWork = TRUE
)

bw <- read_bwg(bw_file, sample_ids = "example")
```

A BigWig-only input uses lazy indexed access by default. Inspect the
track with
[`summary_bwg()`](https://renscq.github.io/bwTools/reference/summary_bwg.md),
[`samples_bwg()`](https://renscq.github.io/bwTools/reference/samples_bwg.md),
[`seqinfo_bwg()`](https://renscq.github.io/bwTools/reference/seqinfo_bwg.md),
[`metadata_bwg()`](https://renscq.github.io/bwTools/reference/metadata_bwg.md),
and
[`zoominfo_bwg()`](https://renscq.github.io/bwTools/reference/zoominfo_bwg.md):

``` r

summary_bwg(bw)
samples_bwg(bw)
seqinfo_bwg(bw)
metadata_bwg(bw)
zoominfo_bwg(bw)
```

Use
[`is_bwg_track()`](https://renscq.github.io/bwTools/reference/is_bwg_track.md)
for a logical contract check and
[`validate_bwg()`](https://renscq.github.io/bwTools/reference/validate_bwg.md)
when an invalid object should raise an error:

``` r

is_bwg_track(bw)
validate_bwg(bw)
```

## Retrieve a region

``` r

region <- retrieve_bwg(
  bw,
  chrom = "chr1",
  start = 500L,
  end = 1600L
)
```

Request `result = "track"` when the result will be passed to another
bwTools operation:

``` r

region_track <- retrieve_bwg(
  bw,
  chrom = "chr1",
  start = 500L,
  end = 1600L,
  result = "track"
)
```

## Calculate statistics

``` r

stats_bwg(
  bw,
  chrom = "chr1",
  start = 1L,
  end = 10000L,
  n_bins = 20L,
  stat = "mean",
  use_zoom = FALSE
)
```

Use `use_zoom = TRUE` for fast large-region summaries when approximate
boundary handling from stored BigWig zoom windows is acceptable.

## Merge tracks

[`merge_bwg()`](https://renscq.github.io/bwTools/reference/merge_bwg.md)
conservatively preserves sample identity with `method = "auto"`:

``` r

a <- retrieve_bwg(bw, "chr1", 500L, 900L, result = "track")
b <- retrieve_bwg(bw, "chr1", 1200L, 1600L, result = "track")
merged <- merge_bwg(a, b)
```

Use `method = "mean"` or `method = "sum"` only when arithmetic
aggregation is intended.

## Write a result

``` r

manifest <- write_bwg(
  region_track,
  outdir = "bwtools-output",
  format = "bigwig",
  overwrite = TRUE
)
```

[`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md)
is the only normal signal-file persistence function. The returned
manifest has one row per selected written sample.

## Constructing a BwgTrack directly

Use
[`bwg_track()`](https://renscq.github.io/bwTools/reference/bwg_track.md)
when signal data already exist in R:

``` r

x <- bwg_track(
  samples = data.frame(sample_id = "sampleA"),
  data = data.frame(
    sample_id = "sampleA",
    chrom = "chr1",
    start = 1L,
    end = 10L,
    value = 2,
    strand = "+"
  ),
  seqinfo = data.frame(chrom = "chr1", length = 1000L)
)
```

The in-memory signal schema is `sample_id`, `chrom`, `start`, `end`,
`value`, and `strand`.
