# Construct a standardized BwgTrack object

Construct a standardized BwgTrack object

## Usage

``` r
bwg_track(samples, data = NULL, seqinfo = NULL, meta = list())
```

## Arguments

- samples:

  Sample metadata table containing `sample_id`. Missing standard columns
  are normalized automatically.

- data:

  Optional in-memory signal table using 1-based closed coordinates.

- seqinfo:

  Optional chromosome metadata table containing `chrom` and `length`.
  When `sample_id` is absent, rows are replicated across samples. When
  supplied, every signal sample/chromosome pair must be represented and
  signal intervals may not exceed a known chromosome length.

- meta:

  Optional metadata list.

## Value

A validated object inheriting from `BwgTrack` and `bwToolsTrack`.
