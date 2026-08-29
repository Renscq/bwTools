# Calculate interval statistics

Calculates coverage-aware statistics over equal-width genomic bins from
a standardized `BwgTrack`. Memory-mode tracks use full-resolution
signal. Lazy BigWig tracks may use stored zoom summaries when
`use_zoom = TRUE`. Exact lazy-BigWig calculations stream overlapping
full-resolution blocks directly into bin accumulators without
materializing an interval table. When all chromosome lengths for a
selected sample are known in `seqinfo`, query ends are clipped to the
known chromosome length and unknown chromosomes are rejected. Queries
starting beyond a known chromosome return a stable empty statistics
table.

## Usage

``` r
stats_bwg(
  x,
  chrom,
  start,
  end,
  n_bins = 1L,
  stat = c("mean", "stdev", "max", "min", "coverage", "sum"),
  sample_ids = NULL,
  use_zoom = TRUE
)
```

## Arguments

- x:

  A validated `BwgTrack` object.

- chrom:

  Chromosome name.

- start:

  Start coordinate, 1-based closed.

- end:

  End coordinate, 1-based closed.

- n_bins:

  Number of equal-width bins.

- stat:

  Statistic: `mean`, `stdev`, `max`, `min`, `coverage`, or `sum`.

- sample_ids:

  Optional sample IDs. Defaults to all samples.

- use_zoom:

  Whether lazy BigWig inputs may use zoom summary levels. Zoom results
  may be approximate when query bins partially overlap stored zoom
  windows. Use `FALSE` for exact values.

## Value

A data.table with one row per sample and genomic bin.
