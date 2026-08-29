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

## Details

Bins are constructed after known chromosome bounds are applied. Exact
lazy BigWig statistics stream indexed full-resolution blocks. Zoom-aware
results expose their source and stored reduction level in `source` and
`resolution`.

## See also

[`retrieve_bwg()`](https://renscq.github.io/bwTools/reference/retrieve_bwg.md),
[`zoominfo_bwg()`](https://renscq.github.io/bwTools/reference/zoominfo_bwg.md)

## Examples

``` r
bw_file <- system.file(
  "extdata", "bwtools_example.bigwig",
  package = "bwTools", mustWork = TRUE
)
x <- read_bwg(bw_file, sample_ids = "example")
stats_bwg(
  x, "chr1", 1L, 1000L,
  n_bins = 10L, stat = "mean", use_zoom = FALSE
)
#>     sample_id  chrom start   end   stat     value covered_bases source
#>        <char> <char> <int> <int> <char>     <num>         <num> <char>
#>  1:   example   chr1     1   100   mean        NA             0   full
#>  2:   example   chr1   101   200   mean        NA             0   full
#>  3:   example   chr1   201   300   mean        NA             0   full
#>  4:   example   chr1   301   400   mean        NA             0   full
#>  5:   example   chr1   401   500   mean  9.100000             1   full
#>  6:   example   chr1   501   600   mean 10.259859           100   full
#>  7:   example   chr1   601   700   mean  9.623916           100   full
#>  8:   example   chr1   701   800   mean  6.963250            20   full
#>  9:   example   chr1   801   900   mean 12.100000            21   full
#> 10:   example   chr1   901  1000   mean 13.622582           100   full
#>     resolution
#>          <num>
#>  1:         NA
#>  2:         NA
#>  3:         NA
#>  4:         NA
#>  5:         NA
#>  6:         NA
#>  7:         NA
#>  8:         NA
#>  9:         NA
#> 10:         NA
```
