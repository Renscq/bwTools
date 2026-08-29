# Retrieve signal intervals

Retrieves one or more genomic regions from a standardized `BwgTrack`.
Query coordinates are 1-based closed. Overlapping or directly adjacent
regions are merged before retrieval, and signal intervals crossing query
boundaries are clipped without rebasing genomic coordinates. When
complete chromosome lengths are available in `seqinfo`, query ends
beyond a chromosome are clipped to the known length and queries starting
beyond the chromosome return an empty result. Unknown chromosomes are
rejected only for selected samples whose `seqinfo` lengths are all
known; length-incomplete metadata remains permissive because it may
describe only chromosomes observed in the signal.

## Usage

``` r
retrieve_bwg(
  x,
  chrom = NULL,
  start = NULL,
  end = NULL,
  regions = NULL,
  sample_ids = NULL,
  strand = NULL,
  result = c("data", "track")
)
```

## Arguments

- x:

  A validated `BwgTrack` object.

- chrom:

  Optional chromosome name for a single-region query.

- start:

  Optional start coordinate for a single-region query.

- end:

  Optional end coordinate for a single-region query.

- regions:

  Optional data.frame or data.table containing `chrom`, `start`, and
  `end` columns. Use either `regions` or `chrom`/`start`/`end`.

- sample_ids:

  Optional sample IDs. Defaults to all samples.

- strand:

  Optional exact strand selector: `+`, `-`, or `*`. `NULL` performs no
  strand filtering. Memory-mode tracks filter row-level signal strand
  values; lazy BigWig tracks use sample-level strand metadata.

- result:

  Result type. `data` returns the standardized signal data.table;
  `track` returns a memory-mode `BwgTrack`.

## Value

A signal data.table or memory-mode `BwgTrack`.

## Details

Lazy BigWig retrieval uses the on-disk R-tree index. Memory-mode
retrieval subsets standardized signal rows. Empty queries return stable
typed outputs; `result = "track"` returns a memory-mode `BwgTrack`
suitable for downstream merge or write operations.

## See also

[`stats_bwg()`](https://renscq.github.io/bwTools/reference/stats_bwg.md),
[`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md)

## Examples

``` r
bw_file <- system.file(
  "extdata", "bwtools_example.bigwig",
  package = "bwTools", mustWork = TRUE
)
x <- read_bwg(bw_file, sample_ids = "example")
retrieve_bwg(x, "chr1", 500L, 1600L)
#>     sample_id  chrom start   end     value strand
#>        <char> <char> <int> <int>     <num> <char>
#>  1:   example   chr1   500   524  9.100000      *
#>  2:   example   chr1   525   549 10.044974      *
#>  3:   example   chr1   550   574 10.752699      *
#>  4:   example   chr1   575   599 11.069216      *
#>  5:   example   chr1   600   624 10.913620      *
#>  6:   example   chr1   625   649 10.295127      *
#>  7:   example   chr1   650   674  9.312052      *
#>  8:   example   chr1   675   699  8.132882      *
#>  9:   example   chr1   700   720  6.963250      *
#> 10:   example   chr1   880   904 12.100000      *
#> 11:   example   chr1   905   929 13.044974      *
#> 12:   example   chr1   930   954 13.752699      *
#> 13:   example   chr1   955   979 14.069216      *
#> 14:   example   chr1   980  1004 13.913620      *
#> 15:   example   chr1  1005  1029 13.295127      *
#> 16:   example   chr1  1030  1040 12.312052      *
#> 17:   example   chr1  1250  1274  8.100000      *
#> 18:   example   chr1  1275  1299  9.044974      *
#> 19:   example   chr1  1300  1324  9.752699      *
#> 20:   example   chr1  1325  1349 10.069216      *
#> 21:   example   chr1  1350  1374  9.913620      *
#> 22:   example   chr1  1375  1399  9.295127      *
#> 23:   example   chr1  1400  1424  8.312052      *
#> 24:   example   chr1  1425  1449  7.132882      *
#> 25:   example   chr1  1450  1474  5.963250      *
#> 26:   example   chr1  1475  1499  5.005459      *
#> 27:   example   chr1  1500  1510  4.418743      *
#>     sample_id  chrom start   end     value strand
#>        <char> <char> <int> <int>     <num> <char>
```
