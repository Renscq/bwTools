# Summarize a BwgTrack object

Summarize a BwgTrack object

## Usage

``` r
summary_bwg(x)
```

## Arguments

- x:

  A validated `BwgTrack` object.

## Value

A one-row data.table describing the object contract and dimensions.

## Details

This is a lightweight object summary and does not calculate genomic
signal statistics. Use
[`stats_bwg()`](https://renscq.github.io/bwTools/reference/stats_bwg.md)
for interval-level quantitative summaries.

## See also

[`stats_bwg()`](https://renscq.github.io/bwTools/reference/stats_bwg.md),
[`samples_bwg()`](https://renscq.github.io/bwTools/reference/samples_bwg.md)

## Examples

``` r
bw_file <- system.file("extdata", "bwtools_example.bigwig", package = "bwTools", mustWork = TRUE)
x <- read_bwg(bw_file, sample_ids = "example")
summary_bwg(x)
#>    samples   mode intervals chromosomes     coordinate schema_version
#>      <int> <char>     <int>       <int>         <char>         <char>
#> 1:       1   lazy        NA           3 1-based closed              2
#>             backend
#>              <char>
#> 1: bwTools-native-R
```
