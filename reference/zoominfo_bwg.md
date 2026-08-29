# Return BigWig zoom-level metadata

Returns stored zoom reduction levels for BigWig-backed samples in a
standardized `BwgTrack`. WIG, bedGraph, and memory-only samples do not
have BigWig zoom levels.

## Usage

``` r
zoominfo_bwg(x, sample_ids = NULL)
```

## Arguments

- x:

  A validated `BwgTrack` object.

- sample_ids:

  Optional sample IDs. Defaults to all samples.

## Value

A data.table containing `sample_id`, `level`, `data_offset`, and
`index_offset`.

## Details

Reduction levels describe stored BigWig zoom summaries. Samples without
a BigWig backing file return no rows rather than an error.

## See also

[`stats_bwg()`](https://renscq.github.io/bwTools/reference/stats_bwg.md),
[`metadata_bwg()`](https://renscq.github.io/bwTools/reference/metadata_bwg.md)

## Examples

``` r
bw_file <- system.file(
  "extdata", "bwtools_example.bigwig",
  package = "bwTools", mustWork = TRUE
)
x <- read_bwg(bw_file, sample_ids = "example")
zoominfo_bwg(x)
#> Empty data.table (0 rows and 4 cols): sample_id,level,data_offset,index_offset
```
