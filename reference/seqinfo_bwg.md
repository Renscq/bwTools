# Return sample-specific chromosome metadata

Return sample-specific chromosome metadata

## Usage

``` r
seqinfo_bwg(x, sample_ids = NULL)
```

## Arguments

- x:

  A validated `BwgTrack` object.

- sample_ids:

  Optional sample IDs. Defaults to all samples.

## Value

A copy of the standardized `sample_id`, `chrom`, and `length` metadata
table.

## Details

BigWig-backed samples usually provide complete chromosome lengths. Text
formats may contain `NA` lengths when the source file does not encode
them. The returned table is a copy.

## See also

[`samples_bwg()`](https://renscq.github.io/bwTools/reference/samples_bwg.md),
[`metadata_bwg()`](https://renscq.github.io/bwTools/reference/metadata_bwg.md)

## Examples

``` r
bw_file <- system.file("extdata", "bwtools_example.bigwig", package = "bwTools", mustWork = TRUE)
x <- read_bwg(bw_file, sample_ids = "example")
seqinfo_bwg(x)
#>    sample_id  chrom length
#>       <char> <char>  <int>
#> 1:   example   chr1  10000
#> 2:   example   chr2   8000
#> 3:   example   chr3   5000
```
