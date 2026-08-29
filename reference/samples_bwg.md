# Return sample metadata

Return sample metadata

## Usage

``` r
samples_bwg(x, sample_ids = NULL)
```

## Arguments

- x:

  A validated `BwgTrack` object.

- sample_ids:

  Optional sample IDs. Defaults to all samples.

## Value

A copy of the standardized sample metadata data.table.

## Details

The returned table is a copy and can be modified without changing `x`.

## See also

[`seqinfo_bwg()`](https://renscq.github.io/bwTools/reference/seqinfo_bwg.md),
[`metadata_bwg()`](https://renscq.github.io/bwTools/reference/metadata_bwg.md)

## Examples

``` r
x <- bwg_track(samples = data.frame(sample_id = "sampleA"), data = data.frame(
  sample_id = "sampleA", chrom = "chr1", start = 1L, end = 2L,
  value = 1, strand = "*"
))
samples_bwg(x)
#>    sample_id   file format strand has_strand
#>       <char> <char> <char> <char>     <lgcl>
#> 1:   sampleA   <NA> memory      *      FALSE
```
