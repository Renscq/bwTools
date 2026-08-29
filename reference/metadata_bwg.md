# Return BwgTrack metadata

Return BwgTrack metadata

## Usage

``` r
metadata_bwg(x)
```

## Arguments

- x:

  A validated `BwgTrack` object.

## Value

The standardized metadata list.

## Details

Standard metadata include `coordinate`, `schema_version`, `backend`, and
`mode`. Treat returned metadata as descriptive object state rather than
a mechanism for mutating the track contract.

## See also

[`samples_bwg()`](https://renscq.github.io/bwTools/reference/samples_bwg.md),
[`seqinfo_bwg()`](https://renscq.github.io/bwTools/reference/seqinfo_bwg.md)

## Examples

``` r
x <- bwg_track(samples = data.frame(sample_id = "sampleA"), data = data.frame(
  sample_id = "sampleA", chrom = "chr1", start = 1L, end = 2L,
  value = 1, strand = "*"
))
metadata_bwg(x)
#> $coordinate
#> [1] "1-based closed"
#> 
#> $schema_version
#> [1] "2"
#> 
#> $backend
#> [1] "bwTools-native-R"
#> 
#> $mode
#> [1] "memory"
#> 
```
