# Validate a BwgTrack object

Performs a strict validation of the public bwTools `BwgTrack` contract.
Valid objects contain standardized `samples`, `data`, `seqinfo`, and
`meta` slots and use 1-based closed in-memory coordinates.

## Usage

``` r
validate_bwg(x)
```

## Arguments

- x:

  A `BwgTrack` object.

## Value

`x`, invisibly. An error is raised when the contract is invalid.

## Details

Validation checks standardized slot schemas, coordinate types, sample
and chromosome consistency, known chromosome bounds, and required
metadata. Package-generated validation failures inherit from
`bwTools_error`.

## See also

[`bwg_track()`](https://renscq.github.io/bwTools/reference/bwg_track.md),
[`is_bwg_track()`](https://renscq.github.io/bwTools/reference/is_bwg_track.md)

## Examples

``` r
x <- bwg_track(
  samples = data.frame(sample_id = "sampleA"),
  data = data.frame(
    sample_id = "sampleA", chrom = "chr1",
    start = 1L, end = 10L, value = 1, strand = "*"
  )
)
validate_bwg(x)
```
