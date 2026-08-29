# Detect genomic signal file format

BigWig detection is based on its four-byte binary signature before any
text decoding is attempted. This keeps format detection safe for binary
files and for installations whose paths are not represented internally
as UTF-8.

## Usage

``` r
detect_bwg_format(file)
```

## Arguments

- file:

  Local signal file path.

## Value

One of `bigwig`, `wig`, or `bedgraph`.

## Details

BigWig is detected from its binary magic signature before text decoding.
WIG and bedGraph are then identified from extension and content. Invalid
or empty files raise a `bwTools_error`.

## See also

[`read_bwg()`](https://renscq.github.io/bwTools/reference/read_bwg.md)

## Examples

``` r
bw_file <- system.file(
  "extdata", "bwtools_example.bigwig",
  package = "bwTools", mustWork = TRUE
)
detect_bwg_format(bw_file)
#> [1] "bigwig"
```
