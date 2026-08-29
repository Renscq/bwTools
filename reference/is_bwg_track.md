# Test whether an object follows the standardized BwgTrack contract

Test whether an object follows the standardized BwgTrack contract

## Usage

``` r
is_bwg_track(x)
```

## Arguments

- x:

  An R object.

## Value

Logical scalar.

## Details

Unlike
[`validate_bwg()`](https://renscq.github.io/bwTools/reference/validate_bwg.md),
this predicate returns `FALSE` for an invalid object rather than raising
an error.

## See also

[`validate_bwg()`](https://renscq.github.io/bwTools/reference/validate_bwg.md),
[`bwg_track()`](https://renscq.github.io/bwTools/reference/bwg_track.md)

## Examples

``` r
is_bwg_track(list())
#> [1] FALSE
```
