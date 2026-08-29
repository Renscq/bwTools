# Return BigWig zoom-level metadata

Return stored BigWig zoom levels for BigWig-backed samples.

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
