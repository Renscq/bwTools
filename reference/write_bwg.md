# Write genomic signal tracks

Write a standardized BwgTrack to BigWig, WIG, or bedGraph files.

## Usage

``` r
write_bwg(
  x,
  outdir,
  format = c("bigwig", "wig", "bedgraph"),
  sample_ids = NULL,
  chrom_sizes = NULL,
  overwrite = FALSE,
  compress = FALSE,
  zoom = TRUE,
  max_zoom_levels = 10L
)
```

## Arguments

- x:

  A validated `BwgTrack` object.

- outdir:

  Output directory.

- format:

  Output format: `bigwig`, `wig`, or `bedgraph`.

- sample_ids:

  Optional sample IDs. Defaults to all samples.

- chrom_sizes:

  Chromosome sizes for BigWig output.

- overwrite:

  Whether existing output files may be replaced.

- compress:

  Whether WIG or bedGraph output should be gzip-compressed.

- zoom:

  Whether BigWig output should include zoom summary levels.

- max_zoom_levels:

  Maximum number of BigWig zoom levels.

## Value

Invisibly returns a manifest data.table with `sample_id`, `file`, and
`format`. Every selected in-memory sample must contain signal;
validation occurs before any output file is written.
