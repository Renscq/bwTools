# Write genomic signal tracks

Write genomic signal tracks

Writes a standardized `BwgTrack` to BigWig, WIG, or bedGraph files. File
persistence is handled only by this function; retrieval, merging, and
statistics never write files.

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

  Chromosome sizes for BigWig output. May be a two-column file or table.
  If omitted, complete chromosome lengths are taken from
  `seqinfo_bwg(x)` when available.

- overwrite:

  Whether existing output files may be replaced.

- compress:

  Whether WIG or bedGraph output should be gzip-compressed. BigWig data
  blocks are always compressed internally by the file format.

- zoom:

  Whether native BigWig output should include zoom summary levels.

- max_zoom_levels:

  Maximum number of BigWig zoom levels to create when `zoom = TRUE`.

## Value

Invisibly returns a data.table with `sample_id`, `file`, and `format`.
Every selected in-memory sample must contain signal; validation occurs
before any output file is written.

## Details

BigWig output uses 0-based half-open coordinates on disk while
preserving the 1-based closed public coordinate contract in memory.
BigWig requires chromosome lengths and writes native zoom summaries by
default. WIG and bedGraph can optionally be gzip-compressed.

## See also

[`read_bwg()`](https://renscq.github.io/bwTools/reference/read_bwg.md),
[`retrieve_bwg()`](https://renscq.github.io/bwTools/reference/retrieve_bwg.md)

## Examples

``` r
if (FALSE) { # \dontrun{
bw_file <- system.file(
  "extdata", "bwtools_example.bigwig",
  package = "bwTools", mustWork = TRUE
)
x <- read_bwg(bw_file, sample_ids = "example", mode = "memory")
write_bwg(
  x,
  outdir = tempfile("bwtools-write-"),
  format = "bigwig"
)
} # }
```
