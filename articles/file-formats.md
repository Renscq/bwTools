# File formats and coordinate conventions

## Supported formats

[`detect_bwg_format()`](https://renscq.github.io/bwTools/reference/detect_bwg_format.md)
returns one of `bigwig`, `wig`, or `bedgraph`. BigWig is detected from
its four-byte binary magic signature before text parsing, so a valid
BigWig does not depend on its file extension.

``` r

library(bwTools)

bw_file <- system.file("extdata", "bwtools_example.bigwig", package = "bwTools")
bg_file <- system.file("extdata", "bwtools_example.bedGraph", package = "bwTools")

detect_bwg_format(bw_file)
detect_bwg_format(bg_file)
```

## Coordinate systems

bwTools exposes one public in-memory coordinate convention:

    1-based closed

A public interval `start = 1, end = 10` contains 10 bases.

On disk:

    BigWig   0-based half-open
    bedGraph 0-based half-open
    WIG      format-specific text representation

The conversion is performed only at the read/write boundary. Do not
manually subtract or add one when calling
[`retrieve_bwg()`](https://renscq.github.io/bwTools/reference/retrieve_bwg.md)
or
[`stats_bwg()`](https://renscq.github.io/bwTools/reference/stats_bwg.md).

## Chromosome names

Chromosome identifiers are preserved exactly. bwTools does not rewrite
`1` to `chr1`, normalize case, or map mitochondrial aliases. Therefore
`1`, `chr1`, `Chr01`, `MT`, and `scaffold_001` are distinct names.

## Sequence metadata

[`seqinfo_bwg()`](https://renscq.github.io/bwTools/reference/seqinfo_bwg.md)
returns:

    sample_id  chrom  length

BigWig stores chromosome lengths, so lazy BigWig tracks normally have
complete sequence metadata. WIG and bedGraph do not encode complete
chromosome sizes, so inferred lengths may be `NA` unless explicit
metadata are supplied through
[`bwg_track()`](https://renscq.github.io/bwTools/reference/bwg_track.md)
or `write_bwg(..., chrom_sizes = ...)`.

## BigWig zoom metadata

[`zoominfo_bwg()`](https://renscq.github.io/bwTools/reference/zoominfo_bwg.md)
returns stored reduction levels for BigWig-backed samples:

``` r

bw <- read_bwg(bw_file, sample_ids = "example")
zoominfo_bwg(bw)
```

WIG, bedGraph, and memory-only tracks have no stored BigWig zoom index
and return a typed zero-row zoom table.

## Format conversion

All conversion follows the same object path:

    source file -> read_bwg() -> BwgTrack -> write_bwg() -> target file

This keeps format-specific coordinate conversion out of retrieval,
statistics, and merging code.
