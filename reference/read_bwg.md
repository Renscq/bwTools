# Read genomic signal tracks

Reads one or more local BigWig, WIG, or bedGraph files into the
standardized bwTools `BwgTrack` contract.

## Usage

``` r
read_bwg(
  files,
  format = c("auto", "bigwig", "wig", "bedgraph"),
  sample_ids = NULL,
  strand = "*",
  mode = c("auto", "memory", "lazy")
)
```

## Arguments

- files:

  One or more local BigWig, WIG, or bedGraph files.

- format:

  Input format. `auto` detects each file independently.

- sample_ids:

  Optional sample IDs. Defaults to input file stems.

- strand:

  Strand labels associated with input files. Values must be `+`, `-`, or
  `*`.

- mode:

  Loading mode. `auto` uses lazy access when all inputs are BigWig and
  memory mode otherwise. `lazy` is supported for BigWig inputs only.

## Value

A validated `BwgTrack` object.

## Details

Public signal coordinates are converted to 1-based closed coordinates on
input. BigWig-only input is lazy in `mode = "auto"`; WIG, bedGraph, and
mixed-format input are materialized in memory. Sample IDs must be unique
within one call.

## See also

[`detect_bwg_format()`](https://renscq.github.io/bwTools/reference/detect_bwg_format.md),
[`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md),
[`bwg_track()`](https://renscq.github.io/bwTools/reference/bwg_track.md)

## Examples

``` r
bw_file <- system.file(
  "extdata", "bwtools_example.bigwig",
  package = "bwTools", mustWork = TRUE
)
x <- read_bwg(bw_file, sample_ids = "example")
summary_bwg(x)
#>    samples   mode intervals chromosomes     coordinate schema_version
#>      <int> <char>     <int>       <int>         <char>         <char>
#> 1:       1   lazy        NA           3 1-based closed              2
#>             backend
#>              <char>
#> 1: bwTools-native-R
```
