# Read genomic signal tracks

Reads local genomic signal files into the standardized bwTools BwgTrack
contract.

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
  memory mode otherwise.

## Value

A validated `BwgTrack` object.
