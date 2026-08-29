# Benchmark a large BigWig workflow

Profiles the native bwTools BigWig workflow using deterministic genomic
regions. The standard suite measures package-metadata-cold and cached
lazy loading, indexed retrieval, zoom statistics, and bounded exact
statistics. Optional operations benchmark region materialization/merge,
native BigWig writing, or full-memory loading.

## Usage

``` r
benchmark_bwg(
  file,
  sample_id = "benchmark",
  regions = NULL,
  chrom = NULL,
  region_sizes = c(10000, 1e+06, 1e+07),
  n_bins = 1000L,
  stat = c("mean", "stdev", "max", "min", "coverage", "sum"),
  iterations = 3L,
  warmup = 1L,
  operations = c("read_lazy", "retrieve", "stats_zoom", "stats_full"),
  full_stats_max_bases = 1e+06,
  write_max_bases = 1e+07,
  write_max_intervals = 2e+06,
  write_zoom_max_intervals = 5e+05,
  write_iterations = 1L,
  write_warmup = 0L,
  verbose = TRUE
)
```

## Arguments

- file:

  One local BigWig file.

- sample_id:

  Sample identifier used for the benchmark track.

- regions:

  Optional explicit genomic regions with `chrom`, `start`, and `end`
  columns. An optional unique `region_id` column is preserved. When
  omitted, centered windows are generated from `region_sizes` on `chrom`
  or on the largest chromosome.

- chrom:

  Optional chromosome used for automatically generated regions.

- region_sizes:

  Region sizes in base pairs when `regions` is omitted.

- n_bins:

  Number of bins used by statistics benchmarks.

- stat:

  Statistic passed to
  [`stats_bwg()`](https://renscq.github.io/bwTools/reference/stats_bwg.md).

- iterations:

  Number of timed iterations per benchmark case.

- warmup:

  Number of untimed warmup iterations per benchmark case.

- operations:

  Operations to run. Supported values are `read_lazy`, `retrieve`,
  `stats_zoom`, `stats_full`, `merge`, `write`, and `read_memory`.

- full_stats_max_bases:

  Maximum region size for `stats_full`. Larger regions are recorded as
  skipped rather than forcing an expensive exact calculation.

- write_max_bases:

  Maximum region size for the optional `write` benchmark. Larger regions
  are recorded as skipped.

- write_max_intervals:

  Maximum number of signal intervals allowed for any writer benchmark
  case. Denser regions are skipped before writing.

- write_zoom_max_intervals:

  Maximum number of signal intervals allowed for `zoom_on` writer
  profiling. `zoom_off` may still run up to `write_max_intervals`. This
  protects the pure-R zoom builder from unexpectedly expensive
  dense-signal benchmarks.

- write_iterations:

  Timed iterations used by writer profiling. Writer benchmarks use an
  independent conservative default because each iteration creates a
  complete BigWig file.

- write_warmup:

  Untimed writer warmup iterations per writer case.

- verbose:

  Whether to print concise benchmark progress messages.

## Value

A `bwToolsBenchmark` object containing `system`, `file`, `config`,
`regions`, `results`, and `summary` components. Memory metrics
distinguish current live R heap after garbage collection from the
approximate maximum R heap observed since `gc(reset = TRUE)`; neither
metric is process RSS.

## Details

This diagnostic function is intentionally file-oriented because its
purpose is to measure file I/O. Normal analysis functions continue to
consume validated `BwgTrack` objects.

`benchmark_bwg()` is deliberately separated from normal analysis. Writer
cases use independent repetition and density guards so profiling does
not accidentally trigger excessive work on dense tracks. Reported heap
metrics are R-heap diagnostics rather than process resident-set size.

## See also

[`read_bwg()`](https://renscq.github.io/bwTools/reference/read_bwg.md),
[`retrieve_bwg()`](https://renscq.github.io/bwTools/reference/retrieve_bwg.md),
[`stats_bwg()`](https://renscq.github.io/bwTools/reference/stats_bwg.md),
[`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md)

## Examples

``` r
if (FALSE) { # \dontrun{
bw_file <- system.file(
  "extdata", "bwtools_example.bigwig",
  package = "bwTools", mustWork = TRUE
)
bench <- benchmark_bwg(
  bw_file,
  sample_id = "example",
  region_sizes = 1000L,
  iterations = 1L,
  warmup = 0L
)
bench$summary
} # }
```
