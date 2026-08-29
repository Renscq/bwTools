# Performance diagnostics and compatibility validation

## Performance diagnostics

[`benchmark_bwg()`](https://renscq.github.io/bwTools/reference/benchmark_bwg.md)
is intended for development, profiling, and regression investigation. It
is not required for routine analysis.

``` r

library(bwTools)

bench <- benchmark_bwg(
  "/data/sample.bigwig",
  sample_id = "sampleA",
  region_sizes = c(1e4, 1e6, 1e7),
  n_bins = 1000L,
  iterations = 3L,
  warmup = 1L
)

bench$summary
```

The standard suite can measure lazy metadata loading, indexed retrieval,
zoom-aware statistics, and bounded exact statistics. Writer profiling
uses independent conservative repetition controls:

``` r

bench_write <- benchmark_bwg(
  "/data/sample.bigwig",
  sample_id = "sampleA",
  region_sizes = c(1e6, 5e6),
  operations = "write",
  write_iterations = 1L,
  write_warmup = 0L,
  verbose = TRUE
)

bench_write$summary
```

Heap metrics are approximate R heap measurements, not operating-system
RSS. The benchmark also does not clear the operating-system file cache.

## External BigWig compatibility

bwTools ships an optional compatibility harness under
`inst/compatibility/`. It can validate against independent
implementations when installed locally:

    pyBigWig / libBigWig
    rtracklayer
    UCSC BigWig utilities

These tools are not runtime dependencies.

Locate the installed runner:

``` r

runner <- system.file(
  "compatibility",
  "run-compatibility.001.R",
  package = "bwTools",
  mustWork = TRUE
)
runner
```

Run it from a shell:

``` bash
Rscript /path/to/run-compatibility.001.R \
  --output bwtools-compatibility
```

The generated `compatibility-report.tsv` uses `PASS`, `SKIP`, and
`FAIL`. Missing optional implementations are `SKIP` unless strict
validation is requested. A detected implementation that produces a real
signal mismatch is `FAIL`.

Compatibility compares signal semantics rather than requiring identical
BigWig binary block segmentation.
