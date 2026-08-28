# bwTools

Current development release: **0.7.2**

`bwTools` provides native-R reading, writing, indexed retrieval, merging, and
statistics for BigWig, WIG, and bedGraph genomic signal tracks.

The 0.6.x series standardized the public API and the `BwgTrack` object contract
so downstream packages can use bwTools without depending on internal list
implementation details. Version 0.7.0 established reproducible large-BigWig
benchmarks; version 0.7.1 introduced block-streamed exact statistics, and version
0.7.2 refines indexed retrieval and benchmark memory diagnostics without changing
the public API.

## Core principles

- BigWig binary I/O is implemented in R.
- No UCSC `bedGraphToBigWig` executable is required.
- No Rcpp or bundled C/C++ runtime is required.
- Public in-memory coordinates are always **1-based closed**.
- BigWig and bedGraph disk coordinates remain **0-based half-open**.
- BigWig files use lazy indexed access by default.
- BigWig format detection checks the binary magic signature before text decoding.
- All normal analysis functions operate on a validated `BwgTrack` object.
- Signal ingestion is centralized in `read_bwg()`; `benchmark_bwg()` is a diagnostic file-I/O exception.
- Only `write_bwg()` writes signal files.
- Performance changes are benchmark-driven rather than assumption-driven.

## Standard public API

| Function | Responsibility | Primary return |
| --- | --- | --- |
| `detect_bwg_format()` | Detect file format | character |
| `read_bwg()` | Read files into bwTools | `BwgTrack` |
| `bwg_track()` | Construct an in-memory/lazy track | `BwgTrack` |
| `validate_bwg()` | Validate the public object contract | invisible `BwgTrack` |
| `is_bwg_track()` | Test the public object contract | logical |
| `samples_bwg()` | Access sample metadata | data.table |
| `seqinfo_bwg()` | Access chromosome metadata | data.table |
| `zoominfo_bwg()` | Access BigWig zoom metadata | data.table |
| `retrieve_bwg()` | Retrieve genomic regions | data.table / `BwgTrack` |
| `merge_bwg()` | Merge or aggregate tracks | `BwgTrack` |
| `metadata_bwg()` | Access object metadata | list |
| `stats_bwg()` | Calculate interval statistics | data.table |
| `summary_bwg()` | Summarize a track | one-row data.table |
| `write_bwg()` | Persist signal tracks | invisible manifest data.table |
| `benchmark_bwg()` | Benchmark large BigWig workflows | `bwToolsBenchmark` |

The old public names `bw_track()` and `bw_detect_format()` are removed in
0.6.1. Use `bwg_track()` and `detect_bwg_format()`.

## Standard argument names

The public API follows these naming rules:

```{text}
files       file paths passed to read_bwg()
file        one BigWig path passed to benchmark_bwg()
x           one BwgTrack object
sample_ids  sample identifiers or sample selection
samples     sample metadata table used by bwg_track()
regions     genomic interval table
chrom       chromosome for a single-region operation
start       1-based closed start
end         1-based closed end
outdir      output directory used only by write_bwg()
```

`sample_names` and selector arguments named `samples` are no longer used.

## BwgTrack contract

A valid object contains exactly four public components:

```{text}
BwgTrack
├── samples
├── data
├── seqinfo
└── meta
```

### samples

Required core columns:

```{text}
sample_id  file  format  strand  has_strand
```

- `sample_id` is unique and non-empty.
- `strand` is `+`, `-`, or `*`.
- `has_strand` is derived from `strand`.
- `format` uses `bigwig`, `wig`, `bedgraph`, or `memory`.

### data

Memory-mode signal uses:

```{text}
sample_id  chrom  start  end  value  strand
```

Coordinates are 1-based closed. Values must be finite numeric values.
The standardized schema stores `start` and `end` as integer columns.

For lazy BigWig tracks, `data` is `NULL`.

### seqinfo

Standard chromosome metadata uses:

```{text}
sample_id  chrom  length
```

Unknown text-format chromosome lengths are represented by `NA`.

### meta

Every object contains:

```{text}
coordinate     = "1-based closed"
schema_version = "2"
backend        = "bwTools-native-R"
mode           = "memory" or "lazy"
```

Downstream packages should use the public accessors whenever possible instead
of relying directly on internal list fields. Use `metadata_bwg()` for operation,
provenance, and schema metadata.

## Bundled example data

```{r}
library(bwTools)

bw_file <- system.file(
  "extdata",
  "bwtools_example.bigwig",
  package = "bwTools",
  mustWork = TRUE
)

bedgraph_file <- system.file(
  "extdata",
  "bwtools_example.bedGraph",
  package = "bwTools",
  mustWork = TRUE
)

chrom_sizes_file <- system.file(
  "extdata",
  "bwtools_example.chrom.sizes",
  package = "bwTools",
  mustWork = TRUE
)
```

## 1. Detect a file format

```{r}
stopifnot(
  identical(detect_bwg_format(bw_file), "bigwig"),
  identical(detect_bwg_format(bedgraph_file), "bedgraph")
)
```

## 2. Read tracks

`mode = "auto"` is the default.

- all-BigWig input -> lazy mode;
- WIG or bedGraph input -> memory mode.

```{r}
bw <- read_bwg(
  bw_file,
  sample_ids = "example"
)

bedgraph <- read_bwg(
  bedgraph_file,
  sample_ids = "example"
)

stopifnot(
  identical(summary_bwg(bw)$mode, "lazy"),
  identical(summary_bwg(bedgraph)$mode, "memory")
)
```

Use explicit loading only when needed:

```{r}
bw_memory <- read_bwg(
  bw_file,
  sample_ids = "example",
  mode = "memory"
)
```

For large RNA-seq and Ribo-seq BigWig tracks, lazy mode is recommended.

## 3. Validate the object contract

```{r}
validate_bwg(bw)
stopifnot(is_bwg_track(bw))
```

Inspect standardized metadata:

```{r}
samples_bwg(bw)
seqinfo_bwg(bw)
summary_bwg(bw)
```

Expected chromosome lengths for the bundled BigWig are:

```{r}
seqinfo <- seqinfo_bwg(bw)

stopifnot(
  identical(seqinfo$chrom, c("chr1", "chr2", "chr3")),
  identical(seqinfo$length, c(10000L, 8000L, 5000L))
)
```

## 4. Retrieve genomic signal

Single region:

```{r}
signal <- retrieve_bwg(
  bw,
  chrom = "chr1",
  start = 450L,
  end = 1600L
)

stopifnot(
  nrow(signal) > 0L,
  all(signal$chrom == "chr1"),
  all(signal$start >= 450L),
  all(signal$end <= 1600L)
)
```

Multiple regions:

```{r}
regions <- data.frame(
  chrom = c("chr1", "chr2"),
  start = c(500L, 450L),
  end = c(900L, 800L)
)

signal_multi <- retrieve_bwg(
  bw,
  regions = regions
)
```

Return a write-ready `BwgTrack` when the retrieved signal needs further package
operations:

```{r}
region_track <- retrieve_bwg(
  bw,
  regions = regions,
  result = "track"
)

validate_bwg(region_track)
```

Sample selection always uses `sample_ids`:

```{r}
signal_one_sample <- retrieve_bwg(
  bw,
  chrom = "chr1",
  start = 500L,
  end = 900L,
  sample_ids = "example"
)
```

## 5. Merge tracks

The default merge is conservative and requires no prior classification of the
input samples or genomic intervals. The following two tracks contain the same
sample on different regions:

```{r}
track_a <- bwg_track(
  samples = data.frame(sample_id = "sampleA"),
  data = data.frame(
    sample_id = "sampleA",
    chrom = "chr1",
    start = 1L,
    end = 100L,
    value = 2,
    strand = "*"
  )
)

track_b <- bwg_track(
  samples = data.frame(sample_id = "sampleA"),
  data = data.frame(
    sample_id = "sampleA",
    chrom = "chr1",
    start = 201L,
    end = 300L,
    value = 4,
    strand = "*"
  )
)

merged <- merge_bwg(track_a, track_b)
stopifnot(nrow(merged$data) == 2L)
```

Default `method = "auto"` behavior:

```{text}
different sample IDs
  -> preserve samples independently

same sample, non-overlapping regions
  -> append directly

same sample, identical records
  -> deduplicate

same sample, overlapping equal values
  -> merge safely

same sample, overlapping different values
  -> error and require explicit mean or sum
```

Explicit arithmetic aggregation uses overlapping tracks. Here the sample ID
is the same and the overlapping values differ, so `method = "auto"` would stop
and require an explicit arithmetic rule:

```{r}
overlap_a <- bwg_track(
  samples = data.frame(sample_id = "sampleA"),
  data = data.frame(
    sample_id = "sampleA", chrom = "chr1",
    start = 1L, end = 100L, value = 2, strand = "*"
  )
)

overlap_b <- bwg_track(
  samples = data.frame(sample_id = "sampleA"),
  data = data.frame(
    sample_id = "sampleA", chrom = "chr1",
    start = 51L, end = 150L, value = 4, strand = "*"
  )
)

merged_mean <- merge_bwg(
  overlap_a,
  overlap_b,
  method = "mean"
)

merged_sum <- merge_bwg(
  overlap_a,
  overlap_b,
  method = "sum"
)
```

For means, the missing-signal rule is explicit:

```{r}
mean_ignore <- merge_bwg(
  overlap_a,
  overlap_b,
  method = "mean",
  missing = "ignore"
)

mean_zero <- merge_bwg(
  overlap_a,
  overlap_b,
  method = "mean",
  missing = "zero"
)
```

## 6. Calculate statistics

All statistics operate on `BwgTrack`; pass files through `read_bwg()` first.

Exact full-resolution calculation:

```{r}
exact_stats <- stats_bwg(
  bw,
  chrom = "chr1",
  start = 1L,
  end = 10000L,
  n_bins = 10L,
  stat = "mean",
  use_zoom = FALSE
)
```

Fast lazy BigWig calculation:

```{r}
fast_stats <- stats_bwg(
  bw,
  chrom = "chr1",
  start = 1L,
  end = 10000L,
  n_bins = 10L,
  stat = "mean",
  use_zoom = TRUE
)
```

Supported statistics:

```{text}
mean
stdev
min
max
coverage
sum
```

For final quantitative analyses, use `use_zoom = FALSE` when exact values are
required.

## 7. Inspect BigWig zoom levels

Zoom metadata are accessed from the track, not directly from a path:

```{r}
zoom_info <- zoominfo_bwg(bw)
print(zoom_info)
```

Expected columns:

```{text}
sample_id  level  data_offset  index_offset
```

## 8. Write tracks

All persistence and format conversion are handled by `write_bwg()`.

```{r}
test_root <- file.path(tempdir(), "bwtools_api_test")
dir.create(test_root, recursive = TRUE, showWarnings = FALSE)

written <- write_bwg(
  bedgraph,
  outdir = test_root,
  format = "bigwig",
  chrom_sizes = chrom_sizes_file,
  overwrite = TRUE
)

stopifnot(file.exists(written$file))
```

Write one selected sample:

```{r}
write_bwg(
  bedgraph,
  outdir = file.path(test_root, "selected"),
  format = "bigwig",
  sample_ids = "example",
  chrom_sizes = chrom_sizes_file,
  overwrite = TRUE
)
```

BigWig writer options remain writer-specific:

```{r}
write_bwg(
  bedgraph,
  outdir = file.path(test_root, "zoom"),
  format = "bigwig",
  chrom_sizes = chrom_sizes_file,
  overwrite = TRUE,
  zoom = TRUE,
  max_zoom_levels = 10L
)
```

They are not duplicated in retrieval or merge functions.

## 9. Construct a BwgTrack manually

External packages can construct the same contract through `bwg_track()`:

```{r}
x <- bwg_track(
  samples = data.frame(
    sample_id = "sampleA",
    strand = "+"
  ),
  data = data.frame(
    sample_id = "sampleA",
    chrom = "chr1",
    start = 101L,
    end = 125L,
    value = 8.5,
    strand = "+"
  ),
  seqinfo = data.frame(
    chrom = "chr1",
    length = 1000000L
  )
)

validate_bwg(x)
```

`bwg_track()` normalizes missing standard fields. For example, the sample table
above receives `file`, `format`, and `has_strand` automatically, while seqinfo
without `sample_id` is expanded to the object's samples.

## 10. Recommended downstream-package integration

Downstream packages should use this pattern:

```{r}
x <- read_bwg(files, sample_ids = sample_ids)
validate_bwg(x)

sample_meta <- samples_bwg(x)
chrom_meta <- seqinfo_bwg(x)
track_meta <- metadata_bwg(x)

signal <- retrieve_bwg(
  x,
  chrom = chrom,
  start = start,
  end = end,
  sample_ids = sample_ids
)
```

Do not depend on private functions such as:

```{text}
bwTools:::bw_metadata()
bwTools:::bw_bigwig_query()
bwTools:::bw_read_index_header()
```

The public contract is the integration boundary.


## 11. Benchmark a large BigWig file

`benchmark_bwg()` is a diagnostic exception to the normal object-only analysis
contract because its purpose is to measure file I/O itself. It accepts one local
BigWig file and returns a structured `bwToolsBenchmark` object; it does not write
benchmark reports automatically.

The default suite is deliberately safe for large files:

```{text}
read_lazy   package-metadata-cold and cached lazy loading
retrieve    indexed retrieval at multiple genomic scales
stats_zoom  zoom-enabled statistics
stats_full  exact statistics only below full_stats_max_bases
```

The default centered benchmark windows are 10 kb, 1 Mb, and 10 Mb on the
largest chromosome. For real performance work, explicit representative regions
are preferred because signal density can vary substantially across a genome.

```{r}
bench <- benchmark_bwg(
  "/data/sample.bigwig",
  sample_id = "sampleA",
  region_sizes = c(1e4, 1e6, 1e7),
  n_bins = 1000L,
  iterations = 3L,
  warmup = 1L
)

bench
bench$summary
```

For a representative dense and sparse region set:

```{r}
benchmark_regions <- data.frame(
  chrom = c("chr1", "chr1", "chr2"),
  start = c(1000001L, 50000001L, 20000001L),
  end = c(1010000L, 51000000L, 30000000L)
)

bench <- benchmark_bwg(
  "/data/sample.bigwig",
  regions = benchmark_regions,
  iterations = 3L
)
```

Optional operations are explicit:

```{r}
bench_io <- benchmark_bwg(
  "/data/sample.bigwig",
  operations = c(
    "read_lazy",
    "retrieve",
    "stats_zoom",
    "stats_full",
    "merge",
    "write"
  )
)
```

`merge` materializes four non-overlapping parts of the largest selected region
before timing `merge_bwg()`. `write` writes only the largest selected region and
benchmarks both `zoom = FALSE` and the default `zoom = TRUE`. This keeps the
benchmark bounded by the chosen genomic regions.

Full-file memory loading is never part of the default suite. It must be selected
explicitly:

```{r}
bench_memory <- benchmark_bwg(
  "/data/sample.bigwig",
  operations = "read_memory",
  iterations = 1L,
  warmup = 0L
)
```

For multi-gigabyte files, `read_memory` may require much more RAM than the
compressed BigWig file size and can terminate the R process if system memory is
insufficient.

### 0.7.2 benchmark-driven retrieval refinement

Version 0.7.2 keeps the block-streamed exact-statistics path introduced in 0.7.1
and changes only the indexed-retrieval materialization strategy and benchmark
memory reporting.

The 0.7.1 real-file benchmark showed that streaming exact statistics substantially
reduced elapsed time, while the growable-vector retrieval path did not improve
large indexed queries consistently. Version 0.7.2 therefore retains the
successful statistics path and replaces repeated vector growth with chunked
primitive vectors:

```{text}
stats_bwg(use_zoom = FALSE)
  BigWig R-tree -> one compressed block -> bin accumulators -> next block

retrieve_bwg()
  BigWig blocks -> primitive vector chunks -> one-shot flatten -> one data.table
```

Each decoded retrieval block contributes only `start`, `end`, and `value`
vectors. These chunks are flattened once after all relevant blocks are decoded.
This avoids both per-block data.table construction and repeated resizing/copying
of large R vectors.

The benchmark memory contract is also more explicit in 0.7.2. Results now
separate current live heap after garbage collection from the maximum heap seen
since the benchmark reset:

```{text}
baseline_heap_mb      live R heap immediately before the timed expression
live_heap_after_mb    live R heap after the expression and garbage collection
peak_heap_mb          maximum R heap observed since gc(reset = TRUE)
peak_heap_delta_mb    peak_heap_mb minus baseline_heap_mb
live_heap_delta_mb    live_heap_after_mb minus baseline_heap_mb
heap_delta_mb         compatibility alias of peak_heap_delta_mb
```

These remain R-heap diagnostics, not operating-system RSS measurements.
`heap_delta_mb` is retained so benchmark-processing scripts written for 0.7.0 or
0.7.1 continue to work.

To compare 0.7.1 and 0.7.2, rerun the same benchmark configuration and inspect:

```{text}
retrieve    median_elapsed_sec / median_peak_heap_delta_mb / median_live_heap_delta_mb
stats_full  should retain the 0.7.1 streaming performance
stats_zoom  should not regress
read_lazy   should not regress
```

Performance is machine-, file-, and signal-density-dependent, so bwTools does
not enforce fixed timing thresholds in package tests. Regression tests require
retrieved signals and exact statistics to remain identical to reference paths.

A larger exact-statistics stress test can still be enabled explicitly by raising
the safety bound:

```{r}
bench_full_large <- benchmark_bwg(
  "/data/sample.bigwig",
  sample_id = "sampleA",
  region_sizes = c(1e6, 1e7),
  n_bins = 1000L,
  iterations = 3L,
  warmup = 1L,
  operations = "stats_full",
  full_stats_max_bases = 1e7
)
```

Keep the default bound for routine diagnostics; raise it only when intentionally
stress-testing the exact full-resolution path.

### Benchmark result contract

The returned object contains:

```{text}
bwToolsBenchmark
├── system   R, platform, OS, bwTools, and data.table versions
├── file     source file size and genome dimensions
├── config   benchmark settings and measurement notes
├── regions  actual benchmark regions
├── results  one row per timed iteration
└── summary  median/min/max summaries by operation and region
```

Important `results` metrics include:

```{text}
elapsed_sec       wall-clock elapsed time
user_sec          R user CPU time
system_sec        R system CPU time
baseline_heap_mb      live R heap immediately before the timed expression
live_heap_after_mb    live R heap after the expression and garbage collection
peak_heap_mb          approximate maximum R heap since the benchmark reset
peak_heap_delta_mb    approximate peak heap increase above baseline
live_heap_delta_mb    retained live heap change relative to baseline
heap_delta_mb         compatibility alias of peak_heap_delta_mb
result_mb         returned R object size
result_rows       returned interval/bin count
mbases_per_sec    genomic-region throughput
file_mb_per_sec   full-memory file-read throughput
output_mb_per_sec native write-output throughput
stat_source       actual statistics source: zoom or full
resolution        actual zoom reduction level when used
status            ok / error / skipped
message           captured diagnostic message
```

The heap metrics are derived from `gc()` and are **not** operating-system RSS.
`live_heap_after_mb` is measured after garbage collection while the benchmark
result remains referenced; `peak_heap_mb` is the maximum R heap observed since
`gc(reset = TRUE)`. Likewise, `metadata_cache_cold` clears only the bwTools
metadata cache; it does not flush the operating-system file cache. These limitations are recorded
in `bench$config` so benchmark reports are not over-interpreted.

To save benchmark tables, use normal table-writing tools explicitly:

```{r}
data.table::fwrite(
  bench$results,
  "bwtools_benchmark.results.tsv",
  sep = "\t"
)

data.table::fwrite(
  bench$summary,
  "bwtools_benchmark.summary.tsv",
  sep = "\t"
)
```

## Return-value contract

| Function | Return |
| --- | --- |
| `read_bwg()` | `BwgTrack` |
| `bwg_track()` | `BwgTrack` |
| `validate_bwg()` | input object invisibly |
| `samples_bwg()` | data.table |
| `seqinfo_bwg()` | data.table |
| `metadata_bwg()` | list |
| `retrieve_bwg(..., result = "data")` | data.table |
| `retrieve_bwg(..., result = "track")` | `BwgTrack` |
| `merge_bwg()` | `BwgTrack` |
| `stats_bwg()` | data.table |
| `zoominfo_bwg()` | data.table |
| `summary_bwg()` | one-row data.table |
| `write_bwg()` | invisible manifest data.table |
| `benchmark_bwg()` | `bwToolsBenchmark` |

## Package regression tests

```{r}
devtools::document()
devtools::test()
devtools::check()
```

Expected development target:

```{text}
FAIL 0
WARN 0
SKIP 0
```

## Development status

- **0.1.x**: native reader and initial BwgTrack contract — implemented.
- **0.2.x**: native BigWig/WIG/bedGraph writer — implemented.
- **0.3.x**: BigWig zoom levels and interval statistics — implemented.
- **0.4.x**: unified indexed genomic retrieval — implemented.
- **0.5.x**: automatic merge and explicit signal aggregation — implemented.
- **0.6.1**: standardized public API and BwgTrack schema v2.
- **0.6.2**: hardened BigWig detection and path-encoding safety.
- **0.6.3**: fixed automatic BigWig dispatch after format detection.
- **0.7.0**: reproducible large-BigWig benchmark suite.
- **0.7.1**: block-streamed exact-statistics optimization.
- **0.7.2**: chunked indexed retrieval and clearer heap diagnostics — current stage.
- **0.8.x**: statistics and zoom refinements where justified by benchmarks.
- **0.9.x**: API freeze, cross-platform checks, and release hardening.
- **1.0.0**: stable public release.
