---
title: "bwTools"
format: gfm
execute:
  echo: true
  eval: false
---

# bwTools

Current development release: **0.8.4**

`bwTools` provides native-R reading, writing, indexed retrieval, merging,
statistics, and performance diagnostics for BigWig, WIG, and bedGraph genomic
signal tracks.

Version 0.8.0 is the **API-stability baseline** before downstream integration.
The optimized BigWig reader, writer, zoom pyramid, indexed retrieval, and
exact/zoom statistics established in 0.7.x are retained. Version 0.8.2 hardens
chromosome boundaries, empty query behavior, and `data`/`seqinfo` consistency.
Version 0.8.3 aligns scalar retrieval interval validation with `stats_bwg()`.
Version 0.8.4 adds reproducible external BigWig compatibility validation without
changing the 15-function public API, BwgTrack schema v2, or runtime dependency
policy.

## Design principles

- BigWig binary I/O is implemented in R.
- No UCSC `bedGraphToBigWig` executable is required.
- No Rcpp or bundled C/C++ runtime is required.
- Public in-memory coordinates are always **1-based closed**.
- BigWig and bedGraph disk coordinates are **0-based half-open**.
- BigWig inputs use lazy indexed access by default.
- Normal analysis functions operate on a validated `BwgTrack`.
- `read_bwg()` is the normal file-ingestion entry point.
- `write_bwg()` is the only normal signal-file persistence function.
- `benchmark_bwg()` is an advanced diagnostic utility, not part of routine
  analysis.

## Public API

| Function | Responsibility | Return |
| --- | --- | --- |
| `detect_bwg_format()` | Detect BigWig/WIG/bedGraph | character |
| `read_bwg()` | Read signal files | `BwgTrack` |
| `bwg_track()` | Construct a standardized track | `BwgTrack` |
| `validate_bwg()` | Validate the object contract | invisible `BwgTrack` |
| `is_bwg_track()` | Test the object contract | logical |
| `samples_bwg()` | Access sample metadata | data.table |
| `seqinfo_bwg()` | Access chromosome metadata | data.table |
| `metadata_bwg()` | Access object metadata | list |
| `zoominfo_bwg()` | Access BigWig zoom metadata | data.table |
| `retrieve_bwg()` | Retrieve genomic regions | data.table / `BwgTrack` |
| `merge_bwg()` | Merge or aggregate tracks | `BwgTrack` |
| `stats_bwg()` | Calculate interval statistics | data.table |
| `summary_bwg()` | Summarize a track | one-row data.table |
| `write_bwg()` | Write BigWig/WIG/bedGraph | invisible manifest |
| `benchmark_bwg()` | Profile large BigWig workflows | `bwToolsBenchmark` |

The parameter names and core return schemas above are regression-tested in
0.8.0 so downstream packages can treat them as the integration boundary.

## Installation

Install the source package with your normal R package workflow. During package
development:

```{r}
devtools::document()
devtools::test()
devtools::check()
```

## Bundled example data

The package contains a small reference dataset:

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

## 1. Detect file format

```{r}
detect_bwg_format(bw_file)
detect_bwg_format(bedgraph_file)
```

Expected formats:

```{text}
bigwig
bedgraph
```

BigWig detection checks the binary magic signature before attempting any text
parsing, so file extensions are not required for valid BigWig files.

## 2. Read signal tracks

### BigWig

```{r}
bw <- read_bwg(
  bw_file,
  sample_ids = "example"
)

bw
summary_bwg(bw)
```

`mode = "auto"` is the default. A BigWig-only input becomes a lazy track:

```{text}
file
  ↓
read_bwg()
  ↓
lazy BwgTrack
  ↓
indexed region access
```

You can explicitly load a small BigWig into memory:

```{r}
bw_memory <- read_bwg(
  bw_file,
  sample_ids = "example",
  mode = "memory"
)
```

### bedGraph and WIG

Text formats are read into memory:

```{r}
bg <- read_bwg(
  bedgraph_file,
  sample_ids = "example"
)

summary_bwg(bg)
```

Mixed BigWig/text inputs use memory mode because lazy mode is currently a
BigWig-only contract.

## 3. BwgTrack object contract

A standardized track contains four components:

```{text}
BwgTrack
├── samples
├── data
├── seqinfo
└── meta
```

Core `samples` columns:

```{text}
sample_id  file  format  strand  has_strand
```

Core memory-mode `data` columns:

```{text}
sample_id  chrom  start  end  value  strand
```

Core `seqinfo` columns:

```{text}
sample_id  chrom  length
```

Required metadata:

```{text}
coordinate     = "1-based closed"
schema_version = "2"
backend        = "bwTools-native-R"
mode           = "memory" or "lazy"
```

Use the public accessors instead of depending on private implementation details:

```{r}
samples_bwg(bw)
seqinfo_bwg(bw)
metadata_bwg(bw)
zoominfo_bwg(bw)
```

Validate an externally constructed object before passing it to downstream code:

```{r}
validate_bwg(bw)
is_bwg_track(bw)
```

## 4. Retrieve genomic regions

For a single region:

```{r}
signal <- retrieve_bwg(
  bw,
  chrom = "chr1",
  start = 500L,
  end = 1600L
)

signal
```

For multiple regions:

```{r}
regions <- data.frame(
  chrom = c("chr1", "chr1", "chr2"),
  start = c(510L, 561L, 450L),
  end = c(560L, 600L, 900L)
)

signal_multi <- retrieve_bwg(
  bw,
  regions = regions
)
```

Overlapping or directly adjacent query regions are merged before retrieval.
Signal intervals crossing query boundaries are clipped, but genomic coordinates
are never rebased.

If the result will be used by another bwTools function, request a `BwgTrack`:

```{r}
region_track <- retrieve_bwg(
  bw,
  regions = regions,
  result = "track"
)
```

This cleanly separates retrieval from persistence:

```{text}
retrieve_bwg()
      ↓
BwgTrack
      ↓
write_bwg()
```

### Query boundary behavior

Query coordinates must be positive 1-based closed integers. When selected
samples have complete known chromosome lengths in `seqinfo`, bwTools uses the
following rules consistently in memory and lazy modes:

```{text}
unknown chromosome             -> error
end beyond chromosome length   -> clip to chromosome length
start beyond chromosome length -> empty typed result
```

For example, a query of `chr1:950-1050` on a chromosome of length 1000 is
interpreted as `chr1:950-1000`. The effective bounded regions are stored in
`metadata_bwg(retrieve_bwg(..., result = "track"))$retrieval_regions`.

A sample is treated as having complete sequence metadata only when all of its
`seqinfo$length` values are known. WIG and bedGraph do not encode chromosome
lengths, so their automatically inferred `seqinfo` contains `length = NA`; a
chromosome not observed in such incomplete metadata is treated as having no
signal rather than automatically being rejected.

Chromosome identifiers are never normalized by bwTools. Names such as `1`,
`chr1`, `Chr01`, `I`, `MT`, and `scaffold_001` remain distinct.

## 5. Interval statistics

Supported statistics:

```{text}
mean
stdev
max
min
coverage
sum
```

Exact full-resolution statistics:

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

Fast zoom-aware statistics:

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

`use_zoom = TRUE` means “allow an appropriate BigWig zoom level”. Small
high-resolution queries may automatically fall back to full-resolution data.
Inspect the output columns:

```{text}
source      full or zoom
resolution  zoom reduction level, or NA for full
```

For final quantitative calculations where exact boundary behavior matters, use
`use_zoom = FALSE`.

When a known chromosome end truncates a query, clipping occurs before `n_bins`
are constructed. Queries beginning entirely beyond a known chromosome return a
typed zero-row statistics table with the normal output columns.

## 6. Merge signal tracks

`merge_bwg()` is the only public merge interface.

### Automatic lossless merge

```{r}
merged <- merge_bwg(track_a, track_b)
```

Default `method = "auto"` follows conservative rules:

- different sample IDs remain different samples;
- non-overlapping segments from the same sample are appended;
- exact duplicate records are deduplicated;
- equal-valued overlaps are normalized;
- conflicting overlapping values in the same sample raise an error instead of
  silently choosing a numerical rule.

### Explicit mean or sum

```{r}
merged_mean <- merge_bwg(
  track_a,
  track_b,
  method = "mean"
)

merged_sum <- merge_bwg(
  track_a,
  track_b,
  method = "sum"
)
```

For means, missing coverage is explicit:

```{r}
merged_ignore <- merge_bwg(
  track_a,
  track_b,
  method = "mean",
  missing = "ignore"
)

merged_zero <- merge_bwg(
  track_a,
  track_b,
  method = "mean",
  missing = "zero"
)
```

Replicates can be grouped with `groups` when arithmetic aggregation is needed.

## 7. Write signal files

All persistence is centralized in `write_bwg()`.

### BigWig

```{r}
manifest_bw <- write_bwg(
  region_track,
  outdir = "output_bigwig",
  format = "bigwig",
  overwrite = TRUE
)

manifest_bw
```

Native BigWig output includes zoom levels by default. Disable them explicitly:

```{r}
manifest_nozoom <- write_bwg(
  region_track,
  outdir = "output_bigwig_nozoom",
  format = "bigwig",
  overwrite = TRUE,
  zoom = FALSE
)
```

`max_zoom_levels = 0L` is also valid and produces no zoom levels.

### WIG

```{r}
manifest_wig <- write_bwg(
  region_track,
  outdir = "output_wig",
  format = "wig",
  overwrite = TRUE
)
```

### bedGraph

```{r}
manifest_bg <- write_bwg(
  region_track,
  outdir = "output_bedgraph",
  format = "bedgraph",
  overwrite = TRUE
)
```

For WIG and bedGraph, enable gzip output with:

```{r}
write_bwg(
  region_track,
  outdir = "output_gz",
  format = "bedgraph",
  overwrite = TRUE,
  compress = TRUE
)
```

The returned manifest always uses:

```{text}
sample_id  file  format
```

## 8. Cross-format workflow

The standard conversion workflow is always:

```{text}
input file
   ↓
read_bwg()
   ↓
BwgTrack
   ↓
write_bwg()
   ↓
BigWig / WIG / bedGraph
```

Do not perform format conversion inside `retrieve_bwg()`, `merge_bwg()`, or
`stats_bwg()`.

## 9. Standardized errors

Errors generated intentionally by the public bwTools API inherit from
`bwTools_error`.

Downstream packages can therefore catch them without matching a particular
message:

```{r}
tryCatch(
  stats_bwg(
    bw,
    chrom = "chr1",
    start = 0L,
    end = 100L
  ),
  bwTools_error = function(e) {
    message("bwTools input error: ", conditionMessage(e))
  }
)
```

Error messages remain human-readable, but downstream integration should prefer
the condition class when control flow depends on an error.

## 10. Advanced performance diagnostics

`benchmark_bwg()` is intended for development and performance diagnosis. It is
not required in normal analysis.

```{r}
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

For writer profiling:

```{r}
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

The benchmark reports approximate R heap behavior, not operating-system RSS.
It does not clear the operating-system file cache.

## 11. External BigWig compatibility validation

Version 0.8.4 ships an optional validation harness for independent BigWig
implementations. These tools are **not** bwTools runtime dependencies.

Supported references are:

```{text}
pyBigWig / libBigWig
rtracklayer
UCSC bigWigInfo, bigWigToBedGraph, and bedGraphToBigWig
```

Locate the installed runner from R:

```{r}
runner <- system.file(
  "compatibility",
  "run-compatibility.001.R",
  package = "bwTools"
)
runner
```

Run it from a shell:

```{bash}
Rscript /path/to/run-compatibility.001.R \
  --output bwtools-compatibility
```

The output is:

```{text}
bwtools-compatibility/compatibility-report.tsv
```

Each available implementation is checked in both directions. Missing reference
software is recorded as `SKIP`; a detected signal mismatch is `FAIL`. Use
`--strict` in a prepared release-validation environment when all supported
references must be present and pass.

The comparison is semantic rather than byte-for-byte: chromosome lengths,
missing-signal positions, and signal values must agree, while external writers
may choose different BigWig block or interval segmentation. See
`inst/COMPATIBILITY.md` for the validation contract.

## 12. API stability rules

Starting with 0.8.0, downstream packages should rely on the documented public
API only.

The following are regression-tested:

- public function names;
- public parameter names and order;
- `BwgTrack` schema v2 core fields and types;
- standard accessor return columns;
- `retrieve_bwg()` core signal columns;
- `stats_bwg()` core statistics columns;
- `write_bwg()` manifest columns;
- BigWig zoom/no-zoom full-resolution round trips;
- WIG and bedGraph compressed/uncompressed round trips;
- package-generated `bwTools_error` conditions.

Functions accessed with `bwTools:::` remain internal implementation details and
are not part of this compatibility promise.

## 13. Quick regression validation

The bundled files can be used for a short installation check:

```{r}
stopifnot(
  identical(detect_bwg_format(bw_file), "bigwig"),
  identical(detect_bwg_format(bedgraph_file), "bedgraph")
)

x <- read_bwg(
  bw_file,
  sample_ids = "example"
)

stopifnot(
  is_bwg_track(x),
  identical(summary_bwg(x)$mode, "lazy"),
  nrow(samples_bwg(x)) == 1L,
  nrow(seqinfo_bwg(x)) > 0L
)

retrieved <- retrieve_bwg(
  x,
  chrom = "chr1",
  start = 450L,
  end = 1600L
)

stopifnot(
  nrow(retrieved) > 0L,
  all(retrieved$chrom == "chr1"),
  all(retrieved$start >= 450L),
  all(retrieved$end <= 1600L)
)
```

Then run the full package test suite:

```{r}
devtools::document()
devtools::test()
devtools::check()
```

## Development status

- **0.1.x** — native reader and initial `BwgTrack` contract.
- **0.2.x** — native BigWig/WIG/bedGraph writer.
- **0.3.x** — BigWig zoom levels and interval statistics.
- **0.4.x** — standardized indexed genomic retrieval.
- **0.5.x** — automatic merge and explicit arithmetic aggregation.
- **0.6.x** — public API and `BwgTrack` schema standardization.
- **0.7.x** — large-file benchmark and reader/writer performance optimization.
- **0.8.x** — API stability, validation, cross-format consistency, and release hardening.
- **0.9.x** — release-candidate checks and API freeze.
- **1.0.0** — stable cross-platform release.

See `inst/API_CONTRACT.md`, `inst/COMPATIBILITY.md`, and `DEVELOPMENT.md` for
the formal downstream integration, external validation, and development
acceptance criteria.
