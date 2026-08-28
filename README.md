# bwTools

Current development release: **0.6.3**

`bwTools` provides native-R reading, writing, indexed retrieval, merging, and
statistics for BigWig, WIG, and bedGraph genomic signal tracks.

The 0.6.x series standardizes the public API and the `BwgTrack` object contract
so downstream packages can use bwTools without depending on internal list
implementation details. Version 0.6.3 keeps BigWig detection binary-safe and fixes automatic
BigWig dispatch before any text decoding is attempted.

## Core principles

- BigWig binary I/O is implemented in R.
- No UCSC `bedGraphToBigWig` executable is required.
- No Rcpp or bundled C/C++ runtime is required.
- Public in-memory coordinates are always **1-based closed**.
- BigWig and bedGraph disk coordinates remain **0-based half-open**.
- BigWig files use lazy indexed access by default.
- BigWig format detection checks the binary magic signature before text decoding.
- All analysis functions operate on a validated `BwgTrack` object.
- Only `read_bwg()` reads signal files.
- Only `write_bwg()` writes signal files.

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

The old public names `bw_track()` and `bw_detect_format()` are removed in
0.6.1. Use `bwg_track()` and `detect_bwg_format()`.

## Standard argument names

The public API follows these naming rules:

```{text}
files       file paths passed to read_bwg()
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
- **0.6.3**: fixed automatic BigWig dispatch after format detection — current stage.
- **0.7.x**: large-file profiling and performance optimization.
- **0.8.x**: statistics and zoom refinements where justified by benchmarks.
- **0.9.x**: API freeze, cross-platform checks, and release hardening.
- **1.0.0**: stable public release.
