# bwTools

Current development release: **0.5.1**

`bwTools` provides native-R I/O and indexed analysis for BigWig, WIG, and
bedGraph genomic signal tracks. It is designed as the signal backend for
GeneTrackR without requiring UCSC command-line programs, Rcpp, libBigWig, or
other compiled BigWig dependencies.

This `README.qmd` is the source document. Render it with:

```{bash}
quarto render README.qmd
```

Code execution is disabled during rendering. The `{r}` chunks below are intended
to be executed interactively in order when validating a development release.

## Design

- BigWig parsing, writing, zoom construction, and statistics are implemented in R.
- BigWig random access uses chromosome B+ tree and R-tree indexes.
- Full-resolution and zoom blocks use zlib compression through base R.
- Public in-memory coordinates are **1-based closed**.
- BigWig and bedGraph coordinates are **0-based half-open** on disk.
- Large BigWig inputs should normally use `mode = "lazy"`.
- Returned objects inherit from `BwgTrack` for GeneTrackR interoperability.

## Bundled example data

The package includes:

| File | Purpose |
| --- | --- |
| `bwtools_example.bigwig` | BigWig with discontinuous gene-like coverage |
| `bwtools_example.bedGraph` | Plain-text reference signal |
| `bwtools_example.chrom.sizes` | Chromosome lengths |

Chromosome lengths are `chr1 = 10000`, `chr2 = 8000`, and `chr3 = 5000`.

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

stopifnot(
  file.exists(bw_file),
  file.exists(bedgraph_file),
  file.exists(chrom_sizes_file)
)
```

For round-trip checks, normalize the signal representation before comparing
files. This keeps the test focused on genomic coordinates and values rather
than row names or storage-specific attributes.

```{r}
canonical_signal <- function(x) {
  columns <- c(
    "sample_id",
    "chrom",
    "start",
    "end",
    "value",
    "strand"
  )

  x <- as.data.frame(x)[, columns, drop = FALSE]
  x$sample_id <- as.character(x$sample_id)
  x$chrom <- as.character(x$chrom)
  x$start <- as.integer(x$start)
  x$end <- as.integer(x$end)
  x$value <- as.numeric(x$value)
  x$strand <- as.character(x$strand)

  x <- x[
    order(
      x$sample_id,
      x$chrom,
      x$start,
      x$end,
      x$strand
    ),
    ,
    drop = FALSE
  ]
  row.names(x) <- NULL
  x
}

assert_signal_equal <- function(expected, observed, tolerance = 1e-4) {
  expected <- canonical_signal(expected)
  observed <- canonical_signal(observed)
  comparison <- all.equal(
    expected,
    observed,
    tolerance = tolerance,
    check.attributes = FALSE
  )

  if (!isTRUE(comparison)) {
    message("Signal round-trip comparison failed:")
    print(comparison)

    if (nrow(expected) != nrow(observed)) {
      message(
        "Row count: expected=",
        nrow(expected),
        ", observed=",
        nrow(observed)
      )
    }

    common_n <- min(nrow(expected), nrow(observed))
    if (common_n > 0L) {
      unequal <- vapply(
        seq_len(common_n),
        function(i) {
          !isTRUE(
            all.equal(
              expected[i, , drop = FALSE],
              observed[i, , drop = FALSE],
              tolerance = tolerance,
              check.attributes = FALSE
            )
          )
        },
        logical(1L)
      )
      first_bad <- which(unequal)[1L]
      if (!is.na(first_bad)) {
        message("First differing row: ", first_bad)
        print(expected[first_bad, , drop = FALSE])
        print(observed[first_bad, , drop = FALSE])
      }
    }

    stop("Signal round-trip validation failed.", call. = FALSE)
  }

  invisible(TRUE)
}
```

## Reader validation

### 1. Detect formats

```{r}
stopifnot(
  identical(bw_detect_format(bw_file), "bigwig"),
  identical(bw_detect_format(bedgraph_file), "bedgraph")
)
```

### 2. Inspect chromosome metadata

```{r}
seqinfo <- seqinfo_bwg(bw_file)
print(seqinfo)

stopifnot(
  identical(seqinfo$chrom, c("chr1", "chr2", "chr3")),
  identical(seqinfo$length, c(10000L, 8000L, 5000L))
)
```

### 3. Read BigWig lazily

```{r}
track_lazy <- read_bwg(
  bw_file,
  format = "bigwig",
  sample_names = "example",
  mode = "lazy"
)

print(track_lazy)
summary_bwg(track_lazy)

stopifnot(
  is_bwg_track(track_lazy),
  is.null(track_lazy$data),
  nrow(track_lazy$samples) == 1L,
  nrow(track_lazy$seqinfo) == 3L
)
```

### 4. Retrieve a region

```{r}
signal_lazy <- retrieve_bwg(
  track_lazy,
  chrom = "chr1",
  start = 450,
  end = 1600
)

print(signal_lazy)

stopifnot(
  nrow(signal_lazy) > 0L,
  all(signal_lazy$chrom == "chr1"),
  all(signal_lazy$start >= 450L),
  all(signal_lazy$end <= 1600L)
)
```

### 5. Compare BigWig with bedGraph

```{r}
track_bedgraph <- read_bwg(
  bedgraph_file,
  format = "bedgraph",
  sample_names = "example",
  mode = "memory"
)

signal_bedgraph <- retrieve_bwg(
  track_bedgraph,
  chrom = "chr1",
  start = 450,
  end = 1600
)

compare_columns <- c(
  "sample_id",
  "chrom",
  "start",
  "end",
  "value",
  "strand"
)

bw_observed <- as.data.frame(signal_lazy)[, compare_columns, drop = FALSE]
bedgraph_observed <- as.data.frame(signal_bedgraph)[, compare_columns, drop = FALSE]

stopifnot(
  isTRUE(
    all.equal(
      bw_observed,
      bedgraph_observed,
      tolerance = 1e-4
    )
  )
)
```

## Writer validation

### 6. Write BigWig natively

```{r}
test_root <- file.path(tempdir(), "bwtools_readme_test")

if (dir.exists(test_root)) {
  unlink(test_root, recursive = TRUE)
}

dir.create(test_root, recursive = TRUE)

written_bigwig <- write_bwg(
  track_bedgraph,
  outdir = test_root,
  format = "bigwig",
  chrom_sizes = chrom_sizes_file,
  overwrite = TRUE
)

print(written_bigwig)

stopifnot(
  nrow(written_bigwig) == 1L,
  file.exists(written_bigwig$file),
  identical(bw_detect_format(written_bigwig$file), "bigwig")
)
```

The native writer creates the BigWig v4 header, chromosome B+ tree,
full-resolution data blocks, full-data R-tree, optional zoom summary blocks,
zoom R-trees, total summary, and terminal magic value entirely in R.

### 7. Validate BigWig round-trip

```{r}
roundtrip_track <- read_bwg(
  written_bigwig$file,
  format = "bigwig",
  sample_names = "example",
  mode = "memory"
)

original_data <- as.data.frame(track_bedgraph$data)
roundtrip_data <- as.data.frame(roundtrip_track$data)

original_data <- original_data[
  order(original_data$chrom, original_data$start, original_data$end),
  ,
  drop = FALSE
]
roundtrip_data <- roundtrip_data[
  order(roundtrip_data$chrom, roundtrip_data$start, roundtrip_data$end),
  ,
  drop = FALSE
]

row.names(original_data) <- NULL
row.names(roundtrip_data) <- NULL

stopifnot(
  isTRUE(
    all.equal(
      original_data,
      roundtrip_data,
      tolerance = 1e-4
    )
  )
)
```

## Zoom-level validation

Do not use the internal `bw_metadata()` structure in user code. Use the public
`zoominfo_bwg()` interface.

### 8. Inspect zoom levels

```{r}
zoom_info <- zoominfo_bwg(written_bigwig$file)
print(zoom_info)

stopifnot(
  nrow(zoom_info) > 0L,
  all(zoom_info$level > 0),
  all(zoom_info$data_offset > 0),
  all(zoom_info$index_offset > zoom_info$data_offset)
)
```

The same information can be obtained from a lazy `BwgTrack`:

```{r}
roundtrip_lazy <- read_bwg(
  written_bigwig$file,
  format = "bigwig",
  sample_names = "example",
  mode = "lazy"
)

zoom_track <- zoominfo_bwg(roundtrip_lazy)
print(zoom_track)

stopifnot(
  nrow(zoom_track) == nrow(zoom_info),
  all(zoom_track$sample_id == "example"),
  identical(zoom_track$level, zoom_info$level)
)
```

### 9. Validate `zoom = FALSE`

```{r}
nozoom_dir <- file.path(test_root, "nozoom")

written_nozoom <- write_bwg(
  track_bedgraph,
  outdir = nozoom_dir,
  format = "bigwig",
  chrom_sizes = chrom_sizes_file,
  overwrite = TRUE,
  zoom = FALSE
)

nozoom_info <- zoominfo_bwg(written_nozoom$file)
print(nozoom_info)

stopifnot(nrow(nozoom_info) == 0L)
```

## Statistics validation

`stats_bwg()` supports `mean`, `stdev`, `max`, `min`, `coverage`, and `sum`.

### 10. Exact full-resolution statistics

```{r}
stats_full <- stats_bwg(
  roundtrip_lazy,
  chrom = "chr1",
  start = 1,
  end = 10000,
  n_bins = 10,
  stat = "mean",
  use_zoom = FALSE
)

print(stats_full)

stopifnot(
  nrow(stats_full) == 10L,
  all(stats_full$source == "full")
)
```

### 11. Zoom-accelerated statistics

```{r}
stats_zoom <- stats_bwg(
  roundtrip_lazy,
  chrom = "chr1",
  start = 1,
  end = 10000,
  n_bins = 10,
  stat = "mean",
  use_zoom = TRUE
)

print(stats_zoom)

stopifnot(
  nrow(stats_zoom) == 10L,
  all(stats_zoom$source %in% c("full", "zoom"))
)
```

For exact scientific calculations, use `use_zoom = FALSE`. Partial overlap with
a stored zoom window is approximate because only summary values are available at
that resolution.

## WIG and bedGraph round-trip

### 12. bedGraph

```{r}
bedgraph_outdir <- file.path(test_root, "bedgraph")

written_bedgraph <- write_bwg(
  track_bedgraph,
  outdir = bedgraph_outdir,
  format = "bedgraph",
  overwrite = TRUE
)

bedgraph_roundtrip <- read_bwg(
  written_bedgraph$file,
  format = "bedgraph",
  sample_names = "example",
  mode = "memory"
)

assert_signal_equal(
  track_bedgraph$data,
  bedgraph_roundtrip$data
)
```

### 13. WIG

```{r}
wig_outdir <- file.path(test_root, "wig")

written_wig <- write_bwg(
  track_bedgraph,
  outdir = wig_outdir,
  format = "wig",
  overwrite = TRUE
)

wig_roundtrip <- read_bwg(
  written_wig$file,
  format = "wig",
  sample_names = "example",
  mode = "memory"
)

assert_signal_equal(
  track_bedgraph$data,
  wig_roundtrip$data
)
```


## Multi-region retrieval

`retrieve_bwg()` is the only genomic extraction interface. It supports both the
original single-region API and a `regions` table for multi-region retrieval.
Extraction never writes files and does not expose writer-specific controls such
as `format`, `compress`, `zoom`, or `max_zoom_levels`.

Overlapping or directly adjacent query regions are merged before retrieval.
Signal intervals crossing a query boundary are clipped while keeping their
original genomic coordinate system.

### 14. Retrieve multiple regions as signal data

```{r}
retrieve_regions <- data.frame(
  chrom = c("chr1", "chr1", "chr2"),
  start = c(510L, 561L, 450L),
  end = c(560L, 600L, 900L)
)

retrieved_data <- retrieve_bwg(
  roundtrip_lazy,
  regions = retrieve_regions
)

print(retrieved_data)

stopifnot(
  nrow(retrieved_data) > 0L,
  all(retrieved_data$chrom %in% c("chr1", "chr2"))
)
```

The two adjacent `chr1` queries are normalized to the genomic union
`chr1:510-600`, preventing duplicate signal records.

### 15. Retrieve a write-ready BwgTrack

Use `result = "track"` when the extracted signal should remain inside the
standard `BwgTrack` object contract:

```{r}
retrieved_track <- retrieve_bwg(
  roundtrip_lazy,
  regions = retrieve_regions,
  result = "track"
)

print(retrieved_track)
print(retrieved_track$meta$retrieval_regions)

stopifnot(
  is_bwg_track(retrieved_track),
  identical(retrieved_track$meta$mode, "memory"),
  isTRUE(retrieved_track$meta$retrieved),
  nrow(retrieved_track$meta$retrieval_regions) == 2L
)
```

The default `result = "data"` is retained for compatibility with existing
GeneTrackR code that expects the signal data.table directly.

### 16. Persist retrieved data only with write_bwg()

File output is deliberately separated from retrieval:

```{r}
retrieved_bigwig_dir <- file.path(
  test_root,
  "retrieved_bigwig"
)

retrieved_bigwig <- write_bwg(
  retrieved_track,
  outdir = retrieved_bigwig_dir,
  format = "bigwig",
  overwrite = TRUE
)

stopifnot(
  nrow(retrieved_bigwig) == 1L,
  file.exists(retrieved_bigwig$file)
)
```

Read the new file back and compare it with the retrieved in-memory track:

```{r}
retrieved_roundtrip <- read_bwg(
  retrieved_bigwig$file,
  format = "bigwig",
  sample_names = "example",
  mode = "memory"
)

assert_signal_equal(
  retrieved_track$data,
  retrieved_roundtrip$data
)
```

This separation keeps responsibilities explicit:

```{text}
read_bwg()      open/load tracks
retrieve_bwg()  extract genomic regions
merge_bwg()     merge or explicitly aggregate tracks
stats_bwg()     calculate interval statistics
write_bwg()     persist or convert tracks
```

## Merge signal tracks

`merge_bwg()` is the single public interface for combining signal tracks.
Users do not need to determine in advance whether inputs contain the same
sample IDs or whether genomic regions overlap.

The default is:

```{text}
merged <- merge_bwg(track1, track2, track3)
```

which is equivalent to:

```{text}
merged <- merge_bwg(
  track1,
  track2,
  track3,
  method = "auto"
)
```

`method = "auto"` is deliberately conservative:

| Input relationship | Automatic behavior |
| --- | --- |
| Different samples | preserve samples independently |
| Same sample, disjoint regions | append records |
| Exact duplicate records | deduplicate automatically |
| Same sample, equal-valued overlaps | merge overlapping signal |
| Same sample, conflicting overlap values | stop and request `mean` or `sum` |

Arithmetic aggregation is therefore never guessed from coordinates alone.
Users only choose `method = "mean"` or `method = "sum"` when changing signal
values is intentional.

### 17. Same sample, different regions

```{r}
merge_a <- bw_track(
  samples = data.frame(
    sample_id = "sampleA",
    file = NA_character_,
    format = "memory",
    strand = "*"
  ),
  data = data.frame(
    sample_id = "sampleA",
    chrom = "chr1",
    start = c(100L, 500L),
    end = c(199L, 599L),
    value = c(2, 4),
    strand = "*"
  ),
  seqinfo = data.frame(
    sample_id = "sampleA",
    chrom = "chr1",
    length = 10000L
  ),
  meta = list(mode = "memory")
)

merge_b <- bw_track(
  samples = data.frame(
    sample_id = "sampleA",
    file = NA_character_,
    format = "memory",
    strand = "*"
  ),
  data = data.frame(
    sample_id = "sampleA",
    chrom = "chr2",
    start = c(200L, 800L),
    end = c(299L, 899L),
    value = c(3, 5),
    strand = "*"
  ),
  seqinfo = data.frame(
    sample_id = "sampleA",
    chrom = "chr2",
    length = 8000L
  ),
  meta = list(mode = "memory")
)

same_sample_merged <- merge_bwg(
  merge_a,
  merge_b
)

stopifnot(
  same_sample_merged$meta$method == "auto",
  nrow(same_sample_merged$samples) == 1L,
  same_sample_merged$samples$sample_id == "sampleA",
  nrow(same_sample_merged$data) == 4L,
  identical(
    unique(same_sample_merged$data$chrom),
    c("chr1", "chr2")
  )
)
```

No method selection was required because the regions are compatible.

### 18. Equal-valued overlaps are resolved automatically

Exact duplicates and same-valued overlaps do not require another option.

```{r}
equal_overlap_a <- bw_track(
  samples = merge_a$samples,
  data = data.frame(
    sample_id = "sampleA",
    chrom = "chr1",
    start = 1L,
    end = 10L,
    value = 2,
    strand = "*"
  ),
  seqinfo = merge_a$seqinfo,
  meta = list(mode = "memory")
)

equal_overlap_b <- bw_track(
  samples = merge_a$samples,
  data = data.frame(
    sample_id = "sampleA",
    chrom = "chr1",
    start = 6L,
    end = 15L,
    value = 2,
    strand = "*"
  ),
  seqinfo = merge_a$seqinfo,
  meta = list(mode = "memory")
)

equal_overlap_merged <- merge_bwg(
  equal_overlap_a,
  equal_overlap_b
)

stopifnot(
  nrow(equal_overlap_merged$data) == 1L,
  equal_overlap_merged$data$start == 1L,
  equal_overlap_merged$data$end == 15L,
  equal_overlap_merged$data$value == 2
)
```

### 19. Conflicting overlaps require an explicit method

Consider two inputs from the same sample:

```{text}
source 1: chr1 1-10  value=2
source 2: chr1 6-15  value=4
```

The overlap contains incompatible values. The default merge stops rather than
silently selecting, averaging, or adding them.

```{r}
conflict_a <- bw_track(
  samples = merge_a$samples,
  data = data.frame(
    sample_id = "sampleA",
    chrom = "chr1",
    start = 1L,
    end = 10L,
    value = 2,
    strand = "*"
  ),
  seqinfo = merge_a$seqinfo,
  meta = list(mode = "memory")
)

conflict_b <- bw_track(
  samples = merge_a$samples,
  data = data.frame(
    sample_id = "sampleA",
    chrom = "chr1",
    start = 6L,
    end = 15L,
    value = 4,
    strand = "*"
  ),
  seqinfo = merge_a$seqinfo,
  meta = list(mode = "memory")
)

stopifnot(
  inherits(
    try(merge_bwg(conflict_a, conflict_b), silent = TRUE),
    "try-error"
  )
)
```

If addition is biologically intended:

```{r}
conflict_sum <- merge_bwg(
  conflict_a,
  conflict_b,
  method = "sum"
)

stopifnot(
  identical(conflict_sum$data$start, c(1L, 6L, 11L)),
  identical(conflict_sum$data$end, c(5L, 10L, 15L)),
  isTRUE(all.equal(conflict_sum$data$value, c(2, 6, 4)))
)
```

If averaging is intended, the default `missing = "ignore"` averages only
members that cover each atomic segment:

```{r}
conflict_mean <- merge_bwg(
  conflict_a,
  conflict_b,
  method = "mean"
)

stopifnot(
  isTRUE(
    all.equal(
      conflict_mean$data$value,
      c(2, 3, 4)
    )
  )
)
```

Use `missing = "zero"` only when uncovered members should contribute zero:

```{r}
conflict_mean_zero <- merge_bwg(
  conflict_a,
  conflict_b,
  method = "mean",
  missing = "zero"
)

stopifnot(
  isTRUE(
    all.equal(
      conflict_mean_zero$data$value,
      c(1, 3, 2)
    )
  )
)
```

### 20. Different samples are preserved automatically

Different sample IDs remain independent in auto mode even when their genomic
coordinates overlap.

```{r}
sample_a_signal <- bw_track(
  samples = data.frame(
    sample_id = "sampleA",
    file = NA_character_,
    format = "memory",
    strand = "*"
  ),
  data = data.frame(
    sample_id = "sampleA",
    chrom = "chr1",
    start = 1L,
    end = 10L,
    value = 2,
    strand = "*"
  ),
  seqinfo = data.frame(
    sample_id = "sampleA",
    chrom = "chr1",
    length = 10000L
  ),
  meta = list(mode = "memory")
)

sample_b_signal <- bw_track(
  samples = data.frame(
    sample_id = "sampleB",
    file = NA_character_,
    format = "memory",
    strand = "*"
  ),
  data = data.frame(
    sample_id = "sampleB",
    chrom = "chr1",
    start = 1L,
    end = 10L,
    value = 4,
    strand = "*"
  ),
  seqinfo = data.frame(
    sample_id = "sampleB",
    chrom = "chr1",
    length = 10000L
  ),
  meta = list(mode = "memory")
)

multi_sample <- merge_bwg(
  sample_a_signal,
  sample_b_signal
)

stopifnot(
  identical(
    multi_sample$samples$sample_id,
    c("sampleA", "sampleB")
  ),
  nrow(multi_sample$data) == 2L
)
```

Coordinates alone therefore never cause different samples to be averaged.

### 21. Explicitly aggregate different samples

When `mean` or `sum` is selected and several sample IDs are present, they are
combined into an output sample named `merged` unless `groups` is supplied.

```{r}
combined_mean <- merge_bwg(
  sample_a_signal,
  sample_b_signal,
  method = "mean"
)

stopifnot(
  combined_mean$samples$sample_id == "merged",
  combined_mean$data$value == 3
)
```

```{r}
combined_sum <- merge_bwg(
  sample_a_signal,
  sample_b_signal,
  method = "sum"
)

stopifnot(
  combined_sum$samples$sample_id == "merged",
  combined_sum$data$value == 6
)
```

### 22. Group replicate samples

Use `groups` when arithmetic merging should occur within biological groups.

```{r}
groups <- data.frame(
  sample_id = c("sampleA", "sampleB"),
  group_id = c("treatment", "treatment")
)

group_mean <- merge_bwg(
  sample_a_signal,
  sample_b_signal,
  method = "mean",
  groups = groups
)

stopifnot(
  group_mean$samples$sample_id == "treatment",
  group_mean$data$value == 3
)
```

A named character vector is also accepted:

```{r}
group_mean_named <- merge_bwg(
  sample_a_signal,
  sample_b_signal,
  method = "mean",
  groups = c(
    sampleA = "treatment",
    sampleB = "treatment"
  )
)

stopifnot(
  group_mean_named$data$value == 3
)
```

### Merge rules

- `method = "auto"` is the default and does not perform arithmetic.
- Different sample IDs remain independent in auto mode regardless of coordinate
  overlap.
- Same-sample disjoint regions are appended automatically.
- Exact duplicates and equal-valued same-sample overlaps are consolidated
  automatically.
- Conflicting same-sample overlaps require explicit `method = "mean"` or
  `method = "sum"`.
- `missing = "ignore"` is the default for means; use `missing = "zero"` only
  when uncovered members should count as zero.
- `groups` is used only with arithmetic methods.
- `strand = "separate"` is the default; opposite strands are not aggregated
  unless `strand = "ignore"` is explicit.
- Chromosome-length conflicts are errors rather than being silently resolved.
- Adjacent equal-valued output segments are merged by default.
- Zero-valued arithmetic segments are omitted by default through
  `drop_zero = TRUE`; auto mode preserves zero-valued input signal.
- `merge_bwg()` always returns a memory-mode `BwgTrack`. File output remains the
  responsibility of `write_bwg()`.

Persist a merged result explicitly:

```{r}
merge_output_dir <- file.path(
  test_root,
  "merge_output"
)

merge_files <- write_bwg(
  same_sample_merged,
  outdir = merge_output_dir,
  format = "bigwig",
  overwrite = TRUE
)

stopifnot(
  nrow(merge_files) == 1L,
  file.exists(merge_files$file)
)
```

## Coordinate contract

A bedGraph record:

```{text}
chr1    499    524    9.1000
```

is 0-based half-open and becomes:

```{text}
chr1    500    524    9.1000
```

in the public 1-based closed in-memory representation. The reverse conversion is
applied when writing BigWig or bedGraph.

## Package regression tests

Run after the manual checks:

```{r}
devtools::document()
devtools::test()
```

Then run:

```{r}
devtools::check()
```

Expected test status:

```{text}
FAIL 0
WARN 0
SKIP 0
```

Any reader, writer, zoom, statistics, object-contract, or text-format failure
should be treated as a regression before continuing development.

## Development status

- **0.1.x**: native reader and stable `BwgTrack` contract — implemented.
- **0.2.x**: native BigWig/WIG/bedGraph writer — implemented.
- **0.3.x**: BigWig zoom levels and interval statistics — implemented.
- **0.4.x**: unified indexed single-/multi-region retrieval — implemented.
- **0.5.x**: unified automatic merge and explicit arithmetic aggregation — current stage.
- **0.6.x**: GeneTrackR backend integration.
- **1.0.0**: stable API and cross-platform release.

See `DEVELOPMENT.md` for detailed acceptance criteria.
