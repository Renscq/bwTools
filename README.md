# bwTools

Current development release: **0.2.0**

`bwTools` is a native-R genomic signal I/O package for BigWig, WIG, and
bedGraph files. It is designed as the signal backend for GeneTrackR without
requiring UCSC command-line programs, Rcpp, libBigWig, or other compiled
BigWig dependencies.

The 0.2.x series adds a native BigWig writer to the 0.1.x indexed reader.
BigWig reading and writing are therefore both implemented directly in R.

## Design

- BigWig binary parsing and writing are implemented in R.
- BigWig random access uses the chromosome B+ tree and full-data R-tree.
- BigWig data blocks use zlib compression through base R.
- Public in-memory coordinates are **1-based closed**.
- BigWig and bedGraph coordinates are **0-based half-open** on disk.
- Large BigWig inputs should normally use `mode = "lazy"`.
- Returned objects inherit from `BwgTrack` for GeneTrackR interoperability.
- BigWig 0.2.x output intentionally omits zoom levels; exact interval access is
  fully supported and zoom summaries are planned for 0.3.x.

## Object contract

A `BwgTrack` contains:

- `samples`: sample metadata and source paths;
- `data`: in-memory signal intervals, or `NULL` in lazy mode;
- `seqinfo`: chromosome metadata;
- `meta`: coordinate convention, backend, mode, and format metadata.

In-memory signal records use:

```text
sample_id  chrom  start  end  value  strand
```

## Bundled example data

The package includes a three-chromosome signal dataset under `inst/extdata/`:

| File | Purpose |
| --- | --- |
| `bwtools_example.bigwig` | BigWig example with discontinuous gene-like coverage |
| `bwtools_example.bedGraph` | Plain-text reference signal |
| `bwtools_example.chrom.sizes` | Chromosome lengths |

Chromosome lengths are `chr1 = 10000`, `chr2 = 8000`, and `chr3 = 5000`.

```r
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

## Read BigWig lazily

```r
track <- read_bwg(
  bw_file,
  format = "bigwig",
  sample_names = "example",
  mode = "lazy"
)

seqinfo_bwg(track)
summary_bwg(track)

signal <- retrieve_bwg(
  track,
  chrom = "chr1",
  start = 450,
  end = 1600
)

signal
```

Lazy mode reads only R-tree blocks overlapping the requested region.

## Read BigWig into memory

```r
track_memory <- read_bwg(
  bw_file,
  format = "bigwig",
  sample_names = "example",
  mode = "memory"
)

track_memory$data
```

For genome-wide RNA-seq or Ribo-seq tracks, lazy mode is preferred.

## Write BigWig natively in R

Read the bundled bedGraph and write it back to BigWig:

```r
track <- read_bwg(
  bedgraph_file,
  format = "bedgraph",
  sample_names = "example",
  mode = "memory"
)

outdir <- file.path(tempdir(), "bwtools_bigwig")

written <- write_bwg(
  track,
  outdir = outdir,
  format = "bigwig",
  chrom_sizes = chrom_sizes_file,
  overwrite = TRUE
)

written
```

No `bedGraphToBigWig` executable is used. The writer directly creates:

1. the BigWig v4 fixed header and total summary;
2. the chromosome B+ tree;
3. zlib-compressed bedGraph-like data blocks;
4. the full-data R-tree index;
5. the terminal BigWig magic value.

The resulting file can immediately be read by `bwTools`:

```r
roundtrip <- read_bwg(
  written$file,
  format = "bigwig",
  sample_names = "example",
  mode = "lazy"
)

retrieve_bwg(
  roundtrip,
  chrom = "chr1",
  start = 450,
  end = 1600
)
```

If a `BwgTrack` already has complete chromosome lengths in `seqinfo`,
`chrom_sizes` can be omitted:

```r
memory_bw <- read_bwg(
  bw_file,
  format = "bigwig",
  sample_names = "example",
  mode = "memory"
)

write_bwg(
  memory_bw,
  outdir = outdir,
  format = "bigwig",
  overwrite = TRUE
)
```

`chrom_sizes` is deliberately required when sequence lengths are unknown. The
writer never infers chromosome length from the last covered signal interval.

## Write bedGraph or WIG

```r
write_bwg(
  track,
  outdir = outdir,
  format = "bedgraph",
  overwrite = TRUE
)

write_bwg(
  track,
  outdir = outdir,
  format = "wig",
  overwrite = TRUE
)
```

Set `compress = TRUE` to create `.gz` bedGraph or WIG files.

## Coordinate conversion

A bedGraph record:

```text
chr1    499    524    9.1000
```

is 0-based half-open and becomes the following 1-based closed in-memory record:

```text
chr1    500    524    9.1000
```

The reverse conversion is applied during writing.

## Writer validation rules

Native BigWig output requires:

- chromosome names present in `chrom_sizes`;
- positive chromosome lengths;
- finite numeric signal values;
- valid 1-based closed intervals;
- intervals not extending beyond chromosome bounds;
- non-overlapping intervals within each chromosome.

Signal records are sorted according to chromosome order in `chrom_sizes` and
then genomic position before binary encoding.

## Automatic format detection

```r
bw_detect_format(bw_file)
bw_detect_format(bedgraph_file)
```

Supported formats are BigWig (`.bw`, `.bigwig`), WIG (`.wig`), and bedGraph
(`.bedgraph`, `.bdg`).

## Development status

- **0.1.x**: native reader and stable `BwgTrack` contract — implemented.
- **0.2.x**: native BigWig/WIG/bedGraph writer — current stage.
- **0.3.x**: BigWig zoom levels and interval statistics.
- **0.4.x**: file-level subset operations.
- **0.5.x**: multi-track merge and explicit arithmetic collapse.
- **0.6.x**: GeneTrackR backend integration.
- **1.0.0**: stable API and cross-platform release.

See `DEVELOPMENT.md` for detailed acceptance criteria.
