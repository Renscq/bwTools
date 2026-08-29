# Retrieve signal intervals

Retrieves one or more genomic regions from a standardized `BwgTrack`.
Query coordinates are 1-based closed. Overlapping or directly adjacent
regions are merged before retrieval, and signal intervals crossing query
boundaries are clipped without rebasing genomic coordinates. When
complete chromosome lengths are available in `seqinfo`, query ends
beyond a chromosome are clipped to the known length and queries starting
beyond the chromosome return an empty result. Unknown chromosomes are
rejected only for selected samples whose `seqinfo` lengths are all
known; length-incomplete metadata remains permissive because it may
describe only chromosomes observed in the signal.

## Usage

``` r
retrieve_bwg(
  x,
  chrom = NULL,
  start = NULL,
  end = NULL,
  regions = NULL,
  sample_ids = NULL,
  strand = NULL,
  result = c("data", "track")
)
```

## Arguments

- x:

  A validated `BwgTrack` object.

- chrom:

  Optional chromosome name for a single-region query.

- start:

  Optional start coordinate for a single-region query.

- end:

  Optional end coordinate for a single-region query.

- regions:

  Optional data.frame or data.table containing `chrom`, `start`, and
  `end` columns. Use either `regions` or `chrom`/`start`/`end`.

- sample_ids:

  Optional sample IDs. Defaults to all samples.

- strand:

  Optional exact strand selector: `+`, `-`, or `*`. `NULL` performs no
  strand filtering. Memory-mode tracks filter row-level signal strand
  values; lazy BigWig tracks use sample-level strand metadata.

- result:

  Result type. `data` returns the standardized signal data.table;
  `track` returns a memory-mode `BwgTrack`.

## Value

A signal data.table or memory-mode `BwgTrack`.
