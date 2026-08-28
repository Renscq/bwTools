# bwTools public API contract

Version: **schema 2 / bwTools 0.7.3**

Downstream packages should treat the following as the stable integration
boundary:

```{text}
read_bwg()
bwg_track()
validate_bwg()
is_bwg_track()
samples_bwg()
seqinfo_bwg()
metadata_bwg()
retrieve_bwg()
merge_bwg()
stats_bwg()
zoominfo_bwg()
summary_bwg()
write_bwg()
detect_bwg_format()
benchmark_bwg()
```

All normal analysis functions consume a validated `BwgTrack`. File paths enter
through `read_bwg()` and signal files leave through `write_bwg()`.

`benchmark_bwg()` is the only diagnostic exception: it accepts one local
BigWig file because its purpose is to measure file-I/O behavior itself. It
returns a `bwToolsBenchmark` object and does not persist benchmark reports.

## Object schema

```{text}
BwgTrack
├── samples: sample_id, file, format, strand, has_strand
├── data: sample_id, chrom, start, end, value, strand OR NULL for lazy mode
├── seqinfo: sample_id, chrom, length
└── meta: coordinate, schema_version, backend, mode
```

In-memory coordinates are always 1-based closed.

Core schema types are fixed after construction:

```{text}
sample_id, file, format, strand, chrom  character
has_strand                            logical
start, end, length                    integer
value                                 numeric
```

Use `bwg_track()` to normalize externally created tables before passing them to
other bwTools functions.

Downstream packages should use public accessors instead of private
`bwTools:::` helpers or direct assumptions about internal implementation.
Use `metadata_bwg()` for schema, operation, and provenance metadata.


## Large-file benchmark contract

`benchmark_bwg()` returns:

```{text}
bwToolsBenchmark
├── system
├── file
├── config
├── regions
├── results
└── summary
```

The benchmark does not alter the `BwgTrack` schema and is not required for
normal downstream integration. Default profiling avoids full-memory BigWig
loading. Exact statistics above `full_stats_max_bases` are recorded as skipped. Optional
writer profiling is bounded by `write_max_bases`, profiles each selected region
with zoom disabled and enabled, and records writer-specific throughput and
output-size diagnostics. Writer timing excludes region retrieval. Summary rows for
`zoom_on` include elapsed-time and output-size overhead relative to the matching
`zoom_off` region. Peak-memory reporting is an approximate R heap maximum from `gc()` rather than
process RSS, and package-metadata-cold runs do not clear the operating-system
file cache.
