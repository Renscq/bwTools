# bwTools public API contract

Version: **schema 2 / bwTools 0.6.2**

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
```

All analysis functions consume a validated `BwgTrack`. File paths enter through
`read_bwg()` and leave through `write_bwg()`.

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
