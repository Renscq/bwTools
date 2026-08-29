# External BigWig compatibility validation

This directory contains an optional development-time compatibility harness for
bwTools. None of the external implementations below are runtime dependencies of
bwTools.

The harness validates the same canonical signal in both directions whenever the
corresponding implementation is available:

| Implementation | bwTools -> external | external -> bwTools |
| --- | --- | --- |
| pyBigWig | open and compare chromosome sizes and signal | write with pyBigWig, read with bwTools |
| rtracklayer | import and compare signal | export BigWig, read with bwTools |
| UCSC utilities | `bigWigInfo` + `bigWigToBedGraph` | `bedGraphToBigWig` + bwTools read |

Unavailable implementations are recorded as `SKIP`. A detected implementation
that disagrees with the canonical signal is recorded as `FAIL` and makes the
runner return a non-zero exit status.

## Run

Locate the installed runner:

```{r}
runner <- system.file(
  "compatibility",
  "run-compatibility.001.R",
  package = "bwTools"
)
runner
```

Then run it from a shell:

```{bash}
Rscript /path/to/run-compatibility.001.R --output bwtools-compatibility
```

Use `--strict` when all three external implementations are expected to be
installed. In strict mode, any unavailable implementation also produces a
non-zero exit status.

```{bash}
Rscript /path/to/run-compatibility.001.R \
  --output bwtools-compatibility \
  --strict
```

The result is written to:

```{text}
bwtools-compatibility/compatibility-report.tsv
```

The canonical inputs are the installed files
`extdata/bwtools_example.bedGraph` and `extdata/bwtools_example.chrom.sizes`.
Signal comparisons are semantic: chromosome lengths, missing-signal positions,
and per-base values are compared, so an external writer is allowed to use a
different internal BigWig block or interval segmentation.
