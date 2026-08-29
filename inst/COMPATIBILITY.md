# BigWig implementation compatibility

Version: **bwTools 0.8.10**

bwTools implements local BigWig binary I/O in R. External implementations are
used only as independent validation references and are not runtime dependencies.

## Compatibility targets

The compatibility harness introduced in 0.8.4 and retained in 0.8.10 supports:

```{text}
pyBigWig / libBigWig
Bioconductor rtracklayer
UCSC bigWigInfo, bigWigToBedGraph, and bedGraphToBigWig
```

The harness checks both directions where the reference implementation supports
writing:

```{text}
bwTools BigWig -> external reader
external writer -> bwTools reader
```

Both zoom-enabled and zoom-disabled bwTools outputs are supplied to external
readers. Compatibility is evaluated by chromosome lengths, missing-signal
positions, and signal values rather than byte-for-byte file identity or internal
block segmentation.

## Runtime dependency policy

External validation software must not appear in bwTools `Imports`, `Depends`,
or normal runtime code. The installed validation runner lives under
`inst/compatibility/` and detects optional tools at execution time.

A missing external implementation is a `SKIP`, not a package test failure. A
present implementation that cannot read a bwTools file, produces a file that
bwTools cannot read, or disagrees on canonical signal is a `FAIL`.

## Reproducible report

Run:

```{bash}
Rscript "$(Rscript -e 'cat(system.file("compatibility", "run-compatibility.001.R", package="bwTools"))')" \
  --output bwtools-compatibility
```

The machine-readable report is:

```{text}
bwtools-compatibility/compatibility-report.tsv
```

Use `--strict` in a prepared validation environment to require pyBigWig,
rtracklayer, and all three UCSC utilities to be present and pass.

The compatibility matrix in release documentation should only claim a PASS
when it is backed by a saved report from an environment where that
implementation was actually available.
