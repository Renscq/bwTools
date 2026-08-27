# bwTools

Current development release: **0.1.2**

`bwTools` is a native-R genomic signal I/O package designed to provide the
BigWig/WIG/bedGraph backend used by GeneTrackR without requiring UCSC command-line
programs or bundled C/C++ libraries.

Version 0.1.0 implements reading only. The BigWig writer, subsetting-to-file,
and merge-to-file stages are developed separately so each binary layer can be
validated before the next one is added.

## Object contract

In-memory signal data use 1-based closed coordinates and the columns:

- `sample_id`
- `chrom`
- `start`
- `end`
- `value`
- `strand`

Objects returned by `read_bwg()` inherit from `BwgTrack`, which allows existing
GeneTrackR signal plotting functions to consume them directly.

## Example

```r
library(bwTools)

track <- read_bwg(
  "sample.bigwig",
  format = "bigwig",
  mode = "lazy"
)

seqinfo_bwg(track)

signal <- retrieve_bwg(
  track,
  chrom = "chr1",
  start = 100000,
  end = 120000
)
```
