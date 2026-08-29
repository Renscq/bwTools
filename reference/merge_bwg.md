# Merge genomic signal tracks

Combines two or more `BwgTrack` objects through one user-facing
interface. The default `method = "auto"` preserves sample identities and
signal values without requiring users to determine whether inputs
contain the same samples or overlapping genomic regions.

## Usage

``` r
merge_bwg(
  ...,
  tracks = NULL,
  method = c("auto", "mean", "sum"),
  groups = NULL,
  missing = c("ignore", "zero"),
  strand = c("separate", "ignore"),
  drop_zero = TRUE,
  merge_adjacent = TRUE
)
```

## Arguments

- ...:

  Two or more `BwgTrack`-compatible objects.

- tracks:

  Optional list of `BwgTrack`-compatible objects. Use either `...` or
  `tracks`.

- method:

  Merge strategy. `auto` performs conservative, non-arithmetic merging;
  `mean` and `sum` perform explicit arithmetic aggregation.

- groups:

  Optional named character vector or table with `sample_id` and
  `group_id` columns. Used only for arithmetic methods.

- missing:

  Mean denominator rule. `ignore` averages only members covering an
  atomic segment; `zero` treats uncovered members as zero. Ignored for
  `auto` and `sum`.

- strand:

  Strand handling rule. `separate` keeps `+`, `-`, and `*` signals
  separate; `ignore` combines strands.

- drop_zero:

  Whether zero-valued arithmetic output segments should be omitted.
  Ignored for `method = "auto"`.

- merge_adjacent:

  Whether directly adjacent equal-valued output segments should be
  merged.

## Value

A memory-mode `BwgTrack` containing merged signal.

## Details

In auto mode, different samples remain separate, disjoint regions from
the same sample are appended, exact duplicates and equal-valued overlaps
are consolidated automatically, and conflicting overlaps from the same
sample stop with an instruction to choose `method = "mean"` or
`method = "sum"`.

Arithmetic methods split partial overlaps at every input boundary.
Without `groups`, one unique sample ID is preserved, while multiple
sample IDs are aggregated into a new sample named `merged`. Use `groups`
to aggregate replicates or treatments into explicit output groups.

File persistence is intentionally excluded. Save the returned `BwgTrack`
explicitly with
[`write_bwg()`](https://renscq.github.io/bwTools/reference/write_bwg.md).
