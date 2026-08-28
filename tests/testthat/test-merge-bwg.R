make_merge_track <- function(
  sample_id,
  chrom,
  start,
  end,
  value,
  strand = "*",
  chrom_length = NULL
) {
  if (length(strand) == 1L) strand <- rep(strand, length(chrom))
  signal <- data.table::data.table(
    sample_id = sample_id,
    chrom = chrom,
    start = as.integer(start),
    end = as.integer(end),
    value = as.numeric(value),
    strand = strand
  )
  sample_strand <- unique(strand)
  if (length(sample_strand) != 1L) sample_strand <- "*"
  samples <- data.table::data.table(
    sample_id = sample_id,
    file = NA_character_,
    format = "memory",
    strand = sample_strand,
    has_strand = sample_strand %in% c("+", "-")
  )
  seqinfo <- NULL
  if (!is.null(chrom_length)) {
    chrom_names <- unique(chrom)
    if (length(chrom_length) == 1L) {
      chrom_length <- rep(chrom_length, length(chrom_names))
    }
    seqinfo <- data.table::data.table(
      sample_id = sample_id,
      chrom = chrom_names,
      length = as.integer(chrom_length)
    )
  }
  bwg_track(samples, signal, seqinfo = seqinfo, meta = list(mode = "memory"))
}

test_that("auto is the default merge strategy", {
  expect_identical(eval(formals(merge_bwg)$method), c("auto", "mean", "sum"))
})

test_that("auto appends disjoint regions from the same sample", {
  a <- make_merge_track(
    "sampleA", "chr1", c(1L, 20L), c(10L, 30L), c(1, 2),
    chrom_length = 100L
  )
  b <- make_merge_track(
    "sampleA", "chr2", c(5L, 40L), c(15L, 50L), c(3, 4),
    chrom_length = 80L
  )

  observed <- merge_bwg(a, b)

  expect_s3_class(observed, "BwgTrack")
  expect_equal(observed$samples$sample_id, "sampleA")
  expect_equal(nrow(observed$data), 4L)
  expect_equal(observed$data$chrom, c("chr1", "chr1", "chr2", "chr2"))
  expect_equal(observed$meta$operation, "merge")
  expect_equal(observed$meta$method, "auto")
})

test_that("auto preserves different samples even when coordinates overlap", {
  a <- make_merge_track("sampleA", "chr1", 1L, 10L, 2, chrom_length = 100L)
  b <- make_merge_track("sampleB", "chr1", 1L, 10L, 4, chrom_length = 100L)

  observed <- merge_bwg(a, b)

  expect_equal(observed$samples$sample_id, c("sampleA", "sampleB"))
  expect_equal(nrow(observed$data), 2L)
  expect_equal(observed$data$value, c(2, 4))
})

test_that("auto deduplicates and joins equal-valued overlaps", {
  a <- make_merge_track("sampleA", "chr1", 1L, 10L, 2, chrom_length = 100L)
  duplicate <- make_merge_track("sampleA", "chr1", 1L, 10L, 2, chrom_length = 100L)
  overlap <- make_merge_track("sampleA", "chr1", 6L, 15L, 2, chrom_length = 100L)

  observed <- merge_bwg(a, duplicate, overlap)

  expect_equal(nrow(observed$data), 1L)
  expect_equal(observed$data$start, 1L)
  expect_equal(observed$data$end, 15L)
  expect_equal(observed$data$value, 2)
})

test_that("auto rejects conflicting same-sample overlaps", {
  a <- make_merge_track("sampleA", "chr1", 1L, 10L, 2, chrom_length = 100L)
  b <- make_merge_track("sampleA", "chr1", 6L, 15L, 4, chrom_length = 100L)

  expect_error(
    merge_bwg(a, b),
    "method = 'mean'.*method = 'sum'"
  )
})

test_that("sum and mean resolve same-sample overlaps", {
  a <- make_merge_track("sampleA", "chr1", 1L, 10L, 2, chrom_length = 100L)
  b <- make_merge_track("sampleA", "chr1", 6L, 15L, 4, chrom_length = 100L)

  summed <- merge_bwg(a, b, method = "sum")
  expect_equal(summed$samples$sample_id, "sampleA")
  expect_equal(summed$data$start, c(1L, 6L, 11L))
  expect_equal(summed$data$end, c(5L, 10L, 15L))
  expect_equal(summed$data$value, c(2, 6, 4))

  mean_ignore <- merge_bwg(a, b, method = "mean")
  expect_equal(mean_ignore$data$value, c(2, 3, 4))

  mean_zero <- merge_bwg(a, b, method = "mean", missing = "zero")
  expect_equal(mean_zero$data$value, c(1, 3, 2))
})

test_that("arithmetic merge combines different samples automatically", {
  a <- make_merge_track("sampleA", "chr1", 1L, 10L, 2, chrom_length = 100L)
  b <- make_merge_track("sampleB", "chr1", 1L, 10L, 4, chrom_length = 100L)

  averaged <- merge_bwg(a, b, method = "mean")
  expect_equal(averaged$samples$sample_id, "merged")
  expect_equal(averaged$data$value, 3)
  expect_equal(averaged$meta$grouping, "all")

  summed <- merge_bwg(a, b, method = "sum")
  expect_equal(summed$samples$sample_id, "merged")
  expect_equal(summed$data$value, 6)
})

test_that("groups support replicate-level arithmetic merging", {
  a <- make_merge_track("A1", "chr1", 1L, 10L, 2, chrom_length = 100L)
  b <- make_merge_track("A2", "chr1", 1L, 10L, 4, chrom_length = 100L)
  c <- make_merge_track("B1", "chr1", 1L, 10L, 10, chrom_length = 100L)
  groups <- data.frame(
    sample_id = c("A1", "A2", "B1"),
    group_id = c("groupA", "groupA", "groupB")
  )

  observed <- merge_bwg(
    a,
    b,
    c,
    method = "mean",
    groups = groups
  )

  expect_equal(observed$samples$sample_id, c("groupA", "groupB"))
  expect_equal(observed$data$value, c(3, 10))
})

test_that("groups are rejected in auto mode", {
  a <- make_merge_track("sampleA", "chr1", 1L, 10L, 2, chrom_length = 100L)
  b <- make_merge_track("sampleB", "chr1", 1L, 10L, 4, chrom_length = 100L)

  expect_error(
    merge_bwg(
      a,
      b,
      groups = c(sampleA = "group", sampleB = "group")
    ),
    "can only be used"
  )
})

test_that("chromosome length conflicts are rejected for every method", {
  a <- make_merge_track("sampleA", "chr1", 1L, 10L, 2, chrom_length = 100L)
  b <- make_merge_track("sampleA", "chr1", 20L, 30L, 4, chrom_length = 200L)
  expect_error(merge_bwg(a, b), "Conflicting chromosome lengths")

  c <- make_merge_track("sampleB", "chr1", 1L, 10L, 4, chrom_length = 200L)
  expect_error(
    merge_bwg(a, c, method = "sum"),
    "Conflicting chromosome lengths"
  )
})

test_that("strand separation prevents cross-strand arithmetic by default", {
  plus <- make_merge_track(
    "plus", "chr1", 1L, 10L, 2,
    strand = "+", chrom_length = 100L
  )
  minus <- make_merge_track(
    "minus", "chr1", 1L, 10L, 4,
    strand = "-", chrom_length = 100L
  )

  separated <- merge_bwg(
    plus,
    minus,
    method = "sum"
  )
  expect_setequal(
    separated$samples$sample_id,
    c("merged_plus", "merged_minus")
  )

  ignored <- merge_bwg(
    plus,
    minus,
    method = "sum",
    strand = "ignore"
  )
  expect_equal(ignored$samples$sample_id, "merged")
  expect_equal(ignored$data$value, 6)
})

test_that("all-zero arithmetic signal remains a valid empty track", {
  a <- make_merge_track("sampleA", "chr1", 1L, 10L, 2, chrom_length = 100L)
  b <- make_merge_track("sampleB", "chr1", 1L, 10L, -2, chrom_length = 100L)

  observed <- merge_bwg(
    a,
    b,
    method = "sum",
    drop_zero = TRUE
  )
  expect_s3_class(observed, "BwgTrack")
  expect_equal(observed$samples$sample_id, "merged")
  expect_equal(nrow(observed$data), 0L)
})

test_that("mean zero counts members without covered intervals", {
  a <- make_merge_track("sampleA", "chr1", 1L, 10L, 2, chrom_length = 100L)
  empty_samples <- data.table::data.table(
    sample_id = "sampleB",
    file = NA_character_,
    format = "memory",
    strand = "*"
  )
  empty_signal <- data.table::data.table(
    sample_id = character(),
    chrom = character(),
    start = integer(),
    end = integer(),
    value = numeric(),
    strand = character()
  )
  empty_seqinfo <- data.table::data.table(
    sample_id = "sampleB",
    chrom = "chr1",
    length = 100L
  )
  b <- bwg_track(
    empty_samples,
    empty_signal,
    seqinfo = empty_seqinfo,
    meta = list(mode = "memory")
  )

  zero_mean <- merge_bwg(
    a,
    b,
    method = "mean",
    missing = "zero"
  )
  ignore_mean <- merge_bwg(
    a,
    b,
    method = "mean",
    missing = "ignore"
  )

  expect_equal(zero_mean$data$value, 1)
  expect_equal(ignore_mean$data$value, 2)
})

test_that("merge exposes no persistence controls", {
  persistence_args <- c(
    "outdir",
    "format",
    "overwrite",
    "compress",
    "zoom",
    "max_zoom_levels"
  )
  expect_length(intersect(names(formals(merge_bwg)), persistence_args), 0L)
})

test_that("collapse_bwg is not part of the public API", {
  expect_false("collapse_bwg" %in% getNamespaceExports("bwTools"))
})
