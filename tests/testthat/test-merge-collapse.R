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
  bw_track(samples, signal, seqinfo = seqinfo, meta = list(mode = "memory"))
}

test_that("same-sample disjoint regions merge directly", {
  a <- make_merge_track("sampleA", "chr1", c(1L, 20L), c(10L, 30L), c(1, 2), chrom_length = 100L)
  b <- make_merge_track("sampleA", "chr2", c(5L, 40L), c(15L, 50L), c(3, 4), chrom_length = 80L)

  observed <- merge_bwg(a, b)

  expect_s3_class(observed, "BwgTrack")
  expect_equal(observed$samples$sample_id, "sampleA")
  expect_equal(nrow(observed$data), 4L)
  expect_equal(observed$data$chrom, c("chr1", "chr1", "chr2", "chr2"))
  expect_equal(observed$meta$operation, "merge")
  expect_equal(observed$meta$input_intervals, 4L)
  expect_equal(observed$meta$output_intervals, 4L)
})

test_that("same-sample overlapping regions collapse by sum or mean", {
  a <- make_merge_track("sampleA", "chr1", 1L, 10L, 2, chrom_length = 100L)
  b <- make_merge_track("sampleA", "chr1", 6L, 15L, 4, chrom_length = 100L)

  summed <- collapse_bwg(a, b, method = "sum", by = "sample")
  expect_equal(summed$data$start, c(1L, 6L, 11L))
  expect_equal(summed$data$end, c(5L, 10L, 15L))
  expect_equal(summed$data$value, c(2, 6, 4))

  mean_zero <- collapse_bwg(
    a,
    b,
    method = "mean",
    by = "sample",
    missing = "zero"
  )
  expect_equal(mean_zero$data$value, c(1, 3, 2))

  mean_ignore <- collapse_bwg(
    a,
    b,
    method = "mean",
    by = "sample",
    missing = "ignore"
  )
  expect_equal(mean_ignore$data$value, c(2, 3, 4))
})

test_that("different samples with independent regions merge without arithmetic", {
  a <- make_merge_track("sampleA", "chr1", 1L, 10L, 2, chrom_length = 100L)
  b <- make_merge_track("sampleB", "chr2", 20L, 30L, 4, chrom_length = 80L)

  observed <- merge_bwg(a, b)

  expect_equal(observed$samples$sample_id, c("sampleA", "sampleB"))
  expect_equal(nrow(observed$data), 2L)
  expect_equal(observed$data$value, c(2, 4))
})

test_that("different samples can collapse into one output sample", {
  a <- make_merge_track("sampleA", "chr1", 1L, 10L, 2, chrom_length = 100L)
  b <- make_merge_track("sampleB", "chr1", 1L, 10L, 4, chrom_length = 100L)

  summed <- collapse_bwg(
    a,
    b,
    method = "sum",
    by = "all",
    output_sample = "combined"
  )
  expect_equal(summed$samples$sample_id, "combined")
  expect_equal(summed$data$value, 6)

  averaged <- collapse_bwg(
    a,
    b,
    method = "mean",
    by = "all",
    output_sample = "combined",
    missing = "zero"
  )
  expect_equal(averaged$data$value, 3)
})

test_that("explicit groups support replicate-level aggregation", {
  a <- make_merge_track("A1", "chr1", 1L, 10L, 2, chrom_length = 100L)
  b <- make_merge_track("A2", "chr1", 1L, 10L, 4, chrom_length = 100L)
  c <- make_merge_track("B1", "chr1", 1L, 10L, 10, chrom_length = 100L)
  groups <- data.frame(
    sample_id = c("A1", "A2", "B1"),
    group_id = c("groupA", "groupA", "groupB")
  )

  observed <- collapse_bwg(
    a,
    b,
    c,
    method = "mean",
    groups = groups,
    missing = "zero"
  )

  expect_equal(observed$samples$sample_id, c("groupA", "groupB"))
  expect_equal(observed$data$value, c(3, 10))
})

test_that("direct merge rejects overlaps and handles duplicates explicitly", {
  a <- make_merge_track("sampleA", "chr1", 1L, 10L, 2, chrom_length = 100L)
  b <- make_merge_track("sampleA", "chr1", 6L, 15L, 4, chrom_length = 100L)
  expect_error(
    merge_bwg(a, b),
    "Use `collapse_bwg\\(\\)`"
  )

  duplicate <- make_merge_track("sampleA", "chr1", 1L, 10L, 2, chrom_length = 100L)
  expect_error(merge_bwg(a, duplicate), "Exact duplicate")
  dropped <- merge_bwg(a, duplicate, duplicates = "drop")
  expect_equal(nrow(dropped$data), 1L)
})

test_that("chromosome length conflicts are rejected", {
  a <- make_merge_track("sampleA", "chr1", 1L, 10L, 2, chrom_length = 100L)
  b <- make_merge_track("sampleA", "chr1", 20L, 30L, 4, chrom_length = 200L)
  expect_error(merge_bwg(a, b), "Conflicting chromosome lengths")

  c <- make_merge_track("sampleB", "chr1", 1L, 10L, 4, chrom_length = 200L)
  expect_error(
    collapse_bwg(a, c, method = "sum", by = "all"),
    "Conflicting chromosome lengths"
  )
})

test_that("strand separation prevents cross-strand arithmetic by default", {
  plus <- make_merge_track("plus", "chr1", 1L, 10L, 2, strand = "+", chrom_length = 100L)
  minus <- make_merge_track("minus", "chr1", 1L, 10L, 4, strand = "-", chrom_length = 100L)

  separated <- collapse_bwg(
    plus,
    minus,
    method = "sum",
    by = "all",
    output_sample = "combined"
  )
  expect_setequal(separated$samples$sample_id, c("combined_plus", "combined_minus"))
  plus_idx <- which(separated$data$sample_id == "combined_plus")
  minus_idx <- which(separated$data$sample_id == "combined_minus")
  expect_equal(separated$data$value[plus_idx], 2)
  expect_equal(separated$data$value[minus_idx], 4)

  ignored <- collapse_bwg(
    plus,
    minus,
    method = "sum",
    by = "all",
    output_sample = "combined",
    strand = "ignore"
  )
  expect_equal(ignored$samples$sample_id, "combined")
  expect_equal(ignored$data$value, 6)
})

test_that("all-zero collapsed signal remains a valid empty track", {
  a <- make_merge_track("sampleA", "chr1", 1L, 10L, 2, chrom_length = 100L)
  b <- make_merge_track("sampleB", "chr1", 1L, 10L, -2, chrom_length = 100L)

  observed <- collapse_bwg(
    a,
    b,
    method = "sum",
    by = "all",
    output_sample = "combined",
    drop_zero = TRUE
  )
  expect_s3_class(observed, "BwgTrack")
  expect_equal(observed$samples$sample_id, "combined")
  expect_equal(nrow(observed$data), 0L)
})

test_that("mean zero counts group members without covered intervals", {
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
  b <- bw_track(
    empty_samples,
    empty_signal,
    seqinfo = empty_seqinfo,
    meta = list(mode = "memory")
  )

  zero_mean <- collapse_bwg(
    a,
    b,
    method = "mean",
    by = "all",
    output_sample = "combined",
    missing = "zero"
  )
  ignore_mean <- collapse_bwg(
    a,
    b,
    method = "mean",
    by = "all",
    output_sample = "combined",
    missing = "ignore"
  )

  expect_equal(zero_mean$data$value, 1)
  expect_equal(ignore_mean$data$value, 2)
})

test_that("merge and collapse expose no persistence controls", {
  persistence_args <- c(
    "outdir",
    "format",
    "overwrite",
    "compress",
    "zoom",
    "max_zoom_levels"
  )
  expect_length(intersect(names(formals(merge_bwg)), persistence_args), 0L)
  expect_length(intersect(names(formals(collapse_bwg)), persistence_args), 0L)
})
