test_that("WIG writing keeps chromosome groups isolated", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = c("chr1", "chr1", "chr2"),
    start = c(10L, 30L, 5L),
    end = c(12L, 31L, 7L),
    value = c(1.5, -2.25, 3.75),
    strand = "*"
  )
  track <- bw_track(
    samples = data.table::data.table(sample_id = "sampleA", strand = "*"),
    data = signal,
    meta = list(mode = "memory")
  )

  outdir <- tempfile(pattern = "bwtools_wig_scope_")
  dir.create(outdir)
  written <- write_bwg(
    track,
    outdir = outdir,
    format = "wig",
    overwrite = TRUE
  )

  lines <- readLines(written$file, warn = FALSE)
  expect_equal(sum(startsWith(lines, "variableStep chrom=chr1 ")), 2L)
  expect_equal(sum(startsWith(lines, "variableStep chrom=chr2 ")), 1L)

  observed <- read_bwg(
    written$file,
    format = "wig",
    sample_names = "sampleA",
    mode = "memory"
  )$data

  expect_equal(nrow(observed), 3L)
  expect_equal(observed$chrom, c("chr1", "chr1", "chr2"))
  expect_equal(observed$start, c(10L, 30L, 5L))
})

test_that("sample-specific seqinfo selects only the requested sample", {
  object <- list(
    seqinfo = data.table::data.table(
      sample_id = c("sampleA", "sampleB"),
      chrom = c("chr1", "chr2"),
      length = c(100L, 200L)
    )
  )

  observed <- bwTools:::bw_resolve_chrom_sizes(
    object,
    chrom_sizes = NULL,
    sample_id = "sampleB"
  )

  expect_equal(observed$chrom, "chr2")
  expect_equal(observed$length, 200)
})

test_that("memory statistics isolate the requested chromosome", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = c("chr1", "chr2"),
    start = c(1L, 1L),
    end = c(10L, 10L),
    value = c(1, 100),
    strand = "*"
  )
  track <- bw_track(
    samples = data.table::data.table(sample_id = "sampleA", strand = "*"),
    data = signal,
    meta = list(mode = "memory")
  )

  observed <- stats_bwg(
    track,
    chrom = "chr1",
    start = 1L,
    end = 10L,
    n_bins = 1L,
    stat = "mean",
    use_zoom = FALSE
  )

  expect_equal(observed$value, 1)
  expect_equal(observed$covered_bases, 10)
})
