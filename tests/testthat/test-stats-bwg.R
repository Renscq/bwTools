test_that("stats_bwg calculates exact coverage-aware statistics", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    start = c(1L, 6L),
    end = c(5L, 10L),
    value = c(2, 4),
    strand = "*"
  )
  track <- bw_track(
    data.table::data.table(sample_id = "sampleA", strand = "*"),
    signal,
    seqinfo = data.table::data.table(sample_id = "sampleA", chrom = "chr1", length = 20L),
    meta = list(mode = "memory")
  )

  mean_result <- stats_bwg(track, "chr1", 1L, 10L, stat = "mean")
  expect_equal(mean_result$value, 3)
  expect_equal(mean_result$covered_bases, 10)
  expect_identical(mean_result$source, "full")

  sum_result <- stats_bwg(track, "chr1", 1L, 10L, stat = "sum")
  expect_equal(sum_result$value, 30)

  coverage_result <- stats_bwg(track, "chr1", 1L, 20L, stat = "coverage")
  expect_equal(coverage_result$value, 0.5)

  binned <- stats_bwg(track, "chr1", 1L, 10L, n_bins = 2L, stat = "mean")
  expect_equal(binned$start, c(1L, 6L))
  expect_equal(binned$end, c(5L, 10L))
  expect_equal(binned$value, c(2, 4))
})

test_that("stats_bwg reports missing bins as NA rather than zero", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    start = 1L,
    end = 5L,
    value = 2,
    strand = "*"
  )
  track <- bw_track(
    data.table::data.table(sample_id = "sampleA", strand = "*"),
    signal,
    meta = list(mode = "memory")
  )
  result <- stats_bwg(track, "chr1", 1L, 10L, n_bins = 2L, stat = "coverage")
  expect_equal(result$value[1L], 1)
  expect_true(is.na(result$value[2L]))
  expect_equal(result$covered_bases, c(5, 0))
})

test_that("stats_bwg can force full-resolution BigWig calculations", {
  bw_file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )
  lazy <- read_bwg(
    bw_file,
    format = "bigwig",
    sample_names = "example",
    mode = "lazy"
  )
  result <- stats_bwg(
    lazy,
    chrom = "chr1",
    start = 1L,
    end = 10000L,
    n_bins = 10L,
    stat = "mean",
    use_zoom = FALSE
  )
  expect_equal(nrow(result), 10L)
  expect_true(all(result$source == "full"))
  expect_true(all(is.na(result$resolution)))
})

test_that("stats_bwg uses zoom levels from native writer for coarse bins", {
  bedgraph_file <- system.file(
    "extdata", "bwtools_example.bedGraph",
    package = "bwTools", mustWork = TRUE
  )
  chrom_sizes_file <- system.file(
    "extdata", "bwtools_example.chrom.sizes",
    package = "bwTools", mustWork = TRUE
  )
  source <- read_bwg(
    bedgraph_file,
    format = "bedgraph",
    sample_names = "example",
    mode = "memory"
  )
  outdir <- tempfile(pattern = "bwtools_stats_zoom_")
  dir.create(outdir)
  written <- write_bwg(
    source,
    outdir = outdir,
    format = "bigwig",
    chrom_sizes = chrom_sizes_file,
    overwrite = TRUE,
    zoom = TRUE
  )
  lazy <- read_bwg(
    written$file,
    format = "bigwig",
    sample_names = "example",
    mode = "lazy"
  )
  result <- stats_bwg(
    lazy,
    chrom = "chr1",
    start = 1L,
    end = 10000L,
    n_bins = 1L,
    stat = "mean",
    use_zoom = TRUE
  )
  expect_identical(result$source, "zoom")
  expect_true(is.finite(result$resolution))

  exact <- stats_bwg(
    lazy,
    chrom = "chr1",
    start = 1L,
    end = 10000L,
    n_bins = 1L,
    stat = "mean",
    use_zoom = FALSE
  )
  expect_equal(result$value, exact$value, tolerance = 1e-4)
})
