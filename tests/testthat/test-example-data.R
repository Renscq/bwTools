test_that("bundled example data are installed and discoverable", {
  bw_file <- system.file(
    "extdata",
    "bwtools_example.bigwig",
    package = "bwTools"
  )
  bedgraph_file <- system.file(
    "extdata",
    "bwtools_example.bedGraph",
    package = "bwTools"
  )
  chrom_sizes_file <- system.file(
    "extdata",
    "bwtools_example.chrom.sizes",
    package = "bwTools"
  )

  expect_true(nzchar(bw_file))
  expect_true(nzchar(bedgraph_file))
  expect_true(nzchar(chrom_sizes_file))
  expect_true(file.exists(bw_file))
  expect_true(file.exists(bedgraph_file))
  expect_true(file.exists(chrom_sizes_file))
})

test_that("bundled BigWig example has the documented chromosome lengths", {
  bw_file <- system.file(
    "extdata",
    "bwtools_example.bigwig",
    package = "bwTools",
    mustWork = TRUE
  )

  track <- read_bwg(bw_file, format = "bigwig", sample_ids = "example", mode = "lazy")
  observed <- seqinfo_bwg(track)
  expected <- data.table::data.table(
    sample_id = "example",
    chrom = c("chr1", "chr2", "chr3"),
    length = c(10000L, 8000L, 5000L)
  )

  expect_equal(observed, expected)
})

test_that("bundled BigWig and bedGraph agree for a representative region", {
  bw_file <- system.file(
    "extdata",
    "bwtools_example.bigwig",
    package = "bwTools",
    mustWork = TRUE
  )
  bedgraph_file <- system.file(
    "extdata",
    "bwtools_example.bedGraph",
    package = "bwTools",
    mustWork = TRUE
  )

  bw_track_obj <- read_bwg(
    bw_file,
    format = "bigwig",
    sample_ids = "example",
    mode = "lazy"
  )
  bedgraph_track_obj <- read_bwg(
    bedgraph_file,
    format = "bedgraph",
    sample_ids = "example",
    mode = "memory"
  )

  bw_signal <- retrieve_bwg(
    bw_track_obj,
    chrom = "chr1",
    start = 450L,
    end = 1600L
  )
  bedgraph_signal <- retrieve_bwg(
    bedgraph_track_obj,
    chrom = "chr1",
    start = 450L,
    end = 1600L
  )

  columns <- c("sample_id", "chrom", "start", "end", "value", "strand")
  bw_observed <- as.data.frame(bw_signal)[, columns, drop = FALSE]
  bedgraph_observed <- as.data.frame(bedgraph_signal)[, columns, drop = FALSE]

  expect_equal(
    bw_observed,
    bedgraph_observed,
    tolerance = 1e-4
  )
})
