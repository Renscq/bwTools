test_that("native BigWig reader returns chromosome metadata", {
  file <- test_path("fixtures", "bwtools_reader_test.bigwig")
  track <- read_bwg(file, format = "bigwig", sample_ids = "sampleA", mode = "lazy")
  observed <- seqinfo_bwg(track)
  expected <- data.table::data.table(
    sample_id = "sampleA",
    chrom = c("chr1", "chr2"),
    length = c(1000L, 500L)
  )
  expect_s3_class(observed, "data.table")
  expect_equal(observed, expected)
})

test_that("native BigWig reader returns 1-based closed intervals", {
  file <- test_path("fixtures", "bwtools_reader_test.bigwig")
  track <- read_bwg(
    file,
    format = "bigwig",
    sample_ids = "sampleA",
    mode = "lazy"
  )
  observed <- retrieve_bwg(track, "chr1", 1L, 200L)
  expect_equal(observed$start, c(10L, 20L, 100L))
  expect_equal(observed$end, c(12L, 20L, 110L))
  expect_equal(observed$value, c(1.5, 2.5, -3.25), tolerance = 1e-6)
})

test_that("native BigWig interval queries clip boundaries", {
  file <- test_path("fixtures", "bwtools_reader_test.bigwig")
  track <- read_bwg(file, format = "bigwig", mode = "lazy")
  observed <- retrieve_bwg(track, "chr1", 11L, 105L)
  expect_equal(observed$start, c(11L, 20L, 100L))
  expect_equal(observed$end, c(12L, 20L, 105L))
})
