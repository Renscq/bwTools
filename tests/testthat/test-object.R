test_that("bw_track creates a GeneTrackR-compatible object", {
  samples <- data.table::data.table(
    sample_id = "sampleA",
    file = "sample.bigwig",
    format = "bigwig",
    strand = "+"
  )
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    start = c(1L, 11L),
    end = c(10L, 20L),
    value = c(1, 2),
    strand = "+"
  )
  x <- bw_track(samples, signal, meta = list(mode = "memory"))
  expect_s3_class(x, "BwgTrack")
  expect_true(is_bwg_track(x))
  expect_s3_class(x$data, "data.table")
  expect_identical(x$meta$coordinate, "1-based closed")
})

test_that("memory retrieval filters and clips the requested region", {
  samples <- data.table::data.table(
    sample_id = "sampleA",
    file = "sample.bedgraph",
    format = "bedgraph",
    strand = "+"
  )
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = c("chr1", "chr1", "chr2"),
    start = c(1L, 20L, 1L),
    end = c(10L, 30L, 10L),
    value = c(1, 2, 3),
    strand = "+"
  )
  x <- bw_track(samples, signal, meta = list(mode = "memory"))
  observed <- retrieve_bwg(x, "chr1", 5L, 25L)
  expect_equal(nrow(observed), 2L)
  expect_equal(observed$start, c(5L, 20L))
  expect_equal(observed$end, c(10L, 25L))
})
