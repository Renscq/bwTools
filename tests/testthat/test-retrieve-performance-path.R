test_that("indexed retrieval preserves bundled BigWig signal after buffer optimization", {
  bw_file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )
  bedgraph_file <- system.file(
    "extdata", "bwtools_example.bedGraph",
    package = "bwTools", mustWork = TRUE
  )

  lazy <- read_bwg(
    bw_file,
    format = "bigwig",
    sample_ids = "example",
    mode = "lazy"
  )
  reference <- read_bwg(
    bedgraph_file,
    format = "bedgraph",
    sample_ids = "example",
    mode = "memory"
  )

  observed <- retrieve_bwg(
    lazy,
    chrom = "chr1",
    start = 450L,
    end = 3600L
  )
  expected <- retrieve_bwg(
    reference,
    chrom = "chr1",
    start = 450L,
    end = 3600L
  )

  expect_equal(observed$sample_id, expected$sample_id)
  expect_equal(observed$chrom, expected$chrom)
  expect_equal(observed$start, expected$start)
  expect_equal(observed$end, expected$end)
  expect_equal(observed$value, expected$value, tolerance = 1e-4)
  expect_equal(observed$strand, expected$strand)
})
