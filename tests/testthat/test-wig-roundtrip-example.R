test_that("bundled example survives WIG round-trip", {
  bedgraph_file <- system.file(
    "extdata",
    "bwtools_example.bedGraph",
    package = "bwTools",
    mustWork = TRUE
  )

  source <- read_bwg(
    bedgraph_file,
    format = "bedgraph",
    sample_names = "example",
    mode = "memory"
  )

  outdir <- tempfile(pattern = "bwtools_wig_example_")
  dir.create(outdir)

  written <- write_bwg(
    source,
    outdir = outdir,
    format = "wig",
    overwrite = TRUE
  )

  observed <- read_bwg(
    written$file,
    format = "wig",
    sample_names = "example",
    mode = "memory"
  )

  columns <- c("sample_id", "chrom", "start", "end", "value", "strand")
  expected <- as.data.frame(source$data)[, columns, drop = FALSE]
  actual <- as.data.frame(observed$data)[, columns, drop = FALSE]

  expected <- expected[
    order(expected$sample_id, expected$chrom, expected$start, expected$end),
    ,
    drop = FALSE
  ]
  actual <- actual[
    order(actual$sample_id, actual$chrom, actual$start, actual$end),
    ,
    drop = FALSE
  ]
  row.names(expected) <- NULL
  row.names(actual) <- NULL

  expect_equal(nrow(actual), 296L)
  expect_equal(actual$sample_id, expected$sample_id)
  expect_equal(actual$chrom, expected$chrom)
  expect_equal(actual$start, expected$start)
  expect_equal(actual$end, expected$end)
  expect_equal(actual$value, expected$value, tolerance = 1e-4)
  expect_equal(actual$strand, expected$strand)
})
