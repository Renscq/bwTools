test_that("bedGraph is converted from 0-based half-open to 1-based closed", {
  file <- tempfile(fileext = ".bedgraph")
  writeLines(c("chr1\t0\t10\t1.5", "chr1\t10\t11\t2.5"), file)
  x <- read_bwg(file, format = "bedgraph", sample_names = "sampleA", mode = "memory")
  expect_equal(x$data$start, c(1L, 11L))
  expect_equal(x$data$end, c(10L, 11L))
  expect_equal(x$data$value, c(1.5, 2.5))
})

test_that("variableStep WIG preserves span", {
  file <- tempfile(fileext = ".wig")
  writeLines(c(
    "variableStep chrom=chr1 span=3",
    "10 1.5",
    "20 2.5"
  ), file)
  x <- read_bwg(file, format = "wig", sample_names = "sampleA", mode = "memory")
  expect_equal(x$data$start, c(10L, 20L))
  expect_equal(x$data$end, c(12L, 22L))
})

test_that("fixedStep WIG uses start, step, and span", {
  file <- tempfile(fileext = ".wig")
  writeLines(c(
    "fixedStep chrom=chr1 start=5 step=3 span=2",
    "1",
    "2",
    "3"
  ), file)
  x <- read_bwg(file, format = "wig", sample_names = "sampleA", mode = "memory")
  expect_equal(x$data$start, c(5L, 8L, 11L))
  expect_equal(x$data$end, c(6L, 9L, 12L))
})
