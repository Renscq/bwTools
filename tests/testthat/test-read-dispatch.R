test_that("automatic BigWig dispatch ignores vapply names", {
  bw_file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )

  detected <- vapply(
    bw_file,
    detect_bwg_format,
    character(1L),
    USE.NAMES = FALSE
  )

  expect_null(names(detected))
  expect_identical(detected[[1L]], "bigwig")

  x <- read_bwg(
    bw_file,
    sample_ids = "example"
  )

  expect_true(is_bwg_track(x))
  expect_identical(samples_bwg(x)$format, "bigwig")
  expect_identical(summary_bwg(x)$mode, "lazy")
})

test_that("format dispatch uses the scalar value rather than attributes", {
  fmt <- stats::setNames("bigwig", "input.bigwig")

  expect_false(identical(fmt, "bigwig"))
  expect_identical(fmt[[1L]], "bigwig")
})
