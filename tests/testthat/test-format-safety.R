test_that("BigWig magic is detected before text decoding", {
  source <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )
  renamed <- tempfile(fileext = ".bin")
  expect_true(file.copy(source, renamed, overwrite = TRUE))

  expect_identical(detect_bwg_format(renamed), "bigwig")

  x <- read_bwg(
    renamed,
    sample_ids = "example",
    mode = "lazy"
  )
  expect_true(is_bwg_track(x))
  expect_identical(samples_bwg(x)$format, "bigwig")
})

test_that("unknown binary input never falls through to text parsing", {
  file <- tempfile(fileext = ".bin")
  con <- base::file(file, open = "wb")
  writeBin(
    as.raw(c(0x01, 0x00, 0xFF, 0x80, rep(0x00, 32))),
    con,
    useBytes = TRUE
  )
  close(con)

  expect_error(
    detect_bwg_format(file),
    "binary file.*valid BigWig signature"
  )
})

test_that("automatic bundled BigWig loading uses binary-safe detection", {
  bw_file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )

  x <- read_bwg(
    bw_file,
    sample_ids = "example"
  )

  expect_true(is_bwg_track(x))
  expect_identical(summary_bwg(x)$mode, "lazy")
  expect_identical(samples_bwg(x)$format, "bigwig")
})
