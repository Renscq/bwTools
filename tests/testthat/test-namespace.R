test_that("standardized public API and S3 registration are available", {
  expected_exports <- c(
    "bwg_track",
    "detect_bwg_format",
    "is_bwg_track",
    "merge_bwg",
    "metadata_bwg",
    "read_bwg",
    "retrieve_bwg",
    "samples_bwg",
    "seqinfo_bwg",
    "summary_bwg",
    "stats_bwg",
    "validate_bwg",
    "write_bwg",
    "zoominfo_bwg"
  )

  expect_setequal(getNamespaceExports("bwTools"), expected_exports)
  expect_true(is.function(getS3method("print", "bwToolsTrack", optional = TRUE)))
})
