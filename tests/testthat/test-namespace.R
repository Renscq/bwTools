test_that("public API and S3 registration are available", {
  expected_exports <- c(
    "bw_detect_format",
    "bw_track",
    "collapse_bwg",
    "merge_bwg",
    "is_bwg_track",
    "read_bwg",
    "retrieve_bwg",
    "seqinfo_bwg",
    "summary_bwg",
    "stats_bwg",
    "write_bwg",
    "zoominfo_bwg"
  )

  expect_true(all(expected_exports %in% getNamespaceExports("bwTools")))
  expect_true(is.function(getS3method("print", "bwToolsTrack", optional = TRUE)))
})
