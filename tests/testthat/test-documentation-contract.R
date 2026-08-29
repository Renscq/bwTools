test_that("installed user documentation structure is stable", {
  doc_dir <- system.file("doc", package = "bwTools")
  expect_true(nzchar(doc_dir) && dir.exists(doc_dir))

  stems <- c(
    "getting-started",
    "reading-writing",
    "retrieve-stats",
    "merge-multisample",
    "file-formats",
    "performance-compatibility"
  )
  expected_html <- file.path(doc_dir, paste0(stems, ".html"))

  expect_true(
    all(file.exists(expected_html)),
    info = paste("Missing installed vignette(s):", paste(
      basename(expected_html[!file.exists(expected_html)]),
      collapse = ", "
    ))
  )
})

test_that("all exported functions have installed help topics", {
  exports <- getNamespaceExports("bwTools")
  missing <- exports[!vapply(
    exports,
    function(fun) length(utils::help(fun, package = "bwTools")) > 0L,
    logical(1L)
  )]

  expect_equal(missing, character())
})

test_that("internal S3 print methods remain registered but unexported", {
  exports <- getNamespaceExports("bwTools")

  expect_true(is.function(utils::getS3method(
    "print",
    "bwToolsBenchmark",
    optional = TRUE
  )))
  expect_true(is.function(utils::getS3method(
    "print",
    "bwToolsTrack",
    optional = TRUE
  )))
  expect_false("print.bwToolsBenchmark" %in% exports)
  expect_false("print.bwToolsTrack" %in% exports)
})
