test_that("external compatibility assets ship without runtime dependencies", {
  root <- normalizePath(
    testthat::test_path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  compat_dir <- file.path(root, "inst", "compatibility")
  runner <- file.path(compat_dir, "run-compatibility.001.R")
  py_bridge <- file.path(compat_dir, "pybigwig-bridge.001.py")
  readme <- file.path(compat_dir, "README.md")
  contract <- file.path(root, "inst", "COMPATIBILITY.md")

  expect_true(dir.exists(compat_dir))
  expect_true(file.exists(runner))
  expect_true(file.exists(py_bridge))
  expect_true(file.exists(readme))
  expect_true(file.exists(contract))

  description <- read.dcf(file.path(root, "DESCRIPTION"))
  runtime <- paste(
    description[1L, intersect(c("Depends", "Imports"), colnames(description))],
    collapse = " "
  )
  expect_false(grepl("pyBigWig|rtracklayer|GenomicRanges|IRanges|GenomeInfoDb", runtime))
})

test_that("compatibility runner defines bidirectional reference checks", {
  root <- normalizePath(
    testthat::test_path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  runner <- readLines(
    file.path(root, "inst", "compatibility", "run-compatibility.001.R"),
    warn = FALSE,
    encoding = "UTF-8"
  )
  runner_text <- paste(runner, collapse = "\n")

  expect_match(runner_text, '"pyBigWig"')
  expect_match(runner_text, '"rtracklayer"')
  expect_match(runner_text, '"UCSC"')
  expect_match(runner_text, '"bidirectional"')
  expect_match(runner_text, '"PASS"')
  expect_match(runner_text, '"SKIP"')
  expect_match(runner_text, '"FAIL"')
  expect_match(runner_text, 'compatibility-report\\.tsv')
})

test_that("compatibility validation uses installed canonical fixtures", {
  root <- normalizePath(
    testthat::test_path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
  runner <- readLines(
    file.path(root, "inst", "compatibility", "run-compatibility.001.R"),
    warn = FALSE,
    encoding = "UTF-8"
  )
  runner_text <- paste(runner, collapse = "\n")

  expect_match(runner_text, "bwtools_example\\.bedGraph")
  expect_match(runner_text, "bwtools_example\\.chrom\\.sizes")
  expect_match(runner_text, "zoom = FALSE", fixed = TRUE)
  expect_match(runner_text, "zoom = TRUE", fixed = TRUE)
  expect_match(runner_text, "missing-signal pattern differs", fixed = TRUE)
})
