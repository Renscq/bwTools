test_that("external compatibility assets ship without runtime dependencies", {
  compat_dir <- system.file("compatibility", package = "bwTools")
  runner <- system.file(
    "compatibility",
    "run-compatibility.001.R",
    package = "bwTools"
  )
  py_bridge <- system.file(
    "compatibility",
    "pybigwig-bridge.001.py",
    package = "bwTools"
  )
  readme <- system.file("compatibility", "README.md", package = "bwTools")
  contract <- system.file("COMPATIBILITY.md", package = "bwTools")

  expect_true(nzchar(compat_dir) && dir.exists(compat_dir))
  expect_true(nzchar(runner) && file.exists(runner))
  expect_true(nzchar(py_bridge) && file.exists(py_bridge))
  expect_true(nzchar(readme) && file.exists(readme))
  expect_true(nzchar(contract) && file.exists(contract))

  description <- utils::packageDescription("bwTools")
  runtime <- paste(
    unlist(description[c("Depends", "Imports")], use.names = FALSE),
    collapse = " "
  )
  expect_false(grepl("pyBigWig|rtracklayer|GenomicRanges|IRanges|GenomeInfoDb", runtime))
})

test_that("compatibility runner defines bidirectional reference checks", {
  runner <- system.file(
    "compatibility",
    "run-compatibility.001.R",
    package = "bwTools"
  )
  expect_true(nzchar(runner) && file.exists(runner))

  runner_text <- paste(
    readLines(runner, warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )

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
  runner <- system.file(
    "compatibility",
    "run-compatibility.001.R",
    package = "bwTools"
  )
  expect_true(nzchar(runner) && file.exists(runner))

  runner_text <- paste(
    readLines(runner, warn = FALSE, encoding = "UTF-8"),
    collapse = "\n"
  )

  expect_match(runner_text, "bwtools_example\\.bedGraph")
  expect_match(runner_text, "bwtools_example\\.chrom\\.sizes")
  expect_match(runner_text, "zoom = FALSE", fixed = TRUE)
  expect_match(runner_text, "zoom = TRUE", fixed = TRUE)
  expect_match(runner_text, "missing-signal pattern differs", fixed = TRUE)
})
