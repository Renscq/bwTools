test_that("release metadata points to canonical project resources", {
  root <- normalizePath(
    testthat::test_path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )

  description <- read.dcf(file.path(root, "DESCRIPTION"))

  expect_match(
    description[[1L, "URL"]],
    "https://github.com/Renscq/bwTools",
    fixed = TRUE
  )
  expect_match(
    description[[1L, "URL"]],
    "https://renscq.github.io/bwTools/",
    fixed = TRUE
  )
  expect_identical(
    description[[1L, "BugReports"]],
    "https://github.com/Renscq/bwTools/issues"
  )
  expect_identical(description[[1L, "NeedsCompilation"]], "no")
})

test_that("shipped README points users to GitHub and pkgdown", {
  root <- normalizePath(
    testthat::test_path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )

  readme <- readLines(file.path(root, "README.md"), warn = FALSE)

  expect_true(any(grepl('pak::pak("Renscq/bwTools")', readme, fixed = TRUE)))
  expect_true(any(grepl("https://renscq.github.io/bwTools/", readme, fixed = TRUE)))
  expect_false(any(grepl("Development status", readme, fixed = TRUE)))
})
