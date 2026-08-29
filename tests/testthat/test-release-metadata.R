test_that("release metadata points to canonical project resources", {
  description <- utils::packageDescription("bwTools")

  expect_match(
    description[["URL"]],
    "https://github.com/Renscq/bwTools",
    fixed = TRUE
  )
  expect_match(
    description[["URL"]],
    "https://renscq.github.io/bwTools/",
    fixed = TRUE
  )
  expect_identical(
    description[["BugReports"]],
    "https://github.com/Renscq/bwTools/issues"
  )
  expect_identical(description[["NeedsCompilation"]], "no")
})
