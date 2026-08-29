test_that("installed user documentation structure is stable", {
  root <- normalizePath(
    testthat::test_path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )

  expected_vignettes <- file.path(
    root,
    "vignettes",
    c(
      "getting-started.Rmd",
      "reading-writing.Rmd",
      "retrieve-stats.Rmd",
      "merge-multisample.Rmd",
      "file-formats.Rmd",
      "performance-compatibility.Rmd"
    )
  )

  expect_true(all(file.exists(expected_vignettes)))

  description <- readLines(file.path(root, "DESCRIPTION"), warn = FALSE)
  expect_true(any(grepl("^VignetteBuilder: knitr$", description)))
  expect_true(any(grepl("knitr", description, fixed = TRUE)))
  expect_true(any(grepl("rmarkdown", description, fixed = TRUE)))

  readme <- readLines(file.path(root, "README.md"), warn = FALSE)
  expect_false(any(grepl("Development status", readme, fixed = TRUE)))
  expect_true(any(grepl("^## Quick start$", readme)))
  expect_true(any(grepl("^## Documentation$", readme)))
  expect_true(any(grepl("^## Public API$", readme)))

  for (file in expected_vignettes) {
    lines <- readLines(file, warn = FALSE)
    expect_true(any(grepl("VignetteIndexEntry", lines, fixed = TRUE)))
    expect_true(any(grepl("VignetteEngine", lines, fixed = TRUE)))
    expect_true(any(grepl("VignetteEncoding", lines, fixed = TRUE)))
  }
})

test_that("all exported functions remain discoverable in shipped documentation", {
  root <- normalizePath(
    testthat::test_path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )

  namespace <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
  exports <- sub("^export\\(([^)]+)\\)$", "\\1", grep("^export\\(", namespace, value = TRUE))

  docs <- c(
    readLines(file.path(root, "README.md"), warn = FALSE),
    unlist(lapply(
      list.files(file.path(root, "vignettes"), pattern = "\\.Rmd$", full.names = TRUE),
      readLines,
      warn = FALSE
    ), use.names = FALSE)
  )

  missing <- exports[!vapply(
    exports,
    function(fun) any(grepl(paste0("`", fun, "\\(\\)`"), docs)),
    logical(1L)
  )]
  expect_equal(missing, character())
})

test_that("every exported function has return and example documentation", {
  root <- normalizePath(
    testthat::test_path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )

  namespace <- readLines(file.path(root, "NAMESPACE"), warn = FALSE)
  exports <- sub("^export\\(([^)]+)\\)$", "\\1", grep("^export\\(", namespace, value = TRUE))
  r_files <- list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)

  missing_return <- character()
  missing_examples <- character()

  for (fun in exports) {
    found <- FALSE
    for (file in r_files) {
      lines <- readLines(file, warn = FALSE)
      idx <- grep(paste0("^", fun, " <- function\\("), lines)
      if (!length(idx)) next
      found <- TRUE
      i <- idx[[1L]] - 1L
      while (i >= 1L && !startsWith(lines[[i]], "#'")) i <- i - 1L
      block <- character()
      while (i >= 1L && startsWith(lines[[i]], "#'")) {
        block <- c(lines[[i]], block)
        i <- i - 1L
      }
      if (!any(grepl("^#' @return", block))) missing_return <- c(missing_return, fun)
      if (!any(grepl("^#' @examples", block))) missing_examples <- c(missing_examples, fun)
      break
    }
    expect_true(found, info = paste("Exported function source not found:", fun))
  }

  expect_equal(missing_return, character())
  expect_equal(missing_examples, character())
})

test_that("S3 print help topics stay internal to pkgdown reference", {
  root <- normalizePath(
    testthat::test_path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )

  topics <- c(
    "print.bwToolsBenchmark.Rd",
    "print.bwToolsTrack.Rd"
  )

  for (topic in topics) {
    lines <- readLines(file.path(root, "man", topic), warn = FALSE)
    expect_true(
      any(grepl("^\\\\keyword\\{internal\\}$", lines)),
      info = paste("Missing internal keyword in", topic)
    )
  }
})

