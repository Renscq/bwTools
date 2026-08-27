test_that("documentation code fences use braced language tags", {
  root <- normalizePath(
    testthat::test_path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )

  files <- list.files(
    root,
    pattern = "\\.(md|qmd)$",
    recursive = TRUE,
    full.names = TRUE
  )

  bad <- character()

  for (file in files) {
    lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
    idx <- grep("^```[[:alnum:]_.+-]+[[:space:]]*$", lines)

    if (length(idx)) {
      bad <- c(
        bad,
        sprintf(
          "%s:%d:%s",
          sub(paste0("^", root, "/?"), "", file),
          idx,
          lines[idx]
        )
      )
    }
  }

  expect_length(bad, 0L, info = paste(bad, collapse = "\n"))
})
