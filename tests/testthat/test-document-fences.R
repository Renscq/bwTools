test_that("shipped documentation code fences use braced language tags", {
  files <- c(
    system.file("API_CONTRACT.md", package = "bwTools"),
    system.file("COMPATIBILITY.md", package = "bwTools"),
    system.file("compatibility", "README.md", package = "bwTools")
  )
  files <- files[nzchar(files)]

  expect_length(files, 3L)
  expect_true(all(file.exists(files)))

  doc_dir <- system.file("doc", package = "bwTools")
  if (nzchar(doc_dir) && dir.exists(doc_dir)) {
    files <- c(
      files,
      list.files(doc_dir, pattern = "\\.Rmd$", full.names = TRUE)
    )
  }

  bad <- character()
  for (file in files) {
    lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
    idx <- grep("^```[[:alnum:]_.+-]+[[:space:]]*$", lines)
    if (length(idx)) {
      bad <- c(
        bad,
        sprintf("%s:%d:%s", basename(file), idx, lines[idx])
      )
    }
  }

  expect_equal(bad, character())
})
