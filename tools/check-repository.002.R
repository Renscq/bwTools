#!/usr/bin/env Rscript

# bwTools repository contract checker
# Script version: dev002
# Date: 2026-08-29

fail <- function(message) {
  stop(message, call. = FALSE)
}

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
required <- c(
  "DESCRIPTION",
  "README.qmd",
  "README.md",
  "_pkgdown.yml",
  ".github/workflows/R-CMD-check.yaml",
  ".github/workflows/pkgdown.yaml"
)

missing <- required[!file.exists(file.path(root, required))]
if (length(missing)) {
  fail(paste("Missing repository files:", paste(missing, collapse = ", ")))
}

read_utf8 <- function(path) {
  readLines(file.path(root, path), warn = FALSE, encoding = "UTF-8")
}

description <- read_utf8("DESCRIPTION")
required_description <- c(
  "Version: 0.8.8",
  "URL: https://github.com/Renscq/bwTools, https://renscq.github.io/bwTools/",
  "BugReports: https://github.com/Renscq/bwTools/issues",
  "Config/Needs/website: pkgdown"
)
for (entry in required_description) {
  if (!any(description == entry)) {
    fail(paste("DESCRIPTION is missing:", entry))
  }
}

readme_qmd <- read_utf8("README.qmd")
readme_md <- read_utf8("README.md")
front_matter <- which(readme_qmd == "---")
if (length(front_matter) >= 2L && front_matter[[1L]] == 1L) {
  derived_readme <- readme_qmd[-seq_len(front_matter[[2L]])]
  while (length(derived_readme) && !nzchar(derived_readme[[1L]])) {
    derived_readme <- derived_readme[-1L]
  }
} else {
  derived_readme <- readme_qmd
}
if (!identical(readme_md, derived_readme)) {
  fail("README.md is not synchronized with the README.qmd body.")
}
if (any(grepl("Development status", readme_qmd, fixed = TRUE))) {
  fail("README must not contain a Development status section.")
}

pkgdown <- read_utf8("_pkgdown.yml")
for (entry in c(
  "url: https://renscq.github.io/bwTools/",
  "    home: https://github.com/Renscq/bwTools/"
)) {
  if (!any(pkgdown == entry)) {
    fail(paste("_pkgdown.yml is missing:", entry))
  }
}

pkgdown_workflow <- read_utf8(".github/workflows/pkgdown.yaml")
required_pkgdown <- c(
  "actions/checkout@v6",
  "r-lib/actions/setup-r@v2",
  "r-lib/actions/setup-r-dependencies@v2",
  "pkgdown::build_site_github_pages",
  "branch: gh-pages",
  "folder: docs"
)
for (entry in required_pkgdown) {
  if (!any(grepl(entry, pkgdown_workflow, fixed = TRUE))) {
    fail(paste("pkgdown workflow is missing:", entry))
  }
}

check_workflow <- read_utf8(".github/workflows/R-CMD-check.yaml")
for (entry in c(
  "macos-latest",
  "windows-latest",
  "ubuntu-latest",
  "r-lib/actions/check-r-package@v2"
)) {
  if (!any(grepl(entry, check_workflow, fixed = TRUE))) {
    fail(paste("R-CMD-check workflow is missing:", entry))
  }
}

rbuildignore <- read_utf8(".Rbuildignore")
for (entry in c(
  "^\\.github$",
  "^_pkgdown\\.yml$",
  "^docs$",
  "^DEVELOPMENT\\.md$",
  "^README\\.qmd$",
  "^tools$"
)) {
  if (!any(rbuildignore == entry)) {
    fail(paste(".Rbuildignore is missing:", entry))
  }
}

print_topics <- c(
  "R/benchmark-bwg.R" = "print.bwToolsBenchmark",
  "R/class-bwg-track.R" = "print.bwToolsTrack"
)
for (path in names(print_topics)) {
  lines <- read_utf8(path)
  fun <- print_topics[[path]]
  idx <- grep(paste0("^", fun, " <- function\\("), lines)
  if (length(idx) != 1L) {
    fail(paste("Could not locate S3 print method:", fun))
  }
  i <- idx[[1L]] - 1L
  block <- character()
  while (i >= 1L && startsWith(lines[[i]], "#'")) {
    block <- c(lines[[i]], block)
    i <- i - 1L
  }
  if (!any(grepl("^#' @keywords internal$", block))) {
    fail(paste("pkgdown-only S3 topic must be internal:", fun))
  }
}

files <- list.files(
  root,
  pattern = "\\.(md|qmd|Rmd)$",
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
if (length(bad)) {
  fail(paste(c("Unbraced language code fences found:", bad), collapse = "\n"))
}

cat("bwTools repository contract: PASS\n")
