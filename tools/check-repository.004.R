#!/usr/bin/env Rscript

# bwTools repository contract checker
# Script version: dev004
# Date: 2026-08-29

fail <- function(message) {
  stop(message, call. = FALSE)
}

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
read_utf8 <- function(path) {
  readLines(file.path(root, path), warn = FALSE, encoding = "UTF-8")
}

required <- c(
  "DESCRIPTION",
  "NAMESPACE",
  "README.qmd",
  "README.md",
  "_pkgdown.yml",
  ".github/workflows/R-CMD-check.yaml",
  ".github/workflows/pkgdown.yaml",
  "inst/API_CONTRACT.md",
  "inst/COMPATIBILITY.md",
  "inst/compatibility/run-compatibility.001.R",
  "inst/compatibility/pybigwig-bridge.001.py",
  "inst/compatibility/README.md",
  "vignettes/getting-started.Rmd",
  "vignettes/reading-writing.Rmd",
  "vignettes/retrieve-stats.Rmd",
  "vignettes/merge-multisample.Rmd",
  "vignettes/file-formats.Rmd",
  "vignettes/performance-compatibility.Rmd"
)
missing <- required[!file.exists(file.path(root, required))]
if (length(missing)) {
  fail(paste("Missing repository files:", paste(missing, collapse = ", ")))
}

description <- read.dcf(file.path(root, "DESCRIPTION"))
if (!identical(description[[1L, "Version"]], "0.8.10")) {
  fail("DESCRIPTION version must be 0.8.10.")
}
for (field in c("URL", "BugReports", "NeedsCompilation")) {
  if (!field %in% colnames(description)) {
    fail(paste("DESCRIPTION is missing field:", field))
  }
}
if (!grepl("https://github.com/Renscq/bwTools", description[[1L, "URL"]], fixed = TRUE) ||
    !grepl("https://renscq.github.io/bwTools/", description[[1L, "URL"]], fixed = TRUE)) {
  fail("DESCRIPTION URL field is not canonical.")
}
if (!identical(description[[1L, "BugReports"]], "https://github.com/Renscq/bwTools/issues")) {
  fail("DESCRIPTION BugReports field is not canonical.")
}
if (!identical(description[[1L, "NeedsCompilation"]], "no")) {
  fail("DESCRIPTION NeedsCompilation must remain 'no'.")
}
runtime <- paste(
  description[1L, intersect(c("Depends", "Imports"), colnames(description))],
  collapse = " "
)
if (grepl("pyBigWig|rtracklayer|GenomicRanges|IRanges|GenomeInfoDb", runtime)) {
  fail("External compatibility references must not become runtime dependencies.")
}

readme_qmd <- read_utf8("README.qmd")
readme_md <- read_utf8("README.md")

if (!any(grepl("README.md is generated from README.qmd", readme_qmd, fixed = TRUE))) {
  fail("README.qmd must state that README.md is generated from it.")
}
if (!any(grepl("README.md is generated from README.qmd", readme_md, fixed = TRUE))) {
  fail("README.md must be a rendered artifact of README.qmd.")
}
if (any(grepl("^```\\{r[},]", readme_md)) ||
    any(grepl("^```\\{\\.text\\}", readme_md))) {
  fail("README.md contains unrendered Quarto/knitr chunk syntax.")
}
if (!any(grepl("^```\\{r\\}", readme_qmd))) {
  fail("README.qmd must keep R chunks in Quarto/knitr {r} form.")
}
if (!any(grepl("^```\\{\\.text\\}", readme_qmd))) {
  fail("README.qmd must use {.text} for non-executable text blocks.")
}
if (any(grepl("Development status", readme_qmd, fixed = TRUE))) {
  fail("README must not contain a Development status section.")
}
for (entry in c(
  'pak::pak("Renscq/bwTools")',
  "https://renscq.github.io/bwTools/"
)) {
  if (!any(grepl(entry, readme_md, fixed = TRUE))) {
    fail(paste("README is missing:", entry))
  }
}

namespace <- read_utf8("NAMESPACE")
exports <- sub(
  "^export\\(([^)]+)\\)$",
  "\\1",
  grep("^export\\(", namespace, value = TRUE)
)
expected_exports <- c(
  "benchmark_bwg",
  "bwg_track",
  "detect_bwg_format",
  "is_bwg_track",
  "merge_bwg",
  "metadata_bwg",
  "read_bwg",
  "retrieve_bwg",
  "samples_bwg",
  "seqinfo_bwg",
  "summary_bwg",
  "stats_bwg",
  "validate_bwg",
  "write_bwg",
  "zoominfo_bwg"
)
if (!setequal(exports, expected_exports) || length(exports) != length(expected_exports)) {
  fail("NAMESPACE public exports do not match the frozen 15-function API.")
}

r_files <- list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE)
for (fun in expected_exports) {
  found <- FALSE
  for (file in r_files) {
    lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
    idx <- grep(paste0("^", fun, " <- function\\("), lines)
    if (!length(idx)) next
    found <- TRUE
    i <- idx[[1L]] - 1L
    block <- character()
    while (i >= 1L && startsWith(lines[[i]], "#'")) {
      block <- c(lines[[i]], block)
      i <- i - 1L
    }
    if (!any(grepl("^#' @return", block))) {
      fail(paste("Exported function is missing @return:", fun))
    }
    if (!any(grepl("^#' @examples", block))) {
      fail(paste("Exported function is missing @examples:", fun))
    }
    break
  }
  if (!found) {
    fail(paste("Could not locate exported function source:", fun))
  }
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

  rd <- file.path(root, "man", paste0(fun, ".Rd"))
  if (!file.exists(rd)) {
    fail(paste("Generated Rd topic is missing:", basename(rd)))
  }
  rd_lines <- readLines(rd, warn = FALSE, encoding = "UTF-8")
  if (!any(grepl("^\\\\keyword\\{internal\\}$", rd_lines))) {
    fail(paste("Generated S3 Rd topic must stay internal:", basename(rd)))
  }
}

compat_runner <- paste(read_utf8("inst/compatibility/run-compatibility.001.R"), collapse = "\n")
for (entry in c(
  '"pyBigWig"',
  '"rtracklayer"',
  '"UCSC"',
  '"bidirectional"',
  '"PASS"',
  '"SKIP"',
  '"FAIL"',
  "compatibility-report.tsv",
  "bwtools_example.bedGraph",
  "bwtools_example.chrom.sizes",
  "zoom = FALSE",
  "zoom = TRUE",
  "missing-signal pattern differs"
)) {
  if (!grepl(entry, compat_runner, fixed = TRUE)) {
    fail(paste("Compatibility runner is missing contract text:", entry))
  }
}

pkgdown_workflow <- read_utf8(".github/workflows/pkgdown.yaml")
required_pkgdown <- c(
  "actions/checkout@v6",
  "r-lib/actions/setup-r@v2",
  "r-lib/actions/setup-r-dependencies@v2",
  "Rscript tools/check-repository.004.R",
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
  "Rscript tools/check-repository.004.R",
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

source_docs <- list.files(
  root,
  pattern = "\\.(qmd|Rmd)$",
  recursive = TRUE,
  full.names = TRUE
)
bad <- character()
for (file in source_docs) {
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
  fail(paste(c(
    "Unbraced language code fences found in source documents:",
    bad
  ), collapse = "\\n"))
}

test_files <- list.files(
  file.path(root, "tests", "testthat"),
  pattern = "\\.R$",
  full.names = TRUE
)
for (file in test_files) {
  lines <- readLines(file, warn = FALSE, encoding = "UTF-8")
  text <- paste(lines, collapse = "\n")
  if (grepl('test_path\\("\\.\\.", "\\.\\."\\)', text)) {
    fail(paste(
      "Installed-package tests must not infer the source root via test_path:",
      basename(file)
    ))
  }
}

cat("bwTools repository contract: PASS\n")
