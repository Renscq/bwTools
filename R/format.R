# Author: Rensc
# Date: 2026-08-28
# Version: dev004
# Function: Detect supported genomic signal formats without decoding binary files as text
# Input: Local file path
# Output: Standardized format name

bw_has_bigwig_magic <- function(x) {
  if (!is.raw(x) || length(x) < 4L) {
    return(FALSE)
  }
  bytes <- as.integer(x[seq_len(4L)])
  identical(bytes, c(38L, 252L, 143L, 136L)) ||
    identical(bytes, c(136L, 143L, 252L, 38L))
}

bw_raw_looks_binary <- function(x) {
  if (!is.raw(x) || length(x) == 0L) {
    return(FALSE)
  }
  bytes <- as.integer(x)
  controls <- bytes < 32L & !bytes %in% c(9L, 10L, 13L)
  any(bytes == 0L) || mean(controls) > 0.01
}

bw_read_format_probe <- function(file, n = 4096L) {
  con <- base::file(file, open = "rb")
  on.exit(close(con), add = TRUE)
  readBin(con, what = "raw", n = as.integer(n), size = 1L)
}

#' Detect genomic signal file format
#'
#' BigWig detection is based on its four-byte binary signature before any text
#' decoding is attempted. This keeps format detection safe for binary files and
#' for installations whose paths are not represented internally as UTF-8.
#'
#' @param file Local signal file path.
#' @return One of `bigwig`, `wig`, or `bedgraph`.
#' @details
#' BigWig is detected from its binary magic signature before text decoding.
#' WIG and bedGraph are then identified from extension and content. Invalid or
#' empty files raise a `bwTools_error`.
#'
#' @examples
#' bw_file <- system.file(
#'   "extdata", "bwtools_example.bigwig",
#'   package = "bwTools", mustWork = TRUE
#' )
#' detect_bwg_format(bw_file)
#' @seealso [read_bwg()]
#' @export
detect_bwg_format <- function(file) {
  file <- bw_validate_local_file(file)

  probe <- bw_read_format_probe(file)
  if (bw_has_bigwig_magic(probe)) {
    return("bigwig")
  }

  display_path <- bw_safe_text(file, fallback = "")
  lower <- tolower(basename(display_path))
  lower <- sub("\\.(gz|bgz)$", "", lower)

  if (grepl("\\.wig$", lower)) {
    return("wig")
  }
  if (grepl("\\.(bedgraph|bdg)$", lower)) {
    return("bedgraph")
  }

  compressed <- grepl("\\.(gz|bgz)$", tolower(display_path))
  if (!compressed && bw_raw_looks_binary(probe)) {
    bw_stop(
      "The input is a binary file but does not contain a valid BigWig signature."
    )
  }

  text_con <- if (compressed) {
    gzfile(file, open = "rt")
  } else {
    base::file(file, open = "rt")
  }
  on.exit(close(text_con), add = TRUE)

  lines <- readLines(text_con, n = 50L, warn = FALSE)
  if (length(lines) > 0L) {
    lines <- bw_safe_text(lines, fallback = "")
    lines <- trimws(lines)
  }
  lines <- lines[nzchar(lines) & !grepl("^(#|browser|track)", lines, ignore.case = TRUE)]
  if (length(lines) == 0L) {
    bw_stop("Unable to detect the signal file format from an empty file.")
  }
  if (grepl("^(fixedStep|variableStep)\\b", lines[1L], ignore.case = TRUE)) {
    return("wig")
  }
  fields <- strsplit(lines[1L], "[[:space:]]+")[[1L]]
  if (length(fields) >= 4L) {
    return("bedgraph")
  }
  bw_stop(paste0("Unable to detect signal format for: ", bw_safe_text(file, fallback = "<path>")))
}
