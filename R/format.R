# Author: Rensc
# Date: 2026-08-28
# Version: dev002
# Function: Detect supported genomic signal formats
# Input: Local file path
# Output: Standardized format name

#' Detect genomic signal file format
#'
#' @param file Local signal file path.
#' @return One of `bigwig`, `wig`, or `bedgraph`.
#' @export
detect_bwg_format <- function(file) {
  file <- bw_validate_local_file(file)
  lower <- tolower(basename(file))
  lower <- sub("\\.(gz|bgz)$", "", lower)

  if (grepl("\\.(bigwig|bw)$", lower)) {
    return("bigwig")
  }
  if (grepl("\\.wig$", lower)) {
    return("wig")
  }
  if (grepl("\\.(bedgraph|bdg)$", lower)) {
    return("bedgraph")
  }

  con <- base::file(file, open = "rb")
  on.exit(close(con), add = TRUE)
  magic <- readBin(con, what = "raw", n = 4L, size = 1L)
  if (length(magic) == 4L) {
    little <- bw_uint(magic, 0L, 4L, "little")
    big <- bw_uint(magic, 0L, 4L, "big")
    if (little == .BWT_BIGWIG_MAGIC || big == .BWT_BIGWIG_MAGIC) {
      return("bigwig")
    }
  }

  text_con <- if (grepl("\\.(gz|bgz)$", tolower(file))) {
    gzfile(file, open = "rt")
  } else {
    base::file(file, open = "rt")
  }
  on.exit(close(text_con), add = TRUE)
  lines <- readLines(text_con, n = 50L, warn = FALSE)
  lines <- trimws(lines)
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
  bw_stop(paste0("Unable to detect signal format for: ", file))
}
