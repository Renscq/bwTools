# Author: Rensc
# Date: 2026-08-27
# Version: dev002
# Function: Provide common validation and path utilities
# Input: User arguments
# Output: Validated values

bw_stop <- function(message) {
  stop(message, call. = FALSE)
}

bw_validate_local_file <- function(file) {
  if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
    bw_stop("`file` must be a single non-missing character string.")
  }
  display_file <- bw_safe_text(file, fallback = "<path>")
  if (grepl("^(https?|ftp)://", display_file, ignore.case = TRUE)) {
    bw_stop("Remote files are not supported by the current native R backend.")
  }
  if (grepl("^file://", display_file, ignore.case = TRUE)) {
    bw_stop("Use a local file path instead of a file:// URI.")
  }
  if (!file.exists(file)) {
    bw_stop(paste0("Signal file does not exist: ", display_file))
  }
  normalizePath(file, winslash = "/", mustWork = TRUE)
}

bw_recycle_argument <- function(x, n, name, default = NULL) {
  if (is.null(x)) {
    x <- default
  }
  if (length(x) == 1L) {
    return(rep(x, n))
  }
  if (length(x) != n) {
    bw_stop(paste0("`", name, "` must have length 1 or match the number of files."))
  }
  x
}

bw_safe_text <- function(x, fallback = "") {
  x <- as.character(x)
  out <- suppressWarnings(iconv(x, from = "", to = "UTF-8", sub = "_"))
  bad <- is.na(out)
  if (any(bad)) {
    out[bad] <- fallback
  }
  out
}

bw_file_stem <- function(file) {
  display_path <- bw_safe_text(file, fallback = "sample")
  name <- basename(display_path)
  name <- sub("\\.(gz|bgz)$", "", name, ignore.case = TRUE)
  name <- sub("\\.(bigwig|bw|wig|bedgraph|bdg)$", "", name, ignore.case = TRUE)
  if (!nzchar(name)) "sample" else name
}

bw_empty_signal <- function(include_sample = FALSE) {
  if (isTRUE(include_sample)) {
    return(data.table::data.table(
      sample_id = character(),
      chrom = character(),
      start = integer(),
      end = integer(),
      value = numeric(),
      strand = character()
    ))
  }
  data.table::data.table(
    chrom = character(),
    start = integer(),
    end = integer(),
    value = numeric()
  )
}
