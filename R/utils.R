# Author: Rensc
# Date: 2026-08-27
# Version: dev001
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
  if (grepl("^(https?|ftp)://", file, ignore.case = TRUE)) {
    bw_stop("Remote files are not supported by the current native R backend.")
  }
  if (grepl("^file://", file, ignore.case = TRUE)) {
    bw_stop("Use a local file path instead of a file:// URI.")
  }
  if (!file.exists(file)) {
    bw_stop(paste0("Signal file does not exist: ", file))
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

bw_file_stem <- function(file) {
  name <- basename(file)
  name <- sub("\\.(gz|bgz)$", "", name, ignore.case = TRUE)
  sub("\\.(bigwig|bw|wig|bedgraph|bdg)$", "", name, ignore.case = TRUE)
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
