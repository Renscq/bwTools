# Author: Rensc
# Date: 2026-08-29
# Version: dev006
# Function: Define and access the standardized BwgTrack object contract
# Input: Sample metadata, signal data, sequence metadata, and package metadata
# Output: Standardized BwgTrack-compatible S3 object

bw_normalize_samples_table <- function(samples, data = NULL) {
  samples <- data.table::copy(data.table::as.data.table(samples))
  if (!"sample_id" %in% names(samples)) {
    bw_stop("`samples` must contain a `sample_id` column.")
  }
  if (nrow(samples) < 1L) {
    bw_stop("`samples` must contain at least one sample.")
  }
  data.table::set(samples, j = "sample_id", value = as.character(samples[["sample_id"]]))
  if (anyNA(samples[["sample_id"]]) || any(!nzchar(samples[["sample_id"]]))) {
    bw_stop("`sample_id` values must be non-missing and non-empty.")
  }
  if (anyDuplicated(samples[["sample_id"]])) {
    bw_stop("`sample_id` values must be unique.")
  }

  if (!"file" %in% names(samples)) {
    data.table::set(samples, j = "file", value = NA_character_)
  }
  data.table::set(samples, j = "file", value = as.character(samples[["file"]]))

  if (!"format" %in% names(samples)) {
    default_format <- if (is.null(data)) NA_character_ else "memory"
    data.table::set(samples, j = "format", value = default_format)
  }
  data.table::set(samples, j = "format", value = tolower(as.character(samples[["format"]])))

  if (!"strand" %in% names(samples)) {
    data.table::set(samples, j = "strand", value = "*")
  }
  data.table::set(samples, j = "strand", value = as.character(samples[["strand"]]))
  if (anyNA(samples[["strand"]]) || any(!samples[["strand"]] %in% c("+", "-", "*"))) {
    bw_stop("`samples$strand` values must be '+', '-', or '*'.")
  }

  data.table::set(
    samples,
    j = "has_strand",
    value = samples[["strand"]] %in% c("+", "-")
  )
  core <- c("sample_id", "file", "format", "strand", "has_strand")
  data.table::setcolorder(samples, c(core, setdiff(names(samples), core)))
  samples[]
}

bw_normalize_signal_table <- function(data, sample_ids) {
  if (is.null(data)) return(NULL)
  data <- data.table::copy(data.table::as.data.table(data))
  required <- c("sample_id", "chrom", "start", "end", "value")
  missing <- setdiff(required, names(data))
  if (length(missing) > 0L) {
    bw_stop(paste0("Signal data are missing required columns: ", paste(missing, collapse = ", "), "."))
  }
  if (!"strand" %in% names(data)) {
    data.table::set(data, j = "strand", value = "*")
  }

  data.table::set(data, j = "sample_id", value = as.character(data[["sample_id"]]))
  data.table::set(data, j = "chrom", value = as.character(data[["chrom"]]))
  start_value <- suppressWarnings(as.numeric(data[["start"]]))
  end_value <- suppressWarnings(as.numeric(data[["end"]]))
  value <- suppressWarnings(as.numeric(data[["value"]]))
  strand <- as.character(data[["strand"]])

  if (anyNA(data[["sample_id"]]) || any(!nzchar(data[["sample_id"]]))) {
    bw_stop("Signal data contain missing or empty sample IDs.")
  }
  unknown <- setdiff(unique(data[["sample_id"]]), sample_ids)
  if (length(unknown) > 0L) {
    bw_stop(paste0("Signal data contain unknown sample IDs: ", paste(unknown, collapse = ", "), "."))
  }
  if (anyNA(data[["chrom"]]) || any(!nzchar(data[["chrom"]]))) {
    bw_stop("Signal data contain missing or empty chromosome names.")
  }
  invalid <- !is.finite(start_value) | !is.finite(end_value) |
    start_value < 1 | end_value < start_value |
    start_value != floor(start_value) | end_value != floor(end_value) |
    start_value > .Machine$integer.max | end_value > .Machine$integer.max
  if (any(invalid)) {
    bw_stop("Signal data contain invalid 1-based closed intervals.")
  }
  if (any(!is.finite(value))) {
    bw_stop("Signal data contain non-finite values.")
  }
  if (anyNA(strand) || any(!strand %in% c("+", "-", "*"))) {
    bw_stop("Signal data contain invalid strand values; use '+', '-', or '*'.")
  }

  data.table::set(data, j = "start", value = as.integer(start_value))
  data.table::set(data, j = "end", value = as.integer(end_value))
  data.table::set(data, j = "value", value = value)
  data.table::set(data, j = "strand", value = strand)
  core <- c("sample_id", "chrom", "start", "end", "value", "strand")
  data.table::setcolorder(data, c(core, setdiff(names(data), core)))
  data.table::setorderv(data, c("sample_id", "chrom", "strand", "start", "end"))
  data.table::setindexv(data, c("sample_id", "chrom", "start", "end"))
  data[]
}

bw_normalize_seqinfo_table <- function(seqinfo, sample_ids) {
  if (is.null(seqinfo)) return(NULL)
  seqinfo <- data.table::copy(data.table::as.data.table(seqinfo))
  required <- c("chrom", "length")
  missing <- setdiff(required, names(seqinfo))
  if (length(missing) > 0L) {
    bw_stop(paste0("`seqinfo` is missing required columns: ", paste(missing, collapse = ", "), "."))
  }

  data.table::set(seqinfo, j = "chrom", value = as.character(seqinfo[["chrom"]]))
  if (anyNA(seqinfo[["chrom"]]) || any(!nzchar(seqinfo[["chrom"]]))) {
    bw_stop("`seqinfo$chrom` values must be non-missing and non-empty.")
  }
  length_value <- suppressWarnings(as.numeric(seqinfo[["length"]]))
  invalid <- !is.na(length_value) & (
    !is.finite(length_value) | length_value < 1 |
      length_value != floor(length_value) |
      length_value > .Machine$integer.max
  )
  if (any(invalid)) {
    bw_stop("`seqinfo$length` must contain positive integer lengths or NA.")
  }
  data.table::set(seqinfo, j = "length", value = as.integer(length_value))

  if (!"sample_id" %in% names(seqinfo)) {
    expanded <- lapply(sample_ids, function(sample_id) {
      x <- data.table::copy(seqinfo)
      data.table::set(x, j = "sample_id", value = sample_id)
      x
    })
    seqinfo <- data.table::rbindlist(expanded, use.names = TRUE, fill = TRUE)
  } else {
    data.table::set(seqinfo, j = "sample_id", value = as.character(seqinfo[["sample_id"]]))
  }

  unknown <- setdiff(unique(seqinfo[["sample_id"]]), sample_ids)
  if (length(unknown) > 0L) {
    bw_stop(paste0("`seqinfo` contains unknown sample IDs: ", paste(unknown, collapse = ", "), "."))
  }

  core <- c("sample_id", "chrom", "length")
  data.table::setcolorder(seqinfo, c(core, setdiff(names(seqinfo), core)))
  data.table::setorderv(seqinfo, c("sample_id", "chrom"))
  seqinfo <- unique(seqinfo)

  key <- paste(seqinfo[["sample_id"]], seqinfo[["chrom"]], sep = "\r")
  for (k in unique(key)) {
    idx <- which(key == k)
    lengths <- unique(seqinfo[["length"]][idx])
    lengths <- lengths[!is.na(lengths)]
    if (length(lengths) > 1L) {
      bw_stop("`seqinfo` contains conflicting lengths for a sample/chromosome pair.")
    }
  }
  seqinfo[]
}

#' Construct a standardized BwgTrack object
#'
#' @param samples Sample metadata table containing `sample_id`. Missing standard
#'   columns are normalized automatically.
#' @param data Optional in-memory signal table using 1-based closed coordinates.
#' @param seqinfo Optional chromosome metadata table containing `chrom` and
#'   `length`. When `sample_id` is absent, rows are replicated across samples.
#'   When supplied, every signal sample/chromosome pair must be represented and
#'   signal intervals may not exceed a known chromosome length.
#' @param meta Optional metadata list.
#' @return A validated object inheriting from `BwgTrack` and `bwToolsTrack`.
#' @details
#' Signal coordinates supplied in `data` must already be 1-based closed.
#' Construction normalizes standard columns, records schema version 2, and
#' validates cross-slot sample/chromosome consistency before returning.
#'
#' @examples
#' x <- bwg_track(
#'   samples = data.frame(sample_id = "sampleA"),
#'   data = data.frame(
#'     sample_id = "sampleA", chrom = "chr1",
#'     start = 1L, end = 10L, value = 2, strand = "+"
#'   ),
#'   seqinfo = data.frame(chrom = "chr1", length = 100L)
#' )
#' x
#' @seealso [validate_bwg()], [read_bwg()]
#' @export
bwg_track <- function(samples, data = NULL, seqinfo = NULL, meta = list()) {
  if (!is.list(meta)) {
    bw_stop("`meta` must be a list.")
  }
  samples <- bw_normalize_samples_table(samples, data = data)
  data <- bw_normalize_signal_table(data, samples[["sample_id"]])
  seqinfo <- bw_normalize_seqinfo_table(seqinfo, samples[["sample_id"]])

  meta$coordinate <- "1-based closed"
  meta$schema_version <- "2"
  meta$backend <- "bwTools-native-R"
  meta$mode <- if (is.null(data)) "lazy" else "memory"

  out <- list(
    samples = samples[],
    data = if (is.null(data)) NULL else data[],
    seqinfo = if (is.null(seqinfo)) NULL else seqinfo[],
    meta = meta
  )
  class(out) <- c("BwgTrack", "bwToolsTrack", "list")
  bw_assert_bwg(out)
  out
}

#' Test whether an object follows the standardized BwgTrack contract
#'
#' @param x An R object.
#' @return Logical scalar.
#' @details
#' Unlike `validate_bwg()`, this predicate returns `FALSE` for an invalid object
#' rather than raising an error.
#'
#' @examples
#' is_bwg_track(list())
#' @seealso [validate_bwg()], [bwg_track()]
#' @export
is_bwg_track <- function(x) {
  length(bw_contract_issues(x)) == 0L
}

#' Return sample metadata
#'
#' @param x A validated `BwgTrack` object.
#' @param sample_ids Optional sample IDs. Defaults to all samples.
#' @return A copy of the standardized sample metadata data.table.
#' @details
#' The returned table is a copy and can be modified without changing `x`.
#'
#' @examples
#' x <- bwg_track(samples = data.frame(sample_id = "sampleA"), data = data.frame(
#'   sample_id = "sampleA", chrom = "chr1", start = 1L, end = 2L,
#'   value = 1, strand = "*"
#' ))
#' samples_bwg(x)
#' @seealso [seqinfo_bwg()], [metadata_bwg()]
#' @export
samples_bwg <- function(x, sample_ids = NULL) {
  bw_assert_bwg(x)
  bw_select_samples(x, sample_ids = sample_ids)
}

#' Return BwgTrack metadata
#'
#' @param x A validated `BwgTrack` object.
#' @return The standardized metadata list.
#' @details
#' Standard metadata include `coordinate`, `schema_version`, `backend`, and
#' `mode`. Treat returned metadata as descriptive object state rather than a
#' mechanism for mutating the track contract.
#'
#' @examples
#' x <- bwg_track(samples = data.frame(sample_id = "sampleA"), data = data.frame(
#'   sample_id = "sampleA", chrom = "chr1", start = 1L, end = 2L,
#'   value = 1, strand = "*"
#' ))
#' metadata_bwg(x)
#' @seealso [samples_bwg()], [seqinfo_bwg()]
#' @export
metadata_bwg <- function(x) {
  bw_assert_bwg(x)
  x$meta
}

#' Print a bwTools signal track
#'
#' @param x A `bwToolsTrack` object.
#' @param ... Additional arguments reserved for future use.
#' @return `x`, invisibly.
#' @export
print.bwToolsTrack <- function(x, ...) {
  bw_assert_bwg(x)
  cat("<bwTools BwgTrack>\n")
  cat("  samples:", nrow(x$samples), "\n")
  cat("  mode:", x$meta$mode %||% if (is.null(x$data)) "lazy" else "memory", "\n")
  cat("  coordinate:", x$meta$coordinate %||% "1-based closed", "\n")
  cat("  schema:", x$meta$schema_version %||% NA_character_, "\n")
  if (!is.null(x$data)) {
    cat("  intervals:", nrow(x$data), "\n")
  }
  invisible(x)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
