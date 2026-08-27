# Author: Rensc
# Date: 2026-08-27
# Version: dev003
# Function: Define the GeneTrackR-compatible genomic signal object
# Input: Sample metadata, signal data, sequence metadata, and package metadata
# Output: BwgTrack-compatible S3 object

#' Construct a BwgTrack-compatible signal object
#'
#' @param samples Sample metadata table. Must contain `sample_id`.
#' @param data Optional in-memory signal table.
#' @param seqinfo Optional chromosome metadata table.
#' @param meta Optional metadata list.
#' @return An object inheriting from `BwgTrack` and `bwToolsTrack`.
#' @export
bw_track <- function(samples, data = NULL, seqinfo = NULL, meta = list()) {
  samples <- data.table::copy(data.table::as.data.table(samples))
  if (!"sample_id" %in% names(samples)) {
    bw_stop("`samples` must contain a `sample_id` column.")
  }
  data.table::set(samples, j = "sample_id", value = as.character(samples[["sample_id"]]))
  if (anyNA(samples$sample_id) || any(!nzchar(samples$sample_id))) {
    bw_stop("`sample_id` values must be non-missing and non-empty.")
  }
  if (anyDuplicated(samples$sample_id)) {
    bw_stop("`sample_id` values must be unique.")
  }

  if (!"strand" %in% names(samples)) {
    data.table::set(samples, j = "strand", value = "*")
  }
  data.table::set(samples, j = "strand", value = as.character(samples[["strand"]]))
  if (!"has_strand" %in% names(samples)) {
    data.table::set(samples, j = "has_strand", value = samples[["strand"]] %in% c("+", "-"))
  }

  if (!is.null(data)) {
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
    data.table::set(data, j = "start", value = as.integer(data[["start"]]))
    data.table::set(data, j = "end", value = as.integer(data[["end"]]))
    data.table::set(data, j = "value", value = as.numeric(data[["value"]]))
    data.table::set(data, j = "strand", value = as.character(data[["strand"]]))
    invalid <- is.na(data$start) | is.na(data$end) | data$start < 1L | data$end < data$start
    if (any(invalid)) {
      bw_stop("Signal data contain invalid 1-based closed intervals.")
    }
    unknown <- setdiff(unique(data$sample_id), samples$sample_id)
    if (length(unknown) > 0L) {
      bw_stop(paste0("Signal data contain unknown sample IDs: ", paste(unknown, collapse = ", "), "."))
    }
    data.table::setorderv(data, c("sample_id", "chrom", "start", "end"))
    data.table::setindexv(data, c("sample_id", "chrom", "start", "end"))
  }

  if (!is.null(seqinfo)) {
    seqinfo <- data.table::as.data.table(seqinfo)
  }

  meta$coordinate <- "1-based closed"
  meta$schema_version <- "1"
  meta$backend <- "bwTools-native-R"

  out <- list(
    samples = samples[],
    data = if (is.null(data)) NULL else data[],
    seqinfo = if (is.null(seqinfo)) NULL else seqinfo[],
    meta = meta
  )
  class(out) <- c("BwgTrack", "bwToolsTrack", "list")
  out
}

#' Test whether an object follows the BwgTrack contract
#'
#' @param x An R object.
#' @return Logical scalar.
#' @export
is_bwg_track <- function(x) {
  inherits(x, "BwgTrack") && is.list(x) && !is.null(x$samples) && !is.null(x$meta)
}

#' Print a bwTools signal track
#'
#' @param x A `bwToolsTrack` object.
#' @param ... Additional arguments reserved for future use.
#' @return `x`, invisibly.
#' @export
print.bwToolsTrack <- function(x, ...) {
  cat("<bwTools BwgTrack>\n")
  cat("  samples:", nrow(x$samples), "\n")
  cat("  mode:", x$meta$mode %||% if (is.null(x$data)) "lazy" else "memory", "\n")
  cat("  coordinate:", x$meta$coordinate %||% "1-based closed", "\n")
  if (!is.null(x$data)) {
    cat("  intervals:", nrow(x$data), "\n")
  }
  invisible(x)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
