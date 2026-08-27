# Author: Rensc
# Date: 2026-08-27
# Version: dev002
# Function: Retrieve genomic intervals from BwgTrack-compatible objects
# Input: BwgTrack and genomic region
# Output: In-memory signal data.table

#' Retrieve signal intervals
#'
#' @param object A `BwgTrack`-compatible object.
#' @param chrom Chromosome name.
#' @param start Start coordinate in 1-based closed coordinates.
#' @param end End coordinate in 1-based closed coordinates.
#' @param samples Optional sample IDs.
#' @param strand Optional strand selector: `+`, `-`, `*`, `both`, or `ignore`.
#' @return A signal data.table in 1-based closed coordinates.
#' @export
retrieve_bwg <- function(object, chrom, start, end, samples = NULL, strand = "ignore") {
  if (!is_bwg_track(object)) {
    bw_stop("`object` must be a BwgTrack-compatible object.")
  }
  query_chrom <- as.character(chrom)[1L]
  query_start <- suppressWarnings(as.integer(start)[1L])
  query_end <- suppressWarnings(as.integer(end)[1L])
  query_strand <- as.character(strand)[1L]

  if (is.na(query_chrom) || !nzchar(query_chrom)) bw_stop("`chrom` must be non-empty.")
  if (is.na(query_start) || query_start < 1L) bw_stop("`start` must be >= 1.")
  if (is.na(query_end) || query_end < query_start) bw_stop("`end` must be >= `start`.")

  sample_tbl <- data.table::copy(data.table::as.data.table(object$samples))
  if (!is.null(samples)) {
    sample_tbl <- sample_tbl[sample_id %in% as.character(samples)]
  }
  if (query_strand %in% c("+", "-")) {
    sample_tbl <- sample_tbl[strand %in% c(query_strand, "*")]
  } else if (identical(query_strand, "*")) {
    sample_tbl <- sample_tbl[strand == "*"]
  } else if (!query_strand %in% c("both", "ignore")) {
    bw_stop("`strand` must be '+', '-', '*', 'both', or 'ignore'.")
  }
  if (nrow(sample_tbl) == 0L) {
    return(bw_empty_signal(include_sample = TRUE))
  }

  if (!is.null(object$data)) {
    dt <- data.table::copy(data.table::as.data.table(object$data))
    dt <- dt[
      sample_id %in% sample_tbl$sample_id &
        chrom == query_chrom &
        end >= query_start &
        start <= query_end
    ]
    if (nrow(dt) == 0L) {
      return(bw_empty_signal(include_sample = TRUE))
    }
    data.table::set(dt, j = "start", value = pmax(as.integer(dt[["start"]]), query_start))
    data.table::set(dt, j = "end", value = pmin(as.integer(dt[["end"]]), query_end))
    data.table::setorderv(dt, c("sample_id", "chrom", "start", "end"))
    return(dt[])
  }

  out <- vector("list", nrow(sample_tbl))
  out_n <- 0L
  for (i in seq_len(nrow(sample_tbl))) {
    if (sample_tbl$format[i] != "bigwig") {
      bw_stop("Lazy retrieval is currently implemented for BigWig sources only.")
    }
    x <- bw_bigwig_query(sample_tbl$file[i], query_chrom, query_start, query_end)
    if (nrow(x) == 0L) next
    data.table::set(x, j = "sample_id", value = sample_tbl$sample_id[i])
    data.table::set(x, j = "strand", value = sample_tbl$strand[i])
    data.table::setcolorder(x, c("sample_id", "chrom", "start", "end", "value", "strand"))
    out_n <- out_n + 1L
    out[[out_n]] <- x
  }
  if (out_n == 0L) {
    return(bw_empty_signal(include_sample = TRUE))
  }
  ans <- data.table::rbindlist(out[seq_len(out_n)], use.names = TRUE)
  data.table::setorderv(ans, c("sample_id", "chrom", "start", "end"))
  ans[]
}

#' Return chromosome metadata
#'
#' @param x A local BigWig path or a `BwgTrack`-compatible object.
#' @return A data.table containing chromosome metadata.
#' @export
seqinfo_bwg <- function(x) {
  if (is_bwg_track(x)) {
    if (!is.null(x$seqinfo)) {
      return(data.table::copy(data.table::as.data.table(x$seqinfo)))
    }
    sample_tbl <- data.table::as.data.table(x$samples)
    if (all(sample_tbl$format == "bigwig")) {
      out <- lapply(seq_len(nrow(sample_tbl)), function(i) {
        z <- bw_bigwig_seqinfo(sample_tbl$file[i])
        data.table::set(z, j = "sample_id", value = sample_tbl$sample_id[i])
        data.table::setcolorder(z, c("sample_id", "chrom", "length"))
        z
      })
      return(data.table::rbindlist(out))
    }
    return(data.table::data.table())
  }
  bw_bigwig_seqinfo(x)
}

#' Summarize a BwgTrack-compatible object
#'
#' @param object A `BwgTrack`-compatible object.
#' @return A named list with object-level counts.
#' @export
summary_bwg <- function(object) {
  if (!is_bwg_track(object)) {
    bw_stop("`object` must be a BwgTrack-compatible object.")
  }
  list(
    samples = nrow(object$samples),
    mode = object$meta$mode %||% if (is.null(object$data)) "lazy" else "memory",
    intervals = if (is.null(object$data)) NA_integer_ else nrow(object$data),
    chromosomes = if (is.null(object$seqinfo)) NA_integer_ else data.table::uniqueN(object$seqinfo$chrom),
    coordinate = object$meta$coordinate %||% "1-based closed",
    backend = object$meta$backend %||% "bwTools-native-R"
  )
}
