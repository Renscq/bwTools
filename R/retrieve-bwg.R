# Author: Rensc
# Date: 2026-08-29
# Version: dev006
# Function: Retrieve one or more genomic regions from BwgTrack-compatible objects
# Input: BwgTrack and genomic region(s)
# Output: In-memory signal data.table or BwgTrack subset

bw_normalize_retrieve_regions <- function(chrom = NULL, start = NULL, end = NULL, regions = NULL) {
  if (is.null(regions)) {
    if (is.null(chrom) || is.null(start) || is.null(end)) {
      bw_stop("Provide either `regions` or all of `chrom`, `start`, and `end`.")
    }
    if (length(chrom) != 1L || length(start) != 1L || length(end) != 1L) {
      bw_stop("`chrom`, `start`, and `end` must be scalar when `regions` is not supplied.")
    }
    regions <- data.table::data.table(
      chrom = chrom,
      start = start,
      end = end
    )
  } else {
    if (!is.null(chrom) || !is.null(start) || !is.null(end)) {
      bw_stop("Use either `regions` or `chrom`/`start`/`end`, not both.")
    }
  }

  x <- data.table::copy(data.table::as.data.table(regions))
  required <- c("chrom", "start", "end")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    bw_stop(paste0(
      "`regions` is missing required columns: ",
      paste(missing, collapse = ", "),
      "."
    ))
  }
  if (nrow(x) < 1L) {
    bw_stop("`regions` must contain at least one genomic interval.")
  }

  x <- x[, required, with = FALSE]
  data.table::set(x, j = "chrom", value = as.character(x[["chrom"]]))
  start_value <- suppressWarnings(as.numeric(x[["start"]]))
  end_value <- suppressWarnings(as.numeric(x[["end"]]))

  invalid <- is.na(x[["chrom"]]) | !nzchar(x[["chrom"]]) |
    !is.finite(start_value) | !is.finite(end_value) |
    start_value < 1 | end_value < start_value |
    start_value != floor(start_value) | end_value != floor(end_value) |
    start_value > .Machine$integer.max | end_value > .Machine$integer.max
  if (any(invalid)) {
    bw_stop(
      "Genomic regions must be valid 1-based closed intervals within the R integer range."
    )
  }

  data.table::set(x, j = "start", value = as.integer(start_value))
  data.table::set(x, j = "end", value = as.integer(end_value))
  chrom_order <- unique(x[["chrom"]])
  data.table::set(x, j = "chrom_order", value = match(x[["chrom"]], chrom_order))
  data.table::setorderv(x, c("chrom_order", "start", "end"))

  out <- vector("list", nrow(x))
  out_n <- 0L
  for (chrom_name in chrom_order) {
    chrom_idx <- which(x[["chrom"]] == chrom_name)
    y <- x[chrom_idx]
    if (nrow(y) < 1L) next

    current_start <- y[["start"]][1L]
    current_end <- y[["end"]][1L]
    if (nrow(y) > 1L) {
      for (i in 2L:nrow(y)) {
        next_start <- y[["start"]][i]
        next_end <- y[["end"]][i]
        touches <- as.numeric(next_start) <= as.numeric(current_end) + 1
        if (touches) {
          current_end <- max(current_end, next_end)
        } else {
          out_n <- out_n + 1L
          out[[out_n]] <- data.table::data.table(
            chrom = chrom_name,
            start = current_start,
            end = current_end
          )
          current_start <- next_start
          current_end <- next_end
        }
      }
    }
    out_n <- out_n + 1L
    out[[out_n]] <- data.table::data.table(
      chrom = chrom_name,
      start = current_start,
      end = current_end
    )
  }

  ans <- data.table::rbindlist(out[seq_len(out_n)], use.names = TRUE)
  data.table::set(ans, j = "start", value = as.integer(ans[["start"]]))
  data.table::set(ans, j = "end", value = as.integer(ans[["end"]]))
  ans[]
}

bw_select_retrieve_samples <- function(x, sample_ids = NULL, strand = NULL) {
  sample_tbl <- bw_select_samples(x, sample_ids = sample_ids)

  if (is.null(strand)) {
    return(sample_tbl[])
  }
  query_strand <- as.character(strand)
  if (length(query_strand) != 1L || is.na(query_strand) ||
      !query_strand %in% c("+", "-", "*")) {
    bw_stop("`strand` must be NULL, '+', '-', or '*'.")
  }
  idx <- which(sample_tbl[["strand"]] == query_strand)
  sample_tbl[idx][]
}

bw_retrieve_memory_regions <- function(object, sample_ids, regions) {
  dt <- data.table::as.data.table(object$data)
  sample_idx <- which(dt[["sample_id"]] %in% sample_ids)
  dt <- dt[sample_idx]
  if (nrow(dt) < 1L) {
    return(bw_empty_signal(include_sample = TRUE))
  }

  out <- vector("list", nrow(regions))
  out_n <- 0L
  for (i in seq_len(nrow(regions))) {
    chrom_name <- regions[["chrom"]][i]
    region_start <- regions[["start"]][i]
    region_end <- regions[["end"]][i]
    idx <- which(
      dt[["chrom"]] == chrom_name &
        dt[["end"]] >= region_start &
        dt[["start"]] <= region_end
    )
    if (length(idx) < 1L) next

    x <- data.table::copy(dt[idx])
    data.table::set(
      x,
      j = "start",
      value = pmax(as.integer(x[["start"]]), region_start)
    )
    data.table::set(
      x,
      j = "end",
      value = pmin(as.integer(x[["end"]]), region_end)
    )
    out_n <- out_n + 1L
    out[[out_n]] <- x
  }

  if (out_n < 1L) {
    return(bw_empty_signal(include_sample = TRUE))
  }
  if (out_n == 1L) {
    return(out[[1L]][])
  }
  ans <- data.table::rbindlist(out[seq_len(out_n)], use.names = TRUE)
  data.table::setorderv(ans, c("sample_id", "chrom", "start", "end"))
  ans[]
}

bw_retrieve_lazy_regions <- function(object, sample_tbl, regions) {
  if (any(sample_tbl[["format"]] != "bigwig")) {
    bw_stop("Lazy retrieval is currently implemented for BigWig sources only.")
  }

  seqinfo <- if (is.null(object$seqinfo)) {
    NULL
  } else {
    data.table::as.data.table(object$seqinfo)
  }

  out <- vector("list", nrow(sample_tbl) * nrow(regions))
  out_n <- 0L
  for (sample_i in seq_len(nrow(sample_tbl))) {
    sample_id <- sample_tbl[["sample_id"]][sample_i]
    available_chrom <- NULL
    if (!is.null(seqinfo) && "chrom" %in% names(seqinfo)) {
      if ("sample_id" %in% names(seqinfo)) {
        seq_idx <- which(seqinfo[["sample_id"]] == sample_id)
        available_chrom <- unique(seqinfo[["chrom"]][seq_idx])
      } else {
        available_chrom <- unique(seqinfo[["chrom"]])
      }
    }

    for (region_i in seq_len(nrow(regions))) {
      chrom_name <- regions[["chrom"]][region_i]
      if (!is.null(available_chrom) && !chrom_name %in% available_chrom) {
        next
      }
      x <- bw_bigwig_query(
        sample_tbl[["file"]][sample_i],
        chrom_name,
        regions[["start"]][region_i],
        regions[["end"]][region_i]
      )
      if (nrow(x) < 1L) next
      data.table::set(x, j = "sample_id", value = sample_id)
      data.table::set(x, j = "strand", value = sample_tbl[["strand"]][sample_i])
      data.table::setcolorder(
        x,
        c("sample_id", "chrom", "start", "end", "value", "strand")
      )
      out_n <- out_n + 1L
      out[[out_n]] <- x
    }
  }

  if (out_n < 1L) {
    return(bw_empty_signal(include_sample = TRUE))
  }
  if (out_n == 1L) {
    return(out[[1L]][])
  }
  ans <- data.table::rbindlist(out[seq_len(out_n)], use.names = TRUE)
  data.table::setorderv(ans, c("sample_id", "chrom", "start", "end"))
  ans[]
}

bw_retrieve_seqinfo <- function(object, sample_ids, regions) {
  if (is.null(object$seqinfo)) {
    return(NULL)
  }
  x <- data.table::copy(data.table::as.data.table(object$seqinfo))
  if (!"chrom" %in% names(x)) {
    return(NULL)
  }
  chrom_idx <- which(x[["chrom"]] %in% regions[["chrom"]])
  x <- x[chrom_idx]
  if ("sample_id" %in% names(x)) {
    sample_idx <- which(x[["sample_id"]] %in% sample_ids)
    x <- x[sample_idx]
  }
  x[]
}

#' Retrieve signal intervals
#'
#' Retrieves one or more genomic regions from a standardized `BwgTrack`.
#' Query coordinates are 1-based closed. Overlapping or directly adjacent
#' regions are merged before retrieval, and signal intervals crossing query
#' boundaries are clipped without rebasing genomic coordinates.
#'
#' @param x A validated `BwgTrack` object.
#' @param chrom Optional chromosome name for a single-region query.
#' @param start Optional start coordinate for a single-region query.
#' @param end Optional end coordinate for a single-region query.
#' @param regions Optional data.frame or data.table containing `chrom`, `start`,
#'   and `end` columns. Use either `regions` or `chrom`/`start`/`end`.
#' @param sample_ids Optional sample IDs. Defaults to all samples.
#' @param strand Optional exact strand selector: `+`, `-`, or `*`. `NULL`
#'   performs no strand filtering.
#' @param result Result type. `data` returns the standardized signal data.table;
#'   `track` returns a memory-mode `BwgTrack`.
#' @return A signal data.table or memory-mode `BwgTrack`.
#' @export
retrieve_bwg <- function(
  x,
  chrom = NULL,
  start = NULL,
  end = NULL,
  regions = NULL,
  sample_ids = NULL,
  strand = NULL,
  result = c("data", "track")
) {
  bw_assert_bwg(x)
  result <- bw_match_arg(result, c("data", "track"), "result")
  normalized_regions <- bw_normalize_retrieve_regions(
    chrom = chrom,
    start = start,
    end = end,
    regions = regions
  )
  sample_tbl <- bw_select_retrieve_samples(
    x,
    sample_ids = sample_ids,
    strand = strand
  )

  data <- if (nrow(sample_tbl) < 1L) {
    bw_empty_signal(include_sample = TRUE)
  } else if (is.null(x$data)) {
    bw_retrieve_lazy_regions(x, sample_tbl, normalized_regions)
  } else {
    bw_retrieve_memory_regions(
      x,
      sample_tbl[["sample_id"]],
      normalized_regions
    )
  }

  if (identical(result, "data")) {
    return(data[])
  }

  seqinfo <- bw_retrieve_seqinfo(
    x,
    sample_tbl[["sample_id"]],
    normalized_regions
  )
  meta <- x$meta
  meta$source_mode <- x$meta$mode %||% if (is.null(x$data)) "lazy" else "memory"
  meta$retrieved <- TRUE
  meta$retrieval_regions <- data.table::copy(normalized_regions)

  bwg_track(
    samples = sample_tbl,
    data = data,
    seqinfo = seqinfo,
    meta = meta
  )
}

#' Return sample-specific chromosome metadata
#'
#' @param x A validated `BwgTrack` object.
#' @param sample_ids Optional sample IDs. Defaults to all samples.
#' @return A copy of the standardized `sample_id`, `chrom`, and `length`
#'   metadata table.
#' @export
seqinfo_bwg <- function(x, sample_ids = NULL) {
  bw_assert_bwg(x)
  sample_tbl <- bw_select_samples(x, sample_ids = sample_ids)
  selected <- sample_tbl[["sample_id"]]

  if (!is.null(x$seqinfo)) {
    seqinfo <- data.table::copy(data.table::as.data.table(x$seqinfo))
    idx <- which(seqinfo[["sample_id"]] %in% selected)
    return(seqinfo[idx][])
  }

  if (nrow(sample_tbl) > 0L && all(sample_tbl[["format"]] == "bigwig")) {
    out <- vector("list", nrow(sample_tbl))
    for (i in seq_len(nrow(sample_tbl))) {
      z <- bw_bigwig_seqinfo(sample_tbl[["file"]][i])
      data.table::set(z, j = "sample_id", value = sample_tbl[["sample_id"]][i])
      data.table::setcolorder(z, c("sample_id", "chrom", "length"))
      out[[i]] <- z
    }
    return(data.table::rbindlist(out, use.names = TRUE, fill = TRUE)[])
  }

  data.table::data.table(
    sample_id = character(),
    chrom = character(),
    length = integer()
  )
}

#' Summarize a BwgTrack object
#'
#' @param x A validated `BwgTrack` object.
#' @return A one-row data.table describing the object contract and dimensions.
#' @export
summary_bwg <- function(x) {
  bw_assert_bwg(x)
  data.table::data.table(
    samples = nrow(x$samples),
    mode = x$meta$mode,
    intervals = if (is.null(x$data)) NA_integer_ else nrow(x$data),
    chromosomes = if (is.null(x$seqinfo)) {
      NA_integer_
    } else {
      data.table::uniqueN(x$seqinfo[["chrom"]])
    },
    coordinate = x$meta$coordinate,
    schema_version = x$meta$schema_version,
    backend = x$meta$backend
  )
}
