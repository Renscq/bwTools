# Author: Rensc
# Date: 2026-08-28
# Version: dev003
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

bw_select_retrieve_samples <- function(object, samples = NULL, strand = "ignore") {
  sample_tbl <- data.table::copy(data.table::as.data.table(object$samples))

  if (!is.null(samples)) {
    selected <- unique(as.character(samples))
    if (length(selected) < 1L || anyNA(selected) || any(!nzchar(selected))) {
      bw_stop("`samples` must contain one or more non-empty sample IDs.")
    }
    missing_samples <- setdiff(selected, sample_tbl[["sample_id"]])
    if (length(missing_samples) > 0L) {
      bw_stop(paste0(
        "Unknown sample IDs: ",
        paste(missing_samples, collapse = ", "),
        "."
      ))
    }
    idx <- match(selected, sample_tbl[["sample_id"]])
    sample_tbl <- sample_tbl[idx]
  }

  query_strand <- as.character(strand)[1L]
  if (is.na(query_strand) || !query_strand %in% c("+", "-", "*", "both", "ignore")) {
    bw_stop("`strand` must be '+', '-', '*', 'both', or 'ignore'.")
  }
  if (query_strand %in% c("+", "-")) {
    idx <- which(sample_tbl[["strand"]] %in% c(query_strand, "*"))
    sample_tbl <- sample_tbl[idx]
  } else if (identical(query_strand, "*")) {
    idx <- which(sample_tbl[["strand"]] == "*")
    sample_tbl <- sample_tbl[idx]
  }

  sample_tbl[]
}

bw_retrieve_memory_regions <- function(object, sample_ids, regions) {
  dt <- data.table::copy(data.table::as.data.table(object$data))
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
#' @description
#' Retrieves one or more genomic regions from a `BwgTrack`-compatible object.
#' Query coordinates are 1-based closed. Overlapping or directly adjacent
#' regions are merged before retrieval, and signal intervals crossing query
#' boundaries are clipped without rebasing genomic coordinates.
#'
#' By default a signal data.table is returned for backward compatibility. Use
#' `result = "track"` when the retrieved data should remain in the standard
#' `BwgTrack` contract, for example before passing it to `write_bwg()`.
#'
#' @param object A `BwgTrack`-compatible object.
#' @param chrom Optional chromosome name for a single-region query.
#' @param start Optional start coordinate for a single-region query.
#' @param end Optional end coordinate for a single-region query.
#' @param samples Optional sample IDs. By default all samples are considered.
#' @param strand Optional strand selector: `+`, `-`, `*`, `both`, or `ignore`.
#' @param regions Optional data.frame or data.table containing `chrom`, `start`,
#'   and `end` columns. Use either `regions` or `chrom`/`start`/`end`.
#' @param result Result type: `data` returns the in-memory signal data.table;
#'   `track` returns a memory-mode `BwgTrack` containing the selected samples,
#'   signal, and chromosome metadata.
#' @return A signal data.table or a memory-mode `BwgTrack`, depending on
#'   `result`.
#' @export
retrieve_bwg <- function(
  object,
  chrom = NULL,
  start = NULL,
  end = NULL,
  samples = NULL,
  strand = "ignore",
  regions = NULL,
  result = c("data", "track")
) {
  if (!is_bwg_track(object)) {
    bw_stop("`object` must be a BwgTrack-compatible object.")
  }
  result <- match.arg(result)
  normalized_regions <- bw_normalize_retrieve_regions(
    chrom = chrom,
    start = start,
    end = end,
    regions = regions
  )
  sample_tbl <- bw_select_retrieve_samples(
    object,
    samples = samples,
    strand = strand
  )

  data <- if (nrow(sample_tbl) < 1L) {
    bw_empty_signal(include_sample = TRUE)
  } else if (is.null(object$data)) {
    bw_retrieve_lazy_regions(object, sample_tbl, normalized_regions)
  } else {
    bw_retrieve_memory_regions(
      object,
      sample_tbl[["sample_id"]],
      normalized_regions
    )
  }

  if (identical(result, "data")) {
    return(data[])
  }

  seqinfo <- bw_retrieve_seqinfo(
    object,
    sample_tbl[["sample_id"]],
    normalized_regions
  )
  meta <- object$meta
  meta$source_mode <- object$meta$mode %||%
    if (is.null(object$data)) "lazy" else "memory"
  meta$mode <- "memory"
  meta$retrieved <- TRUE
  meta$retrieval_regions <- data.table::copy(normalized_regions)

  bw_track(
    samples = sample_tbl,
    data = data,
    seqinfo = seqinfo,
    meta = meta
  )
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
