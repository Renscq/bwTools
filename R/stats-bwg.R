# Author: Rensc
# Date: 2026-08-29
# Version: dev006
# Function: Calculate bounded exact or zoom-accelerated interval statistics
# Input: Validated BwgTrack object and genomic bins
# Output: Per-bin signal statistics in 1-based closed coordinates

bw_empty_stats <- function() {
  data.table::data.table(
    sample_id = character(),
    chrom = character(),
    start = integer(),
    end = integer(),
    stat = character(),
    value = numeric(),
    covered_bases = numeric(),
    source = character(),
    resolution = numeric()
  )
}

bw_stat_bins <- function(start, end, n_bins) {
  start <- suppressWarnings(as.integer(start)[1L])
  end <- suppressWarnings(as.integer(end)[1L])
  n_bins <- suppressWarnings(as.integer(n_bins)[1L])
  if (is.na(start) || start < 1L) bw_stop("`start` must be >= 1.")
  if (is.na(end) || end < start) bw_stop("`end` must be >= `start`.")
  width <- as.numeric(end) - as.numeric(start) + 1
  if (is.na(n_bins) || n_bins < 1L) bw_stop("`n_bins` must be >= 1.")
  if (n_bins > width) {
    bw_stop("`n_bins` cannot exceed the number of bases in the requested interval.")
  }
  query_start0 <- as.numeric(start) - 1
  boundaries <- query_start0 + floor(seq.int(0L, n_bins) * width / n_bins)
  data.table::data.table(
    bin = seq_len(n_bins),
    start0 = boundaries[-length(boundaries)],
    end0 = boundaries[-1L],
    start = as.integer(boundaries[-length(boundaries)] + 1),
    end = as.integer(boundaries[-1L]),
    width = diff(boundaries)
  )
}

bw_finish_stats <- function(
  bins,
  covered,
  sum_data,
  sum_squared,
  min_value,
  max_value,
  stat,
  source,
  resolution
) {
  has_data <- covered > 0

  value <- switch(
    stat,
    mean = {
      out <- rep(NA_real_, length(covered))
      out[has_data] <- sum_data[has_data] / covered[has_data]
      out
    },
    stdev = {
      out <- rep(NA_real_, length(covered))
      one <- covered == 1
      many <- covered > 1
      out[one] <- 0
      if (any(many)) {
        variance <- (
          sum_squared[many] - (sum_data[many]^2 / covered[many])
        ) / (covered[many] - 1)
        variance[variance < 0 & abs(variance) < 1e-10] <- 0
        variance[variance < 0] <- NA_real_
        out[many] <- sqrt(variance)
      }
      out
    },
    max = ifelse(has_data, max_value, NA_real_),
    min = ifelse(has_data, min_value, NA_real_),
    coverage = ifelse(has_data, covered / bins$width, NA_real_),
    sum = ifelse(has_data, sum_data, NA_real_),
    bw_stop("Unsupported statistic.")
  )

  data.table::data.table(
    chrom = NA_character_,
    start = bins$start,
    end = bins$end,
    stat = stat,
    value = as.numeric(value),
    covered_bases = as.numeric(covered),
    source = source,
    resolution = if (is.na(resolution)) NA_real_ else as.numeric(resolution)
  )
}

bw_new_stat_state <- function(n) {
  list(
    covered = numeric(n),
    sum_data = numeric(n),
    sum_squared = numeric(n),
    min_value = rep(Inf, n),
    max_value = rep(-Inf, n)
  )
}

bw_accumulate_signal_vectors <- function(state, start0, end0, value, bins, stat) {
  if (length(start0) < 1L) {
    return(state)
  }

  start0 <- as.numeric(start0)
  end0 <- as.numeric(end0)
  value <- as.numeric(value)
  keep <- is.finite(start0) & is.finite(end0) & is.finite(value) & end0 > start0
  if (!any(keep)) {
    return(state)
  }
  start0 <- start0[keep]
  end0 <- end0[keep]
  value <- value[keep]

  n <- nrow(bins)
  first_bin <- findInterval(start0, bins$end0) + 1L
  last_bin <- findInterval(end0 - 1, bins$end0) + 1L
  keep <- first_bin >= 1L & first_bin <= n & last_bin >= 1L & last_bin <= n
  if (!any(keep)) {
    return(state)
  }
  start0 <- start0[keep]
  end0 <- end0[keep]
  value <- value[keep]
  first_bin <- first_bin[keep]
  last_bin <- last_bin[keep]

  need_sum <- stat %in% c("mean", "stdev", "sum")
  need_squared <- identical(stat, "stdev")
  need_min <- identical(stat, "min")
  need_max <- identical(stat, "max")

  same <- first_bin == last_bin
  same_idx <- which(same)
  if (length(same_idx) > 0L) {
    same_bins <- first_bin[same_idx]
    touched <- unique(same_bins)
    for (bin_id in touched) {
      idx <- same_idx[same_bins == bin_id]
      width <- end0[idx] - start0[idx]
      vals <- value[idx]
      state$covered[bin_id] <- state$covered[bin_id] + sum(width)
      if (need_sum) {
        state$sum_data[bin_id] <- state$sum_data[bin_id] + sum(width * vals)
      }
      if (need_squared) {
        state$sum_squared[bin_id] <- state$sum_squared[bin_id] + sum(width * vals^2)
      }
      if (need_min) {
        state$min_value[bin_id] <- min(state$min_value[bin_id], min(vals))
      }
      if (need_max) {
        state$max_value[bin_id] <- max(state$max_value[bin_id], max(vals))
      }
    }
  }

  cross_idx <- which(!same)
  if (length(cross_idx) > 0L) {
    for (i in cross_idx) {
      for (bin_id in seq.int(first_bin[i], last_bin[i])) {
        segment_start <- max(start0[i], bins$start0[bin_id])
        segment_end <- min(end0[i], bins$end0[bin_id])
        width <- segment_end - segment_start
        if (width <= 0) next
        state$covered[bin_id] <- state$covered[bin_id] + width
        if (need_sum) {
          state$sum_data[bin_id] <- state$sum_data[bin_id] + width * value[i]
        }
        if (need_squared) {
          state$sum_squared[bin_id] <- state$sum_squared[bin_id] + width * value[i]^2
        }
        if (need_min) {
          state$min_value[bin_id] <- min(state$min_value[bin_id], value[i])
        }
        if (need_max) {
          state$max_value[bin_id] <- max(state$max_value[bin_id], value[i])
        }
      }
    }
  }
  state
}

bw_stats_from_signal <- function(signal, chrom, start, end, n_bins, stat) {
  bins <- bw_stat_bins(start, end, n_bins)
  state <- bw_new_stat_state(nrow(bins))

  x <- data.table::as.data.table(signal)
  if (nrow(x) > 0L) {
    start0 <- pmax(as.numeric(x$start) - 1, as.numeric(start) - 1)
    end0 <- pmin(as.numeric(x$end), as.numeric(end))
    state <- bw_accumulate_signal_vectors(
      state = state,
      start0 = start0,
      end0 = end0,
      value = x$value,
      bins = bins,
      stat = stat
    )
  }

  out <- bw_finish_stats(
    bins,
    state$covered,
    state$sum_data,
    state$sum_squared,
    state$min_value,
    state$max_value,
    stat = stat,
    source = "full",
    resolution = NA_integer_
  )
  data.table::set(out, j = "chrom", value = rep(chrom, nrow(out)))
  out[]
}

bw_stats_bigwig_full_stream <- function(
  file,
  chrom,
  start,
  end,
  n_bins,
  stat,
  metadata = NULL
) {
  file <- bw_validate_local_file(file)
  if (is.null(metadata)) {
    metadata <- bw_metadata(file)
  }
  chrom_idx <- match(chrom, metadata$chromosomes$chrom)
  if (is.na(chrom_idx)) {
    bw_stop(paste0("Chromosome `", chrom, "` was not found in the BigWig file."))
  }

  chrom_length <- as.numeric(metadata$chromosomes$length[chrom_idx])
  query_start <- as.numeric(start) - 1
  query_end <- min(as.numeric(end), chrom_length)
  if (query_start >= query_end) {
    return(data.table::data.table())
  }

  start <- as.integer(query_start + 1)
  end <- as.integer(query_end)
  bins <- bw_stat_bins(start, end, n_bins)
  state <- bw_new_stat_state(nrow(bins))

  con <- base::file(file, open = "rb")
  on.exit(close(con), add = TRUE)
  index <- bw_read_index_header(con, metadata$header)
  blocks <- bw_collect_blocks(
    con = con,
    node_offset = index$root_offset,
    endian = metadata$header$endian,
    tid = metadata$chromosomes$tid[chrom_idx],
    query_start = query_start,
    query_end = query_end
  )

  if (nrow(blocks) > 0L) {
    for (i in seq_len(nrow(blocks))) {
      block <- bw_read_exact(con, blocks$offset[i], blocks$size[i])
      if (metadata$header$uncompress_buf_size > 0) {
        block <- bw_decompress_block(block)
      }
      decoded <- bw_decode_block_vectors(
        block,
        metadata$header$endian,
        metadata$chromosomes$tid[chrom_idx],
        query_start,
        query_end
      )
      state <- bw_accumulate_signal_vectors(
        state = state,
        start0 = decoded$start0,
        end0 = decoded$end0,
        value = decoded$value,
        bins = bins,
        stat = stat
      )
    }
  }

  out <- bw_finish_stats(
    bins,
    state$covered,
    state$sum_data,
    state$sum_squared,
    state$min_value,
    state$max_value,
    stat = stat,
    source = "full",
    resolution = NA_integer_
  )
  data.table::set(out, j = "chrom", value = rep(chrom, nrow(out)))
  out[]
}

bw_stats_from_zoom <- function(records, chrom, start, end, n_bins, stat, resolution) {
  bins <- bw_stat_bins(start, end, n_bins)
  n <- nrow(bins)
  covered <- numeric(n)
  sum_data <- numeric(n)
  sum_squared <- numeric(n)
  min_value <- rep(Inf, n)
  max_value <- rep(-Inf, n)

  x <- data.table::as.data.table(records)
  if (nrow(x) > 0L) {
    for (i in seq_len(nrow(x))) {
      record_start <- max(as.numeric(x$start0[i]), as.numeric(start) - 1)
      record_end <- min(as.numeric(x$end0[i]), as.numeric(end))
      record_width <- as.numeric(x$end0[i]) - as.numeric(x$start0[i])
      if (record_end <= record_start || record_width <= 0) next
      pos <- record_start
      while (pos < record_end) {
        j <- findInterval(pos, bins$end0) + 1L
        if (j < 1L) j <- 1L
        if (j > n) break
        segment_end <- min(record_end, bins$end0[j])
        overlap <- segment_end - pos
        if (overlap <= 0) break
        scalar <- overlap / record_width
        covered[j] <- covered[j] + as.numeric(x$valid_count[i]) * scalar
        sum_data[j] <- sum_data[j] + as.numeric(x$sum_data[i]) * scalar
        sum_squared[j] <- sum_squared[j] + as.numeric(x$sum_squared[i]) * scalar
        min_value[j] <- min(min_value[j], as.numeric(x$min_value[i]))
        max_value[j] <- max(max_value[j], as.numeric(x$max_value[i]))
        pos <- segment_end
      }
    }
  }

  out <- bw_finish_stats(
    bins, covered, sum_data, sum_squared, min_value, max_value,
    stat = stat, source = "zoom", resolution = resolution
  )
  data.table::set(out, j = "chrom", value = rep(chrom, nrow(out)))
  out[]
}

bw_stats_bigwig_file <- function(file, sample_id, chrom, start, end, n_bins, stat, use_zoom) {
  metadata <- bw_metadata(file)
  chrom_idx <- match(chrom, metadata$chromosomes$chrom)
  if (is.na(chrom_idx)) {
    bw_stop(paste0("Chromosome `", chrom, "` was not found in BigWig sample `", sample_id, "`."))
  }
  chrom_length <- as.numeric(metadata$chromosomes$length[chrom_idx])
  if (start > chrom_length) {
    return(data.table::data.table())
  }
  end <- min(as.integer(end), as.integer(chrom_length))
  bins <- bw_stat_bins(start, end, n_bins)

  zoom <- if (isTRUE(use_zoom)) bw_choose_zoom_level(metadata, start, end, n_bins) else NULL
  if (!is.null(zoom) && nrow(zoom) > 0L) {
    records <- bw_bigwig_zoom_query(
      file, chrom, start, end,
      level = zoom$level[1L], n_bins = n_bins
    )
    out <- bw_stats_from_zoom(
      records, chrom, start, end, n_bins, stat,
      resolution = zoom$level[1L]
    )
  } else {
    out <- bw_stats_bigwig_full_stream(
      file = file,
      chrom = chrom,
      start = start,
      end = end,
      n_bins = n_bins,
      stat = stat,
      metadata = metadata
    )
  }
  data.table::set(out, j = "sample_id", value = rep(sample_id, nrow(out)))
  data.table::setcolorder(
    out,
    c("sample_id", "chrom", "start", "end", "stat", "value", "covered_bases", "source", "resolution")
  )
  out[]
}

#' Calculate interval statistics
#'
#' Calculates coverage-aware statistics over equal-width genomic bins from a
#' standardized `BwgTrack`. Memory-mode tracks use full-resolution signal.
#' Lazy BigWig tracks may use stored zoom summaries when `use_zoom = TRUE`.
#' Exact lazy-BigWig calculations stream overlapping full-resolution blocks
#' directly into bin accumulators without materializing an interval table. When
#' all chromosome lengths for a selected sample are known in `seqinfo`, query ends
#' are clipped to the known chromosome length and unknown chromosomes are
#' rejected. Queries starting beyond a known chromosome return a stable empty
#' statistics table.
#'
#' @param x A validated `BwgTrack` object.
#' @param chrom Chromosome name.
#' @param start Start coordinate, 1-based closed.
#' @param end End coordinate, 1-based closed.
#' @param n_bins Number of equal-width bins.
#' @param stat Statistic: `mean`, `stdev`, `max`, `min`, `coverage`, or `sum`.
#' @param sample_ids Optional sample IDs. Defaults to all samples.
#' @param use_zoom Whether lazy BigWig inputs may use zoom summary levels.
#'   Zoom results may be approximate when query bins partially overlap stored
#'   zoom windows. Use `FALSE` for exact values.
#' @return A data.table with one row per sample and genomic bin.
#' @export
stats_bwg <- function(
  x,
  chrom,
  start,
  end,
  n_bins = 1L,
  stat = c("mean", "stdev", "max", "min", "coverage", "sum"),
  sample_ids = NULL,
  use_zoom = TRUE
) {
  bw_assert_bwg(x)
  stat <- bw_match_arg(
    stat,
    c("mean", "stdev", "max", "min", "coverage", "sum"),
    "stat"
  )
  use_zoom <- bw_validate_flag(use_zoom, "use_zoom")
  interval <- bw_validate_genomic_interval(chrom, start, end)
  chrom <- interval$chrom
  start <- interval$start
  end <- interval$end
  n_bins <- bw_validate_positive_integer(n_bins, "n_bins", minimum = 1L)

  sample_tbl <- bw_select_samples(x, sample_ids = sample_ids)
  if (nrow(sample_tbl) < 1L) {
    return(bw_empty_stats())
  }

  bounded <- bw_bound_query_interval(
    x,
    sample_ids = sample_tbl[["sample_id"]],
    chrom = chrom,
    start = start,
    end = end
  )
  if (isTRUE(bounded$empty)) {
    return(bw_empty_stats())
  }
  start <- as.integer(bounded$start)
  end <- as.integer(bounded$end)
  bw_stat_bins(start, end, n_bins)

  out <- vector("list", nrow(sample_tbl))
  out_n <- 0L
  if (!is.null(x$data)) {
    dt <- data.table::as.data.table(x$data)
    for (i in seq_len(nrow(sample_tbl))) {
      sid <- sample_tbl[["sample_id"]][i]
      keep <-
        dt[["sample_id"]] == sid &
        dt[["chrom"]] == chrom &
        dt[["end"]] >= start &
        dt[["start"]] <= end
      signal <- dt[which(keep)]
      result <- bw_stats_from_signal(signal, chrom, start, end, n_bins, stat)
      data.table::set(result, j = "sample_id", value = rep(sid, nrow(result)))
      data.table::setcolorder(
        result,
        c(
          "sample_id", "chrom", "start", "end", "stat", "value",
          "covered_bases", "source", "resolution"
        )
      )
      out_n <- out_n + 1L
      out[[out_n]] <- result
    }
  } else {
    for (i in seq_len(nrow(sample_tbl))) {
      if (!identical(sample_tbl[["format"]][i], "bigwig")) {
        bw_stop("Lazy statistics are supported for BigWig-backed samples only.")
      }
      result <- bw_stats_bigwig_file(
        sample_tbl[["file"]][i],
        sample_tbl[["sample_id"]][i],
        chrom,
        start,
        end,
        n_bins,
        stat,
        use_zoom
      )
      if (nrow(result) > 0L) {
        out_n <- out_n + 1L
        out[[out_n]] <- result
      }
    }
  }
  if (out_n < 1L) {
    return(bw_empty_stats())
  }
  data.table::rbindlist(out[seq_len(out_n)], use.names = TRUE, fill = TRUE)[]
}
