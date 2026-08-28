# Author: Rensc
# Date: 2026-08-27
# Version: dev003
# Function: Calculate exact or zoom-accelerated BigWig interval statistics
# Input: Validated BwgTrack object and genomic bins
# Output: Per-bin signal statistics in 1-based closed coordinates

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

bw_finish_stats <- function(bins, covered, sum_data, sum_squared, min_value, max_value, stat, source, resolution) {
  has_data <- covered > 0
  mean_value <- rep(NA_real_, length(covered))
  mean_value[has_data] <- sum_data[has_data] / covered[has_data]

  stdev_value <- rep(NA_real_, length(covered))
  one <- covered == 1
  many <- covered > 1
  stdev_value[one] <- 0
  if (any(many)) {
    variance <- (sum_squared[many] - (sum_data[many]^2 / covered[many])) / (covered[many] - 1)
    variance[variance < 0 & abs(variance) < 1e-10] <- 0
    variance[variance < 0] <- NA_real_
    stdev_value[many] <- sqrt(variance)
  }

  value <- switch(
    stat,
    mean = mean_value,
    stdev = stdev_value,
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

bw_stats_from_signal <- function(signal, chrom, start, end, n_bins, stat) {
  bins <- bw_stat_bins(start, end, n_bins)
  n <- nrow(bins)
  covered <- numeric(n)
  sum_data <- numeric(n)
  sum_squared <- numeric(n)
  min_value <- rep(Inf, n)
  max_value <- rep(-Inf, n)

  x <- data.table::as.data.table(signal)
  if (nrow(x) > 0L) {
    start0 <- pmax(as.numeric(x$start) - 1, as.numeric(start) - 1)
    end0 <- pmin(as.numeric(x$end), as.numeric(end))
    value <- as.numeric(x$value)
    keep <- end0 > start0 & is.finite(value)
    start0 <- start0[keep]
    end0 <- end0[keep]
    value <- value[keep]

    for (i in seq_along(start0)) {
      pos <- start0[i]
      while (pos < end0[i]) {
        j <- findInterval(pos, bins$end0) + 1L
        if (j < 1L) j <- 1L
        if (j > n) break
        segment_end <- min(end0[i], bins$end0[j])
        w <- segment_end - pos
        if (w <= 0) break
        covered[j] <- covered[j] + w
        sum_data[j] <- sum_data[j] + w * value[i]
        sum_squared[j] <- sum_squared[j] + w * value[i]^2
        min_value[j] <- min(min_value[j], value[i])
        max_value[j] <- max(max_value[j], value[i])
        pos <- segment_end
      }
    }
  }

  out <- bw_finish_stats(
    bins, covered, sum_data, sum_squared, min_value, max_value,
    stat = stat, source = "full", resolution = NA_integer_
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
    signal <- bw_bigwig_query(file, chrom, start, end)
    out <- bw_stats_from_signal(signal, chrom, start, end, n_bins, stat)
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
  stat <- match.arg(stat)
  if (!is.logical(use_zoom) || length(use_zoom) != 1L || is.na(use_zoom)) {
    bw_stop("`use_zoom` must be TRUE or FALSE.")
  }
  chrom <- as.character(chrom)
  if (length(chrom) != 1L || is.na(chrom) || !nzchar(chrom)) {
    bw_stop("`chrom` must be a single non-empty chromosome name.")
  }
  start <- suppressWarnings(as.integer(start))
  end <- suppressWarnings(as.integer(end))
  n_bins <- suppressWarnings(as.integer(n_bins))
  if (length(start) != 1L || length(end) != 1L || length(n_bins) != 1L) {
    bw_stop("`start`, `end`, and `n_bins` must be scalar values.")
  }
  bw_stat_bins(start, end, n_bins)

  sample_tbl <- bw_select_samples(x, sample_ids = sample_ids)
  if (nrow(sample_tbl) < 1L) {
    return(data.table::data.table())
  }

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
    return(data.table::data.table())
  }
  data.table::rbindlist(out[seq_len(out_n)], use.names = TRUE, fill = TRUE)[]
}
