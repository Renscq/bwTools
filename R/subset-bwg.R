# Author: Rensc
# Date: 2026-08-28
# Version: dev001
# Function: Subset BwgTrack-compatible signal tracks by genomic intervals
# Input: BwgTrack-compatible object and one or more 1-based closed regions
# Output: In-memory BwgTrack subset or subset files

bw_normalize_subset_regions <- function(regions) {
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
      "`regions` must contain valid 1-based closed intervals within the R integer range."
    )
  }

  data.table::set(x, j = "start", value = as.integer(start_value))
  data.table::set(x, j = "end", value = as.integer(end_value))
  chrom_order <- unique(x[["chrom"]])
  data.table::set(
    x,
    j = "chrom_order",
    value = match(x[["chrom"]], chrom_order)
  )
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

bw_subset_sample_table <- function(object, samples = NULL) {
  sample_tbl <- data.table::copy(data.table::as.data.table(object$samples))
  if (is.null(samples)) {
    return(sample_tbl[])
  }

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
  sample_tbl[idx]
}

bw_subset_seqinfo <- function(object, sample_ids, regions) {
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

bw_subset_memory_data <- function(object, sample_ids, regions) {
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

bw_subset_lazy_data <- function(object, sample_tbl, regions) {
  if (any(sample_tbl[["format"]] != "bigwig")) {
    bw_stop("Lazy subset retrieval is currently supported for BigWig sources only.")
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

bw_subset_manifest <- function(sample_tbl, data, format) {
  counts <- integer(nrow(sample_tbl))
  if (nrow(data) > 0L) {
    sample_match <- match(data[["sample_id"]], sample_tbl[["sample_id"]])
    counts <- tabulate(sample_match, nbins = nrow(sample_tbl))
  }
  data.table::data.table(
    sample_id = sample_tbl[["sample_id"]],
    file = rep(NA_character_, nrow(sample_tbl)),
    format = rep(format, nrow(sample_tbl)),
    status = ifelse(counts > 0L, "ready", "empty"),
    intervals = as.integer(counts)
  )
}

#' Subset genomic signal tracks
#'
#' @description
#' Extracts the union of one or more genomic regions from a `BwgTrack` object.
#' Intervals crossing a requested boundary are clipped to that boundary. For
#' lazy BigWig inputs, only indexed BigWig blocks overlapping the requested
#' regions are decoded, so memory use scales with the selected subset rather
#' than the complete source track.
#'
#' When `outdir = NULL`, an in-memory `BwgTrack` subset is returned. When an
#' output directory is supplied, the subset is written as BigWig, WIG, or
#' bedGraph using `write_bwg()` and a file manifest is returned invisibly.
#'
#' @param object A `BwgTrack`-compatible object.
#' @param regions A data.frame or data.table containing `chrom`, `start`, and
#'   `end` columns in 1-based closed coordinates. Overlapping or directly
#'   adjacent regions are merged before extraction.
#' @param outdir Optional output directory. If `NULL`, return an in-memory
#'   subset without writing files.
#' @param format Output format when `outdir` is supplied: `bigwig`, `wig`, or
#'   `bedgraph`.
#' @param samples Optional sample IDs. By default all samples are retained.
#' @param chrom_sizes Optional chromosome sizes passed to `write_bwg()` for
#'   BigWig output. Complete lengths are otherwise taken from subset `seqinfo`
#'   when available.
#' @param overwrite Whether existing output files may be replaced.
#' @param compress Whether WIG or bedGraph output should be gzip-compressed.
#' @param zoom Whether BigWig output should include zoom summary levels.
#' @param max_zoom_levels Maximum number of BigWig zoom levels to create.
#' @param empty File-output behavior for selected samples with no overlapping
#'   signal: `error` stops before writing, while `skip` records the sample as
#'   empty in the returned manifest and writes only non-empty samples.
#' @return If `outdir = NULL`, a `BwgTrack`-compatible in-memory object.
#'   Otherwise, invisibly returns a data.table with `sample_id`, `file`,
#'   `format`, `status`, and `intervals`.
#' @export
subset_bwg <- function(
  object,
  regions,
  outdir = NULL,
  format = c("bigwig", "wig", "bedgraph"),
  samples = NULL,
  chrom_sizes = NULL,
  overwrite = FALSE,
  compress = FALSE,
  zoom = TRUE,
  max_zoom_levels = 10L,
  empty = c("error", "skip")
) {
  if (!is_bwg_track(object)) {
    bw_stop("`object` must be a BwgTrack-compatible object.")
  }
  format <- match.arg(format)
  empty <- match.arg(empty)
  normalized_regions <- bw_normalize_subset_regions(regions)
  sample_tbl <- bw_subset_sample_table(object, samples)
  sample_ids <- sample_tbl[["sample_id"]]

  data <- if (is.null(object$data)) {
    bw_subset_lazy_data(object, sample_tbl, normalized_regions)
  } else {
    bw_subset_memory_data(object, sample_ids, normalized_regions)
  }
  seqinfo <- bw_subset_seqinfo(object, sample_ids, normalized_regions)

  meta <- object$meta
  meta$source_mode <- object$meta$mode %||%
    if (is.null(object$data)) "lazy" else "memory"
  meta$mode <- "memory"
  meta$subset <- TRUE
  meta$subset_requested_regions <- nrow(data.table::as.data.table(regions))
  meta$subset_regions <- data.table::copy(normalized_regions)

  subset_track <- bw_track(
    samples = sample_tbl,
    data = data,
    seqinfo = seqinfo,
    meta = meta
  )

  if (is.null(outdir)) {
    return(subset_track)
  }
  if (!is.character(outdir) || length(outdir) != 1L || is.na(outdir) || !nzchar(outdir)) {
    bw_stop("`outdir` must be a single non-empty path or NULL.")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outdir <- normalizePath(outdir, winslash = "/", mustWork = TRUE)

  manifest <- bw_subset_manifest(sample_tbl, data, format)
  empty_samples <- manifest[["sample_id"]][manifest[["intervals"]] == 0L]
  if (length(empty_samples) > 0L && identical(empty, "error")) {
    bw_stop(paste0(
      "No signal intervals overlap the requested regions for sample(s): ",
      paste(empty_samples, collapse = ", "),
      ". Use `empty = 'skip'` to write only non-empty samples."
    ))
  }

  active_samples <- manifest[["sample_id"]][manifest[["intervals"]] > 0L]
  if (length(active_samples) > 0L) {
    written <- write_bwg(
      subset_track,
      outdir = outdir,
      format = format,
      samples = active_samples,
      chrom_sizes = chrom_sizes,
      overwrite = overwrite,
      compress = compress,
      zoom = zoom,
      max_zoom_levels = max_zoom_levels
    )
    written_idx <- match(written[["sample_id"]], manifest[["sample_id"]])
    data.table::set(
      manifest,
      i = written_idx,
      j = "file",
      value = written[["file"]]
    )
    data.table::set(
      manifest,
      i = written_idx,
      j = "status",
      value = rep("written", length(written_idx))
    )
  }

  invisible(manifest[])
}
