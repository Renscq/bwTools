# Author: Rensc
# Date: 2026-08-28
# Version: dev001
# Function: Return BigWig zoom-level metadata through a public API
# Input: Local BigWig path or BwgTrack-compatible object
# Output: data.table containing available BigWig zoom levels

bw_empty_zoominfo <- function(include_sample = FALSE) {
  out <- data.table::data.table(
    level = numeric(),
    data_offset = numeric(),
    index_offset = numeric()
  )
  if (isTRUE(include_sample)) {
    data.table::set(out, j = "sample_id", value = character())
    data.table::setcolorder(out, c("sample_id", "level", "data_offset", "index_offset"))
  }
  out
}

#' Return BigWig zoom-level metadata
#'
#' Returns the zoom reduction levels stored in a local BigWig file. For a
#' `BwgTrack` object, zoom metadata are returned for each BigWig-backed sample.
#' WIG and bedGraph sources do not contain BigWig zoom levels.
#'
#' @param x A local BigWig path or a `BwgTrack`-compatible object.
#' @return A data.table containing `level`, `data_offset`, and `index_offset`.
#'   For `BwgTrack` input, a `sample_id` column is added.
#' @export
zoominfo_bwg <- function(x) {
  if (is_bwg_track(x)) {
    sample_tbl <- data.table::copy(data.table::as.data.table(x$samples))
    if (!all(c("sample_id", "file", "format") %in% names(sample_tbl))) {
      bw_stop("`x$samples` does not contain the required BigWig source metadata.")
    }
    sample_tbl <- sample_tbl[format == "bigwig"]
    if (nrow(sample_tbl) < 1L) {
      return(bw_empty_zoominfo(include_sample = TRUE))
    }

    out <- vector("list", nrow(sample_tbl))
    out_n <- 0L
    for (i in seq_len(nrow(sample_tbl))) {
      metadata <- bw_metadata(sample_tbl$file[i])
      z <- data.table::copy(data.table::as.data.table(metadata$zoom_levels))
      if (nrow(z) < 1L) {
        next
      }
      data.table::set(z, j = "sample_id", value = sample_tbl$sample_id[i])
      data.table::setcolorder(z, c("sample_id", "level", "data_offset", "index_offset"))
      out_n <- out_n + 1L
      out[[out_n]] <- z
    }
    if (out_n < 1L) {
      return(bw_empty_zoominfo(include_sample = TRUE))
    }
    return(data.table::rbindlist(out[seq_len(out_n)], use.names = TRUE))
  }

  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    bw_stop("`x` must be a local BigWig path or a BwgTrack-compatible object.")
  }
  file <- bw_validate_local_file(x)
  if (!identical(bw_detect_format(file), "bigwig")) {
    bw_stop("`x` must refer to a BigWig file.")
  }
  metadata <- bw_metadata(file)
  z <- data.table::copy(data.table::as.data.table(metadata$zoom_levels))
  if (nrow(z) < 1L) {
    return(bw_empty_zoominfo(include_sample = FALSE))
  }
  z[]
}
