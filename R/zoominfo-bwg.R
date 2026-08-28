# Author: Rensc
# Date: 2026-08-28
# Version: dev002
# Function: Expose BigWig zoom-level metadata through the standardized BwgTrack API
# Input: BwgTrack object
# Output: Zoom metadata table

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
  out[]
}

#' Return BigWig zoom-level metadata
#'
#' Returns stored zoom reduction levels for BigWig-backed samples in a
#' standardized `BwgTrack`. WIG, bedGraph, and memory-only samples do not have
#' BigWig zoom levels.
#'
#' @param x A validated `BwgTrack` object.
#' @param sample_ids Optional sample IDs. Defaults to all samples.
#' @return A data.table containing `sample_id`, `level`, `data_offset`, and
#'   `index_offset`.
#' @export
zoominfo_bwg <- function(x, sample_ids = NULL) {
  bw_assert_bwg(x)
  sample_tbl <- bw_select_samples(x, sample_ids = sample_ids)
  sample_tbl <- sample_tbl[which(sample_tbl[["format"]] == "bigwig")]
  if (nrow(sample_tbl) < 1L) {
    return(bw_empty_zoominfo(include_sample = TRUE))
  }

  out <- vector("list", nrow(sample_tbl))
  out_n <- 0L
  for (i in seq_len(nrow(sample_tbl))) {
    file <- sample_tbl[["file"]][i]
    if (is.na(file) || !nzchar(file) || !file.exists(file)) {
      next
    }
    metadata <- bw_metadata(file)
    z <- data.table::copy(data.table::as.data.table(metadata$zoom_levels))
    if (nrow(z) < 1L) next
    data.table::set(z, j = "sample_id", value = sample_tbl[["sample_id"]][i])
    data.table::setcolorder(z, c("sample_id", "level", "data_offset", "index_offset"))
    out_n <- out_n + 1L
    out[[out_n]] <- z
  }
  if (out_n < 1L) {
    return(bw_empty_zoominfo(include_sample = TRUE))
  }
  data.table::rbindlist(out[seq_len(out_n)], use.names = TRUE, fill = TRUE)[]
}
