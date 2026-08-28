# Author: Rensc
# Date: 2026-08-28
# Version: dev004
# Function: Read BigWig, WIG, and bedGraph files into standardized BwgTrack objects
# Input: One or more local genomic signal files
# Output: Validated BwgTrack object

#' Read genomic signal tracks
#'
#' Reads one or more local BigWig, WIG, or bedGraph files into the standardized
#' bwTools `BwgTrack` contract.
#'
#' @param files One or more local BigWig, WIG, or bedGraph files.
#' @param format Input format. `auto` detects each file independently.
#' @param sample_ids Optional sample IDs. Defaults to input file stems.
#' @param strand Strand labels associated with input files. Values must be `+`,
#'   `-`, or `*`.
#' @param mode Loading mode. `auto` uses lazy access when all inputs are BigWig
#'   and memory mode otherwise. `lazy` is supported for BigWig inputs only.
#' @return A validated `BwgTrack` object.
#' @export
read_bwg <- function(
  files,
  format = c("auto", "bigwig", "wig", "bedgraph"),
  sample_ids = NULL,
  strand = "*",
  mode = c("auto", "memory", "lazy")
) {
  format <- match.arg(format)
  mode <- match.arg(mode)
  if (!is.character(files) || length(files) < 1L || anyNA(files)) {
    bw_stop("`files` must contain one or more local file paths.")
  }
  files <- vapply(files, bw_validate_local_file, character(1L), USE.NAMES = FALSE)
  n <- length(files)

  if (is.null(sample_ids)) {
    sample_ids <- vapply(files, bw_file_stem, character(1L), USE.NAMES = FALSE)
  } else {
    sample_ids <- bw_recycle_argument(sample_ids, n, "sample_ids")
  }
  sample_ids <- as.character(sample_ids)
  if (anyNA(sample_ids) || any(!nzchar(sample_ids))) {
    bw_stop("`sample_ids` must be non-missing and non-empty.")
  }
  if (anyDuplicated(sample_ids)) {
    bw_stop("`sample_ids` must be unique.")
  }

  strand <- as.character(bw_recycle_argument(strand, n, "strand", default = "*"))
  if (anyNA(strand) || any(!strand %in% c("+", "-", "*"))) {
    bw_stop("`strand` values must be '+', '-', or '*'.")
  }

  formats <- if (identical(format, "auto")) {
    vapply(files, detect_bwg_format, character(1L))
  } else {
    rep(format, n)
  }

  if (identical(mode, "auto")) {
    mode <- if (all(formats == "bigwig")) "lazy" else "memory"
  }
  if (identical(mode, "lazy") && any(formats != "bigwig")) {
    bw_stop("`mode = 'lazy'` is supported for BigWig inputs only.")
  }

  samples <- data.table::data.table(
    sample_id = sample_ids,
    file = files,
    format = formats,
    strand = strand,
    has_strand = strand %in% c("+", "-")
  )

  seqinfo_list <- vector("list", n)
  data_list <- if (identical(mode, "memory")) vector("list", n) else NULL

  for (i in seq_len(n)) {
    fmt <- formats[i]
    x <- NULL
    if (identical(fmt, "bigwig")) {
      seqinfo <- bw_bigwig_seqinfo(files[i])
      if (identical(mode, "memory")) {
        if (anyNA(seqinfo[["length"]])) {
          bw_stop(
            "Full-memory BigWig loading requires known chromosome lengths within the R integer range; use `mode = 'lazy'` otherwise."
          )
        }
        chrom_data <- vector("list", nrow(seqinfo))
        for (j in seq_len(nrow(seqinfo))) {
          chrom_data[[j]] <- bw_bigwig_query(
            files[i],
            seqinfo[["chrom"]][j],
            1L,
            seqinfo[["length"]][j]
          )
        }
        x <- data.table::rbindlist(chrom_data, use.names = TRUE, fill = FALSE)
      }
    } else if (identical(fmt, "bedgraph")) {
      x <- bw_read_bedgraph(files[i])
      seqinfo <- data.table::data.table(
        chrom = unique(x[["chrom"]]),
        length = NA_integer_
      )
    } else {
      x <- bw_read_wig(files[i])
      seqinfo <- data.table::data.table(
        chrom = unique(x[["chrom"]]),
        length = NA_integer_
      )
    }

    data.table::set(seqinfo, j = "sample_id", value = sample_ids[i])
    data.table::setcolorder(seqinfo, c("sample_id", "chrom", "length"))
    seqinfo_list[[i]] <- seqinfo

    if (identical(mode, "memory")) {
      data.table::set(x, j = "sample_id", value = sample_ids[i])
      data.table::set(x, j = "strand", value = strand[i])
      data.table::setcolorder(x, c("sample_id", "chrom", "start", "end", "value", "strand"))
      data_list[[i]] <- x
    }
  }

  seqinfo <- data.table::rbindlist(seqinfo_list, use.names = TRUE, fill = TRUE)
  data <- if (identical(mode, "memory")) {
    data.table::rbindlist(data_list, use.names = TRUE, fill = FALSE)
  } else {
    NULL
  }

  bwg_track(
    samples = samples,
    data = data,
    seqinfo = seqinfo,
    meta = list(
      formats = unique(formats),
      local_only = TRUE
    )
  )
}
