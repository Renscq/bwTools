# Author: Rensc
# Date: 2026-08-27
# Version: dev002
# Function: Read BigWig, WIG, and bedGraph files into BwgTrack-compatible objects
# Input: One or more local genomic signal files
# Output: BwgTrack-compatible object

#' Read genomic signal tracks
#'
#' @param files One or more local BigWig, WIG, or bedGraph files.
#' @param format Input format. Use `auto` to detect each file independently.
#' @param sample_names Optional sample names.
#' @param strand Strand labels associated with input files.
#' @param mode `memory` loads records into memory. `lazy` is currently supported
#'   for BigWig inputs only and performs indexed region retrieval on demand.
#' @return A `BwgTrack`-compatible object.
#' @export
read_bwg <- function(
  files,
  format = c("auto", "bigwig", "wig", "bedgraph"),
  sample_names = NULL,
  strand = "*",
  mode = c("memory", "lazy")
) {
  format <- match.arg(format)
  mode <- match.arg(mode)
  if (!is.character(files) || length(files) < 1L || anyNA(files)) {
    bw_stop("`files` must contain one or more local file paths.")
  }
  files <- vapply(files, bw_validate_local_file, character(1L), USE.NAMES = FALSE)
  n <- length(files)
  sample_names <- bw_recycle_argument(sample_names, n, "sample_names", default = vapply(files, bw_file_stem, character(1L)))
  sample_names <- as.character(sample_names)
  if (anyDuplicated(sample_names)) {
    bw_stop("`sample_names` must be unique.")
  }
  strand <- as.character(bw_recycle_argument(strand, n, "strand", default = "*"))
  if (any(!strand %in% c("+", "-", "*"))) {
    bw_stop("`strand` values must be '+', '-', or '*'.")
  }
  formats <- if (identical(format, "auto")) {
    vapply(files, bw_detect_format, character(1L))
  } else {
    rep(format, n)
  }
  if (mode == "lazy" && any(formats != "bigwig")) {
    bw_stop("bwTools 0.1.x supports lazy access for BigWig files only.")
  }

  samples <- data.table::data.table(
    sample_id = sample_names,
    file = files,
    format = formats,
    strand = strand,
    has_strand = strand %in% c("+", "-")
  )

  seqinfo_list <- vector("list", n)
  data_list <- if (mode == "memory") vector("list", n) else NULL

  for (i in seq_len(n)) {
    fmt <- formats[i]
    if (fmt == "bigwig") {
      seqinfo <- bw_bigwig_seqinfo(files[i])
      if (mode == "memory") {
        if (anyNA(seqinfo$length)) {
          bw_stop(
            "Full-memory BigWig loading currently requires chromosome lengths within the R integer range; use `mode = 'lazy'` for larger chromosomes."
          )
        }
        chrom_data <- vector("list", nrow(seqinfo))
        for (j in seq_len(nrow(seqinfo))) {
          chrom_data[[j]] <- bw_bigwig_query(files[i], seqinfo$chrom[j], 1L, seqinfo$length[j])
        }
        x <- data.table::rbindlist(chrom_data, use.names = TRUE, fill = FALSE)
      }
    } else if (fmt == "bedgraph") {
      x <- bw_read_bedgraph(files[i])
      seqinfo <- data.table::data.table(chrom = unique(x$chrom), length = NA_integer_)
    } else {
      x <- bw_read_wig(files[i])
      seqinfo <- data.table::data.table(chrom = unique(x$chrom), length = NA_integer_)
    }

    data.table::set(seqinfo, j = "sample_id", value = sample_names[i])
    data.table::setcolorder(seqinfo, c("sample_id", "chrom", "length"))
    seqinfo_list[[i]] <- seqinfo

    if (mode == "memory") {
      data.table::set(x, j = "sample_id", value = sample_names[i])
      data.table::set(x, j = "strand", value = strand[i])
      data.table::setcolorder(x, c("sample_id", "chrom", "start", "end", "value", "strand"))
      data_list[[i]] <- x
    }
  }

  seqinfo <- data.table::rbindlist(seqinfo_list, use.names = TRUE, fill = TRUE)
  data <- if (mode == "memory") data.table::rbindlist(data_list, use.names = TRUE, fill = FALSE) else NULL
  bw_track(
    samples = samples,
    data = data,
    seqinfo = seqinfo,
    meta = list(
      mode = mode,
      formats = unique(formats),
      local_only = TRUE
    )
  )
}
