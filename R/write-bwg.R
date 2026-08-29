# Author: Rensc
# Date: 2026-08-27
# Version: dev007
# Function: Write standardized BwgTrack objects to BigWig, WIG, or bedGraph files
# Input: BwgTrack-compatible object, output format, and chromosome sizes
# Output: One signal file per selected sample

bw_check_output_file <- function(file, overwrite = FALSE) {
  if (file.exists(file) && !isTRUE(overwrite)) {
    bw_stop(paste0("File exists: ", file))
  }
  invisible(TRUE)
}

bw_read_chrom_sizes <- function(chrom_sizes) {
  if (is.character(chrom_sizes) && length(chrom_sizes) == 1L && !is.na(chrom_sizes)) {
    file <- bw_validate_local_file(chrom_sizes)
    x <- data.table::fread(
      file,
      header = FALSE,
      select = 1:2,
      col.names = c("chrom", "length"),
      showProgress = FALSE
    )
  } else {
    x <- data.table::copy(data.table::as.data.table(chrom_sizes))
    if (ncol(x) < 2L) {
      bw_stop("`chrom_sizes` must contain at least two columns: chromosome and length.")
    }
    x <- data.table::data.table(
      chrom = as.character(x[[1L]]),
      length = suppressWarnings(as.numeric(x[[2L]]))
    )
  }
  data.table::set(x, j = "chrom", value = as.character(x[["chrom"]]))
  data.table::set(x, j = "length", value = suppressWarnings(as.numeric(x[["length"]])))
  if (nrow(x) < 1L || anyNA(x$chrom) || any(!nzchar(x$chrom))) {
    bw_stop("Chromosome sizes contain missing or empty chromosome names.")
  }
  if (anyDuplicated(x$chrom)) {
    bw_stop("Chromosome sizes contain duplicated chromosome names.")
  }
  invalid_length <- !is.finite(x$length) | x$length < 1 | x$length != floor(x$length) | x$length > 2^32 - 1
  if (any(invalid_length)) {
    bw_stop("Chromosome lengths must be positive integers representable as unsigned 32-bit values.")
  }
  data.table::set(x, j = "tid", value = seq_len(nrow(x)) - 1L)
  data.table::setcolorder(x, c("tid", "chrom", "length"))
  x[]
}

bw_resolve_chrom_sizes <- function(object, chrom_sizes, sample_id) {
  if (!is.null(chrom_sizes)) {
    return(bw_read_chrom_sizes(chrom_sizes))
  }
  if (is.null(object$seqinfo)) {
    bw_stop("`chrom_sizes` is required for BigWig output when `object$seqinfo` is unavailable.")
  }
  x <- data.table::copy(data.table::as.data.table(object$seqinfo))
  if ("sample_id" %in% names(x)) {
    sample_idx <- which(x[["sample_id"]] == sample_id)
    sample_specific <- x[sample_idx]
    if (nrow(sample_specific) > 0L) {
      x <- sample_specific
    } else {
      x <- unique(x[, c("chrom", "length"), with = FALSE])
    }
  }
  if (!all(c("chrom", "length") %in% names(x))) {
    bw_stop("`object$seqinfo` must contain `chrom` and `length` columns for BigWig output.")
  }
  if (anyNA(x$length)) {
    bw_stop("`chrom_sizes` is required because `object$seqinfo` contains unknown chromosome lengths.")
  }
  bw_read_chrom_sizes(x[, c("chrom", "length"), with = FALSE])
}

bw_validate_signal_for_write <- function(signal, chrom_sizes) {
  x <- data.table::copy(data.table::as.data.table(signal))
  required <- c("chrom", "start", "end", "value")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0L) {
    bw_stop(paste0("Signal data are missing required columns: ", paste(missing, collapse = ", "), "."))
  }
  data.table::set(x, j = "chrom", value = as.character(x[["chrom"]]))
  data.table::set(x, j = "start", value = suppressWarnings(as.numeric(x[["start"]])))
  data.table::set(x, j = "end", value = suppressWarnings(as.numeric(x[["end"]])))
  data.table::set(x, j = "value", value = suppressWarnings(as.numeric(x[["value"]])))
  invalid_interval <- is.na(x$chrom) | !nzchar(x$chrom) |
    !is.finite(x$start) | !is.finite(x$end) |
    x$start < 1 | x$end < x$start |
    x$start != floor(x$start) | x$end != floor(x$end) |
    x$end > 2^32 - 1
  if (any(invalid_interval)) {
    bw_stop("Signal data contain invalid 1-based closed intervals for BigWig output.")
  }
  if (any(!is.finite(x$value)) || any(abs(x$value) > .BWT_FLOAT32_MAX)) {
    bw_stop("BigWig output requires finite signal values representable as 32-bit floating-point numbers.")
  }

  tid <- match(x$chrom, chrom_sizes$chrom) - 1L
  if (anyNA(tid)) {
    unknown <- unique(x$chrom[is.na(tid)])
    bw_stop(paste0("Signal data contain chromosomes absent from `chrom_sizes`: ", paste(unknown, collapse = ", "), "."))
  }
  data.table::set(x, j = "tid", value = as.integer(tid))
  chrom_length <- chrom_sizes$length[tid + 1L]
  if (any(x$end > chrom_length)) {
    bad <- which(x$end > chrom_length)[1L]
    bw_stop(paste0(
      "Signal interval exceeds chromosome length: ", x$chrom[bad], ":",
      format(x$start[bad], scientific = FALSE), "-", format(x$end[bad], scientific = FALSE),
      " > ", format(chrom_length[bad], scientific = FALSE), "."
    ))
  }

  data.table::setorderv(x, c("tid", "start", "end"))
  if (nrow(x) > 1L) {
    same_chrom <- x$tid[-1L] == x$tid[-nrow(x)]
    overlaps <- same_chrom & x$start[-1L] <= x$end[-nrow(x)]
    if (any(overlaps)) {
      i <- which(overlaps)[1L] + 1L
      bw_stop(paste0(
        "Overlapping signal intervals are not supported for BigWig writing: ",
        x$chrom[i], ":", format(x$start[i - 1L], scientific = FALSE), "-",
        format(x$end[i - 1L], scientific = FALSE), " overlaps ",
        format(x$start[i], scientific = FALSE), "-", format(x$end[i], scientific = FALSE), "."
      ))
    }
  }
  x[]
}

bw_write_bedgraph_table <- function(signal, file, compress = FALSE) {
  x <- data.table::data.table(
    chrom = as.character(signal$chrom),
    start = as.numeric(signal$start) - 1,
    end = as.numeric(signal$end),
    value = as.numeric(signal$value)
  )
  data.table::fwrite(
    x,
    file,
    sep = "\t",
    col.names = FALSE,
    compress = if (isTRUE(compress)) "gzip" else "auto"
  )
  invisible(file)
}

bw_write_wig_table <- function(signal, file, compress = FALSE) {
  con <- if (isTRUE(compress)) gzfile(file, open = "wt") else base::file(file, open = "wt")
  on.exit(close(con), add = TRUE)
  x <- data.table::copy(data.table::as.data.table(signal))
  span <- as.numeric(x$end) - as.numeric(x$start) + 1
  if (any(!is.finite(span) | span < 1 | span != floor(span))) {
    bw_stop("WIG output requires valid integer interval spans.")
  }
  data.table::set(x, j = "span", value = span)
  data.table::setorderv(x, c("chrom", "start", "end"))

  chrom_levels <- unique(x$chrom)
  for (chrom_name in chrom_levels) {
    chrom_idx <- which(x[["chrom"]] == chrom_name)
    y <- x[chrom_idx]
    if (nrow(y) == 0L) next
    block_start <- c(TRUE, y$span[-1L] != y$span[-nrow(y)])
    block_id <- cumsum(block_start)
    for (id in unique(block_id)) {
      block <- y[which(block_id == id)]
      current_span <- block$span[1L]
      writeLines(paste0("variableStep chrom=", chrom_name, " span=", format(current_span, scientific = FALSE)), con)
      lines <- paste(
        format(block$start, scientific = FALSE, trim = TRUE),
        format(block$value, scientific = FALSE, trim = TRUE, digits = 15),
        sep = "\t"
      )
      writeLines(lines, con)
    }
  }
  invisible(file)
}

bw_write_lazy_copy <- function(sample_tbl, outdir, format, overwrite) {
  out <- vector("list", nrow(sample_tbl))
  for (i in seq_len(nrow(sample_tbl))) {
    sid <- sample_tbl$sample_id[i]
    src <- sample_tbl$file[i]
    src_format <- tolower(sample_tbl$format[i])
    if (!identical(src_format, format)) {
      bw_stop(paste0(
        "Lazy BwgTrack objects can only be copied when the output format matches the source format. ",
        "Sample '", sid, "' is ", src_format, " but requested ", format, "."
      ))
    }
    src <- bw_validate_local_file(src)
    compression_ext <- if (format %in% c("bedgraph", "wig") && grepl("\\.bgz$", src, ignore.case = TRUE)) {
      ".bgz"
    } else if (format %in% c("bedgraph", "wig") && grepl("\\.gz$", src, ignore.case = TRUE)) {
      ".gz"
    } else {
      ""
    }
    ext <- switch(format, bedgraph = ".bedgraph", wig = ".wig", bigwig = ".bigwig")
    dst <- file.path(outdir, paste0(sid, ext, compression_ext))
    bw_check_output_file(dst, overwrite)
    copied <- file.copy(src, dst, overwrite = overwrite)
    if (!isTRUE(copied)) {
      bw_stop(paste0("Failed to copy signal file for sample '", sid, "': ", src))
    }
    out[[i]] <- data.table::data.table(
      sample_id = sid,
      file = normalizePath(dst, winslash = "/", mustWork = TRUE),
      format = format
    )
  }
  data.table::rbindlist(out, use.names = TRUE)
}

#' Write genomic signal tracks
#'
#' @description
#' Write genomic signal tracks
#'
#' Writes a standardized `BwgTrack` to BigWig, WIG, or bedGraph files. File
#' persistence is handled only by this function; retrieval, merging, and
#' statistics never write files.
#'
#' @param x A validated `BwgTrack` object.
#' @param outdir Output directory.
#' @param format Output format: `bigwig`, `wig`, or `bedgraph`.
#' @param sample_ids Optional sample IDs. Defaults to all samples.
#' @param chrom_sizes Chromosome sizes for BigWig output. May be a two-column
#'   file or table. If omitted, complete chromosome lengths are taken from
#'   `seqinfo_bwg(x)` when available.
#' @param overwrite Whether existing output files may be replaced.
#' @param compress Whether WIG or bedGraph output should be gzip-compressed.
#'   BigWig data blocks are always compressed internally by the file format.
#' @param zoom Whether native BigWig output should include zoom summary levels.
#' @param max_zoom_levels Maximum number of BigWig zoom levels to create when
#'   `zoom = TRUE`.
#' @return Invisibly returns a data.table with `sample_id`, `file`, and `format`.
#'   Every selected in-memory sample must contain signal; validation occurs
#'   before any output file is written.
#' @details
#' BigWig output uses 0-based half-open coordinates on disk while preserving
#' the 1-based closed public coordinate contract in memory. BigWig requires
#' chromosome lengths and writes native zoom summaries by default. WIG and
#' bedGraph can optionally be gzip-compressed.
#'
#' @examples
#' \dontrun{
#' bw_file <- system.file(
#'   "extdata", "bwtools_example.bigwig",
#'   package = "bwTools", mustWork = TRUE
#' )
#' x <- read_bwg(bw_file, sample_ids = "example", mode = "memory")
#' write_bwg(
#'   x,
#'   outdir = tempfile("bwtools-write-"),
#'   format = "bigwig"
#' )
#' }
#' @seealso [read_bwg()], [retrieve_bwg()]
#' @export
write_bwg <- function(
  x,
  outdir,
  format = c("bigwig", "wig", "bedgraph"),
  sample_ids = NULL,
  chrom_sizes = NULL,
  overwrite = FALSE,
  compress = FALSE,
  zoom = TRUE,
  max_zoom_levels = 10L
) {
  bw_assert_bwg(x)
  format <- bw_match_arg(format, c("bigwig", "wig", "bedgraph"), "format")
  outdir <- bw_validate_scalar_character(outdir, "outdir")
  overwrite <- bw_validate_flag(overwrite, "overwrite")
  compress <- bw_validate_flag(compress, "compress")
  zoom <- bw_validate_flag(zoom, "zoom")
  max_zoom_levels <- bw_validate_positive_integer(
    max_zoom_levels,
    "max_zoom_levels",
    minimum = 0L
  )
  if (max_zoom_levels > 65535L) {
    bw_stop("`max_zoom_levels` must be between 0 and 65535.")
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outdir <- normalizePath(outdir, winslash = "/", mustWork = TRUE)

  sample_tbl <- bw_select_samples(x, sample_ids = sample_ids)
  if (nrow(sample_tbl) < 1L) {
    bw_stop("No samples were selected for writing.")
  }

  if (is.null(x$data)) {
    return(invisible(bw_write_lazy_copy(sample_tbl, outdir, format, overwrite)))
  }

  dt <- data.table::copy(data.table::as.data.table(x$data))
  dt <- dt[which(dt[["sample_id"]] %in% sample_tbl[["sample_id"]])]
  if (nrow(dt) < 1L) {
    bw_stop("No in-memory signal records are available for the selected samples.")
  }
  available_samples <- unique(as.character(dt[["sample_id"]]))
  missing_signal <- setdiff(
    as.character(sample_tbl[["sample_id"]]),
    available_samples
  )
  if (length(missing_signal) > 0L) {
    bw_stop(paste0(
      "Selected samples contain no in-memory signal records: ",
      paste(missing_signal, collapse = ", "),
      "."
    ))
  }

  out <- vector("list", nrow(sample_tbl))
  for (i in seq_len(nrow(sample_tbl))) {
    sid <- sample_tbl[["sample_id"]][i]
    signal <- dt[which(dt[["sample_id"]] == sid)]

    if (identical(format, "bigwig")) {
      cs <- bw_resolve_chrom_sizes(x, chrom_sizes, sid)
      signal <- bw_validate_signal_for_write(signal, cs)
      file <- file.path(outdir, paste0(sid, ".bigwig"))
      bw_check_output_file(file, overwrite)
      bw_write_bigwig_file(
        signal,
        cs,
        file,
        zoom = zoom,
        max_zoom_levels = max_zoom_levels
      )
    } else if (identical(format, "bedgraph")) {
      file <- file.path(
        outdir,
        paste0(sid, ".bedgraph", if (isTRUE(compress)) ".gz" else "")
      )
      bw_check_output_file(file, overwrite)
      signal <- data.table::copy(signal)
      data.table::setorderv(signal, c("chrom", "start", "end"))
      bw_write_bedgraph_table(signal, file, compress = compress)
    } else {
      file <- file.path(
        outdir,
        paste0(sid, ".wig", if (isTRUE(compress)) ".gz" else "")
      )
      bw_check_output_file(file, overwrite)
      signal <- data.table::copy(signal)
      data.table::setorderv(signal, c("chrom", "start", "end"))
      bw_write_wig_table(signal, file, compress = compress)
    }

    out[[i]] <- data.table::data.table(
      sample_id = sid,
      file = normalizePath(file, winslash = "/", mustWork = TRUE),
      format = format
    )
  }
  ans <- data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
  invisible(ans[])
}
