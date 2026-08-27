# Author: Rensc
# Date: 2026-08-27
# Version: dev002
# Function: Read bedGraph and WIG signal files
# Input: Local text signal file
# Output: Canonical 1-based closed signal table

bw_open_text_connection <- function(file) {
  if (grepl("\\.(gz|bgz)$", tolower(file))) {
    return(gzfile(file, open = "rt"))
  }
  base::file(file, open = "rt")
}

bw_read_bedgraph <- function(file) {
  file <- bw_validate_local_file(file)
  con <- bw_open_text_connection(file)
  on.exit(close(con), add = TRUE)

  out <- list()
  out_n <- 0L
  repeat {
    lines <- readLines(con, n = 100000L, warn = FALSE)
    if (length(lines) == 0L) {
      break
    }
    lines <- trimws(lines)
    keep <- nzchar(lines) & !grepl("^(#|browser|track)", lines, ignore.case = TRUE)
    lines <- lines[keep]
    if (length(lines) == 0L) {
      next
    }

    fields <- strsplit(lines, "[[:space:]]+")
    valid <- lengths(fields) >= 4L
    if (!all(valid)) {
      bw_stop("bedGraph contains a non-empty record with fewer than four fields.")
    }

    chrom <- vapply(fields, `[[`, character(1L), 1L)
    start0 <- suppressWarnings(as.numeric(vapply(fields, `[[`, character(1L), 2L)))
    end0 <- suppressWarnings(as.numeric(vapply(fields, `[[`, character(1L), 3L)))
    value <- suppressWarnings(as.numeric(vapply(fields, `[[`, character(1L), 4L)))

    if (anyNA(start0) || anyNA(end0) || anyNA(value)) {
      bw_stop("bedGraph contains non-numeric start, end, or value fields.")
    }
    if (any(start0 < 0 | end0 <= start0)) {
      bw_stop("bedGraph contains invalid 0-based half-open intervals.")
    }
    if (any(start0 + 1 > .Machine$integer.max) || any(end0 > .Machine$integer.max)) {
      bw_stop("bedGraph coordinates exceed the supported integer range.")
    }

    out_n <- out_n + 1L
    out[[out_n]] <- data.table::data.table(
      chrom = chrom,
      start = as.integer(start0 + 1),
      end = as.integer(end0),
      value = value
    )
  }

  if (out_n == 0L) {
    return(bw_empty_signal())
  }
  ans <- data.table::rbindlist(out[seq_len(out_n)], use.names = TRUE)
  data.table::setorderv(ans, c("chrom", "start", "end"))
  ans[]
}

bw_parse_wig_header <- function(line) {
  fields <- strsplit(trimws(line), "[[:space:]]+")[[1L]]
  mode <- tolower(fields[1L])
  args <- fields[-1L]
  parsed <- list()
  for (arg in args) {
    pair <- strsplit(arg, "=", fixed = TRUE)[[1L]]
    if (length(pair) == 2L) {
      parsed[[tolower(pair[1L])]] <- pair[2L]
    }
  }
  list(mode = mode, args = parsed)
}

bw_read_wig <- function(file) {
  file <- bw_validate_local_file(file)
  con <- bw_open_text_connection(file)
  on.exit(close(con), add = TRUE)

  chrom <- NULL
  mode <- NULL
  span <- 1L
  step <- NA_integer_
  next_start <- NA_integer_
  out <- vector("list", 1024L)
  out_n <- 0L

  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) {
      break
    }
    line <- trimws(line)
    if (!nzchar(line) || grepl("^(#|browser|track)", line, ignore.case = TRUE)) {
      next
    }

    if (grepl("^(fixedStep|variableStep)\\b", line, ignore.case = TRUE)) {
      header <- bw_parse_wig_header(line)
      mode <- header$mode
      chrom <- header$args$chrom
      if (is.null(chrom) || !nzchar(chrom)) {
        bw_stop("WIG step header is missing `chrom`.")
      }
      span <- suppressWarnings(as.integer(header$args$span %||% "1"))
      if (is.na(span) || span < 1L) {
        bw_stop("WIG `span` must be a positive integer.")
      }
      if (identical(mode, "fixedstep")) {
        next_start <- suppressWarnings(as.integer(header$args$start))
        step <- suppressWarnings(as.integer(header$args$step))
        if (is.na(next_start) || next_start < 1L || is.na(step) || step < 1L) {
          bw_stop("fixedStep WIG headers require positive `start` and `step` values.")
        }
      }
      next
    }

    fields <- strsplit(line, "[[:space:]]+")[[1L]]
    if (length(fields) >= 4L) {
      start0 <- suppressWarnings(as.numeric(fields[2L]))
      end0 <- suppressWarnings(as.numeric(fields[3L]))
      value <- suppressWarnings(as.numeric(fields[4L]))
      if (is.na(start0) || is.na(end0) || is.na(value) || start0 < 0 || end0 <= start0) {
        bw_stop("Invalid bedGraph-style record in WIG input.")
      }
      record <- list(chrom = fields[1L], start = start0 + 1, end = end0, value = value)
    } else if (identical(mode, "variablestep")) {
      if (length(fields) < 2L) {
        bw_stop("variableStep WIG records require position and value.")
      }
      pos <- suppressWarnings(as.numeric(fields[1L]))
      value <- suppressWarnings(as.numeric(fields[2L]))
      if (is.na(pos) || pos < 1 || is.na(value)) {
        bw_stop("Invalid variableStep WIG record.")
      }
      record <- list(chrom = chrom, start = pos, end = pos + span - 1L, value = value)
    } else if (identical(mode, "fixedstep")) {
      value <- suppressWarnings(as.numeric(fields[1L]))
      if (is.na(value)) {
        bw_stop("Invalid fixedStep WIG value.")
      }
      record <- list(chrom = chrom, start = next_start, end = next_start + span - 1L, value = value)
      next_start <- next_start + step
    } else {
      bw_stop("WIG data were encountered before a fixedStep or variableStep header.")
    }

    if (record$start > .Machine$integer.max || record$end > .Machine$integer.max) {
      bw_stop("WIG coordinates exceed the supported integer range.")
    }
    out_n <- out_n + 1L
    if (out_n > length(out)) {
      length(out) <- length(out) * 2L
    }
    out[[out_n]] <- record
  }

  if (out_n == 0L) {
    return(bw_empty_signal())
  }
  ans <- data.table::rbindlist(out[seq_len(out_n)])
  data.table::set(ans, j = "start", value = as.integer(ans[["start"]]))
  data.table::set(ans, j = "end", value = as.integer(ans[["end"]]))
  data.table::set(ans, j = "value", value = as.numeric(ans[["value"]]))
  data.table::setorderv(ans, c("chrom", "start", "end"))
  ans[]
}
