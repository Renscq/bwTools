# Author: Rensc
# Date: 2026-08-28
# Version: dev001
# Function: Standardize validation and sample selection for public bwTools APIs
# Input: BwgTrack-compatible objects and sample identifiers
# Output: Validated objects, sample tables, or standardized errors

bw_validate_sample_ids <- function(sample_ids, available = NULL, allow_null = TRUE) {
  if (is.null(sample_ids)) {
    if (isTRUE(allow_null)) return(NULL)
    bw_stop("`sample_ids` must contain at least one sample ID.")
  }

  sample_ids <- unique(as.character(sample_ids))
  if (length(sample_ids) < 1L || anyNA(sample_ids) || any(!nzchar(sample_ids))) {
    bw_stop("`sample_ids` must contain one or more non-empty sample IDs.")
  }
  if (!is.null(available)) {
    missing <- setdiff(sample_ids, as.character(available))
    if (length(missing) > 0L) {
      bw_stop(paste0("Unknown sample IDs: ", paste(missing, collapse = ", "), "."))
    }
  }
  sample_ids
}

bw_select_samples <- function(x, sample_ids = NULL) {
  sample_tbl <- data.table::copy(data.table::as.data.table(x$samples))
  selected <- bw_validate_sample_ids(
    sample_ids,
    available = sample_tbl[["sample_id"]],
    allow_null = TRUE
  )
  if (is.null(selected)) return(sample_tbl[])
  idx <- match(selected, sample_tbl[["sample_id"]])
  sample_tbl[idx][]
}

bw_contract_issues <- function(x) {
  issues <- character()

  if (!is.list(x) || !inherits(x, "BwgTrack")) {
    return("Object must inherit from `BwgTrack` and use the bwTools list contract.")
  }
  required_slots <- c("samples", "data", "seqinfo", "meta")
  missing_slots <- setdiff(required_slots, names(x))
  if (length(missing_slots) > 0L) {
    issues <- c(
      issues,
      paste0("Missing BwgTrack slots: ", paste(missing_slots, collapse = ", "), ".")
    )
    return(issues)
  }

  if (!is.data.frame(x$samples)) {
    issues <- c(issues, "`samples` must be a data.frame or data.table.")
  } else {
    samples <- data.table::as.data.table(x$samples)
    if (nrow(samples) < 1L) {
      issues <- c(issues, "`samples` must contain at least one sample.")
    }
    required <- c("sample_id", "file", "format", "strand", "has_strand")
    missing <- setdiff(required, names(samples))
    if (length(missing) > 0L) {
      issues <- c(
        issues,
        paste0("`samples` is missing required columns: ", paste(missing, collapse = ", "), ".")
      )
    } else {
      if (!is.character(samples[["sample_id"]])) {
        issues <- c(issues, "`samples$sample_id` must be a character column.")
      }
      if (!is.character(samples[["file"]])) {
        issues <- c(issues, "`samples$file` must be a character column.")
      }
      if (!is.character(samples[["format"]])) {
        issues <- c(issues, "`samples$format` must be a character column.")
      }
      if (!is.character(samples[["strand"]])) {
        issues <- c(issues, "`samples$strand` must be a character column.")
      }
      sample_ids <- as.character(samples[["sample_id"]])
      if (anyNA(sample_ids) || any(!nzchar(sample_ids))) {
        issues <- c(issues, "`samples$sample_id` must be non-missing and non-empty.")
      }
      if (anyDuplicated(sample_ids)) {
        issues <- c(issues, "`samples$sample_id` must be unique.")
      }
      strands <- as.character(samples[["strand"]])
      if (anyNA(strands) || any(!strands %in% c("+", "-", "*"))) {
        issues <- c(issues, "`samples$strand` values must be '+', '-', or '*'.")
      }
      if (!is.logical(samples[["has_strand"]]) || anyNA(samples[["has_strand"]])) {
        issues <- c(issues, "`samples$has_strand` must be a non-missing logical column.")
      } else if (any(samples[["has_strand"]] != (strands %in% c("+", "-")))) {
        issues <- c(issues, "`samples$has_strand` is inconsistent with `samples$strand`.")
      }
      formats <- tolower(as.character(samples[["format"]]))
      invalid_format <- !is.na(formats) & !formats %in% c("bigwig", "wig", "bedgraph", "memory")
      if (any(invalid_format)) {
        issues <- c(issues, "`samples$format` contains unsupported format labels.")
      }
    }
  }

  available_samples <- if (is.data.frame(x$samples) && "sample_id" %in% names(x$samples)) {
    as.character(x$samples[["sample_id"]])
  } else {
    character()
  }

  if (!is.null(x$data)) {
    if (!is.data.frame(x$data)) {
      issues <- c(issues, "`data` must be NULL, a data.frame, or a data.table.")
    } else {
      signal <- data.table::as.data.table(x$data)
      required <- c("sample_id", "chrom", "start", "end", "value", "strand")
      missing <- setdiff(required, names(signal))
      if (length(missing) > 0L) {
        issues <- c(
          issues,
          paste0("`data` is missing required columns: ", paste(missing, collapse = ", "), ".")
        )
      } else if (nrow(signal) > 0L) {
        if (!is.character(signal[["sample_id"]])) {
          issues <- c(issues, "`data$sample_id` must be a character column.")
        }
        if (!is.character(signal[["chrom"]])) {
          issues <- c(issues, "`data$chrom` must be a character column.")
        }
        if (!is.integer(signal[["start"]]) || !is.integer(signal[["end"]])) {
          issues <- c(issues, "`data$start` and `data$end` must be integer columns.")
        }
        if (!is.numeric(signal[["value"]])) {
          issues <- c(issues, "`data$value` must be a numeric column.")
        }
        if (!is.character(signal[["strand"]])) {
          issues <- c(issues, "`data$strand` must be a character column.")
        }
        sample_ids <- as.character(signal[["sample_id"]])
        chrom <- as.character(signal[["chrom"]])
        start <- suppressWarnings(as.numeric(signal[["start"]]))
        end <- suppressWarnings(as.numeric(signal[["end"]]))
        value <- suppressWarnings(as.numeric(signal[["value"]]))
        strand <- as.character(signal[["strand"]])
        if (anyNA(sample_ids) || any(!nzchar(sample_ids))) {
          issues <- c(issues, "`data$sample_id` must be non-missing and non-empty.")
        }
        unknown <- setdiff(unique(sample_ids), available_samples)
        if (length(unknown) > 0L) {
          issues <- c(
            issues,
            paste0("`data` contains unknown sample IDs: ", paste(unknown, collapse = ", "), ".")
          )
        }
        if (anyNA(chrom) || any(!nzchar(chrom))) {
          issues <- c(issues, "`data$chrom` must be non-missing and non-empty.")
        }
        invalid_interval <- !is.finite(start) | !is.finite(end) |
          start < 1 | end < start | start != floor(start) | end != floor(end) |
          start > .Machine$integer.max | end > .Machine$integer.max
        if (any(invalid_interval)) {
          issues <- c(issues, "`data` contains invalid 1-based closed intervals.")
        }
        if (any(!is.finite(value))) {
          issues <- c(issues, "`data$value` must contain finite numeric values.")
        }
        if (anyNA(strand) || any(!strand %in% c("+", "-", "*"))) {
          issues <- c(issues, "`data$strand` values must be '+', '-', or '*'.")
        }
      }
    }
  }

  if (!is.null(x$seqinfo)) {
    if (!is.data.frame(x$seqinfo)) {
      issues <- c(issues, "`seqinfo` must be NULL, a data.frame, or a data.table.")
    } else {
      seqinfo <- data.table::as.data.table(x$seqinfo)
      required <- c("sample_id", "chrom", "length")
      missing <- setdiff(required, names(seqinfo))
      if (length(missing) > 0L) {
        issues <- c(
          issues,
          paste0("`seqinfo` is missing required columns: ", paste(missing, collapse = ", "), ".")
        )
      } else if (nrow(seqinfo) > 0L) {
        if (!is.character(seqinfo[["sample_id"]])) {
          issues <- c(issues, "`seqinfo$sample_id` must be a character column.")
        }
        if (!is.character(seqinfo[["chrom"]])) {
          issues <- c(issues, "`seqinfo$chrom` must be a character column.")
        }
        if (!is.integer(seqinfo[["length"]])) {
          issues <- c(issues, "`seqinfo$length` must be an integer column.")
        }
        sample_ids <- as.character(seqinfo[["sample_id"]])
        chrom <- as.character(seqinfo[["chrom"]])
        length_value <- suppressWarnings(as.numeric(seqinfo[["length"]]))
        unknown <- setdiff(unique(sample_ids), available_samples)
        if (length(unknown) > 0L) {
          issues <- c(
            issues,
            paste0("`seqinfo` contains unknown sample IDs: ", paste(unknown, collapse = ", "), ".")
          )
        }
        if (anyNA(chrom) || any(!nzchar(chrom))) {
          issues <- c(issues, "`seqinfo$chrom` must be non-missing and non-empty.")
        }
        invalid_length <- !is.na(length_value) & (
          !is.finite(length_value) | length_value < 1 |
            length_value != floor(length_value) |
            length_value > .Machine$integer.max
        )
        if (any(invalid_length)) {
          issues <- c(issues, "`seqinfo$length` must contain positive integer lengths or NA.")
        }
        key <- paste(sample_ids, chrom, sep = "\r")
        for (k in unique(key)) {
          idx <- which(key == k)
          lengths <- unique(length_value[idx])
          lengths <- lengths[!is.na(lengths)]
          if (length(lengths) > 1L) {
            issues <- c(issues, "`seqinfo` contains conflicting lengths for a sample/chromosome pair.")
            break
          }
        }
      }
    }
  }

  if (!is.list(x$meta)) {
    issues <- c(issues, "`meta` must be a list.")
  } else {
    required_meta <- c("coordinate", "schema_version", "backend", "mode")
    missing_meta <- required_meta[vapply(required_meta, function(name) is.null(x$meta[[name]]), logical(1L))]
    if (length(missing_meta) > 0L) {
      issues <- c(
        issues,
        paste0("`meta` is missing required fields: ", paste(missing_meta, collapse = ", "), ".")
      )
    }
    coordinate <- x$meta$coordinate
    if (!is.null(coordinate) &&
        (!is.character(coordinate) || length(coordinate) != 1L ||
          is.na(coordinate) || !identical(coordinate, "1-based closed"))) {
      issues <- c(issues, "`meta$coordinate` must be the character scalar '1-based closed'.")
    }
    schema_version <- x$meta$schema_version
    if (!is.null(schema_version) &&
        (length(schema_version) != 1L || is.na(schema_version) ||
          !identical(as.character(schema_version), "2"))) {
      issues <- c(issues, "`meta$schema_version` must be the scalar value '2'.")
    }
    backend <- x$meta$backend
    if (!is.null(backend) &&
        (!is.character(backend) || length(backend) != 1L ||
          is.na(backend) || !nzchar(backend))) {
      issues <- c(issues, "`meta$backend` must be a non-empty character scalar.")
    }
    mode <- x$meta$mode
    if (!is.null(mode) &&
        (!is.character(mode) || length(mode) != 1L || is.na(mode) ||
          !mode %in% c("memory", "lazy"))) {
      issues <- c(issues, "`meta$mode` must be the character scalar 'memory' or 'lazy'.")
    }
    if (!is.null(mode)) {
      expected_mode <- if (is.null(x$data)) "lazy" else "memory"
      if (!identical(mode, expected_mode)) {
        issues <- c(issues, "`meta$mode` is inconsistent with the presence of `data`.")
      }
      if (identical(mode, "lazy") && is.data.frame(x$samples)) {
        samples <- data.table::as.data.table(x$samples)
        if (any(is.na(samples[["format"]]) | samples[["format"]] != "bigwig")) {
          issues <- c(issues, "Lazy BwgTrack objects must be backed by BigWig samples only.")
        }
        files <- as.character(samples[["file"]])
        if (anyNA(files) || any(!nzchar(files)) || any(!file.exists(files))) {
          issues <- c(issues, "Lazy BwgTrack samples must reference existing local BigWig files.")
        }
      }
    }
  }

  unique(issues)
}

bw_assert_bwg <- function(x) {
  issues <- bw_contract_issues(x)
  if (length(issues) > 0L) {
    bw_stop(paste(c("Invalid BwgTrack contract:", paste0("- ", issues)), collapse = "\n"))
  }
  invisible(x)
}

#' Validate a BwgTrack object
#'
#' Performs a strict validation of the public bwTools `BwgTrack` contract.
#' Valid objects contain standardized `samples`, `data`, `seqinfo`, and `meta`
#' slots and use 1-based closed in-memory coordinates.
#'
#' @param x A `BwgTrack` object.
#' @return `x`, invisibly. An error is raised when the contract is invalid.
#' @export
validate_bwg <- function(x) {
  bw_assert_bwg(x)
}
