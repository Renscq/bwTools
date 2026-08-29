# Author: Rensc
# Date: 2026-08-29
# Version: dev003
# Function: Standardize validation, sample selection, and genomic query bounds
# Input: BwgTrack-compatible objects, sample identifiers, and genomic intervals
# Output: Validated objects, bounded queries, sample tables, or standardized errors

bw_validate_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    bw_stop(paste0("`", name, "` must be TRUE or FALSE."))
  }
  x
}

bw_match_arg <- function(x, choices, name) {
  choices <- as.character(choices)
  values <- as.character(x)
  if (length(values) > 1L && identical(values, choices)) {
    return(choices[[1L]])
  }
  if (length(values) != 1L || is.na(values) || !values %in% choices) {
    bw_stop(paste0(
      "`", name, "` must be one of: ",
      paste(choices, collapse = ", "),
      "."
    ))
  }
  values[[1L]]
}

bw_validate_scalar_character <- function(
  x,
  name,
  allow_null = FALSE,
  choices = NULL
) {
  if (is.null(x)) {
    if (isTRUE(allow_null)) return(NULL)
    bw_stop(paste0("`", name, "` must be a single non-empty character value."))
  }
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    bw_stop(paste0("`", name, "` must be a single non-empty character value."))
  }
  if (!is.null(choices) && !x %in% choices) {
    bw_stop(paste0(
      "`", name, "` must be one of: ",
      paste(choices, collapse = ", "),
      "."
    ))
  }
  x
}

bw_validate_positive_integer <- function(x, name, minimum = 1L) {
  if (length(x) != 1L || is.na(x)) {
    bw_stop(paste0("`", name, "` must be a single integer >= ", minimum, "."))
  }
  value <- suppressWarnings(as.numeric(x))
  if (!is.finite(value) || value != floor(value) || value < minimum ||
      value > .Machine$integer.max) {
    bw_stop(paste0("`", name, "` must be a single integer >= ", minimum, "."))
  }
  as.integer(value)
}

bw_validate_genomic_interval <- function(chrom, start, end) {
  chrom <- bw_validate_scalar_character(as.character(chrom), "chrom")
  start_value <- bw_validate_positive_integer(start, "start", minimum = 1L)
  end_value <- bw_validate_positive_integer(end, "end", minimum = 1L)
  if (end_value < start_value) {
    bw_stop("`end` must be greater than or equal to `start`.")
  }
  list(chrom = chrom, start = start_value, end = end_value)
}

bw_query_chrom_length <- function(x, sample_ids, chrom) {
  if (is.null(x$seqinfo)) {
    return(NA_integer_)
  }

  seqinfo <- data.table::as.data.table(x$seqinfo)
  required <- c("sample_id", "chrom", "length")
  if (nrow(seqinfo) < 1L || !all(required %in% names(seqinfo))) {
    return(NA_integer_)
  }

  selected <- unique(as.character(sample_ids))
  if (length(selected) < 1L) {
    return(NA_integer_)
  }
  idx <- which(seqinfo[["sample_id"]] %in% selected)
  info <- seqinfo[idx]
  if (nrow(info) < 1L) {
    return(NA_integer_)
  }

  length_value <- suppressWarnings(as.numeric(info[["length"]]))
  info_samples <- unique(as.character(info[["sample_id"]]))
  complete_samples <- info_samples[vapply(
    info_samples,
    function(sample_id) {
      idx <- which(info[["sample_id"]] == sample_id)
      length(idx) > 0L && all(!is.na(length_value[idx]))
    },
    logical(1L)
  )]
  if (length(complete_samples) > 0L) {
    missing <- complete_samples[!vapply(
      complete_samples,
      function(sample_id) {
        any(
          info[["sample_id"]] == sample_id &
            info[["chrom"]] == chrom
        )
      },
      logical(1L)
    )]
    if (length(missing) > 0L) {
      bw_stop(paste0(
        "Chromosome `", chrom,
        "` was not found in complete sequence metadata for selected sample(s): ",
        paste(missing, collapse = ", "),
        "."
      ))
    }
  }

  chrom_idx <- which(info[["chrom"]] == chrom & !is.na(length_value))
  if (length(chrom_idx) < 1L) {
    return(NA_integer_)
  }
  lengths <- unique(length_value[chrom_idx])
  if (length(lengths) > 1L) {
    bw_stop(paste0(
      "Selected samples have conflicting chromosome lengths for `",
      chrom,
      "`."
    ))
  }
  as.integer(lengths[[1L]])
}

bw_bound_query_interval <- function(x, sample_ids, chrom, start, end) {
  chrom_length <- bw_query_chrom_length(x, sample_ids, chrom)
  if (is.na(chrom_length)) {
    return(list(
      chrom = chrom,
      start = start,
      end = end,
      chrom_length = NA_integer_,
      empty = FALSE
    ))
  }
  if (start > chrom_length) {
    return(list(
      chrom = chrom,
      start = start,
      end = end,
      chrom_length = chrom_length,
      empty = TRUE
    ))
  }
  list(
    chrom = chrom,
    start = start,
    end = min(end, chrom_length),
    chrom_length = chrom_length,
    empty = FALSE
  )
}

bw_bound_query_regions <- function(x, sample_ids, regions) {
  regions <- data.table::as.data.table(regions)
  if (nrow(regions) < 1L) {
    return(data.table::data.table(
      chrom = character(),
      start = integer(),
      end = integer()
    ))
  }

  out <- vector("list", nrow(regions))
  out_n <- 0L
  for (i in seq_len(nrow(regions))) {
    bound <- bw_bound_query_interval(
      x,
      sample_ids = sample_ids,
      chrom = regions[["chrom"]][i],
      start = regions[["start"]][i],
      end = regions[["end"]][i]
    )
    if (isTRUE(bound$empty)) next
    out_n <- out_n + 1L
    out[[out_n]] <- data.table::data.table(
      chrom = bound$chrom,
      start = as.integer(bound$start),
      end = as.integer(bound$end)
    )
  }

  if (out_n < 1L) {
    return(data.table::data.table(
      chrom = character(),
      start = integer(),
      end = integer()
    ))
  }
  data.table::rbindlist(out[seq_len(out_n)], use.names = TRUE)[]
}

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

  if (
    is.data.frame(x$data) && is.data.frame(x$seqinfo) &&
      all(c("sample_id", "chrom", "end") %in% names(x$data)) &&
      all(c("sample_id", "chrom", "length") %in% names(x$seqinfo)) &&
      nrow(x$data) > 0L
  ) {
    signal <- data.table::as.data.table(x$data)
    seqinfo <- data.table::as.data.table(x$seqinfo)
    data_sample <- as.character(signal[["sample_id"]])
    data_chrom <- as.character(signal[["chrom"]])
    end_value <- suppressWarnings(as.numeric(signal[["end"]]))

    valid_key <- !is.na(data_sample) & nzchar(data_sample) &
      !is.na(data_chrom) & nzchar(data_chrom)
    if (all(valid_key)) {
      n_signal <- length(data_sample)
      run_start <- if (n_signal == 1L) {
        1L
      } else {
        which(c(
          TRUE,
          data_sample[-1L] != data_sample[-n_signal] |
            data_chrom[-1L] != data_chrom[-n_signal]
        ))
      }
      run_end <- c(run_start[-1L] - 1L, n_signal)
      run_key <- paste(
        data_sample[run_start],
        data_chrom[run_start],
        sep = "\r"
      )
      run_max_end <- vapply(
        seq_along(run_start),
        function(i) {
          values <- end_value[run_start[i]:run_end[i]]
          values <- values[is.finite(values)]
          if (length(values) < 1L) NA_real_ else max(values)
        },
        numeric(1L)
      )

      data_keys <- unique(run_key)
      data_max_end <- vapply(
        data_keys,
        function(key) {
          values <- run_max_end[run_key == key]
          values <- values[is.finite(values)]
          if (length(values) < 1L) NA_real_ else max(values)
        },
        numeric(1L)
      )

      seq_sample <- as.character(seqinfo[["sample_id"]])
      seq_chrom <- as.character(seqinfo[["chrom"]])
      seq_length <- suppressWarnings(as.numeric(seqinfo[["length"]]))
      valid_seq_key <- !is.na(seq_sample) & nzchar(seq_sample) &
        !is.na(seq_chrom) & nzchar(seq_chrom)
      seq_key <- paste(seq_sample, seq_chrom, sep = "\r")
      seq_keys <- unique(seq_key[valid_seq_key])

      missing_key <- data_keys[!data_keys %in% seq_keys]
      if (length(missing_key) > 0L) {
        labels <- sub("\r", "/", missing_key, fixed = TRUE)
        issues <- c(
          issues,
          paste0(
            "`data` contains sample/chromosome pairs absent from `seqinfo`: ",
            paste(labels, collapse = ", "),
            "."
          )
        )
      }

      known_idx <- valid_seq_key & !is.na(seq_length)
      if (any(known_idx)) {
        known_key <- seq_key[known_idx]
        known_length <- seq_length[known_idx]
        match_idx <- match(data_keys, known_key)
        limit <- known_length[match_idx]
        beyond <- !is.na(limit) & is.finite(data_max_end) & data_max_end > limit
        if (any(beyond)) {
          issues <- c(
            issues,
            "`data` contains intervals extending beyond `seqinfo` chromosome lengths."
          )
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
