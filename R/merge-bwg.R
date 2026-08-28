# Author: Rensc
# Date: 2026-08-28
# Version: dev001
# Function: Merge and aggregate BwgTrack-compatible genomic signal objects
# Input: Two or more BwgTrack-compatible objects
# Output: A memory-mode BwgTrack containing merged or aggregated signal

bw_collect_bwg_tracks <- function(..., tracks = NULL) {
  dots <- list(...)
  if (!is.null(tracks)) {
    if (length(dots) > 0L) {
      bw_stop("Use either `...` or `tracks`, not both.")
    }
    dots <- tracks
  }
  if (length(dots) == 1L && is.list(dots[[1L]]) && !is_bwg_track(dots[[1L]])) {
    dots <- dots[[1L]]
  }
  if (length(dots) < 2L) {
    bw_stop("At least two BwgTrack-compatible objects are required.")
  }
  valid <- vapply(dots, is_bwg_track, logical(1L))
  if (!all(valid)) {
    bw_stop("All inputs must be BwgTrack-compatible objects.")
  }
  dots
}

bw_materialize_bwg <- function(object) {
  if (!is.null(object$data)) {
    return(data.table::copy(data.table::as.data.table(object$data)))
  }
  if (is.null(object$seqinfo)) {
    bw_stop("Lazy inputs require chromosome metadata before full-track combination.")
  }

  seqinfo <- data.table::copy(data.table::as.data.table(object$seqinfo))
  if (!all(c("chrom", "length") %in% names(seqinfo))) {
    bw_stop("Lazy inputs require `chrom` and `length` columns in `seqinfo`.")
  }
  known <- !is.na(seqinfo[["length"]]) & is.finite(seqinfo[["length"]])
  if (!all(known)) {
    bw_stop("Lazy inputs require known chromosome lengths before full-track combination.")
  }

  chrom_names <- unique(as.character(seqinfo[["chrom"]]))
  region_list <- lapply(chrom_names, function(chrom_name) {
    idx <- which(seqinfo[["chrom"]] == chrom_name)
    chrom_length <- max(as.numeric(seqinfo[["length"]][idx]))
    if (chrom_length > .Machine$integer.max) {
      bw_stop(
        "Full-track combination currently requires chromosome lengths within the R integer range."
      )
    }
    data.table::data.table(
      chrom = chrom_name,
      start = 1L,
      end = as.integer(chrom_length)
    )
  })
  regions <- data.table::rbindlist(region_list, use.names = TRUE)
  retrieve_bwg(object, regions = regions, result = "data")
}

bw_track_source_table <- function(tracks) {
  out <- vector("list", length(tracks))
  for (i in seq_along(tracks)) {
    samples <- data.table::copy(data.table::as.data.table(tracks[[i]]$samples))
    data.table::set(samples, j = ".source_id", value = paste0("source", i))
    data.table::set(samples, j = ".source_index", value = i)
    out[[i]] <- samples
  }
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}

bw_bind_track_data <- function(tracks) {
  out <- vector("list", length(tracks))
  for (i in seq_along(tracks)) {
    x <- bw_materialize_bwg(tracks[[i]])
    if (nrow(x) > 0L) {
      data.table::set(x, j = ".source_id", value = paste0("source", i))
      data.table::set(x, j = ".source_index", value = i)
    }
    out[[i]] <- x
  }
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}

bw_resolve_sample_strand <- function(source_samples, sample_id) {
  idx <- which(source_samples[["sample_id"]] == sample_id)
  values <- unique(as.character(source_samples[["strand"]][idx]))
  values <- values[!is.na(values) & nzchar(values)]
  if (length(values) < 1L) {
    return("*")
  }
  if (length(values) > 1L) {
    bw_stop(paste0(
      "Sample `", sample_id,
      "` has conflicting strand metadata across inputs: ",
      paste(values, collapse = ", "), "."
    ))
  }
  values[[1L]]
}

bw_build_merged_samples <- function(tracks) {
  source_samples <- bw_track_source_table(tracks)
  sample_ids <- unique(as.character(source_samples[["sample_id"]]))
  out <- lapply(sample_ids, function(sample_id) {
    strand <- bw_resolve_sample_strand(source_samples, sample_id)
    data.table::data.table(
      sample_id = sample_id,
      file = NA_character_,
      format = "memory",
      strand = strand,
      has_strand = strand %in% c("+", "-")
    )
  })
  data.table::rbindlist(out, use.names = TRUE)
}

bw_collect_seqinfo_rows <- function(tracks) {
  out <- vector("list", length(tracks))
  for (i in seq_along(tracks)) {
    x <- tracks[[i]]$seqinfo
    if (is.null(x)) next
    x <- data.table::copy(data.table::as.data.table(x))
    if (!all(c("chrom", "length") %in% names(x))) next
    data.table::set(x, j = "chrom", value = as.character(x[["chrom"]]))
    if (anyNA(x[["chrom"]]) || any(!nzchar(x[["chrom"]]))) {
      bw_stop("Sequence metadata contain missing or empty chromosome names.")
    }
    length_value <- suppressWarnings(as.numeric(x[["length"]]))
    invalid_length <- !is.na(length_value) & (
      !is.finite(length_value) | length_value < 1 | length_value != floor(length_value)
    )
    if (any(invalid_length)) {
      bw_stop("Sequence metadata contain invalid chromosome lengths.")
    }
    data.table::set(x, j = "length", value = length_value)
    if (!"sample_id" %in% names(x)) {
      sample_ids <- as.character(tracks[[i]]$samples[["sample_id"]])
      expanded <- lapply(sample_ids, function(sample_id) {
        y <- data.table::copy(x)
        data.table::set(y, j = "sample_id", value = sample_id)
        y
      })
      x <- data.table::rbindlist(expanded, use.names = TRUE, fill = TRUE)
    } else {
      data.table::set(x, j = "sample_id", value = as.character(x[["sample_id"]]))
      unknown <- setdiff(unique(x[["sample_id"]]), tracks[[i]]$samples[["sample_id"]])
      if (length(unknown) > 0L) {
        bw_stop(paste0(
          "Sequence metadata contain unknown sample IDs: ",
          paste(unknown, collapse = ", "), "."
        ))
      }
    }
    data.table::set(x, j = ".source_id", value = paste0("source", i))
    out[[i]] <- x[, c(".source_id", "sample_id", "chrom", "length"), with = FALSE]
  }
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)
}

bw_validate_seqinfo_lengths <- function(seqinfo, key_cols) {
  if (is.null(seqinfo) || nrow(seqinfo) < 1L) {
    return(invisible(TRUE))
  }
  keys <- unique(seqinfo[, key_cols, with = FALSE])
  for (i in seq_len(nrow(keys))) {
    idx <- rep(TRUE, nrow(seqinfo))
    for (key in key_cols) {
      idx <- idx & seqinfo[[key]] == keys[[key]][i]
    }
    lengths <- unique(as.numeric(seqinfo[["length"]][idx]))
    lengths <- lengths[!is.na(lengths)]
    if (length(lengths) > 1L) {
      label <- paste(
        paste0(key_cols, "=", vapply(key_cols, function(k) keys[[k]][i], character(1L))),
        collapse = ", "
      )
      bw_stop(paste0(
        "Conflicting chromosome lengths detected for ", label, ": ",
        paste(format(lengths, scientific = FALSE), collapse = ", "), "."
      ))
    }
  }
  invisible(TRUE)
}

bw_merge_seqinfo_direct <- function(tracks, sample_ids) {
  seqinfo <- bw_collect_seqinfo_rows(tracks)
  if (nrow(seqinfo) < 1L) return(NULL)
  idx <- which(seqinfo[["sample_id"]] %in% sample_ids)
  seqinfo <- seqinfo[idx]
  bw_validate_seqinfo_lengths(seqinfo, c("chrom"))
  bw_validate_seqinfo_lengths(seqinfo, c("sample_id", "chrom"))

  keys <- unique(seqinfo[, c("sample_id", "chrom"), with = FALSE])
  out <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    idx <- which(
      seqinfo[["sample_id"]] == keys[["sample_id"]][i] &
        seqinfo[["chrom"]] == keys[["chrom"]][i]
    )
    lengths <- unique(as.numeric(seqinfo[["length"]][idx]))
    lengths <- lengths[!is.na(lengths)]
    out[[i]] <- data.table::data.table(
      sample_id = keys[["sample_id"]][i],
      chrom = keys[["chrom"]][i],
      length = if (length(lengths) == 0L) NA_integer_ else as.integer(lengths[[1L]])
    )
  }
  data.table::rbindlist(out, use.names = TRUE)
}

bw_check_direct_merge_overlaps <- function(data) {
  if (nrow(data) < 2L) return(invisible(TRUE))
  x <- data.table::copy(data)
  data.table::setorderv(x, c("sample_id", "chrom", "start", "end"))
  same_group <-
    x[["sample_id"]][-1L] == x[["sample_id"]][-nrow(x)] &
    x[["chrom"]][-1L] == x[["chrom"]][-nrow(x)]
  overlap <- same_group & x[["start"]][-1L] <= x[["end"]][-nrow(x)]
  if (any(overlap)) {
    i <- which(overlap)[1L] + 1L
    bw_stop(paste0(
      "Direct merge found overlapping intervals for sample `",
      x[["sample_id"]][i], "` at ", x[["chrom"]][i], ":",
      x[["start"]][i - 1L], "-", x[["end"]][i - 1L], " and ",
      x[["start"]][i], "-", x[["end"]][i],
      ". Use `collapse_bwg()` for arithmetic aggregation."
    ))
  }
  invisible(TRUE)
}

bw_drop_exact_signal_duplicates <- function(data) {
  key <- c("sample_id", "chrom", "start", "end", "value", "strand")
  data[!duplicated(data[, key, with = FALSE])]
}

#' Directly merge genomic signal tracks
#'
#' @description
#' Combines records from two or more `BwgTrack` objects without changing signal
#' values. Records from different samples may occupy the same genomic region.
#' Records assigned to the same sample must remain non-overlapping after the
#' merge; use `collapse_bwg()` when overlapping values require arithmetic
#' aggregation.
#'
#' Lazy BigWig inputs are materialized by indexed chromosome retrieval because
#' the returned object is a memory-mode `BwgTrack`.
#'
#' @param ... Two or more `BwgTrack`-compatible objects.
#' @param tracks Optional list of `BwgTrack`-compatible objects. Use either
#'   `...` or `tracks`.
#' @param duplicates How exact duplicate signal records should be handled:
#'   `error` rejects them; `drop` removes exact duplicates before overlap
#'   validation.
#' @return A memory-mode `BwgTrack` containing the direct union of input tracks.
#' @export
merge_bwg <- function(
  ...,
  tracks = NULL,
  duplicates = c("error", "drop")
) {
  duplicates <- match.arg(duplicates)
  inputs <- bw_collect_bwg_tracks(..., tracks = tracks)
  data <- bw_bind_track_data(inputs)
  input_intervals <- nrow(data)
  if (nrow(data) < 1L) {
    bw_stop("No signal records are available to merge.")
  }

  data.table::set(data, j = ".source_id", value = NULL)
  data.table::set(data, j = ".source_index", value = NULL)
  exact_key <- c("sample_id", "chrom", "start", "end", "value", "strand")
  duplicated_rows <- duplicated(data[, exact_key, with = FALSE])
  if (any(duplicated_rows)) {
    if (identical(duplicates, "error")) {
      bw_stop(
        "Exact duplicate signal records were found. Use `duplicates = 'drop'` to remove them explicitly."
      )
    }
    data <- bw_drop_exact_signal_duplicates(data)
  }

  bw_check_direct_merge_overlaps(data)
  samples <- bw_build_merged_samples(inputs)
  seqinfo <- bw_merge_seqinfo_direct(inputs, samples[["sample_id"]])
  data.table::setorderv(data, c("sample_id", "chrom", "start", "end"))

  source_summary <- lapply(seq_along(inputs), function(i) {
    list(
      source_id = paste0("source", i),
      samples = as.character(inputs[[i]]$samples[["sample_id"]]),
      mode = inputs[[i]]$meta$mode %||% if (is.null(inputs[[i]]$data)) "lazy" else "memory"
    )
  })

  bw_track(
    samples = samples,
    data = data,
    seqinfo = seqinfo,
    meta = list(
      mode = "memory",
      operation = "merge",
      duplicates = duplicates,
      input_tracks = length(inputs),
      input_intervals = input_intervals,
      output_intervals = nrow(data),
      sources = source_summary
    )
  )
}

bw_normalize_group_map <- function(source_samples, by, groups, output_sample) {
  sample_ids <- unique(as.character(source_samples[["sample_id"]]))
  if (!is.null(groups)) {
    if (is.character(groups) && !is.null(names(groups))) {
      group_map <- data.table::data.table(
        sample_id = names(groups),
        group_id = as.character(groups)
      )
    } else {
      group_map <- data.table::copy(data.table::as.data.table(groups))
      if (!all(c("sample_id", "group_id") %in% names(group_map))) {
        bw_stop("`groups` must contain `sample_id` and `group_id` columns.")
      }
      group_map <- group_map[, c("sample_id", "group_id"), with = FALSE]
      data.table::set(group_map, j = "sample_id", value = as.character(group_map[["sample_id"]]))
      data.table::set(group_map, j = "group_id", value = as.character(group_map[["group_id"]]))
    }
    if (anyDuplicated(group_map[["sample_id"]])) {
      bw_stop("`groups$sample_id` values must be unique.")
    }
    missing_samples <- setdiff(sample_ids, group_map[["sample_id"]])
    if (length(missing_samples) > 0L) {
      bw_stop(paste0(
        "`groups` is missing sample IDs: ",
        paste(missing_samples, collapse = ", "), "."
      ))
    }
    extra_samples <- setdiff(group_map[["sample_id"]], sample_ids)
    if (length(extra_samples) > 0L) {
      bw_stop(paste0(
        "`groups` contains unknown sample IDs: ",
        paste(extra_samples, collapse = ", "), "."
      ))
    }
    if (anyNA(group_map[["group_id"]]) || any(!nzchar(group_map[["group_id"]]))) {
      bw_stop("`groups$group_id` values must be non-missing and non-empty.")
    }
    return(group_map[])
  }

  if (identical(by, "sample")) {
    return(data.table::data.table(sample_id = sample_ids, group_id = sample_ids))
  }

  group_id <- if (is.null(output_sample)) "collapsed" else as.character(output_sample)[1L]
  if (is.na(group_id) || !nzchar(group_id)) {
    bw_stop("`output_sample` must be a non-empty sample ID.")
  }
  data.table::data.table(sample_id = sample_ids, group_id = group_id)
}

bw_validate_member_nonoverlap <- function(data) {
  if (nrow(data) < 2L) return(invisible(TRUE))
  x <- data.table::copy(data)
  data.table::setorderv(x, c(".member_id", "chrom", "start", "end"))
  same_group <-
    x[[".member_id"]][-1L] == x[[".member_id"]][-nrow(x)] &
    x[["chrom"]][-1L] == x[["chrom"]][-nrow(x)]
  overlap <- same_group & x[["start"]][-1L] <= x[["end"]][-nrow(x)]
  if (any(overlap)) {
    i <- which(overlap)[1L] + 1L
    bw_stop(paste0(
      "An individual source/sample contains overlapping intervals at ",
      x[["chrom"]][i], ":", x[["start"]][i - 1L], "-",
      x[["end"]][i - 1L], " and ", x[["start"]][i], "-",
      x[["end"]][i], ". Collapse inputs must be non-overlapping within each source/sample."
    ))
  }
  invisible(TRUE)
}

bw_collapse_group_suffix <- function(group_id, strand_value, use_suffix) {
  if (!isTRUE(use_suffix)) return(group_id)
  suffix <- switch(
    strand_value,
    "+" = "plus",
    "-" = "minus",
    "*" = "unstranded",
    "unstranded"
  )
  paste0(group_id, "_", suffix)
}

bw_collapse_atomic_chrom <- function(
  data,
  method,
  missing,
  n_members,
  drop_zero,
  merge_adjacent
) {
  if (nrow(data) < 1L) {
    return(data.table::data.table(start = integer(), end = integer(), value = numeric()))
  }

  starts <- as.numeric(data[["start"]])
  stops <- as.numeric(data[["end"]]) + 1
  values <- as.numeric(data[["value"]])
  events <- data.table::data.table(
    pos = c(starts, stops),
    delta_sum = c(values, -values),
    delta_count = c(rep.int(1L, length(values)), rep.int(-1L, length(values)))
  )
  events <- events[, list(
    delta_sum = sum(delta_sum),
    delta_count = sum(delta_count)
  ), by = "pos"]
  data.table::setorderv(events, "pos")
  if (nrow(events) < 2L) {
    return(data.table::data.table(start = integer(), end = integer(), value = numeric()))
  }

  current_sum <- 0
  current_count <- 0L
  out <- vector("list", nrow(events) - 1L)
  out_n <- 0L
  for (i in seq_len(nrow(events) - 1L)) {
    current_sum <- current_sum + events[["delta_sum"]][i]
    current_count <- current_count + events[["delta_count"]][i]
    segment_start <- events[["pos"]][i]
    segment_end <- events[["pos"]][i + 1L] - 1
    if (current_count <= 0L || segment_end < segment_start) next

    value <- if (identical(method, "sum")) {
      current_sum
    } else if (identical(missing, "zero")) {
      current_sum / n_members
    } else {
      current_sum / current_count
    }
    if (isTRUE(drop_zero) && isTRUE(all.equal(value, 0, tolerance = 1e-14))) next

    out_n <- out_n + 1L
    out[[out_n]] <- data.table::data.table(
      start = as.integer(segment_start),
      end = as.integer(segment_end),
      value = as.numeric(value)
    )
  }

  if (out_n < 1L) {
    return(data.table::data.table(start = integer(), end = integer(), value = numeric()))
  }
  ans <- data.table::rbindlist(out[seq_len(out_n)], use.names = TRUE)
  if (!isTRUE(merge_adjacent) || nrow(ans) < 2L) return(ans[])

  merged <- vector("list", nrow(ans))
  merged_n <- 1L
  merged[[1L]] <- ans[1L]
  for (i in 2L:nrow(ans)) {
    previous <- merged[[merged_n]]
    adjacent <- as.numeric(previous[["end"]][1L]) + 1 == as.numeric(ans[["start"]][i])
    same_value <- isTRUE(all.equal(
      previous[["value"]][1L],
      ans[["value"]][i],
      tolerance = 1e-12
    ))
    if (adjacent && same_value) {
      previous[["end"]][1L] <- ans[["end"]][i]
      merged[[merged_n]] <- previous
    } else {
      merged_n <- merged_n + 1L
      merged[[merged_n]] <- ans[i]
    }
  }
  data.table::rbindlist(merged[seq_len(merged_n)], use.names = TRUE)
}

bw_build_collapse_seqinfo <- function(tracks, member_map) {
  seqinfo <- bw_collect_seqinfo_rows(tracks)
  if (nrow(seqinfo) < 1L) return(NULL)
  bw_validate_seqinfo_lengths(seqinfo, c("chrom"))
  seqinfo <- merge(
    seqinfo,
    member_map[, c(".source_id", "sample_id", "output_sample"), with = FALSE],
    by = c(".source_id", "sample_id"),
    all.x = FALSE,
    all.y = FALSE,
    sort = FALSE
  )
  if (nrow(seqinfo) < 1L) return(NULL)
  bw_validate_seqinfo_lengths(seqinfo, c("output_sample", "chrom"))

  keys <- unique(seqinfo[, c("output_sample", "chrom"), with = FALSE])
  out <- vector("list", nrow(keys))
  for (i in seq_len(nrow(keys))) {
    idx <- which(
      seqinfo[["output_sample"]] == keys[["output_sample"]][i] &
        seqinfo[["chrom"]] == keys[["chrom"]][i]
    )
    lengths <- unique(as.numeric(seqinfo[["length"]][idx]))
    lengths <- lengths[!is.na(lengths)]
    out[[i]] <- data.table::data.table(
      sample_id = keys[["output_sample"]][i],
      chrom = keys[["chrom"]][i],
      length = if (length(lengths) == 0L) NA_integer_ else as.integer(lengths[[1L]])
    )
  }
  data.table::rbindlist(out, use.names = TRUE)
}

#' Aggregate overlapping genomic signal tracks
#'
#' @description
#' Collapses overlapping signal from two or more `BwgTrack` objects using
#' arithmetic aggregation. Partially overlapping intervals are split at every
#' input boundary before values are calculated, so aggregation does not require
#' identical input interval coordinates.
#'
#' `by = "sample"` aggregates repeated instances of the same sample across
#' input tracks and preserves the sample ID. `by = "all"` aggregates all input
#' samples into one output sample. `groups` provides an explicit sample-to-group
#' mapping for treatment- or replicate-level aggregation.
#'
#' For means, `missing = "zero"` treats uncovered members as zero, whereas
#' `missing = "ignore"` averages only members that cover a given atomic segment.
#'
#' @param ... Two or more `BwgTrack`-compatible objects.
#' @param tracks Optional list of `BwgTrack`-compatible objects. Use either
#'   `...` or `tracks`.
#' @param method Aggregation method: `sum` or `mean`.
#' @param by Default grouping rule when `groups` is not supplied: `sample`
#'   preserves sample IDs; `all` combines all samples.
#' @param groups Optional named character vector or table with `sample_id` and
#'   `group_id` columns.
#' @param output_sample Output sample ID used when `by = "all"`. Defaults to
#'   `collapsed`.
#' @param missing Mean denominator rule: `zero` treats uncovered members as
#'   zero; `ignore` excludes uncovered members. Ignored for `method = "sum"`.
#' @param strand Strand aggregation rule. `separate` keeps `+`, `-`, and `*`
#'   signals separate; `ignore` combines strands.
#' @param drop_zero Whether zero-valued output segments should be omitted.
#' @param merge_adjacent Whether directly adjacent output segments with equal
#'   values should be merged.
#' @return A memory-mode `BwgTrack` containing aggregated signal.
#' @export
collapse_bwg <- function(
  ...,
  tracks = NULL,
  method = c("sum", "mean"),
  by = c("sample", "all"),
  groups = NULL,
  output_sample = NULL,
  missing = c("zero", "ignore"),
  strand = c("separate", "ignore"),
  drop_zero = TRUE,
  merge_adjacent = TRUE
) {
  method <- match.arg(method)
  by <- match.arg(by)
  missing <- match.arg(missing)
  strand <- match.arg(strand)
  inputs <- bw_collect_bwg_tracks(..., tracks = tracks)
  source_samples <- bw_track_source_table(inputs)
  group_map <- bw_normalize_group_map(source_samples, by, groups, output_sample)

  member_map <- unique(source_samples[, c(
    ".source_id", "sample_id", "strand"
  ), with = FALSE])
  data.table::set(
    member_map,
    j = ".member_id",
    value = paste(member_map[[".source_id"]], member_map[["sample_id"]], sep = "::")
  )
  member_group_idx <- match(member_map[["sample_id"]], group_map[["sample_id"]])
  data.table::set(
    member_map,
    j = ".group_id",
    value = group_map[["group_id"]][member_group_idx]
  )
  data.table::set(
    member_map,
    j = ".strand_group",
    value = if (identical(strand, "ignore")) "*" else as.character(member_map[["strand"]])
  )

  strand_counts <- member_map[, list(
    n_strands = data.table::uniqueN(.strand_group)
  ), by = ".group_id"]
  member_strand_idx <- match(member_map[[".group_id"]], strand_counts[[".group_id"]])
  data.table::set(
    member_map,
    j = "n_strands",
    value = strand_counts[["n_strands"]][member_strand_idx]
  )
  data.table::set(
    member_map,
    j = "output_sample",
    value = mapply(
      bw_collapse_group_suffix,
      group_id = member_map[[".group_id"]],
      strand_value = member_map[[".strand_group"]],
      use_suffix = member_map[["n_strands"]] > 1L,
      USE.NAMES = FALSE
    )
  )
  output_key <- unique(
    member_map[, c(".group_id", ".strand_group", "output_sample"), with = FALSE]
  )
  if (anyDuplicated(output_key[["output_sample"]])) {
    bw_stop(
      "Generated output sample IDs are not unique across aggregation groups. Rename group IDs to avoid strand-suffix collisions."
    )
  }
  member_counts <- member_map[, list(
    n_members = data.table::uniqueN(.member_id)
  ), by = "output_sample"]

  data <- bw_bind_track_data(inputs)
  input_intervals <- nrow(data)
  if (nrow(data) < 1L) {
    bw_stop("No signal records are available to collapse.")
  }
  if (any(!is.finite(as.numeric(data[["value"]])))) {
    bw_stop("Arithmetic collapse requires finite signal values.")
  }
  data.table::set(
    data,
    j = ".member_id",
    value = paste(data[[".source_id"]], data[["sample_id"]], sep = "::")
  )
  bw_validate_member_nonoverlap(data)
  member_idx <- match(data[[".member_id"]], member_map[[".member_id"]])
  if (anyNA(member_idx)) {
    bw_stop("Internal grouping error: one or more signal members are not mapped.")
  }
  data.table::set(data, j = ".group_id", value = member_map[[".group_id"]][member_idx])
  data.table::set(
    data,
    j = ".strand_group",
    value = member_map[[".strand_group"]][member_idx]
  )
  data.table::set(
    data,
    j = "output_sample",
    value = member_map[["output_sample"]][member_idx]
  )

  keys <- unique(data[, c("output_sample", ".strand_group", "chrom"), with = FALSE])
  out <- vector("list", nrow(keys))
  out_n <- 0L
  for (i in seq_len(nrow(keys))) {
    idx <- which(
      data[["output_sample"]] == keys[["output_sample"]][i] &
        data[[".strand_group"]] == keys[[".strand_group"]][i] &
        data[["chrom"]] == keys[["chrom"]][i]
    )
    x <- data[idx]
    member_count_idx <- match(keys[["output_sample"]][i], member_counts[["output_sample"]])
    collapsed <- bw_collapse_atomic_chrom(
      x,
      method = method,
      missing = missing,
      n_members = member_counts[["n_members"]][member_count_idx],
      drop_zero = drop_zero,
      merge_adjacent = merge_adjacent
    )
    if (nrow(collapsed) < 1L) next
    data.table::set(collapsed, j = "sample_id", value = keys[["output_sample"]][i])
    data.table::set(collapsed, j = "chrom", value = keys[["chrom"]][i])
    data.table::set(collapsed, j = "strand", value = keys[[".strand_group"]][i])
    data.table::setcolorder(
      collapsed,
      c("sample_id", "chrom", "start", "end", "value", "strand")
    )
    out_n <- out_n + 1L
    out[[out_n]] <- collapsed
  }
  result_data <- if (out_n < 1L) {
    bw_empty_signal(include_sample = TRUE)
  } else {
    x <- data.table::rbindlist(out[seq_len(out_n)], use.names = TRUE)
    data.table::setorderv(x, c("sample_id", "chrom", "start", "end"))
    x[]
  }

  output_map <- unique(member_map[, c(
    ".source_id", "sample_id", ".group_id", ".strand_group", "output_sample"
  ), with = FALSE])
  output_ids <- unique(output_map[["output_sample"]])
  samples <- lapply(output_ids, function(sample_id) {
    idx <- which(output_map[["output_sample"]] == sample_id)
    strand_values <- unique(output_map[[".strand_group"]][idx])
    strand_value <- if (length(strand_values) == 1L) strand_values[[1L]] else "*"
    data.table::data.table(
      sample_id = sample_id,
      file = NA_character_,
      format = "memory",
      strand = strand_value,
      has_strand = strand_value %in% c("+", "-")
    )
  })
  samples <- data.table::rbindlist(samples, use.names = TRUE)
  seqinfo <- bw_build_collapse_seqinfo(inputs, output_map)
  if (!is.null(seqinfo)) {
    seqinfo <- seqinfo[seqinfo[["sample_id"]] %in% samples[["sample_id"]]]
  }

  provenance <- output_map[, list(
    members = list(unique(paste(.source_id, sample_id, sep = "::")))
  ), by = "output_sample"]

  bw_track(
    samples = samples,
    data = result_data,
    seqinfo = seqinfo,
    meta = list(
      mode = "memory",
      operation = "collapse",
      method = method,
      grouping = if (is.null(groups)) by else "groups",
      input_tracks = length(inputs),
      input_intervals = input_intervals,
      output_intervals = nrow(result_data),
      missing = if (identical(method, "mean")) missing else NA_character_,
      strand = strand,
      drop_zero = isTRUE(drop_zero),
      merge_adjacent = isTRUE(merge_adjacent),
      provenance = provenance
    )
  )
}
