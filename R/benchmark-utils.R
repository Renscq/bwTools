# Author: Rensc
# Date: 2026-08-28
# Version: dev001
# Function: Provide internal helpers for large BigWig benchmarking
# Input: Benchmark configuration and measured R expressions
# Output: Normalized regions, timing metrics, and benchmark summaries

.BWT_BENCHMARK_OPERATIONS <- c(
  "read_lazy",
  "retrieve",
  "stats_zoom",
  "stats_full",
  "merge",
  "write",
  "read_memory"
)

bw_benchmark_clear_metadata_cache <- function(file) {
  key <- bw_cache_key(file)
  if (exists(key, envir = .bwtools_cache, inherits = FALSE)) {
    rm(list = key, envir = .bwtools_cache)
  }
  invisible(NULL)
}

bw_benchmark_positive_integer <- function(x, name, allow_zero = FALSE) {
  value <- suppressWarnings(as.numeric(x))
  lower <- if (isTRUE(allow_zero)) 0 else 1
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value < lower || value != floor(value) || value > .Machine$integer.max) {
    comparator <- if (isTRUE(allow_zero)) "non-negative" else "positive"
    bw_stop(paste0("`", name, "` must be a ", comparator, " integer scalar."))
  }
  as.integer(value)
}

bw_benchmark_region_sizes <- function(region_sizes) {
  values <- suppressWarnings(as.numeric(region_sizes))
  if (length(values) < 1L || anyNA(values) || any(!is.finite(values)) ||
      any(values < 1) || any(values != floor(values)) ||
      any(values > .Machine$integer.max)) {
    bw_stop("`region_sizes` must contain positive integer base-pair lengths within the R integer range.")
  }
  unique(as.integer(values))
}

bw_benchmark_default_regions <- function(seqinfo, chrom = NULL, region_sizes) {
  seqinfo <- data.table::copy(data.table::as.data.table(seqinfo))
  known <- which(!is.na(seqinfo[["length"]]) & seqinfo[["length"]] > 0L)
  seqinfo <- seqinfo[known]
  if (nrow(seqinfo) < 1L) {
    bw_stop("Large-file benchmarking requires known chromosome lengths.")
  }

  if (is.null(chrom)) {
    chrom_length_numeric <- as.numeric(seqinfo[["length"]])
    chrom_name <- seqinfo[["chrom"]][which.max(chrom_length_numeric)]
  } else {
    chrom <- as.character(chrom)
    if (length(chrom) != 1L || is.na(chrom) || !nzchar(chrom)) {
      bw_stop("`chrom` must be NULL or a single non-empty chromosome name.")
    }
    idx <- which(seqinfo[["chrom"]] == chrom)
    if (length(idx) < 1L) {
      bw_stop(paste0("Unknown benchmark chromosome: ", chrom, "."))
    }
    chrom_name <- chrom
  }

  idx <- which(seqinfo[["chrom"]] == chrom_name)
  chrom_length <- max(as.numeric(seqinfo[["length"]][idx]), na.rm = TRUE)
  requested <- bw_benchmark_region_sizes(region_sizes)
  actual <- pmin(as.numeric(requested), chrom_length)
  actual <- unique(as.integer(actual))

  out <- vector("list", length(actual))
  for (i in seq_along(actual)) {
    size <- actual[i]
    start <- floor((chrom_length - size) / 2) + 1
    end <- start + size - 1
    out[[i]] <- data.table::data.table(
      region_id = paste0("region_", i),
      chrom = chrom_name,
      start = as.integer(start),
      end = as.integer(end),
      bases = as.integer(size)
    )
  }
  data.table::rbindlist(out, use.names = TRUE)[]
}

bw_benchmark_explicit_regions <- function(regions, seqinfo) {
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
    bw_stop("`regions` must contain at least one benchmark interval.")
  }

  chrom_value <- as.character(x[["chrom"]])
  start_value <- suppressWarnings(as.numeric(x[["start"]]))
  end_value <- suppressWarnings(as.numeric(x[["end"]]))
  invalid <- is.na(chrom_value) | !nzchar(chrom_value) |
    !is.finite(start_value) | !is.finite(end_value) |
    start_value < 1 | end_value < start_value |
    start_value != floor(start_value) | end_value != floor(end_value) |
    start_value > .Machine$integer.max | end_value > .Machine$integer.max
  if (any(invalid)) {
    bw_stop("Benchmark regions must be valid 1-based closed intervals within the R integer range.")
  }
  data.table::set(x, j = "chrom", value = chrom_value)
  data.table::set(x, j = "start", value = as.integer(start_value))
  data.table::set(x, j = "end", value = as.integer(end_value))

  seqinfo <- data.table::copy(data.table::as.data.table(seqinfo))
  chrom_names <- unique(seqinfo[["chrom"]])
  unknown <- setdiff(unique(x[["chrom"]]), chrom_names)
  if (length(unknown) > 0L) {
    bw_stop(paste0(
      "Benchmark regions contain chromosomes absent from the BigWig file: ",
      paste(unknown, collapse = ", "),
      "."
    ))
  }

  for (i in seq_len(nrow(x))) {
    idx <- which(seqinfo[["chrom"]] == x[["chrom"]][i])
    lengths <- as.numeric(seqinfo[["length"]][idx])
    lengths <- lengths[is.finite(lengths)]
    if (length(lengths) < 1L) {
      bw_stop(paste0("Unknown chromosome length for benchmark region: ", x[["chrom"]][i], "."))
    }
    chrom_length <- max(lengths)
    if (x[["start"]][i] > chrom_length) {
      bw_stop(paste0(
        "Benchmark region starts beyond chromosome length: ",
        x[["chrom"]][i], ":", x[["start"]][i], "."
      ))
    }
    if (x[["end"]][i] > chrom_length) {
      x[["end"]][i] <- as.integer(chrom_length)
    }
  }

  region_ids <- if ("region_id" %in% names(x)) {
    as.character(x[["region_id"]])
  } else {
    paste0("region_", seq_len(nrow(x)))
  }
  if (anyNA(region_ids) || any(!nzchar(region_ids)) || anyDuplicated(region_ids)) {
    bw_stop("`regions$region_id` must contain unique non-empty values when supplied.")
  }
  data.table::set(x, j = "region_id", value = region_ids)
  data.table::set(
    x,
    j = "bases",
    value = as.integer(as.numeric(x[["end"]]) - as.numeric(x[["start"]]) + 1)
  )
  x <- x[, c("region_id", "chrom", "start", "end", "bases"), with = FALSE]
  x[]
}

bw_benchmark_gc_used_mb <- function(info) {
  if (!is.matrix(info) || ncol(info) < 2L) {
    return(NA_real_)
  }
  sum(as.numeric(info[, 2L]), na.rm = TRUE)
}

bw_benchmark_gc_peak_mb <- function(info = suppressWarnings(gc())) {
  if (!is.matrix(info) || ncol(info) < 6L) {
    return(NA_real_)
  }
  sum(as.numeric(info[, 6L]), na.rm = TRUE)
}

bw_benchmark_result_rows <- function(x) {
  if (inherits(x, "BwgTrack")) {
    return(if (is.null(x$data)) NA_integer_ else nrow(x$data))
  }
  if (is.data.frame(x)) {
    return(nrow(x))
  }
  NA_integer_
}

bw_benchmark_measure <- function(fun, warmup = 0L) {
  if (warmup > 0L) {
    for (i in seq_len(warmup)) {
      try(fun(), silent = TRUE)
    }
  }

  baseline_gc <- suppressWarnings(gc(reset = TRUE))
  baseline_heap_mb <- bw_benchmark_gc_used_mb(baseline_gc)
  before <- proc.time()
  value <- NULL
  error_message <- NA_character_
  value <- tryCatch(
    fun(),
    error = function(e) {
      error_message <<- conditionMessage(e)
      NULL
    }
  )
  elapsed <- proc.time() - before
  final_gc <- suppressWarnings(gc())
  peak_heap_mb <- bw_benchmark_gc_peak_mb(final_gc)
  heap_delta_mb <- if (
    is.finite(peak_heap_mb) && is.finite(baseline_heap_mb)
  ) max(peak_heap_mb - baseline_heap_mb, 0) else NA_real_

  list(
    value = value,
    status = if (is.na(error_message)) "ok" else "error",
    message = error_message,
    elapsed_sec = unname(elapsed[["elapsed"]]),
    user_sec = unname(elapsed[["user.self"]]),
    system_sec = unname(elapsed[["sys.self"]]),
    baseline_heap_mb = baseline_heap_mb,
    peak_heap_mb = peak_heap_mb,
    heap_delta_mb = heap_delta_mb,
    result_mb = if (is.null(value)) NA_real_ else as.numeric(object.size(value)) / 1024^2,
    result_rows = if (is.null(value)) NA_integer_ else bw_benchmark_result_rows(value)
  )
}

bw_benchmark_result_row <- function(
  measurement,
  operation,
  variant,
  iteration,
  file_mb,
  region = NULL,
  output_mb = NA_real_,
  stat_source = NA_character_,
  resolution = NA_real_,
  note = NA_character_
) {
  bases <- if (is.null(region)) NA_real_ else as.numeric(region[["bases"]])
  elapsed <- measurement$elapsed_sec
  data.table::data.table(
    operation = operation,
    variant = variant,
    region_id = if (is.null(region)) NA_character_ else as.character(region[["region_id"]]),
    chrom = if (is.null(region)) NA_character_ else as.character(region[["chrom"]]),
    start = if (is.null(region)) NA_integer_ else as.integer(region[["start"]]),
    end = if (is.null(region)) NA_integer_ else as.integer(region[["end"]]),
    bases = bases,
    iteration = as.integer(iteration),
    status = measurement$status,
    message = measurement$message,
    elapsed_sec = measurement$elapsed_sec,
    user_sec = measurement$user_sec,
    system_sec = measurement$system_sec,
    baseline_heap_mb = measurement$baseline_heap_mb,
    peak_heap_mb = measurement$peak_heap_mb,
    heap_delta_mb = measurement$heap_delta_mb,
    result_mb = measurement$result_mb,
    result_rows = measurement$result_rows,
    input_file_mb = file_mb,
    output_file_mb = output_mb,
    stat_source = stat_source,
    resolution = resolution,
    file_mb_per_sec = if (
      identical(measurement$status, "ok") && is.finite(file_mb) && elapsed > 0 &&
        operation %in% c("read_memory")
    ) file_mb / elapsed else NA_real_,
    output_mb_per_sec = if (
      identical(measurement$status, "ok") && is.finite(output_mb) && elapsed > 0
    ) output_mb / elapsed else NA_real_,
    mbases_per_sec = if (
      identical(measurement$status, "ok") && is.finite(bases) && elapsed > 0
    ) (bases / 1e6) / elapsed else NA_real_,
    note = note
  )
}

bw_benchmark_skipped_row <- function(operation, variant, file_mb, region = NULL, message) {
  measurement <- list(
    status = "skipped",
    message = message,
    elapsed_sec = NA_real_,
    user_sec = NA_real_,
    system_sec = NA_real_,
    baseline_heap_mb = NA_real_,
    peak_heap_mb = NA_real_,
    heap_delta_mb = NA_real_,
    result_mb = NA_real_,
    result_rows = NA_integer_
  )
  bw_benchmark_result_row(
    measurement,
    operation = operation,
    variant = variant,
    iteration = NA_integer_,
    file_mb = file_mb,
    region = region,
    note = message
  )
}

bw_benchmark_split_region <- function(region, parts = 4L) {
  parts <- min(as.integer(parts), as.integer(region[["bases"]]))
  if (parts < 2L) {
    return(data.table::data.table())
  }
  breaks <- floor(seq(
    as.numeric(region[["start"]]),
    as.numeric(region[["end"]]) + 1,
    length.out = parts + 1L
  ))
  out <- vector("list", parts)
  for (i in seq_len(parts)) {
    start <- as.integer(breaks[i])
    end <- as.integer(breaks[i + 1L] - 1)
    if (i == parts) end <- as.integer(region[["end"]])
    out[[i]] <- data.table::data.table(
      chrom = as.character(region[["chrom"]]),
      start = start,
      end = end
    )
  }
  data.table::rbindlist(out, use.names = TRUE)[]
}

bw_benchmark_summary <- function(results) {
  results <- data.table::copy(data.table::as.data.table(results))
  if (nrow(results) < 1L) {
    return(data.table::data.table())
  }
  keys <- unique(results[, c(
    "operation", "variant", "region_id", "chrom", "start", "end", "bases",
    "stat_source", "resolution"
  ), with = FALSE])
  out <- vector("list", nrow(keys))

  median_or_na <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 1L) NA_real_ else stats::median(x)
  }
  min_or_na <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 1L) NA_real_ else min(x)
  }
  max_or_na <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) < 1L) NA_real_ else max(x)
  }

  for (i in seq_len(nrow(keys))) {
    key <- keys[i]
    same <-
      results[["operation"]] == key[["operation"]] &
      results[["variant"]] == key[["variant"]]
    for (name in c(
      "region_id", "chrom", "start", "end", "bases", "stat_source",
      "resolution"
    )) {
      target <- key[[name]]
      column <- results[[name]]
      if (length(target) == 1L && is.na(target)) {
        same <- same & is.na(column)
      } else {
        same <- same & column == target
      }
    }
    x <- results[which(same)]
    ok <- x[which(x[["status"]] == "ok")]
    out[[i]] <- data.table::data.table(
      operation = key[["operation"]],
      variant = key[["variant"]],
      region_id = key[["region_id"]],
      chrom = key[["chrom"]],
      start = key[["start"]],
      end = key[["end"]],
      bases = key[["bases"]],
      stat_source = key[["stat_source"]],
      resolution = key[["resolution"]],
      iterations = nrow(ok),
      errors = sum(x[["status"]] == "error"),
      skipped = sum(x[["status"]] == "skipped"),
      median_elapsed_sec = median_or_na(ok[["elapsed_sec"]]),
      min_elapsed_sec = min_or_na(ok[["elapsed_sec"]]),
      max_elapsed_sec = max_or_na(ok[["elapsed_sec"]]),
      median_peak_heap_mb = median_or_na(ok[["peak_heap_mb"]]),
      median_heap_delta_mb = median_or_na(ok[["heap_delta_mb"]]),
      median_result_mb = median_or_na(ok[["result_mb"]]),
      median_result_rows = median_or_na(as.numeric(ok[["result_rows"]])),
      median_file_mb_per_sec = median_or_na(ok[["file_mb_per_sec"]]),
      median_output_mb_per_sec = median_or_na(ok[["output_mb_per_sec"]]),
      median_mbases_per_sec = median_or_na(ok[["mbases_per_sec"]])
    )
  }
  data.table::rbindlist(out, use.names = TRUE, fill = TRUE)[]
}

bw_benchmark_system_info <- function() {
  sys <- Sys.info()
  package_version_safe <- function(package) {
    namespace_version <- tryCatch(
      as.character(base::getNamespaceVersion(base::asNamespace(package))),
      error = function(e) NA_character_
    )
    if (!is.na(namespace_version)) {
      return(namespace_version)
    }
    tryCatch(
      as.character(utils::packageVersion(package)),
      error = function(e) NA_character_
    )
  }
  data.table::data.table(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"),
    bwtools_version = package_version_safe("bwTools"),
    r_version = R.version.string,
    platform = R.version$platform,
    os = unname(sys[["sysname"]] %||% NA_character_),
    os_release = unname(sys[["release"]] %||% NA_character_),
    machine = unname(sys[["machine"]] %||% NA_character_),
    host = unname(sys[["nodename"]] %||% NA_character_),
    data_table_version = package_version_safe("data.table")
  )
}

