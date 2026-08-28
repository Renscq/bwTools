# Author: Rensc
# Date: 2026-08-28
# Version: dev001
# Function: Run structured large BigWig performance benchmarks
# Input: One local BigWig file and benchmark configuration
# Output: bwToolsBenchmark object with raw and summarized metrics

#' Benchmark a large BigWig workflow
#'
#' Profiles the native bwTools BigWig workflow using deterministic genomic
#' regions. The standard suite measures package-metadata-cold and cached lazy
#' loading, indexed retrieval, zoom statistics, and bounded exact statistics.
#' Optional operations benchmark region materialization/merge, native BigWig
#' writing, or full-memory loading.
#'
#' This diagnostic function is intentionally file-oriented because its purpose
#' is to measure file I/O. Normal analysis functions continue to consume
#' validated `BwgTrack` objects.
#'
#' @param file One local BigWig file.
#' @param sample_id Sample identifier used for the benchmark track.
#' @param regions Optional explicit genomic regions with `chrom`, `start`, and
#'   `end` columns. An optional unique `region_id` column is preserved. When
#'   omitted, centered windows are generated from
#'   `region_sizes` on `chrom` or on the largest chromosome.
#' @param chrom Optional chromosome used for automatically generated regions.
#' @param region_sizes Region sizes in base pairs when `regions` is omitted.
#' @param n_bins Number of bins used by statistics benchmarks.
#' @param stat Statistic passed to `stats_bwg()`.
#' @param iterations Number of timed iterations per benchmark case.
#' @param warmup Number of untimed warmup iterations per benchmark case.
#' @param operations Operations to run. Supported values are `read_lazy`,
#'   `retrieve`, `stats_zoom`, `stats_full`, `merge`, `write`, and
#'   `read_memory`.
#' @param full_stats_max_bases Maximum region size for `stats_full`. Larger
#'   regions are recorded as skipped rather than forcing an expensive exact
#'   calculation.
#' @return A `bwToolsBenchmark` object containing `system`, `file`, `config`,
#'   `regions`, `results`, and `summary` components. Peak memory is an
#'   approximate R heap maximum from `gc()` rather than operating-system RSS.
#' @export
benchmark_bwg <- function(
  file,
  sample_id = "benchmark",
  regions = NULL,
  chrom = NULL,
  region_sizes = c(1e4, 1e6, 1e7),
  n_bins = 1000L,
  stat = c("mean", "stdev", "max", "min", "coverage", "sum"),
  iterations = 3L,
  warmup = 1L,
  operations = c("read_lazy", "retrieve", "stats_zoom", "stats_full"),
  full_stats_max_bases = 1e6
) {
  file <- bw_validate_local_file(file)
  if (!identical(detect_bwg_format(file), "bigwig")) {
    bw_stop("`benchmark_bwg()` currently supports BigWig input only.")
  }
  sample_id <- as.character(sample_id)
  if (length(sample_id) != 1L || is.na(sample_id) || !nzchar(sample_id)) {
    bw_stop("`sample_id` must be a single non-empty value.")
  }
  stat <- match.arg(stat)
  n_bins <- bw_benchmark_positive_integer(n_bins, "n_bins")
  iterations <- bw_benchmark_positive_integer(iterations, "iterations")
  warmup <- bw_benchmark_positive_integer(warmup, "warmup", allow_zero = TRUE)
  full_stats_max_bases <- bw_benchmark_positive_integer(
    full_stats_max_bases,
    "full_stats_max_bases"
  )

  operations <- unique(as.character(operations))
  if (length(operations) < 1L || anyNA(operations) || any(!nzchar(operations))) {
    bw_stop("`operations` must contain one or more benchmark operation names.")
  }
  unknown_operations <- setdiff(operations, .BWT_BENCHMARK_OPERATIONS)
  if (length(unknown_operations) > 0L) {
    bw_stop(paste0(
      "Unknown benchmark operations: ",
      paste(unknown_operations, collapse = ", "),
      "."
    ))
  }

  file_mb <- as.numeric(file.info(file)$size) / 1024^2
  lazy_track <- read_bwg(
    file,
    sample_ids = sample_id,
    mode = "lazy"
  )
  seqinfo <- seqinfo_bwg(lazy_track, sample_ids = sample_id)
  benchmark_regions <- if (is.null(regions)) {
    bw_benchmark_default_regions(
      seqinfo,
      chrom = chrom,
      region_sizes = region_sizes
    )
  } else {
    if (!is.null(chrom)) {
      bw_stop("`chrom` is used only when `regions` is NULL.")
    }
    bw_benchmark_explicit_regions(regions, seqinfo)
  }

  rows <- list()
  row_n <- 0L
  add_row <- function(x) {
    row_n <<- row_n + 1L
    rows[[row_n]] <<- x
    invisible(NULL)
  }

  if ("read_lazy" %in% operations) {
    for (iteration in seq_len(iterations)) {
      bw_benchmark_clear_metadata_cache(file)
      measurement <- bw_benchmark_measure(
        function() read_bwg(file, sample_ids = sample_id, mode = "lazy"),
        warmup = 0L
      )
      add_row(bw_benchmark_result_row(
        measurement,
        operation = "read_lazy",
        variant = "metadata_cache_cold",
        iteration = iteration,
        file_mb = file_mb,
        note = "Cold refers to the bwTools metadata cache only; operating-system file caches are not cleared."
      ))
    }

    invisible(read_bwg(file, sample_ids = sample_id, mode = "lazy"))
    for (iteration in seq_len(iterations)) {
      measurement <- bw_benchmark_measure(
        function() read_bwg(file, sample_ids = sample_id, mode = "lazy"),
        warmup = warmup
      )
      add_row(bw_benchmark_result_row(
        measurement,
        operation = "read_lazy",
        variant = "metadata_cache_warm",
        iteration = iteration,
        file_mb = file_mb
      ))
    }
  }

  if ("retrieve" %in% operations) {
    for (region_i in seq_len(nrow(benchmark_regions))) {
      region <- benchmark_regions[region_i]
      for (iteration in seq_len(iterations)) {
        measurement <- bw_benchmark_measure(
          function() retrieve_bwg(
            lazy_track,
            chrom = region[["chrom"]],
            start = region[["start"]],
            end = region[["end"]],
            sample_ids = sample_id,
            result = "data"
          ),
          warmup = warmup
        )
        add_row(bw_benchmark_result_row(
          measurement,
          operation = "retrieve",
          variant = "indexed",
          iteration = iteration,
          file_mb = file_mb,
          region = region
        ))
      }
    }
  }

  if ("stats_zoom" %in% operations) {
    for (region_i in seq_len(nrow(benchmark_regions))) {
      region <- benchmark_regions[region_i]
      for (iteration in seq_len(iterations)) {
        measurement <- bw_benchmark_measure(
          function() stats_bwg(
            lazy_track,
            chrom = region[["chrom"]],
            start = region[["start"]],
            end = region[["end"]],
            n_bins = n_bins,
            stat = stat,
            sample_ids = sample_id,
            use_zoom = TRUE
          ),
          warmup = warmup
        )
        stat_source <- NA_character_
        resolution <- NA_real_
        if (!is.null(measurement$value) && is.data.frame(measurement$value)) {
          if ("source" %in% names(measurement$value)) {
            stat_source <- paste(unique(measurement$value[["source"]]), collapse = ",")
          }
          if ("resolution" %in% names(measurement$value)) {
            resolution_values <- unique(as.numeric(measurement$value[["resolution"]]))
            resolution_values <- resolution_values[is.finite(resolution_values)]
            if (length(resolution_values) == 1L) {
              resolution <- resolution_values[[1L]]
            }
          }
        }
        add_row(bw_benchmark_result_row(
          measurement,
          operation = "stats_zoom",
          variant = stat,
          iteration = iteration,
          file_mb = file_mb,
          region = region,
          stat_source = stat_source,
          resolution = resolution
        ))
      }
    }
  }

  if ("stats_full" %in% operations) {
    for (region_i in seq_len(nrow(benchmark_regions))) {
      region <- benchmark_regions[region_i]
      if (as.numeric(region[["bases"]]) > full_stats_max_bases) {
        add_row(bw_benchmark_skipped_row(
          operation = "stats_full",
          variant = stat,
          file_mb = file_mb,
          region = region,
          message = paste0(
            "Region exceeds `full_stats_max_bases` (",
            format(full_stats_max_bases, scientific = FALSE),
            " bp)."
          )
        ))
        next
      }
      for (iteration in seq_len(iterations)) {
        measurement <- bw_benchmark_measure(
          function() stats_bwg(
            lazy_track,
            chrom = region[["chrom"]],
            start = region[["start"]],
            end = region[["end"]],
            n_bins = n_bins,
            stat = stat,
            sample_ids = sample_id,
            use_zoom = FALSE
          ),
          warmup = warmup
        )
        add_row(bw_benchmark_result_row(
          measurement,
          operation = "stats_full",
          variant = stat,
          iteration = iteration,
          file_mb = file_mb,
          region = region,
          stat_source = "full",
          resolution = NA_real_
        ))
      }
    }
  }

  if ("merge" %in% operations) {
    largest <- benchmark_regions[which.max(as.numeric(benchmark_regions[["bases"]]))]
    split_regions <- bw_benchmark_split_region(largest, parts = 4L)
    if (nrow(split_regions) < 2L) {
      add_row(bw_benchmark_skipped_row(
        operation = "merge",
        variant = "auto",
        file_mb = file_mb,
        region = largest,
        message = "The largest benchmark region is too small to split for merge profiling."
      ))
    } else {
      merge_tracks <- lapply(seq_len(nrow(split_regions)), function(i) {
        retrieve_bwg(
          lazy_track,
          chrom = split_regions[["chrom"]][i],
          start = split_regions[["start"]][i],
          end = split_regions[["end"]][i],
          sample_ids = sample_id,
          result = "track"
        )
      })
      for (iteration in seq_len(iterations)) {
        measurement <- bw_benchmark_measure(
          function() merge_bwg(tracks = merge_tracks, method = "auto"),
          warmup = warmup
        )
        add_row(bw_benchmark_result_row(
          measurement,
          operation = "merge",
          variant = "auto_four_parts",
          iteration = iteration,
          file_mb = file_mb,
          region = largest
        ))
      }
    }
  }

  if ("write" %in% operations) {
    largest <- benchmark_regions[which.max(as.numeric(benchmark_regions[["bases"]]))]
    write_track <- retrieve_bwg(
      lazy_track,
      chrom = largest[["chrom"]],
      start = largest[["start"]],
      end = largest[["end"]],
      sample_ids = sample_id,
      result = "track"
    )
    for (zoom_value in c(FALSE, TRUE)) {
      variant <- if (isTRUE(zoom_value)) "zoom_on" else "zoom_off"
      for (iteration in seq_len(iterations)) {
        run_dir <- tempfile(paste0("bwtools_benchmark_write_", variant, "_"))
        dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
        measurement <- bw_benchmark_measure(
          function() write_bwg(
            write_track,
            outdir = run_dir,
            format = "bigwig",
            sample_ids = sample_id,
            overwrite = TRUE,
            zoom = zoom_value
          ),
          warmup = warmup
        )
        output_mb <- NA_real_
        if (!is.null(measurement$value) && is.data.frame(measurement$value) &&
            "file" %in% names(measurement$value) && nrow(measurement$value) > 0L) {
          output_file <- measurement$value[["file"]][1L]
          if (file.exists(output_file)) {
            output_mb <- as.numeric(file.info(output_file)$size) / 1024^2
          }
        }
        add_row(bw_benchmark_result_row(
          measurement,
          operation = "write",
          variant = variant,
          iteration = iteration,
          file_mb = file_mb,
          region = largest,
          output_mb = output_mb,
          note = "Write profiling materializes only the largest selected benchmark region."
        ))
        unlink(run_dir, recursive = TRUE, force = TRUE)
      }
    }
  }

  if ("read_memory" %in% operations) {
    for (iteration in seq_len(iterations)) {
      bw_benchmark_clear_metadata_cache(file)
      measurement <- bw_benchmark_measure(
        function() read_bwg(file, sample_ids = sample_id, mode = "memory"),
        warmup = warmup
      )
      add_row(bw_benchmark_result_row(
        measurement,
        operation = "read_memory",
        variant = "full_file",
        iteration = iteration,
        file_mb = file_mb,
        note = "Full-memory loading may require substantially more memory than the compressed BigWig file size."
      ))
    }
  }

  results <- if (row_n < 1L) {
    data.table::data.table()
  } else {
    data.table::rbindlist(rows[seq_len(row_n)], use.names = TRUE, fill = TRUE)
  }
  summary <- bw_benchmark_summary(results)

  file_info <- data.table::data.table(
    file = file,
    size_mb = file_mb,
    format = "bigwig",
    sample_id = sample_id,
    chromosomes = data.table::uniqueN(seqinfo[["chrom"]]),
    total_genome_bases = sum(as.numeric(seqinfo[["length"]]), na.rm = TRUE)
  )
  config <- list(
    n_bins = n_bins,
    stat = stat,
    iterations = iterations,
    warmup = warmup,
    operations = operations,
    full_stats_max_bases = full_stats_max_bases,
    memory_metric = "approximate R heap max used from gc(reset = TRUE); not process RSS",
    cache_note = "metadata_cache_cold clears only the bwTools metadata cache, not the operating-system file cache"
  )

  out <- list(
    system = bw_benchmark_system_info(),
    file = file_info,
    config = config,
    regions = benchmark_regions[],
    results = results[],
    summary = summary[]
  )
  class(out) <- c("bwToolsBenchmark", "list")
  out
}

#' Print a bwTools benchmark
#'
#' @param x A `bwToolsBenchmark` object.
#' @param ... Additional arguments reserved for future use.
#' @return `x`, invisibly.
#' @export
print.bwToolsBenchmark <- function(x, ...) {
  cat("<bwTools BigWig benchmark>\n")
  if (!is.null(x$file) && nrow(x$file) > 0L) {
    cat("  file:", x$file[["file"]][1L], "\n")
    cat("  size:", format(round(x$file[["size_mb"]][1L], 2), nsmall = 2), "MB\n")
  }
  if (!is.null(x$regions)) {
    cat("  regions:", nrow(x$regions), "\n")
  }
  if (!is.null(x$results)) {
    cat("  benchmark rows:", nrow(x$results), "\n")
    cat("  errors:", sum(x$results[["status"]] == "error"), "\n")
    cat("  skipped:", sum(x$results[["status"]] == "skipped"), "\n")
  }
  if (!is.null(x$summary) && nrow(x$summary) > 0L) {
    cat("\nSummary:\n")
    n_show <- min(nrow(x$summary), 20L)
    print(x$summary[seq_len(n_show)])
    if (nrow(x$summary) > n_show) {
      cat("...", nrow(x$summary) - n_show, "additional summary rows\n")
    }
  }
  invisible(x)
}
