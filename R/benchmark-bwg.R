# Author: Rensc
# Date: 2026-08-29
# Version: dev005
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
#' @param write_max_bases Maximum region size for the optional `write`
#'   benchmark. Larger regions are recorded as skipped.
#' @param write_max_intervals Maximum number of signal intervals allowed for
#'   any writer benchmark case. Denser regions are skipped before writing.
#' @param write_zoom_max_intervals Maximum number of signal intervals allowed
#'   for `zoom_on` writer profiling. `zoom_off` may still run up to
#'   `write_max_intervals`. This protects the pure-R zoom builder from
#'   unexpectedly expensive dense-signal benchmarks.
#' @param write_iterations Timed iterations used by writer profiling. Writer
#'   benchmarks use an independent conservative default because each iteration
#'   creates a complete BigWig file.
#' @param write_warmup Untimed writer warmup iterations per writer case.
#' @param verbose Whether to print concise benchmark progress messages.
#' @return A `bwToolsBenchmark` object containing `system`, `file`, `config`,
#'   `regions`, `results`, and `summary` components. Memory metrics distinguish
#'   current live R heap after garbage collection from the approximate maximum
#'   R heap observed since `gc(reset = TRUE)`; neither metric is process RSS.
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
  full_stats_max_bases = 1e6,
  write_max_bases = 1e7,
  write_max_intervals = 2e6,
  write_zoom_max_intervals = 5e5,
  write_iterations = 1L,
  write_warmup = 0L,
  verbose = TRUE
) {
  file <- bw_validate_local_file(file)
  if (!identical(detect_bwg_format(file), "bigwig")) {
    bw_stop("`benchmark_bwg()` currently supports BigWig input only.")
  }
  sample_id <- as.character(sample_id)
  if (length(sample_id) != 1L || is.na(sample_id) || !nzchar(sample_id)) {
    bw_stop("`sample_id` must be a single non-empty value.")
  }
  stat <- bw_match_arg(
    stat,
    c("mean", "stdev", "max", "min", "coverage", "sum"),
    "stat"
  )
  n_bins <- bw_benchmark_positive_integer(n_bins, "n_bins")
  iterations <- bw_benchmark_positive_integer(iterations, "iterations")
  warmup <- bw_benchmark_positive_integer(warmup, "warmup", allow_zero = TRUE)
  full_stats_max_bases <- bw_benchmark_positive_integer(
    full_stats_max_bases,
    "full_stats_max_bases"
  )
  write_max_bases <- bw_benchmark_positive_integer(
    write_max_bases,
    "write_max_bases"
  )
  write_max_intervals <- bw_benchmark_positive_integer(
    write_max_intervals,
    "write_max_intervals"
  )
  write_zoom_max_intervals <- bw_benchmark_positive_integer(
    write_zoom_max_intervals,
    "write_zoom_max_intervals"
  )
  write_iterations <- bw_benchmark_positive_integer(
    write_iterations,
    "write_iterations"
  )
  write_warmup <- bw_benchmark_positive_integer(
    write_warmup,
    "write_warmup",
    allow_zero = TRUE
  )
  verbose <- bw_validate_flag(verbose, "verbose")

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

    warm_fun <- function() read_bwg(file, sample_ids = sample_id, mode = "lazy")
    bw_benchmark_warmup(warm_fun, warmup)
    for (iteration in seq_len(iterations)) {
      measurement <- bw_benchmark_measure(warm_fun, warmup = 0L)
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
      benchmark_fun <- function() retrieve_bwg(
        lazy_track,
        chrom = region[["chrom"]],
        start = region[["start"]],
        end = region[["end"]],
        sample_ids = sample_id,
        result = "data"
      )
      bw_benchmark_warmup(benchmark_fun, warmup)
      for (iteration in seq_len(iterations)) {
        measurement <- bw_benchmark_measure(benchmark_fun, warmup = 0L)
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
      benchmark_fun <- function() stats_bwg(
        lazy_track,
        chrom = region[["chrom"]],
        start = region[["start"]],
        end = region[["end"]],
        n_bins = n_bins,
        stat = stat,
        sample_ids = sample_id,
        use_zoom = TRUE
      )
      bw_benchmark_warmup(benchmark_fun, warmup)
      for (iteration in seq_len(iterations)) {
        measurement <- bw_benchmark_measure(benchmark_fun, warmup = 0L)
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
      benchmark_fun <- function() stats_bwg(
        lazy_track,
        chrom = region[["chrom"]],
        start = region[["start"]],
        end = region[["end"]],
        n_bins = n_bins,
        stat = stat,
        sample_ids = sample_id,
        use_zoom = FALSE
      )
      bw_benchmark_warmup(benchmark_fun, warmup)
      for (iteration in seq_len(iterations)) {
        measurement <- bw_benchmark_measure(benchmark_fun, warmup = 0L)
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
      benchmark_fun <- function() merge_bwg(tracks = merge_tracks, method = "auto")
      bw_benchmark_warmup(benchmark_fun, warmup)
      for (iteration in seq_len(iterations)) {
        measurement <- bw_benchmark_measure(benchmark_fun, warmup = 0L)
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
    for (region_i in seq_len(nrow(benchmark_regions))) {
      region <- benchmark_regions[region_i]
      if (as.numeric(region[["bases"]]) > write_max_bases) {
        for (variant in c("zoom_off", "zoom_on")) {
          add_row(bw_benchmark_skipped_row(
            operation = "write",
            variant = variant,
            file_mb = file_mb,
            region = region,
            message = paste0(
              "Region exceeds `write_max_bases` (",
              format(write_max_bases, scientific = FALSE),
              " bp)."
            )
          ))
        }
        next
      }

      bw_benchmark_progress(
        verbose,
        "materialize write input ", region[["region_id"]],
        " (", format(region[["bases"]], big.mark = ",", scientific = FALSE), " bp)"
      )
      write_track <- retrieve_bwg(
        lazy_track,
        chrom = region[["chrom"]],
        start = region[["start"]],
        end = region[["end"]],
        sample_ids = sample_id,
        result = "track"
      )
      if (is.null(write_track$data) || nrow(write_track$data) < 1L) {
        for (variant in c("zoom_off", "zoom_on")) {
          add_row(bw_benchmark_skipped_row(
            operation = "write",
            variant = variant,
            file_mb = file_mb,
            region = region,
            message = "No signal intervals are available in this benchmark region."
          ))
        }
        rm(write_track)
        invisible(gc())
        next
      }
      input_metrics <- bw_benchmark_write_input_metrics(write_track)
      if (input_metrics$input_intervals > write_max_intervals) {
        for (variant in c("zoom_off", "zoom_on")) {
          add_row(bw_benchmark_skipped_row(
            operation = "write",
            variant = variant,
            file_mb = file_mb,
            region = region,
            message = paste0(
              "Region contains ", format(input_metrics$input_intervals, big.mark = ","),
              " intervals, exceeding `write_max_intervals` (",
              format(write_max_intervals, big.mark = ","), ")."
            )
          ))
        }
        bw_benchmark_progress(
          verbose,
          "skip ", region[["region_id"]], ": ",
          format(input_metrics$input_intervals, big.mark = ","),
          " intervals exceed write_max_intervals"
        )
        rm(write_track)
        invisible(gc())
        next
      }

      for (zoom_value in c(FALSE, TRUE)) {
        variant <- if (isTRUE(zoom_value)) "zoom_on" else "zoom_off"
        if (isTRUE(zoom_value) &&
            input_metrics$input_intervals > write_zoom_max_intervals) {
          add_row(bw_benchmark_skipped_row(
            operation = "write",
            variant = variant,
            file_mb = file_mb,
            region = region,
            message = paste0(
              "Region contains ", format(input_metrics$input_intervals, big.mark = ","),
              " intervals, exceeding `write_zoom_max_intervals` (",
              format(write_zoom_max_intervals, big.mark = ","), ")."
            )
          ))
          bw_benchmark_progress(
            verbose,
            "skip ", region[["region_id"]], " ", variant, ": ",
            format(input_metrics$input_intervals, big.mark = ","),
            " intervals exceed write_zoom_max_intervals"
          )
          next
        }

        bw_benchmark_progress(
          verbose,
          "write ", region[["region_id"]], " ", variant,
          " (", format(input_metrics$input_intervals, big.mark = ","), " intervals; ",
          write_warmup, " warmup + ", write_iterations, " timed)"
        )

        write_once <- function() {
          run_dir <- tempfile(paste0("bwtools_benchmark_write_", variant, "_"))
          dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
          on.exit(unlink(run_dir, recursive = TRUE, force = TRUE), add = TRUE)
          write_bwg(
            write_track,
            outdir = run_dir,
            format = "bigwig",
            sample_ids = sample_id,
            overwrite = TRUE,
            zoom = zoom_value
          )
        }
        bw_benchmark_warmup(write_once, write_warmup)

        for (iteration in seq_len(write_iterations)) {
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
            warmup = 0L
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
            region = region,
            output_mb = output_mb,
            input_signal_mb = input_metrics$input_signal_mb,
            input_intervals = input_metrics$input_intervals,
            input_payload_mb = input_metrics$input_payload_mb,
            input_covered_bases = input_metrics$input_covered_bases,
            input_data_blocks = input_metrics$input_data_blocks,
            note = paste(
              "Writer timing excludes indexed retrieval; the selected region",
              "is materialized before timing write_bwg()."
            )
          ))
          bw_benchmark_progress(
            verbose,
            "done ", region[["region_id"]], " ", variant,
            " iteration ", iteration, "/", write_iterations,
            " in ", format(round(measurement$elapsed_sec, 2), nsmall = 2), " s"
          )
          unlink(run_dir, recursive = TRUE, force = TRUE)
        }
      }
      rm(write_track)
      invisible(gc())
    }
  }

  if ("read_memory" %in% operations) {
    benchmark_fun <- function() read_bwg(file, sample_ids = sample_id, mode = "memory")
    bw_benchmark_warmup(benchmark_fun, warmup)
    for (iteration in seq_len(iterations)) {
      bw_benchmark_clear_metadata_cache(file)
      measurement <- bw_benchmark_measure(benchmark_fun, warmup = 0L)
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
    write_max_bases = write_max_bases,
    write_max_intervals = write_max_intervals,
    write_zoom_max_intervals = write_zoom_max_intervals,
    write_iterations = write_iterations,
    write_warmup = write_warmup,
    verbose = verbose,
    writer_metric = paste(
      "write timing excludes region retrieval; input_payload_mb estimates the",
      "uncompressed primary type-1 BigWig data payload and excludes indexes/zoom"
    ),
    memory_metric = paste(
      "R heap metrics from gc(): live_heap_after_mb is current used heap after GC;",
      "peak_heap_mb is max used since gc(reset = TRUE); neither is process RSS"
    ),
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
