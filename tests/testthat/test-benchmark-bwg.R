test_that("benchmark_bwg returns a structured benchmark object", {
  file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )

  benchmark <- benchmark_bwg(
    file,
    sample_id = "example",
    region_sizes = c(500L, 2000L),
    n_bins = 10L,
    iterations = 1L,
    warmup = 0L,
    operations = c("read_lazy", "retrieve", "stats_zoom", "stats_full"),
    full_stats_max_bases = 2000L
  )

  expect_s3_class(benchmark, "bwToolsBenchmark")
  expect_named(
    benchmark,
    c("system", "file", "config", "regions", "results", "summary")
  )
  expect_s3_class(benchmark$results, "data.table")
  expect_s3_class(benchmark$summary, "data.table")
  expect_equal(nrow(benchmark$regions), 2L)
  expect_true(all(c(
    "operation", "variant", "region_id", "elapsed_sec", "baseline_heap_mb",
    "live_heap_after_mb", "peak_heap_mb", "peak_heap_delta_mb",
    "live_heap_delta_mb", "heap_delta_mb", "result_mb", "result_rows",
    "status", "message"
  ) %in% names(benchmark$results)))
  expect_true(all(c(
    "read_lazy", "retrieve", "stats_zoom", "stats_full"
  ) %in% unique(benchmark$results$operation)))
  expect_true(all(benchmark$results$status %in% c("ok", "error", "skipped")))
})

test_that("benchmark_bwg records bounded full-stat skips", {
  file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )

  benchmark <- benchmark_bwg(
    file,
    sample_id = "example",
    region_sizes = c(1000L, 2000L),
    iterations = 1L,
    warmup = 0L,
    operations = "stats_full",
    full_stats_max_bases = 500L
  )

  full <- benchmark$results[benchmark$results$operation == "stats_full"]
  expect_true(any(full$status == "ok"))
  expect_true(any(full$status == "skipped"))
  expect_match(
    full$message[full$status == "skipped"][1L],
    "full_stats_max_bases"
  )
})

test_that("benchmark_bwg supports explicit benchmark regions", {
  file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )
  regions <- data.frame(
    chrom = c("chr1", "chr2"),
    start = c(500L, 1000L),
    end = c(900L, 1500L)
  )

  benchmark <- benchmark_bwg(
    file,
    regions = regions,
    iterations = 1L,
    warmup = 0L,
    operations = "retrieve"
  )

  expect_equal(benchmark$regions$chrom, c("chr1", "chr2"))
  expect_equal(benchmark$regions$bases, c(401L, 501L))
  expect_equal(
    benchmark$results$operation,
    rep("retrieve", 2L)
  )
})

test_that("optional merge and write benchmark bounded selected regions", {
  file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )

  benchmark <- benchmark_bwg(
    file,
    sample_id = "example",
    region_sizes = c(1000L, 2000L),
    iterations = 1L,
    warmup = 0L,
    operations = c("merge", "write")
  )

  expect_true("merge" %in% benchmark$results$operation)
  expect_true("write" %in% benchmark$results$operation)
  write_rows <- benchmark$results[benchmark$results$operation == "write"]
  expect_equal(nrow(write_rows), 4L)
  expect_setequal(write_rows$variant, c("zoom_on", "zoom_off"))
  expect_setequal(write_rows$region_id, c("region_1", "region_2"))
  expect_true(all(write_rows$status == "ok"))
  expect_true(all(write_rows$output_file_mb > 0))
  expect_true(all(write_rows$input_signal_mb > 0))
  expect_true(all(write_rows$input_intervals > 0L))
  expect_true(all(write_rows$input_payload_mb > 0))
  expect_true(all(write_rows$input_data_blocks > 0L))
  expect_true(all(is.na(write_rows$intervals_per_sec) | write_rows$intervals_per_sec >= 0))
  expect_true(all(is.na(write_rows$input_payload_mb_per_sec) | write_rows$input_payload_mb_per_sec >= 0))

  write_summary <- benchmark$summary[benchmark$summary$operation == "write"]
  expect_true(all(c(
    "median_input_signal_mb",
    "median_input_intervals",
    "median_input_payload_mb",
    "median_output_file_mb",
    "median_output_to_signal_ratio",
    "median_payload_to_output_ratio",
    "median_intervals_per_sec",
    "zoom_elapsed_overhead_sec",
    "zoom_elapsed_overhead_pct",
    "zoom_size_overhead_mb",
    "zoom_size_overhead_pct"
  ) %in% names(write_summary)))
  zoom_on <- write_summary[write_summary$variant == "zoom_on"]
  expect_true(all(is.finite(zoom_on$zoom_elapsed_overhead_sec)))
  expect_true(all(is.finite(zoom_on$zoom_size_overhead_mb)))
})

test_that("read_memory profiling is explicit and works on bundled data", {
  file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )

  benchmark <- benchmark_bwg(
    file,
    iterations = 1L,
    warmup = 0L,
    operations = "read_memory"
  )

  expect_equal(benchmark$results$operation, "read_memory")
  expect_equal(benchmark$results$status, "ok")
  expect_true(benchmark$results$result_rows > 0L)
  expect_true(
    is.na(benchmark$results$file_mb_per_sec) ||
      benchmark$results$file_mb_per_sec >= 0
  )
})

test_that("benchmark_bwg rejects non-BigWig input", {
  file <- system.file(
    "extdata", "bwtools_example.bedGraph",
    package = "bwTools", mustWork = TRUE
  )
  expect_error(
    benchmark_bwg(file, iterations = 1L, warmup = 0L),
    "supports BigWig input only"
  )
})

test_that("explicit benchmark regions remain separate cases", {
  file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )
  regions <- data.frame(
    region_id = c("left", "right"),
    chrom = c("chr1", "chr1"),
    start = c(500L, 600L),
    end = c(599L, 699L)
  )

  benchmark <- benchmark_bwg(
    file,
    regions = regions,
    iterations = 1L,
    warmup = 0L,
    operations = "retrieve"
  )

  expect_equal(nrow(benchmark$regions), 2L)
  expect_equal(benchmark$regions$region_id, c("left", "right"))
  expect_equal(benchmark$regions$start, c(500L, 600L))
  expect_equal(benchmark$regions$end, c(599L, 699L))
})


test_that("benchmark memory metrics separate live and peak R heap", {
  file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )

  benchmark <- benchmark_bwg(
    file,
    sample_id = "example",
    region_sizes = 500L,
    iterations = 1L,
    warmup = 0L,
    operations = "retrieve"
  )

  result <- benchmark$results[1L]
  expect_true(is.finite(result$baseline_heap_mb))
  expect_true(is.finite(result$live_heap_after_mb))
  expect_true(is.finite(result$peak_heap_mb))
  expect_true(is.finite(result$peak_heap_delta_mb))
  expect_true(is.finite(result$live_heap_delta_mb))
  expect_equal(result$heap_delta_mb, result$peak_heap_delta_mb)
  expect_true(result$peak_heap_mb >= result$live_heap_after_mb)

  expect_true(all(c(
    "median_live_heap_after_mb",
    "median_peak_heap_mb",
    "median_peak_heap_delta_mb",
    "median_live_heap_delta_mb",
    "median_heap_delta_mb"
  ) %in% names(benchmark$summary)))
})


test_that("write benchmark respects write_max_bases", {
  file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )

  benchmark <- benchmark_bwg(
    file,
    sample_id = "example",
    region_sizes = c(100L, 2000L),
    iterations = 1L,
    warmup = 0L,
    operations = "write",
    write_max_bases = 1500L
  )

  write_rows <- benchmark$results[benchmark$results$operation == "write"]
  small <- write_rows[write_rows$region_id == "region_1"]
  large <- write_rows[write_rows$region_id == "region_2"]
  expect_true(all(small$status == "ok"))
  expect_true(all(large$status == "skipped"))
  expect_true(all(grepl("write_max_bases", large$message, fixed = TRUE)))
})


test_that("writer benchmark uses independent conservative repetitions", {
  file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )

  benchmark <- benchmark_bwg(
    file,
    sample_id = "example",
    region_sizes = 1000L,
    iterations = 3L,
    warmup = 1L,
    operations = "write",
    write_iterations = 1L,
    write_warmup = 0L,
    verbose = FALSE
  )

  write_rows <- benchmark$results[benchmark$results$operation == "write"]
  expect_equal(nrow(write_rows), 2L)
  expect_setequal(write_rows$variant, c("zoom_off", "zoom_on"))
  expect_true(all(write_rows$iteration == 1L))
})

test_that("zoom writer benchmark can be skipped independently by interval density", {
  file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )

  benchmark <- benchmark_bwg(
    file,
    sample_id = "example",
    region_sizes = 2000L,
    operations = "write",
    write_max_intervals = 100000L,
    write_zoom_max_intervals = 1L,
    write_iterations = 1L,
    write_warmup = 0L,
    verbose = FALSE
  )

  write_rows <- benchmark$results[benchmark$results$operation == "write"]
  zoom_off <- write_rows[write_rows$variant == "zoom_off"]
  zoom_on <- write_rows[write_rows$variant == "zoom_on"]
  expect_equal(zoom_off$status, "ok")
  expect_equal(zoom_on$status, "skipped")
  expect_match(zoom_on$message, "write_zoom_max_intervals", fixed = TRUE)
})

test_that("benchmark warmup helper runs the requested count once", {
  calls <- 0L
  fun <- function() {
    calls <<- calls + 1L
    invisible(NULL)
  }
  bwTools:::bw_benchmark_warmup(fun, 2L)
  expect_equal(calls, 2L)
})
