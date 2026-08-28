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
    region_sizes = c(100L, 2000L),
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
    region_sizes = 1000L,
    iterations = 1L,
    warmup = 0L,
    operations = c("merge", "write")
  )

  expect_true("merge" %in% benchmark$results$operation)
  expect_true("write" %in% benchmark$results$operation)
  write_rows <- benchmark$results[benchmark$results$operation == "write"]
  expect_setequal(write_rows$variant, c("zoom_on", "zoom_off"))
  expect_true(all(write_rows$status == "ok"))
  expect_true(all(write_rows$output_file_mb > 0))
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
