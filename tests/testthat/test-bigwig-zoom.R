test_that("native BigWig writer creates readable zoom levels", {
  bedgraph_file <- system.file(
    "extdata", "bwtools_example.bedGraph",
    package = "bwTools", mustWork = TRUE
  )
  chrom_sizes_file <- system.file(
    "extdata", "bwtools_example.chrom.sizes",
    package = "bwTools", mustWork = TRUE
  )
  source <- read_bwg(
    bedgraph_file,
    format = "bedgraph",
    sample_ids = "example",
    mode = "memory"
  )
  outdir <- tempfile(pattern = "bwtools_zoom_")
  dir.create(outdir)
  written <- write_bwg(
    source,
    outdir = outdir,
    format = "bigwig",
    chrom_sizes = chrom_sizes_file,
    overwrite = TRUE,
    zoom = TRUE
  )

  written_track <- read_bwg(
    written$file,
    format = "bigwig",
    sample_ids = "example",
    mode = "lazy"
  )
  zoom_info <- zoominfo_bwg(written_track)
  expect_gt(nrow(zoom_info), 0L)
  expect_true(all(zoom_info$level > 0L))
  expect_true(all(zoom_info$data_offset > 0))
  expect_true(all(zoom_info$index_offset > zoom_info$data_offset))

  zoom <- zoom_info[1L]
  records <- bwTools:::bw_bigwig_zoom_query(
    written$file,
    chrom = "chr1",
    start = 1L,
    end = 10000L,
    level = zoom$level
  )
  expect_gt(nrow(records), 0L)
  expect_true(all(records$valid_count > 0))
  expect_true(all(records$end0 > records$start0))
})

test_that("BigWig zoom levels can be disabled", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    start = c(10L, 100L),
    end = c(20L, 120L),
    value = c(1, 2),
    strand = "*"
  )
  track <- bwg_track(
    data.table::data.table(sample_id = "sampleA", strand = "*"),
    signal,
    meta = list(mode = "memory")
  )
  outdir <- tempfile(pattern = "bwtools_nozoom_")
  dir.create(outdir)
  written <- write_bwg(
    track,
    outdir = outdir,
    format = "bigwig",
    chrom_sizes = data.frame(chrom = "chr1", length = 1000L),
    overwrite = TRUE,
    zoom = FALSE
  )
  written_track <- read_bwg(
    written$file,
    format = "bigwig",
    sample_ids = "sampleA",
    mode = "lazy"
  )
  zoom_info <- zoominfo_bwg(written_track)
  expect_equal(nrow(zoom_info), 0L)
})

test_that("zoom summary records preserve whole-window aggregates", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    start = c(1L, 11L, 31L),
    end = c(10L, 20L, 40L),
    value = c(1, 3, 2),
    strand = "*"
  )
  chrom_sizes <- bwTools:::bw_read_chrom_sizes(
    data.frame(chrom = "chr1", length = 100L)
  )
  signal2 <- bwTools:::bw_validate_signal_for_write(signal, chrom_sizes)
  records <- bwTools:::bw_build_zoom_records(signal2, chrom_sizes, level = 20L)

  expect_equal(records$start0, c(0, 20))
  expect_equal(records$end0, c(20, 40))
  expect_equal(records$valid_count, c(20, 10))
  expect_equal(records$min_value, c(1, 2))
  expect_equal(records$max_value, c(3, 2))
  expect_equal(records$sum_data, c(40, 20))
  expect_equal(records$sum_squared, c(100, 40))
})

test_that("hierarchical zoom aggregation matches direct full-signal summaries", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = c(rep("chr1", 6L), rep("chr2", 3L)),
    start = c(1L, 12L, 19L, 41L, 88L, 131L, 5L, 37L, 81L),
    end = c(8L, 18L, 35L, 67L, 105L, 170L, 20L, 60L, 120L),
    value = c(1, 2, -1, 4, 3, 0.5, 7, -2, 5),
    strand = "*"
  )
  chrom_sizes <- bwTools:::bw_read_chrom_sizes(
    data.frame(chrom = c("chr1", "chr2"), length = c(200L, 150L))
  )
  signal2 <- bwTools:::bw_validate_signal_for_write(signal, chrom_sizes)

  level1 <- bwTools:::bw_build_zoom_records(signal2, chrom_sizes, level = 20L)
  hierarchical <- bwTools:::bw_aggregate_zoom_records(
    level1,
    chrom_sizes,
    level = 80L
  )
  direct <- bwTools:::bw_build_zoom_records(signal2, chrom_sizes, level = 80L)

  cols <- c(
    "tid", "start0", "end0", "valid_count",
    "min_value", "max_value", "sum_data", "sum_squared"
  )
  expect_equal(hierarchical[, ..cols], direct[, ..cols], tolerance = 1e-10)
})

test_that("hierarchical zoom pyramid preserves direct summaries across levels", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    start = seq.int(1L, 9991L, by = 10L),
    end = pmin(seq.int(1L, 9991L, by = 10L) + 6L, 10000L),
    value = sin(seq_len(1000L) / 11),
    strand = "*"
  )
  chrom_sizes <- bwTools:::bw_read_chrom_sizes(
    data.frame(chrom = "chr1", length = 10000L)
  )
  signal2 <- bwTools:::bw_validate_signal_for_write(signal, chrom_sizes)
  levels <- c(20, 80, 320, 1280)

  pyramid <- bwTools:::bw_build_zoom_pyramid(signal2, chrom_sizes, levels)
  expect_length(pyramid, length(levels))

  cols <- c(
    "tid", "start0", "end0", "valid_count",
    "min_value", "max_value", "sum_data", "sum_squared"
  )
  for (i in seq_along(levels)) {
    direct <- bwTools:::bw_build_zoom_records(signal2, chrom_sizes, levels[i])
    expect_equal(pyramid[[i]][, ..cols], direct[, ..cols], tolerance = 1e-10)
  }
})

test_that("first zoom level splits boundary-crossing intervals exactly", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    start = c(18L, 45L),
    end = c(26L, 53L),
    value = c(2, -3),
    strand = "*"
  )
  chrom_sizes <- bwTools:::bw_read_chrom_sizes(
    data.frame(chrom = "chr1", length = 53L)
  )
  signal2 <- bwTools:::bw_validate_signal_for_write(signal, chrom_sizes)
  records <- bwTools:::bw_build_zoom_records(signal2, chrom_sizes, level = 20L)

  expect_equal(records$start0, c(0, 20, 40))
  expect_equal(records$end0, c(20, 40, 53))
  expect_equal(records$valid_count, c(3, 6, 9))
  expect_equal(records$min_value, c(2, 2, -3))
  expect_equal(records$max_value, c(2, 2, -3))
  expect_equal(records$sum_data, c(6, 12, -27))
  expect_equal(records$sum_squared, c(12, 24, 81))
})

test_that("written hierarchical zoom levels match direct summaries", {
  bedgraph_file <- system.file(
    "extdata", "bwtools_example.bedGraph",
    package = "bwTools", mustWork = TRUE
  )
  chrom_sizes_file <- system.file(
    "extdata", "bwtools_example.chrom.sizes",
    package = "bwTools", mustWork = TRUE
  )
  source <- read_bwg(
    bedgraph_file,
    format = "bedgraph",
    sample_ids = "example",
    mode = "memory"
  )
  chrom_sizes <- bwTools:::bw_read_chrom_sizes(chrom_sizes_file)
  signal <- bwTools:::bw_validate_signal_for_write(source$data, chrom_sizes)

  outdir <- tempfile(pattern = "bwtools_hierarchical_zoom_")
  dir.create(outdir)
  written <- write_bwg(
    source,
    outdir = outdir,
    format = "bigwig",
    chrom_sizes = chrom_sizes_file,
    overwrite = TRUE,
    zoom = TRUE
  )
  written_track <- read_bwg(
    written$file,
    format = "bigwig",
    sample_ids = "example",
    mode = "lazy"
  )
  zoom_info <- zoominfo_bwg(written_track)
  expect_gt(nrow(zoom_info), 1L)

  for (level in zoom_info$level) {
    observed <- bwTools:::bw_bigwig_zoom_query(
      written$file,
      chrom = "chr1",
      start = 1L,
      end = 10000L,
      level = level
    )
    expected <- bwTools:::bw_build_zoom_records(
      signal,
      chrom_sizes,
      level = level
    )[tid == 0]

    expect_equal(observed$start0, expected$start0)
    expect_equal(observed$end0, expected$end0)
    expect_equal(observed$valid_count, expected$valid_count)
    expect_equal(observed$min_value, expected$min_value, tolerance = 1e-5)
    expect_equal(observed$max_value, expected$max_value, tolerance = 1e-5)
    expect_equal(observed$sum_data, expected$sum_data, tolerance = 1e-4)
    expect_equal(observed$sum_squared, expected$sum_squared, tolerance = 1e-3)
  }
})
