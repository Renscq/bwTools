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
    sample_names = "example",
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

  zoom_info <- zoominfo_bwg(written$file)
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
  track <- bw_track(
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
  zoom_info <- zoominfo_bwg(written$file)
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
