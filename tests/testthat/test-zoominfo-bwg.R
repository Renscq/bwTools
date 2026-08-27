test_that("zoominfo_bwg reports native writer zoom levels", {
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
  outdir <- tempfile(pattern = "bwtools_zoominfo_")
  dir.create(outdir)
  written <- write_bwg(
    source,
    outdir = outdir,
    format = "bigwig",
    chrom_sizes = chrom_sizes_file,
    overwrite = TRUE,
    zoom = TRUE
  )

  path_info <- zoominfo_bwg(written$file)
  expect_gt(nrow(path_info), 0L)
  expect_named(path_info, c("level", "data_offset", "index_offset"))
  expect_true(all(path_info$level > 0))
  expect_true(all(path_info$data_offset > 0))
  expect_true(all(path_info$index_offset > path_info$data_offset))

  track <- read_bwg(
    written$file,
    format = "bigwig",
    sample_names = "example",
    mode = "lazy"
  )
  object_info <- zoominfo_bwg(track)
  expect_named(object_info, c("sample_id", "level", "data_offset", "index_offset"))
  expect_true(all(object_info$sample_id == "example"))
  expect_equal(object_info$level, path_info$level)
})

test_that("zoominfo_bwg returns an empty table when zoom is disabled", {
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
  outdir <- tempfile(pattern = "bwtools_zoominfo_nozoom_")
  dir.create(outdir)
  written <- write_bwg(
    track,
    outdir = outdir,
    format = "bigwig",
    chrom_sizes = data.frame(chrom = "chr1", length = 1000L),
    overwrite = TRUE,
    zoom = FALSE
  )

  info <- zoominfo_bwg(written$file)
  expect_equal(nrow(info), 0L)
  expect_named(info, c("level", "data_offset", "index_offset"))
})
