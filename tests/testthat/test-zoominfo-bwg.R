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
    sample_ids = "example",
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

  track <- read_bwg(
    written$file,
    format = "bigwig",
    sample_ids = "example",
    mode = "lazy"
  )
  object_info <- zoominfo_bwg(track)
  expect_gt(nrow(object_info), 0L)
  expect_named(object_info, c("sample_id", "level", "data_offset", "index_offset"))
  expect_true(all(object_info$sample_id == "example"))
  expect_true(all(object_info$level > 0))
  expect_true(all(object_info$data_offset > 0))
  expect_true(all(object_info$index_offset > object_info$data_offset))
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
  track <- bwg_track(
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

  written_track <- read_bwg(
    written$file,
    format = "bigwig",
    sample_ids = "sampleA",
    mode = "lazy"
  )
  info <- zoominfo_bwg(written_track)
  expect_equal(nrow(info), 0L)
  expect_named(info, c("sample_id", "level", "data_offset", "index_offset"))
})
