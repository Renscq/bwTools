test_that("bwg_track normalizes the public object schema", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    start = 10,
    end = 20,
    value = 2
  )
  x <- bwg_track(
    samples = data.frame(sample_id = "sampleA"),
    data = signal,
    seqinfo = data.frame(chrom = "chr1", length = 1000L)
  )

  expect_true(is_bwg_track(x))
  expect_invisible(validate_bwg(x))
  expect_named(
    samples_bwg(x),
    c("sample_id", "file", "format", "strand", "has_strand")
  )
  expect_equal(samples_bwg(x)$format, "memory")
  expect_equal(samples_bwg(x)$strand, "*")
  expect_named(
    x$data,
    c("sample_id", "chrom", "start", "end", "value", "strand")
  )
  expect_equal(seqinfo_bwg(x)$sample_id, "sampleA")
  expect_identical(metadata_bwg(x)$coordinate, "1-based closed")
  expect_identical(x$meta$coordinate, "1-based closed")
  expect_identical(x$meta$schema_version, "2")
  expect_identical(x$meta$mode, "memory")
})

test_that("validate_bwg rejects shallow or inconsistent objects", {
  shallow <- structure(
    list(samples = data.frame(sample_id = "sampleA"), meta = list()),
    class = "BwgTrack"
  )
  expect_false(is_bwg_track(shallow))
  expect_error(validate_bwg(shallow), "Invalid BwgTrack contract")

  x <- bwg_track(
    samples = data.frame(sample_id = "sampleA"),
    data = data.frame(
      sample_id = "sampleA",
      chrom = "chr1",
      start = 1L,
      end = 10L,
      value = 1,
      strand = "*"
    )
  )
  x$meta$coordinate <- "0-based half-open"
  expect_false(is_bwg_track(x))
  expect_error(validate_bwg(x), "meta\\$coordinate")


  y <- bwg_track(
    samples = data.frame(sample_id = "sampleA"),
    data = data.frame(
      sample_id = "sampleA",
      chrom = "chr1",
      start = 1L,
      end = 10L,
      value = 1,
      strand = "*"
    )
  )
  y$data$start <- as.numeric(y$data$start)
  expect_false(is_bwg_track(y))
  expect_error(validate_bwg(y), "data\\$start")
})

test_that("sample accessors use sample_ids consistently", {
  signal <- data.table::data.table(
    sample_id = c("sampleA", "sampleB"),
    chrom = "chr1",
    start = c(1L, 11L),
    end = c(10L, 20L),
    value = c(1, 2),
    strand = "*"
  )
  x <- bwg_track(
    samples = data.frame(sample_id = c("sampleA", "sampleB")),
    data = signal,
    seqinfo = data.frame(chrom = "chr1", length = 100L)
  )

  expect_equal(samples_bwg(x, sample_ids = "sampleB")$sample_id, "sampleB")
  expect_equal(seqinfo_bwg(x, sample_ids = "sampleB")$sample_id, "sampleB")
  expect_error(samples_bwg(x, sample_ids = "missing"), "Unknown sample IDs")
  expect_error(seqinfo_bwg(x, sample_ids = "missing"), "Unknown sample IDs")
})

test_that("read_bwg auto mode is safe for BigWig and text input", {
  bw_file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )
  bedgraph_file <- system.file(
    "extdata", "bwtools_example.bedGraph",
    package = "bwTools", mustWork = TRUE
  )

  expect_identical(detect_bwg_format(bw_file), "bigwig")
  expect_identical(detect_bwg_format(bedgraph_file), "bedgraph")

  bw <- read_bwg(bw_file, sample_ids = "bw")
  text <- read_bwg(bedgraph_file, sample_ids = "bg")
  mixed <- read_bwg(
    c(bw_file, bedgraph_file),
    sample_ids = c("bw", "bg")
  )

  expect_identical(summary_bwg(bw)$mode, "lazy")
  expect_identical(summary_bwg(text)$mode, "memory")
  expect_identical(summary_bwg(mixed)$mode, "memory")
  expect_true(is_bwg_track(mixed))
})

test_that("main analysis APIs require BwgTrack input", {
  bw_file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )
  expect_error(seqinfo_bwg(bw_file), "Invalid BwgTrack contract")
  expect_error(zoominfo_bwg(bw_file), "Invalid BwgTrack contract")
  expect_error(
    stats_bwg(bw_file, "chr1", 1L, 100L),
    "Invalid BwgTrack contract"
  )
})


test_that("public data-table returns use stable core columns", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    start = 10L,
    end = 20L,
    value = 2,
    strand = "*"
  )
  x <- bwg_track(
    samples = data.frame(sample_id = "sampleA"),
    data = signal,
    seqinfo = data.frame(chrom = "chr1", length = 100L)
  )

  retrieved <- retrieve_bwg(x, "chr1", 1L, 50L)
  stats <- stats_bwg(
    x,
    chrom = "chr1",
    start = 1L,
    end = 50L,
    use_zoom = FALSE
  )
  summary <- summary_bwg(x)
  zoom <- zoominfo_bwg(x)

  expect_named(
    retrieved,
    c("sample_id", "chrom", "start", "end", "value", "strand")
  )
  expect_named(
    stats,
    c(
      "sample_id", "chrom", "start", "end", "stat", "value",
      "covered_bases", "source", "resolution"
    )
  )
  expect_named(
    summary,
    c(
      "samples", "mode", "intervals", "chromosomes", "coordinate",
      "schema_version", "backend"
    )
  )
  expect_named(zoom, c("sample_id", "level", "data_offset", "index_offset"))
})

test_that("merge validates a malformed single BwgTrack without flattening it", {
  malformed <- structure(
    list(
      samples = data.frame(sample_id = "sampleA"),
      data = NULL,
      seqinfo = NULL,
      meta = list()
    ),
    class = c("BwgTrack", "bwToolsTrack", "list")
  )
  valid <- bwg_track(
    samples = data.frame(sample_id = "sampleB"),
    data = data.frame(
      sample_id = "sampleB", chrom = "chr1", start = 1L, end = 2L,
      value = 1, strand = "*"
    )
  )

  expect_error(
    merge_bwg(malformed, valid),
    "Input track 1 does not satisfy the BwgTrack contract"
  )
})
