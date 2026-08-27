test_that("subset_bwg clips and merges overlapping or adjacent regions", {
  bedgraph_file <- system.file(
    "extdata",
    "bwtools_example.bedGraph",
    package = "bwTools",
    mustWork = TRUE
  )
  source <- read_bwg(
    bedgraph_file,
    format = "bedgraph",
    sample_names = "example",
    mode = "memory"
  )
  regions <- data.frame(
    chrom = c("chr1", "chr1"),
    start = c(510L, 561L),
    end = c(560L, 600L)
  )

  observed <- subset_bwg(source, regions)
  expected <- retrieve_bwg(source, "chr1", 510L, 600L)

  expect_true(is_bwg_track(observed))
  expect_identical(observed$meta$mode, "memory")
  expect_true(isTRUE(observed$meta$subset))
  expect_equal(nrow(observed$meta$subset_regions), 1L)
  expect_equal(observed$meta$subset_regions$start, 510L)
  expect_equal(observed$meta$subset_regions$end, 600L)
  expect_equal(observed$data, expected, tolerance = 1e-10)
  expect_equal(observed$data$start[1L], 510L)
  expect_equal(observed$data$end[nrow(observed$data)], 600L)
})

test_that("subset_bwg uses indexed lazy BigWig retrieval", {
  bw_file <- system.file(
    "extdata",
    "bwtools_example.bigwig",
    package = "bwTools",
    mustWork = TRUE
  )
  source <- read_bwg(
    bw_file,
    format = "bigwig",
    sample_names = "example",
    mode = "lazy"
  )
  regions <- data.frame(
    chrom = c("chr1", "chr2"),
    start = c(510L, 450L),
    end = c(600L, 525L)
  )

  observed <- subset_bwg(source, regions)

  expect_true(is_bwg_track(observed))
  expect_false(is.null(observed$data))
  expect_true(all(observed$data$chrom %in% c("chr1", "chr2")))
  expect_equal(unique(observed$seqinfo$chrom), c("chr1", "chr2"))
  expect_equal(unique(observed$seqinfo$length), c(10000L, 8000L))
  expect_true(all(observed$data$start >= c(510L, 450L)[match(observed$data$chrom, c("chr1", "chr2"))]))
  expect_true(all(observed$data$end <= c(600L, 525L)[match(observed$data$chrom, c("chr1", "chr2"))]))
})

test_that("subset_bwg writes BigWig subsets that round-trip exactly", {
  bw_file <- system.file(
    "extdata",
    "bwtools_example.bigwig",
    package = "bwTools",
    mustWork = TRUE
  )
  source <- read_bwg(
    bw_file,
    format = "bigwig",
    sample_names = "example",
    mode = "lazy"
  )
  regions <- data.frame(
    chrom = c("chr1", "chr2"),
    start = c(510L, 450L),
    end = c(1500L, 900L)
  )
  expected <- subset_bwg(source, regions)

  outdir <- tempfile(pattern = "bwtools_subset_bigwig_")
  dir.create(outdir)
  manifest <- subset_bwg(
    source,
    regions,
    outdir = outdir,
    format = "bigwig",
    overwrite = TRUE
  )

  expect_equal(manifest$status, "written")
  expect_true(file.exists(manifest$file))
  expect_equal(manifest$intervals, nrow(expected$data))

  observed <- read_bwg(
    manifest$file,
    format = "bigwig",
    sample_names = "example",
    mode = "memory"
  )
  expect_equal(observed$data, expected$data, tolerance = 1e-4)
  expect_equal(
    seqinfo_bwg(observed)$length,
    c(10000L, 8000L)
  )
})

test_that("subset_bwg supports WIG and bedGraph subset output", {
  bedgraph_file <- system.file(
    "extdata",
    "bwtools_example.bedGraph",
    package = "bwTools",
    mustWork = TRUE
  )
  source <- read_bwg(
    bedgraph_file,
    format = "bedgraph",
    sample_names = "example",
    mode = "memory"
  )
  regions <- data.frame(
    chrom = c("chr1", "chr3"),
    start = c(510L, 700L),
    end = c(1500L, 1200L)
  )
  expected <- subset_bwg(source, regions)$data

  for (format in c("wig", "bedgraph")) {
    outdir <- tempfile(pattern = paste0("bwtools_subset_", format, "_"))
    dir.create(outdir)
    manifest <- subset_bwg(
      source,
      regions,
      outdir = outdir,
      format = format,
      overwrite = TRUE
    )
    observed <- read_bwg(
      manifest$file,
      format = format,
      sample_names = "example",
      mode = "memory"
    )$data
    expect_equal(observed, expected, tolerance = 1e-4)
  }
})

test_that("subset_bwg keeps empty in-memory results and makes file behavior explicit", {
  bedgraph_file <- system.file(
    "extdata",
    "bwtools_example.bedGraph",
    package = "bwTools",
    mustWork = TRUE
  )
  source <- read_bwg(
    bedgraph_file,
    format = "bedgraph",
    sample_names = "example",
    mode = "memory"
  )
  regions <- data.frame(chrom = "chr1", start = 2000L, end = 2100L)

  empty_track <- subset_bwg(source, regions)
  expect_true(is_bwg_track(empty_track))
  expect_equal(nrow(empty_track$data), 0L)

  outdir <- tempfile(pattern = "bwtools_subset_empty_")
  dir.create(outdir)
  expect_error(
    subset_bwg(
      source,
      regions,
      outdir = outdir,
      format = "bedgraph"
    ),
    "No signal intervals overlap"
  )

  manifest <- subset_bwg(
    source,
    regions,
    outdir = outdir,
    format = "bedgraph",
    empty = "skip"
  )
  expect_equal(manifest$status, "empty")
  expect_equal(manifest$intervals, 0L)
  expect_true(is.na(manifest$file))
})

test_that("subset_bwg preserves selected samples without silent loss", {
  samples <- data.table::data.table(
    sample_id = c("sampleA", "sampleB"),
    strand = "*"
  )
  signal <- data.table::data.table(
    sample_id = c("sampleA", "sampleB"),
    chrom = c("chr1", "chr1"),
    start = c(10L, 100L),
    end = c(20L, 110L),
    value = c(1, 2),
    strand = "*"
  )
  source <- bw_track(samples, signal, meta = list(mode = "memory"))
  regions <- data.frame(chrom = "chr1", start = 1L, end = 50L)

  observed <- subset_bwg(source, regions)
  expect_equal(observed$samples$sample_id, c("sampleA", "sampleB"))
  expect_equal(unique(observed$data$sample_id), "sampleA")

  selected <- subset_bwg(source, regions, samples = "sampleA")
  expect_equal(selected$samples$sample_id, "sampleA")
  expect_equal(unique(selected$data$sample_id), "sampleA")

  outdir <- tempfile(pattern = "bwtools_subset_partial_empty_")
  dir.create(outdir)
  manifest <- subset_bwg(
    source,
    regions,
    outdir = outdir,
    format = "bedgraph",
    empty = "skip"
  )
  expect_equal(manifest$sample_id, c("sampleA", "sampleB"))
  expect_equal(manifest$status, c("written", "empty"))
  expect_true(file.exists(manifest$file[1L]))
  expect_true(is.na(manifest$file[2L]))
})

test_that("subset_bwg validates region coordinates", {
  track <- bw_track(
    samples = data.table::data.table(sample_id = "sampleA", strand = "*"),
    data = data.table::data.table(
      sample_id = "sampleA",
      chrom = "chr1",
      start = 10L,
      end = 20L,
      value = 1,
      strand = "*"
    ),
    meta = list(mode = "memory")
  )

  expect_error(
    subset_bwg(track, data.frame(chrom = "chr1", start = 0L, end = 10L)),
    "valid 1-based closed intervals"
  )
  expect_error(
    subset_bwg(track, data.frame(chrom = "chr1", start = 20L, end = 10L)),
    "valid 1-based closed intervals"
  )
  expect_error(
    subset_bwg(track, data.frame(chrom = "chr1", start = 1L)),
    "missing required columns"
  )
})
