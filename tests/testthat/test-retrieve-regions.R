test_that("retrieve_bwg preserves the single-region data API", {
  samples <- data.frame(
    sample_id = "sampleA",
    file = NA_character_,
    format = "bedgraph",
    strand = "*"
  )
  signal <- data.frame(
    sample_id = "sampleA",
    chrom = "chr1",
    start = c(100L, 150L, 300L),
    end = c(120L, 180L, 320L),
    value = c(1, 2, 3),
    strand = "*"
  )
  track <- bwg_track(samples, signal, meta = list(mode = "memory"))

  observed <- retrieve_bwg(track, "chr1", 110L, 160L)

  expect_s3_class(observed, "data.table")
  expect_equal(observed$start, c(110L, 150L))
  expect_equal(observed$end, c(120L, 160L))
  expect_equal(observed$value, c(1, 2))
})

test_that("retrieve_bwg handles multiple regions as a genomic union", {
  samples <- data.frame(
    sample_id = "sampleA",
    file = NA_character_,
    format = "bedgraph",
    strand = "*"
  )
  signal <- data.frame(
    sample_id = "sampleA",
    chrom = "chr1",
    start = c(500L, 550L, 700L),
    end = c(524L, 574L, 724L),
    value = c(1, 2, 3),
    strand = "*"
  )
  track <- bwg_track(samples, signal, meta = list(mode = "memory"))
  regions <- data.frame(
    chrom = c("chr1", "chr1"),
    start = c(510L, 561L),
    end = c(560L, 600L)
  )

  observed <- retrieve_bwg(track, regions = regions)

  expect_equal(nrow(observed), 2L)
  expect_equal(observed$start, c(510L, 550L))
  expect_equal(observed$end, c(524L, 574L))
  expect_equal(observed$value, c(1, 2))
})

test_that("retrieve_bwg can return a write-ready BwgTrack", {
  samples <- data.frame(
    sample_id = c("sampleA", "sampleB"),
    file = c(NA_character_, NA_character_),
    format = c("bedgraph", "bedgraph"),
    strand = c("*", "*")
  )
  signal <- data.frame(
    sample_id = c("sampleA", "sampleB"),
    chrom = c("chr1", "chr2"),
    start = c(100L, 200L),
    end = c(120L, 220L),
    value = c(1, 2),
    strand = c("*", "*")
  )
  seqinfo <- data.frame(
    sample_id = c("sampleA", "sampleA", "sampleB", "sampleB"),
    chrom = c("chr1", "chr2", "chr1", "chr2"),
    length = c(1000L, 800L, 1000L, 800L)
  )
  track <- bwg_track(
    samples,
    signal,
    seqinfo = seqinfo,
    meta = list(mode = "memory")
  )
  regions <- data.frame(chrom = "chr1", start = 110L, end = 130L)

  observed <- retrieve_bwg(
    track,
    regions = regions,
    sample_ids = "sampleA",
    result = "track"
  )

  expect_true(is_bwg_track(observed))
  expect_identical(observed$meta$mode, "memory")
  expect_true(isTRUE(observed$meta$retrieved))
  expect_equal(nrow(observed$meta$retrieval_regions), 1L)
  expect_equal(observed$samples$sample_id, "sampleA")
  expect_equal(observed$data$start, 110L)
  expect_equal(observed$data$end, 120L)
  expect_equal(unique(observed$seqinfo$chrom), "chr1")
})

test_that("retrieve_bwg performs indexed multi-region lazy BigWig retrieval", {
  file <- system.file(
    "extdata",
    "bwtools_example.bigwig",
    package = "bwTools",
    mustWork = TRUE
  )
  track <- read_bwg(
    file,
    format = "bigwig",
    sample_ids = "example",
    mode = "lazy"
  )
  regions <- data.frame(
    chrom = c("chr1", "chr2"),
    start = c(510L, 450L),
    end = c(600L, 900L)
  )

  observed <- retrieve_bwg(track, regions = regions)

  expect_gt(nrow(observed), 0L)
  expect_true(all(observed$chrom %in% c("chr1", "chr2")))
  expect_true(all(
    (observed$chrom == "chr1" & observed$start >= 510L & observed$end <= 600L) |
      (observed$chrom == "chr2" & observed$start >= 450L & observed$end <= 900L)
  ))
})

test_that("retrieve_bwg validates query mode and sample IDs", {
  samples <- data.frame(
    sample_id = "sampleA",
    file = NA_character_,
    format = "bedgraph",
    strand = "*"
  )
  signal <- data.frame(
    sample_id = "sampleA",
    chrom = "chr1",
    start = 10L,
    end = 20L,
    value = 1,
    strand = "*"
  )
  track <- bwg_track(samples, signal, meta = list(mode = "memory"))

  expect_error(
    retrieve_bwg(
      track,
      chrom = "chr1",
      start = 1L,
      end = 10L,
      regions = data.frame(chrom = "chr1", start = 1L, end = 10L)
    ),
    "either `regions` or `chrom`/`start`/`end`"
  )
  expect_error(
    retrieve_bwg(track, regions = data.frame(chrom = "chr1", start = 0L, end = 10L)),
    "valid 1-based closed intervals"
  )
  expect_error(
    retrieve_bwg(track, "chr1", 1L, 10L, sample_ids = "missing"),
    "Unknown sample IDs"
  )
})

test_that("retrieved BwgTrack persists through write_bwg without retrieval-side output", {
  file <- system.file(
    "extdata",
    "bwtools_example.bigwig",
    package = "bwTools",
    mustWork = TRUE
  )
  source <- read_bwg(
    file,
    format = "bigwig",
    sample_ids = "example",
    mode = "lazy"
  )
  regions <- data.frame(
    chrom = c("chr1", "chr2"),
    start = c(510L, 450L),
    end = c(600L, 900L)
  )
  retrieved <- retrieve_bwg(
    source,
    regions = regions,
    result = "track"
  )

  outdir <- file.path(tempdir(), "bwtools_retrieve_write")
  if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)
  written <- write_bwg(
    retrieved,
    outdir = outdir,
    format = "bigwig",
    overwrite = TRUE
  )
  roundtrip <- read_bwg(
    written$file,
    format = "bigwig",
    sample_ids = "example",
    mode = "memory"
  )

  expected <- as.data.frame(retrieved$data)
  observed <- as.data.frame(roundtrip$data)
  row.names(expected) <- NULL
  row.names(observed) <- NULL
  expect_equal(observed, expected, tolerance = 1e-4)
})
