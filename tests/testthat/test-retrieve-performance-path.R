test_that("indexed retrieval preserves bundled BigWig signal after chunked-vector optimization", {
  bw_file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )
  bedgraph_file <- system.file(
    "extdata", "bwtools_example.bedGraph",
    package = "bwTools", mustWork = TRUE
  )

  lazy <- read_bwg(
    bw_file,
    format = "bigwig",
    sample_ids = "example",
    mode = "lazy"
  )
  reference <- read_bwg(
    bedgraph_file,
    format = "bedgraph",
    sample_ids = "example",
    mode = "memory"
  )

  observed <- retrieve_bwg(
    lazy,
    chrom = "chr1",
    start = 450L,
    end = 3600L
  )
  expected <- retrieve_bwg(
    reference,
    chrom = "chr1",
    start = 450L,
    end = 3600L
  )

  expect_equal(observed$sample_id, expected$sample_id)
  expect_equal(observed$chrom, expected$chrom)
  expect_equal(observed$start, expected$start)
  expect_equal(observed$end, expected$end)
  expect_equal(observed$value, expected$value, tolerance = 1e-4)
  expect_equal(observed$strand, expected$strand)
})


test_that("chunked indexed retrieval concatenates multiple BigWig data blocks", {
  n <- 4000L
  starts <- seq.int(1L, by = 2L, length.out = n)
  signal <- data.table::data.table(
    sample_id = rep("sampleA", n),
    chrom = rep("chr1", n),
    start = starts,
    end = starts,
    value = as.numeric(seq_len(n)),
    strand = rep("*", n)
  )
  track <- bwg_track(
    samples = data.frame(sample_id = "sampleA"),
    data = signal,
    seqinfo = data.frame(chrom = "chr1", length = 10000L)
  )

  outdir <- tempfile("bwtools_chunked_query_")
  dir.create(outdir)
  written <- write_bwg(
    track,
    outdir = outdir,
    format = "bigwig",
    overwrite = TRUE,
    zoom = FALSE
  )
  lazy <- read_bwg(
    written$file,
    sample_ids = "sampleA",
    mode = "lazy"
  )
  observed <- retrieve_bwg(
    lazy,
    chrom = "chr1",
    start = 1L,
    end = 10000L
  )

  expect_equal(nrow(observed), n)
  expect_equal(observed$start, signal$start)
  expect_equal(observed$end, signal$end)
  expect_equal(observed$value, signal$value, tolerance = 1e-5)
})
