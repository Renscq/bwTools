test_that("native BigWig writer round-trips bundled example signal", {
  bedgraph_file <- system.file(
    "extdata",
    "bwtools_example.bedGraph",
    package = "bwTools",
    mustWork = TRUE
  )
  chrom_sizes_file <- system.file(
    "extdata",
    "bwtools_example.chrom.sizes",
    package = "bwTools",
    mustWork = TRUE
  )

  source <- read_bwg(
    bedgraph_file,
    format = "bedgraph",
    sample_names = "example",
    mode = "memory"
  )
  outdir <- tempfile(pattern = "bwtools_writer_")
  dir.create(outdir)
  written <- write_bwg(
    source,
    outdir = outdir,
    format = "bigwig",
    chrom_sizes = chrom_sizes_file,
    overwrite = TRUE
  )

  expect_equal(nrow(written), 1L)
  expect_true(file.exists(written$file))
  expect_identical(written$format, "bigwig")

  observed <- read_bwg(
    written$file,
    format = "bigwig",
    sample_names = "example",
    mode = "lazy"
  )
  expected_seqinfo <- data.table::data.table(
    sample_id = "example",
    chrom = c("chr1", "chr2", "chr3"),
    length = c(10000L, 8000L, 5000L)
  )
  expect_equal(seqinfo_bwg(observed), expected_seqinfo)

  for (chrom in c("chr1", "chr2", "chr3")) {
    chrom_length <- expected_seqinfo$length[match(chrom, expected_seqinfo$chrom)]
    expected <- retrieve_bwg(source, chrom, 1L, chrom_length)
    actual <- retrieve_bwg(observed, chrom, 1L, chrom_length)
    columns <- c("sample_id", "chrom", "start", "end", "value", "strand")
    expect_equal(
      as.data.frame(actual)[, columns, drop = FALSE],
      as.data.frame(expected)[, columns, drop = FALSE],
      tolerance = 1e-4
    )
  }
})

test_that("BigWig writer can use complete BwgTrack seqinfo without chrom_sizes", {
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
    mode = "memory"
  )
  outdir <- tempfile(pattern = "bwtools_writer_seqinfo_")
  dir.create(outdir)
  written <- write_bwg(
    source,
    outdir = outdir,
    format = "bigwig",
    overwrite = TRUE
  )
  expect_true(file.exists(written$file))
  expect_equal(seqinfo_bwg(written$file)$length, c(10000L, 8000L, 5000L))
})

test_that("BigWig writer rejects overlaps and chromosome overflow", {
  samples <- data.table::data.table(sample_id = "sampleA", strand = "*")
  overlap_signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = c("chr1", "chr1"),
    start = c(10L, 15L),
    end = c(20L, 25L),
    value = c(1, 2),
    strand = "*"
  )
  overlap_track <- bw_track(samples, overlap_signal, meta = list(mode = "memory"))
  expect_error(
    write_bwg(
      overlap_track,
      outdir = tempfile(pattern = "bwtools_overlap_"),
      format = "bigwig",
      chrom_sizes = data.frame(chrom = "chr1", length = 100L)
    ),
    "Overlapping signal intervals"
  )

  overflow_signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = "chr1",
    start = 95L,
    end = 110L,
    value = 1,
    strand = "*"
  )
  overflow_track <- bw_track(samples, overflow_signal, meta = list(mode = "memory"))
  expect_error(
    write_bwg(
      overflow_track,
      outdir = tempfile(pattern = "bwtools_overflow_"),
      format = "bigwig",
      chrom_sizes = data.frame(chrom = "chr1", length = 100L)
    ),
    "exceeds chromosome length"
  )
})

test_that("native writer preserves bedGraph and WIG signal intervals", {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = c("chr1", "chr1", "chr2"),
    start = c(10L, 30L, 5L),
    end = c(12L, 31L, 7L),
    value = c(1.5, -2.25, 3.75),
    strand = "*"
  )
  track <- bw_track(
    samples = data.table::data.table(sample_id = "sampleA", strand = "*"),
    data = signal,
    meta = list(mode = "memory")
  )
  expected <- signal[, c("sample_id", "chrom", "start", "end", "value", "strand"), with = FALSE]
  data.table::setorderv(expected, c("sample_id", "chrom", "start", "end"))

  for (format in c("bedgraph", "wig")) {
    outdir <- tempfile(pattern = paste0("bwtools_", format, "_"))
    dir.create(outdir)
    written <- write_bwg(track, outdir = outdir, format = format, overwrite = TRUE)
    observed <- read_bwg(
      written$file,
      format = format,
      sample_names = "sampleA",
      mode = "memory"
    )$data
    expect_equal(observed, expected, tolerance = 1e-10)
  }
})

test_that("lazy BigWig tracks can be copied without decoding", {
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
  outdir <- tempfile(pattern = "bwtools_lazy_copy_")
  dir.create(outdir)
  written <- write_bwg(
    source,
    outdir = outdir,
    format = "bigwig",
    overwrite = TRUE
  )
  expect_true(file.exists(written$file))
  expect_identical(unname(tools::md5sum(written$file)), unname(tools::md5sum(bw_file)))
})

test_that("native writer supports multi-level chromosome and R-tree indexes", {
  n_chrom <- 300L
  chrom <- sprintf("chr%03d", seq_len(n_chrom))
  samples <- data.table::data.table(sample_id = "sampleA", strand = "*")
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = chrom,
    start = 11L,
    end = 20L,
    value = as.numeric((seq_len(n_chrom) %% 7L) - 3L),
    strand = "*"
  )
  track <- bw_track(samples, signal, meta = list(mode = "memory"))
  chrom_sizes <- data.frame(chrom = chrom, length = rep(1000L, n_chrom))
  outdir <- tempfile(pattern = "bwtools_multilevel_")
  dir.create(outdir)

  written <- write_bwg(
    track,
    outdir = outdir,
    format = "bigwig",
    chrom_sizes = chrom_sizes,
    overwrite = TRUE
  )
  observed_seqinfo <- seqinfo_bwg(written$file)
  expect_equal(nrow(observed_seqinfo), n_chrom)
  expect_equal(observed_seqinfo$chrom[c(1L, n_chrom)], chrom[c(1L, n_chrom)])
  expect_equal(observed_seqinfo$length[c(1L, n_chrom)], c(1000L, 1000L))

  observed_first <- read_bwg(
    written$file,
    format = "bigwig",
    sample_names = "sampleA",
    mode = "lazy"
  )
  expect_equal(retrieve_bwg(observed_first, chrom[1L], 1L, 1000L)$value, signal$value[1L], tolerance = 1e-6)
  expect_equal(retrieve_bwg(observed_first, chrom[n_chrom], 1L, 1000L)$value, signal$value[n_chrom], tolerance = 1e-6)
})
