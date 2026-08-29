test_that("query coordinates use strict 1-based closed integer validation", {
  x <- bwg_track(
    samples = data.frame(sample_id = "sampleA"),
    data = data.frame(
      sample_id = "sampleA",
      chrom = "chr1",
      start = 1L,
      end = 10L,
      value = 1,
      strand = "*"
    ),
    seqinfo = data.frame(chrom = "chr1", length = 100L)
  )

  invalid_starts <- list(0L, -1L, 1.5, NA_real_, Inf, .Machine$integer.max + 1)
  for (value in invalid_starts) {
    expect_error(
      retrieve_bwg(x, "chr1", value, 10L),
      class = "bwTools_error"
    )
    expect_error(
      stats_bwg(x, "chr1", value, 10L),
      class = "bwTools_error"
    )
  }

  point <- retrieve_bwg(x, "chr1", 5L, 5L)
  expect_equal(point$start, 5L)
  expect_equal(point$end, 5L)

  expect_error(
    retrieve_bwg(x, "chr1", 10L, 9L),
    "`end`"
  )
  expect_error(
    stats_bwg(x, "chr1", 10L, 9L),
    "`end`"
  )
})

test_that("retrieve_bwg clips complete-reference queries and preserves empty schemas", {
  x <- bwg_track(
    samples = data.frame(sample_id = "sampleA"),
    data = data.frame(
      sample_id = "sampleA",
      chrom = "chr1",
      start = 90L,
      end = 100L,
      value = 2,
      strand = "*"
    ),
    seqinfo = data.frame(chrom = "chr1", length = 100L)
  )

  clipped <- retrieve_bwg(
    x,
    chrom = "chr1",
    start = 95L,
    end = 120L,
    result = "track"
  )
  expect_true(is_bwg_track(clipped))
  expect_equal(clipped$data$start, 95L)
  expect_equal(clipped$data$end, 100L)
  expect_equal(clipped$meta$retrieval_regions$start, 95L)
  expect_equal(clipped$meta$retrieval_regions$end, 100L)

  empty <- retrieve_bwg(x, "chr1", 101L, 120L)
  expect_named(
    empty,
    c("sample_id", "chrom", "start", "end", "value", "strand")
  )
  expect_equal(nrow(empty), 0L)
  expect_type(empty$start, "integer")
  expect_type(empty$value, "double")

  empty_track <- retrieve_bwg(
    x,
    chrom = "chr1",
    start = 101L,
    end = 120L,
    result = "track"
  )
  expect_true(is_bwg_track(empty_track))
  expect_equal(nrow(empty_track$data), 0L)
  expect_equal(nrow(empty_track$meta$retrieval_regions), 0L)
})

test_that("complete seqinfo rejects unknown chromosomes but incomplete text metadata does not", {
  complete <- bwg_track(
    samples = data.frame(sample_id = "sampleA"),
    data = data.frame(
      sample_id = "sampleA",
      chrom = "chr1",
      start = 1L,
      end = 10L,
      value = 1,
      strand = "*"
    ),
    seqinfo = data.frame(chrom = "chr1", length = 100L)
  )

  expect_error(
    retrieve_bwg(complete, "chr2", 1L, 10L),
    "Chromosome `chr2`",
    class = "bwTools_error"
  )
  expect_error(
    stats_bwg(complete, "chr2", 1L, 10L),
    "Chromosome `chr2`",
    class = "bwTools_error"
  )

  incomplete <- bwg_track(
    samples = data.frame(sample_id = "sampleA"),
    data = data.frame(
      sample_id = "sampleA",
      chrom = "chr1",
      start = 1L,
      end = 10L,
      value = 1,
      strand = "*"
    ),
    seqinfo = data.frame(chrom = "chr1", length = NA_integer_)
  )

  retrieved <- retrieve_bwg(incomplete, "chr2", 1L, 10L)
  expect_equal(nrow(retrieved), 0L)

  stats <- stats_bwg(
    incomplete,
    chrom = "chr2",
    start = 1L,
    end = 10L,
    n_bins = 2L,
    use_zoom = FALSE
  )
  expect_equal(nrow(stats), 2L)
  expect_true(all(is.na(stats$value)))
  expect_equal(stats$covered_bases, c(0, 0))

  mixed <- bwg_track(
    samples = data.frame(sample_id = "sampleA"),
    data = data.frame(
      sample_id = "sampleA",
      chrom = "chr1",
      start = 1L,
      end = 10L,
      value = 1,
      strand = "*"
    ),
    seqinfo = data.frame(
      chrom = c("chr1", "chr2"),
      length = c(100L, NA_integer_)
    )
  )
  expect_equal(nrow(retrieve_bwg(mixed, "chr3", 1L, 10L)), 0L)
})

test_that("stats_bwg applies chromosome-end clipping consistently in memory and lazy modes", {
  source <- bwg_track(
    samples = data.frame(sample_id = "sampleA"),
    data = data.frame(
      sample_id = "sampleA",
      chrom = "chr1",
      start = 90L,
      end = 100L,
      value = 2,
      strand = "*"
    ),
    seqinfo = data.frame(chrom = "chr1", length = 100L)
  )

  outdir <- tempfile("bwtools-boundary-stats-")
  written <- write_bwg(
    source,
    outdir = outdir,
    format = "bigwig",
    overwrite = TRUE,
    zoom = FALSE
  )
  lazy <- read_bwg(
    written$file,
    format = "bigwig",
    sample_ids = "sampleA",
    mode = "lazy"
  )

  memory_stats <- stats_bwg(
    source,
    chrom = "chr1",
    start = 95L,
    end = 120L,
    n_bins = 2L,
    stat = "mean",
    use_zoom = FALSE
  )
  lazy_stats <- stats_bwg(
    lazy,
    chrom = "chr1",
    start = 95L,
    end = 120L,
    n_bins = 2L,
    stat = "mean",
    use_zoom = FALSE
  )

  expect_equal(memory_stats$start, c(95L, 98L))
  expect_equal(memory_stats$end, c(97L, 100L))
  expect_equal(memory_stats$value, c(2, 2))
  expect_equal(lazy_stats, memory_stats, tolerance = 1e-6)

  expect_error(
    stats_bwg(
      source,
      chrom = "chr1",
      start = 95L,
      end = 120L,
      n_bins = 7L,
      use_zoom = FALSE
    ),
    "`n_bins` cannot exceed"
  )

  expected_empty_names <- c(
    "sample_id", "chrom", "start", "end", "stat", "value",
    "covered_bases", "source", "resolution"
  )
  memory_empty <- stats_bwg(
    source,
    chrom = "chr1",
    start = 101L,
    end = 120L,
    use_zoom = FALSE
  )
  lazy_empty <- stats_bwg(
    lazy,
    chrom = "chr1",
    start = 101L,
    end = 120L,
    use_zoom = FALSE
  )
  expect_named(memory_empty, expected_empty_names)
  expect_named(lazy_empty, expected_empty_names)
  expect_equal(nrow(memory_empty), 0L)
  expect_equal(nrow(lazy_empty), 0L)
})

test_that("BwgTrack validation enforces data and seqinfo coordinate consistency", {
  samples <- data.frame(sample_id = "sampleA")

  expect_error(
    bwg_track(
      samples = samples,
      data = data.frame(
        sample_id = "sampleA", chrom = "chr2", start = 1L, end = 10L,
        value = 1, strand = "*"
      ),
      seqinfo = data.frame(chrom = "chr1", length = 100L)
    ),
    "absent from `seqinfo`",
    class = "bwTools_error"
  )

  expect_error(
    bwg_track(
      samples = samples,
      data = data.frame(
        sample_id = "sampleA", chrom = "chr1", start = 95L, end = 101L,
        value = 1, strand = "*"
      ),
      seqinfo = data.frame(chrom = "chr1", length = 100L)
    ),
    "beyond `seqinfo` chromosome lengths",
    class = "bwTools_error"
  )

  valid_unknown_length <- bwg_track(
    samples = samples,
    data = data.frame(
      sample_id = "sampleA", chrom = "chr1", start = 95L, end = 101L,
      value = 1, strand = "*"
    ),
    seqinfo = data.frame(chrom = "chr1", length = NA_integer_)
  )
  expect_true(is_bwg_track(valid_unknown_length))

  tampered <- bwg_track(
    samples = samples,
    data = data.frame(
      sample_id = "sampleA", chrom = "chr1", start = 95L, end = 100L,
      value = 1, strand = "*"
    ),
    seqinfo = data.frame(chrom = "chr1", length = 100L)
  )
  tampered$data$end <- 101L
  expect_false(is_bwg_track(tampered))
  expect_error(
    validate_bwg(tampered),
    "beyond `seqinfo` chromosome lengths",
    class = "bwTools_error"
  )
})

test_that("empty in-memory signal has stable retrieve and stats behavior", {
  empty_signal <- data.frame(
    sample_id = character(),
    chrom = character(),
    start = integer(),
    end = integer(),
    value = numeric(),
    strand = character()
  )
  x <- bwg_track(
    samples = data.frame(sample_id = "sampleA"),
    data = empty_signal,
    seqinfo = data.frame(chrom = "chr1", length = 100L)
  )

  retrieved <- retrieve_bwg(x, "chr1", 1L, 20L)
  expect_equal(nrow(retrieved), 0L)
  expect_named(
    retrieved,
    c("sample_id", "chrom", "start", "end", "value", "strand")
  )

  stats <- stats_bwg(
    x,
    chrom = "chr1",
    start = 1L,
    end = 20L,
    n_bins = 2L,
    stat = "coverage",
    use_zoom = FALSE
  )
  expect_equal(nrow(stats), 2L)
  expect_true(all(is.na(stats$value)))
  expect_equal(stats$covered_bases, c(0, 0))

  expect_error(
    write_bwg(x, tempfile("bwtools-empty-write-"), format = "bigwig"),
    "No in-memory signal records"
  )
  expect_error(
    merge_bwg(x, x),
    "No signal records are available to merge"
  )
})

test_that("signal values reject non-finite input and preserve finite signed values", {
  samples <- data.frame(sample_id = "sampleA")
  for (value in list(NA_real_, NaN, Inf, -Inf)) {
    expect_error(
      bwg_track(
        samples = samples,
        data = data.frame(
          sample_id = "sampleA", chrom = "chr1", start = 1L, end = 1L,
          value = value, strand = "*"
        )
      ),
      "non-finite",
      class = "bwTools_error"
    )
  }

  signed <- bwg_track(
    samples = samples,
    data = data.frame(
      sample_id = "sampleA",
      chrom = "chr1",
      start = c(1L, 6L, 11L),
      end = c(5L, 10L, 15L),
      value = c(-3, 0, 2.5),
      strand = "*"
    ),
    seqinfo = data.frame(chrom = "chr1", length = 20L)
  )
  outdir <- tempfile("bwtools-signed-values-")
  written <- write_bwg(
    signed,
    outdir = outdir,
    format = "bigwig",
    overwrite = TRUE,
    zoom = FALSE
  )
  observed <- read_bwg(
    written$file,
    format = "bigwig",
    sample_ids = "sampleA",
    mode = "memory"
  )
  expect_equal(observed$data$value, c(-3, 0, 2.5), tolerance = 1e-6)

  too_large <- bwg_track(
    samples = samples,
    data = data.frame(
      sample_id = "sampleA", chrom = "chr1", start = 1L, end = 1L,
      value = .Machine$double.xmax, strand = "*"
    ),
    seqinfo = data.frame(chrom = "chr1", length = 20L)
  )
  expect_error(
    write_bwg(
      too_large,
      outdir = tempfile("bwtools-float32-limit-"),
      format = "bigwig",
      overwrite = TRUE,
      zoom = FALSE
    ),
    "32-bit floating-point"
  )
})

test_that("chromosome identifiers are preserved exactly across supported formats", {
  chroms <- c("1", "chr1", "Chr01", "I", "chrM", "MT", "scaffold_001")
  x <- bwg_track(
    samples = data.frame(sample_id = "sampleA"),
    data = data.frame(
      sample_id = rep("sampleA", length(chroms)),
      chrom = chroms,
      start = rep(1L, length(chroms)),
      end = rep(1L, length(chroms)),
      value = seq_along(chroms),
      strand = rep("*", length(chroms))
    ),
    seqinfo = data.frame(chrom = chroms, length = rep(10L, length(chroms)))
  )

  for (format in c("bigwig", "wig", "bedgraph")) {
    outdir <- tempfile(paste0("bwtools-chrom-names-", format, "-"))
    written <- write_bwg(
      x,
      outdir = outdir,
      format = format,
      overwrite = TRUE,
      zoom = FALSE
    )
    observed <- read_bwg(
      written$file,
      format = format,
      sample_ids = "sampleA",
      mode = "memory"
    )
    expect_setequal(unique(observed$data$chrom), chroms)
  }
})

test_that("the in-memory coordinate contract is safe at the R integer boundary", {
  max_int <- .Machine$integer.max
  x <- bwg_track(
    samples = data.frame(sample_id = "sampleA"),
    data = data.frame(
      sample_id = "sampleA",
      chrom = "chrMax",
      start = max_int,
      end = max_int,
      value = 1,
      strand = "*"
    ),
    seqinfo = data.frame(chrom = "chrMax", length = max_int)
  )
  expect_true(is_bwg_track(x))

  outdir <- tempfile("bwtools-max-coordinate-")
  written <- write_bwg(
    x,
    outdir = outdir,
    format = "bigwig",
    overwrite = TRUE,
    zoom = FALSE
  )
  lazy <- read_bwg(
    written$file,
    format = "bigwig",
    sample_ids = "sampleA",
    mode = "lazy"
  )
  observed <- retrieve_bwg(
    lazy,
    chrom = "chrMax",
    start = max_int,
    end = max_int
  )
  expect_equal(observed$start, max_int)
  expect_equal(observed$end, max_int)
  expect_equal(observed$value, 1, tolerance = 1e-6)

  expect_error(
    bwg_track(
      samples = data.frame(sample_id = "sampleA"),
      data = data.frame(
        sample_id = "sampleA",
        chrom = "chrMax",
        start = as.numeric(max_int) + 1,
        end = as.numeric(max_int) + 1,
        value = 1,
        strand = "*"
      )
    ),
    "invalid 1-based closed intervals",
    class = "bwTools_error"
  )
})

test_that("unsorted, overlapping, adjacent, and duplicate regions form a genomic union", {
  x <- bwg_track(
    samples = data.frame(sample_id = "sampleA"),
    data = data.frame(
      sample_id = "sampleA",
      chrom = "chr1",
      start = c(95L, 150L),
      end = c(105L, 160L),
      value = c(1, 2),
      strand = "*"
    ),
    seqinfo = data.frame(chrom = "chr1", length = 200L)
  )
  regions <- data.frame(
    chrom = rep("chr1", 5L),
    start = c(151L, 101L, 90L, 95L, 151L),
    end = c(170L, 120L, 100L, 99L, 170L)
  )

  observed <- retrieve_bwg(x, regions = regions, result = "track")
  expect_equal(
    observed$meta$retrieval_regions,
    data.table::data.table(
      chrom = c("chr1", "chr1"),
      start = c(90L, 151L),
      end = c(120L, 170L)
    )
  )
  expect_equal(observed$data$start, c(95L, 151L))
  expect_equal(observed$data$end, c(105L, 160L))
})

test_that("multi-sample queries reject conflicting known chromosome lengths", {
  x <- bwg_track(
    samples = data.frame(sample_id = c("sampleA", "sampleB")),
    data = data.frame(
      sample_id = c("sampleA", "sampleB"),
      chrom = c("chr1", "chr1"),
      start = c(1L, 1L),
      end = c(10L, 10L),
      value = c(1, 2),
      strand = c("*", "*")
    ),
    seqinfo = data.frame(
      sample_id = c("sampleA", "sampleB"),
      chrom = c("chr1", "chr1"),
      length = c(100L, 120L)
    )
  )

  expect_error(
    retrieve_bwg(x, "chr1", 1L, 10L),
    "conflicting chromosome lengths",
    class = "bwTools_error"
  )
  expect_error(
    stats_bwg(x, "chr1", 1L, 10L),
    "conflicting chromosome lengths",
    class = "bwTools_error"
  )

  single <- retrieve_bwg(
    x,
    "chr1",
    1L,
    10L,
    sample_ids = "sampleA"
  )
  expect_equal(single$sample_id, "sampleA")
})
