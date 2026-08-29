test_that("mixed BigWig and bedGraph inputs preserve sample metadata", {
  bw_file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )
  bg_file <- system.file(
    "extdata", "bwtools_example.bedGraph",
    package = "bwTools", mustWork = TRUE
  )

  x <- read_bwg(
    c(bw_file, bg_file),
    sample_ids = c("rna_plus", "rna_minus"),
    strand = c("+", "-")
  )

  expect_s3_class(x, "BwgTrack")
  expect_identical(summary_bwg(x)$mode, "memory")
  expect_identical(samples_bwg(x)$sample_id, c("rna_plus", "rna_minus"))
  expect_identical(samples_bwg(x)$format, c("bigwig", "bedgraph"))
  expect_identical(samples_bwg(x)$strand, c("+", "-"))
  expect_true(all(c("rna_plus", "rna_minus") %in% unique(x$data$sample_id)))

  seqinfo <- seqinfo_bwg(x)
  expect_true(any(seqinfo$sample_id == "rna_plus" & !is.na(seqinfo$length)))
  expect_true(any(seqinfo$sample_id == "rna_minus" & is.na(seqinfo$length)))
})

test_that("BigWig-only multi-sample input remains lazy", {
  bw_file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )

  x <- read_bwg(
    c(bw_file, bw_file),
    sample_ids = c("plus", "minus"),
    strand = c("+", "-")
  )

  expect_identical(summary_bwg(x)$mode, "lazy")
  expect_null(x$data)
  expect_identical(samples_bwg(x)$sample_id, c("plus", "minus"))

  plus <- retrieve_bwg(
    x,
    chrom = "chr1",
    start = 450L,
    end = 1600L,
    strand = "+"
  )
  expect_true(nrow(plus) > 0L)
  expect_identical(unique(plus$sample_id), "plus")
  expect_true(all(plus$strand == "+"))
})

test_that("duplicate sample IDs are rejected during multi-file reads", {
  bw_file <- system.file(
    "extdata", "bwtools_example.bigwig",
    package = "bwTools", mustWork = TRUE
  )

  expect_error(
    read_bwg(
      c(bw_file, bw_file),
      sample_ids = c("sample", "sample")
    ),
    "unique",
    class = "bwTools_error"
  )
})

test_that("strand retrieval uses row-level signal for mixed-strand memory samples", {
  plus <- data.table::data.table(
    sample_id = "sample",
    chrom = "chr1",
    start = 1L,
    end = 10L,
    value = 2,
    strand = "+"
  )
  minus <- data.table::data.table(
    sample_id = "sample",
    chrom = "chr1",
    start = 20L,
    end = 30L,
    value = 4,
    strand = "-"
  )
  samples <- data.table::data.table(
    sample_id = "sample",
    file = NA_character_,
    format = "memory",
    strand = "*"
  )
  seqinfo <- data.table::data.table(
    sample_id = "sample",
    chrom = "chr1",
    length = 100L
  )
  x <- bwg_track(
    samples = samples,
    data = data.table::rbindlist(list(plus, minus)),
    seqinfo = seqinfo
  )

  plus_data <- retrieve_bwg(
    x,
    chrom = "chr1",
    start = 1L,
    end = 100L,
    strand = "+"
  )
  minus_data <- retrieve_bwg(
    x,
    chrom = "chr1",
    start = 1L,
    end = 100L,
    strand = "-"
  )
  unstranded <- retrieve_bwg(
    x,
    chrom = "chr1",
    start = 1L,
    end = 100L,
    strand = "*"
  )

  expect_equal(plus_data$start, 1L)
  expect_true(all(plus_data$strand == "+"))
  expect_equal(minus_data$start, 20L)
  expect_true(all(minus_data$strand == "-"))
  expect_equal(nrow(unstranded), 0L)

  plus_track <- retrieve_bwg(
    x,
    chrom = "chr1",
    start = 1L,
    end = 100L,
    strand = "+",
    result = "track"
  )
  expect_identical(samples_bwg(plus_track)$sample_id, "sample")
  expect_identical(samples_bwg(plus_track)$strand, "+")
  expect_true(samples_bwg(plus_track)$has_strand)
  expect_true(all(plus_track$data$strand == "+"))
})

test_that("auto merge mixed strands remains retrievable by exact strand", {
  make_track <- function(strand, start, value) {
    bwg_track(
      samples = data.table::data.table(
        sample_id = "sample",
        file = NA_character_,
        format = "memory",
        strand = strand
      ),
      data = data.table::data.table(
        sample_id = "sample",
        chrom = "chr1",
        start = as.integer(start),
        end = as.integer(start + 4L),
        value = as.numeric(value),
        strand = strand
      ),
      seqinfo = data.table::data.table(
        sample_id = "sample",
        chrom = "chr1",
        length = 100L
      )
    )
  }

  x <- merge_bwg(
    make_track("+", 1L, 2),
    make_track("-", 20L, 4)
  )

  expect_identical(samples_bwg(x)$strand, "*")
  expect_setequal(unique(x$data$strand), c("+", "-"))

  plus <- retrieve_bwg(x, "chr1", 1L, 100L, strand = "+")
  minus <- retrieve_bwg(x, "chr1", 1L, 100L, strand = "-")

  expect_equal(plus$value, 2)
  expect_true(all(plus$strand == "+"))
  expect_equal(minus$value, 4)
  expect_true(all(minus$strand == "-"))
})

test_that("multi-sample writer returns one manifest row per selected sample", {
  bg_file <- system.file(
    "extdata", "bwtools_example.bedGraph",
    package = "bwTools", mustWork = TRUE
  )
  chrom_sizes <- system.file(
    "extdata", "bwtools_example.chrom.sizes",
    package = "bwTools", mustWork = TRUE
  )

  x <- read_bwg(
    c(bg_file, bg_file),
    sample_ids = c("sample_a", "sample_b"),
    strand = c("+", "-")
  )
  outdir <- tempfile("bwtools-multisample-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE, force = TRUE), add = TRUE)

  manifest <- write_bwg(
    x,
    outdir = outdir,
    format = "bigwig",
    chrom_sizes = chrom_sizes,
    zoom = FALSE
  )

  expect_identical(manifest$sample_id, c("sample_a", "sample_b"))
  expect_equal(nrow(manifest), 2L)
  expect_true(all(file.exists(manifest$file)))
  expect_true(all(manifest$format == "bigwig"))

  observed <- read_bwg(
    manifest$file,
    sample_ids = manifest$sample_id,
    strand = c("+", "-")
  )
  expect_identical(summary_bwg(observed)$mode, "lazy")
  expect_identical(samples_bwg(observed)$sample_id, c("sample_a", "sample_b"))
})

test_that("writer rejects selected memory samples with no signal before writing", {
  samples <- data.table::data.table(
    sample_id = c("with_signal", "empty"),
    file = NA_character_,
    format = "memory",
    strand = "*"
  )
  signal <- data.table::data.table(
    sample_id = "with_signal",
    chrom = "chr1",
    start = 1L,
    end = 10L,
    value = 1,
    strand = "*"
  )
  seqinfo <- data.table::data.table(
    sample_id = c("with_signal", "empty"),
    chrom = "chr1",
    length = 100L
  )
  x <- bwg_track(samples, signal, seqinfo = seqinfo)
  outdir <- tempfile("bwtools-empty-sample-")
  dir.create(outdir)
  on.exit(unlink(outdir, recursive = TRUE, force = TRUE), add = TRUE)

  expect_error(
    write_bwg(
      x,
      outdir = outdir,
      format = "bedgraph"
    ),
    "empty",
    class = "bwTools_error"
  )
  expect_length(list.files(outdir), 0L)

  manifest <- write_bwg(
    x,
    outdir = outdir,
    format = "bedgraph",
    sample_ids = "with_signal"
  )
  expect_identical(manifest$sample_id, "with_signal")
  expect_equal(nrow(manifest), 1L)
  expect_true(file.exists(manifest$file))
})
