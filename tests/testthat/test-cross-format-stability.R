bw_test_canonical_signal <- function(x) {
  out <- data.table::copy(data.table::as.data.table(x))
  out <- out[, c("sample_id", "chrom", "start", "end", "value", "strand"), with = FALSE]
  data.table::set(out, j = "sample_id", value = as.character(out[["sample_id"]]))
  data.table::set(out, j = "chrom", value = as.character(out[["chrom"]]))
  data.table::set(out, j = "start", value = as.integer(out[["start"]]))
  data.table::set(out, j = "end", value = as.integer(out[["end"]]))
  data.table::set(out, j = "value", value = as.numeric(out[["value"]]))
  data.table::set(out, j = "strand", value = as.character(out[["strand"]]))
  data.table::setorderv(out, c("sample_id", "chrom", "start", "end", "value"))
  out <- as.data.frame(out)
  row.names(out) <- NULL
  out
}

bw_test_track <- function() {
  signal <- data.table::data.table(
    sample_id = "sampleA",
    chrom = c("chr1", "chr1", "chr1", "chr2", "chr2"),
    start = c(1L, 25L, 100L, 5L, 70L),
    end = c(10L, 39L, 125L, 20L, 90L),
    value = c(0, -2.5, 3.25, 1.5, -0.125),
    strand = "*"
  )
  bwg_track(
    samples = data.frame(sample_id = "sampleA"),
    data = signal,
    seqinfo = data.frame(
      chrom = c("chr1", "chr2"),
      length = c(500L, 300L)
    )
  )
}

test_that("BigWig zoom settings preserve full-resolution signal", {
  x <- bw_test_track()
  expected <- bw_test_canonical_signal(x$data)

  out_nozoom <- tempfile("bwtools-bigwig-nozoom-")
  out_zoom <- tempfile("bwtools-bigwig-zoom-")

  manifest_nozoom <- write_bwg(
    x,
    out_nozoom,
    format = "bigwig",
    overwrite = TRUE,
    zoom = FALSE
  )
  manifest_zoom <- write_bwg(
    x,
    out_zoom,
    format = "bigwig",
    overwrite = TRUE,
    zoom = TRUE
  )

  nozoom <- read_bwg(
    manifest_nozoom$file,
    sample_ids = "sampleA",
    mode = "memory"
  )
  zoom <- read_bwg(
    manifest_zoom$file,
    sample_ids = "sampleA",
    mode = "memory"
  )

  expect_equal(bw_test_canonical_signal(nozoom$data), expected, tolerance = 1e-6)
  expect_equal(bw_test_canonical_signal(zoom$data), expected, tolerance = 1e-6)
  expect_equal(nrow(zoominfo_bwg(nozoom)), 0L)
  expect_gt(nrow(zoominfo_bwg(zoom)), 0L)
})

test_that("WIG and bedGraph round trips preserve standardized signal", {
  x <- bw_test_track()
  expected <- bw_test_canonical_signal(x$data)

  for (format in c("wig", "bedgraph")) {
    for (compress in c(FALSE, TRUE)) {
      outdir <- tempfile(paste0("bwtools-", format, "-"))
      manifest <- write_bwg(
        x,
        outdir,
        format = format,
        overwrite = TRUE,
        compress = compress
      )
      observed <- read_bwg(
        manifest$file,
        sample_ids = "sampleA",
        mode = "memory"
      )
      expect_equal(
        bw_test_canonical_signal(observed$data),
        expected,
        tolerance = 1e-6
      )
    }
  }
})

test_that("cross-format output returns a stable manifest", {
  x <- bw_test_track()
  outdir <- tempfile("bwtools-manifest-")

  manifest <- write_bwg(
    x,
    outdir,
    format = "bedgraph",
    overwrite = TRUE
  )

  expect_s3_class(manifest, "data.table")
  expect_named(manifest, c("sample_id", "file", "format"))
  expect_identical(manifest$sample_id, "sampleA")
  expect_identical(manifest$format, "bedgraph")
  expect_true(file.exists(manifest$file))
})
