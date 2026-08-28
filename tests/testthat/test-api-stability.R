test_that("public API signatures remain stable", {
  expected <- list(
    benchmark_bwg = c(
      "file", "sample_id", "regions", "chrom", "region_sizes", "n_bins",
      "stat", "iterations", "warmup", "operations", "full_stats_max_bases",
      "write_max_bases", "write_max_intervals", "write_zoom_max_intervals",
      "write_iterations", "write_warmup", "verbose"
    ),
    bwg_track = c("samples", "data", "seqinfo", "meta"),
    detect_bwg_format = "file",
    is_bwg_track = "x",
    merge_bwg = c(
      "...", "tracks", "method", "groups", "missing", "strand",
      "drop_zero", "merge_adjacent"
    ),
    metadata_bwg = "x",
    read_bwg = c("files", "format", "sample_ids", "strand", "mode"),
    retrieve_bwg = c(
      "x", "chrom", "start", "end", "regions", "sample_ids", "strand",
      "result"
    ),
    samples_bwg = c("x", "sample_ids"),
    seqinfo_bwg = c("x", "sample_ids"),
    summary_bwg = "x",
    stats_bwg = c(
      "x", "chrom", "start", "end", "n_bins", "stat", "sample_ids",
      "use_zoom"
    ),
    validate_bwg = "x",
    write_bwg = c(
      "x", "outdir", "format", "sample_ids", "chrom_sizes", "overwrite",
      "compress", "zoom", "max_zoom_levels"
    ),
    zoominfo_bwg = c("x", "sample_ids")
  )

  for (name in names(expected)) {
    fun <- getExportedValue("bwTools", name)
    expect_identical(names(formals(fun)), expected[[name]])
  }
})

test_that("public errors inherit from bwTools_error", {
  err <- tryCatch(
    read_bwg(character()),
    error = identity
  )

  expect_s3_class(err, "bwTools_error")
  expect_match(conditionMessage(err), "`files`")

  err <- tryCatch(
    read_bwg("unused", format = "unknown"),
    error = identity
  )
  expect_s3_class(err, "bwTools_error")
  expect_match(conditionMessage(err), "`format`")

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
  err <- tryCatch(
    stats_bwg(x, chrom = "chr1", start = 1.5, end = 10L),
    error = identity
  )
  expect_s3_class(err, "bwTools_error")
  expect_match(conditionMessage(err), "`start`")
})

test_that("logical public arguments are validated consistently", {
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

  expect_error(
    stats_bwg(x, "chr1", 1L, 10L, use_zoom = 1),
    "`use_zoom` must be TRUE or FALSE",
    class = "bwTools_error"
  )
  expect_error(
    merge_bwg(x, x, drop_zero = 1),
    "`drop_zero` must be TRUE or FALSE",
    class = "bwTools_error"
  )
  expect_error(
    merge_bwg(x, x, merge_adjacent = NA),
    "`merge_adjacent` must be TRUE or FALSE",
    class = "bwTools_error"
  )

  outdir <- tempfile("bwtools-invalid-write-")
  expect_error(
    write_bwg(x, outdir, format = "bedgraph", overwrite = 1),
    "`overwrite` must be TRUE or FALSE",
    class = "bwTools_error"
  )
})
