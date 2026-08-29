#!/usr/bin/env Rscript
# Author: Rensc
# Date: 2026-08-29
# Version: dev001
# Function: Run optional cross-implementation BigWig compatibility validation
# Input: Installed bwTools plus optional pyBigWig, rtracklayer, and UCSC utilities
# Output: compatibility-report.tsv and temporary validation artifacts

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list(output = file.path(tempdir(), "bwtools-compatibility"), strict = FALSE)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (identical(arg, "--strict")) {
      out$strict <- TRUE
      i <- i + 1L
      next
    }
    if (identical(arg, "--output")) {
      if (i == length(args)) stop("`--output` requires a directory path.")
      out$output <- args[[i + 1L]]
      i <- i + 2L
      next
    }
    if (startsWith(arg, "--output=")) {
      out$output <- sub("^--output=", "", arg)
      i <- i + 1L
      next
    }
    stop("Unknown argument: ", arg)
  }
  out
}

config <- parse_args(args)
dir.create(config$output, recursive = TRUE, showWarnings = FALSE)
output_dir <- normalizePath(config$output, winslash = "/", mustWork = TRUE)

if (!requireNamespace("bwTools", quietly = TRUE)) {
  stop("bwTools must be installed before running compatibility validation.")
}
if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("data.table is required by bwTools and must be installed.")
}

package_root <- system.file(package = "bwTools")
compat_dir <- system.file("compatibility", package = "bwTools")
bedgraph_file <- system.file("extdata", "bwtools_example.bedGraph", package = "bwTools")
chrom_sizes_file <- system.file("extdata", "bwtools_example.chrom.sizes", package = "bwTools")
python_bridge <- file.path(compat_dir, "pybigwig-bridge.001.py")

required_files <- c(package_root, compat_dir, bedgraph_file, chrom_sizes_file, python_bridge)
if (any(!nzchar(required_files)) || any(!file.exists(required_files))) {
  stop("Installed bwTools compatibility assets are incomplete.")
}

report <- list()
add_result <- function(implementation, direction, status, version = NA_character_, detail = "", artifact = NA_character_) {
  report[[length(report) + 1L]] <<- data.table::data.table(
    implementation = as.character(implementation),
    direction = as.character(direction),
    status = as.character(status),
    version = as.character(version),
    detail = as.character(detail),
    artifact = as.character(artifact)
  )
  invisible(NULL)
}

run_process <- function(command, args = character()) {
  stdout <- tempfile("bwtools-compat-stdout-")
  stderr <- tempfile("bwtools-compat-stderr-")
  on.exit(unlink(c(stdout, stderr)), add = TRUE)
  status <- suppressWarnings(system2(command, args = args, stdout = stdout, stderr = stderr))
  list(
    status = if (is.null(status)) 0L else as.integer(status),
    stdout = if (file.exists(stdout)) readLines(stdout, warn = FALSE) else character(),
    stderr = if (file.exists(stderr)) readLines(stderr, warn = FALSE) else character()
  )
}

find_python_with_pybigwig <- function() {
  candidates <- unique(c(Sys.which("python"), Sys.which("python3")))
  candidates <- unname(candidates[nzchar(candidates)])
  for (candidate in candidates) {
    check <- run_process(candidate, c(shQuote(python_bridge), "check"))
    if (identical(check$status, 0L)) {
      version <- if (length(check$stdout)) tail(check$stdout, 1L) else "unknown"
      return(list(command = candidate, version = version))
    }
  }
  NULL
}

canonical_track <- bwTools::read_bwg(
  bedgraph_file,
  format = "bedgraph",
  sample_ids = "reference",
  mode = "memory"
)
chrom_sizes <- data.table::fread(
  chrom_sizes_file,
  header = FALSE,
  col.names = c("chrom", "length"),
  showProgress = FALSE
)

canonical_dense <- function(signal, chrom_sizes) {
  signal <- data.table::as.data.table(signal)
  out <- setNames(vector("list", nrow(chrom_sizes)), chrom_sizes$chrom)
  for (i in seq_len(nrow(chrom_sizes))) {
    values <- rep(NA_real_, chrom_sizes$length[[i]])
    idx <- which(signal[["chrom"]] == chrom_sizes$chrom[[i]])
    if (length(idx)) {
      for (j in idx) {
        values[signal[["start"]][[j]]:signal[["end"]][[j]]] <- signal[["value"]][[j]]
      }
    }
    out[[i]] <- values
  }
  out
}

expected_dense <- canonical_dense(canonical_track$data, chrom_sizes)

assert_dense_equal <- function(observed, label, tolerance = 1e-5) {
  observed_dense <- canonical_dense(observed, chrom_sizes)
  if (!identical(names(observed_dense), names(expected_dense))) {
    stop(label, ": chromosome names differ.")
  }
  for (chrom in names(expected_dense)) {
    x <- observed_dense[[chrom]]
    y <- expected_dense[[chrom]]
    if (!identical(is.na(x), is.na(y))) {
      stop(label, ": missing-signal pattern differs on ", chrom, ".")
    }
    keep <- !is.na(y)
    if (any(abs(x[keep] - y[keep]) > tolerance)) {
      stop(label, ": signal values differ on ", chrom, ".")
    }
  }
  invisible(TRUE)
}

bwtools_outdir <- file.path(output_dir, "bwtools")
nozoom_outdir <- file.path(bwtools_outdir, "nozoom")
zoom_outdir <- file.path(bwtools_outdir, "zoom")
dir.create(nozoom_outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(zoom_outdir, recursive = TRUE, showWarnings = FALSE)
manifest_nozoom <- bwTools::write_bwg(
  canonical_track,
  nozoom_outdir,
  format = "bigwig",
  chrom_sizes = chrom_sizes,
  overwrite = TRUE,
  zoom = FALSE
)
nozoom_file <- manifest_nozoom$file[[1L]]

manifest_zoom <- bwTools::write_bwg(
  canonical_track,
  zoom_outdir,
  format = "bigwig",
  chrom_sizes = chrom_sizes,
  overwrite = TRUE,
  zoom = TRUE
)
zoom_file <- manifest_zoom$file[[1L]]

for (path in c(nozoom_file, zoom_file)) {
  observed <- bwTools::read_bwg(path, sample_ids = "reference", mode = "memory")
  assert_dense_equal(observed$data, basename(path))
}
add_result(
  "bwTools",
  "self-roundtrip",
  "PASS",
  as.character(utils::packageVersion("bwTools")),
  "Canonical bedGraph signal preserved with zoom disabled and enabled.",
  bwtools_outdir
)

# pyBigWig: validate bwTools files, then generate an external file for bwTools.
python <- find_python_with_pybigwig()
if (is.null(python)) {
  add_result("pyBigWig", "bidirectional", "SKIP", detail = "No Python interpreter with importable pyBigWig was found.")
} else {
  py_failed <- FALSE
  py_detail <- character()
  for (path in c(nozoom_file, zoom_file)) {
    result <- run_process(
      python$command,
      c(shQuote(python_bridge), "validate", shQuote(path), shQuote(bedgraph_file), shQuote(chrom_sizes_file))
    )
    if (!identical(result$status, 0L)) {
      py_failed <- TRUE
      py_detail <- c(py_detail, result$stderr, result$stdout)
    }
  }
  py_output <- file.path(output_dir, "pybigwig-reference.bigwig")
  if (!py_failed) {
    result <- run_process(
      python$command,
      c(shQuote(python_bridge), "write", shQuote(py_output), shQuote(bedgraph_file), shQuote(chrom_sizes_file))
    )
    if (!identical(result$status, 0L)) {
      py_failed <- TRUE
      py_detail <- c(py_detail, result$stderr, result$stdout)
    }
  }
  if (!py_failed) {
    tryCatch({
      observed <- bwTools::read_bwg(py_output, sample_ids = "reference", mode = "memory")
      assert_dense_equal(observed$data, "pyBigWig -> bwTools")
    }, error = function(e) {
      py_failed <<- TRUE
      py_detail <<- c(py_detail, conditionMessage(e))
    })
  }
  add_result(
    "pyBigWig",
    "bidirectional",
    if (py_failed) "FAIL" else "PASS",
    python$version,
    paste(py_detail, collapse = " | "),
    if (file.exists(py_output)) py_output else NA_character_
  )
}

# rtracklayer: import bwTools output, then export a reference BigWig for bwTools.
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  add_result("rtracklayer", "bidirectional", "SKIP", detail = "rtracklayer is not installed.")
} else {
  rt_failed <- FALSE
  rt_detail <- character()
  rt_output <- file.path(output_dir, "rtracklayer-reference.bigwig")
  tryCatch({
    for (path in c(nozoom_file, zoom_file)) {
      gr <- rtracklayer::import(path, format = "BigWig")
      gr_df <- as.data.frame(gr)
      observed <- data.table::data.table(
        chrom = as.character(gr_df$seqnames),
        start = as.integer(gr_df$start),
        end = as.integer(gr_df$end),
        value = as.numeric(gr_df$score)
      )
      assert_dense_equal(observed, paste0("bwTools -> rtracklayer: ", basename(path)))
    }

    signal <- data.table::as.data.table(canonical_track$data)
    gr <- GenomicRanges::GRanges(
      seqnames = signal$chrom,
      ranges = IRanges::IRanges(start = signal$start, end = signal$end),
      score = signal$value
    )
    GenomeInfoDb::seqlengths(gr) <- stats::setNames(chrom_sizes$length, chrom_sizes$chrom)
    rtracklayer::export(gr, rt_output, format = "BigWig")
    observed <- bwTools::read_bwg(rt_output, sample_ids = "reference", mode = "memory")
    assert_dense_equal(observed$data, "rtracklayer -> bwTools")
  }, error = function(e) {
    rt_failed <<- TRUE
    rt_detail <<- c(rt_detail, conditionMessage(e))
  })
  add_result(
    "rtracklayer",
    "bidirectional",
    if (rt_failed) "FAIL" else "PASS",
    as.character(utils::packageVersion("rtracklayer")),
    paste(rt_detail, collapse = " | "),
    if (file.exists(rt_output)) rt_output else NA_character_
  )
}

# UCSC utilities: validate via bigWigInfo/bigWigToBedGraph and generate via bedGraphToBigWig.
ucsc_tools <- c(
  bigWigInfo = Sys.which("bigWigInfo"),
  bigWigToBedGraph = Sys.which("bigWigToBedGraph"),
  bedGraphToBigWig = Sys.which("bedGraphToBigWig")
)
if (any(!nzchar(ucsc_tools))) {
  missing <- names(ucsc_tools)[!nzchar(ucsc_tools)]
  add_result(
    "UCSC",
    "bidirectional",
    "SKIP",
    detail = paste0("Missing utilities: ", paste(missing, collapse = ", "), ".")
  )
} else {
  ucsc_failed <- FALSE
  ucsc_detail <- character()
  ucsc_output <- file.path(output_dir, "ucsc-reference.bigwig")
  tryCatch({
    for (path in c(nozoom_file, zoom_file)) {
      info <- run_process(ucsc_tools[["bigWigInfo"]], shQuote(path))
      if (!identical(info$status, 0L)) {
        stop("bigWigInfo rejected ", basename(path), ": ", paste(c(info$stderr, info$stdout), collapse = " | "))
      }
      bg_out <- tempfile("bwtools-ucsc-", fileext = ".bedGraph")
      converted <- run_process(ucsc_tools[["bigWigToBedGraph"]], c(shQuote(path), shQuote(bg_out)))
      if (!identical(converted$status, 0L)) {
        stop("bigWigToBedGraph failed for ", basename(path), ": ", paste(c(converted$stderr, converted$stdout), collapse = " | "))
      }
      observed <- bwTools::read_bwg(bg_out, format = "bedgraph", sample_ids = "reference", mode = "memory")
      assert_dense_equal(observed$data, paste0("bwTools -> UCSC: ", basename(path)))
    }

    generated <- run_process(
      ucsc_tools[["bedGraphToBigWig"]],
      c(shQuote(bedgraph_file), shQuote(chrom_sizes_file), shQuote(ucsc_output))
    )
    if (!identical(generated$status, 0L)) {
      stop("bedGraphToBigWig failed: ", paste(c(generated$stderr, generated$stdout), collapse = " | "))
    }
    observed <- bwTools::read_bwg(ucsc_output, sample_ids = "reference", mode = "memory")
    assert_dense_equal(observed$data, "UCSC -> bwTools")
  }, error = function(e) {
    ucsc_failed <<- TRUE
    ucsc_detail <<- c(ucsc_detail, conditionMessage(e))
  })
  version_result <- run_process(ucsc_tools[["bigWigInfo"]], "-version")
  version <- paste(c(version_result$stdout, version_result$stderr), collapse = " ")
  add_result(
    "UCSC",
    "bidirectional",
    if (ucsc_failed) "FAIL" else "PASS",
    if (nzchar(version)) version else "unknown",
    paste(ucsc_detail, collapse = " | "),
    if (file.exists(ucsc_output)) ucsc_output else NA_character_
  )
}

report_dt <- data.table::rbindlist(report, use.names = TRUE, fill = TRUE)
report_file <- file.path(output_dir, "compatibility-report.tsv")
data.table::fwrite(report_dt, report_file, sep = "\t")
print(report_dt)
cat("\nCompatibility report: ", report_file, "\n", sep = "")

external <- report_dt[implementation != "bwTools"]
if (any(external$status == "FAIL")) {
  quit(save = "no", status = 1L)
}
if (isTRUE(config$strict) && (nrow(external) < 3L || any(external$status != "PASS"))) {
  quit(save = "no", status = 2L)
}
quit(save = "no", status = 0L)
