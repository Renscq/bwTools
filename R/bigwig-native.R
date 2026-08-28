# Author: Rensc
# Date: 2026-08-29
# Version: dev004
# Function: Read local BigWig files using native R binary I/O
# Input: Local BigWig file paths and genomic intervals
# Output: Chromosome metadata and 1-based closed signal intervals

bw_cache_key <- function(file) {
  info <- file.info(file)
  paste(file, as.numeric(info$size[1L]), as.numeric(info$mtime[1L]), sep = "|")
}

bw_read_exact <- function(con, offset, n) {
  offset <- as.numeric(offset)[1L]
  n <- as.numeric(n)[1L]
  if (!is.finite(offset) || offset < 0 || offset > .BWT_MAX_EXACT_UINT64) {
    bw_stop("Invalid or unsupported BigWig file offset.")
  }
  if (!is.finite(n) || n < 0 || n > .Machine$integer.max) {
    bw_stop("Invalid or unsupported BigWig read size.")
  }
  seek(con, where = offset, origin = "start")
  out <- readBin(con, what = "raw", n = as.integer(n), size = 1L)
  if (length(out) != as.integer(n)) {
    bw_stop(paste0("Unexpected end of BigWig file at offset ", format(offset, scientific = FALSE), "."))
  }
  out
}

bw_uint <- function(x, offset, size, endian = "little") {
  offset <- as.integer(offset)[1L]
  size <- as.integer(size)[1L]
  idx <- seq.int(offset + 1L, offset + size)
  if (offset < 0L || size < 1L || max(idx) > length(x)) {
    bw_stop("Binary integer read exceeds the available buffer.")
  }
  bytes <- as.numeric(x[idx])
  if (identical(endian, "big")) {
    bytes <- rev(bytes)
  }
  sum(bytes * 256^(seq_along(bytes) - 1L))
}

bw_u8 <- function(x, offset) {
  as.integer(bw_uint(x, offset, 1L, "little"))
}

bw_u16 <- function(x, offset, endian) {
  as.integer(bw_uint(x, offset, 2L, endian))
}

bw_u32 <- function(x, offset, endian) {
  bw_uint(x, offset, 4L, endian)
}

bw_u64 <- function(x, offset, endian) {
  out <- bw_uint(x, offset, 8L, endian)
  if (!is.finite(out) || out > .BWT_MAX_EXACT_UINT64) {
    bw_stop("A 64-bit BigWig offset exceeds the exact integer range supported by native R numeric values.")
  }
  out
}

bw_uint_vector <- function(x, offsets, size, endian = "little") {
  offsets <- as.numeric(offsets)
  size <- as.integer(size)[1L]
  if (length(offsets) == 0L) {
    return(numeric())
  }
  if (size < 1L || any(!is.finite(offsets)) || any(offsets < 0)) {
    bw_stop("Invalid vectorized binary integer read.")
  }
  if (any(offsets + size > length(x))) {
    bw_stop("Vectorized binary integer read exceeds the available buffer.")
  }
  idx <- rep(offsets, each = size) + rep(seq.int(0L, size - 1L), times = length(offsets)) + 1
  bytes <- matrix(as.numeric(x[idx]), nrow = length(offsets), ncol = size, byrow = TRUE)
  if (identical(endian, "big")) {
    bytes <- bytes[, seq.int(size, 1L), drop = FALSE]
  }
  drop(bytes %*% 256^(seq_len(size) - 1L))
}

bw_u32_vector <- function(x, offsets, endian) {
  bw_uint_vector(x, offsets, 4L, endian)
}

bw_float32_vector <- function(x, offsets, endian) {
  bits <- bw_u32_vector(x, offsets, endian)
  sign_bit <- floor(bits / 2^31)
  exponent <- floor((bits %% 2^31) / 2^23)
  fraction <- bits %% 2^23
  sign_value <- ifelse(sign_bit == 0, 1, -1)
  out <- sign_value * (1 + fraction / 2^23) * 2^(exponent - 127)
  subnormal <- exponent == 0
  out[subnormal] <- sign_value[subnormal] * (fraction[subnormal] / 2^23) * 2^-126
  infinite <- exponent == 255 & fraction == 0
  nan_value <- exponent == 255 & fraction != 0
  out[infinite] <- sign_value[infinite] * Inf
  out[nan_value] <- NaN
  out
}



bw_double <- function(x, offset, endian) {
  offset <- as.integer(offset)[1L]
  idx <- seq.int(offset + 1L, offset + 8L)
  if (offset < 0L || max(idx) > length(x)) {
    bw_stop("Binary double read exceeds the available buffer.")
  }
  con <- rawConnection(x[idx], open = "rb")
  on.exit(close(con), add = TRUE)
  readBin(con, what = numeric(), n = 1L, size = 8L, endian = endian)
}

bw_fixed_string <- function(x) {
  zero <- which(as.integer(x) == 0L)
  if (length(zero) > 0L) {
    end <- zero[1L] - 1L
    x <- if (end > 0L) x[seq_len(end)] else raw()
  }
  if (length(x) == 0L) {
    return("")
  }
  rawToChar(x)
}

bw_detect_endian <- function(header_raw) {
  if (bw_u32(header_raw, 0L, "little") == .BWT_BIGWIG_MAGIC) {
    return("little")
  }
  if (bw_u32(header_raw, 0L, "big") == .BWT_BIGWIG_MAGIC) {
    return("big")
  }
  bw_stop("The input file is not a valid BigWig file.")
}

bw_read_header <- function(con) {
  x <- bw_read_exact(con, 0, .BWT_HEADER_SIZE)
  endian <- bw_detect_endian(x)
  list(
    endian = endian,
    version = bw_u16(x, 4L, endian),
    n_levels = bw_u16(x, 6L, endian),
    chrom_tree_offset = bw_u64(x, 8L, endian),
    data_offset = bw_u64(x, 16L, endian),
    index_offset = bw_u64(x, 24L, endian),
    field_count = bw_u16(x, 32L, endian),
    defined_field_count = bw_u16(x, 34L, endian),
    sql_offset = bw_u64(x, 36L, endian),
    summary_offset = bw_u64(x, 44L, endian),
    uncompress_buf_size = bw_u32(x, 52L, endian),
    extension_offset = bw_u64(x, 56L, endian)
  )
}

bw_read_zoom_headers <- function(con, header) {
  n_levels <- as.integer(header$n_levels)
  if (n_levels < 1L) {
    return(data.table::data.table(
      level = integer(),
      data_offset = numeric(),
      index_offset = numeric()
    ))
  }
  x <- bw_read_exact(con, .BWT_HEADER_SIZE, as.numeric(n_levels) * 24L)
  offsets <- seq.int(0L, by = 24L, length.out = n_levels)
  data.table::data.table(
    level = as.numeric(vapply(offsets, function(z) bw_u32(x, z, header$endian), numeric(1L))),
    data_offset = vapply(offsets, function(z) bw_u64(x, z + 8L, header$endian), numeric(1L)),
    index_offset = vapply(offsets, function(z) bw_u64(x, z + 16L, header$endian), numeric(1L))
  )
}

bw_read_total_summary <- function(con, header) {
  if (header$summary_offset <= 0) {
    return(NULL)
  }
  x <- bw_read_exact(con, header$summary_offset, .BWT_SUMMARY_SIZE)
  list(
    n_bases_covered = bw_u64(x, 0L, header$endian),
    min_value = bw_double(x, 8L, header$endian),
    max_value = bw_double(x, 16L, header$endian),
    sum_data = bw_double(x, 24L, header$endian),
    sum_squared = bw_double(x, 32L, header$endian)
  )
}

bw_read_chrom_node <- function(con, offset, key_size, endian, chrom, chrom_length, depth = 0L) {
  if (depth > 128L) {
    bw_stop("Chromosome B+ tree exceeds the supported recursion depth.")
  }
  node_header <- bw_read_exact(con, offset, 4L)
  is_leaf <- bw_u8(node_header, 0L) != 0L
  n_children <- bw_u16(node_header, 2L, endian)
  entry_size <- key_size + 8L
  if (n_children == 0L) {
    return(list(chrom = chrom, length = chrom_length))
  }
  entries <- bw_read_exact(con, offset + 4, as.numeric(n_children) * entry_size)
  for (i in seq_len(n_children)) {
    base <- (i - 1L) * entry_size
    key_raw <- entries[seq.int(base + 1L, base + key_size)]
    if (is_leaf) {
      tid <- bw_u32(entries, base + key_size, endian)
      seq_len_value <- bw_u32(entries, base + key_size + 4L, endian)
      r_index <- as.integer(tid + 1)
      if (r_index < 1L || r_index > length(chrom)) {
        bw_stop("Invalid chromosome ID in BigWig chromosome tree.")
      }
      chrom[r_index] <- bw_fixed_string(key_raw)
      chrom_length[r_index] <- seq_len_value
    } else {
      child_offset <- bw_u64(entries, base + key_size, endian)
      child <- bw_read_chrom_node(con, child_offset, key_size, endian, chrom, chrom_length, depth + 1L)
      chrom <- child$chrom
      chrom_length <- child$length
    }
  }
  list(chrom = chrom, length = chrom_length)
}

bw_read_chrom_tree <- function(con, header) {
  x <- bw_read_exact(con, header$chrom_tree_offset, 32L)
  endian <- header$endian
  if (bw_u32(x, 0L, endian) != .BWT_CIRTREE_MAGIC) {
    bw_stop("Invalid BigWig chromosome B+ tree magic number.")
  }
  key_size <- bw_u32(x, 8L, endian)
  value_size <- bw_u32(x, 12L, endian)
  item_count <- bw_u64(x, 16L, endian)
  if (key_size < 1 || key_size > .Machine$integer.max) {
    bw_stop("Invalid chromosome key size in BigWig file.")
  }
  if (value_size != 8) {
    bw_stop("Unsupported chromosome value size in BigWig file.")
  }
  if (item_count > .Machine$integer.max) {
    bw_stop("Unsupported number of chromosomes in BigWig file.")
  }
  n <- as.integer(item_count)
  out <- bw_read_chrom_node(
    con = con,
    offset = header$chrom_tree_offset + 32,
    key_size = as.integer(key_size),
    endian = endian,
    chrom = rep(NA_character_, n),
    chrom_length = rep(NA_real_, n)
  )
  if (anyNA(out$chrom) || any(!nzchar(out$chrom))) {
    bw_stop("Failed to resolve all chromosomes from the BigWig chromosome tree.")
  }
  data.table::data.table(tid = seq_len(n) - 1L, chrom = out$chrom, length = out$length)
}

bw_metadata <- function(file, use_cache = TRUE) {
  file <- bw_validate_local_file(file)
  key <- bw_cache_key(file)
  if (isTRUE(use_cache) && exists(key, envir = .bwtools_cache, inherits = FALSE)) {
    return(get(key, envir = .bwtools_cache, inherits = FALSE))
  }
  con <- base::file(file, open = "rb")
  on.exit(close(con), add = TRUE)
  header <- bw_read_header(con)
  zoom_levels <- bw_read_zoom_headers(con, header)
  total_summary <- bw_read_total_summary(con, header)
  chromosomes <- bw_read_chrom_tree(con, header)
  out <- list(
    file = file,
    header = header,
    zoom_levels = zoom_levels,
    total_summary = total_summary,
    chromosomes = chromosomes
  )
  if (isTRUE(use_cache)) {
    assign(key, out, envir = .bwtools_cache)
  }
  out
}

bw_read_index_header <- function(con, header, index_offset = header$index_offset) {
  index_offset <- as.numeric(index_offset)[1L]
  if (!is.finite(index_offset) || index_offset <= 0) {
    bw_stop("The BigWig file does not contain the requested R-tree index.")
  }
  x <- bw_read_exact(con, index_offset, .BWT_INDEX_HEADER_SIZE)
  endian <- header$endian
  if (bw_u32(x, 0L, endian) != .BWT_INDEX_MAGIC) {
    bw_stop("Invalid BigWig R-tree index magic number.")
  }
  list(
    block_size = bw_u32(x, 4L, endian),
    n_items = bw_u64(x, 8L, endian),
    chrom_start = bw_u32(x, 16L, endian),
    base_start = bw_u32(x, 20L, endian),
    chrom_end = bw_u32(x, 24L, endian),
    base_end = bw_u32(x, 28L, endian),
    index_size = bw_u64(x, 32L, endian),
    items_per_slot = bw_u32(x, 40L, endian),
    root_offset = index_offset + .BWT_INDEX_HEADER_SIZE
  )
}

bw_index_entry_overlaps <- function(chrom_start, base_start, chrom_end, base_end, tid, query_start, query_end) {
  if (tid < chrom_start || tid > chrom_end) {
    return(FALSE)
  }
  if (chrom_start != chrom_end) {
    if (tid == chrom_start && base_start >= query_end) return(FALSE)
    if (tid == chrom_end && base_end <= query_start) return(FALSE)
    return(TRUE)
  }
  !(base_start >= query_end || base_end <= query_start)
}

bw_collect_blocks <- function(con, node_offset, endian, tid, query_start, query_end, depth = 0L) {
  if (depth > 128L) {
    bw_stop("BigWig R-tree exceeds the supported recursion depth.")
  }
  node_header <- bw_read_exact(con, node_offset, 4L)
  is_leaf <- bw_u8(node_header, 0L) != 0L
  n_children <- bw_u16(node_header, 2L, endian)
  if (n_children == 0L) {
    return(data.table::data.table(offset = numeric(), size = numeric()))
  }
  entry_size <- if (is_leaf) 32L else 24L
  entries <- bw_read_exact(con, node_offset + 4, as.numeric(n_children) * entry_size)
  out <- vector("list", n_children)
  out_n <- 0L
  for (i in seq_len(n_children)) {
    base <- (i - 1L) * entry_size
    chrom_start <- bw_u32(entries, base, endian)
    base_start <- bw_u32(entries, base + 4L, endian)
    chrom_end <- bw_u32(entries, base + 8L, endian)
    base_end <- bw_u32(entries, base + 12L, endian)
    data_offset <- bw_u64(entries, base + 16L, endian)
    if (tid < chrom_start) break
    if (!bw_index_entry_overlaps(chrom_start, base_start, chrom_end, base_end, tid, query_start, query_end)) next
    if (is_leaf) {
      data_size <- bw_u64(entries, base + 24L, endian)
      out_n <- out_n + 1L
      out[[out_n]] <- data.table::data.table(offset = data_offset, size = data_size)
    } else {
      child <- bw_collect_blocks(con, data_offset, endian, tid, query_start, query_end, depth + 1L)
      if (nrow(child) > 0L) {
        out_n <- out_n + 1L
        out[[out_n]] <- child
      }
    }
  }
  if (out_n == 0L) {
    return(data.table::data.table(offset = numeric(), size = numeric()))
  }
  ans <- data.table::rbindlist(out[seq_len(out_n)])
  ans <- unique(ans, by = c("offset", "size"))
  ans[]
}

bw_decompress_block <- function(x) {
  tryCatch(
    base::memDecompress(x, type = "gzip"),
    error = function(e) bw_stop(paste0("Failed to decompress a zlib-compressed BigWig data block: ", conditionMessage(e)))
  )
}

bw_decode_block_vectors <- function(x, endian, tid, query_start, query_end) {
  if (length(x) < .BWT_DATA_HEADER_SIZE) {
    bw_stop("Truncated BigWig data block header.")
  }
  block_tid <- bw_u32(x, 0L, endian)
  block_start <- bw_u32(x, 4L, endian)
  block_step <- bw_u32(x, 12L, endian)
  block_span <- bw_u32(x, 16L, endian)
  block_type <- bw_u8(x, 20L)
  n_items <- bw_u16(x, 22L, endian)
  if (block_tid != tid || n_items == 0L) {
    return(list(start0 = numeric(), end0 = numeric(), value = numeric()))
  }

  item_size <- switch(
    as.character(block_type),
    "1" = 12L,
    "2" = 8L,
    "3" = 4L,
    bw_stop("Unsupported BigWig data block type.")
  )
  required_size <- .BWT_DATA_HEADER_SIZE + as.numeric(n_items) * item_size
  if (required_size > length(x)) {
    bw_stop("Truncated BigWig data block payload.")
  }

  item_offsets <- .BWT_DATA_HEADER_SIZE + seq.int(0, n_items - 1L) * item_size
  if (block_type == 1L) {
    start0 <- bw_u32_vector(x, item_offsets, endian)
    end0 <- bw_u32_vector(x, item_offsets + 4L, endian)
    value <- bw_float32_vector(x, item_offsets + 8L, endian)
  } else if (block_type == 2L) {
    start0 <- bw_u32_vector(x, item_offsets, endian)
    end0 <- start0 + block_span
    value <- bw_float32_vector(x, item_offsets + 4L, endian)
  } else {
    start0 <- block_start + seq.int(0, n_items - 1L) * block_step
    end0 <- start0 + block_span
    value <- bw_float32_vector(x, item_offsets, endian)
  }

  keep <- end0 > query_start & start0 < query_end
  if (!any(keep)) {
    return(list(start0 = numeric(), end0 = numeric(), value = numeric()))
  }
  start0 <- pmax(start0[keep], query_start)
  end0 <- pmin(end0[keep], query_end)
  value <- as.numeric(value[keep])
  valid <- end0 > start0
  if (!any(valid)) {
    return(list(start0 = numeric(), end0 = numeric(), value = numeric()))
  }
  list(
    start0 = as.numeric(start0[valid]),
    end0 = as.numeric(end0[valid]),
    value = value[valid]
  )
}

bw_decode_block <- function(x, endian, tid, query_start, query_end, chrom) {
  decoded <- bw_decode_block_vectors(
    x = x,
    endian = endian,
    tid = tid,
    query_start = query_start,
    query_end = query_end
  )
  if (length(decoded$start0) < 1L) {
    return(bw_empty_signal())
  }

  start1 <- decoded$start0 + 1
  end1 <- decoded$end0
  if (any(start1 > .Machine$integer.max) || any(end1 > .Machine$integer.max)) {
    bw_stop("Queried coordinates exceed the supported integer coordinate range.")
  }
  data.table::data.table(
    chrom = rep(chrom, length(start1)),
    start = as.integer(start1),
    end = as.integer(end1),
    value = decoded$value
  )
}

bw_decode_zoom_block <- function(x, endian, tid, query_start, query_end, chrom) {
  if (length(x) %% .BWT_ZOOM_RECORD_SIZE != 0L) {
    bw_stop("Truncated or malformed BigWig zoom data block.")
  }
  n <- length(x) %/% .BWT_ZOOM_RECORD_SIZE
  if (n < 1L) {
    return(data.table::data.table())
  }
  offsets <- seq.int(0L, by = .BWT_ZOOM_RECORD_SIZE, length.out = n)
  record_tid <- bw_u32_vector(x, offsets, endian)
  start0 <- bw_u32_vector(x, offsets + 4L, endian)
  end0 <- bw_u32_vector(x, offsets + 8L, endian)
  valid_count <- bw_u32_vector(x, offsets + 12L, endian)
  min_value <- bw_float32_vector(x, offsets + 16L, endian)
  max_value <- bw_float32_vector(x, offsets + 20L, endian)
  sum_data <- bw_float32_vector(x, offsets + 24L, endian)
  sum_squared <- bw_float32_vector(x, offsets + 28L, endian)
  keep <- record_tid == tid & end0 > query_start & start0 < query_end
  if (!any(keep)) {
    return(data.table::data.table())
  }
  data.table::data.table(
    chrom = rep(chrom, sum(keep)),
    start0 = as.numeric(start0[keep]),
    end0 = as.numeric(end0[keep]),
    valid_count = as.numeric(valid_count[keep]),
    min_value = as.numeric(min_value[keep]),
    max_value = as.numeric(max_value[keep]),
    sum_data = as.numeric(sum_data[keep]),
    sum_squared = as.numeric(sum_squared[keep])
  )
}

bw_choose_zoom_level <- function(metadata, start, end, n_bins) {
  z <- metadata$zoom_levels
  if (is.null(z) || nrow(z) < 1L) {
    return(NULL)
  }
  bases_per_bin <- (as.numeric(end) - as.numeric(start) + 1) / as.numeric(n_bins)
  target <- bases_per_bin / 2
  eligible <- z[z[["level"]] <= target]
  if (nrow(eligible) < 1L) {
    return(NULL)
  }
  eligible[which.max(eligible[["level"]])]
}

bw_bigwig_zoom_query <- function(file, chrom, start, end, level = NULL, n_bins = 1L) {
  file <- bw_validate_local_file(file)
  metadata <- bw_metadata(file)
  chrom_idx <- match(chrom, metadata$chromosomes$chrom)
  if (is.na(chrom_idx)) {
    bw_stop(paste0("Chromosome `", chrom, "` was not found in the BigWig file."))
  }
  query_start <- as.numeric(start) - 1
  query_end <- min(as.numeric(end), metadata$chromosomes$length[chrom_idx])
  if (query_start >= query_end) {
    return(data.table::data.table())
  }
  if (is.null(level)) {
    zoom <- bw_choose_zoom_level(metadata, start, end, n_bins)
  } else {
    level <- as.numeric(level)[1L]
    zoom_idx <- match(level, metadata$zoom_levels$level)
    zoom <- if (is.na(zoom_idx)) NULL else metadata$zoom_levels[zoom_idx]
  }
  if (is.null(zoom) || nrow(zoom) < 1L) {
    return(data.table::data.table())
  }

  con <- base::file(file, open = "rb")
  on.exit(close(con), add = TRUE)
  index <- bw_read_index_header(con, metadata$header, zoom$index_offset[1L])
  blocks <- bw_collect_blocks(
    con = con,
    node_offset = index$root_offset,
    endian = metadata$header$endian,
    tid = metadata$chromosomes$tid[chrom_idx],
    query_start = query_start,
    query_end = query_end
  )
  if (nrow(blocks) < 1L) {
    return(data.table::data.table())
  }
  out <- vector("list", nrow(blocks))
  out_n <- 0L
  for (i in seq_len(nrow(blocks))) {
    block <- bw_read_exact(con, blocks$offset[i], blocks$size[i])
    if (metadata$header$uncompress_buf_size > 0) {
      block <- bw_decompress_block(block)
    }
    decoded <- bw_decode_zoom_block(
      block, metadata$header$endian, metadata$chromosomes$tid[chrom_idx],
      query_start, query_end, chrom
    )
    if (nrow(decoded) > 0L) {
      out_n <- out_n + 1L
      out[[out_n]] <- decoded
    }
  }
  if (out_n < 1L) {
    return(data.table::data.table())
  }
  ans <- data.table::rbindlist(out[seq_len(out_n)], use.names = TRUE)
  data.table::setorderv(ans, c("start0", "end0"))
  attr(ans, "zoom_level") <- as.numeric(zoom$level[1L])
  ans[]
}

bw_bigwig_seqinfo <- function(file) {
  metadata <- bw_metadata(file)
  x <- metadata$chromosomes
  lengths <- rep(NA_integer_, nrow(x))
  safe <- is.finite(x$length) & x$length <= .Machine$integer.max
  lengths[safe] <- as.integer(x$length[safe])
  data.table::data.table(chrom = as.character(x$chrom), length = lengths)
}

bw_bigwig_query <- function(file, chrom, start, end) {
  file <- bw_validate_local_file(file)
  if (!is.character(chrom) || length(chrom) != 1L || is.na(chrom) || !nzchar(chrom)) {
    bw_stop("`chrom` must be a single non-missing character string.")
  }
  start <- suppressWarnings(as.integer(start)[1L])
  end <- suppressWarnings(as.integer(end)[1L])
  if (is.na(start) || start < 1L) bw_stop("`start` must be >= 1 in 1-based closed coordinates.")
  if (is.na(end) || end < start) bw_stop("`end` must be >= `start`.")

  metadata <- bw_metadata(file)
  chrom_idx <- match(chrom, metadata$chromosomes$chrom)
  if (is.na(chrom_idx)) {
    bw_stop(paste0("Chromosome `", chrom, "` was not found in the BigWig file."))
  }
  chrom_length <- metadata$chromosomes$length[chrom_idx]
  query_start <- as.numeric(start - 1L)
  query_end <- min(as.numeric(end), chrom_length)
  if (query_start >= chrom_length || query_start >= query_end) {
    return(bw_empty_signal())
  }

  con <- base::file(file, open = "rb")
  on.exit(close(con), add = TRUE)
  index <- bw_read_index_header(con, metadata$header)
  blocks <- bw_collect_blocks(
    con = con,
    node_offset = index$root_offset,
    endian = metadata$header$endian,
    tid = metadata$chromosomes$tid[chrom_idx],
    query_start = query_start,
    query_end = query_end
  )
  if (nrow(blocks) == 0L) {
    return(bw_empty_signal())
  }

  start_chunks <- vector("list", nrow(blocks))
  end_chunks <- vector("list", nrow(blocks))
  value_chunks <- vector("list", nrow(blocks))
  chunk_n <- 0L
  interval_n <- 0
  needs_sort <- FALSE
  previous_start <- -Inf
  previous_end <- -Inf

  for (i in seq_len(nrow(blocks))) {
    block <- bw_read_exact(con, blocks$offset[i], blocks$size[i])
    if (metadata$header$uncompress_buf_size > 0) {
      block <- bw_decompress_block(block)
    }
    decoded <- bw_decode_block_vectors(
      block,
      metadata$header$endian,
      metadata$chromosomes$tid[chrom_idx],
      query_start,
      query_end
    )
    n_decoded <- length(decoded$start0)
    if (n_decoded < 1L) next

    start1 <- decoded$start0 + 1
    end1 <- decoded$end0
    if (any(start1 > .Machine$integer.max) || any(end1 > .Machine$integer.max)) {
      bw_stop("Queried coordinates exceed the supported integer coordinate range.")
    }

    if (n_decoded > 1L) {
      start_diff <- diff(start1)
      end_diff <- diff(end1)
      if (any(start_diff < 0 | (start_diff == 0 & end_diff < 0))) {
        needs_sort <- TRUE
      }
    }
    if (
      start1[1L] < previous_start ||
      (start1[1L] == previous_start && end1[1L] < previous_end)
    ) {
      needs_sort <- TRUE
    }
    previous_start <- start1[n_decoded]
    previous_end <- end1[n_decoded]

    interval_n <- interval_n + as.numeric(n_decoded)
    if (interval_n > .Machine$integer.max) {
      bw_stop("Queried BigWig interval count exceeds the supported R vector length.")
    }

    chunk_n <- chunk_n + 1L
    start_chunks[[chunk_n]] <- as.integer(start1)
    end_chunks[[chunk_n]] <- as.integer(end1)
    value_chunks[[chunk_n]] <- decoded$value
  }

  if (chunk_n < 1L) {
    return(bw_empty_signal())
  }

  start_value <- unlist(start_chunks, use.names = FALSE)
  end_value <- unlist(end_chunks, use.names = FALSE)
  signal_value <- unlist(value_chunks, use.names = FALSE)
  interval_n <- length(start_value)

  ans <- data.table::data.table(
    chrom = rep(chrom, interval_n),
    start = start_value,
    end = end_value,
    value = signal_value
  )
  if (isTRUE(needs_sort)) {
    data.table::setorderv(ans, c("start", "end"))
  }
  ans[]
}

