# Author: Rensc
# Date: 2026-08-27
# Version: dev003
# Function: Write local BigWig files using native R binary I/O
# Input: Canonical 1-based closed signal intervals and chromosome sizes
# Output: Local BigWig files without external executables or compiled code

bw_pack_uint <- function(value, size) {
  value <- as.numeric(value)[1L]
  size <- as.integer(size)[1L]
  max_value <- 256^size - 1
  if (!is.finite(value) || value < 0 || value != floor(value) || value > max_value || value > .BWT_MAX_EXACT_UINT64) {
    bw_stop("Unsigned integer value cannot be represented exactly in the requested binary width.")
  }
  as.raw(floor(value / 256^(seq.int(0L, size - 1L))) %% 256)
}

bw_pack_u16 <- function(value) {
  bw_pack_uint(value, 2L)
}

bw_pack_u32 <- function(value) {
  bw_pack_uint(value, 4L)
}

bw_pack_u64 <- function(value) {
  bw_pack_uint(value, 8L)
}

bw_pack_u32_vector <- function(values) {
  values <- as.numeric(values)
  if (length(values) == 0L) {
    return(raw())
  }
  if (any(!is.finite(values)) || any(values < 0) || any(values != floor(values)) || any(values > 2^32 - 1)) {
    bw_stop("One or more values cannot be represented as unsigned 32-bit integers.")
  }
  bytes <- rbind(
    values %% 256,
    floor(values / 256) %% 256,
    floor(values / 256^2) %% 256,
    floor(values / 256^3) %% 256
  )
  as.raw(as.vector(bytes))
}

bw_pack_float32_vector <- function(values) {
  values <- as.numeric(values)
  if (any(!is.finite(values)) || any(abs(values) > .BWT_FLOAT32_MAX)) {
    bw_stop("BigWig values must be finite and representable as 32-bit floating-point values.")
  }
  writeBin(values, raw(), size = 4L, endian = "little")
}

bw_pack_double <- function(value) {
  value <- as.numeric(value)[1L]
  if (!is.finite(value)) {
    bw_stop("BigWig summary values must be finite.")
  }
  writeBin(value, raw(), size = 8L, endian = "little")
}

bw_pack_fixed_string <- function(value, key_size) {
  raw_value <- charToRaw(enc2utf8(as.character(value)[1L]))
  key_size <- as.integer(key_size)[1L]
  if (length(raw_value) > key_size) {
    bw_stop("Chromosome name exceeds the configured BigWig chromosome key size.")
  }
  c(raw_value, raw(key_size - length(raw_value)))
}

bw_write_raw <- function(con, x) {
  if (!is.raw(x)) {
    x <- as.raw(x)
  }
  writeBin(x, con, size = 1L)
  invisible(length(x))
}

bw_tell <- function(con) {
  as.numeric(seek(con, where = NA, origin = "start", rw = "write"))
}

bw_build_chrom_tree <- function(chrom_sizes, absolute_offset, block_size = .BWT_CHROM_TREE_BLOCK_SIZE) {
  chrom_sizes <- data.table::copy(data.table::as.data.table(chrom_sizes))
  n <- nrow(chrom_sizes)
  if (n < 1L) {
    bw_stop("At least one chromosome is required for BigWig output.")
  }
  key_size <- max(nchar(chrom_sizes$chrom, type = "bytes"))
  if (!is.finite(key_size) || key_size < 1L) {
    bw_stop("Chromosome names must be non-empty strings.")
  }

  entries <- lapply(seq_len(n), function(i) {
    list(
      key = chrom_sizes$chrom[i],
      tid = as.numeric(chrom_sizes$tid[i]),
      length = as.numeric(chrom_sizes$length[i])
    )
  })
  ord <- order(vapply(entries, function(x) x$key, character(1L)), method = "radix")
  entries <- entries[ord]

  make_leaf <- function(x) {
    list(is_leaf = TRUE, entries = x, children = NULL, first_key = x[[1L]]$key, offset = NA_real_)
  }
  leaves <- split(entries, ceiling(seq_along(entries) / block_size))
  level <- lapply(leaves, make_leaf)

  while (length(level) > 1L) {
    groups <- split(level, ceiling(seq_along(level) / block_size))
    level <- lapply(groups, function(children) {
      list(
        is_leaf = FALSE,
        entries = NULL,
        children = children,
        first_key = children[[1L]]$first_key,
        offset = NA_real_
      )
    })
  }
  root <- level[[1L]]

  node_size <- function(node) {
    n_children <- if (node$is_leaf) length(node$entries) else length(node$children)
    4 + n_children * (key_size + 8)
  }

  assign_offsets <- function(node, offset) {
    node$offset <- offset
    next_offset <- offset + node_size(node)
    if (!node$is_leaf) {
      for (i in seq_along(node$children)) {
        assigned <- assign_offsets(node$children[[i]], next_offset)
        node$children[[i]] <- assigned$node
        next_offset <- assigned$next_offset
      }
    }
    list(node = node, next_offset = next_offset)
  }

  header_size <- 32
  assigned <- assign_offsets(root, absolute_offset + header_size)
  root <- assigned$node

  serialize_node <- function(node) {
    n_children <- if (node$is_leaf) length(node$entries) else length(node$children)
    out <- c(as.raw(if (node$is_leaf) 1L else 0L), as.raw(0L), bw_pack_u16(n_children))
    if (node$is_leaf) {
      for (entry in node$entries) {
        out <- c(
          out,
          bw_pack_fixed_string(entry$key, key_size),
          bw_pack_u32(entry$tid),
          bw_pack_u32(entry$length)
        )
      }
    } else {
      for (child in node$children) {
        out <- c(
          out,
          bw_pack_fixed_string(child$first_key, key_size),
          bw_pack_u64(child$offset)
        )
      }
      for (child in node$children) {
        out <- c(out, serialize_node(child))
      }
    }
    out
  }

  body <- serialize_node(root)
  header <- c(
    bw_pack_u32(.BWT_CIRTREE_MAGIC),
    bw_pack_u32(block_size),
    bw_pack_u32(key_size),
    bw_pack_u32(8L),
    bw_pack_u64(n),
    bw_pack_u64(0L)
  )
  raw <- c(header, body)
  expected <- assigned$next_offset - absolute_offset
  if (length(raw) != expected) {
    bw_stop("Internal chromosome B+ tree size mismatch during BigWig writing.")
  }
  raw
}

bw_encode_bedgraph_block <- function(tid, start0, end0, value) {
  n <- length(start0)
  if (n < 1L || length(end0) != n || length(value) != n) {
    bw_stop("Invalid BigWig data block input.")
  }
  if (n > 65535L) {
    bw_stop("A BigWig data block cannot contain more than 65535 entries.")
  }
  header <- c(
    bw_pack_u32(tid),
    bw_pack_u32(start0[1L]),
    bw_pack_u32(end0[n]),
    bw_pack_u32(0L),
    bw_pack_u32(0L),
    as.raw(1L),
    as.raw(0L),
    bw_pack_u16(n)
  )
  start_raw <- matrix(bw_pack_u32_vector(start0), nrow = 4L)
  end_raw <- matrix(bw_pack_u32_vector(end0), nrow = 4L)
  value_raw <- matrix(bw_pack_float32_vector(value), nrow = 4L)
  payload <- as.raw(as.vector(rbind(start_raw, end_raw, value_raw)))
  c(header, payload)
}

bw_build_rtree <- function(
  index_entries,
  absolute_offset,
  block_size = .BWT_RTREE_BLOCK_SIZE,
  items_per_slot = 1L
) {
  index_entries <- data.table::copy(data.table::as.data.table(index_entries))
  if (nrow(index_entries) < 1L) {
    bw_stop("Cannot build a BigWig R-tree without data blocks.")
  }

  make_leaf <- function(idx) {
    list(
      is_leaf = TRUE,
      entries = index_entries[idx],
      children = NULL,
      offset = NA_real_,
      bounds = list(
        chrom_start = index_entries$chrom_start[idx[1L]],
        base_start = index_entries$base_start[idx[1L]],
        chrom_end = index_entries$chrom_end[idx[length(idx)]],
        base_end = index_entries$base_end[idx[length(idx)]]
      )
    )
  }

  leaf_groups <- split(seq_len(nrow(index_entries)), ceiling(seq_len(nrow(index_entries)) / block_size))
  level <- lapply(leaf_groups, make_leaf)

  make_parent <- function(children) {
    first <- children[[1L]]$bounds
    last <- children[[length(children)]]$bounds
    list(
      is_leaf = FALSE,
      entries = NULL,
      children = children,
      offset = NA_real_,
      bounds = list(
        chrom_start = first$chrom_start,
        base_start = first$base_start,
        chrom_end = last$chrom_end,
        base_end = last$base_end
      )
    )
  }

  while (length(level) > 1L) {
    groups <- split(level, ceiling(seq_along(level) / block_size))
    level <- lapply(groups, make_parent)
  }
  root <- level[[1L]]

  node_size <- function(node) {
    n_children <- if (node$is_leaf) nrow(node$entries) else length(node$children)
    4 + n_children * if (node$is_leaf) 32 else 24
  }

  assign_offsets <- function(node, offset) {
    node$offset <- offset
    next_offset <- offset + node_size(node)
    if (!node$is_leaf) {
      for (i in seq_along(node$children)) {
        assigned <- assign_offsets(node$children[[i]], next_offset)
        node$children[[i]] <- assigned$node
        next_offset <- assigned$next_offset
      }
    }
    list(node = node, next_offset = next_offset)
  }

  root_offset <- absolute_offset + .BWT_INDEX_HEADER_SIZE
  assigned <- assign_offsets(root, root_offset)
  root <- assigned$node

  serialize_node <- function(node) {
    n_children <- if (node$is_leaf) nrow(node$entries) else length(node$children)
    out <- c(as.raw(if (node$is_leaf) 1L else 0L), as.raw(0L), bw_pack_u16(n_children))
    if (node$is_leaf) {
      for (i in seq_len(nrow(node$entries))) {
        x <- node$entries[i]
        out <- c(
          out,
          bw_pack_u32(x$chrom_start),
          bw_pack_u32(x$base_start),
          bw_pack_u32(x$chrom_end),
          bw_pack_u32(x$base_end),
          bw_pack_u64(x$data_offset),
          bw_pack_u64(x$data_size)
        )
      }
    } else {
      for (child in node$children) {
        b <- child$bounds
        out <- c(
          out,
          bw_pack_u32(b$chrom_start),
          bw_pack_u32(b$base_start),
          bw_pack_u32(b$chrom_end),
          bw_pack_u32(b$base_end),
          bw_pack_u64(child$offset)
        )
      }
      for (child in node$children) {
        out <- c(out, serialize_node(child))
      }
    }
    out
  }

  body <- serialize_node(root)
  actual_tree_size <- assigned$next_offset - root_offset
  if (length(body) != actual_tree_size) {
    bw_stop("Internal R-tree size mismatch during BigWig writing.")
  }
  index_size <- if (root$is_leaf) {
    4 + 24 * nrow(root$entries)
  } else {
    actual_tree_size
  }
  b <- root$bounds
  header <- c(
    bw_pack_u32(.BWT_INDEX_MAGIC),
    bw_pack_u32(block_size),
    bw_pack_u64(nrow(index_entries)),
    bw_pack_u32(b$chrom_start),
    bw_pack_u32(b$base_start),
    bw_pack_u32(b$chrom_end),
    bw_pack_u32(b$base_end),
    bw_pack_u64(index_size),
    bw_pack_u32(items_per_slot),
    bw_pack_u32(0L)
  )
  c(header, body)
}

bw_bigwig_summary <- function(signal) {
  width <- as.numeric(signal$end) - as.numeric(signal$start) + 1
  value <- as.numeric(signal$value)
  list(
    n_bases_covered = sum(width),
    min_value = min(value),
    max_value = max(value),
    sum_data = sum(width * value),
    sum_squared = sum(width * value^2)
  )
}

bw_build_header <- function(
  chrom_tree_offset,
  data_offset,
  index_offset,
  summary_offset,
  uncompress_buf_size,
  n_levels = 0L
) {
  n_levels <- as.integer(n_levels)[1L]
  if (is.na(n_levels) || n_levels < 0L || n_levels > 65535L) {
    bw_stop("`n_levels` must be between 0 and 65535.")
  }
  c(
    bw_pack_u32(.BWT_BIGWIG_MAGIC),
    bw_pack_u16(4L),
    bw_pack_u16(n_levels),
    bw_pack_u64(chrom_tree_offset),
    bw_pack_u64(data_offset),
    bw_pack_u64(index_offset),
    bw_pack_u16(0L),
    bw_pack_u16(0L),
    bw_pack_u64(0L),
    bw_pack_u64(summary_offset),
    bw_pack_u32(uncompress_buf_size),
    bw_pack_u64(0L)
  )
}

bw_build_zoom_headers_raw <- function(zoom_headers) {
  zoom_headers <- data.table::as.data.table(zoom_headers)
  if (nrow(zoom_headers) < 1L) {
    return(raw())
  }
  out <- vector("list", nrow(zoom_headers))
  for (i in seq_len(nrow(zoom_headers))) {
    out[[i]] <- c(
      bw_pack_u32(zoom_headers$level[i]),
      bw_pack_u32(0L),
      bw_pack_u64(zoom_headers$data_offset[i]),
      bw_pack_u64(zoom_headers$index_offset[i])
    )
  }
  do.call(c, out)
}

bw_build_summary_raw <- function(summary) {
  c(
    bw_pack_u64(summary$n_bases_covered),
    bw_pack_double(summary$min_value),
    bw_pack_double(summary$max_value),
    bw_pack_double(summary$sum_data),
    bw_pack_double(summary$sum_squared)
  )
}

bw_select_zoom_levels <- function(signal, chrom_sizes, max_zoom_levels = .BWT_DEFAULT_MAX_ZOOM_LEVELS) {
  max_zoom_levels <- suppressWarnings(as.integer(max_zoom_levels)[1L])
  if (is.na(max_zoom_levels) || max_zoom_levels < 0L || max_zoom_levels > 65535L) {
    bw_stop("`max_zoom_levels` must be an integer between 0 and 65535.")
  }
  if (max_zoom_levels == 0L || nrow(signal) < 1L) {
    return(integer())
  }
  max_chrom <- max(as.numeric(chrom_sizes$length))
  if (!is.finite(max_chrom) || max_chrom < 1) {
    return(integer())
  }
  mean_width <- floor(mean(as.numeric(signal$end) - as.numeric(signal$start) + 1))
  first_level <- max(10, mean_width * 16)
  first_level <- min(first_level, max_chrom)
  if (first_level < 1 || first_level > 2^32 - 1) {
    return(integer())
  }
  levels <- numeric()
  current <- as.numeric(first_level)
  while (length(levels) < max_zoom_levels && current <= max_chrom && current <= 2^32 - 1) {
    levels <- c(levels, current)
    if (current > (2^32 - 1) / 4) break
    next_level <- current * 4
    if (next_level <= current) break
    current <- next_level
  }
  as.numeric(unique(levels))
}

bw_zoom_chrom_lengths <- function(chrom_sizes) {
  chrom_sizes <- data.table::as.data.table(chrom_sizes)
  stats::setNames(
    as.numeric(chrom_sizes$length),
    as.character(chrom_sizes$tid)
  )
}

bw_finish_zoom_records <- function(records, chrom_sizes, level) {
  records <- data.table::as.data.table(records)
  if (nrow(records) < 1L) {
    return(data.table::data.table())
  }

  chrom_length <- bw_zoom_chrom_lengths(chrom_sizes)
  end0 <- pmin(
    as.numeric(records$start0) + as.numeric(level),
    unname(chrom_length[as.character(records$tid)])
  )
  if (any(!is.finite(end0))) {
    bw_stop("Failed to resolve chromosome length while constructing BigWig zoom levels.")
  }
  data.table::set(records, j = "end0", value = as.numeric(end0))
  data.table::setcolorder(
    records,
    c(
      "tid", "start0", "end0", "valid_count",
      "min_value", "max_value", "sum_data", "sum_squared"
    )
  )
  data.table::setorderv(records, c("tid", "start0", "end0"))
  records[]
}

bw_split_crossing_zoom_intervals <- function(
  tid,
  start0,
  end0,
  value,
  chrom_length,
  level
) {
  n <- length(start0)
  if (n < 1L) {
    return(data.table::data.table())
  }

  tid_parts <- vector("list", n)
  start_parts <- vector("list", n)
  width_parts <- vector("list", n)
  value_parts <- vector("list", n)

  for (i in seq_len(n)) {
    pos <- as.numeric(start0[i])
    interval_end <- as.numeric(end0[i])
    chr_len <- as.numeric(chrom_length[i])
    if (!is.finite(chr_len)) {
      bw_stop("Failed to resolve chromosome length while constructing BigWig zoom levels.")
    }

    first_bin <- floor(pos / level) * level
    last_bin <- floor((interval_end - 1) / level) * level
    bin_starts <- seq(from = first_bin, to = last_bin, by = level)
    bin_ends <- pmin(bin_starts + level, chr_len)
    overlap_start <- pmax(pos, bin_starts)
    overlap_end <- pmin(interval_end, bin_ends)
    width <- overlap_end - overlap_start
    keep <- width > 0
    if (!any(keep)) next

    n_keep <- sum(keep)
    tid_parts[[i]] <- rep(as.numeric(tid[i]), n_keep)
    start_parts[[i]] <- as.numeric(bin_starts[keep])
    width_parts[[i]] <- as.numeric(width[keep])
    value_parts[[i]] <- rep(as.numeric(value[i]), n_keep)
  }

  keep_parts <- lengths(start_parts) > 0L
  if (!any(keep_parts)) {
    return(data.table::data.table())
  }
  data.table::data.table(
    tid = unlist(tid_parts[keep_parts], use.names = FALSE),
    start0 = unlist(start_parts[keep_parts], use.names = FALSE),
    width = unlist(width_parts[keep_parts], use.names = FALSE),
    value = unlist(value_parts[keep_parts], use.names = FALSE)
  )
}

bw_build_zoom_records <- function(signal, chrom_sizes, level) {
  level <- as.numeric(level)[1L]
  if (!is.finite(level) || level < 1 || level != floor(level)) {
    bw_stop("Zoom reduction level must be a positive integer.")
  }

  signal <- data.table::as.data.table(signal)
  chrom_sizes <- data.table::as.data.table(chrom_sizes)
  if (nrow(signal) < 1L) {
    return(data.table::data.table())
  }

  tid <- as.numeric(signal$tid)
  start0 <- as.numeric(signal$start) - 1
  end0 <- as.numeric(signal$end)
  value <- as.numeric(signal$value)
  width <- end0 - start0
  start_bin <- floor(start0 / level) * level
  end_bin <- floor((end0 - 1) / level) * level
  same_bin <- start_bin == end_bin

  pieces <- vector("list", 2L)
  if (any(same_bin)) {
    pieces[[1L]] <- data.table::data.table(
      tid = tid[same_bin],
      start0 = start_bin[same_bin],
      width = width[same_bin],
      value = value[same_bin]
    )
  } else {
    pieces[[1L]] <- data.table::data.table()
  }

  if (any(!same_bin)) {
    chrom_length <- bw_zoom_chrom_lengths(chrom_sizes)
    pieces[[2L]] <- bw_split_crossing_zoom_intervals(
      tid = tid[!same_bin],
      start0 = start0[!same_bin],
      end0 = end0[!same_bin],
      value = value[!same_bin],
      chrom_length = unname(chrom_length[as.character(tid[!same_bin])]),
      level = level
    )
  } else {
    pieces[[2L]] <- data.table::data.table()
  }

  pieces <- pieces[vapply(pieces, nrow, integer(1L)) > 0L]
  if (length(pieces) < 1L) {
    return(data.table::data.table())
  }
  x <- data.table::rbindlist(pieces, use.names = TRUE)
  data.table::set(
    x,
    j = "sum_data",
    value = as.numeric(x$width) * as.numeric(x$value)
  )
  data.table::set(
    x,
    j = "sum_squared",
    value = as.numeric(x$width) * as.numeric(x$value)^2
  )

  ans <- x[
    ,
    .(
      valid_count = sum(as.numeric(width)),
      min_value = min(as.numeric(value)),
      max_value = max(as.numeric(value)),
      sum_data = sum(as.numeric(sum_data)),
      sum_squared = sum(as.numeric(sum_squared))
    ),
    by = .(tid, start0)
  ]
  bw_finish_zoom_records(ans, chrom_sizes, level)
}

bw_aggregate_zoom_records <- function(records, chrom_sizes, level) {
  level <- as.numeric(level)[1L]
  if (!is.finite(level) || level < 1 || level != floor(level)) {
    bw_stop("Zoom reduction level must be a positive integer.")
  }

  records <- data.table::as.data.table(records)
  chrom_sizes <- data.table::as.data.table(chrom_sizes)
  if (nrow(records) < 1L) {
    return(data.table::data.table())
  }

  parent_start <- floor(as.numeric(records$start0) / level) * level
  chrom_length <- bw_zoom_chrom_lengths(chrom_sizes)
  parent_end <- pmin(
    parent_start + level,
    unname(chrom_length[as.character(records$tid)])
  )
  if (any(!is.finite(parent_end))) {
    bw_stop("Failed to resolve chromosome length while aggregating BigWig zoom levels.")
  }
  if (any(as.numeric(records$end0) > parent_end)) {
    bw_stop(
      "Previous BigWig zoom records cross the requested parent zoom window."
    )
  }

  x <- data.table::data.table(
    tid = as.numeric(records$tid),
    start0 = as.numeric(parent_start),
    valid_count = as.numeric(records$valid_count),
    min_value = as.numeric(records$min_value),
    max_value = as.numeric(records$max_value),
    sum_data = as.numeric(records$sum_data),
    sum_squared = as.numeric(records$sum_squared)
  )
  ans <- x[
    ,
    .(
      valid_count = sum(valid_count),
      min_value = min(min_value),
      max_value = max(max_value),
      sum_data = sum(sum_data),
      sum_squared = sum(sum_squared)
    ),
    by = .(tid, start0)
  ]
  bw_finish_zoom_records(ans, chrom_sizes, level)
}

bw_build_zoom_pyramid <- function(signal, chrom_sizes, levels) {
  levels <- as.numeric(levels)
  if (length(levels) < 1L) {
    return(list())
  }
  if (any(!is.finite(levels)) || any(levels < 1) || any(levels != floor(levels))) {
    bw_stop("Zoom reduction levels must be positive integers.")
  }
  if (length(levels) > 1L && any(levels[-1L] <= levels[-length(levels)])) {
    bw_stop("Zoom reduction levels must be strictly increasing.")
  }

  out <- vector("list", length(levels))
  out[[1L]] <- bw_build_zoom_records(signal, chrom_sizes, levels[1L])
  if (length(levels) > 1L) {
    for (i in 2:length(levels)) {
      out[[i]] <- bw_aggregate_zoom_records(
        out[[i - 1L]],
        chrom_sizes,
        levels[i]
      )
    }
  }
  out
}

bw_encode_zoom_block <- function(records) {
  records <- data.table::as.data.table(records)
  n <- nrow(records)
  if (n < 1L) return(raw())
  tid_raw <- matrix(bw_pack_u32_vector(records$tid), nrow = 4L)
  start_raw <- matrix(bw_pack_u32_vector(records$start0), nrow = 4L)
  end_raw <- matrix(bw_pack_u32_vector(records$end0), nrow = 4L)
  count_raw <- matrix(bw_pack_u32_vector(records$valid_count), nrow = 4L)
  min_raw <- matrix(bw_pack_float32_vector(records$min_value), nrow = 4L)
  max_raw <- matrix(bw_pack_float32_vector(records$max_value), nrow = 4L)
  sum_raw <- matrix(bw_pack_float32_vector(records$sum_data), nrow = 4L)
  sumsq_raw <- matrix(bw_pack_float32_vector(records$sum_squared), nrow = 4L)
  as.raw(as.vector(rbind(
    tid_raw, start_raw, end_raw, count_raw,
    min_raw, max_raw, sum_raw, sumsq_raw
  )))
}

bw_write_zoom_level <- function(con, records, level) {
  records <- data.table::copy(data.table::as.data.table(records))
  if (nrow(records) < 1L) {
    bw_stop("Cannot write an empty BigWig zoom level.")
  }
  block_limit <- floor(.BWT_WRITE_BUFFER_SIZE / .BWT_ZOOM_RECORD_SIZE)
  chunks <- split(seq_len(nrow(records)), ceiling(seq_len(nrow(records)) / block_limit))
  data_offset <- bw_tell(con)
  bw_write_raw(con, bw_pack_u32(length(chunks)))

  index_parts <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    x <- records[chunks[[i]]]
    raw_block <- bw_encode_zoom_block(x)
    if (length(raw_block) > .BWT_WRITE_BUFFER_SIZE) {
      bw_stop("Internal BigWig zoom block exceeds the configured uncompressed buffer size.")
    }
    compressed <- base::memCompress(raw_block, type = "gzip")
    block_offset <- bw_tell(con)
    bw_write_raw(con, compressed)
    index_parts[[i]] <- data.table::data.table(
      chrom_start = x$tid[1L],
      base_start = x$start0[1L],
      chrom_end = x$tid[nrow(x)],
      base_end = x$end0[nrow(x)],
      data_offset = block_offset,
      data_size = length(compressed)
    )
  }
  index_entries <- data.table::rbindlist(index_parts, use.names = TRUE)
  index_offset <- bw_tell(con)
  bw_write_raw(con, bw_build_rtree(
    index_entries,
    absolute_offset = index_offset,
    items_per_slot = block_limit
  ))
  list(
    level = as.numeric(level),
    data_offset = as.numeric(data_offset),
    index_offset = as.numeric(index_offset)
  )
}

bw_write_bigwig_file <- function(
  signal,
  chrom_sizes,
  file,
  zoom = TRUE,
  max_zoom_levels = .BWT_DEFAULT_MAX_ZOOM_LEVELS
) {
  signal <- data.table::copy(data.table::as.data.table(signal))
  chrom_sizes <- data.table::copy(data.table::as.data.table(chrom_sizes))
  if (nrow(signal) < 1L) {
    bw_stop("BigWig output requires at least one signal interval.")
  }

  zoom_levels <- if (isTRUE(zoom)) {
    bw_select_zoom_levels(signal, chrom_sizes, max_zoom_levels = max_zoom_levels)
  } else {
    integer()
  }
  n_levels <- length(zoom_levels)
  summary <- bw_bigwig_summary(signal)
  summary_offset <- .BWT_HEADER_SIZE + n_levels * 24L
  chrom_tree_offset <- summary_offset + .BWT_SUMMARY_SIZE
  chrom_tree_raw <- bw_build_chrom_tree(chrom_sizes, absolute_offset = chrom_tree_offset)
  data_offset <- chrom_tree_offset + length(chrom_tree_raw)

  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = paste0(".", basename(file), "."), tmpdir = dirname(file), fileext = ".tmp")
  con <- base::file(tmp, open = "w+b")
  ok <- FALSE
  on.exit({
    if (!is.null(con)) try(close(con), silent = TRUE)
    if (!ok && file.exists(tmp)) unlink(tmp)
  }, add = TRUE)

  placeholder_size <- .BWT_HEADER_SIZE + n_levels * 24L + .BWT_SUMMARY_SIZE
  bw_write_raw(con, raw(placeholder_size))
  bw_write_raw(con, chrom_tree_raw)
  if (bw_tell(con) != data_offset) {
    bw_stop("Internal BigWig data offset mismatch before writing signal blocks.")
  }

  block_limit <- floor((.BWT_WRITE_BUFFER_SIZE - .BWT_DATA_HEADER_SIZE) / 12L)
  chrom_groups <- split(seq_len(nrow(signal)), signal$tid, drop = TRUE)
  n_blocks <- sum(vapply(chrom_groups, function(idx) ceiling(length(idx) / block_limit), numeric(1L)))
  bw_write_raw(con, bw_pack_u64(n_blocks))

  index_parts <- vector("list", n_blocks)
  block_i <- 0L
  for (idx in chrom_groups) {
    chunks <- split(idx, ceiling(seq_along(idx) / block_limit))
    for (chunk_idx in chunks) {
      x <- signal[chunk_idx]
      block_raw <- bw_encode_bedgraph_block(
        tid = x$tid[1L],
        start0 = as.numeric(x$start) - 1,
        end0 = as.numeric(x$end),
        value = x$value
      )
      if (length(block_raw) > .BWT_WRITE_BUFFER_SIZE) {
        bw_stop("Internal BigWig data block exceeds the configured uncompressed buffer size.")
      }
      compressed <- base::memCompress(block_raw, type = "gzip")
      block_offset <- bw_tell(con)
      bw_write_raw(con, compressed)
      block_i <- block_i + 1L
      index_parts[[block_i]] <- data.table::data.table(
        chrom_start = as.numeric(x$tid[1L]),
        base_start = as.numeric(x$start[1L]) - 1,
        chrom_end = as.numeric(x$tid[1L]),
        base_end = as.numeric(x$end[nrow(x)]),
        data_offset = block_offset,
        data_size = length(compressed)
      )
    }
  }
  if (block_i != n_blocks) {
    bw_stop("Internal BigWig block count mismatch.")
  }

  index_entries <- data.table::rbindlist(index_parts, use.names = TRUE)
  index_offset <- bw_tell(con)
  bw_write_raw(con, bw_build_rtree(index_entries, absolute_offset = index_offset))

  zoom_headers <- vector("list", n_levels)
  if (n_levels > 0L) {
    previous_records <- NULL
    for (i in seq_along(zoom_levels)) {
      records <- if (i == 1L) {
        bw_build_zoom_records(signal, chrom_sizes, zoom_levels[i])
      } else {
        bw_aggregate_zoom_records(
          previous_records,
          chrom_sizes,
          zoom_levels[i]
        )
      }
      zoom_headers[[i]] <- bw_write_zoom_level(con, records, zoom_levels[i])
      previous_records <- records
      rm(records)
    }
    rm(previous_records)
  }
  zoom_headers_dt <- if (n_levels > 0L) {
    data.table::rbindlist(zoom_headers, use.names = TRUE)
  } else {
    data.table::data.table(level = integer(), data_offset = numeric(), index_offset = numeric())
  }

  bw_write_raw(con, bw_pack_u32(.BWT_BIGWIG_MAGIC))

  header_raw <- bw_build_header(
    chrom_tree_offset = chrom_tree_offset,
    data_offset = data_offset,
    index_offset = index_offset,
    summary_offset = summary_offset,
    uncompress_buf_size = .BWT_WRITE_BUFFER_SIZE,
    n_levels = n_levels
  )
  zoom_header_raw <- bw_build_zoom_headers_raw(zoom_headers_dt)
  summary_raw <- bw_build_summary_raw(summary)
  if (length(header_raw) != .BWT_HEADER_SIZE ||
      length(zoom_header_raw) != n_levels * 24L ||
      length(summary_raw) != .BWT_SUMMARY_SIZE) {
    bw_stop("Internal BigWig header, zoom header, or summary size mismatch.")
  }
  seek(con, where = 0, origin = "start", rw = "write")
  bw_write_raw(con, header_raw)
  if (n_levels > 0L) {
    seek(con, where = .BWT_HEADER_SIZE, origin = "start", rw = "write")
    bw_write_raw(con, zoom_header_raw)
  }
  seek(con, where = summary_offset, origin = "start", rw = "write")
  bw_write_raw(con, summary_raw)
  close(con)
  con <- NULL

  if (file.exists(file)) {
    unlink(file)
  }
  if (!file.rename(tmp, file)) {
    copied <- file.copy(tmp, file, overwrite = TRUE)
    unlink(tmp)
    if (!isTRUE(copied)) {
      bw_stop(paste0("Failed to move completed BigWig file into place: ", file))
    }
  }
  ok <- TRUE
  normalizePath(file, winslash = "/", mustWork = TRUE)
}

