# Author: Rensc
# Date: 2026-08-27
# Version: dev006
# Function: Define package constants and internal caches
# Input: None
# Output: Internal constants and environments

.BWT_BIGWIG_MAGIC <- 0x888FFC26
.BWT_CIRTREE_MAGIC <- 0x78CA8C91
.BWT_INDEX_MAGIC <- 0x2468ACE0
.BWT_MAX_EXACT_UINT64 <- 2^53 - 1
.BWT_HEADER_SIZE <- 64L
.BWT_INDEX_HEADER_SIZE <- 48L
.BWT_DATA_HEADER_SIZE <- 24L
.BWT_SUMMARY_SIZE <- 40L
.BWT_WRITE_BUFFER_SIZE <- 32768L
.BWT_CHROM_TREE_BLOCK_SIZE <- 256L
.BWT_RTREE_BLOCK_SIZE <- 64L
.BWT_DEFAULT_MAX_ZOOM_LEVELS <- 10L
.BWT_ZOOM_RECORD_SIZE <- 32L
.BWT_ZOOM_BLOCK_HEADER_SIZE <- 4L
.BWT_FLOAT32_MAX <- 3.402823466e38

.bwtools_cache <- new.env(parent = emptyenv())

# Tell data.table that bwTools intentionally uses data.table query semantics.
.datatable.aware <- TRUE

#' bwTools: Native R tools for genomic signal tracks
#'
#' Native R tools for reading, writing, retrieving, merging, and summarizing
#' BigWig, WIG, and bedGraph signal tracks. Public analysis functions operate on
#' the standardized BwgTrack schema v2 and use 1-based closed coordinates.
#'
#' @keywords internal
"_PACKAGE"
