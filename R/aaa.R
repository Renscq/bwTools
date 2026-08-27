# Author: Rensc
# Date: 2026-08-27
# Version: dev001
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

.bwtools_cache <- new.env(parent = emptyenv())
