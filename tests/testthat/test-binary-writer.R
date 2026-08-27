test_that("native unsigned integer packers round-trip reader helpers", {
  values32 <- c(0, 1, 255, 256, 65535, 2^31, 2^32 - 1)
  raw32 <- bwTools:::bw_pack_u32_vector(values32)
  observed32 <- bwTools:::bw_u32_vector(raw32, seq(0, by = 4, length.out = length(values32)), "little")
  expect_equal(observed32, values32)

  values64 <- c(0, 1, 2^32, 2^40)
  for (value in values64) {
    raw64 <- bwTools:::bw_pack_u64(value)
    expect_equal(bwTools:::bw_u64(raw64, 0L, "little"), value)
  }
})

test_that("native float32 writer is compatible with the reader decoder", {
  values <- c(1.5, -3.25, 0, 100.125)
  raw_values <- bwTools:::bw_pack_float32_vector(values)
  observed <- bwTools:::bw_float32_vector(
    raw_values,
    seq(0, by = 4, length.out = length(values)),
    "little"
  )
  expect_equal(observed, values, tolerance = 1e-6)
})

test_that("native zoom record encoder is compatible with the reader decoder", {
  records <- data.table::data.table(
    tid = c(0, 0),
    start0 = c(100, 200),
    end0 = c(150, 260),
    valid_count = c(50, 40),
    min_value = c(1.5, -2),
    max_value = c(3.5, 4),
    sum_data = c(100, 80),
    sum_squared = c(250, 200)
  )
  raw_records <- bwTools:::bw_encode_zoom_block(records)
  observed <- bwTools:::bw_decode_zoom_block(
    raw_records,
    endian = "little",
    tid = 0,
    query_start = 0,
    query_end = 1000,
    chrom = "chr1"
  )
  expect_equal(observed$start0, records$start0)
  expect_equal(observed$end0, records$end0)
  expect_equal(observed$valid_count, records$valid_count)
  expect_equal(observed$min_value, records$min_value, tolerance = 1e-6)
  expect_equal(observed$max_value, records$max_value, tolerance = 1e-6)
  expect_equal(observed$sum_data, records$sum_data, tolerance = 1e-6)
  expect_equal(observed$sum_squared, records$sum_squared, tolerance = 1e-6)
})
