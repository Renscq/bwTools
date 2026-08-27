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
