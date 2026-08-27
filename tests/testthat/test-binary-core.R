test_that("binary integer helpers decode little and big endian values", {
  x <- as.raw(c(0x26, 0xFC, 0x8F, 0x88))
  expect_equal(bwTools:::bw_uint(x, 0L, 4L, "little"), 0x888FFC26)
  expect_equal(bwTools:::bw_uint(rev(x), 0L, 4L, "big"), 0x888FFC26)
})

test_that("IEEE754 float32 decoder handles positive and negative values", {
  x <- writeBin(c(1.5, -3.25), raw(), size = 4L, endian = "little")
  observed <- bwTools:::bw_float32_vector(x, c(0, 4), "little")
  expect_equal(observed, c(1.5, -3.25), tolerance = 1e-7)
})
