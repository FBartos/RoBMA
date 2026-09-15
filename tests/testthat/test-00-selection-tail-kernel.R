context("(00) Standard-normal tail kernel")

# The selection kernels evaluate the standard-normal upper tail millions of
# times per density line, through a batched approximation instead of the
# library's scalar erfc. These tests certify that approximation against
# `pnorm()` over a dense argument grid, including the far tail where the
# library routine itself loses relative accuracy, and assert that the batched
# result does not depend on how the arguments are grouped.

test_that("the tail kernel agrees with pnorm across the supported range", {

  grid <- sort(unique(c(
    seq(-40, 40, by = 0.001),
    seq(-1, 1, by = 1e-5),
    seq(0.6, 0.75, by = 1e-6),          # the central/middle boundary
    seq(5.5, 5.8, by = 1e-6),           # the middle/far boundary
    c(0, -0.5, 0.5, 8, 20, 30, 37, 37.5)
  )))
  kernel    <- .selnorm_normal_upper_tail(grid)
  reference <- stats::pnorm(grid, lower.tail = FALSE)

  expect_true(all(is.finite(kernel)))
  expect_true(all(kernel >= 0 & kernel <= 1))

  # Relative agreement wherever the reference is representable. Cody's ranges
  # are shared with pnorm, so this is double rounding, not approximation error.
  comparable <- reference > 0 & abs(grid) <= 37.5
  expect_lt(max(abs(kernel[comparable] - reference[comparable]) /
                  reference[comparable]), 1e-13)

  # Beyond the range both algorithms report an exact zero.
  expect_true(all(kernel[grid > 37.5193] == 0))
  expect_identical(.selnorm_normal_upper_tail(0), 0.5)
})


test_that("the tail kernel matches an arbitrary-precision reference", {

  # Reference values computed to 20 significant digits with an
  # arbitrary-precision implementation of the complementary error function.
  exact <- c(
    `2`  = 0.0227501319481792072,
    `5`  = 2.8665157187919391167e-7,
    `8`  = 6.2209605742717841235e-16,
    `20` = 2.7536241186062336951e-89,
    `37` = 5.7255712225245768227e-300
  )
  kernel <- .selnorm_normal_upper_tail(as.numeric(names(exact)))
  expect_lt(max(abs(kernel / exact - 1)), 1e-15)
})


test_that("a tail value does not depend on its position in the batch", {

  grid <- sort(unique(c(seq(-12, 12, by = 0.0007), 0, 0.67448975,
                        5.656854249492380195)))
  reference <- .selnorm_normal_upper_tail(grid)

  for (offset in 1:5) {
    shifted <- .selnorm_normal_upper_tail(c(rep(0, offset), grid))
    expect_identical(shifted[-seq_len(offset)], reference)
  }

  # Length-one requests take the same path as a full batch.
  probe <- sample(seq_along(grid), 40L)
  expect_identical(vapply(grid[probe], .selnorm_normal_upper_tail, numeric(1L)),
                   reference[probe])
})


test_that("the tail kernel handles non-finite and extreme arguments", {

  expect_identical(.selnorm_normal_upper_tail(Inf), 0)
  expect_identical(.selnorm_normal_upper_tail(-Inf), 1)
  expect_true(is.na(.selnorm_normal_upper_tail(NaN)))
  expect_true(is.na(.selnorm_normal_upper_tail(NA_real_)))
  expect_identical(.selnorm_normal_upper_tail(numeric(0)), numeric(0))
})
