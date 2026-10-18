test_that("multilevel sufficient-statistic underflow is unavailable rather than a zero effect", {

  expect_error(.evaluate.brma.multilevel_blup.norm(matrix(0, 1L, 1L),
    matrix(0, 1L, 1L), matrix(2^511, 1L, 1L), 2^-563, 2^1022, 1L),
    "Multilevel Gaussian conditional arithmetic is numerically unavailable.", fixed = TRUE)
  zero <- .evaluate.brma.multilevel_blup.norm(matrix(0, 1L, 1L),
    matrix(0, 1L, 1L), matrix(2^511, 1L, 1L), 0, 2^1022, 1L)
  expect_equal(zero$cluster, matrix(0, 1L, 1L), tolerance = 0)
  expect_equal(zero$estimate, matrix(0, 1L, 1L), tolerance = 0)
})


test_that("multilevel total conditional means use the stable weighted expression", {

  value <- .evaluate.brma.multilevel_blup.norm(matrix(0, 1L, 2L),
    matrix(1, 1L, 2L), matrix(1, 1L, 2L), c(1, -2^56), c(0, 1), c(1L, 1L))
  # The first fitted response has exactly zero sampling noise, so its total
  # latent conditional mean is y1 even when its two contributions nearly cancel.
  expect_equal(value$total[1L, 1L], 1, tolerance = 0)
  expect_equal(value$total[1L, 2L], .5 * (-2^56) + .5 * value$cluster[1L, 2L], tolerance = 0)
})
