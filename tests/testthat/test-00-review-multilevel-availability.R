test_that("multilevel sufficient-statistic underflow is unavailable rather than a zero effect", {

  expect_error(.evaluate.brma.multilevel_blup.norm(matrix(0, 1L, 1L),
    matrix(0, 1L, 1L), matrix(2^511, 1L, 1L), 2^-563, 2^1022, 1L),
    "Multilevel Gaussian conditional arithmetic is numerically unavailable.", fixed = TRUE)
  zero <- .evaluate.brma.multilevel_blup.norm(matrix(0, 1L, 1L),
    matrix(0, 1L, 1L), matrix(2^511, 1L, 1L), 0, 2^1022, 1L)
  expect_equal(zero$cluster, matrix(0, 1L, 1L), tolerance = 0)
  expect_equal(zero$estimate, matrix(0, 1L, 1L), tolerance = 0)
})
