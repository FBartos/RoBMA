test_that("exact recovery excludes the integer Ledermann boundary", {

  expect_identical(vapply(c(3L, 6L, 10L), .covariance_ledermann_rank, integer(1)),
    c(0L, 2L, 5L))
  small <- diag(3) + tcrossprod(c(.1, .2, .3))
  expect_null(.covariance_minimum_rank_factor(small))
  loading <- matrix(c(.1, .2, -.3, .4, .15, -.2,
    .3, -.2, .1, -.15, .4, .25, -.2, .1, .3, .2, -.1, .4), 6L, 3L)
  expect_null(.covariance_minimum_rank_factor(tcrossprod(loading) + diag(.2, 6L)))
  signed <- matrix(c(.04, -.012, -.012, .05), 2L)
  recovered <- .covariance_minimum_rank_factor(signed)
  expect_identical(recovered$rank, 1L)
  expect_lt(max(abs(tcrossprod(recovered$loading) + diag(recovered$diagonal) - signed)), 1e-15)
})
