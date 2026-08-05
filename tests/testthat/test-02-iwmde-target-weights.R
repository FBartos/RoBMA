context("IWMDE linear target weights")


test_that("linear target weights are validated before exact zeros are removed", {

  tiny <- .Machine$double.eps / 2
  expect_identical(
    .iwmde_linear_weights(c(mu = tiny, tau = 0)),
    c(mu = tiny)
  )
  expect_identical(
    .iwmde_linear_weights(matrix(
      c(1, 0),
      nrow     = 1L,
      dimnames = list(NULL, c("mu", "tau"))
    )),
    c(mu = 1)
  )

  for (invalid in list(NA_real_, NaN, Inf, -Inf)) {
    expect_error(
      .iwmde_linear_weights(c(mu = 0, tau = invalid)),
      "must all be finite"
    )
  }
  expect_error(
    .iwmde_linear_weights(stats::setNames(c(1, 2), c("mu", "mu"))),
    "duplicate parameter name"
  )
  expect_error(.iwmde_linear_weights(c(1, 2)), "fully named")
  expect_error(.iwmde_linear_weights(c(mu = "one")), "must be numeric")
  expect_error(
    .iwmde_linear_weights(matrix(
      1:4,
      nrow     = 2L,
      dimnames = list(NULL, c("mu", "tau"))
    )),
    "exactly one row"
  )
})
