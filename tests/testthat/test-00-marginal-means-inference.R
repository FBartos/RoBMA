context("Marginal-means inference")


.fixed_marginal_posterior <- function(value = 0, n = 20L) {

  posterior <- structure(
    rep(value, n),
    class = c("marginal_posterior", "numeric")
  )
  attr(posterior, "posterior_atoms") <- BayesTools::posterior_atom_attribute(
    data.frame(x = value, mass = 1)
  )
  return(posterior)
}


test_that("BF-free marginal means retain structurally fixed cells", {

  posterior <- list(zero = .fixed_marginal_posterior())
  inference <- .marginal_means_inclusion_bf(
    posterior       = posterior,
    null_hypothesis = 0,
    compute         = FALSE
  )

  expect_named(inference, "zero")
  expect_true(is.na(inference[["zero"]]))
})


test_that("structurally fixed marginal cells have unavailable BFs", {

  posterior <- list(zero = .fixed_marginal_posterior())
  testthat::local_mocked_bindings(
    Savage_Dickey_BF = function(...) {
      stop(
        paste0(
          "The posterior contains a declared point mass at the exact null ",
          "hypothesis value. The ordinary Savage-Dickey density ratio is invalid."
        ),
        call. = FALSE
      )
    },
    .package = "BayesTools"
  )
  inference <- .marginal_means_inclusion_bf(
    posterior       = posterior,
    null_hypothesis = 0,
    compute         = TRUE
  )

  expect_named(inference, "zero")
  expect_true(is.na(inference[["zero"]]))
  expect_match(
    attr(inference[["zero"]], "warnings", exact = TRUE),
    "structurally fixed"
  )
})
