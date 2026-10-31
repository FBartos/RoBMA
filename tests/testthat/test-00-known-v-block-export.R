context("(00) known-V block export keeps interleaved rows")

# The sampling covariance of four estimates in two dependent pairs whose rows
# are interleaved: V[1, 3] and V[2, 4] are the nonzero covariances.
.known_v_interleaved <- function() {

  V <- diag(c(.1, .2, .3, .4))
  V[1, 3] <- V[3, 1] <- .03
  V[2, 4] <- V[4, 2] <- .04

  V
}


test_that("the compact export of interleaved blocks is the matrix", {

  V         <- .known_v_interleaved()
  canonical <- .known_v_canonicalize(V, warn_singular = FALSE)
  expect_identical(lapply(.known_v_blocks(canonical), `[[`, "index"),
                   list(c(1L, 3L), c(2L, 4L)))

  converted <- .known_v_as_input(canonical)
  expect_true(is.matrix(converted))
  expect_identical(unname(converted), V)
  prepared <- .known_v_prepare(converted, rep(TRUE, 4L), "auto",
                               known_v_is_scale = TRUE, warn_singular = FALSE)
  expect_identical(unname(.known_v_materialize(prepared)), V)
  expect_identical(lapply(.known_v_blocks(prepared), `[[`, "index"),
                   list(c(1L, 3L), c(2L, 4L)))
})


test_that("blocks over contiguous rows keep their compact export", {

  V <- diag(c(.1, .2, .3, .4, .5))
  V[1, 2] <- V[2, 1] <- .03
  V[3, 4] <- V[4, 3] <- .04
  canonical <- .known_v_canonicalize(V, warn_singular = FALSE)
  converted <- .known_v_as_input(canonical)
  expect_true(is.list(converted) && !is.matrix(converted))
  expect_identical(lapply(converted, dim), list(c(2L, 2L), c(2L, 2L), c(1L, 1L)))
  prepared <- .known_v_prepare(converted, rep(TRUE, 5L), "auto",
                               known_v_is_scale = TRUE, warn_singular = FALSE)
  expect_identical(unname(.known_v_materialize(prepared)), V)

  # A singleton row between the rows of one block is an interleaved layout.
  V <- diag(c(.1, .2, .3))
  V[1, 3] <- V[3, 1] <- .03
  converted <- .known_v_as_input(.known_v_canonicalize(V, warn_singular = FALSE))
  expect_identical(unname(converted), V)
})


test_that("a fitted model with interleaved known-V rows keeps its covariance", {

  V   <- .known_v_interleaved()
  dat <- data.frame(yi = c(.1, .2, .3, .4), estimate = factor(1:4))
  R   <- diag(1:4)
  dimnames(R) <- list(levels(dat$estimate), levels(dat$estimate))
  object <- brma.mv(
    yi = yi, V = V, random = ~ 1 | estimate, R = R, Rscale = "none",
    data = dat, known_v_parameterization = "auto", measure = "GEN",
    prior_effect = BayesTools::prior("normal", list(mean = 0, sd = 1)),
    prior_heterogeneity = BayesTools::prior_random(
      sd = BayesTools::prior("normal", list(mean = 0, sd = 1), list(0, Inf))),
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE
  )
  retained <- .data_known_v_data(object[["data"]])

  expect_identical(unname(.known_v_materialize(retained)), V)
  expect_identical(lapply(.known_v_blocks(retained), `[[`, "index"),
                   list(c(1L, 3L), c(2L, 4L)))
})
