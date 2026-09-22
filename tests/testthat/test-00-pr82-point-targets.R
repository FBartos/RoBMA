test_that("all-point linear targets use one declared-metadata arithmetic path", {

  samples <- matrix(c(.1, .2), ncol = 1L, dimnames = list(NULL, "other"))
  context <- .iwmde_context_ensure_caches(list(
    posterior_samples = samples,
    flat_prior_list = list(a = BayesTools::prior("point", list(1)),
                          b = BayesTools::prior("point", list(1)),
                          c = BayesTools::prior("point", list(1)))
  ))
  weights <- c(a = .1, b = .2, c = .3)
  replacement <- list(type = "linear", weights = weights)
  target <- .iwmde_linear_values(context, samples, weights)[[1L]]
  # Matrix arithmetic gives the atom's exact binary64 representation. The
  # neighboring literal is a distinct user request and must not be snapped.
  expect_false(identical(target, .6))
  component <- .iwmde_linear_components(context, replacement)
  expect_identical(component[["point_location"]], rep(target, 2L))
  expect_identical(component[["point_masses"]][["x"]], target)
  state <- list(row = samples[1L, ], row_index = 1L,
                 parameters = list(mu = 0), baseline_log_prior = 0)
  expect_identical(.iwmde_linear_value_row(context, state[["row"]], weights), target)
  expect_true(.iwmde_replace_linear_row(context, state, target, replacement)[["valid"]])
  expect_false(.iwmde_replace_linear_row(context, state, .6, replacement)[["valid"]])
  expect_true(.iwmde_replace_linear_parameters(context, state, target, replacement)[["valid"]])
  candidates <- .iwmde_build_replacement_samples(context, "contrast", c(target, .6),
                                                  list(state), replacement)
  expect_identical(candidates[["valid"]], c(TRUE, FALSE))
  expect_identical(.iwmde_predictor_linear_log_prior_delta(
    context, c(target, .6), list(state), replacement), c(0, -Inf))
})

test_that("point-mixture target locations follow persisted branch priors", {

  mixture <- BayesTools::prior_mixture(list(
    BayesTools::prior("point", list(1)), BayesTools::prior("point", list(2))
  ), is_null = c(TRUE, FALSE))
  samples <- cbind(a_indicator = c(1, 2), other = c(.1, .2))
  context <- .iwmde_context_ensure_caches(list(
    posterior_samples = samples, indicator_names = "a_indicator",
    flat_prior_list = list(a = mixture, b = BayesTools::prior("point", list(1)),
                          c = BayesTools::prior("point", list(1)))
  ))
  weights <- c(a = .1, b = .2, c = .3)
  expected <- as.numeric(rbind(c(1, 1, 1), c(2, 1, 1)) %*% weights)
  expect_identical(.iwmde_linear_values(context, samples, weights), expected)
  expect_identical(.iwmde_linear_components(context, list(weights = weights))[["point_location"]], expected)
  for (row in 1:2) {
    expect_identical(.iwmde_linear_value_row(context, samples[row, ], weights), expected[[row]])
  }
})
