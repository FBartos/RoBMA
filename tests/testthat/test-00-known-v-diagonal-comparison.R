context("A diagonal known V is the ordinary data of its standard errors")

test_that("an ordinary fit and an identical diagonal known-V fit compare", {

  withr::local_seed(2631)
  k       <- 8L
  yi      <- stats::rnorm(k, .3, .2)
  sei     <- stats::runif(k, .1, .2)
  prior_effect <- BayesTools::prior("normal", list(mean = 0, sd = 1))

  ordinary <- brma(
    yi = yi, sei = sei, measure = "GEN", prior_unit_information_sd = 1,
    prior_effect = prior_effect, prior_heterogeneity = NULL,
    only_priors = TRUE, silent = TRUE
  )
  known_v <- function(V, ...) {
    brma.mv(
      yi = yi, V = V, measure = "GEN", prior_unit_information_sd = 1,
      prior_effect = prior_effect, only_priors = TRUE, silent = TRUE, ...
    )
  }
  diagonal <- known_v(diag(sei^2))
  correlated_V <- diag(sei^2)
  correlated_V[1L, 2L] <- correlated_V[2L, 1L] <- .5 * sei[1L] * sei[2L]
  correlated <- known_v(correlated_V, known_v_parameterization = "block_mvn")

  # no correlated block: the equivalent ordinary data, however V is supplied
  expect_true(.is_data_known_v(diagonal[["data"]]))
  expect_identical(.get_outcome_hash(diagonal), .get_outcome_hash(ordinary))
  for (V in list(sei^2, lapply(sei^2, function(x) matrix(x, 1L, 1L)))) {
    expect_identical(.get_outcome_hash(known_v(V)), .get_outcome_hash(ordinary))
  }
  from_variances <- brma(
    yi = yi, vi = sei^2, measure = "GEN", prior_unit_information_sd = 1,
    prior_effect = prior_effect, prior_heterogeneity = NULL,
    only_priors = TRUE, silent = TRUE
  )
  expect_identical(.get_outcome_hash(from_variances), .get_outcome_hash(ordinary))

  # a correlated block, or other data, is a different target
  expect_false(identical(.get_outcome_hash(correlated), .get_outcome_hash(ordinary)))
  expect_false(identical(.get_outcome_hash(correlated), .get_outcome_hash(diagonal)))
  other_sei <- brma(
    yi = yi, sei = sei * 1.01, measure = "GEN", prior_unit_information_sd = 1,
    prior_effect = prior_effect, prior_heterogeneity = NULL,
    only_priors = TRUE, silent = TRUE
  )
  expect_false(identical(.get_outcome_hash(other_sei), .get_outcome_hash(diagonal)))

  # marginal likelihood comparisons
  stored_marglik <- function(logml) {
    BayesTools:::.bt_marglik_from_upstream(
      structure(
        list(logml = logml, niter = 4L, method = "normal", repetitions = 1L),
        class = "bridge"
      ),
      maxiter = 100, nonfinite = "error",
      chain_metadata = list(count = 2L, draws = 512L, draws_per_chain = c(256L, 256L))
    )
  }
  with_marglik <- function(object, logml) {
    structure(list(data = object[["data"]], marglik = stored_marglik(logml)),
              class = "brma")
  }
  model_ordinary   <- with_marglik(ordinary, 0)
  model_diagonal   <- with_marglik(diagonal, log(2))
  model_correlated <- with_marglik(correlated, log(2))

  expect_equal(bf(model_ordinary, model_diagonal)[["bf"]], .5, tolerance = 1e-12)
  expect_equal(unname(post_prob(model_ordinary, model_diagonal)), c(1, 2) / 3,
               tolerance = 1e-12)
  expect_error(bf(model_ordinary, model_correlated), "same outcome data")
  expect_error(post_prob(model_ordinary, model_correlated), "same outcome data")

  # LOO comparisons: the same posterior draws give the same pointwise scores
  precision <- 1 + sum(1 / sei^2)
  draws     <- stats::rnorm(400L, sum(yi / sei^2) / precision, sqrt(1 / precision))
  replay <- function(object) {
    review_attach_replay(review_compile_object(object), cbind(mu = draws))
  }
  # (eight estimates leave some importance weights heavy-tailed; the scores
  # compared are the same computation for both fits)
  add_loo_quietly <- function(object) suppressWarnings(add_loo(replay(object)))
  loo_ordinary <- add_loo_quietly(ordinary)
  loo_diagonal <- add_loo_quietly(diagonal)
  comparison   <- loo_compare(loo_ordinary, loo_diagonal)
  expect_equal(unname(comparison[, "elpd_diff"]), c(0, 0), tolerance = 1e-10)
  expect_error(
    loo_compare(loo_ordinary, add_loo_quietly(correlated)),
    "different data, unit, or retained-context targets"
  )
})
