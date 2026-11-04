context("Cached fits for the random-effect prior plots")

source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_refit_if_cached("plot-prior-random")

# Two fitted brma.mv models for plot_prior() and print_prior() of random-effect
# selectors: a total SD split by a Dirichlet(2, 3) allocation between two named
# blocks, and a block whose SD is modelled on the scale.

test_that("brma.mv fits with an SD allocation and a modelled SD are prepared for the prior plots", {

  withr::local_seed(4)
  data <- data.frame(
    study = rep(1:6, each = 2),
    crit  = rep(c("a", "b"), 6),
    sei   = stats::runif(12, 0.08, 0.30)
  )
  data[["yi"]] <- 0.3 + rep(stats::rnorm(6, 0, 0.15), each = 2) +
    stats::rnorm(12, 0, data[["sei"]])
  allocation <- BayesTools::random_variance_allocation(
    name    = "heterogeneity",
    terms   = c(Study = "Study", Crit = "Crit"),
    sd      = BayesTools::prior("normal", list(0, 0.4), list(0, Inf)),
    weights = BayesTools::prior("dirichlet", list(alpha = c(2, 3)))
  )
  fit_allocated <- suppressWarnings(brma.mv(
    yi = yi, V = diag(data[["sei"]]^2), data = data, measure = "GEN",
    random = list(Study = ~ 1 | study, Crit = ~ 1 | crit),
    prior_heterogeneity = BayesTools::prior_random(allocation = list(allocation)),
    prior_unit_information_sd = 1,
    chains = 1, burnin = 300, adapt = 300, sample = 500, seed = 1,
    silent = TRUE, parallel = FALSE
  ))

  withr::local_seed(2)
  scale_data <- data.frame(
    study = rep(1:6, each = 2),
    x     = stats::rnorm(12),
    sei   = stats::runif(12, 0.08, 0.30)
  )
  scale_data[["yi"]] <- 0.3 + rep(stats::rnorm(6, 0, 0.15), each = 2) +
    stats::rnorm(12, 0, scale_data[["sei"]])
  fit_scale <- suppressWarnings(brma.mv(
    yi = yi, V = diag(scale_data[["sei"]]^2), data = scale_data, measure = "GEN",
    random = ~ 1 | study, scale = ~ x, prior_unit_information_sd = 1,
    chains = 1, burnin = 300, adapt = 300, sample = 500, seed = 1,
    silent = TRUE, parallel = FALSE
  ))

  save_fit("plot_prior_random_allocated", fit_allocated)
  save_fit("plot_prior_random_scale", fit_scale)

  expect_true(is_cached_fit_valid("plot_prior_random_allocated"))
  expect_true(is_cached_fit_valid("plot_prior_random_scale"))
})
