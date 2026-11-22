test_that("automatic Gaussian sampling uses the compiled standardized moderator", {

  skip_on_cran()
  previous <- RoBMA.options()
  on.exit(do.call(RoBMA.options, previous), add = TRUE)
  RoBMA.options(sampler.joint = "auto", sampler.blocks = TRUE,
                sampler.probabilities = character(), selection.sampler = "slice.coarse")
  dat <- data.frame(
    yi = c(-.4, .8, .1, -.3, .6, .4), vi = c(.1, .04, .2, .3, .1, .2),
    x = c(11, 13, 21, 15, 19, 23), group = factor(rep(c("a", "b"), 3L))
  )
  settings <- list(yi = quote(yi), vi = quote(vi), mods = ~ x, data = dat, measure = "SMD",
                   chains = 1L, adapt = 500L, burnin = 50L, sample = 100L, seed = 706L,
                   convergence_checks = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL))
  fit <- function(...) do.call(BMA, utils::modifyList(settings, list(...)))[["fit"]]
  report <- function(backend) attr(backend, "sampler_declarations")[["robma"]][["joint"]]
  assigned <- function(backend) backend[["samplers"]][["Sampler"]] == "RoBMA::GaussianJoint"

  auto <- fit()
  expect_equal(sum(assigned(auto)), 6L)
  expect_identical(report(auto)[["policy"]], "auto")
  expect_identical(report(auto)[["sampler"]], "RoBMA::GaussianJoint")
  expect_identical(report(auto)[["status"]], "active")

  unscaled <- fit(standardize_continuous_predictors = FALSE)
  factor_design <- fit(mods = ~ group)
  for (backend in list(unscaled, factor_design)) {
    expect_false(any(assigned(backend)))
    expect_identical(report(backend)[["status"]], "graph_declined")
  }
  known_v <- settings
  known_v[["vi"]] <- NULL
  known_v[["V"]] <- diag(dat[["vi"]])
  expect_false(any(assigned(do.call(BMA.mv, known_v)[["fit"]])))

  RoBMA.options(sampler.joint = "off")
  off <- fit(standardize_continuous_predictors = FALSE)
  expect_false(any(assigned(off)))
  expect_identical(off[["mcmc"]], unscaled[["mcmc"]])

  RoBMA.options(sampler.joint = "gaussian")
  expect_equal(sum(assigned(fit(standardize_continuous_predictors = FALSE))), 6L)
  RoBMA.options(selection.sampler = "jags.default")
  stock <- fit()
  expect_false(any(assigned(stock)))
  expect_identical(report(stock)[["status"]], "disabled")
  RoBMA.options(sampler.joint = "auto")
  expect_false(any(assigned(fit())))

  RoBMA.options(selection.sampler = "slice", sampler.blocks = FALSE)
  expect_equal(sum(assigned(fit())), 6L)
})
