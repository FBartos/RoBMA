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

test_that("poor Gaussian proposal coverage leaves the ordinary samplers and RNG draws intact", {

  skip_on_cran()
  previous <- RoBMA.options()
  on.exit(do.call(RoBMA.options, previous), add = TRUE)
  RoBMA.options(sampler.probabilities = character(), selection.sampler = "slice.coarse")
  fixture <- function(amplitude, se, half_rows) {
    x <- seq_len(half_rows)
    y <- rep(c(-amplitude, amplitude), length.out = half_rows)
    data.frame(y = c(y, rev(y)), se = se, x = c(-rev(x), x))
  }
  fit <- function(dat, policy) {
    RoBMA.options(sampler.joint = policy)
    BMA(yi = y, sei = se, mods = ~ x, data = dat, measure = "SMD",
        prior_effect = BayesTools::prior("normal", list(0, .7071068)),
        prior_mods = BayesTools::prior("normal", list(0, .3535534)),
        prior_heterogeneity = BayesTools::prior("normal", list(0, .02), list(0, Inf)),
        chains = 1L, adapt = 500L, burnin = 50L, sample = 100L, seed = 708L,
        convergence_checks = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL))[["fit"]]
  }

  tail <- fixture(.55, .005, 4L)
  first_rank_bin <- fixture(5.5e-6, 1.1e-6, 24L)
  for (dat in list(tail, first_rank_bin)) {
    base <- fit(dat, "off")
    for (policy in c("auto", "gaussian")) {
      declined <- fit(dat, policy)
      declarations <- attr(declined, "sampler_declarations", exact = TRUE)
      marker <- declarations[["joint"]][["markers"]][[1L]]
      expect_identical(marker[["tag"]], "robma_gaussian_joint")
      expect_identical(declarations[["robma"]][["joint"]][["status"]], "native_declined")
      expect_identical(declined[["samplers"]], base[["samplers"]])
      expect_identical(declined[["mcmc"]], base[["mcmc"]])
      expect_identical(declined[["end_state"]], base[["end_state"]])
    }
  }
})

test_that("automatic Gaussian sampling budgets retained estimates before compilation", {

  skip_on_cran()
  previous <- RoBMA.options()
  on.exit(do.call(RoBMA.options, previous), add = TRUE)
  RoBMA.options(max_cores = 2L, sampler.probabilities = character(), selection.sampler = "slice.coarse")
  rows <- seq_len(801L)
  dat <- data.frame(y = .14 + .05 * sin(rows) + .001 * rows, se = .2, x = rows)
  fit <- function(policy, keep = rep(TRUE, nrow(dat)), standardized = TRUE, parallel = FALSE) {
    RoBMA.options(sampler.joint = policy)
    BMA(yi = y, sei = se, mods = ~ x, data = dat, subset = keep, measure = "SMD",
        standardize_continuous_predictors = standardized,
        chains = if (parallel) 2L else 1L, parallel = parallel,
        adapt = 500L, burnin = 50L, sample = 100L, seed = 165L,
        convergence_checks = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL))[["fit"]]
  }
  declarations <- function(backend) attr(backend, "sampler_declarations", exact = TRUE)
  report <- function(backend) declarations(backend)[["robma"]][["joint"]]
  tag <- function(backend) declarations(backend)[["joint"]][["markers"]][[1L]][["tag"]]

  boundary <- fit("auto", keep = rows <= 800L)
  expect_identical(tag(boundary), "robma_gaussian_joint")
  expect_true(report(boundary)[["status"]] %in% c("active", "native_declined"))

  bypassed <- fit("auto", parallel = TRUE)
  ordinary <- fit("off", parallel = TRUE)
  expect_null(tag(bypassed))
  expect_identical(report(bypassed)[["status"]], "budget_bypassed")
  expect_match(report(bypassed)[["reason"]], "801 retained estimates.*limit of 800")
  expect_identical(report(bypassed)[["sampler"]], "RoBMA::GaussianJoint")
  expect_length(report(bypassed)[["nodes"]], 0L)
  expect_identical(bypassed[["samplers"]], ordinary[["samplers"]])
  expect_identical(bypassed[["mcmc"]], ordinary[["mcmc"]])
  expect_identical(bypassed[["end_state"]], ordinary[["end_state"]])

  explicit <- fit("gaussian")
  expect_identical(tag(explicit), "robma_gaussian_joint")
  expect_true(report(explicit)[["status"]] %in% c("active", "native_declined"))

  unsupported <- fit("auto", standardized = FALSE)
  expect_null(tag(unsupported))
  expect_identical(report(unsupported)[["status"]], "graph_declined")
})
