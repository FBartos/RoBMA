context("qCMDE/IWMDE classify inclusion indicators by their declared role")

# A normal meta-regression on one continuous predictor 'predictor' with 400
# synthetic draws; the heterogeneity has the prior 'tau_prior'.
.review_predictor_replay <- function(predictor, tau_prior) {

  dat <- data.frame(yi = c(.1, .25, .3, .4, .2, .15, .05, .35), sei = .2)
  dat[[predictor]] <- c(-1, 0, 1, 2, 3, 1, -2, 0)
  object <- brma(
    yi = dat[["yi"]], sei = dat[["sei"]],
    mods = stats::as.formula(paste("~", predictor)), data = dat, measure = "GEN",
    prior_unit_information_sd = 1,
    prior_effect = BayesTools::prior("normal", list(0, 1)),
    prior_mods = BayesTools::prior("normal", list(0, 1)),
    prior_heterogeneity = tau_prior, only_priors = TRUE, silent = TRUE
  )
  setup <- review_compile_object(object)
  withr::local_seed(11)
  n <- 400L
  samples <- cbind(
    stats::rnorm(n, .1, .1), stats::rnorm(n, .05, .1),
    abs(stats::rnorm(n, .15, .08)) + .02
  )
  colnames(samples) <- c(
    review_coefficient_name(setup, "intercept"),
    review_coefficient_name(setup, predictor), "tau"
  )

  review_attach_replay(setup, samples)
}

.review_indicator_bf <- function(object, statement) {

  hypothesis(
    object, statement, columns = "all", density_method = "qCMDE",
    density_control = list(n_points = 20, samples = 100), seed = 1
  )
}

test_that("qCMDE hypotheses run for a continuous predictor named like an indicator", {

  tau_prior <- BayesTools::prior("normal", list(0, 1), truncation = list(0, Inf))
  object  <- .review_predictor_replay("treatment_indicator", tau_prior)
  control <- .review_predictor_replay("treatment_flag", tau_prior)

  # point, linear and non-linear targets run like those of the renamed predictor
  for (statement in c("treatment_indicator = 0", "intercept = 0", "tau = 0.15")) {
    renamed  <- sub("treatment_indicator", "treatment_flag", statement)
    expected <- .review_indicator_bf(control, renamed)
    observed <- tryCatch(.review_indicator_bf(object, statement), error = identity)
    expect_s3_class(observed, "data.frame")
    if (is.data.frame(observed)) {
      expect_true(is.finite(observed[["BF"]][[1L]]), info = statement)
      expect_equal(observed[["BF"]], expected[["BF"]], tolerance = 1e-8, info = statement)
    }
  }
})

test_that("a declared coefficient named like an indicator is no indicator", {

  tau_prior <- BayesTools::prior("normal", list(0, 1), truncation = list(0, Inf))
  context <- .iwmde_context(.review_predictor_replay("treatment_indicator", tau_prior))

  expect_identical(context[["coefficient_names"]],
                   c("mu_intercept", "mu_treatment_indicator"))
  expect_identical(context[["indicator_names"]], character())
  expect_false(.iwmde_parameter_is_indicator("mu_treatment_indicator", context))
  expect_identical(
    .iwmde_parameter_spec(context, "mu_treatment_indicator")[["status"]], "ok"
  )
})

test_that("model-inclusion indicators stay indicators beside a coefficient named like one", {

  dat <- data.frame(
    yi = c(.1, .2, .3, .4, .2, .1), sei = .2,
    x = c(-1, 0, 1, 2, 0, 1), treatment_indicator = c(-1, 0, 1, 2, 3, 1)
  )
  object <- RoBMA(
    yi = yi, sei = sei, mods = ~ x + treatment_indicator, data = dat,
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE
  )
  setup <- review_compile_object(object)
  priors <- list()
  for (compiled in setup[["compiled"]]) priors <- c(priors, compiled[["prior_list"]])
  priors <- c(priors, .create_fit_priors(object[["data"]], object[["priors"]]))
  columns <- c(
    review_prior_columns(priors), "mu_x_indicator", "mu_treatment_indicator_indicator",
    "tau_indicator", "bias_indicator"
  )
  samples <- matrix(.5, 20L, length(columns), dimnames = list(NULL, columns))
  samples[, grep("_indicator$", columns)[-match("mu_treatment_indicator", columns, nomatch = 0L)]] <-
    rep(c(1, 2), 10)
  object <- review_attach_replay(setup, samples)
  context <- .iwmde_context(object)

  expect_true("mu_treatment_indicator" %in% context[["coefficient_names"]])
  expect_setequal(
    context[["indicator_names"]],
    c("mu_x_indicator", "mu_treatment_indicator_indicator", "tau_indicator",
      "bias_indicator")
  )
  for (name in context[["indicator_names"]]) {
    expect_true(.iwmde_parameter_is_indicator(name, context), info = name)
  }
  expect_false(.iwmde_parameter_is_indicator("mu_treatment_indicator", context))
})

test_that("contexts without parameter-map metadata keep the suffix rule", {

  samples <- cbind(theta_indicator = c(1, 2), mu_treatment_indicator = c(.1, .2), tau = 1)

  expect_identical(
    .iwmde_indicator_names(samples),
    c("theta_indicator", "mu_treatment_indicator")
  )
  expect_identical(
    .iwmde_indicator_names(samples, coefficient_names = "mu_treatment_indicator"),
    "theta_indicator"
  )
  expect_true(.iwmde_parameter_is_indicator("theta_indicator", list(indicator_names = "theta_indicator")))
  expect_true(.iwmde_parameter_is_indicator("bias_indicator", list()))
  expect_false(.iwmde_parameter_is_indicator("tau", list()))
  expect_null(.iwmde_declared_coefficient_names(list()))
})
