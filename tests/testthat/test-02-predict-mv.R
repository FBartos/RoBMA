context("Predict brma.mv")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))

fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)

.brma_mv_prior_object <- function(random = FALSE) {

  dat <- data.frame(
    yi     = c(0.08, 0.13, 0.18, 0.20),
    study  = rep(c("s1", "s2"), each = 2),
    effect = rep(c("a", "b"), 2),
    x      = c(0, 1, 0, 1)
  )
  V <- matrix(0, nrow = 4, ncol = 4)
  block <- matrix(c(0.04, 0.018, 0.018, 0.05), nrow = 2)
  V[1:2, 1:2] <- block
  V[3:4, 3:4] <- block

  args <- list(
    yi                        = quote(yi),
    V                         = V,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  if (random) {
    args[["random"]] <- quote(~ 1 | study / effect)
  }

  do.call(brma.mv, args)
}

.brma_mv_random_sd_names <- function(object) {

  random_effects <- object[["formula_design"]][["mu"]][["random_effects"]]
  stats::setNames(
    vapply(random_effects, function(term) {
      term[["sd_parameter_names"]][1L]
    }, character(1)),
    vapply(random_effects, function(term) {
      term[["block_name"]]
    }, character(1))
  )
}

test_that("brma.mv newdata uses fitted continuous-predictor scaling metadata", {

  object <- brma.mv(
    yi                        = yi,
    V                         = diag(c(0.04, 0.05, 0.06)),
    mods                      = ~ x,
    data                      = data.frame(
      yi = c(0.10, 0.20, 0.30),
      x  = c(10, 20, 30)
    ),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  posterior_samples <- matrix(
    c(1, 2, 0.10, 0.20, 0, 0),
    nrow     = 2L,
    dimnames = list(NULL, c("mu_intercept", "mu_x", "tau"))
  )
  newdata <- data.frame(x = c(10, 20, 30, 40))

  prediction <- predict(
    object,
    newdata            = newdata,
    type               = "terms",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )
  expected <- cbind(
    posterior_samples[, "mu_intercept"] + 10 * posterior_samples[, "mu_x"],
    posterior_samples[, "mu_intercept"] + 20 * posterior_samples[, "mu_x"],
    posterior_samples[, "mu_intercept"] + 30 * posterior_samples[, "mu_x"],
    posterior_samples[, "mu_intercept"] + 40 * posterior_samples[, "mu_x"]
  )

  expect_equal(unname(as.matrix(prediction)), expected, tolerance = 1e-12)
  expect_equal(
    attr(prediction, "data")[["mods"]][["x"]],
    newdata[["x"]]
  )
})


test_that("type cluster is rejected for brma.mv models", {

  nonrandom_object <- .brma_mv_prior_object(random = FALSE)
  random_object    <- .brma_mv_prior_object(random = TRUE)

  expect_error(
    predict(nonrandom_object, type = "cluster"),
    "cluster.*multilevel|multilevel.*cluster"
  )
  expect_error(
    predict(random_object, type = "cluster"),
    "cluster.*multilevel|multilevel.*cluster"
  )
})

test_that("known-V newdata response predictions require and validate V_new", {

  object <- .brma_mv_prior_object(random = FALSE)
  newdata <- data.frame(row = 1:2)
  V_new   <- matrix(c(0.04, 0.01, 0.01, 0.05), nrow = 2)
  posterior_samples <- matrix(
    c(0.1, 0.2, 0.05, 0.06),
    nrow = 2,
    dimnames = list(NULL, c("mu", "tau"))
  )

  expect_error(
    predict(object, newdata = newdata, type = "response"),
    "V_new"
  )
  expect_error(
    predict(
      object,
      newdata            = data.frame(sei = c(0.10, 0.10)),
      V_new              = V_new,
      type               = "response",
      .posterior_samples = posterior_samples
    ),
    "diag\\(V_new\\)"
  )
  expect_error(
    predict(
      object,
      newdata            = newdata,
      V_new              = matrix(0.04, nrow = 1, ncol = 1),
      type               = "response",
      .posterior_samples = posterior_samples
    ),
    "dimensions.*V_new"
  )
  expect_error(
    predict(
      object,
      newdata            = newdata,
      V_new              = matrix(c(0.04, 0.02, 0.01, 0.05), nrow = 2),
      type               = "response",
      .posterior_samples = posterior_samples
    ),
    "symmetric"
  )
  expect_error(
    predict(
      object,
      newdata            = newdata,
      V_new              = matrix(c(0.04, 0.08, 0.08, 0.05), nrow = 2),
      type               = "response",
      .posterior_samples = posterior_samples
    ),
    "positive semidefinite"
  )

  set.seed(1)
  response <- predict(
    object,
    newdata            = newdata,
    V_new              = V_new,
    type               = "response",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )
  expect_brma_samples_matrix(response, 2, "known-V V_new response")
  expect_equal(attr(response, "data")[["outcome"]][["sei"]], sqrt(diag(V_new)))
  response_target <- attr(response, "brma_mv_prediction_target")
  expect_equal(response_target[["type"]], "response")
  expect_true(response_target[["known_v"]])
  expect_false(response_target[["random_formula"]])
  expect_true(response_target[["v_new"]])
  expect_equal(response_target[["mean_target"]], "fixed_location")
  expect_equal(response_target[["covariance_target"]], "V_new_plus_heterogeneity")

  set.seed(1)
  response_list <- predict(
    object,
    newdata            = newdata,
    V_new              = list(matrix(0.04, 1), matrix(0.05, 1)),
    type               = "response",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )
  expect_brma_samples_matrix(response_list, 2, "known-V block-list V_new response")
})


test_that("known-V prediction compares supplied sampling scales row-wise", {

  expect_silent(.predict_known_v_newdata_check_variance(
    supplied = c(1e-300, 1e300),
    expected = c(1e-300, 1e300),
    label    = "vi"
  ))
  expect_error(
    .predict_known_v_newdata_check_variance(
      supplied = c(1e-300 * (1 + 1e-10), 1e300),
      expected = c(1e-300, 1e300),
      label    = "vi"
    ),
    "must match diag\\(V_new\\)"
  )

  actual <- .predict_known_v_newdata_add_vi(
    newdata = data.frame(sei = c(1e-150, 1e150)),
    V_new   = diag(c(1e-300, 1e300))
  )
  expect_identical(actual[["vi"]], c(1e-300, 1e300))
  expect_false("sei" %in% names(actual))
})


test_that("known-V newdata response predictions preserve V_new covariance", {

  object <- .brma_mv_prior_object(random = FALSE)
  S      <- 20000L
  V_new  <- matrix(c(0.04, 0.018, 0.018, 0.05), nrow = 2)
  posterior_samples <- cbind(
    mu  = rep(0, S),
    tau = rep(0, S)
  )

  set.seed(11)
  response <- predict(
    object,
    newdata            = data.frame(row = 1:2),
    V_new              = V_new,
    type               = "response",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )

  expect_equal(
    unname(stats::cov(as.matrix(response))),
    unname(V_new),
    tolerance = 0.004
  )
})

test_that("known-V weightfunction response prediction does not fall back to diagonal selection RNG", {

  object <- .brma_mv_prior_object(random = FALSE)
  object[["priors"]][["outcome"]][["bias"]] <- BayesTools::prior_weightfunction(
    "two-sided",
    steps   = c(0.05),
    weights = BayesTools::wf_fixed(c(1, 1))
  )

  expect_error(
    predict(object, type = "response", bias_adjusted = FALSE),
    "known-V covariance"
  )
})

test_that("random-formula brma.mv predict uses cached fit", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  fit_brma  <- fits[[name]]
  n_studies <- nobs(fit_brma)

  terms <- predict(fit_brma, type = "terms", quiet = TRUE)
  scale <- predict(fit_brma, type = "terms.scale", quiet = TRUE)

  expect_brma_samples_matrix(terms, n_studies, "random brma.mv terms")
  expect_type(scale, "list")
  expect_equal(names(scale), "Component 1")
  expect_brma_samples_matrix(scale[["Component 1"]], n_studies, "random brma.mv terms.scale")
  expect_equal(
    colnames(scale[["Component 1"]]),
    paste0("tau[", seq_len(n_studies), "]")
  )
})

test_that("3-level location-scale brma.mv fixtures expose targeted SD components", {

  model_names <- c(
    total  = "brma.mv_block_mvn_3lvl_scale_total",
    top    = "brma.mv_block_mvn_3lvl_scale_top",
    bottom = "brma.mv_block_mvn_3lvl_scale_bottom"
  )
  skip_if_missing_fits(unname(model_names))

  expected_components <- list(
    total  = "total",
    top    = c("study", "effect"),
    bottom = c("study", "effect")
  )
  expected_parameters <- list(
    total  = c(tau = "log_tau"),
    top    = c(study = "log_tau_study", effect = "log_tau_effect"),
    bottom = c(study = "log_tau_study", effect = "log_tau_effect")
  )
  expected_rows <- list(
    total  = c("(log_tau) x"),
    top    = c("(log_tau_study) x", "(log_tau_effect) intercept"),
    bottom = c("(log_tau_study) intercept", "(log_tau_effect) x")
  )
  unexpected_rows <- list(
    total  = character(),
    top    = c("(log_tau_effect) x"),
    bottom = c("(log_tau_study) x")
  )

  for (target in names(model_names)) {
    fit_brma <- fits[[model_names[[target]]]]
    n_studies <- nobs(fit_brma)

    terms_scale <- predict(fit_brma, type = "terms.scale", quiet = TRUE)
    new_scale   <- predict(
      fit_brma,
      newdata = data.frame(x = c(0, 1)),
      type    = "terms.scale",
      quiet   = TRUE
    )
    scale_fitted <- fitted(fit_brma, component = "scale")

    expect_equal(names(terms_scale), expected_components[[target]], info = target)
    expect_equal(names(new_scale), expected_components[[target]], info = target)
    expect_equal(names(scale_fitted), expected_components[[target]], info = target)
    expect_equal(
      .data_scale_formula_parameters(fit_brma[["data"]]),
      expected_parameters[[target]],
      info = target
    )

    for (row in expected_rows[[target]]) {
      expect_true(
        any(grepl(row, rownames(fit_brma[["summary"]]), fixed = TRUE)),
        info = paste(target, row)
      )
    }
    for (row in unexpected_rows[[target]]) {
      expect_false(
        any(grepl(row, rownames(fit_brma[["summary"]]), fixed = TRUE)),
        info = paste(target, row)
      )
    }

    for (component in expected_components[[target]]) {
      expect_brma_samples_matrix(
        terms_scale[[component]],
        n_studies,
        paste(target, component, "terms.scale")
      )
      expect_brma_samples_matrix(
        new_scale[[component]],
        2L,
        paste(target, component, "newdata terms.scale")
      )
      expect_equal(
        colnames(terms_scale[[component]]),
        paste0("tau[", seq_len(n_studies), "]"),
        info = paste(target, component)
      )
      expect_equal(
        unname(scale_fitted[[component]]),
        unname(colMeans(as.matrix(terms_scale[[component]]))),
        tolerance = 1e-12,
        info      = paste(target, component)
      )
    }
  }
})

test_that("random-formula brma.mv fitted scale summarizes list-valued predictions", {

  object <- .brma_mv_prior_object(random = TRUE)
  n_studies <- nobs(object)

  location_samples <- matrix(
    seq_len(2 * n_studies),
    nrow = 2,
    dimnames = list(NULL, paste0("mu[", seq_len(n_studies), "]"))
  )
  scale_samples <- matrix(
    seq_len(2 * n_studies) / 10,
    nrow = 2,
    dimnames = list(NULL, paste0("tau[", seq_len(n_studies), "]"))
  )

  testthat::local_mocked_bindings(
    predict.brma = function(object, newdata = NULL, type = "response",
                            conditional = FALSE, quiet = TRUE, ...) {

      if (identical(type, "terms.scale")) {
        return(list("Component 1" = scale_samples))
      }
      if (identical(type, "terms")) {
        return(location_samples)
      }
      stop("Unexpected prediction type in fitted() regression test.",
           call. = FALSE)
    },
    .package = "RoBMA"
  )

  scale <- fitted(object, component = "scale")
  all   <- fitted(object, component = "all")

  expect_type(scale, "list")
  expect_equal(names(scale), "Component 1")
  expect_equal(scale[["Component 1"]], .diagnostic_set_names(colMeans(scale_samples), object))

  expect_type(all[["scale"]], "list")
  expect_equal(all[["location"]], .diagnostic_set_names(colMeans(location_samples), object))
  expect_equal(all[["scale"]][["Component 1"]], scale[["Component 1"]])
})

test_that("cached brma.mv fitted wrappers agree with prediction targets", {

  model_names <- c(
    "brma.mv_block_mvn_fixed_random_null",
    "brma.mv_block_mvn_random_scale"
  )
  skip_if_missing_fits(model_names)

  for (name in model_names) {
    fit_brma <- fits[[name]]
    terms    <- as.matrix(predict(fit_brma, type = "terms", quiet = TRUE))
    estimate <- as.matrix(predict(fit_brma, type = "estimate", quiet = TRUE))

    expect_equal(
      unname(fitted(fit_brma)),
      unname(colMeans(terms)),
      tolerance = 1e-12,
      info      = paste(name, "location fitted")
    )
    expect_equal(
      unname(fitted(fit_brma, conditioning_depth = "estimate")),
      unname(colMeans(estimate)),
      tolerance = 1e-12,
      info      = paste(name, "estimate-depth fitted")
    )
    expect_error(
      fitted(fit_brma, unit = "cluster"),
      "cluster",
      info = name
    )
  }

  fit_scale     <- fits[["brma.mv_block_mvn_random_scale"]]
  scale_samples <- predict(fit_scale, type = "terms.scale", quiet = TRUE)
  scale_fitted  <- fitted(fit_scale, component = "scale")
  all_fitted    <- fitted(fit_scale, component = "all")

  expect_type(scale_fitted, "list")
  expect_equal(names(scale_fitted), names(scale_samples))
  expect_equal(names(all_fitted), c("location", "scale"))
  expect_equal(all_fitted[["scale"]], scale_fitted)
  for (component in names(scale_samples)) {
    expect_equal(
      unname(scale_fitted[[component]]),
      unname(colMeans(as.matrix(scale_samples[[component]]))),
      tolerance = 1e-12,
      info      = component
    )
  }
  expect_error(
    fitted(fit_scale, component = "scale", transform = exp),
    "only available for location fitted values",
    fixed = TRUE
  )
})

test_that("random-formula brma.mv estimate and response predictions use cached fit", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  fit_brma  <- fits[[name]]
  n_studies <- nobs(fit_brma)

  estimate <- predict(fit_brma, type = "estimate", quiet = TRUE)
  set.seed(1)
  response <- predict(fit_brma, type = "response", quiet = TRUE)

  expect_brma_samples_matrix(estimate, n_studies, "random brma.mv estimate")
  expect_brma_samples_matrix(response, n_studies, "random brma.mv response")
  expect_equal(attr(estimate, "title"), "Conditional True Effects:")
})

test_that("known-V random response draws combine generative components", {

  object            <- .brma_mv_prior_object(random = TRUE)
  posterior_samples <- matrix(0, nrow = 2L, ncol = 1L)
  mu_samples        <- matrix(c(.1, .2, .3, .4), nrow = 2L, byrow = TRUE)
  random_noise      <- matrix(c(.5, .6, .7, .8), nrow = 2L, byrow = TRUE)
  sampling_noise    <- matrix(c(.01, .02, .03, .04), nrow = 2L, byrow = TRUE)
  known_V           <- .known_v_newdata_prepare(diag(c(.04, .05)), k = 2L)
  captured          <- NULL

  testthat::local_mocked_bindings(
    .predict_brma_mv_new_effect_random_draws = function(object, data,
                                                        posterior_samples,
                                                        max_bytes = NULL) {
      random_noise
    },
    .known_v_sampling_noise = function(known_V, S, K) {
      captured <<- list(known_V = known_V, S = S, K = K)
      sampling_noise
    },
    .package = "RoBMA"
  )

  out <- .predict_brma_mv_known_v_response_draws(
    object            = object,
    data              = list(outcome = data.frame(yi = c(0, 0))),
    known_V           = known_V,
    mu_samples        = mu_samples,
    posterior_samples = posterior_samples,
    same_data         = FALSE
  )

  expect_equal(
    as.numeric(out),
    as.numeric(mu_samples + random_noise + sampling_noise)
  )
  expect_identical(captured[["known_V"]], known_V)
  expect_identical(captured[c("S", "K")], list(S = 2L, K = 2L))
  generation <- attr(out, "brma_mv_response_generation")
  expect_identical(
    generation[c("method", "same_data")],
    list(method = "component_generative", same_data = FALSE)
  )
  expect_identical(generation[["chunk_size"]], 2L)
  expect_identical(generation[["n_chunks"]], 1L)
  expect_lte(generation[["working_peak"]], generation[["max_bytes"]])
})

test_that("known-V response generation bounds component workspace by chunks", {

  S                    <- 5L
  K                    <- 2L
  object               <- .brma_mv_prior_object(random = TRUE)
  posterior_samples    <- matrix(seq_len(S), ncol = 1L)
  mu_samples           <- matrix(seq_len(S * K), nrow = S, ncol = K)
  known_V              <- .known_v_newdata_prepare(diag(c(.04, .05)), k = K)
  max_bytes            <- .predict_brma_mv_response_peak_bytes(2L, K)
  random_chunk_sizes   <- integer()
  sampling_chunk_sizes <- integer()

  testthat::local_mocked_bindings(
    .predict_brma_mv_new_effect_random_draws = function(object, data,
                                                        posterior_samples,
                                                        max_bytes = NULL) {
      random_chunk_sizes <<- c(random_chunk_sizes, nrow(posterior_samples))
      matrix(0, nrow = nrow(posterior_samples), ncol = K)
    },
    .known_v_sampling_noise = function(known_V, S, K) {
      sampling_chunk_sizes <<- c(sampling_chunk_sizes, S)
      matrix(0, nrow = S, ncol = K)
    },
    .package = "RoBMA"
  )

  out <- .predict_brma_mv_known_v_response_draws(
    object            = object,
    data              = list(outcome = data.frame(yi = numeric(K))),
    known_V           = known_V,
    mu_samples        = mu_samples,
    posterior_samples = posterior_samples,
    same_data         = FALSE,
    max_bytes         = max_bytes
  )

  expect_equal(as.numeric(out), as.numeric(mu_samples))
  expect_identical(random_chunk_sizes, c(2L, 2L, 1L))
  expect_identical(sampling_chunk_sizes, c(2L, 2L, 1L))
  expect_identical(attr(out, "brma_mv_response_generation")[["chunk_size"]], 2L)
  expect_identical(attr(out, "brma_mv_response_generation")[["n_chunks"]], 3L)
  expect_lte(
    attr(out, "brma_mv_response_generation")[["working_peak"]],
    max_bytes
  )
})


test_that("random-formula brma.mv newdata estimate predictions use BayesTools targets", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  newdata_scale  <- data.frame(
    x = c(0, 1)
  )
  newdata_random <- data.frame(
    study  = c("s1", "s2"),
    effect = c("a", "b"),
    x      = c(0, 1)
  )

  scale <- predict(fit_brma, newdata = newdata_scale, type = "terms.scale", quiet = TRUE)
  expect_type(scale, "list")
  expect_equal(names(scale), "Component 1")
  expect_brma_samples_matrix(scale[["Component 1"]], nrow(newdata_scale), "random brma.mv newdata terms.scale")

  estimate <- predict(fit_brma, newdata = newdata_random, type = "estimate", quiet = TRUE)
  expect_brma_samples_matrix(
    estimate,
    nrow(newdata_random),
    "random brma.mv newdata estimate"
  )
  expect_error(
    predict(fit_brma, newdata = newdata_random, type = "response", quiet = TRUE),
    "V_new"
  )
})


test_that("explicit random-formula rows use marginal new-effect targets", {

  object <- brma.mv(
    yi                         = yi,
    V                          = diag(c(.04, .05)),
    data                       = data.frame(
      yi    = c(.1, .2),
      study = c("s1", "s2")
    ),
    random                     = ~ 1 | study,
    marginalize_estimate_level = FALSE,
    measure                    = "GEN",
    prior_unit_information_sd  = 1,
    only_priors                = TRUE
  )
  new_data <- .prepare_newdata(
    object         = object,
    newdata        = data.frame(.prediction_row = 1:2),
    type           = "estimate",
    bias_adjusted  = TRUE,
    include_scale  = TRUE,
    include_random = TRUE
  )
  posterior_samples <- matrix(0, nrow = 3L, ncol = 1L)
  calls             <- list()

  testthat::local_mocked_bindings(
    .evaluate.brma.random_effects = function(fit, data, priors,
                                             posterior_samples = NULL,
                                             same_data = TRUE,
                                             required = FALSE,
                                             formula_target = "conditional",
                                             blocks = NULL,
                                             object = NULL) {
      calls[[length(calls) + 1L]] <<- list(
        same_data      = same_data,
        formula_target = formula_target,
        blocks         = blocks
      )
      matrix(.25, nrow = nrow(posterior_samples), ncol = nrow(data[["outcome"]]))
    },
    .package = "RoBMA"
  )

  draws <- .predict_brma_mv_new_effect_random_draws(
    object            = object,
    data              = new_data,
    posterior_samples = posterior_samples
  )

  expect_equal(draws, matrix(.25, nrow = 3L, ncol = 2L))
  expect_length(calls, 1L)
  expect_false(calls[[1L]][["same_data"]])
  expect_identical(calls[[1L]][["formula_target"]], "marginal")
  expect_identical(calls[[1L]][["blocks"]], "study")
  expect_equal(length(unique(new_data[["location"]][["study"]])), 2L)
})


test_that("new-effect grouping synthesis does not hide missing design values", {

  object <- brma.mv(
    yi                        = yi,
    V                         = diag(c(.04, .05, .06)),
    data                      = data.frame(
      yi    = c(.1, .2, .3),
      x     = c(0, 1, 2),
      study = c("s1", "s2", "s3")
    ),
    random                    = list(
      grouped_by_x = ~ 1 | x,
      slope        = ~ us(0 + x | study)
    ),
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )

  expect_error(
    .prepare_newdata(
      object         = object,
      newdata        = data.frame(.prediction_row = 1:2),
      type           = "estimate",
      bias_adjusted  = TRUE,
      include_scale  = TRUE,
      include_random = TRUE
    ),
    "random-effect variables. Missing: x",
    fixed = TRUE
  )
})


test_that("random-formula brma.mv newdata estimate includes marginalized random draws", {

  object <- .brma_mv_prior_object(random = TRUE)
  sd_names <- .brma_mv_random_sd_names(object)
  posterior_samples <- matrix(
    0,
    nrow = 2,
    ncol = 1 + length(sd_names),
    dimnames = list(NULL, c("mu", unname(sd_names)))
  )
  posterior_samples[, "mu"] <- c(.10, .20)
  posterior_samples[, unname(sd_names)] <- .05
  newdata <- data.frame(
    study  = c("s1", "s2"),
    effect = c("a", "b"),
    x      = c(0, 1),
    yi     = c(.10, .15),
    vi     = c(.04, .05)
  )
  contribution <- matrix(.25, nrow = nrow(posterior_samples), ncol = nrow(newdata))

  testthat::local_mocked_bindings(
    .evaluate.brma.random_effects = function(fit, data, priors,
                                             posterior_samples = NULL,
                                             same_data = TRUE,
                                             required = FALSE,
                                             formula_target = "conditional",
                                             blocks = NULL,
                                             object = NULL) {
      matrix(0, nrow = nrow(posterior_samples), ncol = nrow(data[["outcome"]]))
    },
    .predict_known_v_marginalized_random_draws = function(object, data,
                                                          posterior_samples) {
      contribution
    },
    .package = "RoBMA"
  )

  estimate <- predict(
    object,
    newdata            = newdata,
    type               = "estimate",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )
  terms <- predict(
    object,
    newdata            = newdata,
    type               = "terms",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )

  expect_equal(
    unname(as.matrix(estimate) - as.matrix(terms)),
    unname(contribution),
    tolerance = 1e-12
  )
})


test_that("marginalized newdata random draws share repeated group levels", {

  object <- .brma_mv_prior_object(random = TRUE)
  term   <- .data_marginalized_random_effects(object[["data"]])[[1L]]
  sd_name <- term[["sd_parameter_names"]][[1L]]
  posterior_samples <- matrix(
    c(2, 3),
    ncol     = 1L,
    dimnames = list(NULL, sd_name)
  )
  new_data <- object[["data"]]
  new_data[["outcome"]] <- data.frame(
    yi   = c(0, 0, 0),
    sei  = c(0, 0, 0),
    slab = paste0("new", 1:3)
  )
  new_data[["location"]] <- data.frame(
    study  = c("s1", "s1", "s2"),
    effect = c("a", "a", "b"),
    x      = c(0, 0, 1)
  )

  set.seed(91)
  draws <- .predict_known_v_marginalized_random_draws(
    object            = object,
    data              = new_data,
    posterior_samples = posterior_samples
  )

  expect_equal(draws[, 1], draws[, 2], tolerance = 1e-12)
  expect_false(isTRUE(all.equal(draws[, 1], draws[, 3])))
})


test_that("marginalized newdata estimate draws match sampled random covariance", {

  dat <- data.frame(
    yi    = c(.10, .20),
    study = c("s1", "s2")
  )
  V <- diag(rep(.04, 2L))
  marginalized <- brma.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | study,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  sampled <- brma.mv(
    yi                         = yi,
    V                          = V,
    random                     = ~ 1 | study,
    data                       = dat,
    measure                    = "GEN",
    prior_unit_information_sd  = 1,
    marginalize_estimate_level = FALSE,
    only_priors                = TRUE
  )
  sampled_design <- .fitted_formula_design(sampled, "mu", required = TRUE)
  sd_name        <- sampled_design[["random_effects"]][[1L]][["sd_parameter_names"]][[1L]]
  S              <- 20000L
  posterior_samples <- matrix(
    rep(.5, S),
    ncol     = 1L,
    dimnames = list(NULL, sd_name)
  )
  new_location <- data.frame(study = c("s1", "s1", "s2"))
  new_data <- marginalized[["data"]]
  new_data[["outcome"]] <- data.frame(
    yi   = c(0, 0, 0),
    sei  = c(0, 0, 0),
    slab = paste0("new", 1:3)
  )
  new_data[["location"]] <- new_location

  formula_fit <- .posterior_formula_fit(
    fit               = sampled[["fit"]],
    posterior_samples = posterior_samples,
    formula_design    = TRUE
  )
  attr(formula_fit, "formula_design") <- list(mu = sampled_design)
  location_priors <- attr(sampled[["fit"]], "prior_list")
  if (is.null(location_priors)) {
    location_priors <- sampled[["priors"]][["location"]]
  }
  sampled_vcov <- BayesTools::random_effects_marginal_vcov(
    fit               = formula_fit,
    parameter         = "mu",
    data              = new_location,
    posterior_samples = posterior_samples,
    prior_list        = location_priors,
    blocks            = "study",
    new_levels        = "sample"
  )[["samples"]]
  expected <- apply(sampled_vcov, c(2L, 3L), mean)

  set.seed(202)
  draws <- .predict_known_v_marginalized_random_draws(
    object            = marginalized,
    data              = new_data,
    posterior_samples = posterior_samples
  )

  expect_equal(unname(stats::cov(draws)), unname(expected), tolerance = .015)
})


test_that("conditional newdata random effects split known-R new-level policy by block", {

  posterior_samples <- matrix(
    0,
    nrow     = 2L,
    ncol     = 1L,
    dimnames = list(NULL, "mu")
  )
  fit <- coda::as.mcmc(posterior_samples)
  attr(fit, "prior_list") <- list()
  data <- list(
    outcome  = data.frame(yi = c(.10, .20, .30)),
    location = data.frame(study = c("s1", "s2", "s3"))
  )
  attr(data, "random") <- TRUE
  formula_design <- list(
    random_effects = list(
      list(
        block_name       = "known",
        compile_mode     = "sampled",
        group_covariance = structure(list(), class = "random_group_covariance")
      ),
      list(
        block_name   = "ordinary",
        compile_mode = "sampled"
      )
    )
  )
  object <- list(
    fit    = fit,
    data   = data,
    priors = list(location = list())
  )
  calls <- list()

  testthat::local_mocked_bindings(
    .fitted_formula_design = function(object, parameter, required = FALSE) {
      formula_design
    },
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    JAGS_predict_formula = function(fit, parameter, formula = NULL, data = NULL,
                                    prior_list = NULL,
                                    formula_target = "conditional",
                                    blocks = NULL, new_levels = NULL, ...) {

      calls[[length(calls) + 1L]] <<- list(
        blocks     = blocks,
        new_levels = new_levels
      )
      value <- if (identical(blocks, "known")) 1 else 2
      list(random = matrix(value, nrow = 3L, ncol = 2L))
    },
    .package = "BayesTools"
  )

  out <- .evaluate.brma.random_effects(
    fit               = fit,
    data              = data,
    priors            = object[["priors"]],
    posterior_samples = posterior_samples,
    same_data         = FALSE,
    formula_target    = "conditional",
    object            = object
  )

  expect_equal(
    vapply(calls, `[[`, character(1), "blocks"),
    c("known", "ordinary")
  )
  expect_equal(
    vapply(calls, `[[`, character(1), "new_levels"),
    c("error", "sample")
  )
  expect_equal(out, matrix(3, nrow = 2L, ncol = 3L))
})


test_that("known-R brma.mv prediction rejects new group levels", {

  name <- "brma.mv_block_mvn_known_R"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  newdata <- data.frame(
    x     = c(0, 1),
    z     = c(-1, 0),
    study = c("s1", "s2")
  )

  estimate <- predict(
    fit_brma,
    newdata = newdata,
    type    = "estimate",
    quiet   = TRUE
  )
  expect_brma_samples_matrix(
    estimate,
    nrow(newdata),
    "known-R brma.mv newdata estimate"
  )
  newdata[["study"]][2L] <- "s4"
  expect_error(
    predict(fit_brma, newdata = newdata, type = "estimate", quiet = TRUE),
    "R_new"
  )
})


test_that("marginal covariance random effects are returned as random draws", {

  covariance <- array(0, dim = c(2L, 2L, 2L))
  covariance[1L, , ] <- diag(c(.04, .09))
  covariance[2L, , ] <- diag(c(.01, .16))
  draws <- matrix(c(.1, .2, .3, .4), nrow = 2L, byrow = TRUE)
  captured <- NULL

  testthat::local_mocked_bindings(
    .outcome_rng.norm_known_v_covariance = function(mu_samples,
                                                    covariance_samples) {
      captured <<- list(mu = mu_samples, covariance = covariance_samples)
      draws
    },
    .package = "RoBMA"
  )

  out <- .random_effects_from_marginal_vcov(
    prediction = list(vcov = list(samples = covariance)),
    S          = 2L,
    K          = 2L
  )

  expect_equal(captured[["mu"]], matrix(0, nrow = 2L, ncol = 2L))
  expect_equal(captured[["covariance"]], covariance)
  expect_equal(out, t(draws))
})


test_that("random-formula brma.mv newdata response predictions use V_new covariance", {

  object <- .brma_mv_prior_object(random = TRUE)
  sd_names <- .brma_mv_random_sd_names(object)
  posterior_samples <- matrix(
    0,
    nrow = 2,
    ncol = 1 + length(sd_names),
    dimnames = list(NULL, c("mu", unname(sd_names)))
  )
  posterior_samples[, "mu"] <- c(0.1, 0.2)
  posterior_samples[, sd_names[["study"]]]        <- c(0.25, 0.30)
  posterior_samples[, sd_names[["effect_study"]]] <- c(0.05, 0.06)

  newdata <- data.frame(
    study  = c("s3", "s3"),
    effect = c("a", "c"),
    x      = c(0, 1)
  )
  V_new <- matrix(c(0.04, 0.012, 0.012, 0.05), nrow = 2)
  random_noise   <- matrix(c(.30, .40, .50, .60), nrow = 2L, byrow = TRUE)
  sampling_noise <- matrix(c(.01, .02, .03, .04), nrow = 2L, byrow = TRUE)
  captured       <- NULL

  testthat::local_mocked_bindings(
    .predict_brma_mv_new_effect_random_draws = function(object, data,
                                                        posterior_samples,
                                                        max_bytes = NULL) {
      random_noise
    },
    .known_v_sampling_noise = function(known_V, S, K) {
      captured <<- list(known_V = known_V, S = S, K = K)
      sampling_noise
    },
    .package = "RoBMA"
  )

  response <- predict(
    object,
    newdata            = newdata,
    V_new              = V_new,
    type               = "response",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )

  expect_brma_samples_matrix(response, nrow(newdata), "random brma.mv V_new response")
  expect_equal(
    unname(as.matrix(response)),
    matrix(c(.41, .52, .73, .84), nrow = 2L, byrow = TRUE),
    tolerance = 1e-12
  )
  expect_equal(.known_v_materialize(captured[["known_V"]]), V_new)
  expect_identical(captured[c("S", "K")], list(S = 2L, K = 2L))
  expect_equal(attr(response, "data")[["outcome"]][["sei"]], sqrt(diag(V_new)))
  response_target <- attr(response, "brma_mv_prediction_target")
  expect_equal(response_target[["type"]], "response")
  expect_true(response_target[["known_v"]])
  expect_true(response_target[["random_formula"]])
  expect_true(response_target[["v_new"]])
  expect_false(response_target[["random_effects_in_mean"]])
  expect_equal(response_target[["mean_target"]], "fixed_location")
  expect_equal(
    response_target[["covariance_target"]],
    "V_new_plus_marginal_random_effect_generation"
  )
  expect_equal(
    response_target[["new_levels"]],
    "marginal_new_effects"
  )
})

test_that("random-formula brma.mv V_new response matches generative covariance oracle", {

  object    <- .brma_mv_prior_object(random = TRUE)
  study_sd  <- 0.25
  effect_sd <- 0.05
  S         <- 20000L
  posterior_samples <- matrix(
    0,
    nrow     = S,
    ncol     = 1L,
    dimnames = list(NULL, "mu")
  )

  newdata <- data.frame(
    study  = c("s3", "s3"),
    effect = c("a", "b"),
    x      = c(0, 1)
  )
  V_new    <- matrix(c(0.04, 0.012, 0.012, 0.05), nrow = 2)
  new_data <- .prepare_newdata(
    object         = object,
    newdata        = .predict_known_v_newdata_add_vi(newdata, V_new),
    type           = "response",
    bias_adjusted  = FALSE,
    include_scale  = TRUE,
    include_random = TRUE
  )

  testthat::local_mocked_bindings(
    .predict_brma_mv_new_effect_random_draws = function(object, data,
                                                        posterior_samples,
                                                        max_bytes = NULL) {
      study  <- stats::rnorm(S, sd = study_sd)
      effect <- matrix(stats::rnorm(S * 2L, sd = effect_sd), nrow = S)
      effect + matrix(study, nrow = S, ncol = 2L)
    },
    .package = "RoBMA"
  )

  set.seed(82)
  response <- .predict_brma_mv_known_v_response_draws(
    object            = object,
    data              = new_data,
    known_V           = .known_v_newdata_prepare(V_new, k = nrow(newdata)),
    mu_samples        = matrix(0, nrow = S, ncol = 2L),
    posterior_samples = posterior_samples,
    same_data         = FALSE
  )

  expected <- V_new +
    matrix(study_sd^2, nrow = 2L, ncol = 2L) +
    diag(rep(effect_sd^2, 2L), nrow = 2L)
  expect_equal(unname(stats::cov(response)), unname(expected), tolerance = .01)
})


test_that("random-formula brma.mv V_new response evaluates row-varying marginalized scale", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  posterior_samples <- .get_posterior_samples(fit_brma[["fit"]])
  posterior_samples <- posterior_samples[seq_len(min(10L, nrow(posterior_samples))), , drop = FALSE]
  newdata <- data.frame(
    study  = c("s_new", "s_new"),
    effect = c("e1", "e2"),
    x      = c(0, 1)
  )
  V_new <- matrix(c(0.04, 0.010, 0.010, 0.05), nrow = 2)

  set.seed(1)
  response <- predict(
    fit_brma,
    newdata            = newdata,
    V_new              = V_new,
    type               = "response",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )

  expect_brma_samples_matrix(response, nrow(newdata), "random-scale brma.mv V_new response")
  expect_true(all(is.finite(as.matrix(response))))
  expect_equal(
    attr(response, "brma_mv_prediction_target")[["covariance_target"]],
    "V_new_plus_marginal_random_effect_generation"
  )
})

test_that("3-level component-scale brma.mv V_new response evaluates named SD sources", {

  model_names <- c(
    top    = "brma.mv_block_mvn_3lvl_scale_top",
    bottom = "brma.mv_block_mvn_3lvl_scale_bottom"
  )
  skip_if_missing_fits(unname(model_names))

  newdata <- data.frame(
    study  = c("s_new", "s_new"),
    effect = c("e1", "e2"),
    x      = c(0, 1)
  )
  V_new <- matrix(c(0.04, 0.010, 0.010, 0.05), nrow = 2)

  for (target in names(model_names)) {
    fit_brma <- fits[[model_names[[target]]]]
    posterior_samples <- .get_posterior_samples(fit_brma[["fit"]])
    posterior_samples <- posterior_samples[seq_len(min(10L, nrow(posterior_samples))), , drop = FALSE]

    set.seed(1)
    response <- predict(
      fit_brma,
      newdata            = newdata,
      V_new              = V_new,
      type               = "response",
      quiet              = TRUE,
      .posterior_samples = posterior_samples
    )
    response_target <- attr(response, "brma_mv_prediction_target")

    expect_brma_samples_matrix(response, nrow(newdata), paste(target, "V_new response"))
    expect_true(all(is.finite(as.matrix(response))), info = target)
    expect_equal(
      response_target[["covariance_target"]],
      "V_new_plus_marginal_random_effect_generation",
      info = target
    )
  }
})


test_that("known-R new-effect prediction chunks covariance draws to memory budget", {

  name <- "brma.mv_block_mvn_known_R"
  skip_if_missing_fits(name)

  object            <- fits[[name]]
  posterior_samples <- .get_posterior_samples(object[["fit"]])
  posterior_samples <- posterior_samples[
    seq_len(min(4L, nrow(posterior_samples))),
    ,
    drop = FALSE
  ]
  newdata <- data.frame(
    x     = c(0, 1),
    z     = c(-1, 0),
    study = c("s1", "s2")
  )
  new_data <- .prepare_newdata(
    object         = object,
    newdata        = newdata,
    type           = "estimate",
    bias_adjusted  = TRUE,
    include_scale  = TRUE,
    include_random = TRUE
  )
  original_evaluate <- .evaluate.brma.random_effects
  chunk_sizes       <- integer()

  testthat::local_mocked_bindings(
    .evaluate.brma.random_effects = function(fit, data, priors,
                                             posterior_samples = NULL, ...) {
      chunk_sizes <<- c(chunk_sizes, nrow(posterior_samples))
      original_evaluate(
        fit               = fit,
        data              = data,
        priors            = priors,
        posterior_samples = posterior_samples,
        ...
      )
    },
    .package = "RoBMA"
  )

  draws <- .predict_brma_mv_new_effect_random_draws(
    object            = object,
    data              = new_data,
    posterior_samples = posterior_samples,
    max_bytes         = .known_v_covariance_peak_bytes(1L, nrow(newdata))
  )

  expect_equal(dim(draws), c(nrow(posterior_samples), nrow(newdata)))
  expect_true(all(is.finite(draws)))
  expect_true(length(chunk_sizes) > 1L)
  expect_true(all(chunk_sizes <= 1L))
})

test_that("random-formula brma.mv estimate separates fixed and random effects", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  terms    <- as.matrix(predict(fit_brma, type = "terms", quiet = TRUE))
  estimate <- as.matrix(predict(fit_brma, type = "estimate", quiet = TRUE))
  random_components <- ranef(fit_brma)
  random            <- Reduce(`+`, lapply(random_components[["location"]], as.matrix))

  expect_equal(
    unname(estimate - terms),
    unname(random),
    tolerance = 1e-12
  )
})

test_that("random-formula brma.mv ranef decomposes random blocks", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  out      <- ranef(fit_brma)

  expect_type(out, "list")
  expect_named(out, "location")
  expect_true(all(c("effect_study", "study") %in% names(out[["location"]])))
  for (component in c("effect_study", "study")) {
    expect_brma_samples_matrix(
      out[["location"]][[component]],
      nobs(fit_brma),
      paste(name, component, "ranef")
    )
  }

  terms    <- as.matrix(predict(fit_brma, type = "terms", quiet = TRUE))
  estimate <- as.matrix(predict(fit_brma, type = "estimate", quiet = TRUE))
  total    <- Reduce(`+`, lapply(out[["location"]], as.matrix))
  total_component <- ranef(fit_brma, component = "total")
  study_component <- ranef(fit_brma, component = "study")
  location_component <- ranef(fit_brma, component = "location")

  expect_equal(
    unname(total),
    unname(estimate - terms),
    tolerance = 1e-10
  )
  expect_equal(
    unname(as.matrix(total_component)),
    unname(total),
    tolerance = 1e-10
  )
  expect_equal(
    unname(as.matrix(study_component)),
    unname(as.matrix(out[["location"]][["study"]])),
    tolerance = 1e-12
  )
  expect_equal(names(location_component), names(out[["location"]]))
})

test_that("marginalized brma.mv ranef uses Gaussian BLUP component", {

  dat <- data.frame(
    yi     = c(0.10, 0.20, 0.30, 0.40),
    effect = paste0("e", 1:4)
  )
  V <- matrix(
    c(
      0.04, 0.01, 0.00, 0.00,
      0.01, 0.05, 0.00, 0.00,
      0.00, 0.00, 0.04, 0.012,
      0.00, 0.00, 0.012, 0.05
    ),
    nrow = 4
  )
  object <- brma.mv(
    yi                        = yi,
    V                         = V,
    data                      = dat,
    random                    = ~ 1 | effect,
    known_v_parameterization  = "block_mvn",
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  posterior_samples <- matrix(
    c(
      0.05, 0.20,
      0.06, 0.30
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(NULL, c("mu", "mu__xREx__effect_intercept"))
  )

  terms    <- as.matrix(predict(
    object,
    type               = "terms",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  ))
  estimate <- as.matrix(predict(
    object,
    type               = "estimate",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  ))
  out <- ranef(
    object,
    .posterior_samples = posterior_samples
  )
  out_list <- ranef(
    object,
    simplify          = FALSE,
    .posterior_samples = posterior_samples
  )

  tau <- matrix(
    posterior_samples[, "mu__xREx__effect_intercept"],
    nrow = nrow(posterior_samples),
    ncol = nrow(dat)
  )
  expected <- .evaluate.brma.known_v_blup.norm(
    mu_samples = terms,
    tau_within = tau,
    yi         = dat[["yi"]],
    known_V    = .data_known_v_data(object[["data"]])
  ) - terms

  expect_brma_samples_matrix(out, nrow(dat), "marginalized ranef")
  expect_type(out_list, "list")
  expect_named(out_list, "location")
  expect_equal(names(out_list[["location"]]), "effect")
  expect_equal(
    unname(as.matrix(out)),
    unname(expected),
    tolerance = 1e-12
  )
  expect_equal(
    unname(as.matrix(out_list[["location"]][["effect"]])),
    unname(expected),
    tolerance = 1e-12
  )
  expect_equal(
    unname(estimate - terms),
    unname(expected),
    tolerance = 1e-12
  )
})


test_that("known-R row multipliers reach random prediction consumers", {

  dat <- data.frame(
    yi       = c(0.10, 0.20, 0.30),
    estimate = factor(c("e1", "e2", "e3"))
  )
  known_R <- diag(c(4, 9, 16))
  dimnames(known_R) <- list(levels(dat[["estimate"]]),
                            levels(dat[["estimate"]]))
  object <- brma.mv(
    yi                        = yi,
    V                         = diag(0.01, nrow(dat)),
    random                    = ~ 1 | estimate,
    R                         = known_R,
    Rscale                    = "none",
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  term <- .data_marginalized_random_effects(object[["data"]])[[1L]]
  sd_name <- term[["sd_parameter_names"]][[1L]]
  posterior_samples <- matrix(
    c(
      0, 0.20,
      0, 0.30
    ),
    nrow     = 2L,
    byrow    = TRUE,
    dimnames = list(NULL, c("mu", sd_name))
  )

  row_sd <- tcrossprod(
    posterior_samples[, sd_name],
    sqrt(c(4, 9, 16))
  )
  row_variance <- row_sd^2
  expected_blup <- row_variance / (row_variance + 0.01) *
    matrix(dat[["yi"]], nrow = 2L, ncol = 3L, byrow = TRUE)

  scale <- predict(
    object,
    type               = "terms.scale",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )
  pooled <- pooled_heterogeneity(
    object,
    .posterior_samples = posterior_samples
  )
  estimate <- predict(
    object,
    type               = "estimate",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )
  random_effect <- ranef(
    object,
    component          = "estimate",
    .posterior_samples = posterior_samples
  )

  expect_equal(unname(as.matrix(scale[["estimate"]])), row_sd,
               tolerance = 1e-12)
  expect_equal(
    unname(as.matrix(pooled)),
    matrix(sqrt(rowMeans(row_sd^2)), ncol = 1L),
    tolerance = 1e-12
  )
  expect_equal(unname(as.matrix(estimate)), expected_blup,
               tolerance = 1e-12)
  expect_equal(unname(as.matrix(random_effect)), expected_blup,
               tolerance = 1e-12)

  set.seed(1905)
  random_draws <- .predict_known_v_marginalized_random_draws(
    object            = object,
    data              = object[["data"]],
    posterior_samples = posterior_samples
  )
  set.seed(1905)
  expected_draws <- row_sd * matrix(
    stats::rnorm(length(row_sd)),
    nrow = nrow(row_sd),
    ncol = ncol(row_sd)
  )
  expect_equal(random_draws, expected_draws, tolerance = 1e-12)

  expect_error(
    predict(
      object,
      newdata            = data.frame(row = 1:2),
      type               = "terms.scale",
      quiet              = TRUE,
      .posterior_samples = posterior_samples
    ),
    "Known-R.*grouping"
  )
  mapped_scale <- predict(
    object,
    newdata = data.frame(
      estimate = factor(c("e1", "e3"), levels = levels(dat[["estimate"]]))
    ),
    type               = "terms.scale",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )
  expect_equal(
    unname(as.matrix(mapped_scale[["estimate"]])),
    row_sd[, c(1L, 3L), drop = FALSE],
    tolerance = 1e-12
  )
})


test_that("random-slope terms.scale replays the prediction design", {

  dat <- data.frame(
    yi    = c(0.10, 0.20, 0.30, 0.40),
    x     = c(0, 1, 0, 1),
    study = c("s1", "s1", "s2", "s2")
  )
  object <- brma.mv(
    yi                        = yi,
    V                         = diag(0.04, nrow(dat)),
    random                    = ~ diag(1 + x | study),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  term <- .fitted_formula_design(
    object,
    "mu",
    required = TRUE
  )[["random_effects"]][[1L]]
  posterior_samples <- matrix(
    c(
      0.20, 0.30,
      0.25, 0.35
    ),
    nrow     = 2L,
    byrow    = TRUE,
    dimnames = list(NULL, term[["sd_parameter_names"]])
  )
  new_x <- c(0, 2)
  new_z <- (new_x - mean(dat[["x"]])) / stats::sd(dat[["x"]])
  expected <- sqrt(
    posterior_samples[, 1L]^2 +
      tcrossprod(posterior_samples[, 2L]^2, new_z^2)
  )

  synthetic_groups <- predict(
    object,
    newdata            = data.frame(x = new_x),
    type               = "terms.scale",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )[["study"]]
  same_group <- predict(
    object,
    newdata            = data.frame(x = new_x, study = c("s1", "s1")),
    type               = "terms.scale",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )[["study"]]
  different_groups <- predict(
    object,
    newdata            = data.frame(x = new_x, study = c("s1", "s2")),
    type               = "terms.scale",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )[["study"]]

  expect_equal(unname(as.matrix(synthetic_groups)), expected,
               tolerance = 1e-12)
  expect_equal(unname(as.matrix(same_group)), expected,
               tolerance = 1e-12)
  expect_equal(unname(as.matrix(different_groups)), expected,
               tolerance = 1e-12)
  expect_error(
    predict(
      object,
      newdata            = data.frame(study = c("s1", "s2")),
      type               = "terms.scale",
      quiet              = TRUE,
      .posterior_samples = posterior_samples
    ),
    "random-effect.*x"
  )

  slope_only <- brma.mv(
    yi                        = yi,
    V                         = diag(0.04, nrow(dat)),
    random                    = ~ diag(0 + x | study),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  slope_term <- .fitted_formula_design(
    slope_only,
    "mu",
    required = TRUE
  )[["random_effects"]][[1L]]
  slope_samples <- matrix(
    c(0.30, 0.40),
    ncol     = 1L,
    dimnames = list(NULL, slope_term[["sd_parameter_names"]])
  )
  fitted_z <- (dat[["x"]] - mean(dat[["x"]])) / stats::sd(dat[["x"]])
  expected_fitted <- tcrossprod(slope_samples[, 1L], abs(fitted_z))
  fitted_scale <- predict(
    slope_only,
    type               = "terms.scale",
    quiet              = TRUE,
    .posterior_samples = slope_samples
  )[["study"]]

  expect_equal(unname(as.matrix(fitted_scale)), expected_fitted,
               tolerance = 1e-12)
})


test_that("scale-bound random slopes use row-effective prediction SDs", {

  dat <- data.frame(
    yi    = c(0.10, 0.20, 0.30, 0.40),
    x     = c(0, 1, 0, 1),
    w     = c(-1, -1, 1, 1),
    study = c("s1", "s1", "s2", "s2")
  )
  object <- brma.mv(
    yi                        = yi,
    V                         = diag(0.04, nrow(dat)),
    random                    = ~ diag(0 + x | study),
    scale                     = ~ w,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  posterior_samples <- matrix(
    c(0.20, log(2)),
    nrow     = 1L,
    dimnames = list(NULL, c("log_tau_intercept", "log_tau_w"))
  )
  newdata <- data.frame(
    x = c(0.5, 1.5),
    w = c(0, 1)
  )
  new_x <- (newdata[["x"]] - mean(dat[["x"]])) / stats::sd(dat[["x"]])
  expected <- matrix(
    0.20 * exp(log(2) * newdata[["w"]]) * abs(new_x),
    nrow = 1L
  )

  scale <- predict(
    object,
    newdata            = newdata,
    type               = "terms.scale",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )

  expect_equal(names(scale), "study")
  expect_equal(unname(as.matrix(scale[["study"]])), expected,
               tolerance = 1e-12)
  expect_error(
    predict(
      object,
      newdata            = data.frame(w = c(0, 1)),
      type               = "terms.scale",
      quiet              = TRUE,
      .posterior_samples = posterior_samples
    ),
    "random-effect.*x"
  )
})


test_that("scale-bound random factors replay their prediction design", {

  dat <- data.frame(
    yi    = seq(0.1, 0.6, by = 0.1),
    arm   = factor(rep(c("A", "B", "C"), 2L)),
    w     = rep(c(0, 1), each = 3L),
    study = rep(c("s1", "s2"), each = 3L)
  )
  object <- brma.mv(
    yi                             = yi,
    V                              = diag(0.04, nrow(dat)),
    random                         = ~ diag(1 + arm | study),
    scale                          = ~ w,
    data                           = dat,
    measure                        = "GEN",
    set_contrast_factor_predictors = "treatment",
    prior_unit_information_sd      = 1,
    only_priors                    = TRUE
  )
  term <- .fitted_formula_design(
    object,
    "mu",
    required = TRUE
  )[["random_effects"]][[1L]]
  allocation_name <- term[["sd_binding"]][["allocations"]][[1L]][[
    "weight_name"
  ]]
  posterior_samples <- matrix(
    c(0.25, 0, rep(1 / 3, 3L)),
    nrow = 1L,
    dimnames = list(
      NULL,
      c(
        "log_tau_intercept",
        "log_tau_w",
        paste0(allocation_name, "[", 1:3, "]")
      )
    )
  )
  newdata <- data.frame(
    arm = factor(c("A", "B", "C"), levels = levels(dat[["arm"]])),
    w   = 0
  )
  expected <- matrix(0.25 * c(1, sqrt(2), sqrt(2)), nrow = 1L)

  scale <- predict(
    object,
    newdata            = newdata,
    type               = "terms.scale",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )

  expect_equal(unname(as.matrix(scale[["study"]])), expected,
               tolerance = 1e-12)
})


test_that("terms.scale keeps sampled and scale-bound random blocks", {

  dat <- data.frame(
    yi   = c(0.10, 0.20, 0.30, 0.40, 0.15, 0.25),
    vi   = c(0.04, 0.05, 0.06, 0.07, 0.02, 0.03),
    type = factor(
      c("RCT", "RCT", "RCT", "cohort", "cohort", "cohort"),
      levels = c("RCT", "cohort")
    )
  )
  dat[["study"]]   <- factor(seq_len(nrow(dat)))
  dat[["cohort"]]  <- as.numeric(dat[["type"]] == "cohort")
  dat[["bias_id"]] <- factor("cohort_bias")
  object <- brma.mv(
    yi = yi,
    V  = vi,
    random = list(
      coh_bias    = ~ diag(0 + cohort | bias_id),
      ran_effects = ~ 1 | study
    ),
    scale = list(ran_effects = ~ type),
    data                           = dat,
    known_v_parameterization       = "block_mvn",
    measure                        = "GEN",
    prior_unit_information_sd      = 1,
    only_priors                    = TRUE
  )
  terms <- .fitted_formula_design(
    object,
    "mu",
    required = TRUE
  )[["random_effects"]]
  coh_bias_term <- terms[[which(vapply(
    terms,
    function(term) identical(term[["block_name"]], "coh_bias"),
    logical(1)
  ))]]
  posterior_samples <- matrix(
    c(0.40, 0.20, log(2)),
    nrow = 1L,
    dimnames = list(
      NULL,
      c(
        coh_bias_term[["sd_parameter_names"]],
        "log_tau_ran_effects_intercept",
        "log_tau_ran_effects_type"
      )
    )
  )
  cohort_z <- (dat[["cohort"]] - mean(dat[["cohort"]])) /
    stats::sd(dat[["cohort"]])
  expected_coh_bias <- matrix(0.40 * abs(cohort_z), nrow = 1L)
  expected_ran_effects <- matrix(
    c(rep(0.20, 3L), rep(0.40, 3L)),
    nrow = 1L
  )

  scale <- predict(
    object,
    type               = "terms.scale",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )

  expect_equal(names(scale), c("coh_bias", "ran_effects"))
  expect_equal(unname(as.matrix(scale[["coh_bias"]])), expected_coh_bias,
               tolerance = 1e-12)
  expect_equal(unname(as.matrix(scale[["ran_effects"]])),
               expected_ran_effects, tolerance = 1e-12)
})


test_that("random-formula brma.mv terms.scale returns random SD components", {

  object  <- .brma_mv_prior_object(random = TRUE)
  newdata <- data.frame(
    study  = c("s1", "s2"),
    effect = c("a", "b"),
    x      = c(0, 1)
  )
  sd_names          <- .brma_mv_random_sd_names(object)
  posterior_samples <- matrix(
    c(0.1, 0.2, 0.30, 0.40, 0.05, 0.06),
    nrow     = 2,
    dimnames = list(NULL, c("mu", unname(sd_names)))
  )

  scale <- predict(
    object,
    type               = "terms.scale",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )

  expect_type(scale, "list")
  expect_equal(names(scale), names(sd_names))
  for (component in names(scale)) {
    expect_brma_samples_matrix(scale[[component]], nobs(object), component)
    expect_equal(
      colnames(scale[[component]]),
      paste0("tau[", seq_len(nobs(object)), "]")
    )
  }

  expect_error(
    predict(object, newdata = newdata, type = "estimate", quiet = TRUE),
    "Posterior samples are required"
  )
})

test_that("known-V estimate predictions use the full covariance BLUP", {

  names <- c("brma.mv_latent", "brma.mv_whitened", "brma.mv_block_mvn")
  skip_if_missing_fits(names)

  for (name in names) {
    fit_brma <- fits[[name]]

    theta <- as.matrix(predict(fit_brma, type = "estimate", quiet = TRUE))
    mu    <- as.matrix(predict(fit_brma, type = "terms", quiet = TRUE))
    tau   <- as.matrix(predict(fit_brma, type = "terms.scale", quiet = TRUE))

    expected <- .evaluate.brma.known_v_blup.norm(
      mu_samples = mu,
      tau_within = tau,
      yi         = fit_brma[["data"]][["outcome"]][["yi"]],
      known_V    = .data_known_v_data(fit_brma[["data"]])
    )

    expect_equal(
      unname(theta),
      unname(expected),
      tolerance = 1e-12,
      info = paste(name, ": known-V estimate uses full covariance BLUP")
    )
  }
})
