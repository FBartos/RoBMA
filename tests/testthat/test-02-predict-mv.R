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
    aggregate_scale <- predict(
      fit_brma,
      newdata = TRUE,
      type    = "terms.scale",
      quiet   = TRUE
    )
    scale_fitted <- fitted(fit_brma, component = "scale")

    expect_equal(names(terms_scale), expected_components[[target]], info = target)
    expect_equal(names(new_scale), expected_components[[target]], info = target)
    expect_equal(names(aggregate_scale), expected_components[[target]], info = target)
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
      expect_brma_samples_matrix(
        aggregate_scale[[component]],
        1L,
        paste(target, component, "aggregate terms.scale")
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

test_that("random-formula brma.mv existing-data response uses conditional random effects in mean", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  posterior_samples <- .get_posterior_samples(fit_brma[["fit"]])
  posterior_samples <- posterior_samples[seq_len(min(25L, nrow(posterior_samples))), , drop = FALSE]
  captured <- NULL

  testthat::local_mocked_bindings(
    .outcome_rng.norm_known_v = function(mu_samples, tau_within, known_V) {
      captured <<- list(mu = mu_samples, tau = tau_within)
      mu_samples
    },
    .outcome_rng.norm_known_v_covariance = function(mu_samples, covariance_samples) {
      stop("Existing-data random-formula response should not use V_new covariance RNG.",
           call. = FALSE)
    },
    .package = "RoBMA"
  )

  expected <- as.matrix(predict(
    fit_brma,
    type               = "estimate",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  ))
  response <- as.matrix(predict(
    fit_brma,
    type               = "response",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  ))

  fixed_only <- as.matrix(predict(
    fit_brma,
    type               = "terms",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  ))

  expect_equal(unname(response), unname(expected), tolerance = 1e-12)
  expect_false(isTRUE(all.equal(unname(response), unname(fixed_only), tolerance = 1e-8)))
  expect_true(any(captured[["tau"]] > 0))
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

  aggregate_estimate <- predict(fit_brma, newdata = TRUE, type = "effect", quiet = TRUE)
  expect_brma_samples_matrix(
    aggregate_estimate,
    1,
    "random brma.mv aggregate effect"
  )

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
  aggregate <- predict(
    fit_brma,
    newdata = TRUE,
    type    = "effect",
    quiet   = TRUE
  )
  expect_brma_samples_matrix(
    aggregate,
    1,
    "known-R brma.mv aggregate effect"
  )

  newdata[["study"]][2L] <- "s4"
  expect_error(
    predict(fit_brma, newdata = newdata, type = "estimate", quiet = TRUE),
    "known group covariance"
  )
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
  captured <- NULL

  testthat::local_mocked_bindings(
    .outcome_rng.norm_known_v_covariance = function(mu_samples, covariance_samples) {
      captured <<- list(mu = mu_samples, covariance = covariance_samples)
      mu_samples
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
    matrix(c(0.1, 0.1, 0.2, 0.2), nrow = 2, byrow = TRUE),
    tolerance = 1e-12
  )
  expect_equal(unname(captured[["mu"]]), unname(as.matrix(response)), tolerance = 1e-12)
  expect_true(any(abs(captured[["covariance"]][, 1L, 2L] - V_new[1L, 2L]) > 0))
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
    "V_new_plus_marginal_random_effect_covariance"
  )
  expect_equal(
    response_target[["new_levels"]],
    "marginal_random_effect_covariance"
  )
})

test_that("random-formula brma.mv V_new response covariance matches ZGZ oracle", {

  object   <- .brma_mv_prior_object(random = TRUE)
  sd_names <- .brma_mv_random_sd_names(object)
  study_sd <- c(0.25, 0.30)
  effect_sd <- c(0.05, 0.06)
  posterior_samples <- matrix(
    0,
    nrow = 2,
    ncol = 1 + length(sd_names),
    dimnames = list(NULL, c("mu", unname(sd_names)))
  )
  posterior_samples[, "mu"] <- c(0.1, 0.2)
  posterior_samples[, sd_names[["study"]]]        <- study_sd
  posterior_samples[, sd_names[["effect_study"]]] <- effect_sd

  newdata <- data.frame(
    study  = c("s3", "s3"),
    effect = c("a", "b"),
    x      = c(0, 1)
  )
  V_new      <- matrix(c(0.04, 0.012, 0.012, 0.05), nrow = 2)
  known_V_new <- .known_v_newdata_prepare(V_new, k = nrow(newdata))
  new_data   <- .prepare_newdata(
    object         = object,
    newdata        = .predict_known_v_newdata_add_vi(newdata, V_new),
    type           = "response",
    bias_adjusted  = FALSE,
    include_scale  = TRUE,
    include_random = TRUE
  )

  covariance_samples <- .predict_known_v_newdata_response_covariance(
    object            = object,
    data              = new_data,
    known_V_new       = known_V_new,
    posterior_samples = posterior_samples
  )

  for (s in seq_len(nrow(posterior_samples))) {
    expected <- V_new +
      matrix(study_sd[s]^2, nrow = 2, ncol = 2) +
      diag(rep(effect_sd[s]^2, 2), nrow = 2, ncol = 2)
    expect_equal(
      unname(covariance_samples[s, , ]),
      unname(expected),
      tolerance = 1e-12
    )
  }
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
  captured <- NULL

  testthat::local_mocked_bindings(
    .outcome_rng.norm_known_v_covariance = function(mu_samples, covariance_samples) {
      captured <<- covariance_samples
      mu_samples
    },
    .package = "RoBMA"
  )

  response <- predict(
    fit_brma,
    newdata            = newdata,
    V_new              = V_new,
    type               = "response",
    quiet              = TRUE,
    .posterior_samples = posterior_samples
  )

  expect_brma_samples_matrix(response, nrow(newdata), "random-scale brma.mv V_new response")
  expect_equal(dim(captured), c(nrow(posterior_samples), nrow(newdata), nrow(newdata)))
  expect_true(all(captured[, 1L, 1L] > V_new[1L, 1L]))
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
      "V_new_plus_marginal_random_effect_covariance",
      info = target
    )
  }
})


test_that("random-formula brma.mv V_new response respects covariance memory budget", {

  object   <- .brma_mv_prior_object(random = TRUE)
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
    effect = c("a", "b"),
    x      = c(0, 1)
  )
  V_new <- matrix(c(0.04, 0.012, 0.012, 0.05), nrow = 2)
  old_options <- options(
    RoBMA.known_v_covariance_max_bytes = .known_v_covariance_bytes(1L, nrow(newdata))
  )
  on.exit(options(old_options), add = TRUE)

  expect_error(
    predict(
      object,
      newdata            = newdata,
      V_new              = V_new,
      type               = "response",
      quiet              = TRUE,
      .posterior_samples = posterior_samples
    ),
    "exceeding the configured budget"
  )
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
