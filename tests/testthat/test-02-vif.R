context("VIF")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-metafor.R"))

test_that("VIF rejects data-only objects", {

  dat <- data.frame(
    yi  = c(0.10, 0.20, 0.30, 0.40),
    sei = c(0.10, 0.12, 0.11, 0.15),
    x   = c(0, 1, 0, 1)
  )

  object <- brma.norm(
    yi        = yi,
    sei       = sei,
    mods      = ~ x,
    data      = dat,
    only_data = TRUE
  )

  expect_error(
    vif(object),
    "data-only"
  )
})

.vif_full_gls_vcov_oracle <- function(X, covariance_samples) {

  K <- nrow(X)
  P <- ncol(X)

  vcov_sum <- matrix(0, nrow = P, ncol = P)
  for (s in seq_len(dim(covariance_samples)[1L])) {
    M <- covariance_samples[s, , ]
    M <- (M + t(M)) / 2
    W <- solve(M)
    vcov_sum <- vcov_sum + solve(crossprod(X, W %*% X))
  }

  vcov <- vcov_sum / dim(covariance_samples)[1L]
  dimnames(vcov) <- list(colnames(X), colnames(X))

  return(vcov)
}

.vif_known_v_table_oracle <- function(object, covariance_samples) {

  X            <- .get_model_matrix(object)
  tables <- lapply(seq_len(dim(covariance_samples)[1L]), function(s) {
    covariance <- solve(crossprod(X, solve(covariance_samples[s, , ], X)))
    .vif_table_from_vcov(covariance, attr(X, "assign"), object, X)
  })
  result <- tables[[1L]]
  for(column in c("GVIF", "GVIF^(1/(2*df))")) {
    result[[column]] <- rowMeans(vapply(tables, `[[`, numeric(nrow(result)), column))
  }
  result
}

test_that("the live factor-plan VIF backend matches an independent dense oracle", {

  dat <- data.frame(yi = c(.1, .2, -.1, .3),
    x = c(-1, 0, 1, -1), z = c(0, 1, 1, 2))
  V <- matrix(c(.08, .02, 0, 0, .02, .09, 0, 0,
    0, 0, .07, .03, 0, 0, .03, .10), 4L)
  setup <- review_compile_object(brma.mv(yi = yi, V = V, mods = ~ x + z,
    data = dat, measure = "GEN", prior_unit_information_sd = 1,
    standardize_continuous_predictors = FALSE, only_priors = TRUE, silent = TRUE))
  samples <- matrix(0, 2L, 3L, dimnames = list(NULL, c("intercept", "x", "z")))
  colnames(samples)[1:3] <- vapply(c("intercept", "x", "z"),
    function(term) review_coefficient_name(setup, term), character(1))
  object <- review_attach_replay(setup, samples)
  BayesTools::JAGS_validate_fit_contract(object$fit, requires = c("name_encoding",
    "formula_name_map", "formula_design", "parameter_map", "draw_geometry"))
  X <- .get_model_matrix(object)
  covariance <- array(0, c(2L, 4L, 4L))
  for (j in 1:2) covariance[j, , ] <- V
  actual <- .vif_vcov_from_factor_plan(object, X, samples,
    .data_known_v_data(object$data))
  expect_equal(actual, .vif_full_gls_vcov_oracle(X, covariance), tolerance = 1e-12)
  diagonal <- covariance
  for (j in 1:2) diagonal[j, , ] <- diag(diag(covariance[j, , ]))
  expect_false(isTRUE(all.equal(actual, .vif_full_gls_vcov_oracle(X, diagonal),
    tolerance = 1e-8)))
})


test_that("the diagonal covariance helper retains its mean-covariance default", {

  X <- cbind(
    intercept = 1,
    x         = c(-2, -1, 0, 1, 2),
    z         = c(0, 1, 0, 2, 1)
  )
  vi          <- c(0.02, 0.05, 0.03, 0.08, 0.04)
  weights     <- rep(1, length(vi))
  tau_samples <- c(0, 1.2)
  tau_within  <- matrix(tau_samples, nrow = 2, ncol = length(vi))

  expected <- Reduce(`+`, lapply(tau_samples, function(tau) {
    W <- diag(weights / (vi + tau^2))
    solve(crossprod(X, W %*% X))
  })) / length(tau_samples)
  actual <- .vif_vcov_from_tau_samples(
    X                  = X,
    vi                 = vi,
    weights            = weights,
    tau_within_samples = tau_within
  )
  plug_in <- solve(crossprod(
    X,
    diag(weights / (vi + mean(tau_samples)^2)) %*% X
  ))

  dimnames(expected) <- dimnames(actual)
  dimnames(plug_in)  <- dimnames(actual)
  expect_equal(actual, expected, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(actual, plug_in, tolerance = 1e-8)))
})

test_that("diagonal and known-V covariance helpers agree before applying a functional", {

  X <- cbind(
    intercept = 1,
    x         = c(-2, -1, 0, 1, 2),
    z         = c(0, 1, 0, 2, 1)
  )
  vi          <- c(0.02, 0.05, 0.03, 0.08, 0.04)
  weights     <- rep(1, length(vi))
  tau_samples <- c(0, 0.15, 0.5)
  tau_within  <- matrix(tau_samples, nrow = 3, ncol = length(vi))
  covariance_samples <- array(0, dim = c(3, length(vi), length(vi)))
  for (s in seq_along(tau_samples)) {
    covariance_samples[s, , ] <- diag((vi + tau_samples[s]^2) / weights)
  }

  ordinary <- .vif_vcov_from_tau_samples(
    X                  = X,
    vi                 = vi,
    weights            = weights,
    tau_within_samples = tau_within
  )
  known_v <- .vif_full_gls_vcov_oracle(X, covariance_samples)

  expect_equal(ordinary, known_v, tolerance = 1e-12)
})

test_that("scalar heterogeneity draws remain compact for large designs", {

  S <- 15000L
  K <- 10000L
  tau <- matrix(seq(0, 1, length.out = S), ncol = 1)

  normalized <- .vif_tau_matrix(tau, K, "tau")

  expect_equal(dim(normalized), c(S, 1L))
  expect_lt(as.numeric(object.size(normalized)), 2 * S * 8)
})

test_that("compact cluster heterogeneity broadcasts across rows", {

  X <- cbind(
    intercept = 1,
    x         = c(-1, 0, 1, 2)
  )
  vi      <- c(0.02, 0.05, 0.03, 0.08)
  cluster <- c(1, 1, 2, 2)
  within  <- matrix(c(0.1, 0.2), ncol = 1)
  between <- matrix(c(0.3, 0.4), ncol = 1)

  compact <- .vif_vcov_from_tau_samples(
    X                   = X,
    vi                  = vi,
    weights             = rep(1, length(vi)),
    tau_within_samples  = within,
    tau_between_samples = between,
    cluster             = cluster
  )
  expanded <- .vif_vcov_from_tau_samples(
    X                   = X,
    vi                  = vi,
    weights             = rep(1, length(vi)),
    tau_within_samples  = matrix(within[, 1], nrow = 2, ncol = length(vi)),
    tau_between_samples = matrix(between[, 1], nrow = 2, ncol = length(vi)),
    cluster             = cluster
  )

  expect_equal(compact, expanded, tolerance = 1e-12)
})

test_that("known-V VIF uses sampling covariance alone without scalar tau", {

  object <- structure(
    list(
      summary = matrix(
        c(0.1, 0.2),
        nrow     = 2,
        dimnames = list(c("(mu) intercept", "(mu) x"), "Mean")
      ),
      data = structure(list(), scale = FALSE)
    ),
    class = c("brma.mv", "brma.norm", "brma")
  )
  V <- matrix(c(0.04, 0.01, 0.01, 0.05), nrow = 2)

  posterior_samples <- matrix(
    c(0.1, 0.2, 0.3, 0.4),
    nrow = 2,
    dimnames = list(NULL, c("mu", "x"))
  )
  known_V <- .known_v_newdata_prepare(V, k = nrow(V))
  covariance_samples <- .known_v_add_base_covariance(
    known_V            = known_V,
    covariance_samples = .known_v_diagonal_extra_covariance_samples(
      object            = object,
      posterior_samples = posterior_samples,
      K                 = 2
    )
  )

  expect_equal(covariance_samples[1, , ], V)
})

test_that("fixed location coefficient extraction requires declared location columns", {

  samples <- matrix(seq_len(6), nrow = 2,
    dimnames = list(NULL, c("(mu) intercept", "(mu) x", "(tau) intercept")))
  filtered <- .diagnostic_fixed_location_coefficient_samples(samples)
  expect_equal(colnames(filtered), c("(mu) intercept", "(mu) x"))
  for (invalid in list(samples[, 3, drop = FALSE], matrix(seq_len(4), 2L))) {
    expect_error(.diagnostic_fixed_location_coefficient_samples(invalid),
      "Fixed location coefficient columns could not be identified.", fixed = TRUE)
  }
})

.vif_draw_replay <- function(kind = "simple") {

  data <- data.frame(yi = c(.2, -.1, .4, .1, .6, -.2),
    sei = c(.1, .3, .2, .4, .15, .5), x = c(-2, -1, 0, 1, 2, .5),
    z = c(0, 1, 0, 2, 1, -.5), cluster = rep(1:3, each = 2))
  if (kind == "mv_diagonal") data$cluster <- seq_len(nrow(data))
  args <- list(yi = quote(yi), mods = ~ x + z, data = data, measure = "GEN",
    prior_unit_information_sd = 1, standardize_continuous_predictors = FALSE,
    only_priors = TRUE, silent = TRUE)
  if (kind == "factor") {
    data$f <- factor(rep(letters[1:3], each = 2L))
    args$data <- data
    args$mods <- ~ x + f
  }
  if (kind == "ensemble") {
    args[c("prior_effect_null", "prior_heterogeneity_null", "prior_mods_null")] <- list(NULL, NULL, NULL)
    args$prior_bias <- list(
      BayesTools::prior_PET("normal", list(0, 1), list(0, Inf)),
      BayesTools::prior_PEESE("normal", list(0, 1), list(0, Inf)))
  }
  if (kind %in% c("mv", "mv_diagonal")) {
    args$V <- diag(data$sei^2)
    args$random <- ~ 1 | cluster
  } else {
    args$sei <- quote(sei)
  }
  if (kind == "scale") args$scale <- ~ x
  if (kind == "cluster") args$cluster <- quote(cluster)
  constructor <- switch(kind, mv = brma.mv, mv_diagonal = brma.mv, ensemble = RoBMA, brma)
  setup <- review_compile_object(do.call(constructor, args))
  priors <- list()
  for (compiled in setup$compiled) priors <- c(priors, compiled$prior_list)
  priors <- c(priors, .create_fit_priors(setup$object$data, setup$object$priors))
  columns <- unique(unlist(lapply(names(priors), function(name) {
    BayesTools:::.prior_linear_prior_columns(name, priors[[name]])
  })))
  if (kind == "ensemble") columns <- unique(c(columns, "bias_indicator", "PET", "PEESE"))
  tau <- rep(c(.05, .2, .8, 1.2, 2), 3L)
  slope <- rep(c(-.3, .2, .1, -.1, .4), 3L)
  samples <- matrix(0, length(tau), length(columns), dimnames = list(NULL, columns))
  samples[, grepl("tau|xREx.*intercept$", columns)] <- tau
  if ("log_tau_x" %in% columns) samples[, "log_tau_x"] <- slope
  if ("rho" %in% columns) samples[, "rho"] <- .4
  if (kind == "ensemble") {
    samples[, grepl("indicator", columns)] <- 1
    samples[, "bias_indicator"] <- match(rep(1:3, length.out = length(tau)),
      .bias_design_member_classes(setup$object$priors))
  }
  object <- review_attach_replay(setup, samples)
  class(object) <- setdiff(class(object), "only_priors.brma")
  covariance <- array(0, c(length(tau), nrow(data), nrow(data)))
  for (s in seq_along(tau)) {
    within <- if (kind == "scale") tau[s] * exp(slope[s] * data$x) else tau[s]
    covariance[s, , ] <- diag(data$sei^2 + switch(kind,
      cluster = .6 * within^2, mv = 0, mv_diagonal = 0, within^2), nrow(data))
    if (kind %in% c("cluster", "mv", "mv_diagonal")) for (rows in split(seq_len(nrow(data)), data$cluster)) {
      covariance[s, rows, rows] <- covariance[s, rows, rows] +
        tau[s]^2 * if (kind == "cluster") .4 else 1
    }
  }
  list(object = object, covariance = covariance, data = data)
}

test_that("public VIF summarizes nonlinear draw-wise diagnostics in every covariance class", {

  for (kind in c("simple", "scale", "cluster", "mv")) {
    setup <- .vif_draw_replay(kind)
    object <- setup$object
    X <- .get_model_matrix(object)
    expected <- t(vapply(seq_len(dim(setup$covariance)[1L]), function(s) {
      covariance <- solve(crossprod(X, solve(setup$covariance[s, , ], X)))
      diag(solve(stats::cov2cor(covariance)[-1L, -1L]))
    }, numeric(2L)))
    actual <- vif(object, posterior_correlation = FALSE,
                  probs = c(.025, .975), include_samples = TRUE)
    expect_equal(actual$vif$GVIF, unname(colMeans(expected)), tolerance = 1e-11, info = kind)
    expect_equal(actual$vif[["GVIF^(1/(2*df))"]], unname(colMeans(sqrt(expected))),
                 tolerance = 1e-11, info = kind)
    expect_equal(actual$vif$CI_0.025, unname(apply(expected, 2, stats::quantile, .025)),
                 tolerance = 1e-11, info = kind)
    expect_equal(actual$vif$CI_0.975, unname(apply(expected, 2, stats::quantile, .975)),
                 tolerance = 1e-11, info = kind)
    expect_equal(actual$samples$GVIF, as.vector(expected), tolerance = 1e-11, info = kind)
    expect_equal(actual$samples$draw, rep(1:15, 2L), info = kind)
    expect_equal(as.data.frame(actual), data.frame(actual, check.names = FALSE))
    expect_true(all(c("CI_0.025", "CI_0.975") %in% names(data.frame(actual))))
    expect_output(print(actual), "CI_0.025", fixed = TRUE)
    averaged <- .vif_known_v_table_oracle(object, setup$covariance)
    expect_equal(actual$vif$GVIF, averaged$GVIF, tolerance = 1e-11)
    old <- .vif_table_from_vcov(.vif_full_gls_vcov_oracle(X, setup$covariance),
                               attr(X, "assign"), object, X)
    expect_gt(max(abs(actual$vif$GVIF - old$GVIF)), 1e-5, label = kind)
    thinned <- expect_warning(vif(object, posterior_correlation = FALSE,
      max_samples = 10, include_samples = TRUE), "deterministically thinned")
    rows <- .thin_sample_rows(15, 10)
    expect_equal(thinned$samples$draw, rep(rows, 2L))
    expect_equal(thinned$vif$GVIF, unname(colMeans(expected[rows, , drop = FALSE])), tolerance = 1e-11)
  }
})

test_that("posterior VIF options validate probabilities without hiding numerical missingness", {

  object <- .vif_draw_replay()$object
  expect_error(vif(object, probs = c(.5, .5)), "'probs' must not contain duplicate probabilities.", fixed = TRUE)
  expect_error(vif(object, probs = NA_real_), "probs")
  expect_error(vif(object, probs = 1.1), "probs")
  expect_error(vif(object, include_samples = NA), "include_samples")
  template <- data.frame(term = "x", df = 1, GVIF = 0, "GVIF^(1/(2*df))" = 0, check.names = FALSE)
  summary <- .vif_summarize_draws(matrix(c(1, NaN), ncol = 1), template, 1:2, probs = .5)
  expect_true(is.nan(summary$GVIF))
  expect_true(is.na(summary$CI_0.5))
})

test_that("posterior VIF is invariant to equivalent simple and known-V representations", {

  ordinary <- vif(.vif_draw_replay("simple")$object, posterior_correlation = FALSE,
    probs = c(.025, .975), include_samples = TRUE)
  known_v <- vif(.vif_draw_replay("mv_diagonal")$object, posterior_correlation = FALSE,
    probs = c(.025, .975), include_samples = TRUE)
  expect_equal(ordinary$vif, known_v$vif, tolerance = 1e-11)
  expect_equal(ordinary$samples, known_v$samples, tolerance = 1e-11)
})



test_that("posterior GVIF preserves factor groups and interval meanings", {

  setup <- .vif_draw_replay("factor")
  object <- setup$object
  X <- .get_model_matrix(object)
  expected <- t(vapply(seq_len(dim(setup$covariance)[1L]), function(s) {
    covariance <- solve(crossprod(X, solve(setup$covariance[s, , ], X)))
    R <- stats::cov2cor(covariance)[-1L, -1L]
    c(det(R[2:3, 2:3]) / det(R), det(R[2:3, 2:3]) / det(R))
  }, numeric(2L)))
  result <- vif(object, posterior_correlation = FALSE, probs = c(.1, .9), include_samples = TRUE)
  expect_identical(result$vif$df, c(1L, 2L))
  expect_equal(result$vif$GVIF, colMeans(expected), tolerance = 1e-11)
  expect_equal(result$vif[["GVIF^(1/(2*df))"]],
    c(mean(sqrt(expected[, 1])), mean(expected[, 2]^.25)), tolerance = 1e-11)
  expect_equal(result$vif$CI_0.1, unname(apply(expected, 2, stats::quantile, .1)), tolerance = 1e-11)
  expect_equal(result$vif$CI_0.9, unname(apply(expected, 2, stats::quantile, .9)), tolerance = 1e-11)
  expect_equal(as.data.frame(result), data.frame(result, check.names = FALSE))
  expect_true(all(vapply(result$samples[c("draw", "df", "GVIF", "GVIF^(1/(2*df))")], is.numeric, logical(1))))
})

test_that("posterior VIF aligns model-specific draws and conditions bias intervals", {

  setup <- .vif_draw_replay("ensemble")
  object <- setup$object
  data <- setup$data
  samples <- .get_posterior_samples(object$fit)
  classes <- .bias_design_member_classes(object$priors)[samples[, "bias_indicator"]]
  selected <- .thin_sample_rows(nrow(samples), 10)
  expected <- matrix(NA_real_, length(selected), 4L)
  for(i in seq_along(selected)){
    row <- selected[i]
    X <- cbind(1, data$x, data$z)
    if(classes[row] == 2L) X <- cbind(X, data$sei)
    if(classes[row] == 3L) X <- cbind(X, data$sei^2)
    covariance <- solve(crossprod(X, solve(setup$covariance[row, , ], X)))
    values <- diag(solve(stats::cov2cor(covariance)[-1L, -1L]))
    columns <- c(1:2, if(classes[row] != 1L) classes[row] + 1L)
    expected[i, columns] <- values
  }
  result <- expect_warning(vif(object, posterior_correlation = FALSE, max_samples = 10,
    probs = c(.25, .75), include_samples = TRUE), "deterministically thinned")
  expect_identical(result$vif$term, c("x", "z", "PET", "PEESE"))
  expect_equal(result$samples$draw, rep(selected, 4L))
  expect_equal(result$samples$GVIF, as.vector(expected), tolerance = 1e-10)
  expect_equal(result$vif$GVIF, unname(colMeans(expected, na.rm = TRUE)), tolerance = 1e-10)
  expect_equal(result$vif[["GVIF^(1/(2*df))"]], unname(colMeans(sqrt(expected), na.rm = TRUE)), tolerance = 1e-10)
  expect_equal(result$vif$CI_0.25, unname(apply(expected, 2, stats::quantile, .25, na.rm = TRUE)), tolerance = 1e-10)
  expect_equal(result$vif$CI_0.75, unname(apply(expected, 2, stats::quantile, .75, na.rm = TRUE)), tolerance = 1e-10)
  samples[classes == 3L, "bias_indicator"] <- match(1L, .bias_design_member_classes(object$priors))
  object$fit <- BayesTools::JAGS_with_draws(object$fit, samples)
  absent <- vif(object, posterior_correlation = FALSE, probs = .5, include_samples = TRUE)
  expect_true(all(is.na(absent$samples$GVIF[absent$samples$term == "PEESE"])))
  expect_true(all(is.na(as.matrix(absent$vif[absent$vif$term == "PEESE", c("GVIF", "GVIF^(1/(2*df))", "CI_0.5")]))))
  expect_equal(nrow(absent$samples), nrow(samples) * 4L)
})

test_that("posterior VIF sample rows distinguish colliding display labels", {

  template <- data.frame(term = c("PET", "PET"), df = 1L,
    GVIF = NA_real_, "GVIF^(1/(2*df))" = NA_real_, check.names = FALSE)
  values <- matrix(c(1, 4, NA, 2), nrow = 2L)
  included <- matrix(c(TRUE, TRUE, FALSE, TRUE), nrow = 2L)
  result <- .vif_summarize_draws(values, template, c(2L, 9L),
    probs = .5, include_samples = TRUE, included = included)
  samples <- attr(result, "posterior_samples")
  expect_identical(samples$term_index, c(1L, 1L, 2L, 2L))
  expect_identical(samples$draw, c(2L, 9L, 2L, 9L))
  expect_identical(samples$term, rep("PET", 4L))
  expect_equal(result$GVIF, c(2.5, 2))
  expect_equal(samples$GVIF[samples$term_index == 2L], c(NA_real_, 2))
})

skip_if_no_fits()
skip_if_not_installed("metafor")

fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

for_each_case(vif_cases(), function(case) {
  test_that_case("VIF matches metafor", case, {
    expect_vif_matches_metafor(case)
  })
})

test_that("VIF returns optional posterior correlation diagnostics", {

  name <- "bcg_meta-regression"
  skip_if_missing_fits(name)

  fit_brma      <- fits[[name]]
  vif_with_post <- vif(fit_brma)
  post_cor      <- vif_with_post[["posterior_correlation"]]

  expect_s3_class(vif_with_post, "vif.brma")
  expect_vif_table(vif_with_post[["vif"]], 2)
  expect_type(post_cor, "double")
  expect_equal(nrow(post_cor), ncol(post_cor))
  expect_equal(unname(diag(post_cor)), rep(1, nrow(post_cor)), tolerance = 1e-12)
  expect_equal(post_cor, t(post_cor), tolerance = 1e-12)
  expect_true(all(is.finite(post_cor)))

  vif_without_post <- vif(fit_brma, posterior_correlation = FALSE)
  expect_null(vif_without_post[["posterior_correlation"]])
})

test_that("VIF rejects models without moderators", {

  model_names <- c(
    "bcg_meta-analysis",
    "konstantopoulos2011_3lvl",
    "nielweise2008_glmm",
    "dat.lehmann2018-PET",
    "dat.lehmann2018-3PSM",
    "brma.mv_block_mvn",
    "brma.mv_block_mvn_random_scale"
  )
  skip_if_missing_fits(model_names)

  for (name in model_names) {
    expect_error(
      vif(fits[[name]]),
      "only meaningful for models with moderators",
      info = paste("vif rejects model without moderators for", name)
    )
  }
})

test_that("ordinary covariance helper averages posterior GLS covariances", {

  name <- "bcg_meta-regression"
  skip_if_missing_fits(name)

  object <- fits[[name]]
  X      <- .get_model_matrix(object)
  K      <- nrow(X)
  tau <- .evaluate.brma.tau(
    fit               = object[["fit"]],
    scale_data        = object[["data"]][["scale"]],
    scale_formula     = NULL,
    scale_priors      = object[["priors"]][["scale"]],
    is_scale          = FALSE,
    is_multilevel     = FALSE,
    K                 = K,
    fixed_tau         = .fixed_tau_prior_value(object[["priors"]]),
    fixed_rho         = .fixed_rho_prior_value(object[["priors"]])
  )
  vi <- .outcome_data_vi(object)
  weights <- .outcome_data_weights(object)
  if (is.null(weights)) weights <- rep(1, K)
  expected <- Reduce(`+`, lapply(seq_len(nrow(tau[["tau_within"]])), function(s) {
    precision <- weights / (vi + tau[["tau_within"]][s, ]^2)
    solve(crossprod(X, X * precision))
  })) / nrow(tau[["tau_within"]])
  dimnames(expected) <- list(colnames(X), colnames(X))

  expect_equal(.vif_vcov_brma(object, X), expected, tolerance = 1e-12)
})

test_that("product-space covariance helpers average posterior covariance", {

  names <- c(
    "dat.lehmann2018_BMA.norm_mods",
    "dat.lehmann2018_RoBMA_mods",
    "bcg_BMA.glmm_3lvl_location_scale"
  )
  skip_if_missing_fits(names)

  for (name in names) {
    object        <- fits[[name]]
    X             <- .get_model_matrix(object)
    vi            <- .outcome_data_vi(object)
    weights       <- .outcome_data_weights(object)
    is_scale      <- .is_scale(object)
    is_multilevel <- .is_multilevel(object)
    tau <- .evaluate.brma.tau(
      fit               = object[["fit"]],
      scale_data        = object[["data"]][["scale"]],
      scale_formula     = if (is_scale) .create_fit_formula_list(
        data = object[["data"]],
        "scale"
      ) else NULL,
      scale_priors      = object[["priors"]][["scale"]],
      is_scale          = is_scale,
      is_multilevel     = is_multilevel,
      K                 = nrow(X),
      fixed_tau         = .fixed_tau_prior_value(object[["priors"]]),
      fixed_rho         = .fixed_rho_prior_value(object[["priors"]])
    )
    block_indices <- if (is_multilevel) {
      split(
        seq_len(nrow(X)),
        object[["data"]][["outcome"]][["cluster"]]
      )
    } else {
      list()
    }
    tau_within <- tau[["tau_within"]]
    expected <- Reduce(`+`, lapply(seq_len(nrow(tau_within)), function(s) {
      M <- diag(
        (vi + tau_within[s, ]^2) / weights,
        nrow = length(vi)
      )
      for (index in block_indices) {
        M[index, index] <- M[index, index] + tcrossprod(
          tau[["tau_between"]][s, index]
        )
      }
      solve(crossprod(X, solve(M, X)))
    })) / nrow(tau_within)
    dimnames(expected) <- list(colnames(X), colnames(X))

    expect_equal(
      .vif_vcov_brma(object, X),
      expected,
      tolerance = 1e-10,
      info      = name
    )
  }
})

test_that("fixed heterogeneity VIF matches metafor covariance", {

  X <- cbind(
    intercept = 1,
    x         = c(-2, -1, 0, 1, 2),
    z         = c(0, 1, 0, 2, 1)
  )
  vi      <- c(0.02, 0.05, 0.03, 0.08, 0.04)
  tau     <- 0.25
  weights <- rep(1, length(vi))
  actual <- .vif_vcov_from_tau_samples(
    X                  = X,
    vi                 = vi,
    weights            = weights,
    tau_within_samples = matrix(tau, nrow = 8, ncol = length(vi))
  )
  fit_metafor <- metafor::rma(
    yi        = c(-0.2, 0.1, 0.3, 0.05, 0.4),
    vi        = vi,
    mods      = X[, -1, drop = FALSE],
    intercept = TRUE,
    tau2      = tau^2
  )
  expected <- fit_metafor[["vb"]]
  dimnames(expected) <- dimnames(actual)

  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("equivalent brma and diagonal-V brma.mv fits have VIF parity", {

  names <- c("vif_parity_brma", "vif_parity_brma_mv")
  skip_if_missing_fits(names)

  ordinary <- vif(
    fits[[names[[1]]]],
    posterior_correlation = FALSE
  )[["vif"]]
  known_v <- vif(
    fits[[names[[2]]]],
    posterior_correlation = FALSE
  )[["vif"]]

  expect_equal(ordinary[["term"]], known_v[["term"]])
  expect_equal(ordinary[["df"]], known_v[["df"]])
  expect_equal(ordinary[["GVIF"]], known_v[["GVIF"]], tolerance = 0.05)
  expect_equal(
    ordinary[["GVIF^(1/(2*df))"]],
    known_v[["GVIF^(1/(2*df))"]],
    tolerance = 0.03
  )
})

test_that("VIF supports brma.mv known-V GLS meta-regressions", {

  name <- "brma.mv_block_mvn_mods"
  skip_if_missing_fits(name)

  object  <- fits[[name]]
  covariance <- .known_v_marginal_covariance_samples(object)
  expected   <- .vif_known_v_table_oracle(
    object             = object,
    covariance_samples = covariance
  )
  actual <- vif(object, posterior_correlation = FALSE)[["vif"]]

  expect_vif_table(actual, 2, info = name)
  expect_equal(actual, expected, tolerance = 1e-12)
})

test_that("VIF supports brma.mv random-formula marginal GLS covariance", {

  name <- "brma.mv_block_mvn_random_mods_scale"
  skip_if_missing_fits(name)

  object      <- fits[[name]]
  known_V     <- .data_known_v_data(object[["data"]])
  # The scale formula supplies row-specific external SDs through RoBMA's
  # adapter; these derived values are not monitored columns in the raw fit.
  random_vcov <- .brma_mv_random_effects_marginal_vcov(
    object            = object,
    posterior_samples = .get_posterior_samples(object[["fit"]])
  )
  covariance_samples <- .known_v_marginal_covariance_samples(object)

  expected <- .vif_known_v_table_oracle(
    object             = object,
    covariance_samples = covariance_samples
  )
  result <- vif(object)
  actual <- result[["vif"]]
  post_cor <- result[["posterior_correlation"]]

  expect_vif_table(actual, 2, info = name)
  expect_equal(actual, expected, tolerance = 1e-10)
  expect_equal(colnames(post_cor), c("(mu) intercept", "(mu) x", "(mu) z"))
  expect_false(any(grepl("tau", colnames(post_cor), fixed = TRUE)))
  expect_false(any(grepl("tau2_prop(", colnames(post_cor), fixed = TRUE)))
  expect_equal(
    covariance_samples,
    .known_v_add_base_covariance(
      known_V            = known_V,
      covariance_samples = random_vcov[["samples"]]
    )
  )
  expect_equal(random_vcov[["metadata"]][["n_rows"]], nobs(object))
  expect_true(length(random_vcov[["metadata"]][["included_blocks"]]) > 0)
  expect_equal(nrow(random_vcov[["metadata"]][["skipped_blocks"]]), 0)
  expect_true(any(abs(random_vcov[["samples"]]) > 0))
})

test_that("VIF supports brma.mv known-R marginal GLS covariance", {

  name <- "brma.mv_block_mvn_known_R"
  skip_if_missing_fits(name)

  object      <- fits[[name]]
  known_V     <- .data_known_v_data(object[["data"]])
  random_vcov <- BayesTools::random_effects_marginal_vcov(
    fit       = object[["fit"]],
    parameter = "mu"
  )
  covariance_samples <- .known_v_marginal_covariance_samples(object)

  expected <- .vif_known_v_table_oracle(
    object             = object,
    covariance_samples = covariance_samples
  )
  actual <- vif(object, posterior_correlation = FALSE)[["vif"]]
  group_covariance <- random_vcov[["metadata"]][["blocks"]][["study"]][["group_covariance"]]

  expect_vif_table(actual, 2, info = name)
  expect_equal(actual, expected, tolerance = 1e-10)
  expect_equal(
    covariance_samples,
    .known_v_add_base_covariance(
      known_V            = known_V,
      covariance_samples = random_vcov[["samples"]]
    )
  )
  expect_equal(random_vcov[["metadata"]][["included_blocks"]], "study")
  expect_equal(group_covariance[["scale"]], "none")
  expect_equal(unname(diag(group_covariance[["kernel"]])), c(4, 9, 16))
  expect_true(any(abs(random_vcov[["samples"]]) > 0))
})

test_that("VIF avoids dense known-V covariance allocation", {

  name <- "brma.mv_block_mvn_random_mods_scale"
  skip_if_missing_fits(name)

  object       <- fits[[name]]
  expected     <- vif(object, posterior_correlation = FALSE)[["vif"]]
  old_options  <- options(
    RoBMA.known_v_covariance_max_bytes = 1
  )
  on.exit(options(old_options), add = TRUE)

  result   <- vif(object, posterior_correlation = FALSE)
  actual   <- result[["vif"]]
  expect_equal(actual, expected, tolerance = 1e-10)
  expect_null(attr(result, "known_v_diagnostic"))
})

test_that("VIF exposes deterministic known-V draw thinning", {

  name <- "brma.mv_block_mvn_random_mods_scale"
  skip_if_missing_fits(name)

  object <- fits[[name]]

  result <- expect_warning(
    vif(object, posterior_correlation = FALSE, max_samples = 10),
    "deterministically thinned"
  )
  metadata <- attr(result, "known_v_diagnostic")

  expect_vif_table(result[["vif"]], 2, info = name)
  expect_equal(metadata[["n_used_samples"]], 10L)
  expect_true(metadata[["thinned"]])
})

test_that("Full known-V covariance helper respects allocation budget", {

  name <- "brma.mv_block_mvn_random_mods_scale"
  skip_if_missing_fits(name)

  object            <- fits[[name]]
  posterior_samples <- .get_posterior_samples(object[["fit"]])
  K                 <- nobs(object)
  one_draw_mem      <- .known_v_covariance_bytes(1L, K)

  expect_error(
    .known_v_marginal_covariance_samples(
      object            = object,
      posterior_samples = posterior_samples[seq_len(2L), , drop = FALSE],
      max_bytes         = one_draw_mem
    ),
    "exceeding the configured budget"
  )
})
