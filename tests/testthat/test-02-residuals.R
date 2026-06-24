context("Residuals")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))
source(testthat::test_path("helper-test-matrix.R"))
source(testthat::test_path("helper-metafor.R"))

REFERENCE_DIR <<- testthat::test_path("..", "results", "residuals")

skip_if_no_fits()
skip_if_not_installed("metafor")

fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)

for_each_case(residual_metafor_cases(), function(case) {
  test_that_case("Residuals match metafor", case, {
    expect_residuals_match_metafor(case)
  })
})

test_that("Residuals for equivalent interaction parameterizations match", {

  model_names <- c("bcg_meta-regression3", "bcg_meta-regression3b")
  skip_if_missing_fits(model_names)

  fit_brma1 <- fits[["bcg_meta-regression3"]]
  fit_brma2 <- fits[["bcg_meta-regression3b"]]

  expect_equal(residuals(fit_brma1), residuals(fit_brma2), tolerance = 0.10,
               info = "outcome residuals match")
  expect_equal(
    residuals(fit_brma1, type = "pearson"),
    residuals(fit_brma2, type = "pearson"),
    tolerance = 0.10,
    info      = "pearson residuals match"
  )
  expect_equal(rstandard(fit_brma1), rstandard(fit_brma2), tolerance = 0.10,
               info = "rstandard residuals match")
})

test_that("Extended rstudent residuals for interaction parameterizations match", {

  skip_if_not_full_diagnostics(
    "Interaction parameterization LOO-PIT checks are redundant with core residual coverage."
  )
  model_names <- c("bcg_meta-regression4", "bcg_meta-regression4b")
  skip_if_missing_fits(model_names)

  brma_rstudent1 <- suppressWarnings(rstudent(fits[["bcg_meta-regression4"]]))
  brma_rstudent2 <- suppressWarnings(rstudent(fits[["bcg_meta-regression4b"]]))

  expect_true(cor(brma_rstudent1[, "z"], brma_rstudent2[, "z"]) > 0.90,
              info = "rstudent residuals match")
})

test_that("Extended rstudent residuals align with metafor oracles", {

  skip_if_not_full_diagnostics(
    "Default residual coverage keeps one representative LOO-PIT check per family."
  )

  extended <- list(
    list(
      name      = "bcg_meta-regression4",
      resid_cor = 0.90,
      se_cor    = NULL,
      z_cor     = 0.90
    ),
    list(
      name      = "konstantopoulos2011_3lvl2",
      resid_cor = 0.60,
      se_cor    = 0.80,
      z_cor     = 0.60
    ),
    list(
      name      = "dat.lehmann2018-PETreg",
      equal     = TRUE
    )
  )
  skip_if_missing_fits(vapply(extended, `[[`, "", "name"))

  for (case in extended) {
    name <- case[["name"]]
    if (isTRUE(case[["equal"]])) {
      .expect_residual_table_equal(
        name,
        suppressWarnings(rstudent(fits[[name]])),
        rstudent(info[[name]][["metafor"]]),
        tolerance   = 0.05,
        z_tolerance = 0.10,
        label       = "rstudent"
      )
    } else {
      .expect_residual_rstudent_rank(
        name,
        suppressWarnings(rstudent(fits[[name]])),
        rstudent(info[[name]][["metafor"]]),
        resid_cor = case[["resid_cor"]],
        se_cor    = case[["se_cor"]],
        z_cor     = case[["z_cor"]]
      )
    }
  }
})

test_that("Extended rstudent residuals without direct metafor oracle are calibrated", {

  skip_if_not_full_diagnostics(
    "These checks have weaker or family-specific oracles and duplicate default branches."
  )

  scale_name <- "bangertdrowns2004_location-scale"
  sel_name   <- "dat.lehmann2018-3PSMreg"
  glmm_name  <- "bcg_glmm_reg"
  skip_if_missing_fits(c(scale_name, sel_name, glmm_name))

  scale_rstudent <- rstudent(fits[[scale_name]])
  scale_standard <- rstandard(info[[scale_name]][["metafor"]])
  expect_equal(mean(scale_rstudent$z), 0, 0.05)
  expect_equal(sd(scale_rstudent$z), 1, 0.05)
  expect_true(cor(scale_rstudent$z, scale_standard$z, method = "spearman") > 0.8)

  sel_resid    <- residuals(fits[[sel_name]])
  sel_rstudent <- suppressWarnings(rstudent(fits[[sel_name]]))
  expect_true(sel_rstudent$z[4] > 4)
  expect_true(all(abs(sel_rstudent$z[-4]) < 4))
  expect_true(cor(sel_resid, sel_rstudent$z) > 0.9)

  glmm_resid    <- residuals(fits[[glmm_name]])
  glmm_rstudent <- suppressWarnings(rstudent(fits[[glmm_name]], type = "estimate"))
  expect_true(cor(glmm_resid, glmm_rstudent$z) > 0.9)
})

test_that("Residuals for BMA.norm fits are internally consistent", {

  name <- "dat.lehmann2018_BMA.norm"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  n        <- nobs(fit_brma)

  resid_marginal <- residuals(fit_brma, conditioning_depth = "marginal")
  resid_estimate <- residuals(fit_brma, conditioning_depth = "estimate")

  expect_residual_vector(resid_marginal, n)
  expect_residual_vector(resid_estimate, n)
  expect_error(residuals(fit_brma, conditioning_depth = "cluster"),
               "is only available for multilevel models")
  expect_true(cor(resid_marginal, resid_estimate) > 0.9)
  expect_true(sd(resid_marginal) > sd(resid_estimate))

  pearson_marginal <- residuals(
    fit_brma,
    type               = "pearson",
    conditioning_depth = "marginal"
  )
  pearson_estimate <- residuals(
    fit_brma,
    type               = "pearson",
    conditioning_depth = "estimate"
  )

  expect_residual_vector(pearson_marginal, n)
  expect_residual_vector(pearson_estimate, n)
  expect_true(cor(pearson_marginal, pearson_estimate) > 0.9)

  brma_rstandard <- rstandard(fit_brma)
  expect_residual_table(brma_rstandard, n)
  expect_equal(residuals(fit_brma, type = "rstandard"), brma_rstandard[["z"]],
               tolerance = 1e-12)

  brma_rstudent <- suppressWarnings(rstudent(fit_brma))
  expect_residual_table(brma_rstudent, n)
  expect_equal(suppressWarnings(residuals(fit_brma, type = "rstudent")),
               brma_rstudent[["z"]], tolerance = 1e-12)
  expect_error(
    residuals(fit_brma, type = "rstudent", conditioning_depth = "marginal"),
    "do not set 'conditioning_depth'"
  )
})

test_that("Residuals for BMA.glmm fits are internally consistent", {

  model_names <- c("bcg_BMA.glmm", "nielweise2008_BMA.glmm", "bcg_BMA.glmm_3lvl_location_scale")
  skip_if_missing_fits(model_names)

  for (name in c("bcg_BMA.glmm", "nielweise2008_BMA.glmm")) {
    fit_brma <- fits[[name]]
    n        <- nobs(fit_brma)

    expect_residual_vector(residuals(fit_brma), n)
    expect_error(residuals(fit_brma, type = "pearson"), "normal outcome models")
    expect_error(rstandard(fit_brma), "normal outcome models")
    expect_residual_table(suppressWarnings(rstudent(fit_brma)), n)
  }

  fit_3lvl <- fits[["bcg_BMA.glmm_3lvl_location_scale"]]
  n_3lvl   <- nrow(fit_3lvl[["data"]][["outcome"]])

  expect_null(fit_3lvl[["marglik"]])
  expect_residual_vector(residuals(fit_3lvl, conditioning_depth = "cluster"),
                         n_3lvl)
  expect_error(rstudent(fit_3lvl, unit = "cluster"), "Cluster-unit rstudent residuals")
})

test_that("Residuals for RoBMA fits are internally consistent", {

  model_names <- c("dat.lehmann2018_RoBMA", "dat.lehmann2018_RoBMA_3lvl_mods_scale")
  skip_if_missing_fits(model_names)

  fit_brma <- fits[["dat.lehmann2018_RoBMA"]]
  n        <- nobs(fit_brma)

  expect_residual_vector(residuals(fit_brma), n)
  expect_error(residuals(fit_brma, type = "pearson"), "selection models")
  expect_error(rstandard(fit_brma), "selection models")

  brma_rstudent <- suppressWarnings(rstudent(fit_brma))
  expect_residual_table(brma_rstudent, n)
  expect_equal(suppressWarnings(residuals(fit_brma, type = "rstudent")),
               brma_rstudent[["z"]], tolerance = 1e-12)

  fit_3lvl <- fits[["dat.lehmann2018_RoBMA_3lvl_mods_scale"]]
  n_3lvl   <- nrow(fit_3lvl[["data"]][["outcome"]])

  expect_null(fit_3lvl[["marglik"]])
  expect_residual_vector(residuals(fit_3lvl, conditioning_depth = "cluster"),
                         n_3lvl)
})


.known_v_first_block_cdf_oracle <- function(fit_brma) {

  setup         <- .estimate_likelihood_setup.brma(fit_brma)
  known_V       <- .data_known_v_data(setup[["data"]])
  block         <- known_V[["block_indices"]][[1L]]
  s             <- 1L
  extra         <- .known_v_extra_variance_from_setup(setup)
  covariance    <- known_V[["V"]][block, block, drop = FALSE] +
    diag(extra[s, block], nrow = length(block))
  yi            <- setup[["yi"]][block]
  mu            <- setup[["mu"]][s, block]
  lower_tail    <- TRUE

  if (identical(setup[["effect_direction"]], "negative")) {
    yi         <- -yi
    mu         <- -mu
    lower_tail <- FALSE
  }

  if (length(block) == 1L) {
    residual <- yi - mu
    variance <- covariance[1L, 1L]
  } else {
    precision <- chol2inv(chol(covariance))
    residual  <- as.vector(precision %*% (yi - mu)) / diag(precision)
    variance  <- 1 / diag(precision)
  }

  stats::pnorm(
    residual,
    mean       = 0,
    sd         = sqrt(variance),
    lower.tail = lower_tail
  )
}


.known_v_rstandard_gls_oracle <- function(fit_brma) {

  setup         <- .estimate_likelihood_setup.brma(fit_brma)
  covariance_samples <- .known_v_marginal_covariance_samples(
    object            = fit_brma,
    posterior_samples = setup[["posterior_samples"]]
  )
  X             <- .get_model_matrix(fit_brma)
  yi            <- setup[["yi"]]
  K             <- setup[["K"]]
  S             <- setup[["S"]]
  I_K           <- diag(K)
  resid_samples <- matrix(NA_real_, nrow = S, ncol = K)
  se_samples    <- matrix(NA_real_, nrow = S, ncol = K)
  z_samples     <- matrix(NA_real_, nrow = S, ncol = K)

  for (s in seq_len(S)) {
    covariance <- covariance_samples[s, , ]
    W          <- chol2inv(chol(covariance))
    WX         <- W %*% X
    XtWX_inv   <- .hat_solve_crossprod(crossprod(X, WX))
    H          <- X %*% XtWX_inv %*% t(WX)
    beta_hat   <- as.vector(XtWX_inv %*% crossprod(X, as.vector(W %*% yi)))
    residual   <- yi - as.vector(X %*% beta_hat)
    se         <- sqrt(pmax(diag((I_K - H) %*% covariance %*% t(I_K - H)), 0))

    resid_samples[s, ] <- residual
    se_samples[s, ]    <- se
    z_samples[s, ]     <- residual / se
  }

  data.frame(
    resid = colMeans(resid_samples),
    se    = colMeans(se_samples),
    z     = colMeans(z_samples)
  )
}


.known_v_rstandard_estimate_oracle <- function(fit_brma) {

  setup          <- .estimate_likelihood_setup.brma(fit_brma)
  known_V        <- .data_known_v_data(setup[["data"]])
  extra          <- .known_v_extra_variance_from_setup(setup)
  X              <- .get_model_matrix(fit_brma)
  yi             <- setup[["yi"]]
  K              <- setup[["K"]]
  S              <- setup[["S"]]
  offset_samples <- setup[["mu_random"]]
  if (is.null(offset_samples)) {
    offset_samples <- matrix(0, nrow = S, ncol = K)
  }

  resid_samples <- matrix(NA_real_, nrow = S, ncol = K)
  se_samples    <- matrix(NA_real_, nrow = S, ncol = K)
  z_samples     <- matrix(NA_real_, nrow = S, ncol = K)

  for (s in seq_len(S)) {
    covariance <- known_V[["V"]] + diag(extra[s, ], nrow = K)
    W          <- chol2inv(chol(covariance))
    WX         <- W %*% X
    XtWX_inv   <- .hat_solve_crossprod(crossprod(X, WX))
    y_offset   <- yi - offset_samples[s, ]
    beta_hat   <- as.vector(XtWX_inv %*% crossprod(X, as.vector(W %*% y_offset)))
    raw_resid  <- y_offset - as.vector(X %*% beta_hat)
    residual   <- as.vector(known_V[["V"]] %*% W %*% raw_resid)
    C          <- W - WX %*% XtWX_inv %*% t(WX)
    se         <- sqrt(pmax(diag(known_V[["V"]] %*% C %*% known_V[["V"]]), 0))

    resid_samples[s, ] <- residual
    se_samples[s, ]    <- se
    z_samples[s, ]     <- residual / se
  }

  data.frame(
    resid = colMeans(resid_samples),
    se    = colMeans(se_samples),
    z     = colMeans(z_samples)
  )
}


test_that("brma.mv known-V residual diagnostics are internally consistent", {

  mv_names <- c(
    "brma.mv_latent",
    "brma.mv_whitened",
    "brma.mv_block_mvn",
    "brma.mv_block_mvn_random_scale",
    "brma.mv_block_mvn_known_R"
  )
  skip_if_missing_fits(mv_names)

  for (name in mv_names) {
    fit_brma  <- fits[[name]]
    n         <- nobs(fit_brma)
    is_random <- .is_random(fit_brma)

    expect_residual_vector(residuals(fit_brma), n, info = name)
    expect_residual_vector(
      residuals(fit_brma, conditioning_depth = "estimate"),
      n,
      info = name
    )

    pearson <- residuals(fit_brma, type = "pearson")
    standard <- rstandard(fit_brma)

    expect_residual_vector(pearson, n, info = name)
    expect_residual_table(standard, n, info = name)

    if (is_random) {
      pearson_estimate <- residuals(
        fit_brma,
        type               = "pearson",
        conditioning_depth = "estimate"
      )
      standard_estimate <- rstandard(fit_brma, conditioning_depth = "estimate")

      expect_residual_vector(pearson_estimate, n, info = name)
      expect_residual_table(standard_estimate, n, info = name)
    }

    student <- suppressWarnings(rstudent(fit_brma))
    expect_residual_table(student, n, info = name)
    expect_equal(
      suppressWarnings(residuals(fit_brma, type = "LOO-PIT")),
      student[["z"]],
      tolerance = 1e-12,
      info      = name
    )
  }
})

test_that("v14 brma.mv residual diagnostics return finite estimate-unit output", {

  mv_names <- c(
    "brma.mv_v14_konstantopoulos2011_cs",
    "brma.mv_v14_assink2016_nested",
    "brma.mv_v14_ishak2007_har",
    "brma.mv_v14_begg1989_study_treatment"
  )
  skip_if_missing_fits(mv_names)

  for (name in mv_names) {
    fit_brma <- fits[[name]]
    n        <- nobs(fit_brma)

    expect_residual_vector(residuals(fit_brma), n, info = name)
    expect_residual_vector(
      residuals(fit_brma, conditioning_depth = "estimate"),
      n,
      info = name
    )
    expect_residual_vector(
      residuals(fit_brma, type = "pearson"),
      n,
      info = name
    )
    expect_residual_table(rstandard(fit_brma), n, info = name)
    expect_residual_table(
      rstandard(fit_brma, conditioning_depth = "estimate"),
      n,
      info = name
    )
    expect_residual_table(suppressWarnings(rstudent(fit_brma)), n, info = name)

    fit_missing <- fit_brma
    fit_missing[["loo"]] <- NULL
    expect_error(
      rstudent(fit_missing),
      "LOO has not been computed",
      info = name
    )
  }
})


test_that("brma.mv known-V residual diagnostics match Schur and GLS oracles", {

  name <- "brma.mv_block_mvn"
  skip_if_missing_fits(name)

  fit_brma    <- fits[[name]]
  setup       <- .estimate_likelihood_setup.brma(fit_brma)
  cdf_matrix  <- .cdf_lik_estimate.brma(fit_brma)
  known_V     <- .data_known_v_data(setup[["data"]])
  first_block <- known_V[["block_indices"]][[1L]]

  expect_equal(
    unname(cdf_matrix[1L, first_block]),
    unname(.known_v_first_block_cdf_oracle(fit_brma)),
    tolerance = 1e-12
  )

  expect_equal(
    unname(.known_v_pearson_residual_se_samples(fit_brma, "marginal")),
    unname(sqrt(
      matrix(
        diag(known_V[["V"]]),
        nrow  = setup[["S"]],
        ncol  = setup[["K"]],
        byrow = TRUE
      ) + .known_v_extra_variance_from_setup(setup)
    )),
    tolerance = 1e-12
  )
  expect_equal(
    unname(.known_v_pearson_residual_se_samples(fit_brma, "estimate")),
    unname(sqrt(matrix(
      diag(known_V[["V"]]),
      nrow  = setup[["S"]],
      ncol  = setup[["K"]],
      byrow = TRUE
    ))),
    tolerance = 1e-12
  )

  expect_equal(
    unname(as.matrix(rstandard(fit_brma))),
    unname(as.matrix(.known_v_rstandard_gls_oracle(fit_brma))),
    tolerance = 1e-10
  )
  expect_equal(
    unname(as.matrix(rstandard(fit_brma, conditioning_depth = "estimate"))),
    unname(as.matrix(.known_v_rstandard_estimate_oracle(fit_brma))),
    tolerance = 1e-10
  )
})


test_that("brma.mv random known-V estimate residuals use BLUP and sampled offsets", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  setup    <- .estimate_likelihood_setup.brma(fit_brma)
  yi_mat   <- matrix(
    setup[["yi"]],
    nrow  = setup[["S"]],
    ncol  = setup[["K"]],
    byrow = TRUE
  )
  expected_resid <- colMeans(yi_mat - .known_v_estimate_blup_from_setup(setup))
  known_V        <- .data_known_v_data(setup[["data"]])
  covariance_samples <- .known_v_marginal_covariance_samples(
    object            = fit_brma,
    posterior_samples = setup[["posterior_samples"]]
  )
  marginal_se <- matrix(NA_real_, nrow = setup[["S"]], ncol = setup[["K"]])
  for (s in seq_len(setup[["S"]])) {
    marginal_se[s, ] <- sqrt(diag(covariance_samples[s, , ]))
  }

  expect_equal(
    unname(residuals(fit_brma, conditioning_depth = "estimate")),
    unname(expected_resid),
    tolerance = 1e-12
  )
  expect_equal(
    unname(.known_v_pearson_residual_se_samples(fit_brma, "marginal")),
    unname(marginal_se),
    tolerance = 1e-12
  )
  expect_equal(
    unname(.known_v_pearson_residual_se_samples(fit_brma, "estimate")),
    unname(sqrt(matrix(
      diag(known_V[["V"]]),
      nrow  = setup[["S"]],
      ncol  = setup[["K"]],
      byrow = TRUE
    ))),
    tolerance = 1e-12
  )
  expect_equal(
    unname(as.matrix(rstandard(fit_brma, conditioning_depth = "estimate"))),
    unname(as.matrix(.known_v_rstandard_estimate_oracle(fit_brma))),
    tolerance = 1e-10
  )
  expect_equal(
    unname(as.matrix(rstandard(fit_brma))),
    unname(as.matrix(.known_v_rstandard_gls_oracle(fit_brma))),
    tolerance = 1e-10
  )
})


test_that("Known-V Pearson residuals chunk marginal covariance without changing results", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  fit_brma     <- fits[[name]]
  expected     <- residuals(fit_brma, type = "pearson")
  K            <- nobs(fit_brma)
  one_draw_mem <- .known_v_covariance_bytes(1L, K)
  old_options  <- options(
    RoBMA.known_v_covariance_max_bytes = 2 * one_draw_mem
  )
  on.exit(options(old_options), add = TRUE)

  actual   <- residuals(fit_brma, type = "pearson")
  metadata <- attr(actual, "known_v_diagnostic")
  attr(actual, "known_v_diagnostic") <- NULL

  expect_equal(unname(actual), unname(expected), tolerance = 1e-10)
  expect_true(metadata[["n_chunks"]] > 1L)
  expect_equal(metadata[["n_used_samples"]], metadata[["n_posterior_samples"]])
})


test_that("legacy known-V CDF wrapper uses estimate-unit target", {

  name <- "brma.mv_block_mvn_random_scale"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]

  cdf_legacy <- .cdf.brma(fit_brma, conditioning_depth = "estimate")
  cdf_target <- .cdf_lik_estimate.brma(fit_brma)

  expect_equal(unname(cdf_legacy), unname(cdf_target))
  expect_error(
    .cdf.brma(fit_brma, conditioning_depth = "marginal"),
    "conditioning_depth = 'estimate'"
  )
})


test_that("brma.mv known-V residual consumers accept residual tables", {

  name <- "brma.mv_block_mvn"
  skip_if_missing_fits(name)

  fit_brma <- fits[[name]]
  qq_data  <- qqnorm(fit_brma, as_data = TRUE)
  fun_data <- funnel(fit_brma, residual = TRUE, as_data = TRUE)

  expect_true(is.list(qq_data))
  expect_true(is.list(fun_data))
  expect_true("points" %in% names(qq_data))
  expect_true("points" %in% names(fun_data))
})
