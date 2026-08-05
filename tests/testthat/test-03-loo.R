# ============================================================================ #
#
# Test LOO-PSIS functionality for brma objects
#
# ============================================================================
context("loo methods for brma objects")

# The LOO-PSIS functionality is necessary for the residuals and funnel plot
# functionality. Due to the computational costs (and the possibility to test)
# against other metafor output) it is primary tested therein.

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-contracts.R"))
skip_on_cran()
skip_if_no_fits()
skip_if_not_installed("loo")

# list cached fits lazily
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)
info      <- lazy_infos(fit_names, validate = FALSE)


.expected_known_v_estimate_target <- function(object) {

  known_V <- .data_known_v_data(object[["data"]])

  if (length(.known_v_correlated_blocks(known_V)) > 0L) {
    return("known_v_estimate")
  }

  return("factorized_estimate")
}


test_that("log_lik returns finite pointwise matrices for cached fits", {

  for (name in c("bcg_meta-analysis", "bcg_meta-regression")) {
    fit_brma <- fits[[name]]
    log_lik  <- log_lik(fit_brma)

    expect_true(is.matrix(log_lik), info = name)
    expect_equal(ncol(log_lik), nobs(fit_brma), info = name)
    expect_true(all(is.finite(log_lik)), info = name)
    expect_true(stats::var(rowSums(log_lik)) > 0, info = name)
  }
})

test_that("logLik is absent and information criteria fail explicitly", {

  fit_name <- "brma.mv_block_mvn"
  skip_if_missing_fits(fit_name)
  fit_brma <- fits[[fit_name]]
  draws    <- log_lik(fit_brma)

  expect_true("log_lik" %in% getNamespaceExports("RoBMA"))
  expect_true(is.matrix(draws))
  expect_equal(ncol(draws), nobs(fit_brma))
  expect_null(getS3method("logLik", "brma", optional = TRUE))
  expect_error(stats::logLik(fit_brma), "no applicable method")
  expect_error(stats::AIC(fit_brma), "AIC\\(\\) is not defined for brma")
  expect_error(stats::BIC(fit_brma), "BIC\\(\\) is not defined for brma")
})

test_that("loo and WAIC expose cached diagnostics and missing-cache errors", {

  for (name in c("bcg_meta-analysis", "bcg_meta-regression")) {
    fit_brma <- fits[[name]]
    loo_brma <- suppressWarnings(loo(fit_brma))

    expect_s3_class(loo_brma, "loo")
    expect_true(all(c("elpd_loo", "p_loo", "looic") %in% rownames(loo_brma[["estimates"]])),
                info = name)
    expect_true(all(is.finite(loo_brma[["estimates"]][, "Estimate"])), info = name)

    fit_missing <- fit_brma
    fit_missing[["loo"]] <- NULL
    expect_error(loo(fit_missing), "LOO has not been computed", info = name)
    expect_error(waic(fit_missing), "WAIC has not been computed", info = name)

    fit_waic <- suppressWarnings(add_waic(fit_missing))
    waic_brma <- waic(fit_waic)
    waic_loo  <- loo::waic(fit_waic)

    expect_s3_class(waic_brma, "waic")
    expect_s3_class(waic_loo, "waic")
    expect_equal(waic_loo, waic_brma)
    expect_true(all(c("elpd_waic", "p_waic", "waic") %in% rownames(waic_brma[["estimates"]])),
                info = name)
    expect_true(all(is.finite(waic_brma[["estimates"]][, "Estimate"])), info = name)
  }
})

test_that("loo_compare compares two brma models", {

  # get two brma fits
  fit_brma  <- fits[["bcg_meta-analysis"]]
  fit_brma2 <- fits[["bcg_meta-regression"]]

  # compare
  out <- suppressWarnings(loo_compare(fit_brma, fit_brma2))

  # check structure
  expect_identical(class(out), c("compare.loo", "matrix", "array"))
  expect_equal(nrow(out), 2)
  expect_true(all(c("elpd_diff", "se_diff") %in% colnames(out)))
})

test_that("loo_compare accepts loo objects", {

  # get two brma fits
  fit_brma  <- fits[["bcg_meta-analysis"]]
  fit_brma2 <- fits[["bcg_meta-regression"]]

  # get two loo brma fits
  loo_brma  <- loo(fits[["bcg_meta-analysis"]])
  loo_brma2 <- suppressWarnings(loo(fits[["bcg_meta-regression"]]))

  # compare
  out <- loo_compare(loo_brma, loo_brma2)

  # check structure
  expect_identical(class(out), c("compare.loo", "matrix", "array"))
  expect_equal(nrow(out), 2)
  expect_true(all(c("elpd_diff", "se_diff") %in% colnames(out)))
})

test_that("cluster-unit LOO has a comparable joint-target label", {

  fit_name <- "konstantopoulos2011_3lvl"
  skip_if_missing_fits(fit_name)

  fit_brma <- suppressWarnings(add_loo(fits[[fit_name]], unit = "cluster"))
  loo_brma <- loo(fit_brma, unit = "cluster")
  target   <- attr(loo_brma, "RoBMA_target", exact = TRUE)
  loo_copy <- loo_brma
  out      <- loo_compare(loo_brma, loo_copy)

  expect_equal(target[["target"]], "cluster_joint")
  expect_equal(nrow(out), 2L)
})

test_that("cluster-unit GLMM likelihood requires certified nested quadrature", {

  fit_name <- "bcg_glmm_3lvl_scale"
  skip_if_missing_fits(fit_name)
  fit_brma <- fits[[fit_name]]

  expect_error(
    log_lik(fit_brma, unit = "cluster"),
    "Cluster-unit GLMM log-likelihood is unavailable"
  )
  expect_error(
    add_loo(fit_brma, unit = "cluster"),
    "Cluster-unit GLMM log-likelihood is unavailable"
  )
  expect_error(
    add_waic(fit_brma, unit = "cluster"),
    "Cluster-unit GLMM log-likelihood is unavailable"
  )
})


test_that("fitted GLMM methods support point nuisance priors", {

  fit_names <- c("bcg_glmm", "nielweise2008_glmm")
  skip_if_missing_fits(fit_names)

  fit_bin <- fits[[fit_names[[1L]]]]
  fit_bin[["priors"]][["outcome"]][["pi"]] <- BayesTools::prior_factor(
    "point", list(location = 0.5), contrast = "independent"
  )
  fit_bin[["loo"]]  <- NULL
  fit_bin[["waic"]] <- NULL

  fit_pois <- fits[[fit_names[[2L]]]]
  fit_pois[["priors"]][["outcome"]][["phi"]] <- BayesTools::prior_factor(
    "point", list(location = -1), contrast = "independent"
  )
  fit_pois[["loo"]]  <- NULL
  fit_pois[["waic"]] <- NULL

  for (fit_brma in list(fit_bin, fit_pois)) {
    log_lik <- log_lik(fit_brma)
    expect_true(is.matrix(log_lik))
    expect_true(all(is.finite(log_lik)))

    fit_loo <- suppressWarnings(add_loo(fit_brma))
    expect_s3_class(loo(fit_loo), "loo")

    fit_waic <- suppressWarnings(add_waic(fit_brma))
    expect_s3_class(waic(fit_waic), "waic")
  }
})


test_that("loo_compare rejects fewer than two models", {

  # get one brma fit
  fit_brma <- fits[["bcg_meta-analysis"]]

  expect_error(loo_compare(fit_brma, "At least two models"))
})

test_that("log_lik, LOO, weights, diagnostics, and WAIC are available for product-space fits", {

  product_names <- c(
    "dat.lehmann2018_BMA.norm",
    "bcg_BMA.glmm",
    "dat.lehmann2018_RoBMA"
  )
  skip_if_missing_fits(product_names)

  for (name in product_names) {

    fit_brma <- fits[[name]]

    log_lik <- log_lik(fit_brma)
    expect_true(is.matrix(log_lik), info = name)
    expect_equal(ncol(log_lik), nobs(fit_brma), info = name)
    expect_true(all(is.finite(log_lik)), info = name)

    loo_result <- loo(fit_brma)
    expect_s3_class(loo_result, "loo")

    weights <- loo_weights(fit_brma)
    expect_true(is.matrix(weights), info = name)
    expect_equal(dim(weights), dim(log_lik), info = name)
    expect_equal(colSums(weights), rep(1, ncol(weights)), tolerance = 1e-10,
                 info = name)

    expect_no_error(suppressWarnings(check_loo(fit_brma)))

    if (identical(name, "dat.lehmann2018_BMA.norm")) {
      fit_waic <- fit_brma
      fit_waic[["waic"]] <- NULL
      fit_waic   <- suppressWarnings(add_waic(fit_waic))
      waic_result <- waic(fit_waic)
    } else {
      waic_result <- suppressWarnings(loo::waic(log_lik))
    }
    expect_s3_class(waic_result, "waic")
  }
})

test_that("LOO comparison and model weights support product-space fits", {

  product_names <- c("dat.lehmann2018_BMA.norm", "dat.lehmann2018_RoBMA")
  skip_if_missing_fits(product_names)

  out <- suppressWarnings(loo_compare(
    fits[["dat.lehmann2018_BMA.norm"]],
    fits[["dat.lehmann2018_RoBMA"]]
  ))

  expect_true(is.matrix(out))
  expect_equal(nrow(out), 2)
  expect_true("elpd_diff" %in% colnames(out))
  expect_true("se_diff" %in% colnames(out))

  model_weights <- suppressWarnings(loo_model_weights(
    fits[["dat.lehmann2018_BMA.norm"]],
    fits[["dat.lehmann2018_RoBMA"]]
  ))
  expect_length(model_weights, 2L)
  expect_equal(sum(model_weights), 1, tolerance = 1e-12)
  expect_true(all(model_weights >= 0 & model_weights <= 1))
})


test_that("brma.mv known-V fits expose conditional estimate-unit LOO and WAIC", {

  mv_names <- "brma.mv_block_mvn"
  if (is_certification_profile()) {
    mv_names <- c(
      "brma.mv_latent",
      "brma.mv_whitened",
      mv_names,
      "brma.mv_block_mvn_fixed_random_null",
      "brma.mv_block_mvn_known_R",
      "brma.mv_block_mvn_random_scale",
      "brma.mv_block_mvn_3lvl_scale_total",
      "brma.mv_block_mvn_3lvl_scale_top",
      "brma.mv_block_mvn_3lvl_scale_bottom"
    )
    mv_names <- intersect(mv_names, active_fit_catalog()[["name"]])
    if (length(mv_names) == 0L) {
      testthat::skip("No known-V LOO fixtures are active in this case.")
    }
  }
  skip_if_missing_fits(mv_names)

  for (name in mv_names) {
    fit_brma <- fits[[name]]
    is_known_r <- identical(name, "brma.mv_block_mvn_known_R")
    if (is_known_r) {
      fit_brma <- suppressWarnings(add_loo(fit_brma))
    }
    log_lik  <- log_lik(fit_brma)
    target   <- attr(log_lik, "RoBMA_target", exact = TRUE)

    expect_equal(ncol(log_lik), nobs(fit_brma), info = name)
    expect_true(all(is.finite(log_lik)), info = name)
    expect_equal(
      target[["target"]],
      .expected_known_v_estimate_target(fit_brma),
      info = name
    )
    expect_true(isTRUE(target[["known_v"]]), info = name)
    expect_true(isTRUE(target[["known_v_estimate_backend"]]), info = name)
    expect_equal(
      target[["known_v_parameterization"]],
      .data_known_v_data(fit_brma[["data"]])[["parameterization"]],
      info = name
    )
    expect_equal(
      target[["known_v_parameterization_requested"]],
      .data_known_v_data(fit_brma[["data"]])[["parameterization_requested"]],
      info = name
    )
    expect_true(length(target[["dependency_component_sizes"]]) > 0L, info = name)
    expect_equal(isTRUE(target[["known_r"]]), is_known_r, info = name)
    if (is_known_r) {
      expect_equal(target[["known_r_blocks"]], "study", info = name)
      expect_true(grepl("conditions on sampled fitted random effects",
                        target[["known_r_semantics"]], fixed = TRUE),
                  info = name)
      expect_equal(target[["random_effects"]], "conditioned", info = name)
    } else {
      expect_equal(length(target[["known_r_blocks"]]), 0L, info = name)
      expect_null(target[["known_r_semantics"]], info = name)
    }

    loo_result <- loo(fit_brma)
    loo_target <- attr(loo_result, "RoBMA_target", exact = TRUE)
    expect_s3_class(loo_result, "loo")
    expect_equal(loo_target[["target"]], target[["target"]], info = name)
    expect_true(isTRUE(loo_target[["known_v_estimate_backend"]]), info = name)
    expect_equal(loo_target[["known_r"]], target[["known_r"]], info = name)
    expect_equal(loo_target[["known_r_blocks"]], target[["known_r_blocks"]],
                 info = name)
    expect_equal(loo_target[["known_r_semantics"]],
                 target[["known_r_semantics"]], info = name)
    expect_equal(
      loo_target[["dependency_component_sizes"]],
      target[["dependency_component_sizes"]],
      info = name
    )

    fit_waic   <- suppressWarnings(add_waic(fit_brma))
    waic_result <- waic(fit_waic)
    waic_target <- attr(waic_result, "RoBMA_target", exact = TRUE)
    expect_s3_class(waic_result, "waic")
    expect_equal(waic_target[["target"]], target[["target"]], info = name)
    expect_true(isTRUE(waic_target[["known_v_estimate_backend"]]), info = name)
    expect_equal(waic_target[["known_r"]], target[["known_r"]], info = name)
    expect_equal(waic_target[["known_r_blocks"]], target[["known_r_blocks"]],
                 info = name)
    expect_equal(waic_target[["known_r_semantics"]],
                 target[["known_r_semantics"]], info = name)
    expect_equal(
      waic_target[["dependency_component_sizes"]],
      target[["dependency_component_sizes"]],
      info = name
    )
  }
})

test_that("v14 brma.mv metafor fixtures cache usable estimate-unit LOO", {

  skip_if_not_certification("This case exercises the high-draw v14 fixtures.")
  mv_names <- c(
    "brma.mv_v14_konstantopoulos2011_cs",
    "brma.mv_v14_assink2016_nested",
    "brma.mv_v14_ishak2007_har",
    "brma.mv_v14_begg1989_study_treatment"
  )
  mv_names <- intersect(mv_names, active_fit_catalog()[["name"]])
  if (length(mv_names) == 0L) {
    testthat::skip("No v14 LOO fixtures are active in this case.")
  }
  skip_if_missing_fits(mv_names)

  for (name in mv_names) {
    fit_brma <- fits[[name]]
    log_lik  <- log_lik(fit_brma)
    target   <- attr(log_lik, "RoBMA_target", exact = TRUE)

    expect_equal(ncol(log_lik), nobs(fit_brma), info = name)
    expect_true(all(is.finite(log_lik)), info = name)
    expect_equal(
      target[["target"]],
      .expected_known_v_estimate_target(fit_brma),
      info = name
    )
    expect_true(isTRUE(target[["known_v_estimate_backend"]]), info = name)

    loo_result <- suppressWarnings(loo(fit_brma))
    loo_target <- attr(loo_result, "RoBMA_target", exact = TRUE)
    weights    <- suppressWarnings(loo_weights(fit_brma))

    expect_s3_class(loo_result, "loo")
    expect_true(all(c("elpd_loo", "p_loo", "looic") %in%
      rownames(loo_result[["estimates"]])), info = name)
    expect_true(all(is.finite(loo_result[["estimates"]][, "Estimate"])),
                info = name)
    expect_equal(loo_target[["target"]], target[["target"]], info = name)
    expect_true(isTRUE(loo_target[["known_v_estimate_backend"]]), info = name)
    expect_true(is.matrix(weights), info = name)
    expect_equal(dim(weights), dim(log_lik), info = name)
    expect_equal(colSums(weights), rep(1, ncol(weights)), tolerance = 1e-10,
                 info = name)
  }
})

test_that("brma.mv estimate-unit LOO rejects unsupported and missing targets", {

  mv_name      <- "brma.mv_block_mvn"
  regular_name <- "bcg_meta-analysis"
  skip_if_missing_fits(c(mv_name, regular_name))

  fit_brma <- fits[[mv_name]]

  expect_error(
    log_lik(fit_brma, unit = "cluster"),
    "cluster"
  )
  expect_error(
    add_loo(fit_brma, unit = "cluster"),
    "cluster"
  )

  fit_missing <- fit_brma
  fit_missing[["loo"]] <- NULL
  expect_error(
    loo(fit_missing),
    "LOO has not been computed"
  )

  expect_error(
    suppressWarnings(loo_compare(fit_brma, fits[[regular_name]])),
    "same outcome target|same data|target"
  )
})


test_that("loo_compare accepts brma.mv backends with the same known-V target", {

  mv_names <- c("brma.mv_latent", "brma.mv_whitened", "brma.mv_block_mvn")
  skip_if_missing_fits(mv_names)

  out <- suppressWarnings(loo_compare(
    fits[["brma.mv_latent"]],
    fits[["brma.mv_whitened"]],
    fits[["brma.mv_block_mvn"]]
  ))

  expect_true(is.matrix(out))
  expect_equal(nrow(out), 3)
  expect_true("elpd_diff" %in% colnames(out))
})

test_that("LOO and WAIC compare random-effect representations on one target", {

  representation_names <- c(
    none         = "brma.mv_block_mvn",
    marginalized = "brma.mv_block_mvn_random"
  )
  if (is_certification_profile()) {
    representation_names <- c(
      representation_names,
      sampled = "brma.mv_block_mvn_known_R",
      mixed   = "brma.mv_block_mvn_random_scale"
    )
  }
  skip_if_missing_fits(unname(representation_names))

  representation_fits <- lapply(representation_names, function(name) {
    fit_brma <- suppressWarnings(add_loo(fits[[name]]))
    suppressWarnings(add_waic(fit_brma))
  })
  targets <- lapply(representation_fits, function(fit_brma) {
    attr(log_lik(fit_brma), "RoBMA_target", exact = TRUE)
  })

  expected_random_effects <- c(none = "none", marginalized = "none")
  expected_estimate_level <- c(none = "none", marginalized = "marginalized")
  if (is_certification_profile()) {
    expected_random_effects <- c(
      expected_random_effects,
      sampled = "conditioned",
      mixed   = "conditioned"
    )
    expected_estimate_level <- c(
      expected_estimate_level,
      sampled = "none",
      mixed   = "marginalized"
    )
  }

  expect_equal(
    vapply(targets, function(target) target[["random_effects"]], character(1)),
    expected_random_effects
  )
  expect_equal(
    vapply(targets, function(target) {
      target[["estimate_level_random"]]
    }, character(1)),
    expected_estimate_level
  )

  comparison_fields <- c("data_hash", "unit", "conditioning_depth", "target")
  for (field in comparison_fields) {
    values <- vapply(targets, function(target) target[[field]], character(1))
    expect_equal(length(unique(values)), 1L, info = field)
  }

  loo_comparison <- suppressWarnings(do.call(
    loo_compare,
    unname(lapply(representation_fits, loo))
  ))
  waic_comparison <- suppressWarnings(do.call(
    loo_compare,
    unname(lapply(representation_fits, waic))
  ))

  expect_equal(nrow(loo_comparison), length(representation_names))
  expect_equal(nrow(waic_comparison), length(representation_names))
})

test_that("matched sampled and marginalized effects have equivalent scores", {

  skip_if_not_certification(
    "Matched sampled-effect scores require certification-only fits."
  )

  representation_names <- c(
    marginalized = "brma.mv_block_mvn_random",
    sampled      = "brma.mv_block_mvn_random_sampled"
  )
  skip_if_missing_fits(unname(representation_names))

  representation_fits <- lapply(representation_names, function(name) {
    fit_brma <- suppressWarnings(add_loo(fits[[name]]))
    suppressWarnings(add_waic(fit_brma))
  })
  loo_comparison <- suppressWarnings(do.call(
    loo_compare,
    unname(lapply(representation_fits, loo))
  ))
  waic_comparison <- suppressWarnings(do.call(
    loo_compare,
    unname(lapply(representation_fits, waic))
  ))

  expect_lt(
    abs(loo_comparison[2L, "elpd_diff"]),
    0.25
  )
  expect_lt(
    abs(waic_comparison[2L, "elpd_diff"]),
    0.25
  )
})

test_that("marginalized normal score equals the sampled-effect convolution", {

  yi  <- c(-0.4, 0.2, 0.9)
  mu  <- c(-0.1, 0.3, 0.5)
  tau <- c(0.2, 0.5, 0.8)
  sei <- c(0.1, 0.25, 0.4)

  integrated <- vapply(seq_along(yi), function(i) {
    stats::integrate(
      function(random_effect) {
        stats::dnorm(
          yi[[i]],
          mean = mu[[i]] + random_effect,
          sd   = sei[[i]]
        ) * stats::dnorm(random_effect, mean = 0, sd = tau[[i]])
      },
      lower       = -Inf,
      upper       = Inf,
      subdivisions = 1000L,
      rel.tol     = 1e-11
    )[["value"]]
  }, numeric(1))
  marginalized <- exp(.outcome_pdf.norm(
    yi         = yi,
    mu_samples = matrix(mu, nrow = 1L),
    tau_within = matrix(tau, nrow = 1L),
    sei        = sei
  )[1L, ])

  expect_equal(marginalized, integrated, tolerance = 1e-9)
})

test_that("sampled local effects expose their observed PSIS reliability risk", {

  skip_if_not_certification(
    "The sampled-effect reliability contrast uses a certification fixture."
  )

  marginalized <- fits[["brma.mv_block_mvn_random"]]
  sampled      <- fits[["brma.mv_block_mvn_random_sampled"]]
  marginalized_k <- loo::pareto_k_values(loo(marginalized))
  sampled_k      <- loo::pareto_k_values(loo(sampled))

  expect_false(any(marginalized_k > 0.7))
  expect_true(any(sampled_k > 0.7))
  expect_warning(
    check_loo(sampled),
    "Some Pareto k values are high"
  )
})

test_that("LOO and WAIC cannot be mixed in one comparison table", {

  fit_brma <- fits[["brma.mv_block_mvn_random"]]
  fit_brma <- suppressWarnings(add_waic(fit_brma))

  expect_error(
    loo_compare(loo(fit_brma), waic(fit_brma)),
    "cannot be compared in the same table"
  )
})

test_that("cached LOO and WAIC reject stale outcome and known-V targets", {

  mv_name <- "brma.mv_block_mvn"
  skip_if_missing_fits(mv_name)
  fit_mv <- fits[[mv_name]]
  fit_mv <- suppressWarnings(add_loo(fit_mv))
  fit_mv <- suppressWarnings(add_waic(fit_mv))

  stale_outcome <- fit_mv
  stale_outcome[["data"]][["outcome"]][["yi"]][1L] <-
    stale_outcome[["data"]][["outcome"]][["yi"]][1L] + 0.01
  expect_error(loo(stale_outcome), "current outcome data")
  expect_error(waic(stale_outcome), "current outcome data")

  stale_V <- fit_mv
  known_V <- .data_known_v_data(stale_V[["data"]])
  if (!is.null(known_V[["blocks"]])) {
    known_V[["blocks"]][[1L]][["covariance"]][1L, 2L] <-
      known_V[["blocks"]][[1L]][["covariance"]][1L, 2L] + 0.001
    known_V[["blocks"]][[1L]][["covariance"]][2L, 1L] <-
      known_V[["blocks"]][[1L]][["covariance"]][1L, 2L]
  } else if (!is.null(known_V[["V"]])) {
    known_V[["V"]][1L, 2L] <- known_V[["V"]][1L, 2L] + 0.001
    known_V[["V"]][2L, 1L] <- known_V[["V"]][1L, 2L]
  } else {
    known_V[["diagonal"]][1L] <- known_V[["diagonal"]][1L] + 0.01
  }
  attr(stale_V[["data"]], "known_V_data") <- known_V

  expect_error(loo(stale_V), "current outcome data")
  expect_error(waic(stale_V), "current outcome data")
})

test_that("cached target checks ignore provenance but reject target changes", {

  fit_brma <- fits[["brma.mv_block_mvn"]]
  fit_brma <- suppressWarnings(add_loo(fit_brma))
  fit_brma <- suppressWarnings(add_waic(fit_brma))

  provenance_only <- fit_brma
  loo_target <- attr(
    provenance_only[["loo"]][["estimate"]],
    "RoBMA_target",
    exact = TRUE
  )
  waic_target <- attr(
    provenance_only[["waic"]][["estimate"]],
    "RoBMA_target",
    exact = TRUE
  )
  loo_target[["random_effects"]]           <- "conditioned"
  waic_target[["estimate_level_random"]] <- "marginalized"
  attr(provenance_only[["loo"]][["estimate"]], "RoBMA_target") <- loo_target
  attr(provenance_only[["waic"]][["estimate"]], "RoBMA_target") <- waic_target

  expect_no_error(loo(provenance_only))
  expect_no_error(waic(provenance_only))

  changed_target <- fit_brma
  target <- attr(
    changed_target[["loo"]][["estimate"]],
    "RoBMA_target",
    exact = TRUE
  )
  target[["target"]] <- "different_estimate_target"
  attr(changed_target[["loo"]][["estimate"]], "RoBMA_target") <- target
  expect_error(loo(changed_target), "current likelihood target")

  incomplete_target <- fit_brma
  target <- attr(
    incomplete_target[["waic"]][["estimate"]],
    "RoBMA_target",
    exact = TRUE
  )
  target[["conditioning_depth"]] <- NULL
  attr(incomplete_target[["waic"]][["estimate"]], "RoBMA_target") <- target
  expect_error(waic(incomplete_target), "incomplete RoBMA target metadata")
})


# ---------------------------------------------------------------------------- #
# Outcome type specific tests
# ---------------------------------------------------------------------------- #

test_that(".outcome_pdf.norm computes correct log-likelihood", {

  set.seed(123)

  yi  <- c(0.1, 0.2, 0.3)
  sei <- c(0.1, 0.1, 0.1)
  K   <- length(yi)
  S   <- 10

  mu_samples <- matrix(0.15, nrow = S, ncol = K)
  tau_within <- matrix(0.05, nrow = S, ncol = K)

  log_lik <- .outcome_pdf.norm(yi, mu_samples, tau_within, sei)

  expect_equal(dim(log_lik), c(S, K))

  # compute expected log-likelihood directly for first observation
  total_sd <- sqrt(0.05^2 + 0.1^2)
  expected_ll <- dnorm(0.1, mean = 0.15, sd = total_sd, log = TRUE)
  expect_equal(log_lik[1, 1], expected_ll, tolerance = 1e-10)
})

test_that(".outcome_cdf.norm keeps negative-direction tail precision", {

  yi         <- 10
  sei        <- 1
  mu_samples <- matrix(0, nrow = 1, ncol = 1)
  tau_within <- matrix(0, nrow = 1, ncol = 1)

  cdf_vals <- .outcome_cdf.norm(
    yi         = yi,
    mu_samples = mu_samples,
    tau_within = tau_within,
    sei        = sei,
    lower.tail = FALSE
  )

  expect_equal(cdf_vals[1, 1], stats::pnorm(yi, lower.tail = FALSE))
  expect_gt(cdf_vals[1, 1], 0)
})

test_that(".outcome_cdf.selnorm forwards lower.tail through step kernel", {

  yi         <- 3
  sei        <- 1
  mu_samples <- matrix(0, nrow = 1, ncol = 1)
  tau_within <- matrix(0, nrow = 1, ncol = 1)
  prior      <- BayesTools::prior_weightfunction(
    "one-sided", c(.05), BayesTools::wf_fixed(c(1, 1))
  )
  spec <- .selection_spec(
    priors           = list(outcome = list(bias = prior)),
    yi               = yi,
    sei              = sei,
    effect_direction = "positive"
  )
  selection_context <- c(spec, list(
    omega       = matrix(c(1, 1), nrow = 1),
    alpha       = 0,
    phack_kind  = 0L,
    kernel_mode = SELKERNEL_STEP
  ))

  cdf_vals <- .outcome_cdf.selnorm(
    yi                = yi,
    mu_samples        = mu_samples,
    tau_within        = tau_within,
    sei               = sei,
    selection_context = selection_context,
    lower.tail        = FALSE
  )

  expect_equal(cdf_vals[1, 1], stats::pnorm(yi, lower.tail = FALSE),
               tolerance = 1e-12)
  expect_gt(cdf_vals[1, 1], 0)
})

test_that("GLMM nuisance grids integrate over prior support", {

  prior_pi <- BayesTools::prior("beta", list(1, 1))
  pi_grid  <- .glmm_binom_logit_pi_grid(
    ai       = c(0L, 50L),
    ci       = c(50L, 0L),
    n1i      = c(100L, 100L),
    n2i      = c(100L, 100L),
    prior_pi = prior_pi,
    n_pi     = 17
  )

  expect_equal(pi_grid[["grid"]][, 1], pi_grid[["grid"]][, 2])
  expect_equal(sum(exp(pi_grid[["log_weights"]][, 1])), 1, tolerance = 1e-12)

  prior_phi <- BayesTools::prior("normal", list(0, 1))
  phi_grid  <- .glmm_pois_log_phi_grid(
    x1i       = c(0L, 100L),
    x2i       = c(100L, 0L),
    t1i       = c(10, 1000),
    t2i       = c(1000, 10),
    prior_phi = prior_phi,
    n_phi     = 17
  )

  expect_equal(phi_grid[["grid"]][, 1], phi_grid[["grid"]][, 2])
  expect_equal(sum(exp(phi_grid[["log_weights"]][, 1])), 1, tolerance = 1e-12)
})

test_that(".outcome_pdf.binom computes correct log-likelihood", {

  ai  <- c(10, 20)
  ci  <- c(15, 25)
  n1i <- c(100, 100)
  n2i <- c(100, 100)
  K   <- length(ai)
  S   <- 5

  mu_samples <- matrix(0.0, nrow = S, ncol = K) # log-OR = 0
  tau_within <- matrix(0.0, nrow = S, ncol = K)
  prior_pi   <- BayesTools::prior("beta", list(1, 1))

  log_lik <- .outcome_pdf.binom(
    ai         = ai,
    ci         = ci,
    n1i        = n1i,
    n2i        = n2i,
    mu_samples = mu_samples,
    tau_within = tau_within,
    prior_pi   = prior_pi,
    n_theta    = 3,
    n_pi       = 120
  )

  expect_equal(dim(log_lik), c(S, K))

  expected <- lchoose(n1i[1], ai[1]) + lchoose(n2i[1], ci[1]) +
    lbeta(ai[1] + ci[1] + 1, n1i[1] + n2i[1] - ai[1] - ci[1] + 1)

  expect_equal(log_lik[1, 1], expected, tolerance = 1e-8)
})

test_that(".outcome_pdf.binom handles boundary cell studies", {

  prior_pi   <- BayesTools::prior("beta", list(1, 1))
  mu_samples <- matrix(0, nrow = 2, ncol = 4)
  tau_within <- matrix(0, nrow = 2, ncol = 4)

  log_lik <- .outcome_pdf.binom(
    ai         = c(10, 0, 10, 0),
    ci         = c(10, 0, 0, 10),
    n1i        = c(10, 10, 10, 10),
    n2i        = c(10, 10, 10, 10),
    mu_samples = mu_samples,
    tau_within = tau_within,
    prior_pi   = prior_pi
  )

  expect_equal(dim(log_lik), c(2, 4))
  expect_true(all(is.finite(log_lik)))

  ai       <- c(10, 0, 10, 0)
  ci       <- c(10, 0, 0, 10)
  n1i      <- rep(10, 4)
  n2i      <- rep(10, 4)
  expected <- lchoose(n1i, ai) + lchoose(n2i, ci) +
    lbeta(ai + ci + 1, n1i + n2i - ai - ci + 1)

  expect_equal(log_lik[1, ], expected, tolerance = 1e-8)
  expect_equal(log_lik[2, ], expected, tolerance = 1e-8)
})

test_that(".outcome_pdf.binom matches R reference", {

  prior_pi <- BayesTools::prior("beta", list(1.5, 2.5))

  ai  <- c(3L, 0L, 8L, 12L)
  ci  <- c(2L, 4L, 0L, 10L)
  n1i <- c(20L, 18L, 15L, 30L)
  n2i <- c(22L, 17L, 16L, 31L)

  mu_samples <- matrix(
    c(-0.4, 0.1, 0.7, 0.2,
      -0.2, 0.3, 0.5, 0.0,
       0.1, 0.5, 0.3, -0.2),
    nrow = 3,
    ncol = 4
  )
  tau_within <- matrix(
    c(0.0, 0.1, 0.4, 0.2,
      0.2, 0.3, 0.1, 0.5,
      0.4, 0.2, 0.0, 0.3),
    nrow = 3,
    ncol = 4
  )

  expect_equal(
    .outcome_pdf.binom(ai, ci, n1i, n2i, mu_samples, tau_within, prior_pi, n_theta = 7, n_pi = 9),
    .outcome_pdf.binom_r(ai, ci, n1i, n2i, mu_samples, tau_within, prior_pi, n_theta = 7, n_pi = 9),
    tolerance = 1e-10
  )

  weights <- c(1, 0.5, 1.25, 2)
  weighted <- .outcome_pdf.binom(
    ai, ci, n1i, n2i, mu_samples, tau_within, prior_pi,
    weights = weights, n_theta = 7, n_pi = 9
  )

  expect_equal(
    weighted,
    .outcome_pdf.binom_r(
      ai, ci, n1i, n2i, mu_samples, tau_within, prior_pi,
      weights = weights, n_theta = 7, n_pi = 9
    ),
    tolerance = 1e-10
  )
  expect_false(isTRUE(all.equal(
    weighted,
    sweep(.outcome_pdf.binom(
      ai, ci, n1i, n2i, mu_samples, tau_within, prior_pi,
      n_theta = 7, n_pi = 9
    ), 2, weights, "*"),
    tolerance = 1e-8
  )))
})

test_that(".outcome_pdf.pois computes correct log-likelihood", {

  x1i <- c(10, 20)
  x2i <- c(15, 25)
  t1i <- c(100, 100)
  t2i <- c(100, 100)
  K   <- length(x1i)
  S   <- 5

  mu_samples <- matrix(0.0, nrow = S, ncol = K) # log-IRR = 0
  tau_within <- matrix(0.0, nrow = S, ncol = K)
  prior_phi  <- BayesTools::prior("normal", list(0, 1))

  log_lik <- .outcome_pdf.pois(
    x1i        = x1i,
    x2i        = x2i,
    t1i        = t1i,
    t2i        = t2i,
    mu_samples = mu_samples,
    tau_within = tau_within,
    prior_phi  = prior_phi,
    n_theta    = 3,
    n_phi      = 80
  )

  expect_equal(dim(log_lik), c(S, K))

  expected <- log(stats::integrate(
    function(phi) {
      lambda <- t1i[1] * exp(phi)
      exp(stats::dpois(x1i[1], lambda = lambda, log = TRUE) +
          stats::dpois(x2i[1], lambda = lambda, log = TRUE) +
          stats::dnorm(phi, log = TRUE))
    },
    lower   = -Inf,
    upper   = Inf,
    rel.tol = 1e-10
  )[["value"]])

  expect_equal(log_lik[1, 1], expected, tolerance = 1e-5)
})

test_that(".outcome_pdf.pois matches R reference", {

  prior_phi <- BayesTools::prior("normal", list(-1, 1.5))

  x1i <- c(0L, 3L, 10L, 12L)
  x2i <- c(1L, 0L, 8L, 9L)
  t1i <- c(12, 30, 45, 50)
  t2i <- c(10, 28, 43, 48)

  mu_samples <- matrix(
    c(-0.4, 0.1, 0.7, 0.2,
      -0.2, 0.3, 0.5, 0.0,
       0.1, 0.5, 0.3, -0.2),
    nrow = 3,
    ncol = 4
  )
  tau_within <- matrix(
    c(0.0, 0.1, 0.4, 0.2,
      0.2, 0.3, 0.1, 0.5,
      0.4, 0.2, 0.0, 0.3),
    nrow = 3,
    ncol = 4
  )

  expect_equal(
    .outcome_pdf.pois(x1i, x2i, t1i, t2i, mu_samples, tau_within, prior_phi, n_theta = 7, n_phi = 9),
    .outcome_pdf.pois_r(x1i, x2i, t1i, t2i, mu_samples, tau_within, prior_phi, n_theta = 7, n_phi = 9),
    tolerance = 1e-10
  )

  weights <- c(1, 0.5, 1.25, 2)
  weighted <- .outcome_pdf.pois(
    x1i, x2i, t1i, t2i, mu_samples, tau_within, prior_phi,
    weights = weights, n_theta = 7, n_phi = 9
  )

  expect_equal(
    weighted,
    .outcome_pdf.pois_r(
      x1i, x2i, t1i, t2i, mu_samples, tau_within, prior_phi,
      weights = weights, n_theta = 7, n_phi = 9
    ),
    tolerance = 1e-10
  )
  expect_false(isTRUE(all.equal(
    weighted,
    sweep(.outcome_pdf.pois(
      x1i, x2i, t1i, t2i, mu_samples, tau_within, prior_phi,
      n_theta = 7, n_phi = 9
    ), 2, weights, "*"),
    tolerance = 1e-8
  )))
})

test_that("native GLMM cluster likelihood matches R composition", {

  skip_if_not(.has_native_glmm_cluster(), "Native GLMM cluster kernels unavailable.")

  set.seed(2024)
  S <- 5
  K <- 5

  setup <- list(
    mu          = matrix(rnorm(S * K, 0, 0.25), nrow = S, ncol = K),
    tau_within  = matrix(runif(S * K, 0.05, 0.25), nrow = S, ncol = K),
    tau_between = matrix(runif(S * K, 0.02, 0.18), nrow = S, ncol = K),
    cluster     = list(a = c(1L, 3L), b = c(2L, 4L, 5L)),
    weights     = c(1, 0.5, 1.25, 2, 0.75)
  )

  bin_data <- list(outcome = data.frame(
    ai  = c(3L, 0L, 8L, 12L, 2L),
    ci  = c(2L, 4L, 0L, 10L, 1L),
    n1i = c(20L, 18L, 15L, 30L, 22L),
    n2i = c(22L, 17L, 16L, 31L, 21L)
  ))
  bin_priors <- list(outcome = list(
    pi = BayesTools::prior("beta", list(1.5, 2.5))
  ))

  expect_equal(
    .log_lik_cluster_glmm_native(
      setup        = setup,
      data         = bin_data,
      priors       = bin_priors,
      outcome_type = "bin",
      n_theta      = 5,
      n_gamma      = 5,
      n_pi         = 7
    ),
    .log_lik_cluster_glmm_r(
      setup        = setup,
      data         = bin_data,
      priors       = bin_priors,
      outcome_type = "bin",
      n_theta      = 5,
      n_gamma      = 5,
      n_pi         = 7
    ),
    tolerance = 1e-10
  )

  pois_data <- list(outcome = data.frame(
    x1i = c(0L, 3L, 10L, 12L, 1L),
    x2i = c(1L, 0L, 8L, 9L, 2L),
    t1i = c(12, 30, 45, 50, 25),
    t2i = c(10, 28, 43, 48, 24)
  ))
  pois_priors <- list(outcome = list(
    phi = BayesTools::prior("normal", list(-1, 1.5))
  ))

  expect_equal(
    .log_lik_cluster_glmm_native(
      setup        = setup,
      data         = pois_data,
      priors       = pois_priors,
      outcome_type = "pois",
      n_theta      = 5,
      n_gamma      = 5,
      n_phi        = 7
    ),
    .log_lik_cluster_glmm_r(
      setup        = setup,
      data         = pois_data,
      priors       = pois_priors,
      outcome_type = "pois",
      n_theta      = 5,
      n_gamma      = 5,
      n_phi        = 7
    ),
    tolerance = 1e-10
  )
})


# ---------------------------------------------------------------------------- #
# loo_weights and check_loo S3 tests
# ---------------------------------------------------------------------------- #

test_that("loo_weights and check_loo return stable diagnostics", {

  fit_brma <- fits[["bcg_meta-analysis"]]

  # check loo_weights
  weights <- loo_weights(fit_brma)
  expect_true(is.matrix(weights))
  expect_equal(dim(weights), dim(log_lik(fit_brma)))
  expect_equal(colSums(weights), rep(1, ncol(weights)), tolerance = 1e-10)

  # check check_loo (should not throw anything for this clean fit)
  expect_silent(check_loo(fit_brma))

  # simulate bad k
  fit_bad <- fit_brma
  fit_bad[["loo"]][["estimate"]][["diagnostics"]][["pareto_k"]][1] <- 0.8
  expect_warning(check_loo(fit_bad), "Some Pareto k values are high")
})

