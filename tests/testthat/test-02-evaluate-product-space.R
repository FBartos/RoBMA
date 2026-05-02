context("Evaluate product-space posterior structure")

source(testthat::test_path("common-functions.R"))

skip_if_no_fits()
fit_names <- list_fits()
fits      <- lazy_fits(fit_names, validate = FALSE)

test_that("RoBMA mixed-bias branch PDF, CDF, and RNG use selected-normal branches", {

  name <- "dat.lehmann2018_RoBMA"
  skip_if_missing_fits(name)

  object            <- fits[[name]]
  setup             <- .estimate_likelihood_setup.brma(object)
  posterior_samples <- setup[["posterior_samples"]]
  bias_indicator    <- .extract_bias_indicator(object, posterior_samples = posterior_samples)
  use_normal        <- .extract_use_normal(object, posterior_samples = posterior_samples)

  normal_rows   <- head(which(use_normal), 100)
  weighted_rows <- head(which(!use_normal), 100)
  if (length(normal_rows) == 0L || length(weighted_rows) == 0L) {
    skip("Cached RoBMA fit lacks both normal and weighted publication-bias branches.")
  }

  selected_rows <- sort(c(normal_rows, weighted_rows))
  selection_context <- .selection_context(
    object            = object,
    posterior_samples = posterior_samples
  )
  selected_context <- .selection_context_subset_rows(selection_context, selected_rows)

  log_lik <- .outcome_pdf.selnorm(
    yi                = setup[["yi"]],
    mu_samples        = setup[["mu"]][selected_rows, , drop = FALSE],
    tau_within        = setup[["tau_within"]][selected_rows, , drop = FALSE],
    sei               = setup[["sei"]],
    selection_context = selected_context
  )
  cdf_vals <- .outcome_cdf.selnorm(
    yi                = setup[["yi"]],
    mu_samples        = setup[["mu"]][selected_rows, , drop = FALSE],
    tau_within        = setup[["tau_within"]][selected_rows, , drop = FALSE],
    sei               = setup[["sei"]],
    selection_context = selected_context
  )

  expect_equal(dim(log_lik), c(length(selected_rows), setup[["K"]]))
  expect_equal(dim(cdf_vals), c(length(selected_rows), setup[["K"]]))
  expect_true(all(is.finite(log_lik)))
  expect_true(all(is.finite(cdf_vals)))
  expect_true(all(cdf_vals > 0 & cdf_vals < 1))

  normal_branch   <- sort(unique(bias_indicator[normal_rows]))[1]
  weighted_branch <- sort(unique(bias_indicator[weighted_rows]))[1]
  normal_branch_rows   <- normal_rows[bias_indicator[normal_rows] == normal_branch]
  weighted_branch_rows <- weighted_rows[bias_indicator[weighted_rows] == weighted_branch]

  normal_pdf <- .outcome_pdf.norm(
    yi         = if (setup[["effect_direction"]] == "negative") -setup[["yi"]] else setup[["yi"]],
    mu_samples = if (setup[["effect_direction"]] == "negative") {
      -setup[["mu"]][normal_branch_rows, , drop = FALSE]
    } else {
      setup[["mu"]][normal_branch_rows, , drop = FALSE]
    },
    tau_within = setup[["tau_within"]][normal_branch_rows, , drop = FALSE],
    sei        = setup[["sei"]]
  )
  normal_cdf <- .outcome_cdf.norm(
    yi         = if (setup[["effect_direction"]] == "negative") -setup[["yi"]] else setup[["yi"]],
    mu_samples = if (setup[["effect_direction"]] == "negative") {
      -setup[["mu"]][normal_branch_rows, , drop = FALSE]
    } else {
      setup[["mu"]][normal_branch_rows, , drop = FALSE]
    },
    tau_within = setup[["tau_within"]][normal_branch_rows, , drop = FALSE],
    sei        = setup[["sei"]],
    lower.tail = if (setup[["effect_direction"]] == "negative") FALSE else TRUE
  )

  expect_equal(
    log_lik[match(normal_branch_rows, selected_rows), , drop = FALSE],
    normal_pdf,
    tolerance = 1e-10
  )
  expect_equal(
    cdf_vals[match(normal_branch_rows, selected_rows), , drop = FALSE],
    normal_cdf,
    tolerance = 1e-10
  )

  weighted_pdf_norm <- .outcome_pdf.norm(
    yi         = if (setup[["effect_direction"]] == "negative") -setup[["yi"]] else setup[["yi"]],
    mu_samples = if (setup[["effect_direction"]] == "negative") {
      -setup[["mu"]][weighted_branch_rows, , drop = FALSE]
    } else {
      setup[["mu"]][weighted_branch_rows, , drop = FALSE]
    },
    tau_within = setup[["tau_within"]][weighted_branch_rows, , drop = FALSE],
    sei        = setup[["sei"]]
  )
  expect_false(isTRUE(all.equal(
    log_lik[match(weighted_branch_rows, selected_rows), , drop = FALSE],
    weighted_pdf_norm,
    tolerance = 1e-8
  )))

  set.seed(11)
  rng_normal <- .outcome_rng.selnorm(
    mu_samples        = setup[["mu"]][normal_branch_rows, , drop = FALSE],
    tau_within        = setup[["tau_within"]][normal_branch_rows, , drop = FALSE],
    sei               = setup[["sei"]],
    selection_context = .selection_context_subset_rows(selection_context, normal_branch_rows)
  )
  set.seed(11)
  rng_normal_expected <- .outcome_rng.norm(
    mu_samples = setup[["mu"]][normal_branch_rows, , drop = FALSE],
    tau_within = setup[["tau_within"]][normal_branch_rows, , drop = FALSE],
    sei        = setup[["sei"]]
  )
  expect_equal(dim(rng_normal), dim(rng_normal_expected))
  expect_equal(as.vector(rng_normal), as.vector(rng_normal_expected),
               tolerance = 1e-12)

  set.seed(12)
  rng_mixed <- .outcome_rng.selnorm(
    mu_samples        = setup[["mu"]][selected_rows, , drop = FALSE],
    tau_within        = setup[["tau_within"]][selected_rows, , drop = FALSE],
    sei               = setup[["sei"]],
    selection_context = selected_context
  )

  expect_equal(dim(rng_mixed), c(length(selected_rows), setup[["K"]]))
  expect_true(all(is.finite(rng_mixed)))
})

.expect_indicator_for_prior <- function(posterior_samples, prior, column, info) {

  if (is.null(prior) || !BayesTools::is.prior.mixture(prior)) {
    return(invisible(TRUE))
  }

  expect_true(column %in% colnames(posterior_samples),
              info = paste(info, "missing posterior indicator", column))
  if (!column %in% colnames(posterior_samples)) {
    return(invisible(FALSE))
  }

  values <- posterior_samples[, column]

  expect_true(all(is.finite(values)),
              info = paste(info, column, "contains non-finite values"))
  expect_true(all(values == as.integer(values)),
              info = paste(info, column, "contains non-integer values"))
  expect_true(all(values >= 1 & values <= length(prior)),
              info = paste(info, column, "outside prior component range"))
}

test_that("product-space posterior indicators have expected columns and valid ranges", {

  product_names <- c(
    "dat.lehmann2018_BMA.norm",
    "dat.lehmann2018_BMA.norm_mods",
    "dat.lehmann2018_BMA.norm_scale",
    "bcg_BMA.glmm",
    "bcg_BMA.glmm_3lvl_location_scale",
    "dat.lehmann2018_RoBMA",
    "dat.lehmann2018_RoBMA_3lvl_mods_scale"
  )
  skip_if_missing_fits(product_names)

  for (name in product_names) {

    object            <- fits[[name]]
    posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
    priors            <- object[["priors"]]

    .expect_indicator_for_prior(
      posterior_samples = posterior_samples,
      prior             = priors[["outcome"]][["mu"]],
      column            = "mu_indicator",
      info              = name
    )
    .expect_indicator_for_prior(
      posterior_samples = posterior_samples,
      prior             = priors[["outcome"]][["tau"]],
      column            = "tau_indicator",
      info              = name
    )
    .expect_indicator_for_prior(
      posterior_samples = posterior_samples,
      prior             = priors[["outcome"]][["rho"]],
      column            = "rho_indicator",
      info              = name
    )
    .expect_indicator_for_prior(
      posterior_samples = posterior_samples,
      prior             = priors[["outcome"]][["bias"]],
      column            = "bias_indicator",
      info              = name
    )

    for (term in names(priors[["mods"]])) {
      column <- paste0(
        BayesTools::JAGS_parameter_names(term, formula_parameter = "mu"),
        "_indicator"
      )
      .expect_indicator_for_prior(
        posterior_samples = posterior_samples,
        prior             = priors[["mods"]][[term]],
        column            = column,
        info              = name
      )
    }

    for (term in names(priors[["scale"]])) {
      column <- paste0(
        BayesTools::JAGS_parameter_names(term, formula_parameter = "log_tau"),
        "_indicator"
      )
      .expect_indicator_for_prior(
        posterior_samples = posterior_samples,
        prior             = priors[["scale"]][[term]],
        column            = column,
        info              = name
      )
    }
  }
})

test_that("RoBMA inactive bias branches expose zeros or neutral weights", {

  product_names <- c("dat.lehmann2018_RoBMA", "dat.lehmann2018_RoBMA_3lvl_mods_scale")
  skip_if_missing_fits(product_names)

  for (name in product_names) {

    object            <- fits[[name]]
    posterior_samples <- suppressWarnings(coda::as.mcmc(object[["fit"]]))
    priors_bias       <- object[["priors"]][["outcome"]][["bias"]]
    if (is.null(priors_bias) || !"bias_indicator" %in% colnames(posterior_samples)) {
      next
    }
    if (!BayesTools::is.prior.mixture(priors_bias)) {
      priors_bias <- list(priors_bias)
    }

    bias_indicator      <- as.integer(posterior_samples[, "bias_indicator"])
    branch_is_PET       <- vapply(priors_bias, BayesTools::is.prior.PET, logical(1))
    branch_is_PEESE     <- vapply(priors_bias, BayesTools::is.prior.PEESE, logical(1))
    branch_is_weightfun <- vapply(priors_bias, BayesTools::is.prior.weightfunction, logical(1))

    if ("PET" %in% colnames(posterior_samples)) {
      inactive <- which(!branch_is_PET[bias_indicator])
      expect_equal(as.numeric(posterior_samples[inactive, "PET"]), rep(0, length(inactive)),
                   tolerance = 1e-12, info = paste(name, "inactive PET branch"))
    }
    if ("PEESE" %in% colnames(posterior_samples)) {
      inactive <- which(!branch_is_PEESE[bias_indicator])
      expect_equal(as.numeric(posterior_samples[inactive, "PEESE"]), rep(0, length(inactive)),
                   tolerance = 1e-12, info = paste(name, "inactive PEESE branch"))
    }

    omega_cols <- grep("^omega(\\[|$)", colnames(posterior_samples), value = TRUE)
    if (length(omega_cols) > 0L) {
      inactive <- which(!branch_is_weightfun[bias_indicator])
      omega_inactive <- as.matrix(posterior_samples[inactive, omega_cols, drop = FALSE])
      expect_equal(dim(omega_inactive), c(length(inactive), length(omega_cols)),
                   info = paste(name, "inactive weightfunction branch"))
      expect_equal(as.vector(omega_inactive), rep(1, length(omega_inactive)),
                   tolerance = 1e-12,
                   info = paste(name, "inactive weightfunction branch"))
    }
  }
})
