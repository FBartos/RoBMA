context("Summary heterogeneity helpers")

.heterogeneity_correlation_object <- function(random, posterior_samples) {

  dat <- data.frame(
    yi    = c(.1, .2, .3, .4),
    time  = c(1, 2, 1, 2),
    study = c("s1", "s1", "s2", "s2"),
    site  = c("a", "b", "a", "b")
  )
  object <- brma.mv(
    yi                        = yi,
    vi                        = rep(.04, 4),
    random                    = random,
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  design <- RoBMA:::.fitted_formula_design(object, "mu", required = TRUE)
  fit    <- coda::mcmc.list(coda::mcmc(posterior_samples))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list")     <- design[["prior_list"]]
  attr(fit, "formula_design") <- list(mu = design)
  attr(fit, "parameter_map") <- BayesTools:::.bt_build_parameter_map(
    columns        = colnames(posterior_samples),
    prior_list     = design[["prior_list"]],
    formula_design = list(mu = design)
  )
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  object[["fit"]] <- BayesTools:::.bt_attach_fit_contract(fit)

  return(object)
}


test_that("homogeneous structured heterogeneity includes public correlations", {

  posterior_samples <- cbind(
    mu_intercept       = c(.1, .2, .3, .4),
    mu__xREx__study_sd  = c(.2, .4, .6, .8),
    mu__xREx__study_rho = c(-.3, .1, .4, .8)
  )
  formulas <- list(~ ar(time | study), ~ cs(time | study), ~ car(time | study))
  for (random in formulas) {
    draws <- posterior_samples
    if (identical(random, formulas[[3L]])) {
      draws[, "mu__xREx__study_rho"] <- c(.1, .2, .4, .8)
    }
    object       <- .heterogeneity_correlation_object(random, draws)
    probs        <- c(.1, .9)
    result       <- summary_heterogeneity(object, probs = probs)
    estimates    <- result[["estimates"]]
    expected_cor <- draws[, "mu__xREx__study_rho"]

    expect_identical(rownames(estimates), c("sd", "var", "cor"))
    expect_equal(estimates["cor", "Mean"], mean(expected_cor))
    expect_equal(estimates["cor", "Median"], median(expected_cor))
    expect_equal(
      as.numeric(estimates["cor", c("0.1", "0.9")]),
      unname(quantile(expected_cor, probs))
    )
    expect_equal(estimates["sd", "Mean"], mean(draws[, "mu__xREx__study_sd"]))
    expect_equal(estimates["var", "Mean"], mean(draws[, "mu__xREx__study_sd"]^2))
    expect_identical(
      summary_heterogeneity(object, component = "study", probs = probs),
      result
    )
    expect_equal(
      summary_heterogeneity(object, component = "total", probs = probs)[["estimates"]][
        "cor", "Mean"
      ],
      mean(expected_cor)
    )
    expect_true(any(grepl("^cor +", capture.output(print(result)))))
    frame <- as.data.frame(result)
    expect_identical(frame[["component"]], rep("study", 3L))
    expect_identical(frame[["parameter"]], c("sd", "var", "cor"))
    expect_identical(tail(names(frame), 2L), c("CI_0.1", "CI_0.9"))
    expect_identical(data.frame(result), frame)
  }
})


test_that("correlations stay with their component and out of additive totals", {

  draws <- cbind(
    mu_intercept                              = c(.1, .2, .3, .4),
    mu__xREx__study_sd                         = sqrt(.25) * c(1, 2, 3, 4),
    mu__xREx__study_rho                        = c(-.3, .1, .4, .8),
    mu__xREx__site_sd                          = sqrt(.75) * c(1, 2, 3, 4),
    mu__xREx__site_rho                         = c(.2, .3, .4, .5),
    mu__xRE_ALLOCx_heterogeneity__allocation_sd = c(1, 2, 3, 4),
    "mu__xRE_ALLOCx_heterogeneity__weight[1]"   = rep(.25, 4),
    "mu__xRE_ALLOCx_heterogeneity__weight[2]"   = rep(.75, 4)
  )
  object <- .heterogeneity_correlation_object(
    list(study = ~ ar(time | study), site = ~ cs(time | site)),
    draws
  )
  result <- summary_heterogeneity(object)
  for (block in c("study", "site")) {
    label     <- paste0(block, ": cor")
    estimates <- result[[block]][["estimates"]]
    expect_identical(rownames(estimates), c("sd", "var", label))
    expect_equal(
      estimates[label, "Mean"],
      mean(draws[, paste0("mu__xREx__", block, "_rho")])
    )
    expect_identical(
      summary_heterogeneity(object, component = block),
      result[[block]]
    )
  }
  expect_identical(
    rownames(summary_heterogeneity(object, component = "total")[["estimates"]]),
    c("sd_total", "var_total")
  )
  frame <- as.data.frame(result)
  expect_identical(data.frame(result), frame)
  expect_identical(
    frame[["component"]][grepl(": cor$", frame[["parameter"]])],
    c("study", "site")
  )
})


test_that("allocation correlations are not duplicated or invented", {

  draws <- cbind(
    mu_intercept                              = c(.1, .2, .3, .4),
    mu__xREx__study_rho                        = c(-.3, .1, .4, .8),
    mu__xRE_ALLOCx_heterogeneity__allocation_sd = c(1, 2, 3, 4),
    "mu__xRE_ALLOCx_heterogeneity__weight[1]"   = rep(.25, 4),
    "mu__xRE_ALLOCx_heterogeneity__weight[2]"   = rep(.75, 4)
  )
  object <- .heterogeneity_correlation_object(~ har(time | study), draws)
  result <- summary_heterogeneity(object)
  expect_equal(sum(rownames(result[["estimates"]]) == "cor"), 1L)
  expect_equal(
    result[["estimates"]]["cor", "Mean"],
    mean(draws[, "mu__xREx__study_rho"])
  )
  expect_equal(result[["estimates"]]["sd_common", "Mean"], 2.5)
  expect_identical(summary_heterogeneity(object, component = "study"), result)

  independent <- .heterogeneity_correlation_object(
    ~ 1 | study,
    cbind(mu_intercept = c(.1, .2, .3, .4),
          mu__xREx__study_intercept = c(.2, .4, .6, .8))
  )
  expect_identical(
    rownames(summary_heterogeneity(independent)[["estimates"]]),
    c("sd", "var")
  )
})


test_that("heterogeneity summaries coerce to component-aware data frames", {

  make_summary <- function(component, offset) {
    estimates <- BayesTools::ensemble_estimates_table(
      samples = list(tau = offset + 1:20, tau2 = (offset + 1:20)^2),
      parameters = c("tau", "tau2"),
      probs = c(.025, .975)
    )
    structure(
      list(estimates = estimates, component = component),
      class = "summary_heterogeneity.brma"
    )
  }
  study <- make_summary("study", 0)
  site  <- make_summary("site", 1)
  summaries <- structure(
    list(study = study, site = site),
    class = c("summary_heterogeneity.brma_list", "list")
  )

  study_frame <- as.data.frame(study)
  long_frame  <- as.data.frame(summaries)
  list_frames <- as.data.frame(summaries, format = "list")

  expect_identical(names(study_frame)[1:2], c("component", "parameter"))
  expect_identical(study_frame[["component"]], rep("study", 2L))
  expect_identical(study_frame[["parameter"]], c("tau", "tau2"))
  expect_identical(names(long_frame)[1:2], c("component", "parameter"))
  expect_setequal(long_frame[["component"]], c("study", "site"))
  expect_true(all(vapply(list_frames, is.data.frame, logical(1))))
  expect_identical(data.frame(study), data.frame(study_frame))
  expect_identical(data.frame(summaries), data.frame(long_frame))
})

test_that("rho allocation retains endpoints and rejects invalid values", {

  tau <- matrix(.3, nrow = 3L, ncol = 2L)
  out <- RoBMA:::.heterogeneity_components(
    tau_total     = tau,
    rho           = c(0, .5, 1),
    is_multilevel = TRUE
  )

  expect_equal(out[["rho"]], c(0, .5, 1))
  expect_equal(out[["tau_within"]][1, ], c(.3, .3))
  expect_equal(out[["tau_between"]][1, ], c(0, 0))
  expect_equal(out[["tau_within"]][3, ], c(0, 0))
  expect_equal(out[["tau_between"]][3, ], c(.3, .3))
  expect_equal(
    out[["tau_within"]]^2 + out[["tau_between"]]^2,
    tau^2
  )

  expect_error(
    RoBMA:::.heterogeneity_components(
      tau[1, , drop = FALSE],
      -1e-12,
      TRUE
    ),
    "within \\[0, 1\\]"
  )
  expect_error(
    RoBMA:::.heterogeneity_components(
      tau[1, , drop = FALSE],
      1 + 1e-12,
      TRUE
    ),
    "within \\[0, 1\\]"
  )

  posterior <- matrix(
    c(.2, .4),
    ncol     = 1L,
    dimnames = list(NULL, "tau")
  )
  fixed <- RoBMA:::.evaluate.brma.tau(
    fit               = NULL,
    scale_data        = NULL,
    scale_formula     = NULL,
    scale_priors      = NULL,
    is_scale          = FALSE,
    is_multilevel     = TRUE,
    K                 = 1L,
    posterior_samples = posterior,
    fixed_rho         = 1
  )
  expect_equal(fixed[["tau_within"]], matrix(0, nrow = 2L))
  expect_equal(as.numeric(fixed[["tau_between"]]), c(.2, .4))
})


test_that("fixed rho summaries use evaluated allocation samples", {

  observed_rho <- NULL
  testthat::local_mocked_bindings(
    .outcome_data_vi = function(object) c(1, 1),
    .get_model_matrix = function(object) matrix(1, nrow = 2L, ncol = 1L),
    .get_posterior_samples = function(fit) matrix(
      c(.2, .4),
      ncol     = 1L,
      dimnames = list(NULL, "tau")
    ),
    .evaluate.brma.tau = function(...) list(
      tau_within  = matrix(0, nrow = 2L, ncol = 2L),
      tau_between = matrix(c(.2, .2, .4, .4), nrow = 2L),
      rho         = c(1, 1)
    ),
    .summary_heterogeneity_samples = function(
        tau_within_samples, tau_between_samples, rho_samples, ...) {
      observed_rho <<- rho_samples
      list(rho = rho_samples)
    },
    .package = "RoBMA"
  )

  data <- list(scale = NULL)
  attr(data, "cluster") <- TRUE
  attr(data, "scale")   <- FALSE
  object <- structure(list(
    fit    = NULL,
    data   = data,
    priors = list(
      outcome = list(
        rho = BayesTools::prior(
          "point",
          parameters = list(location = 1)
        )
      ),
      scale = NULL
    )
  ), class = "brma")

  summary_heterogeneity.brma(object)

  expect_identical(observed_rho, c(1, 1))
})

test_that("scale heterogeneity summaries aggregate variances before tau", {

  tau_within <- matrix(
    c(1, 3,
      2, 4),
    nrow = 2,
    byrow = TRUE
  )

  samples <- RoBMA:::.summary_heterogeneity_samples(
    tau_within_samples  = tau_within,
    tau_between_samples = matrix(0, nrow = 2, ncol = 2),
    v_tilde             = 1,
    is_multilevel       = FALSE
  )

  expected_tau2 <- rowMeans(tau_within^2)
  old_tau2      <- rowMeans(tau_within)^2

  expect_equal(samples[["tau2"]], expected_tau2)
  expect_equal(samples[["tau"]], sqrt(expected_tau2))
  expect_false(isTRUE(all.equal(samples[["tau2"]], old_tau2)))
  expect_equal(samples[["I2"]], rowMeans(100 * tau_within^2 / (tau_within^2 + 1)))
  expect_equal(samples[["H2"]], rowMeans(tau_within^2 + 1))
})

test_that("multilevel scale heterogeneity partitions variance and I2", {

  tau_within <- matrix(
    c(1, 2,
      3, 4),
    nrow = 2,
    byrow = TRUE
  )
  tau_between <- matrix(
    c(2, 1,
      1, 3),
    nrow = 2,
    byrow = TRUE
  )

  samples <- RoBMA:::.summary_heterogeneity_samples(
    tau_within_samples  = tau_within,
    tau_between_samples = tau_between,
    rho_samples         = c(0.8, 0.2),
    v_tilde             = 2,
    is_multilevel       = TRUE
  )

  sigma2_within  <- tau_within^2
  sigma2_between <- tau_between^2
  sigma2_total   <- sigma2_within + sigma2_between
  denominator    <- sigma2_total + 2

  expect_equal(samples[["tau2 [within]"]], rowMeans(sigma2_within))
  expect_equal(samples[["tau2 [between]"]], rowMeans(sigma2_between))
  expect_equal(samples[["tau2"]], rowMeans(sigma2_total))
  expect_equal(samples[["rho"]], c(0.8, 0.2))
  expect_equal(samples[["I2"]], rowMeans(100 * sigma2_total / denominator))
  expect_equal(samples[["I2 [within]"]], rowMeans(100 * sigma2_within / denominator))
  expect_equal(samples[["I2 [between]"]], rowMeans(100 * sigma2_between / denominator))
  expect_equal(samples[["I2"]], samples[["I2 [within]"]] + samples[["I2 [between]"]])
})


test_that("heterogeneity summaries reject invalid rho without projection", {

  tau <- matrix(0, nrow = 2L, ncol = 1L)
  unidentified <- RoBMA:::.summary_heterogeneity_samples(
    tau_within_samples  = tau,
    tau_between_samples = tau,
    v_tilde             = 1,
    is_multilevel       = TRUE
  )

  expect_true(all(is.na(unidentified[["rho"]])))
  expect_error(
    RoBMA:::.summary_heterogeneity_samples(
      tau_within_samples  = tau,
      tau_between_samples = tau,
      rho_samples         = c(-1e-12, 1),
      v_tilde             = 1,
      is_multilevel       = TRUE
    ),
    "within \\[0, 1\\]"
  )
})
