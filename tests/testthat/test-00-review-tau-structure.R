.review_tau_random_formula_fit <- function() {

  data <- data.frame(yi = c(.1, .2, -.1, .3), vi = rep(.04, 4L),
    study = factor(c("a", "a", "b", "b")))
  object <- brma.mv(yi = yi, V = diag(data$vi), data = data, random = ~ 1 | study,
    measure = "GEN",
    prior_effect = BayesTools::prior("normal", list(0, 1)),
    prior_heterogeneity = BayesTools::prior("normal", list(0, 1),
      truncation = list(lower = 0, upper = Inf)),
    only_priors = TRUE)
  .object_formula_fit(object)
}


test_that("missing tau follows declared formula structure instead of column spelling", {

  formula_fit <- .review_tau_random_formula_fit()
  design <- BayesTools::JAGS_formula_design(formula_fit)
  expect_gt(length(design[["mu"]][["random_effects"]]), 0L)
  clean <- matrix(.2, 3L, 1L, dimnames = list(NULL, "mu_intercept"))
  spoofed <- clean
  colnames(spoofed) <- "unrelated__xRE_name"
  for (pooled in c(FALSE, TRUE)) {
    evaluate <- if (pooled) .evaluate.brma.pooled_tau else .evaluate.brma.tau
    call <- list(fit = formula_fit, scale_data = NULL, scale_formula = NULL,
      scale_priors = NULL, is_scale = FALSE, is_multilevel = FALSE,
      posterior_samples = clean)
    if (!pooled) call[["K"]] <- 2L
    declared <- tryCatch(do.call(evaluate, call), error = identity)
    expect_false(inherits(declared, "error"))
    if (!inherits(declared, "error")) {
      expect_identical(declared[["tau_total"]], matrix(0, 3L, if (pooled) 1L else 2L))
    }
    call["fit"] <- list(NULL)
    call[["posterior_samples"]] <- spoofed
    expect_error(do.call(evaluate, call), "Missing posterior tau columns.", fixed = TRUE)
  }
})


test_that("explicit random structure overrides the formula fallback in both directions", {

  formula_fit <- .review_tau_random_formula_fit()
  samples <- matrix(.2, 3L, 1L, dimnames = list(NULL, "unrelated__xRE_name"))
  for (pooled in c(FALSE, TRUE)) {
    evaluate <- if (pooled) .evaluate.brma.pooled_tau else .evaluate.brma.tau
    call <- list(fit = formula_fit, scale_data = NULL, scale_formula = NULL,
      scale_priors = NULL, is_scale = FALSE, is_multilevel = FALSE,
      posterior_samples = samples, is_random = FALSE)
    if (!pooled) call[["K"]] <- 2L
    expect_error(do.call(evaluate, call), "Missing posterior tau columns.", fixed = TRUE)
    call["fit"] <- list(NULL)
    call[["is_random"]] <- TRUE
    result <- tryCatch(do.call(evaluate, call), error = identity)
    expect_false(inherits(result, "error"))
    if (!inherits(result, "error")) {
      expect_identical(result[["tau_total"]], matrix(0, 3L, if (pooled) 1L else 2L))
    }
  }
})
