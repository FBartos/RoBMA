test_that("mean-centered selection intercepts use retained source metadata", {

  dat <- data.frame(yi = c(.1, .2, -.1, .05),
    study = factor(rep(c("a", "b"), each = 2L)), esid = factor(rep(1:2, 2L)))
  create <- function(mode) {

    bselmodel.mv(yi = dat$yi, V = diag(.04, 4L), random = ~ 1 | study / esid,
      data = dat, measure = "SMD", only_priors = TRUE,
      prior_effect = BayesTools::prior("t", list(0, .7, 3)),
      prior_heterogeneity = BayesTools::prior_random(
        sd = BayesTools::prior("point", list(.25)),
        study = BayesTools::random_block(parameterization = "mean_centered")),
      selection = selection_model(other_random_effects = mode, group = "study"))
  }
  retained <- create("condition")
  sources <- .data_selection_model(retained$data)$sources$random
  study_source <- sources[[which(vapply(sources, function(source) source$name == "study", logical(1L)))]]
  expect_true(study_source$retained)
  terms <- retained$formula_design$mu$random_effects
  study <- terms[[which(vapply(terms, function(term) term$block_name == "study", logical(1L)))]]
  expect_identical(study$parameterization_resolved, "mean_centered")
  condition <- tryCatch(create("integrate"), error = identity)
  expect_identical(conditionMessage(condition),
    "Mean-centered parameterization is unavailable for integrated selection source 'study'.")
  expect_null(conditionCall(condition))
})
