test_that("automatic estimate marginalization honors sampled mean-centered requests", {

  for (repeated in c(FALSE, TRUE)) {
    dat <- data.frame(
      yi = if (repeated) c(0.1, 0.2, 0.3, 0.4) else c(0.1, 0.2, 0.3),
      study = factor(if (repeated) c("a", "a", "b", "b") else c("a", "b", "c"))
    )
    for (mode in c("mean_centered", "centered", "noncentered")) {
      expect_silent(object <- brma.mv(
        yi = yi, V = diag(0.01, nrow(dat)), random = ~ 1 | study,
        data = dat, measure = "GEN", prior_unit_information_sd = 1,
        prior_effect = BayesTools::prior("normal", list(0, 1)),
        prior_heterogeneity = BayesTools::prior_random(
          sd = BayesTools::prior("gamma", list(shape = 2, rate = 2)),
          study = BayesTools::random_block(parameterization = mode)),
        standardize_continuous_predictors = FALSE, only_priors = TRUE, silent = TRUE
      ))
      term <- .fitted_formula_design(object, "mu", required = TRUE)[["random_effects"]][[1L]]
      if (mode == "mean_centered") {
        expect_identical(term[["compile_mode"]], "sampled")
        expect_identical(term[["parameterization_resolved"]], "mean_centered")
        expect_type(term[["mean_translation"]], "list")
        expect_false(.is_marginalized_estimate_level_candidate(term, nrow(dat)))
      } else {
        expect_identical(term[["compile_mode"]], if (repeated) "sampled" else "marginalized")
      }
    }
  }
})


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
    paste0("Mean-centered parameterization is unavailable for integrated selection source 'study'. ",
      "Use 'random_block(parameterization = \"centered\")' or ",
      "'random_block(parameterization = \"noncentered\")' in 'prior_heterogeneity'."))
  expect_null(conditionCall(condition))
})
