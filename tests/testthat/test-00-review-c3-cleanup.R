test_that("prediction validates probabilities and explicit row admission", {

  object <- brma(yi = c(.1, -.2), sei = c(.3, .4), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE)
  samples <- cbind(mu = c(.1, .2), tau = .3)
  for (probs in list(NA_real_, NaN, Inf, -.01, 1.01, "bad")) {
    expect_error(predict(object, type = "terms", probs = probs, quiet = TRUE,
      .posterior_samples = samples), "probs", fixed = TRUE)
  }
  for (probs in list(c(0, 1), c(.1, .5, .9), NULL)) {
    expect_s3_class(predict(object, type = "terms", probs = probs, quiet = TRUE,
      .posterior_samples = samples), "brma_samples")
  }
  expect_error(predict(object, type = "cluster", quiet = TRUE,
    .posterior_samples = samples), paste0(
    "type = 'cluster' is only available for multilevel (3-level) models. ",
    "Use type = 'terms' for non-multilevel models."
  ), fixed = TRUE)
  for (newdata in list(data.frame(yi = numeric(), sei = numeric()),
                       list(yi = numeric(), sei = numeric()))) {
    expect_error(predict(object, newdata = newdata, type = "terms", quiet = TRUE,
      .posterior_samples = samples), "'newdata' must contain at least one row.", fixed = TRUE)
  }
  for (newdata in list(list(c(.1, .2)), list(yi = c(.1, .2), c(.3, .4)))) {
    expect_error(predict(object, newdata = newdata, type = "terms", quiet = TRUE,
      .posterior_samples = samples), "'newdata' list must have a name for every column.", fixed = TRUE)
  }
  expect_s3_class(predict(object, newdata = list(yi = c(.1, .2)), type = "terms",
    quiet = TRUE, .posterior_samples = samples), "brma_samples")
})


test_that("best-weight newdata identifies its missing publication-group column", {

  dat <- data.frame(yi = c(.4, -.2), sei = c(.8, 1.1), paper = c("a", "a"))
  object <- bselmodel(yi = yi, sei = sei, data = dat, measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE,
    prior_bias = BayesTools::prior_weightfunction("one-sided", .5,
      BayesTools::wf_fixed(c(1, .2)), model = BayesTools::selection_model(
        group = "paper", weight_rule = "best")))
  samples <- cbind(mu = c(.1, .2), tau = .3)
  expect_error(predict(object, newdata = dat[c("yi", "sei")], type = "response",
    quiet = TRUE, .posterior_samples = samples), paste0(
    "The 'newdata' must contain the publication-group column 'paper' ",
    "required by 'weight_rule = \"best\"'."
  ), fixed = TRUE)
})
