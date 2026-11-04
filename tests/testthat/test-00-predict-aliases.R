test_that("prediction abbreviations preserve their canonical conditioning target", {

  object <- brma(yi = c(-1, 2), sei = c(1, 1), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE)
  posterior <- cbind(mu = c(0, 0, 0), tau = c(1, 1, 1))
  aliases <- c(b = "blup", es = "estimate", ef = "effect",
    o = "outcome", r = "response", m = "marginal")
  for (alias in names(aliases)) {
    canonical <- aliases[[alias]]
    withr::local_seed(927)
    reference <- predict(object, type = canonical, quiet = TRUE,
      .posterior_samples = posterior)
    set.seed(927)
    actual <- predict(object, type = alias, quiet = TRUE,
      .posterior_samples = posterior)
    expect_equal(actual, reference, info = canonical)
  }
  expect_equal(as.numeric(predict(object, type = "b", quiet = TRUE,
    .posterior_samples = posterior)), rep(c(-.5, 1), each = 3L))
  for (alias in c("b", "blup")) {
    error <- tryCatch(predict(object, type = alias,
      conditioning_depth = "marginal", quiet = TRUE,
      .posterior_samples = posterior), error = identity)
    expect_s3_class(error, "error")
    expect_identical(if (inherits(error, "error")) conditionMessage(error) else NA_character_,
      "type = 'blup' fixes conditioning_depth = 'estimate'.")
  }
  for (alias in c("cl", "cluster")) {
    error <- tryCatch(predict(object, type = alias,
      conditioning_depth = "estimate", quiet = TRUE,
      .posterior_samples = posterior), error = identity)
    expect_s3_class(error, "error")
    expect_identical(if (inherits(error, "error")) conditionMessage(error) else NA_character_,
      "type = 'cluster' fixes conditioning_depth = 'cluster'.")
  }
  multilevel <- brma(yi = c(-1, 2), sei = c(1, 1), cluster = c("a", "a"),
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE)
  multilevel_posterior <- cbind(posterior, rho = rep(.5, 3L), `gamma[1]` = rep(0, 3L))
  set.seed(927)
  reference <- predict(multilevel, type = "cluster", quiet = TRUE,
    .posterior_samples = multilevel_posterior)
  set.seed(927)
  actual <- predict(multilevel, type = "cl", quiet = TRUE,
    .posterior_samples = multilevel_posterior)
  expect_equal(actual, reference)
})
