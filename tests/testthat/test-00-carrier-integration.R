.carrier_integration_fit_cache <- new.env(parent = emptyenv())

.carrier_integration_fits <- function(){

  if(!is.null(.carrier_integration_fit_cache$fits)) return(.carrier_integration_fit_cache$fits)
  set.seed(1)
  k <- 48L
  data <- data.frame(g1 = factor(rep(c(5, 10, 20), length.out = k), levels = c(5, 10, 20)),
    g2 = factor(rep(1:4, length.out = k)[sample.int(k)], levels = 1:4), sei = stats::runif(k, .1, .3))
  data$yi <- stats::rnorm(k, c(0, .2, .4)[as.integer(data$g1)] + c(0, .1, .3, .5)[as.integer(data$g2)], data$sei)
  fit <- function(mods, contrast = "treatment", prior_effect = NULL, prior_mods = NULL, averaged = FALSE){
    args <- list(yi = quote(yi), sei = quote(sei), mods = mods, data = data, measure = "SMD",
      set_contrast_factor_predictors = contrast, chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 1, silent = TRUE)
    if(!is.null(prior_effect)) args$prior_effect <- prior_effect
    if(!is.null(prior_mods)) args$prior_mods <- prior_mods
    expect_warning(do.call(if(averaged) BMA.norm else brma, args), "Only one chain was run", fixed = TRUE)
  }
  fits <- list(treatment = fit(~g1-1), meandif = fit(~g1-1, "meandif"),
    nonzero = fit(~g1, prior_effect = BayesTools::prior("point", list(1))), averaged = fit(~g1, averaged = TRUE),
    ordinary = fit(~g1, "meandif"), gated = fit(~g1, prior_effect = BayesTools::prior("point", list(1)),
      prior_mods = list(g1 = BayesTools::prior_ordered(BayesTools::prior_spike_and_slab(BayesTools::prior("gamma", list(2, 2))), allocation = c(.25, .75)))))
  .carrier_integration_fit_cache$fits <- fits
  fits
}

test_that("expanded public factor groups retain every statement-major result row", {
  skip_on_cran()
  fit <- .carrier_integration_fits()$ordinary
  statements <- c("g1 > 0", "mu > 0", "g1 > 0")
  separate <- lapply(statements, function(statement) hypothesis(fit, statement, density_method = "KDE", columns = "all", seed = 1, standardized_coefficients = TRUE))
  combined <- hypothesis(fit, statements, density_method = "KDE", columns = "all", seed = 1, standardized_coefficients = TRUE)
  expect_identical(nrow(combined), 7L)
  expected <- unlist(lapply(separate, function(result) attr(result, "raw_log_BF")), use.names = FALSE)
  expect_identical(attr(combined, "raw_log_BF"), expected)
  expect_length(attr(combined, "numerical_diagnostics"), 7L)
  expect_length(attr(combined, "prior_numerical_diagnostics"), 7L)
  expect_length(attr(combined, "hypothesis_ast")$statements, 3L)
  for(logBF in c(FALSE, TRUE)) for(BF01 in c(FALSE, TRUE)){
    viewed <- update(combined, logBF = logBF, BF01 = BF01)
    expect_identical(attr(viewed, "raw_log_BF"), expected)
    expected_view <- unlist(lapply(separate, function(result) as.numeric(update(result, logBF = logBF, BF01 = BF01)$BF)), use.names = FALSE)
    expect_identical(as.numeric(viewed$BF), expected_view)
  }
  plain <- .output_table_as_long_data_frame(combined, "location")
  expect_identical(nrow(plain), 7L)
  expect_null(attr(plain$BF, "canonical_log_BF"))
})
