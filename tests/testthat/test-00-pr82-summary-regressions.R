test_that("unpaired random inclusion gates retain usable labels", {

  quantities <- data.frame(
    canonical_name = c("inclusion(a)", "inclusion(b)", "inclusion(c)"),
    role = rep("random_inclusion", 3L)
  )
  quantities$extraction_key <- I(list(
    list(source_parameter = "missing"), list(), list(source_parameter = "paired")
  ))
  local_mocked_bindings(
    parameter_catalog = function(...) list(quantities = quantities),
    .package = "BayesTools"
  )
  local_mocked_bindings(
    .random_inclusion_sd_names = function(...) c(paired = "tau"),
    .package = "RoBMA"
  )
  expect_identical(.summary_random_inclusion_labels(
    list(), quantities$canonical_name, quantities$canonical_name
  ), c("a", "b", "tau"))
})

test_that("summary sections preserve colliding labels and values", {

  tables <- list(data.frame(Mean = 1, row.names = "tau"),
                 data.frame(Mean = 2, row.names = "tau"))
  result <- .summary_brma_combine_tables(tables, "Estimates")
  expect_identical(result$Mean, c(1, 2))
  expect_identical(rownames(result), c("tau", "tau.1"))
})

test_that("integer and double zero spikes have the same model title", {

  for (zero in list(0L, 0)) {
    object <- brma(yi = c(0.1, 0.2), sei = c(0.2, 0.3), measure = "SMD",
      prior_heterogeneity = prior("point", list(location = zero)),
      only_priors = TRUE, silent = TRUE)
    expect_match(.summary.brma_model_names(object), "Fixed-Effect", fixed = TRUE)
  }
})
