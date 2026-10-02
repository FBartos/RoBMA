context("Ordered prior style forwarding")


test_that("ordered priors leave matched styles to BayesTools and retain explicit styles", {

  expect_identical(.set_dots_prior(NULL), list(col = "grey60", lty = 1, col.fill = "#B3B3B34C"))
  expect_identical(.set_dots_prior(NULL, ordered = TRUE), list(col.fill = "#B3B3B34C"))
  custom <- list(col = c("blue", "red"), lty = c(3, 4), col.fill = "white", lwd = 2)
  expect_identical(.set_dots_prior(custom, ordered = TRUE), custom)
  expect_identical(.set_dots_prior(list(linetype = "dotted"), ordered = TRUE),
    list(linetype = "dotted", col.fill = "#B3B3B34C"))
  style <- .set_dots_prior_marginal_means(NULL, 1L, 2L)
  expect_identical(style[["col"]], c("black", "black"))
  expect_identical(style[["lty"]], c(2, 2))
})


test_that("ordered style selection uses declared coefficient dependencies", {

  quantity <- function(dependencies) {
    result <- data.frame(quantity_id = "target")
    result[["extraction_key"]] <- I(list(list(type = "factor_level", dependencies = dependencies)))
    result
  }
  testthat::local_mocked_bindings(JAGS_ordered_parameter_spec = function(...) {
    list(ordered = list(coefficient_names = c("ordered[1]", "ordered[2]")))
  }, .package = "BayesTools")
  object <- list(fit = list())
  for (dependencies in list("ordered[1]", c("ordinary", "ordered[2]"))) {
    entry <- list(selection = list(quantities = quantity(dependencies)))
    expect_true(.plot_brma_ordered_target(object, entry))
  }
  entry <- list(selection = list(quantities = quantity("ordinary")))
  expect_false(.plot_brma_ordered_target(object, entry))
  entry <- list(selection = list(quantities = quantity(character())), member_quantity_ids = "parent")
  parent <- quantity(c("ordered[1]", "ordered[2]"))
  parent[["quantity_id"]] <- "parent"
  local_mocked_bindings(.brma_parameter_catalog_metadata = function(...) {
    list(catalog = list(quantities = parent))
  }, .package = "RoBMA")
  expect_true(.plot_brma_ordered_target(object, entry))
})


test_that("plot dispatch forwards matched ordered styles and ordinary prior defaults", {

  quantities <- data.frame(quantity_id = "target")
  quantities[["extraction_key"]] <- I(list(list(dependencies = "ordered[1]")))
  entry <- list(component = "mods", selection = list(quantities = quantities))
  object <- structure(list(fit = list()), class = "brma")
  captured <- NULL
  local_mocked_bindings(
    .check_and_select_plot_parameter = function(...) "selected",
    .brma_parameter_select_entry = function(...) entry,
    .plot_output_setup = function(...) NULL,
    .as_mixed_posteriors_parameters = function(...) "selected",
    .brma_as_mixed_posteriors = function(...) list(selected = seq_len(20)),
    .plot_brma_has_posterior_density = function(...) FALSE,
    .package = "RoBMA"
  )
  local_mocked_bindings(
    JAGS_ordered_parameter_spec = function(...) {
      list(ordered = list(coefficient_names = "ordered[1]"))
    },
    plot_posterior = function(...) {
      captured <<- list(...)
      invisible(NULL)
    }, .package = "BayesTools"
  )
  plot(object, prior = TRUE, par_name = "selected")
  expect_identical(captured[["dots_prior"]], list(col.fill = "#B3B3B34C"))
  custom <- list(col = "blue", linetype = "dotdash", col.fill = "white")
  plot(object, prior = TRUE, dots_prior = custom, par_name = "selected")
  expect_identical(captured[["dots_prior"]], custom)
  entry[["selection"]][["quantities"]][["extraction_key"]] <- I(list(list(dependencies = "ordinary")))
  plot(object, prior = TRUE, par_name = "selected")
  expect_identical(captured[["dots_prior"]],
    list(col = "grey60", lty = 1, col.fill = "#B3B3B34C"))
})
