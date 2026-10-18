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


test_that("selected ordered prior styles resolve scalar defaults and alias conflicts", {

  expect_identical(.set_dots_prior(NULL, ordered = TRUE, posterior_col = "purple"),
    list(col = "purple", lty = 2, col.fill = "#B3B3B34C"))
  aliases <- list(color = "blue", colour = "red", lty = 1, linetype = "dotted")
  expect_identical(.set_dots_prior(aliases, ordered = TRUE, posterior_col = "purple"),
    list(color = "blue", colour = "red", lty = "dotted", linetype = "dotted",
      col = "blue", col.fill = "#B3B3B34C"))
  explicit <- list(col = "green", color = "blue", colour = "red", lty = 1,
    linetype = "dotdash")
  expected <- explicit
  expected[["lty"]] <- "dotdash"
  expected[["col.fill"]] <- "#B3B3B34C"
  expect_identical(.set_dots_prior(explicit, ordered = TRUE, posterior_col = "purple"),
    expected)
  expect_identical(.set_dots_prior(list(colour = "red"), ordered = TRUE,
    posterior_col = "purple"), list(colour = "red", col = "red", lty = 2,
    col.fill = "#B3B3B34C"))
})


test_that("ordered style selection uses declared coefficient dependencies", {

  quantity <- function(dependencies, term = "f") {
    result <- data.frame(quantity_id = "target", term = term, formula_parameter = "mu")
    result[["extraction_key"]] <- I(list(list(type = "factor_level", dependencies = dependencies)))
    result
  }
  owners <- "mu_f"
  calls <- 0L
  testthat::local_mocked_bindings(JAGS_ordered_parameter_spec = function(fit, parameters) {
    calls <<- calls + 1L
    expect_identical(parameters, owners)
    list(f = list(coefficient_names = c("mu_f[1]", "mu_f[2]")),
      h = list(coefficient_names = c("mu_h[1]", "mu_h[2]")))
  }, .package = "BayesTools")
  fit <- list()
  ordered <- BayesTools::prior_ordered(BayesTools::prior("normal", list(0, 1)))
  attr(fit, "prior_list") <- list(mu_f = ordered, mu_h = ordered,
    mu_x = BayesTools::prior("normal", list(0, 1)))
  object <- list(fit = fit)
  for (dependencies in list("mu_f[1]", c("mu_x", "mu_f[2]"))) {
    entry <- list(selection = list(quantities = quantity(dependencies)))
    expect_true(.plot_brma_ordered_target(object, entry))
  }
  entry <- list(selection = list(quantities = quantity("mu_x", term = "x")))
  expect_false(.plot_brma_ordered_target(object, entry))
  expect_identical(calls, 2L)
  malformed <- quantity("mu_f[1]")
  malformed[["formula_parameter"]] <- NULL
  expect_false(.plot_brma_ordered_target(object, list(selection = list(quantities = malformed))))
  expect_identical(calls, 2L)
  entry <- list(selection = list(quantities = quantity(character(), term = "")),
    member_quantity_ids = c("parent_f", "parent_h"))
  parent <- rbind(quantity(c("mu_f[1]", "mu_f[2]")), quantity("mu_h[1]", term = "h"))
  parent[["quantity_id"]] <- c("parent_f", "parent_h")
  local_mocked_bindings(.brma_parameter_catalog_metadata = function(...) {
    list(catalog = list(quantities = parent))
  }, .package = "RoBMA")
  owners <- c("mu_f", "mu_h")
  expect_true(.plot_brma_ordered_target(object, entry))
  expect_identical(calls, 3L)
})


test_that("ordered plot ownership validates only selected bound terms", {

  normal <- BayesTools::prior("normal", list(0, 1))
  ordered <- BayesTools::prior_ordered(normal)
  data <- data.frame(f = base::ordered(c("base", "middle", "high", "base")),
    h = base::ordered(c("base", "middle", "high", "base")), x = seq_len(4))
  compiled <- BayesTools::JAGS_formula(~ f + h + x, "mu", data,
    list(intercept = normal, f = ordered, h = ordered, x = normal))
  columns <- c("mu_intercept", "mu_f[1]", "mu_f[2]", "mu_h[1]", "mu_h[2]", "mu_x")
  values <- matrix(1, 4L, length(columns), dimnames = list(NULL, columns))
  fit <- structure(list(mcmc = coda::mcmc.list(coda::mcmc(values)), sample = 4L),
    class = c("runjags", "BayesTools_fit", "list"))
  attr(fit, "prior_list") <- compiled[["prior_list"]]
  attr(fit, "formula_design") <- list(mu = compiled[["formula_design"]])
  attr(fit, "parameter_map") <- BayesTools:::.bt_build_parameter_map(
    columns = columns, prior_list = compiled[["prior_list"]],
    formula_design = list(mu = compiled[["formula_design"]]))
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  fit <- BayesTools:::.bt_attach_fit_contract(fit)
  quantities <- BayesTools::parameter_catalog(fit)[["quantities"]]
  priors <- attr(fit, "prior_list", exact = TRUE)
  metadata <- attr(priors[["mu_h"]], "ordered_metadata", exact = TRUE)
  metadata[["numeric_literals"]] <- NULL
  attr(priors[["mu_h"]], "ordered_metadata") <- metadata
  attr(fit, "prior_list") <- priors
  object <- list(fit = fit)

  entry <- function(term) {
    list(selection = list(quantities = quantities[quantities[["term"]] == term, , drop = FALSE]))
  }
  expect_false(.plot_brma_ordered_target(object, entry("x")))
  expect_true(.plot_brma_ordered_target(object, entry("f")))
  error <- tryCatch(.plot_brma_ordered_target(object, entry("h")), error = identity)
  expect_s3_class(error, "BayesTools_ordered_metadata_unavailable")
  expect_s3_class(error, "BayesTools_refit_required")
})


test_that("plot dispatch forwards matched ordered styles and ordinary prior defaults", {

  quantities <- data.frame(quantity_id = "target", term = "f", formula_parameter = "mu")
  quantities[["extraction_key"]] <- I(list(list(dependencies = "mu_f[1]")))
  entry <- list(parameter = "selected", component = "mods",
    selection = list(quantities = quantities))
  fit <- list()
  attr(fit, "prior_list") <- list(mu_f = BayesTools::prior_ordered(
    BayesTools::prior("normal", list(0, 1))), mu_x = BayesTools::prior("normal", list(0, 1)))
  object <- structure(list(fit = fit), class = "brma")
  captured <- NULL
  local_mocked_bindings(
    .check_and_select_plot_parameter = function(...) entry,
    .brma_parameter_select_entry = function(...) entry,
    .plot_output_setup = function(...) NULL,
    .as_mixed_posteriors_parameters = function(...) "selected",
    .brma_as_mixed_posteriors = function(...) list(selected = seq_len(20)),
    .plot_brma_has_posterior_density = function(...) FALSE,
    .package = "RoBMA"
  )
  local_mocked_bindings(
    JAGS_ordered_parameter_spec = function(fit, parameters) {
      expect_identical(parameters, "mu_f")
      list(ordered = list(coefficient_names = "mu_f[1]"))
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
  entry[["selection"]][["quantities"]][["term"]] <- "x"
  entry[["selection"]][["quantities"]][["extraction_key"]] <- I(list(list(dependencies = "mu_x")))
  plot(object, prior = TRUE, par_name = "selected")
  expect_identical(captured[["dots_prior"]],
    list(col = "grey60", lty = 1, col.fill = "#B3B3B34C"))
})
