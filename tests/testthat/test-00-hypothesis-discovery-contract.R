test_that("hypothesis discovery reports fitted transform routes", {

  direct_entry <- list(
    role              = "fixed_coefficient",
    status            = "estimated",
    fixed_value       = NA_real_,
    formula_parameter = NA_character_
  )
  direct <- .hypothesis_quantities_brma_route(
    object = structure(list(), class = c("brma.glmm", "brma")),
    row    = data.frame(parameter = "mu", component = "mods"),
    entry  = direct_entry
  )
  expect_identical(direct[["point_test_methods"]], "KDE, qCMDE")
  expect_false(direct[["bracket"]])

  grouped_entry <- direct_entry
  grouped_entry[["role"]] <- "formula_coefficient_group"
  grouped <- .hypothesis_quantities_brma_route(
    object = structure(list(), class = "brma"),
    row    = data.frame(parameter = "mu_group", component = "mods"),
    entry  = grouped_entry
  )
  expect_true(grouped[["bracket"]])

  transformed_entry <- direct_entry
  transformed_entry[["formula_parameter"]] <- "log_tau"
  testthat::local_mocked_bindings(
    .hypothesis_brma_formula_coefficient_target = function(...) list(),
    .hypothesis_brma_formula_transform_route = function(target) {
      list(type = "exp_affine")
    },
    .package = "RoBMA"
  )
  transformed <- .hypothesis_quantities_brma_route(
    object = structure(list(), class = "brma"),
    row    = data.frame(parameter = "log_tau_x", component = "scale"),
    entry  = transformed_entry
  )
  expect_identical(transformed[["point_test_methods"]], "KDE")
  expect_identical(transformed[["direction_test_methods"]], "KDE")

  fixed_entry <- direct_entry
  fixed_entry[["status"]]      <- "fixed"
  fixed_entry[["fixed_value"]] <- 0
  fixed <- .hypothesis_quantities_brma_route(
    object = structure(list(), class = "brma"),
    row    = data.frame(parameter = "mu", component = "mods"),
    entry  = fixed_entry
  )
  expect_false(fixed[["point_test"]])
  expect_false(fixed[["direction_test"]])
})


test_that("marginal-means discovery reserves brackets for grouped levels", {

  object <- structure(list(
    term_map = data.frame(
      term      = c("intercept", "group"),
      parameter = c("mu_intercept", "mu_group")
    ),
    inference = list(conditional = list(
      mu_intercept = c(-1, 0, 1),
      mu_group     = list(a = c(-1, 0), b = c(0, 1))
    )),
    source_object = structure(list(), class = c("brma.glmm", "brma"))
  ), class = "marginal_means.brma")

  out <- hypothesis_quantities(object)
  scalar <- out[out[["parameter"]] == "mu_intercept", , drop = FALSE]
  grouped <- out[out[["parameter"]] == "mu_group", , drop = FALSE]

  expect_true(all(is.na(scalar[["bracket"]])))
  expect_true(all(grouped[["bracket"]] == "mu_group[level]"))
  expect_identical(unique(out[["point_test_methods"]]), "KDE, qCMDE")
})


test_that("hypothesis ordinate metadata is limited to requested values", {

  ordinate <- function(value) list(
    value      = value,
    ordinate   = value + 1,
    diagnostics = list(marker = value)
  )
  posterior <- stats::rnorm(10)
  attr(posterior, "posterior_ordinate") <- structure(
    list(status = "ok", ordinates = list(ordinate(0), ordinate(1))),
    class = c("BayesTools_posterior_ordinates", "list")
  )
  refs <- data.frame(level = NA_character_, value = 1)

  out <- .hypothesis_brma_keep_requested_ordinates(posterior, refs)
  entries <- .iwmde_posterior_ordinate_entries(
    attr(out, "posterior_ordinate", exact = TRUE)
  )

  expect_length(entries, 1L)
  expect_equal(entries[[1L]][["value"]], 1)
})
