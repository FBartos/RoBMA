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
  expect_identical(
    direct[["reason"]],
    .iwmde_capability(
      object         = structure(list(), class = c("brma.glmm", "brma")),
      density_method = "IWMDE"
    )[["reason"]]
  )

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
  expect_match(transformed[["reason"]], "atom-free, unconditional")

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


test_that("hypothesis discovery shares the runtime qCMDE/IWMDE capability", {

  data <- structure(list(), random = TRUE)
  object <- structure(
    list(data = data),
    class = c("brma.mv", "brma")
  )
  testthat::local_mocked_bindings(
    .brma_parameter_catalog = function(object) {

      data.frame(
        alias      = "mu",
        parameter = "mu",
        component = "mods",
        term       = "intercept",
        stringsAsFactors = FALSE
      )
    },
    .package = "RoBMA"
  )

  capability <- .iwmde_capability(
    object         = object,
    density_method = "qCMDE"
  )
  out <- hypothesis_quantities(object)

  expect_false(capability[["available"]])
  expect_identical(out[["point_test_methods"]], "KDE")
  expect_identical(out[["reason"]], capability[["reason"]])
  expect_error(
    .check_iwmde_available(object, "qCMDE/IWMDE hypothesis()"),
    capability[["reason"]],
    fixed = TRUE
  )

  known_v_object <- object
  attr(known_v_object[["data"]], "known_V") <- TRUE
  expect_true(.iwmde_capability(
    object         = known_v_object,
    density_method = "qCMDE"
  )[["available"]])
  expect_invisible(
    .check_iwmde_available(known_v_object, "qCMDE/IWMDE hypothesis()")
  )
})


test_that("marginal discovery shares its source-model IWMDE capability", {

  data <- structure(list(), random = TRUE)
  source_object <- structure(
    list(data = data),
    class = c("brma.mv", "brma")
  )
  object <- structure(list(
    term_map = data.frame(
      term      = "intercept",
      parameter = "mu_intercept"
    ),
    inference = list(conditional = list(
      mu_intercept = c(-1, 0, 1)
    )),
    source_object = source_object
  ), class = "marginal_means.brma")

  capability <- .iwmde_capability(
    object         = source_object,
    density_method = "IWMDE"
  )
  out <- hypothesis_quantities(object)

  expect_identical(out[["point_test_methods"]], "KDE")
  expect_identical(out[["reason"]], capability[["reason"]])
})


test_that("marginal-means discovery reports selectable levels", {

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
  expect_identical(grouped[["bracket"]], c("mu_group[a]", "mu_group[b]"))
  expect_identical(unique(out[["point_test_methods"]]), "KDE, qCMDE")
})


test_that("marginal-means discovery rejects a fixed no-intercept scalar", {

  object <- structure(list(
    term_map = data.frame(
      term      = "intercept",
      parameter = "mu_intercept"
    ),
    inference = list(conditional = list(
      mu_intercept = list(intercept = matrix(0, nrow = 1L, ncol = 20L))
    )),
    source_object = structure(list(), class = "brma")
  ), class = "marginal_means.brma")

  out <- hypothesis_quantities(object)

  expect_identical(out[["bracket"]], "mu_intercept[intercept]")
  expect_false(out[["point_test"]])
  expect_false(out[["direction_test"]])
  expect_identical(out[["point_test_methods"]], "")
  expect_identical(out[["direction_test_methods"]], "")
  expect_match(out[["reason"]], "fixed by the fitted model")
})


test_that("marginal-means discovery keeps sampled siblings of a fixed level", {

  object <- structure(list(
    term_map = data.frame(
      term      = "x",
      parameter = "mu_x"
    ),
    inference = list(conditional = list(
      mu_x = list(
        `-1SD` = c(-1, -2, -3),
        `0SD`  = c(0, 0, 0),
        `1SD`  = c(1, 2, 3)
      )
    )),
    source_object = structure(list(), class = "brma")
  ), class = "marginal_means.brma")

  out <- hypothesis_quantities(object)
  fixed <- out[["bracket"]] == "mu_x[0SD]"

  expect_identical(
    out[["bracket"]],
    c("mu_x[-1SD]", "mu_x[0SD]", "mu_x[1SD]")
  )
  expect_identical(out[["point_test"]], c(TRUE, FALSE, TRUE))
  expect_identical(out[["direction_test"]], c(TRUE, FALSE, TRUE))
  expect_identical(
    out[["point_test_methods"]][!fixed],
    rep("KDE, qCMDE, IWMDE", 2L)
  )
  expect_identical(out[["point_test_methods"]][fixed], "")
  expect_match(out[["reason"]][fixed], "fixed by the fitted model")
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
