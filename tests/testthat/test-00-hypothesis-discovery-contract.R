test_that("hypothesis plans refuse methods by the runtime qCMDE/IWMDE capability", {

  # A point plan whose targets are eligible: every method is available up to
  # the capability of the fitted model.
  plan <- list(point = TRUE, refusal = NULL, targets = list(), method_refusals = list())
  glmm <- structure(list(), class = c("brma.glmm", "brma"))
  plan <- .hypothesis_plan_finish(plan, glmm)
  capability <- .iwmde_capability(object = glmm, density_method = "IWMDE")

  expect_null(.hypothesis_plan_status(plan, "KDE"))
  expect_null(.hypothesis_plan_status(plan, "qCMDE"))
  expect_null(.hypothesis_plan_status(plan, "normal"))
  expect_identical(
    .hypothesis_plan_status(plan, "IWMDE"),
    list(
      reason = capability[["reason"]],
      class  = c("RoBMA_hypothesis_method", "RoBMA_hypothesis_unavailable")
    )
  )
  expect_error(
    .hypothesis_plan_check(plan, "IWMDE"),
    capability[["reason"]],
    fixed = TRUE,
    class = "RoBMA_hypothesis_method"
  )
  rendered <- .hypothesis_quantities_render_plans(
    list(point = list(plan), region = list(plan), contrast = NULL),
    bracket = FALSE
  )
  expect_identical(rendered[["point_test_methods"]], "KDE, qCMDE")
  expect_identical(rendered[["reason"]], capability[["reason"]])
  expect_true(is.na(rendered[["contrast_test"]]))

  # A refusal of the whole plan (a fixed quantity) applies to point and
  # region statements and every method.
  fixed <- .hypothesis_plan_finish(
    list(point = FALSE, refusal = .hypothesis_plan_fixed_refusal("mu"),
         targets = list(), method_refusals = list()),
    glmm
  )
  rendered <- .hypothesis_quantities_render_plans(
    list(point = list(fixed), region = list(fixed), contrast = NULL),
    bracket = FALSE
  )
  expect_false(rendered[["point_test"]])
  expect_false(rendered[["direction_test"]])
  expect_identical(rendered[["point_test_methods"]], "")
  expect_match(rendered[["reason"]], "fixed by the fitted model", fixed = TRUE)
})


test_that("hypothesis discovery shares the runtime qCMDE/IWMDE capability", {

  object <- single_sd_random_object(BayesTools::prior("gamma", list(2, 2)))
  capability <- .iwmde_capability(
    object         = object,
    density_method = "qCMDE"
  )
  out <- hypothesis_quantities(object)
  intercept <- out[out[["parameter"]] == "mu_intercept", , drop = FALSE]

  expect_false(capability[["available"]])
  expect_identical(unique(intercept[["point_test_methods"]]), "KDE")
  expect_identical(unique(intercept[["reason"]]), capability[["reason"]])
  expect_error(
    hypothesis(object, "mu_intercept = 0", density_method = "qCMDE"),
    capability[["reason"]],
    fixed = TRUE,
    class = "RoBMA_hypothesis_method"
  )
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

  # Objects without a fit cannot be tested.
  unfitted <- hypothesis_quantities(structure(list(), class = "brma"))
  expect_false(any(unfitted[["point_test"]]))
  expect_identical(unique(unfitted[["reason"]]), "'hypothesis' requires a fitted brma object.")
})


test_that("random discovery uses the authoritative likelihood-aware target", {

  object <- single_sd_random_object(BayesTools::prior("gamma", list(2, 2)))
  attr(object[["data"]], "known_V") <- TRUE
  parameter <- "(mu) tau(intercept)"
  plan_of <- function() .hypothesis_plans(object, paste0("`", parameter, "` = 0.3"))[[1L]]

  testthat::local_mocked_bindings(
    .brma_random_parameter_density_target = function(object, parameter, ...) {
      list(reason = "unsupported")
    },
    .package = "RoBMA"
  )
  unsupported <- plan_of()
  expect_null(.hypothesis_plan_status(unsupported, "KDE"))
  for (method in c("qCMDE", "IWMDE")) {
    expect_identical(
      .hypothesis_plan_status(unsupported, method),
      list(reason = "unsupported",
           class  = c("RoBMA_hypothesis_target", "RoBMA_hypothesis_unavailable"))
    )
  }

  testthat::local_mocked_bindings(
    .brma_random_parameter_density_target = function(object, parameter, ...) {
      list(parameter = "rho_raw", parameter_spec = list(type = "primitive"))
    },
    .package = "RoBMA"
  )
  supported <- plan_of()
  for (method in c("KDE", "qCMDE", "IWMDE")) {
    expect_null(.hypothesis_plan_status(supported, method), info = method)
  }
})


# A marginal mean of a fixture: draws with an exact normal prior density and
# declared atoms (none, or a point mass carrying all mass for a fixed mean).
.discovery_marginal_mean <- function(values, fixed = FALSE) {

  with_draw_metadata(
    values,
    class         = c("marginal_posterior.simple", "marginal_posterior"),
    prior_density = BayesTools::prior("normal", list(0, 1)),
    atoms         = if (fixed) {
      BayesTools::posterior_atom_attribute(point_masses = data.frame(x = 0, mass = 1))
    } else {
      BayesTools::posterior_atom_attribute()
    }
  )
}


.discovery_marginal_object <- function(term_map, means, source_object) {

  structure(list(
    term_map      = term_map,
    inference     = list(averaged = means, conditional = means),
    source_object = source_object
  ), class = "marginal_means.brma")
}


test_that("marginal discovery shares its source-model IWMDE capability", {

  source_object <- structure(
    list(data = structure(list(), random = TRUE)),
    class = c("brma.mv", "brma")
  )
  object <- .discovery_marginal_object(
    term_map      = data.frame(term = "intercept", parameter = "mu_intercept"),
    means         = list(mu_intercept = .discovery_marginal_mean(c(-1, 0, 1))),
    source_object = source_object
  )

  out <- hypothesis_quantities(object)

  # A source object without a fit cannot compute qCMDE/IWMDE ordinates.
  expect_identical(out[["point_test_methods"]], "KDE")
  expect_match(out[["reason"]], "does not contain the source fitted brma", fixed = TRUE)

  source_object[["fit"]] <- list(fitted = TRUE)
  object[["source_object"]] <- source_object
  out <- hypothesis_quantities(object)
  capability <- .iwmde_capability(object = source_object, density_method = "IWMDE")
  expect_identical(out[["point_test_methods"]], "KDE")
  expect_identical(out[["reason"]], capability[["reason"]])
})


test_that("marginal-means discovery reports selectable levels", {

  object <- .discovery_marginal_object(
    term_map = data.frame(
      term      = c("intercept", "group"),
      parameter = c("mu_intercept", "mu_group")
    ),
    means = list(
      mu_intercept = .discovery_marginal_mean(c(-1, 0, 1)),
      mu_group     = list(
        a = .discovery_marginal_mean(c(-1, 0)),
        b = .discovery_marginal_mean(c(0, 1))
      )
    ),
    source_object = structure(list(fit = list(fitted = TRUE)), class = c("brma.glmm", "brma"))
  )

  out <- hypothesis_quantities(object)
  scalar <- out[out[["parameter"]] == "mu_intercept", , drop = FALSE]
  grouped <- out[out[["parameter"]] == "mu_group", , drop = FALSE]

  expect_true(all(is.na(scalar[["bracket"]])))
  expect_identical(grouped[["bracket"]], c("mu_group[a]", "mu_group[b]"))
  expect_identical(unique(out[["point_test_methods"]]), "KDE, qCMDE")
  expect_true(is.na(scalar[["contrast_test"]]))
  expect_identical(grouped[["contrast_test"]], c(FALSE, FALSE))
})


test_that("marginal-means discovery rejects a fixed no-intercept scalar", {

  object <- .discovery_marginal_object(
    term_map = data.frame(term = "intercept", parameter = "mu_intercept"),
    means    = list(mu_intercept = list(
      intercept = .discovery_marginal_mean(rep(0, 20L), fixed = TRUE)
    )),
    source_object = structure(list(fit = list(fitted = TRUE)), class = "brma")
  )

  out <- hypothesis_quantities(object)

  expect_identical(out[["bracket"]], "mu_intercept[intercept]")
  expect_false(out[["point_test"]])
  expect_false(out[["direction_test"]])
  expect_identical(out[["point_test_methods"]], "")
  expect_match(out[["reason"]], "fixed by the fitted model")
  expect_error(
    hypothesis(object, "mu_intercept[intercept] > 0"),
    class = "RoBMA_hypothesis_fixed"
  )
})


test_that("marginal-means discovery keeps sampled siblings of a fixed level", {

  # The fixed level is declared by its atoms, not by constant draws: the
  # sampled levels' constant draws do not fix them.
  object <- .discovery_marginal_object(
    term_map = data.frame(term = "x", parameter = "mu_x"),
    means    = list(mu_x = list(
      `-1SD` = .discovery_marginal_mean(c(-1, -2, -3)),
      `0SD`  = .discovery_marginal_mean(c(0, 0, 0), fixed = TRUE),
      `1SD`  = .discovery_marginal_mean(c(1, 1, 1))
    )),
    source_object = structure(list(fit = list(fitted = TRUE)), class = "brma")
  )

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

  ordinate <- function(value) BayesTools::posterior_ordinate_attribute(
    value          = value,
    ordinate       = value + 1,
    method         = "qCMDE",
    density_method = "qCMDE",
    diagnostics    = list(marker = value)
  )
  posterior <- stats::rnorm(10)
  BayesTools::posterior_metadata(posterior, "posterior_ordinate") <-
    BayesTools::posterior_ordinate_append(ordinate(0), ordinate(1))
  refs <- data.frame(level = NA_character_, value = 1)

  out <- .hypothesis_brma_keep_requested_ordinates(posterior, refs)
  entries <- .iwmde_posterior_ordinate_entries(
    BayesTools::posterior_metadata(out, "posterior_ordinate")
  )

  expect_length(entries, 1L)
  expect_equal(entries[[1L]][["value"]], 1)
})
