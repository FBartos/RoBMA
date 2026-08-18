.marginal_means_route_test_parameter <- function(a, b, condition_keys = NULL) {

  levels <- list(
    A = structure(
      a,
      class           = c("marginal_posterior.simple", "numeric"),
      linear_weights  = c(a = 1, b = 0),
      posterior_atoms = BayesTools::posterior_atom_attribute()
    ),
    B = structure(
      b,
      class           = c("marginal_posterior.simple", "numeric"),
      linear_weights  = c(a = 0, b = 1),
      posterior_atoms = BayesTools::posterior_atom_attribute()
    )
  )
  if (!is.null(condition_keys)) {
    attr(levels[["A"]], "condition_key") <- condition_keys[[1L]]
    attr(levels[["B"]], "condition_key") <- condition_keys[[2L]]
  }
  class(levels) <- c(
    "marginal_posterior.factor",
    "marginal_posterior",
    "list"
  )
  attr(levels, "parameter") <- "mu_alloc"
  attr(levels, "prior_density_context") <- BayesTools:::.prior_density_context(
    prior_list = list(
      a = BayesTools::prior("normal", list(mean = 0, sd = 1)),
      b = BayesTools::prior("normal", list(mean = 0, sd = 1))
    ),
    column_names = c("a", "b"),
    n_grid       = 128
  )

  return(levels)
}


.marginal_means_route_test_object <- function(model_averaged = TRUE) {

  averaged <- .marginal_means_route_test_parameter(
    a = c(rep(1, 75), rep(-1, 25)),
    b = rep(0, 100)
  )
  conditional <- .marginal_means_route_test_parameter(
    a              = rep(1, 100),
    b              = rep(0, 100),
    condition_keys = c("event-A", "event-B")
  )
  inference <- structure(
    list(
      averaged    = list(mu_alloc = averaged),
      conditional = list(mu_alloc = conditional),
      inference   = list()
    ),
    class = c("marginal_inference", "list")
  )
  source_class <- if (model_averaged) c("RoBMA", "brma") else "brma"
  object <- list(
    inference      = inference,
    term_map       = data.frame(
      term             = "alloc",
      parameter        = "mu_alloc",
      label            = "alloc",
      stringsAsFactors = FALSE
    ),
    density_method = "KDE",
    model_averaged = model_averaged,
    source_object  = structure(list(), class = source_class)
  )
  class(object) <- "marginal_means.brma"

  return(object)
}


test_that("cross-level region hypotheses use direct averaged odds", {

  object <- .marginal_means_route_test_object(model_averaged = TRUE)
  out <- hypothesis(
    object,
    "alloc[A] > alloc[B]",
    columns = "all",
    seed    = 913
  )

  posterior_odds <- (75 / 100) / (25 / 100)
  prior_odds     <- 0.5 / 0.5
  direct_bf      <- posterior_odds / prior_odds

  expect_equal(out[["posterior"]], posterior_odds, tolerance = 1e-12)
  expect_equal(out[["prior"]], prior_odds, tolerance = 0.04)
  expect_equal(as.numeric(attr(out, "raw_BF")), direct_bf, tolerance = 0.12)
  expect_equal(out[["method"]], "prior-posterior odds")
})


test_that("cross-level region counts ignore density sample budgets", {

  object <- .marginal_means_route_test_object(model_averaged = TRUE)
  hypothesis_text <- "alloc[A] > alloc[B]"

  kde <- hypothesis(
    object,
    hypothesis_text,
    columns = "all",
    seed    = 913
  )
  iwmde <- hypothesis(
    object,
    hypothesis_text,
    density_method  = "IWMDE",
    density_control = list(samples = 20),
    columns         = "all",
    seed            = 913
  )
  qcmde <- hypothesis(
    object,
    hypothesis_text,
    density_method  = "qCMDE",
    density_control = list(samples = 20),
    columns         = "all",
    seed            = 913
  )

  expect_identical(iwmde, kde)
  expect_identical(qcmde, kde)
  expect_equal(kde[["posterior"]], 3, tolerance = 1e-12)
})


test_that("marginal-means point routes fail closed when events are incoherent", {

  object <- .marginal_means_route_test_object(model_averaged = TRUE)

  expect_error(
    hypothesis(object, "alloc[A] = 0 vs alloc[A] > 0"),
    "cannot mix point and region"
  )
  expect_error(
    hypothesis(object, "alloc[A] = 0 vs alloc[B] != 0"),
    "spanning multiple marginal-means levels"
  )
  expect_error(
    hypothesis(object, c("alloc[A] = 0", "alloc[A] > 0")),
    "cannot mix point-null and region statements"
  )
})


test_that("single-model point hypotheses use averaged marginal draws", {

  object <- .marginal_means_route_test_object(model_averaged = FALSE)
  captured <- NULL
  testthat::local_mocked_bindings(
    hypothesis_BF = function(posterior, ...) {

      captured <<- posterior[["conditional"]]
      return("ok")
    },
    .package = "BayesTools"
  )

  expect_equal(hypothesis(object, "alloc[A] = 0"), "ok")
  expect_identical(captured, object[["inference"]][["averaged"]])
})


test_that("single-model marginal hypotheses share one averaged posterior", {

  object <- .marginal_means_route_test_object(model_averaged = FALSE)

  mixed <- .hypothesis_marginal_means_route(
    object     = object,
    hypothesis = BayesTools::hypothesis_parse(
      "mu_alloc[A] > 0 vs mu_alloc[A] = 0"
    ),
    parameter  = "mu_alloc"
  )
  cross_level_point <- .hypothesis_marginal_means_route(
    object     = object,
    hypothesis = BayesTools::hypothesis_parse(
      "mu_alloc[A] - mu_alloc[B] = 0"
    ),
    parameter  = "mu_alloc"
  )
  mixed_statements <- .hypothesis_marginal_means_route(
    object     = object,
    hypothesis = BayesTools::hypothesis_parse(c(
      "mu_alloc[A] = 0",
      "mu_alloc[A] > 0"
    )),
    parameter  = "mu_alloc"
  )

  expect_identical(mixed, list(route = "mixed", inference_type = "averaged"))
  expect_identical(
    cross_level_point,
    list(route = "point", inference_type = "averaged")
  )
  expect_identical(
    mixed_statements,
    list(route = "mixed", inference_type = "averaged")
  )
})


test_that("cross-level point densities ignore child-level ordinates", {

  object <- .marginal_means_route_test_object(model_averaged = FALSE)
  hypothesis_text <- "alloc[A] - alloc[B] = 0"
  baseline <- suppressWarnings(hypothesis(
    object,
    hypothesis_text,
    columns = "all",
    seed    = 819
  ))

  levels <- object[["inference"]][["averaged"]][["mu_alloc"]]
  attr(levels[["A"]], "posterior_ordinate") <- list(
    value = 0, ordinate = 1e-12, method = "bogus-child-A"
  )
  attr(levels[["B"]], "posterior_ordinate") <- list(
    value = 0, ordinate = 1e12, method = "bogus-child-B"
  )
  object[["inference"]][["averaged"]][["mu_alloc"]] <- levels
  perturbed <- suppressWarnings(hypothesis(
    object,
    hypothesis_text,
    columns = "all",
    seed    = 819
  ))

  expect_identical(perturbed, baseline)
  expect_identical(baseline[["method"]], "Savage-Dickey")
})


test_that("hypothesis labels retain expanded factor levels", {

  display <- BayesTools::hypothesis_parse("Preregistered > 0")
  transformed <- BayesTools::hypothesis_parse("mu_Preregistered > 0")
  out <- data.frame(
    Alternative = rep("mu_Preregistered > 0", 2L),
    Null        = rep("mu_Preregistered <= 0", 2L),
    BF          = c(2, 3),
    row.names   = c(
      "mu_Preregistered[Not Pre-Registered]",
      "mu_Preregistered[Pre-Registered]"
    ),
    check.names = FALSE
  )
  attr(out, "hypothesis_ast") <- transformed

  restored <- .hypothesis_brma_restore_hypothesis_labels(
    out             = out,
    hypothesis      = display,
    parameter_label = "Preregistered"
  )

  expect_identical(
    restored[["Alternative"]],
    rep("Preregistered > 0", 2L)
  )
  expect_identical(restored[["Null"]], rep("Preregistered <= 0", 2L))
  expect_identical(
    rownames(restored),
    c(
      "Preregistered[Not Pre-Registered]",
      "Preregistered[Pre-Registered]"
    )
  )
  expect_identical(attr(restored, "hypothesis_ast", exact = TRUE), display)
})
