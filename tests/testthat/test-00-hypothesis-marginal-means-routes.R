.marginal_means_route_test_parameter <- function(a, b, condition_keys = NULL) {

  # Each level is one N(0, 1) coefficient of the joint prior context.
  levels <- list(
    A = with_draw_metadata(
      a,
      class          = c("marginal_posterior.simple", "numeric"),
      linear_weights = c(a = 1, b = 0),
      prior_density  = BayesTools::prior("normal", list(mean = 0, sd = 1)),
      atoms          = BayesTools::posterior_atom_attribute()
    ),
    B = with_draw_metadata(
      b,
      class          = c("marginal_posterior.simple", "numeric"),
      linear_weights = c(a = 0, b = 1),
      prior_density  = BayesTools::prior("normal", list(mean = 0, sd = 1)),
      atoms          = BayesTools::posterior_atom_attribute()
    )
  )
  if (!is.null(condition_keys)) {
    BayesTools::posterior_metadata(levels[["A"]], "condition") <-
      list(condition_key = condition_keys[[1L]])
    BayesTools::posterior_metadata(levels[["B"]], "condition") <-
      list(condition_key = condition_keys[[2L]])
  }
  class(levels) <- c(
    "marginal_posterior.factor",
    "marginal_posterior",
    "list"
  )
  attr(levels, "parameter") <- "mu_alloc"
  BayesTools::posterior_metadata(levels, "prior_context") <-
    BayesTools:::.prior_density_context(
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

  # Refused within one statement by its plan, and across the statements of
  # one request by the request's refusal; both with the statement class.
  expect_error(
    hypothesis(object, "alloc[A] = 0 vs alloc[A] > 0"),
    "cannot mix point and region",
    class = "RoBMA_hypothesis_statement"
  )
  expect_error(
    hypothesis(object, "alloc[A] = 0 vs alloc[B] != 0"),
    "spanning multiple marginal-means levels",
    class = "RoBMA_hypothesis_statement"
  )
  expect_error(
    hypothesis(object, c("alloc[A] = 0", "alloc[A] > 0")),
    "cannot mix point-null and region statements",
    class = "RoBMA_hypothesis_statement"
  )
  expect_error(
    hypothesis(object, c("alloc[A] = 0", "alloc[B] = 0")),
    "spanning multiple marginal-means levels",
    class = "RoBMA_hypothesis_statement"
  )
  # These are statements to restate, not unavailable tests: the statement
  # class only, as for a nonlinear expression of the levels.
  single <- .marginal_means_route_test_object(model_averaged = FALSE)
  requests <- list(
    list(object, "alloc[A] = 0 vs alloc[A] > 0"),
    list(object, "alloc[A] = 0 vs alloc[B] != 0"),
    list(object, c("alloc[A] = 0", "alloc[A] > 0")),
    list(object, c("alloc[A] = 0", "alloc[B] = 0")),
    list(single, "alloc[A] * alloc[B] = 0")
  )
  for (request in requests) {
    error <- tryCatch(hypothesis(request[[1L]], request[[2L]]), error = identity)
    expect_identical(
      class(error),
      c("RoBMA_hypothesis_statement", "error", "condition"),
      info = paste(request[[2L]], collapse = "; ")
    )
  }
  expect_error(
    hypothesis(single, "alloc[A] * alloc[B] = 0"),
    "A linear target must be a linear combination of levels of 'mu_alloc' and numbers.",
    fixed = TRUE
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
  plans <- function(statements) {
    lapply(statements, function(statement) {
      .hypothesis_plan_marginal_means(object, statement, parameter = "mu_alloc")
    })
  }

  mixed             <- plans("mu_alloc[A] > 0 vs mu_alloc[A] = 0")
  cross_level_point <- plans("mu_alloc[A] - mu_alloc[B] = 0")
  mixed_statements  <- plans(c("mu_alloc[A] = 0", "mu_alloc[A] > 0"))

  for (request in list(mixed, cross_level_point, mixed_statements)) {
    for (plan in request) {
      expect_identical(plan[["inference_type"]], "averaged")
      expect_null(plan[["refusal"]])
    }
    expect_null(.hypothesis_plan_marginal_means_request_refusal(
      plans = request, object = object, parameter = "mu_alloc"
    ))
  }
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
  BayesTools::posterior_metadata(levels[["A"]], "posterior_ordinate") <-
    BayesTools::posterior_ordinate_attribute(
      value = 0, ordinate = 1e-12, method = "bogus-child-A",
      density_method = "qCMDE"
    )
  BayesTools::posterior_metadata(levels[["B"]], "posterior_ordinate") <-
    BayesTools::posterior_ordinate_attribute(
      value = 0, ordinate = 1e12, method = "bogus-child-B",
      density_method = "qCMDE"
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


test_that("the plans of a marginal-means object are built once and shared by its methods", {

  .hypothesis_plan_cache_clear()
  withr::defer(.hypothesis_plan_cache_clear())
  object   <- .marginal_means_route_test_object(model_averaged = FALSE)
  original <- .hypothesis_plan_marginal_means_build
  builds   <- 0L
  testthat::local_mocked_bindings(
    .hypothesis_plan_marginal_means_build = function(...) {
      builds <<- builds + 1L
      original(...)
    },
    .package = "RoBMA"
  )

  first  <- hypothesis(object, "alloc[A] > alloc[B]", columns = "all", seed = 913)
  expect_identical(builds, 1L)
  again  <- hypothesis(object, "alloc[A] > alloc[B]", columns = "all", seed = 913)
  expect_identical(builds, 1L)
  expect_identical(again, first)
  # Another statement is another plan; the parameter the statement selects, or
  # names, is what the plan is of.
  hypothesis(object, "alloc[A] > 0", columns = "all", seed = 913)
  expect_identical(builds, 2L)
  hypothesis(object, "alloc[A] > 0", parameter = "mu_alloc", columns = "all",
             seed = 913)
  expect_identical(builds, 2L)
  hypothesis(object, "alloc[B] > 0", columns = "all", seed = 913)
  expect_identical(builds, 3L)
  # hypothesis_quantities() plans the statements of every level, and plans them
  # once for the next call.
  quantities <- hypothesis_quantities(object)
  after      <- builds
  expect_gt(after, 3L)
  expect_identical(hypothesis_quantities(object), quantities)
  expect_identical(builds, after)
  # Another object, even one that differs in a label only, has plans of its own.
  perturbed <- object
  perturbed[["term_map"]][["label"]] <- "another label"
  hypothesis(perturbed, "alloc[A] > alloc[B]", columns = "all", seed = 913)
  expect_identical(builds, after + 1L)

  # The plans of the cache are those of a call without one.
  testthat::local_mocked_bindings(
    .hypothesis_plan_cache = function(object = NULL) new.env(parent = emptyenv()),
    .package = "RoBMA"
  )
  expect_identical(
    hypothesis(object, "alloc[A] > alloc[B]", columns = "all", seed = 913),
    first
  )
  expect_identical(hypothesis_quantities(object), quantities)
})
