.marginal_means_route_test_parameter <- function(a, b, condition_keys = NULL,
                                                 parameter = "mu_alloc") {

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
  attr(levels, "parameter") <- parameter
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


# The object above with a second marginal-means parameter, 'other', that has
# the same levels: 40% of the averaged draws of its level A lie above zero
# (75% for 'alloc').
.marginal_means_route_two_parameter_object <- function(model_averaged = TRUE) {

  object <- .marginal_means_route_test_object(model_averaged = model_averaged)
  object[["inference"]][["averaged"]][["mu_other"]] <-
    .marginal_means_route_test_parameter(
      a         = c(rep(1, 40), rep(-1, 60)),
      b         = rep(0, 100),
      parameter = "mu_other"
    )
  object[["inference"]][["conditional"]][["mu_other"]] <-
    .marginal_means_route_test_parameter(
      a              = rep(1, 100),
      b              = rep(0, 100),
      condition_keys = c("event-A", "event-B"),
      parameter      = "mu_other"
    )
  object[["term_map"]] <- rbind(
    object[["term_map"]],
    data.frame(
      term             = "other",
      parameter        = "mu_other",
      label            = "other",
      stringsAsFactors = FALSE
    )
  )

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
    class = "BayesTools_hypothesis_statement"
  )
  expect_error(
    hypothesis(object, "alloc[A] = 0 vs alloc[B] != 0"),
    "spanning multiple marginal-means levels",
    class = "BayesTools_hypothesis_statement"
  )
  expect_error(
    hypothesis(object, c("alloc[A] = 0", "alloc[A] > 0")),
    "cannot mix point-null and region statements",
    class = "BayesTools_hypothesis_statement"
  )
  expect_error(
    hypothesis(object, c("alloc[A] = 0", "alloc[B] = 0")),
    "spanning multiple marginal-means levels",
    class = "BayesTools_hypothesis_statement"
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
      c("BayesTools_hypothesis_statement", "error", "condition"),
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
  # The engine evaluates a group with the draw-based method directly.
  testthat::local_mocked_bindings(
    hypothesis_BF.default = function(posterior, ...) {

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
      engine_plan_marginal_means(object, statement, parameter = "mu_alloc")
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
    expect_null(BayesTools:::.hypothesis_model_plan_marginal_means_request_refusal(
      plans = request, provider = engine_provider(object)
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


test_that("the plans of a marginal-means object are built once and shared by its methods", {

  BayesTools::hypothesis_cache_clear()
  withr::defer(BayesTools::hypothesis_cache_clear())
  object   <- .marginal_means_route_test_object(model_averaged = FALSE)
  original <- BayesTools:::.hypothesis_model_plan_marginal_means_build
  builds   <- 0L
  testthat::local_mocked_bindings(
    .hypothesis_model_plan_marginal_means_build = function(...) {
      builds <<- builds + 1L
      original(...)
    },
    .package = "BayesTools"
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
    .hypothesis_model_plan_cache = function(provider = NULL) new.env(parent = emptyenv()),
    .package = "BayesTools"
  )
  expect_identical(
    hypothesis(object, "alloc[A] > alloc[B]", columns = "all", seed = 913),
    first
  )
  expect_identical(hypothesis_quantities(object), quantities)
})

test_that("marginal-means statements on several parameters are tested in one call", {

  BayesTools::hypothesis_cache_clear()
  withr::defer(BayesTools::hypothesis_cache_clear())
  call <- function(object, statements, ...) {
    hypothesis(object, statements, columns = "all", seed = 913, ...)
  }
  original <- BayesTools:::.hypothesis_model_plan_marginal_means_build
  builds   <- 0L
  testthat::local_mocked_bindings(
    .hypothesis_model_plan_marginal_means_build = function(...) {
      builds <<- builds + 1L
      original(...)
    },
    .package = "BayesTools"
  )

  # Region statements are available to single-model and model-averaged
  # objects. The posterior odds are the counts of the draws, 75% and 40%
  # above zero.
  for (model_averaged in c(FALSE, TRUE)) {
    object <- .marginal_means_route_two_parameter_object(model_averaged)
    forward <- call(object, c("alloc[A] > 0", "other[A] > 0"))
    expect_identical(rownames(forward), c("mu_alloc", "mu_other"))
    expect_equal(forward[["posterior"]], c(3, 2 / 3), tolerance = 1e-12)
    alone <- list(
      call(object, "alloc[A] > 0"),
      call(object, "other[A] > 0")
    )
    for (column in c("prior", "posterior")) {
      expect_equal(
        forward[[column]],
        c(alone[[1L]][[column]], alone[[2L]][[column]]),
        tolerance = 1e-12
      )
    }
    expect_equal(
      attr(forward, "raw_log_BF"),
      c(attr(alone[[1L]], "raw_log_BF"), attr(alone[[2L]], "raw_log_BF")),
      tolerance = 1e-12
    )

    # The rows keep the order of the statements, and the statement numbers.
    reversed <- call(object, c("other[A] > 0", "alloc[A] > 0"))
    expect_identical(rownames(reversed), c("mu_other", "mu_alloc"))
    expect_equal(reversed[["posterior"]], c(2 / 3, 3), tolerance = 1e-12)
    mixed <- call(object, c("alloc[A] > 0", "other[A] > 0", "alloc[A] > -0.5"))
    expect_identical(rownames(mixed), c("mu_alloc (1)", "mu_other", "mu_alloc (3)"))
    expect_equal(mixed[["posterior"]][c(1L, 2L)], c(3, 2 / 3), tolerance = 1e-12)
    expect_equal(
      attr(mixed, "raw_log_BF")[[3L]],
      attr(call(object, "alloc[A] > -0.5"), "raw_log_BF"),
      tolerance = 1e-12
    )
  }

  # The plans of a request are those of its statements, built for the
  # parameter each statement names: a plan is shared with the call that names
  # the parameter itself.
  object <- .marginal_means_route_two_parameter_object(TRUE)
  BayesTools::hypothesis_cache_clear()
  builds <- 0L
  call(object, c("alloc[A] > 0", "other[A] > 0"))
  expect_identical(builds, 2L)
  call(object, c("other[A] > 0", "alloc[A] > 0"))
  expect_identical(builds, 2L)
  call(object, "other[A] > 0", parameter = "mu_other")
  expect_identical(builds, 2L)
})


test_that("a parameter that does not hold for a statement is refused with the statement and parameters", {

  mismatch <- c(
    "BayesTools_hypothesis_parameter_mismatch", "BayesTools_hypothesis_statement",
    "error", "condition"
  )
  for (model_averaged in c(FALSE, TRUE)) {
    object <- .marginal_means_route_two_parameter_object(model_averaged)
    requests <- list(
      list("other[A] > 0", "alloc"),
      list("other[A] = 0", "mu_alloc"),
      list(c("alloc[A] > 0", "other[A] > 0"), "alloc"),
      list("alloc[A] > other[A]", "alloc"),
      list("alloc[A] > 0", "other")
    )
    for (request in requests) {
      error <- tryCatch(
        hypothesis(object, request[[1L]], parameter = request[[2L]]),
        error = identity
      )
      expect_identical(
        class(error), mismatch,
        info = paste(request[[1L]], request[[2L]])
      )
    }
    error <- tryCatch(
      hypothesis(object, c("alloc[A] > 0", "other[A] > 0"), parameter = "alloc"),
      error = identity
    )
    expect_identical(
      conditionMessage(error),
      paste0(
        "Statement '`other[A]` > 0' refers to the marginal means of ",
        "'mu_other', but 'parameter' selects 'mu_alloc'. A statement must ",
        "refer to the marginal means of one parameter; omit 'parameter' to ",
        "test statements on different parameters in one call."
      )
    )
    # Statements of the selected parameter are as before.
    expect_identical(
      hypothesis(object, "alloc[A] > 0", parameter = "alloc", seed = 913),
      hypothesis(object, "alloc[A] > 0", seed = 913)
    )
  }
})


test_that("a statement on the marginal means of two parameters is refused without 'parameter' advice", {

  for (model_averaged in c(FALSE, TRUE)) {
    object <- .marginal_means_route_two_parameter_object(model_averaged)
    for (statement in c("alloc[A] > other[A]", "alloc[A] > other[A] vs alloc[A] = other[A]")) {
      error <- tryCatch(hypothesis(object, statement), error = identity)
      expect_identical(
        class(error),
        c("BayesTools_hypothesis_ambiguous", "BayesTools_hypothesis_statement",
          "error", "condition"),
        info = statement
      )
      expect_match(
        conditionMessage(error),
        "references multiple marginal-means parameters (mu_alloc, mu_other), but one statement must refer to the marginal means of one parameter",
        fixed = TRUE
      )
      expect_false(grepl("Specify 'parameter'", conditionMessage(error), fixed = TRUE))

      # With 'parameter' the statement still names two parameters, and
      # omitting 'parameter' does not help: the message is the same, under
      # the class of the mismatch with the selected parameter.
      for (parameter in c("alloc", "mu_other")) {
        selected <- tryCatch(
          hypothesis(object, statement, parameter = parameter),
          error = identity
        )
        expect_identical(
          class(selected),
          c("BayesTools_hypothesis_parameter_mismatch", "BayesTools_hypothesis_statement",
            "error", "condition"),
          info = paste(statement, parameter)
        )
        expect_identical(
          conditionMessage(selected), conditionMessage(error),
          info = paste(statement, parameter)
        )
        expect_false(grepl("omit 'parameter'", conditionMessage(selected), fixed = TRUE))
      }
    }
  }
})


test_that("an unknown name in a marginal-means request is named in the refusal", {

  object <- .marginal_means_route_two_parameter_object(model_averaged = FALSE)
  for (statements in list(
    "zz[A] > 0", c("alloc[A] > 0", "zz[A] > 0"), c("zz[A] > 0", "alloc[A] > 0")
  )) {
    error <- tryCatch(hypothesis(object, statements), error = identity)
    expect_identical(
      class(error),
      c("BayesTools_hypothesis_statement", "BayesTools_parameter_not_found",
        "BayesTools_parameter_resolution_error", "error", "condition"),
      info = paste(statements, collapse = " | ")
    )
    expect_identical(error[["alias"]], "zz")
    expect_match(
      conditionMessage(error), "The name 'zz' is not a marginal-means quantity.",
      fixed = TRUE
    )
    expect_true(all(
      c("alloc", "other", "mu_alloc", "mu_other") %in% error[["available"]]
    ))
    expect_false("zz" %in% error[["available"]])
  }

  # several unknown names of one statement are all named
  error <- tryCatch(hypothesis(object, "yy[A] > zz[A]"), error = identity)
  expect_identical(error[["alias"]], "yy")
  expect_match(
    conditionMessage(error),
    "The names 'yy', 'zz' are not marginal-means quantities.",
    fixed = TRUE
  )
  expect_s3_class(error, "BayesTools_parameter_not_found")
})


test_that("model-averaged requests count point levels across parameters", {

  object <- .marginal_means_route_two_parameter_object(model_averaged = TRUE)
  single <- .marginal_means_route_two_parameter_object(model_averaged = FALSE)
  refusal <- function(object, statements) {
    plans <- lapply(statements, function(statement) {
      engine_plan_marginal_means(object, statement)
    })
    BayesTools:::.hypothesis_model_plan_marginal_means_request_refusal(
      plans = plans, provider = engine_provider(object)
    )
  }

  # The same level label of two parameters is two levels.
  levels <- c("alloc[A] = 0", "other[A] = 0")
  expect_identical(
    refusal(object, levels)[["reason"]],
    paste0(
      "Point hypotheses spanning multiple marginal-means levels are not ",
      "supported because the levels need not share one conditioning event."
    )
  )
  expect_error(
    hypothesis(object, levels),
    "spanning multiple marginal-means levels",
    class = "BayesTools_hypothesis_statement"
  )
  # One level of one parameter, and region statements on several parameters,
  # are not refused; points and regions do not mix across parameters either.
  expect_null(refusal(object, c("alloc[A] = 0", "alloc[A] = 0.5")))
  expect_null(refusal(object, c("alloc[A] > 0", "other[A] > 0", "alloc[B] > 0")))
  expect_identical(
    refusal(object, c("alloc[A] = 0", "other[A] > 0"))[["class"]],
    "BayesTools_hypothesis_statement"
  )
  expect_error(
    hypothesis(object, c("alloc[A] = 0", "other[A] > 0")),
    "cannot mix point-null and region statements",
    class = "BayesTools_hypothesis_statement"
  )
  # Single-model objects use one averaged posterior for all of them.
  for (statements in list(levels, c("alloc[A] = 0", "other[A] > 0"))) {
    expect_null(refusal(single, statements))
  }
})

test_that("the settings of a density method are validated after the statements of marginal means", {

  object <- .marginal_means_route_test_object(model_averaged = FALSE)
  # settings of a method that has none
  expect_error(
    hypothesis(object, "alloc[A] > 0", density_method = "KDE",
               density_control = list(foo = 1)),
    "'density_control' is only used when 'density_method' is 'qCMDE' or 'IWMDE'.",
    fixed = TRUE
  )
  # and unrecognized settings of an ordinate method, also for a region
  expect_error(
    hypothesis(object, "alloc[A] > 0", density_method = "qCMDE",
               density_control = list(foo = 1)),
    "'density_control' contains unrecognized setting(s): 'foo'.",
    fixed = TRUE
  )
  # a statement the plans refuse is reported before the settings
  expect_error(
    hypothesis(object, "alloc[A] * alloc[B] = 0", density_method = "KDE",
               density_control = list(foo = 1)),
    class = "BayesTools_hypothesis_statement"
  )
})
