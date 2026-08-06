.marginal_means_route_test_parameter <- function(a, b, condition_keys = NULL) {

  levels <- list(
    A = structure(
      a,
      class          = c("marginal_posterior.simple", "numeric"),
      linear_weights = c(a = 1, b = 0)
    ),
    B = structure(
      b,
      class          = c("marginal_posterior.simple", "numeric"),
      linear_weights = c(a = 0, b = 1)
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
