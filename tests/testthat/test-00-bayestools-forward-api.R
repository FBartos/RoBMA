context("BayesTools forward API integration")

if (!exists("fit_cache_paths", mode = "function")) {
  source(testthat::test_path("common-functions.R"))
}

.read_cached_fit_with_formula_design <- function(name, parameter = "mu") {

  paths <- unique(c(
    fit_cache_paths(name)[["fit"]],
    testthat::test_path("test_files", "fits", paste0(name, ".RDS"))
  ))

  for (path in paths[file.exists(paths)]) {
    object <- suppressWarnings(try(readRDS(path), silent = TRUE))
    if (inherits(object, "try-error") || is.null(object[["fit"]])) {
      next
    }

    design <- BayesTools::JAGS_formula_design(object[["fit"]], parameter)
    if (!is.null(design)) {
      return(list(object = object, design = design))
    }
  }

  skip(paste0("Cached fit '", name, "' with formula design metadata is unavailable."))
}

test_that("BayesTools forward API guard passes for the active namespace", {

  expect_true(isTRUE(.check_bayestools_forward_api()))
  expect_true(all(
    c("fit", "parameter") %in% names(formals(BayesTools::JAGS_formula_design))
  ))
  expect_true("interpret_records" %in% getNamespaceExports("BayesTools"))
  expect_true(all(
    c("legend", "legend_title", "legend_labels", "legend_position", "density_method") %in%
      names(formals(BayesTools::plot_marginal))
  ))
  expect_true("density_method" %in% names(formals(BayesTools::plot_posterior)))
  expect_true("density_method" %in% names(formals(BayesTools::Savage_Dickey_BF)))
})

test_that(".get_model_matrix uses fitted BayesTools design metadata", {

  outcome <- data.frame(
    yi  = c(0.10, 0.20, 0.30),
    sei = c(0.20, 0.30, 0.40)
  )
  data <- list(
    outcome = outcome,
    mods    = data.frame(x = c(1, 2, 3))
  )
  attr(data, "mods")             <- TRUE
  attr(data, "effect_direction") <- "negative"

  X_design <- cbind(
    mu_intercept    = 1,
    x__xXx__factorB = c(0, 1, 0)
  )
  attr(X_design, "assign") <- c(0L, 1L)

  design <- list(
    model_matrix     = X_design,
    column_names     = colnames(X_design),
    raw_column_names = c("(Intercept)", "x:factorB"),
    assign           = c(0L, 1L),
    rank             = 2L,
    aliased          = stats::setNames(c(FALSE, FALSE), colnames(X_design))
  )
  fit <- structure(list(), class = "BayesTools_fit")
  attr(fit, "formula_design") <- list(mu = design)

  object <- list(
    data   = data,
    fit    = fit,
    priors = list(
      outcome = list(
        bias = BayesTools::prior_PET("normal", parameters = list(mean = 0, sd = 1))
      )
    )
  )

  X <- .get_model_matrix(object)

  expect_equal(as.numeric(X[, colnames(X_design)]), as.numeric(X_design))
  expect_equal(X[, "PET"], -outcome[["sei"]])
  expect_equal(attr(X, "assign"), c(0L, 1L, 2L))
  expect_equal(attr(X, "rank"), 3L)
  expect_equal(
    attr(X, "raw_colnames"),
    c(mu_intercept = "(Intercept)", x__xXx__factorB = "x:factorB", PET = "PET")
  )
})

test_that(".get_model_matrix preserves BayesTools design metadata from real fits", {

  cached <- .read_cached_fit_with_formula_design("bcg_meta-regression3")
  object <- cached[["object"]]
  design <- cached[["design"]]
  X      <- .get_model_matrix(object)

  expect_s3_class(design, "BayesTools_formula_design")
  expect_equal(dim(X), dim(design[["model_matrix"]]))
  expect_equal(as.numeric(X), as.numeric(design[["model_matrix"]]))
  expect_equal(colnames(X), design[["column_names"]])
  expect_equal(attr(X, "assign"), design[["assign"]])
  expect_equal(attr(X, "rank"), design[["rank"]])
  expect_equal(
    attr(X, "raw_colnames"),
    stats::setNames(design[["raw_column_names"]], design[["column_names"]])
  )
  expect_true(any(grepl("__xXx__", colnames(X), fixed = TRUE)))
  expect_true(any(grepl(":", attr(X, "raw_colnames"), fixed = TRUE)))
})

test_that("formula design accessors support current non-fitted objects", {

  dat <- data.frame(
    yi  = c(.10, .20, .15, .30),
    sei = c(.10, .12, .11, .15),
    x   = c(-1, 0, 1, 2),
    z   = factor(c("a", "b", "a", "b"))
  )

  prior_object <- BMA(
    yi           = yi,
    sei          = sei,
    mods         = ~ x + z,
    data         = dat,
    measure      = "SMD",
    only_priors  = TRUE
  )

  prior_design <- .fitted_formula_design(prior_object, "mu", required = TRUE)

  expect_s3_class(prior_design, "BayesTools_formula_design")
  expect_true(all(c("intercept", "x", "z") %in% prior_design[["model_terms"]]))
  expect_true("mu_x" %in% names(prior_design[["prior_list"]]))

  data_object <- brma.norm(
    yi        = yi,
    sei       = sei,
    mods      = ~ x + z,
    data      = dat,
    only_data = TRUE
  )

  data_design <- .fitted_formula_design(data_object, "mu", required = TRUE)

  expect_s3_class(data_design, "BayesTools_formula_design")
  expect_equal(nrow(data_design[["model_matrix"]]), nrow(dat))
  expect_equal(data_design[["assign"]], attr(data_design[["model_matrix"]], "assign"))
  expect_true(all(c("intercept", "x", "z") %in% data_design[["model_terms"]]))
})

test_that("Weightfunction observed p-values use selection mapping sign", {

  outcome <- data.frame(
    yi  = c(0.10, -0.20, 0.30),
    sei = c(0.20, 0.25, 0.50)
  )
  data <- list(outcome = outcome)
  attr(data, "effect_direction") <- "negative"

  object <- list(
    data   = data,
    priors = list(
      outcome = list(
        bias = BayesTools::prior_weightfunction(
          "one-sided",
          c(0.05),
          weights = BayesTools::wf_fixed(c(1, 0.5))
        )
      )
    )
  )

  expect_equal(
    .weightfunction_observed_p_values(object),
    stats::pnorm(-outcome[["yi"]] / outcome[["sei"]], lower.tail = FALSE)
  )
})
