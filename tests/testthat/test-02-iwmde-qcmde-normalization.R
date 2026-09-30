context("qCMDE normalization on fitted models")

source(testthat::test_path("common-functions.R"))
source(testthat::test_path("helper-iwmde.R"))


# Fitted models whose rows cannot report a Gaussian kernel: a normal selection
# model in which every row carries the weight function, and binomial GLMMs.
.qcmde_merge_fit_names <- function() {

  utils::head(intersect(
    c("dat.lehmann2018-3PSM", "bcg_glmm", "nielweise2008_glmm"),
    list_fits()
  ), 2L)
}


# The context, the rows and row states, the replacement, and the normalization
# grid of the location target of a fitted model, for the first `n_rows` active
# rows.
.qcmde_fit_inputs <- function(fit_name, n_rows = 6L, n_points = 15L) {

  context   <- .iwmde_context(load_fit(fit_name, validate = FALSE))
  columns   <- colnames(context[["posterior_samples"]])
  parameter <- intersect(c("mu", "mu_intercept"), columns)[[1L]]
  spec      <- .iwmde_parameter_spec(context, parameter, NULL)
  values    <- .iwmde_parameter_values(context, parameter, spec)
  component <- .iwmde_parameter_components(context, parameter, spec)
  active    <- component[["active"]] & is.finite(values)
  rows      <- utils::head(which(active), n_rows)
  states    <- .iwmde_row_states(context, rows, parameter, spec)
  keep      <- vapply(states, function(state) {
    is.finite(state[["baseline_log_q"]])
  }, logical(1))
  rows      <- rows[keep]
  states    <- states[keep]
  support   <- .iwmde_parameter_support(context, parameter, rows, spec)
  transform <- .iwmde_parameter_transform(support)
  grid      <- .iwmde_normalization_grid(
    values               = values[active],
    display_grid         = numeric(),
    support              = support,
    transform            = transform,
    normalization_points = n_points,
    normalization_prob   = .99
  )
  quantiles <- as.numeric(stats::quantile(
    values[active], probs = seq(.2, .8, length.out = 7L), names = FALSE,
    type = 8
  ))

  list(
    context     = context,
    parameter   = parameter,
    rows        = rows,
    row_states  = states,
    replacement = .iwmde_replacement_spec(context, parameter, spec),
    transform   = transform,
    grid        = grid,
    values      = quantiles
  )
}


.qcmde_plain <- function(x) {

  matrix(as.numeric(x), nrow(x), ncol(x))
}


test_that("the joint density of a call is that of the same values in separate calls", {

  fit_names <- .qcmde_merge_fit_names()
  skip_if(length(fit_names) == 0L, "No cached selection or GLMM fixture is active.")

  for (fit_name in fit_names) {
    input <- .qcmde_fit_inputs(fit_name)
    first  <- input[["values"]][1:3]
    second <- input[["values"]][4:7]
    log_q  <- function(values) {
      .qcmde_plain(.iwmde_log_q_grid(
        input[["context"]], input[["parameter"]], values,
        input[["row_states"]], input[["replacement"]]
      ))
    }
    joint <- log_q(c(first, second))
    expect_true(all(is.finite(joint)), info = fit_name)
    expect_identical(joint[1:3, , drop = FALSE], log_q(first), info = fit_name)
    expect_identical(joint[4:7, , drop = FALSE], log_q(second), info = fit_name)
    expect_identical(joint[c(6L, 2L, 7L), , drop = FALSE],
                     log_q(input[["values"]][c(6L, 2L, 7L)]), info = fit_name)
  }
})


test_that("an estimate-only pass evaluated in one call assembles what separate calls do", {

  fit_names <- .qcmde_merge_fit_names()
  skip_if(length(fit_names) == 0L, "No cached selection or GLMM fixture is active.")

  original <- .iwmde_log_q_grid
  for (fit_name in fit_names) {
    input <- .qcmde_fit_inputs(fit_name)
    ordinary <- .iwmde_qcmde_ordinary_rows(input[["row_states"]])
    expect_true(.iwmde_qcmde_estimate_only(
      input[["context"]], input[["row_states"]], ordinary
    ), info = fit_name)

    calls <- 0L
    pass  <- function(estimate_only) {
      calls <<- 0L
      testthat::with_mocked_bindings(
        .iwmde_qcmde_normalization_pass(
          context            = input[["context"]],
          parameter          = input[["parameter"]],
          display_grid       = input[["values"]][3:4],
          normalization_grid = input[["grid"]],
          transform          = input[["transform"]],
          normalization_prob = .99,
          row_states         = input[["row_states"]],
          replacement        = input[["replacement"]],
          estimator_rows     = input[["rows"]]
        ),
        .iwmde_log_q_grid = function(...) {
          calls <<- calls + 1L
          original(...)
        },
        .iwmde_qcmde_estimate_only = function(context, row_states, rows) {
          estimate_only
        },
        .package = "RoBMA"
      )
    }
    merged     <- pass(TRUE)
    merged_calls <- calls
    sequential <- pass(FALSE)
    sequential_calls <- calls

    expect_lt(merged_calls, sequential_calls, label = fit_name)
    merged[["log_q_display"]]     <- .qcmde_plain(merged[["log_q_display"]])
    sequential[["log_q_display"]] <- .qcmde_plain(sequential[["log_q_display"]])
    expect_identical(merged, sequential, info = fit_name)
  }
})


test_that("a cached normalization gives the statements the results of the uncached path", {

  fit_names <- .qcmde_merge_fit_names()
  skip_if(length(fit_names) == 0L, "No cached selection or GLMM fixture is active.")

  original   <- .iwmde_qcmde_normalization_pass
  statements <- c("mu = 0", "mu = 0.15", "mu = -0.05")
  for (fit_name in fit_names) {
    fit     <- load_fit(fit_name, validate = FALSE)
    context <- .iwmde_context(fit)
    env     <- .iwmde_qcmde_cache_env(context)
    expect_true(is.environment(env), info = fit_name)
    passes <- 0L
    run <- function(cache) {
      passes <<- 0L
      rm(list = ls(env, all.names = TRUE), envir = env)
      testthat::with_mocked_bindings(
        {
          results <- lapply(statements, function(statement) {
            hypothesis(fit, statement, density_method = "qCMDE")
          })
          # The density curve has the display grid of a plot and the row budget
          # of a density, so its normalization is its own; a second request for
          # the curve reuses it.
          curve <- function() {
            .iwmde_estimate(
              context         = .iwmde_context(fit),
              parameter       = "mu",
              density_method  = "qCMDE",
              density_control = NULL,
              outputs         = "density"
            )[["density"]]
          }
          curves <- list(curve(), curve())
          list(results = results, curves = curves, passes = passes)
        },
        .iwmde_qcmde_normalization_pass = function(...) {
          passes <<- passes + 1L
          original(...)
        },
        .iwmde_qcmde_cache_env = if (cache) {
          function(context) env
        } else {
          function(context) NULL
        },
        .package = "RoBMA"
      )
    }
    cached   <- run(TRUE)
    uncached <- run(FALSE)

    # One pass serves the statements, and one the requests of the density curve.
    expect_identical(cached$passes, 2L, info = fit_name)
    expect_identical(uncached$passes, length(statements) + 2L, info = fit_name)
    for (i in seq_along(statements)) {
      expect_identical(as.data.frame(cached$results[[i]]),
                       as.data.frame(uncached$results[[i]]),
                       info = paste(fit_name, statements[[i]]))
      expect_identical(density_diagnostics(cached$results[[i]]),
                       density_diagnostics(uncached$results[[i]]),
                       info = paste(fit_name, statements[[i]]))
    }
    for (i in 1:2) {
      expect_identical(cached$curves[[i]][["y"]], uncached$curves[[i]][["y"]],
                       info = fit_name)
      expect_identical(cached$curves[[i]], uncached$curves[[i]], info = fit_name)
    }
  }
})


test_that("the key of a normalization names the fit, target, rows and control, and not the value tested", {

  fit_names <- .qcmde_merge_fit_names()
  skip_if(length(fit_names) < 2L, "Two cached selection or GLMM fixtures are needed.")

  key <- function(context, value = 0, parameter = "mu", control = NULL) {
    plan <- .iwmde_plan(
      context         = context,
      parameter       = parameter,
      density_method  = "qCMDE",
      density_control = control,
      outputs         = "ordinate",
      values          = value
    )
    .iwmde_plan_normalization_cache(context, plan, plan[["grids"]][["display_grid"]])
  }
  first  <- .iwmde_context(load_fit(fit_names[[1L]], validate = FALSE))
  second <- .iwmde_context(load_fit(fit_names[[2L]], validate = FALSE))
  base   <- key(first)

  expect_false(base[["display_in_key"]])
  expect_identical(key(first, .15)[["key"]], base[["key"]])
  expect_identical(key(first, -.3)[["key"]], base[["key"]])
  # Another fit, another target of the fit, and another control are others.
  expect_false(identical(key(second)[["key"]], base[["key"]]))
  scale <- intersect(c("tau", "log_tau_intercept"),
                     colnames(first[["posterior_samples"]]))
  if (length(scale) > 0L) {
    expect_false(identical(
      key(first, .3, parameter = scale[[1L]])[["key"]], base[["key"]]
    ))
  }
  for (control in list(list(normalization_prob = .99),
                       list(normalization_points = 40L),
                       list(samples = 50L))) {
    expect_false(identical(key(first, control = control)[["key"]], base[["key"]]),
                 info = paste(names(control), control[[1L]]))
  }
  # One cache serves the session; the keys keep the fits apart.
  expect_identical(base[["cache"]], key(first, .2)[["cache"]])
  expect_identical(base[["cache"]], key(second)[["cache"]])
})
