test_that("parallel density chunks preserve rows, budgets, and cleanup", {

  object <- list(marker = "complete fitted object", fit = list(draws = 1))
  attr(object[["fit"]], "runtime_state") <- list(as.raw(1:4))
  samples <- matrix(c(11, 13, 17, 19, 23, 29), ncol = 1L,
    dimnames = list(NULL, "draw_id"))
  z <- c(-1, 0, 2)
  control <- set_selection_likelihood_control(seed = 17L)
  withr::local_options(RoBMA.known_v_covariance_max_bytes = 123456,
    BayesTools.random_effects_memory_limit_bytes = 234567)
  events <- character()
  sampled_rows <- list()
  initialized_limits <- numeric()
  fail_worker <- FALSE
  setup <- function(context){

    events <<- c(events, paste(context$phase, context$role))
    invisible(NULL)
  }
  testthat::local_mocked_bindings(
    .selection_runtime_setup = function() setup,
    .zplot_selection_marginal = function(object, posterior_samples, z_sequence,
        z_threshold, conditioning_depth, integration_control, extrapolate_only){
      expect_identical(object, list(marker = "complete fitted object", fit = list(draws = 1)))
      expect_identical(z_sequence, z)
      expect_null(z_threshold)
      expect_identical(conditioning_depth, "marginal")
      expect_identical(integration_control, control)
      expect_identical(getOption("RoBMA.known_v_covariance_max_bytes"), 123456)
      expect_identical(getOption("BayesTools.random_effects_memory_limit_bytes"), 234567)
      expect_false(extrapolate_only)
      ids <- posterior_samples[, "draw_id"]
      sampled_rows[[length(sampled_rows) + 1L]] <<- ids
      if(fail_worker && 19 %in% ids) stop("posterior draw 19 failed")
      density <- outer(ids, z, "+")
      list(fitted = density, extrapolated = density / 2,
        weights = ids, EDR = NULL)
    },
    .package = "RoBMA"
  )
  testthat::local_mocked_bindings(
    .JAGS_require_packages = function(required_packages, cl, operation){
      expect_identical(operation, "Parallel zplot density computation")
      expect_identical(required_packages, c("BayesTools", "RoBMA"))
      events <<- c(events, "parity")
    },
    .package = "BayesTools"
  )
  testthat::local_mocked_bindings(
    makePSOCKcluster = function(cores, ...) as.list(seq_len(cores)),
    stopCluster = function(cl) events <<- c(events, "stop"),
    clusterCall = function(cl, fun, ...){
      old <- options(RoBMA.known_v_covariance_max_bytes = 1,
        BayesTools.random_effects_memory_limit_bytes = 2)
      on.exit(options(old), add = TRUE)
      lapply(cl, function(node){
        options(RoBMA.known_v_covariance_max_bytes = 1,
          BayesTools.random_effects_memory_limit_bytes = 2)
        value <- fun(...)
        initialized_limits <<- c(initialized_limits, getOption("RoBMA.known_v_covariance_max_bytes"),
          getOption("BayesTools.random_effects_memory_limit_bytes"))
        value
      })
    },
    clusterApply = function(cl, x, fun, ...) lapply(x, fun, ...),
    .package = "parallel"
  )
  set.seed(123)
  seed <- .Random.seed
  expect_message(result <- .zplot_selection_marginal_parallel(object, samples, z,
    "marginal", control, cores = 2L),
    "Computing zplot densities using 2 parallel workers.", fixed = TRUE)
  expect_identical(unlist(sampled_rows, use.names = FALSE), as.numeric(samples))
  expected <- outer(as.numeric(samples), z, "+")
  expect_identical(result$fitted, expected)
  expect_identical(result$extrapolated, expected / 2)
  expect_identical(result$weights, as.numeric(samples))
  expect_identical(events, c("parity", "start coordinator", "start worker",
    "start worker", "stop", "finish coordinator"))
  expect_identical(.Random.seed, seed)
  expect_identical(attr(object[["fit"]], "runtime_state"), list(as.raw(1:4)))
  expect_identical(initialized_limits, c(123456, 234567, 123456, 234567))

  events <- character()
  fail_worker <- TRUE
  expect_error(suppressMessages(.zplot_selection_marginal_parallel(object, samples,
    z, "marginal", control, cores = 2L)), "posterior draw 19 failed", fixed = TRUE)
  expect_identical(tail(events, 2L), c("stop", "finish coordinator"))
})


test_that("density dispatch thins once and keeps summaries and reference curves serial", {

  samples <- matrix(c(11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67), ncol = 1L)
  object <- structure(list(fit = list()), class = c("brma.mv", "brma"))
  calls <- list()
  evaluate <- function(route, posterior_samples, z, cores = NULL){

    calls[[length(calls) + 1L]] <<- list(route = route, samples = posterior_samples, cores = cores)
    density <- outer(as.numeric(posterior_samples), z, "+")
    list(fitted = density, extrapolated = density / 2,
      weights = rep(1, nrow(posterior_samples)), EDR = rep(.1, nrow(posterior_samples)))
  }
  testthat::local_mocked_bindings(
    .get_posterior_samples = function(...) samples,
    .is_weightfunction = function(...) TRUE,
    .effect_direction = function(...) "positive",
    .zplot_requires_selection_marginal = function(...) TRUE,
    .zplot_selection_marginal = function(object, posterior_samples, z_sequence, z_threshold,
        conditioning_depth, integration_control, extrapolate_only = FALSE){
      evaluate("serial", posterior_samples, if(is.null(z_threshold)) z_sequence else z_threshold)
    },
    .zplot_selection_marginal_parallel = function(object, posterior_samples, z_sequence,
        conditioning_depth, integration_control, cores){
      evaluate("parallel", posterior_samples, z_sequence, cores)
    },
    .package = "RoBMA"
  )
  z <- c(-1, 0, 2)
  .zplot_fun.brma(object, z_sequence = z, max_samples = 10, parallel = TRUE, cores = 2L)
  .zplot_density_pair(object, z, max_samples = 10, parallel = TRUE, cores = 2L)
  expect_identical(vapply(calls, `[[`, character(1L), "route"), rep("parallel", 2L))
  expected <- matrix(c(11, 17, 19, 29, 31, 41, 43, 53, 59, 67), ncol = 1L)
  expect_true(all(vapply(calls, function(call) identical(call$samples, expected), logical(1L))))
  expect_identical(vapply(calls, `[[`, integer(1L), "cores"), c(2L, 2L))

  calls <- list()
  .zplot_fun.brma(object, z_threshold = 1.96, max_samples = 10, parallel = TRUE, cores = 2L)
  .zplot_fun.brma(object, z_sequence = z, max_samples = 10, extrapolate = TRUE,
    parallel = TRUE, cores = 2L)
  .zplot_fun.brma(object, z_sequence = z, max_samples = 10, parallel = TRUE, cores = 20L)
  .zplot_density_pair(object, z, max_samples = 10, parallel = TRUE, cores = 1L)
  expect_identical(vapply(calls, `[[`, character(1L), "route"), rep("serial", 4L))
  expect_identical(formals(plot.zplot_brma)$parallel, FALSE)
  expect_identical(formals(lines.zplot_brma)$parallel, FALSE)
  expect_error(plot.zplot_brma(object, parallel = NA), "'parallel'", fixed = TRUE)
  expect_error(lines.zplot_brma(object, cores = 0, as_data = TRUE), "'cores'", fixed = TRUE)
})


test_that("generic QMC z projections retain their design across draw partitions", {

  sei <- c(.8, 1, 1.2)
  prior <- BayesTools::prior_weightfunction("one-sided", .025,
    BayesTools::wf_fixed(c(1, .2)))
  context <- .selection_spec(list(outcome = list(bias = prior)), numeric(3L), sei, "positive")
  context$omega <- matrix(rep(c(1, .2), 3L), 3L, byrow = TRUE)
  context$alpha <- numeric(3L)
  context$phack_kind <- integer(3L)
  context$kernel_mode <- rep(SELKERNEL_STEP, 3L)
  context$vector_rule <- integer(3L)
  context$use_normal <- rep(FALSE, 3L)
  context <- BayesTools::selection_context_validate(context, n_samples = 3L)
  mean <- matrix(c(-.5, .2, .8, .3, -.4, .1, 0, .5, -.5), 3L, byrow = TRUE)
  covariance <- matrix(c(1, .3, -.2, .3, 1, .15, -.2, .15, 1), 3L)
  packed <- matrix(covariance[lower.tri(covariance, diag = TRUE)],
    3L, 6L, byrow = TRUE)
  z <- c(-2, .3, 2.2)
  control <- set_selection_likelihood_control(points_per_scramble = 32L,
    max_points_per_scramble = 4096L, scrambles = 8L, seed = 17L)
  designs <- new.env(parent = emptyenv())
  set.seed(123)
  seed <- .Random.seed
  full <- .zplot_joint_block(z, mean, packed, sei, context, FALSE, control, designs)
  pieces <- lapply(list(1L, 2:3), function(rows){
    .zplot_joint_block(z, mean[rows, , drop = FALSE], packed[rows, , drop = FALSE], sei,
      BayesTools::selection_context_subset_rows(context, rows), FALSE, control)
  })
  expect_gt(length(ls(designs)), 1L)
  expect_identical(do.call(rbind, lapply(pieces, `[[`, "density")), full$density)
  expect_identical(unlist(lapply(pieces, `[[`, "relative_mcse"), use.names = FALSE),
    full$relative_mcse)
  expect_identical(unlist(lapply(pieces, `[[`, "log_density"), use.names = FALSE),
    full$log_density)

  # In projection mode the historical log_density field already contains
  # -log(A), rather than a Gaussian likelihood. A distant, finite dummy
  # observation therefore cannot introduce Gaussian cancellation into either
  # the returned normalizer or its refinement comparison.
  native <- get(".Call", envir = environment(.zplot_joint_block))
  project <- .zplot_joint_block
  environment(project) <- new.env(parent = environment(project))
  environment(project)$.Call <- function(.NAME, ..., PACKAGE = NULL){

    arguments <- list(...)
    arguments[[1L]] <- rep(1e8, length(arguments[[1L]]))
    do.call(native, c(list(.NAME), arguments, list(PACKAGE = PACKAGE)))
  }
  shifted <- project(z, mean, packed, sei, context, FALSE, control)
  expect_identical(shifted, full)
  expect_identical(.Random.seed, seed)
})
