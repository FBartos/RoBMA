.review_grid_diagnostic_transport <- function(method, mcse, change, kind, consumer = "grid",
    cache_change = NULL) {

  # Controlled diagnostic transport only: retain the actual native numerical
  # values and all original wrappers/grid/cache, replace reported diagnostics
  # immediately after .Call, and restore the exact wrapper afterwards.
  name <- paste0(".selection_joint_", if (method == "rank_one") "cluster" else "factor", "_loglik_block")
  namespace <- asNamespace("RoBMA")
  original <- get(name, namespace)
  pieces <- as.list(body(original))
  at <- which(vapply(pieces, function(x) is.call(x) && identical(x[[1L]], as.name("<-")) &&
    identical(x[[2L]], as.name("result")) && is.call(x[[3L]]) && identical(x[[3L]][[1L]], as.name(".Call")), logical(1)))
  stopifnot(length(at) == 1L)
  trace(name, where = namespace, at = at + 1L, print = FALSE,
    tracer = substitute({ result[["relative_mcse"]][] <- MCSE; result[["relative_change"]][] <- CHANGE },
      list(MCSE = mcse, CHANGE = change)))
  on.exit({
    untrace(name, where = namespace)
    stopifnot(identical(body(get(name, namespace)), body(original)), identical(formals(get(name, namespace)), formals(original)))
  }, add = TRUE)
  yi <- c(.8, 1.2); sei <- c(.4, .6); mean <- c(-.1, .2)
  residual <- c(.2, .25)
  loading <- if (method == "rank_one") matrix(c(.8, .8), 2L) else matrix(c(.8, .8, .3, -.3), 2L)
  covariance <- diag(residual^2) + tcrossprod(loading)
  prior <- BayesTools::prior_weightfunction("one-sided", .05, BayesTools::wf_fixed(c(1, .1)))
  object <- bselmodel.mv(yi = yi, V = covariance, measure = "GEN", prior_unit_information_sd = 1,
    prior_bias = prior, only_priors = TRUE, silent = TRUE)
  queries <- seq(-.1, .1, length.out = if (kind == "anchor") 81L else 16L)
  queries <- c(queries, queries[1:2])
  S <- length(queries)
  shared <- .selection_normalizer_grid(list(data = object$data), queries, 1L)
  if (kind == "fallback") shared$geometry_limit <- 0
  metadata <- list(sign = 1L, rows = 1:2, block = 1L,
    state = list(shared = shared, values = queries, state_index = rep(1L, S), rows = 1L,
      mean = matrix(mean, 1L), basis = matrix(if (kind == "constant") c(0, 0) else c(1, .5), 1L), current = 0))
  context <- .selection_spec(list(outcome = list(bias = prior)), yi, sei, "positive")
  context$omega <- matrix(c(1, .1), S, 2L, byrow = TRUE)
  context$kernel_mode <- rep(1L, S); context$vector_rule <- integer(S)
  context$alpha <- numeric(S); context$phack_kind <- integer(S); context$use_normal <- rep(FALSE, S)
  # Diagnostic transport needs finite real native values. Rank-one product
  # fallback now refuses this deliberately tiny budget, so retain its accepted
  # ordinary ladder; higher-rank randomized coverage keeps the short schedule.
  quadrature <- .selection_joint_cluster_quadrature_rules(if (method == "rank_one")
    SELNORM_CLUSTER_QUADRATURE_ORDERS else c(1L, 3L, 5L))
  rank <- ncol(loading)
  plan <- list(points_per_scramble = 8L, max_points_per_scramble = 8L,
    factor_points_per_proposal = 8L, factor_max_points_per_proposal = 8L,
    scrambles = 2L, relative_tolerance = .005, seed = 1L,
    quadrature = quadrature, factor_quadrature = quadrature, factor_ranks = rank,
    block_methods = method,
    design_keys = paste0("factor_", rank),
    designs = setNames(list(BayesTools::selection_qmc_design(2L * rank, 8L, 2L, seed = 1L)), paste0("factor_", rank)))
  means <- matrix(mean, S, 2L, byrow = TRUE)
  if (kind != "constant") means <- means + outer(queries, c(1, .5))
  call <- function() {
    if (consumer == "zplot") return(.zplot_normalizer_block(yi, means, NULL, sei, context,
      plan, plan, 1L, new.env(parent = emptyenv()),
      list(residual_sd = matrix(residual, S, 2L, byrow = TRUE),
        loading = matrix(as.numeric(loading), S, length(loading), byrow = TRUE))))
    args <- list(yi, means, matrix(residual, S, 2L, byrow = TRUE),
      matrix(as.numeric(loading), S, length(loading), byrow = TRUE), sei, context, plan)
    if (method == "factor") args$block_index <- 1L
    if (consumer == "raw") args$return_normalizer <- TRUE else args$normalizer_grid <- metadata
    do.call(get(name, namespace), args)
  }
  first <- tryCatch(call(), error = identity)
  if (!is.null(cache_change)) {
    for (entry in as.list(shared$cases)) entry$relative_change[entry$have] <- cache_change
  }
  list(first = first, repeat_value = tryCatch(call(), error = identity), shared = shared)
}

test_that("F3-3 every direct factor grid obeys native diagnostics including cached rows", {

  remedy <- paste0(". Diagnostics cover deterministic refinement or randomized integration; ",
    "rank-one product fallback uses a deterministic work budget. ",
    "Increase 'max_points_per_scramble' or 'scrambles' in ",
    "'set_selection_likelihood_control()'. Pass the control as 'selection_control' when fitting, ",
    "'integration_control' in zplot(), or 'density_control$integration_control' for supported post-fit densities and hypotheses.")
  for (method in c("rank_one", "factor")) for (kind in c("constant", "fallback", "anchor")) {
    for (metrics in list(c(0, .02), c(0, NaN), c(0, Inf), c(.02, 0))) {
      result <- .review_grid_diagnostic_transport(method, metrics[1L], metrics[2L], kind)
      message <- paste0("Selection ", if (method == "rank_one") "cluster" else "factor",
        " normalizer was rejected by diagnostics: relative Monte Carlo standard error was ",
        format(metrics[1L], digits = 4), " and refinement/design relative change was ", format(metrics[2L], digits = 4),
        " (accepted when 2.36 times the standard error and the change are within 'relative_tolerance')", remedy)
      for (value in result[c("first", "repeat_value")]) {
        expect_s3_class(value, "error")
        if (inherits(value, "error")) expect_identical(conditionMessage(value), message)
      }
      expect_true(all(vapply(as.list(result$shared$cases), function(entry) !any(entry$have), logical(1))))
    }
    passing <- .review_grid_diagnostic_transport(method, 0, 0, kind)
    expect_true(is.numeric(passing$first) && all(is.finite(passing$first)))
    expect_identical(passing$first, passing$repeat_value)
    if (kind == "anchor") {
      expect_gt(passing$shared$stats[["interpolated_points"]], 0)
      # Native acceptance alone does not certify interpolation: existing
      # anchors additionally require MCSE exactly zero and finite nonnegative
      # eta = 2 * nested change no larger than the original tolerance.
      for (metrics in list(c(.001, 0), c(0, .004))) {
        direct <- .review_grid_diagnostic_transport(method, metrics[1L], metrics[2L], kind)
        expect_true(is.numeric(direct$first) && all(is.finite(direct$first)))
        expect_identical(direct$shared$stats[["interpolated_points"]], 0)
      }
    }
  }
})

test_that("F3-3 cached direct rows are revalidated independently of fresh storage", {

  for (method in c("rank_one", "factor")) {
    result <- .review_grid_diagnostic_transport(method, 0, 0, "constant", cache_change = .02)
    expect_true(is.numeric(result$first) && all(is.finite(result$first)))
    expect_s3_class(result$repeat_value, "error")
    expect_match(conditionMessage(result$repeat_value), "refinement/design relative change was 0.02", fixed = TRUE)
    expect_equal(result$shared$stats[["constant_evaluations"]], 1)
  }
})

test_that("F3-3 raw normalizers retain downstream zplot diagnostic ownership", {

  for (method in c("rank_one", "factor")) {
    raw <- .review_grid_diagnostic_transport(method, 0, .02, "constant", consumer = "raw")
    expect_type(raw$first, "list")
    expect_false(inherits(raw$first, "error"))
    if (!inherits(raw$first, "error")) {
      expect_true(all(raw$first$relative_mcse == 0))
      expect_true(all(raw$first$relative_change == .02))
      expect_true(all(is.finite(raw$first$log_normalizer)))
    }
    zplot <- .review_grid_diagnostic_transport(method, 0, .02, "constant", consumer = "zplot")
    expect_s3_class(zplot$first, "error")
    expect_identical(conditionMessage(zplot$first), paste0(
      "Zplot marginal integration was rejected by diagnostics: relative integration error was 0.02",
      " (accepted when 2.36 times the standard error and the change are within 'relative_tolerance'). ",
      "Diagnostics cover deterministic refinement or randomized integration; ",
      "rank-one product fallback uses a deterministic work budget. ",
      "Increase 'max_points_per_scramble' in 'integration_control = set_selection_likelihood_control()'."))
    passing <- .review_grid_diagnostic_transport(method, 0, 0, "constant", consumer = "zplot")
    expect_false(inherits(passing$first, "error"))
    expect_true(all(is.finite(passing$first$log_density)))
    expect_true(all(passing$first$relative_mcse == 0))
  }
})

test_that("the original tiny rank-one grid schedule remains an explicit refusal", {

  # Same normalizer inputs and controls as the original diagnostic transport
  # fixture. Its baseline QMC result already failed diagnostics (about24.5%
  # error); production quadrature is used above only to supply accepted values.
  block <- .test_cert_block(list(
    .test_cert_parallel_group(1L, .8, .2, -.1, .4),
    .test_cert_parallel_group(1L, .8, .25, .2, .6)), c(1, .1),
    c(stats::qnorm(.05, lower.tail = FALSE), -Inf))
  result <- .test_cert_native_cluster(block, initial = 8L, points = 8L, scrambles = 2L,
    rules = .selection_joint_cluster_quadrature_rules(c(1L, 3L, 5L)))
  expect_true(is.nan(result$log_normalizer))
  expect_identical(result$log_density, -Inf)
  expect_identical(result$relative_mcse, Inf)
  expect_identical(result$relative_change, Inf)
  expect_false(.test_cert_deterministic(result))
  expect_error(.selection_joint_factor_check_diagnostics(result,
    list(relative_tolerance = .005), "cluster"), "deterministic work budget", fixed = TRUE)
})
