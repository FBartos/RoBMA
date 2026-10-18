.review_grid_diagnostic_transport <- function(method, mcse, change, kind) {

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
  quadrature <- .selection_joint_cluster_quadrature_rules(c(1L, 3L, 5L))
  rank <- ncol(loading)
  plan <- list(points_per_scramble = 8L, max_points_per_scramble = 8L,
    factor_points_per_proposal = 8L, factor_max_points_per_proposal = 8L,
    scrambles = 2L, relative_tolerance = .005, seed = 1L,
    quadrature = quadrature, factor_quadrature = quadrature, factor_ranks = rank,
    design_keys = paste0("factor_", rank),
    designs = setNames(list(BayesTools::selection_qmc_design(2L * rank, 8L, 2L, seed = 1L)), paste0("factor_", rank)))
  means <- matrix(mean, S, 2L, byrow = TRUE)
  if (kind != "constant") means <- means + outer(queries, c(1, .5))
  call <- function() {
    args <- list(yi, means, matrix(residual, S, 2L, byrow = TRUE),
      matrix(as.numeric(loading), S, length(loading), byrow = TRUE), sei, context, plan)
    if (method == "factor") args$block_index <- 1L
    args$normalizer_grid <- metadata
    do.call(get(name, namespace), args)
  }
  list(first = tryCatch(call(), error = identity), repeat_value = tryCatch(call(), error = identity),
    shared = shared)
}

test_that("F3-3 every direct factor grid obeys native diagnostics including cached rows", {

  remedy <- paste0(". Increase 'max_points_per_scramble' or 'scrambles' in ",
    "'set_selection_likelihood_control()'. Pass the control as 'selection_control' when fitting, ",
    "'integration_control' in zplot(), or 'density_control$integration_control' for supported post-fit densities and hypotheses.")
  for (method in c("rank_one", "factor")) for (kind in c("constant", "fallback", "anchor")) {
    for (metrics in list(c(0, .02), c(0, NaN), c(0, Inf), c(.02, 0))) {
      result <- .review_grid_diagnostic_transport(method, metrics[1L], metrics[2L], kind)
      message <- paste0("Selection ", if (method == "rank_one") "cluster" else "factor",
        " normalizer was rejected by diagnostics: relative Monte Carlo standard error was ",
        format(metrics[1L], digits = 4), " and nested-design relative change was ", format(metrics[2L], digits = 4), remedy)
      for (value in result[c("first", "repeat_value")]) {
        expect_s3_class(value, "error")
        if (inherits(value, "error")) expect_identical(conditionMessage(value), message)
      }
    }
    passing <- .review_grid_diagnostic_transport(method, 0, 0, kind)
    expect_true(is.numeric(passing$first) && all(is.finite(passing$first)))
    expect_identical(passing$first, passing$repeat_value)
  }
})
