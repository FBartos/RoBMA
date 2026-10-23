context("Selection QMC stopping multiplier")
skip_on_cran()

# A randomized (QMC) selection normalizer is accepted only when
#   max(z * relative_mcse, relative_change) <= relative_tolerance,   z = 2.36,
# by the native stop of the factor and cluster fallback and by every judge of a
# returned QMC diagnostic: the JAGS wrappers of the factor, cluster, dense
# covariance and conditioned-sampling normalizers and their R counterparts. The
# deferred geometries below reach the fallback, and a single stage (initial =
# maximum budget) fixes its diagnostics, so a tolerance can be placed in the band
#   max(mcse, change) <= tolerance < z * mcse
# where the one-standard-error rule accepted and the scaled rule must not. Every rejection
# message reports the raw values and names the rule in a clause (selnorm-qmc-z.h and
# .selection_qmc_criterion()); the tests pin its wording.

.test_qmc_native <- function(case, seed, initial, points, tolerance) {

  if (case[["kind"]] == "factor") {
    .test_cert_native_factor(case[["block"]], seed = seed, initial = initial,
                             points = points, tolerance = tolerance)
  } else {
    .test_cert_native_cluster(case[["block"]], seed = seed, initial = initial,
                              points = points, tolerance = tolerance)
  }
}

# The single-stage diagnostics of a case and the tolerance band they define.
.test_qmc_band <- function(name, seed) {

  case <- .test_cert_proposal_cases()[[name]]
  one  <- .test_qmc_native(case, seed, 256L, 256L, .005)
  z    <- .selection_qmc_z()
  list(
    case = case, seed = seed, one = one,
    mcse = one[["relative_mcse"]], change = one[["relative_change"]],
    low  = max(one[["relative_mcse"]], one[["relative_change"]]),
    high = z * one[["relative_mcse"]],
    inside = sqrt(max(one[["relative_mcse"]], one[["relative_change"]]) * z * one[["relative_mcse"]])
  )
}

# Fixtures whose single-stage diagnostics leave a band: a depth-one forest, a rank-two
# tensor and two rank-one clusters (seeds of the fixed QMC designs).
.test_qmc_band_specs <- list(
  list(name = "crossing_20", seed = 3L),
  list(name = "tensor_77",   seed = 1L),
  list(name = "cluster_733", seed = 1L),
  list(name = "cluster_189", seed = 2L)
)

test_that("the stopping multiplier is one native constant, the t(7) quantile for eight scrambles", {

  z <- .selection_qmc_z()
  expect_identical(z, 2.36)
  expect_lt(abs(z - stats::qt(.975, 7)), .01)
})

test_that("the native stop continues while the scaled error exceeds the tolerance", {

  for (spec in .test_qmc_band_specs) {
    info <- paste(spec[["name"]], "seed", spec[["seed"]])
    band <- .test_qmc_band(spec[["name"]], spec[["seed"]])
    # fixture: the randomized fallback, with a band where one standard error would stop
    expect_false(.test_cert_deterministic(band[["one"]]), info = info)
    expect_gt(band[["mcse"]], 0, label = info)
    expect_lt(band[["low"]], band[["high"]], label = info)

    # at a tolerance inside the band the first stage is the same fallback result
    single <- .test_qmc_native(band[["case"]], spec[["seed"]], 256L, 256L, band[["inside"]])
    expect_false(.test_cert_deterministic(single), info = info)
    expect_identical(single[["relative_mcse"]], band[["mcse"]], info = info)
    expect_identical(single[["log_normalizer"]], band[["one"]][["log_normalizer"]], info = info)

    # the one-standard-error rule would stop at the first stage there; the scaled
    # rule goes on to more points (a different value with a different spread)
    long <- .test_qmc_native(band[["case"]], spec[["seed"]], 256L, 4096L, band[["inside"]])
    expect_false(.test_cert_deterministic(long), info = info)
    expect_false(identical(long[["relative_mcse"]], band[["mcse"]]), info = info)
    expect_false(identical(long[["log_normalizer"]], band[["one"]][["log_normalizer"]]), info = info)

    # at the boundary z * mcse = tolerance the stage is accepted, with the raw MCSE
    at <- .test_qmc_native(band[["case"]], spec[["seed"]], 256L, 4096L, band[["high"]])
    expect_false(.test_cert_deterministic(at), info = info)
    expect_identical(at[["relative_mcse"]], band[["mcse"]], info = info)
    expect_identical(at[["relative_change"]], band[["change"]], info = info)
    expect_identical(at[["log_normalizer"]], band[["one"]][["log_normalizer"]], info = info)

    # a hair below the boundary it is not
    below <- .test_qmc_native(band[["case"]], spec[["seed"]], 256L, 4096L, band[["high"]] * (1 - 1e-9))
    expect_false(identical(below[["relative_mcse"]], band[["mcse"]]), info = info)
  }
})

test_that("a budget-exhausted fallback is judged by the scaled error in R", {

  z <- .selection_qmc_z()
  mcse <- .003
  change <- .001
  tolerance <- .005
  expect_lte(mcse, tolerance)
  expect_gt(z * mcse, tolerance)
  band <- list(relative_mcse = mcse, relative_change = change)
  plan <- list(relative_tolerance = tolerance)
  accepted <- list(relative_mcse = .002, relative_change = .004)
  expect_lte(z * accepted[["relative_mcse"]], tolerance)
  boundary <- list(relative_mcse = .0025, relative_change = 0)
  boundary_plan <- list(relative_tolerance = z * boundary[["relative_mcse"]])

  # the reported values stay the raw ones, followed by the rule that rejected them
  raw <- function(error, mcse_text) expect_match(conditionMessage(error), mcse_text, fixed = TRUE)
  with_change <- " (accepted when 2.36 times the standard error and the change are within 'relative_tolerance')"
  error_only  <- " (accepted when 2.36 times the standard error is within 'relative_tolerance')"
  expect_identical(.selection_qmc_criterion(), with_change)
  expect_identical(.selection_qmc_criterion(change = FALSE), error_only)

  factor_judge <- .selection_joint_factor_check_diagnostics
  for (method in c("factor", "cluster")) {
    error <- tryCatch(factor_judge(band, plan, method), error = identity)
    expect_s3_class(error, "error")
    raw(error, paste0("relative Monte Carlo standard error was 0.003 and nested-design relative change was 0.001",
                      with_change, ". Increase 'max_points_per_scramble' or 'scrambles' in "))
    expect_null(factor_judge(accepted, plan, method))
    expect_null(factor_judge(boundary, boundary_plan, method))
    boundary_plan_below <- list(relative_tolerance = boundary_plan[["relative_tolerance"]] * (1 - 1e-9))
    expect_error(factor_judge(boundary, boundary_plan_below, method), "rejected by diagnostics")
    # the nested-design change keeps its own, unscaled, bound
    expect_error(factor_judge(list(relative_mcse = 0, relative_change = .006), plan, method),
                 "rejected by diagnostics")
  }
  # a vector of rows is judged row by row
  expect_error(factor_judge(list(relative_mcse = c(.001, .003), relative_change = c(0, 0)), plan, "factor"),
               "relative Monte Carlo standard error was 0.003")

  error <- tryCatch(.selection_conditioned_sampling_diagnostics(band, plan), error = identity)
  expect_s3_class(error, "error")
  raw(error, paste0("relative Monte Carlo standard error was 0.003 and nested-design relative change was 0.001",
                    with_change, ". Increase 'max_points_per_scramble' or 'scrambles' with "))
  expect_null(.selection_conditioned_sampling_diagnostics(accepted, plan))
  expect_null(.selection_conditioned_sampling_diagnostics(boundary, boundary_plan))
  expect_error(.selection_conditioned_sampling_diagnostics(
    list(relative_mcse = 0, relative_change = .006), plan), "rejected by diagnostics")

  error <- tryCatch(.selection_joint_dense_loglik_check_mcse(mcse, plan), error = identity)
  expect_s3_class(error, "error")
  raw(error, paste0("relative Monte Carlo standard error was 0.003", error_only,
                    ". Increase 'points_per_scramble' or 'scrambles' in "))
  expect_null(.selection_joint_dense_loglik_check_mcse(.002, plan))
  expect_null(.selection_joint_dense_loglik_check_mcse(boundary[["relative_mcse"]], boundary_plan))
  expect_error(.selection_joint_dense_loglik_check_mcse(c(.001, mcse, Inf), plan), "0.003")
})

test_that("the refinement of an event integral keeps going while the scaled error exceeds the tolerance", {

  expect_lte(.selection_qmc_z() * .0021, .005)
  plan <- list(points_per_scramble = 4L, max_points_per_scramble = 8L,
    relative_tolerance = .005, scrambles = 8L, seed = 1L, designs = list())
  calls <- 0L
  # the first budget leaves mcse .003 (z * .003 > .005), the refined one .002 (z * .002 <= .005)
  result <- .selection_joint_checked_event(function(current, rows) {
    calls <<- calls + 1L
    mcse <- if (current[["points_per_scramble"]] == 4L) .003 else .002
    list(log_mass = -.1, relative_mcse = mcse, relative_quadrature_error = 0)
  }, plan, "Test event", "Increase the integration budget.")
  expect_identical(calls, 2L)
  expect_identical(result[["relative_mcse"]], .002)
  # at the budget limit the band is rejected with the raw value
  plan[["points_per_scramble"]] <- 8L
  expect_error(.selection_joint_checked_event(function(current, rows) {
    list(log_mass = -.1, relative_mcse = .003, relative_quadrature_error = 0)
  }, plan, "Test event", "Increase the integration budget."),
    paste0("Test event was rejected by diagnostics: relative integration MCSE was 0.003",
           " (accepted when 2.36 times the standard error is within 'relative_tolerance')",
           ". Increase the integration budget."),
    fixed = TRUE)
  expect_silent(.selection_joint_checked_event(function(current, rows) {
    list(log_mass = -.1, relative_mcse = .0021, relative_quadrature_error = 0)
  }, plan, "Test event", "Increase the integration budget."))
})

test_that("the normalizer-only zplot block judges the scaled error", {

  prior <- BayesTools::prior_weightfunction("one-sided", .025, BayesTools::wf_fixed(c(1, .2)))
  spec <- .selection_spec(list(outcome = list(bias = prior)), c(0, 0), c(1, 1),
    effect_direction = "positive", signed_data = FALSE)
  spec$omega       <- matrix(c(1, .2), 1)
  spec$alpha       <- 0
  spec$phack_kind  <- 0L
  spec$kernel_mode <- SELKERNEL_STEP
  spec$vector_rule <- rep(0L, nrow(spec$omega))
  spec$use_normal  <- FALSE
  mcse <- .003
  expect_gt(.selection_qmc_z() * mcse, .005)
  factors <- list(residual_sd = matrix(1, 1L, 2L), loading = matrix(0, 1L, 2L))
  control <- set_selection_likelihood_control(points_per_scramble = 8L, max_points_per_scramble = 8L)
  block <- function(method, relative_mcse) {
    testthat::local_mocked_bindings(
      .selection_joint_cluster_loglik_block = function(...) {
        list(log_density = 0, relative_mcse = relative_mcse, relative_change = 0, log_normalizer = log(.8))
      },
      .selection_joint_dense_loglik_block = function(...) {
        list(log_density = 0, relative_mcse = relative_mcse, log_normalizer = log(.8))
      },
      .package = "RoBMA"
    )
    tryCatch(.zplot_normalizer_block(
      c(0, 0), matrix(0, 1L, 2L), matrix(c(1, .7, 1), 1L), c(1, 1), spec, control,
      list(block_methods = method, points_per_scramble = 8L, designs = list(`2` = numeric())),
      1L, new.env(), if (method == "rank_one") factors else NULL), error = identity)
  }
  for (method in c("rank_one", "dense")) {
    failure <- block(method, mcse)
    expect_s3_class(failure, "error")
    # the message reports the raw error, as every judge does
    expect_match(conditionMessage(failure),
      paste0("relative integration error was ", format(mcse, digits = 4),
             " (accepted when 2.36 times the standard error and the change are within 'relative_tolerance')",
             ". Increase"),
      fixed = TRUE, info = method)
    accepted <- block(method, .002)
    expect_false(inherits(accepted, "error"), info = method)
    expect_identical(accepted[["log_density"]], -log(.8), info = method)
  }
})

# ---------------------------------------------------------------------------- #
# JAGS wrappers
# ---------------------------------------------------------------------------- #

# The package-owned message of a JAGS node evaluated when the model is created,
# NULL when the node accepts. The parent beta keeps the node from being constant.
.test_qmc_jags <- function(syntax, data) {

  connection <- textConnection(syntax)
  on.exit(close(connection), add = TRUE)
  condition <- tryCatch(rjags::jags.model(
    connection, data = data, n.chains = 1L, n.adapt = 0L, quiet = TRUE,
    inits = list(beta = 0, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 1L)
  ), error = identity)
  if (!inherits(condition, "error")) return(NULL)
  trimws(sub("^RUNTIME ERROR:\n", "", conditionMessage(condition)))
}

.test_qmc_jags_finite <- function(x) ifelse(is.finite(x), x, sign(x) * 1e300)

.test_qmc_jags_step <- function(kind, band, tolerance, points = 256L) {

  block <- band[["case"]][["block"]]
  rank  <- block[["rank"]]
  rules <- if (kind == "factor") .selection_joint_factor_quadrature_rules() else
    .selection_joint_cluster_quadrature_rules(SELNORM_CLUSTER_QUADRATURE_ORDERS)
  qmc <- BayesTools::selection_qmc_design(
    dimensions = 2L * rank, points = points, scrambles = 8L, seed = band[["seed"]])
  data <- list(
    K = block[["rows"]], R = rank, P = points, D = 2L * rank,
    y = rep(0, block[["rows"]]), offset = block[["mean"]], sd = block[["sd"]],
    loading = if (kind == "factor") matrix(as.double(block[["loading"]]), block[["rows"]], rank) else
      as.double(block[["loading"]]),
    sei = block[["sei"]], omega = block[["omega"]],
    z_lower = c(block[["z_lower"]][1L], -1e300), z_upper = c(1e300, block[["z_lower"]][1L]),
    obs_bin = rep(1L, block[["rows"]]),
    nodes = rules[["nodes"]], log_weights = rules[["log_weights"]],
    orders = as.double(rules[["orders"]]), qmc = qmc
  )
  if (kind != "factor") data[["R"]] <- NULL
  syntax <- paste0(
    "model {\n beta ~ dnorm(0,1)\n for (i in 1:K) { mu[i] <- beta + offset[i] }\n",
    "y[1:K] ~ dselnorm_", kind, "_step(mu[1:K],sd[1:K],",
    if (kind == "factor") "loading[1:K,1:R]" else "loading[1:K]",
    ",sei[1:K],omega[1:2],z_lower[1:2],z_upper[1:2],obs_bin[1:K],1,1,1,",
    "nodes,log_weights,orders,qmc[1:8,1:P,1:D],", points, ",", points, ",8,",
    format(tolerance, digits = 17), ",0) }"
  )
  .test_qmc_jags(syntax, data)
}

test_that("the factor and cluster JAGS wrappers judge the scaled error of an exhausted budget", {

  skip_if_not_installed("rjags")
  for (spec in list(list(kind = "factor", name = "crossing_20", seed = 3L),
                    list(kind = "cluster", name = "cluster_189", seed = 2L))) {
    info <- paste(spec[["kind"]], spec[["name"]])
    band <- .test_qmc_band(spec[["name"]], spec[["seed"]])
    expect_lt(band[["low"]], band[["high"]], label = info)
    # the native stage at the same tolerance is the fallback result, so the wrapper
    # judges exactly these diagnostics
    expect_identical(.test_qmc_native(band[["case"]], spec[["seed"]], 256L, 256L,
                                      band[["inside"]])[["relative_mcse"]], band[["mcse"]], info = info)
    rejected <- .test_qmc_jags_step(spec[["kind"]], band, band[["inside"]])
    expect_identical(rejected, paste0(
      "Selection ", spec[["kind"]], " normalizer was rejected by diagnostics: ",
      "relative Monte Carlo standard error was ", sprintf("%g", band[["mcse"]]),
      " and nested-design relative change was ", sprintf("%g", band[["change"]]),
      " (accepted when 2.36 times the standard error and the change are within 'relative_tolerance')",
      ". Increase 'max_points_per_scramble' or 'scrambles' in 'selection_control'."), info = info)
    expect_null(.test_qmc_jags_step(spec[["kind"]], band, band[["high"]] * (1 + 1e-9)), info = info)
  }
})

test_that("the dense covariance JAGS wrapper judges the scaled error of its fixed design", {

  skip_if_not_installed("rjags")
  quadrature <- .selection_joint_cluster_quadrature_rules(c(1L, 3L, 5L))
  points <- 64L
  qmc <- BayesTools::selection_qmc_design(dimensions = 4L, points = points, scrambles = 8L, seed = 1L)
  cutoff <- stats::qnorm(.05, lower.tail = FALSE)
  residual_sd <- c(.2, .25)
  loading <- matrix(c(.8, .8, .3, -.3), 2L)
  covariance <- diag(residual_sd^2) + tcrossprod(loading)
  native <- function(tolerance) .Call("RoBMA_selnorm_mnorm_step_loglik_batch",
    c(.8, 1.2), matrix(c(-.1, .2), 1L), matrix(covariance[lower.tri(covariance, diag = TRUE)], 1L),
    c(.4, .6), matrix(c(1, .1), 1L), c(cutoff, -Inf), c(Inf, cutoff), c(1L, 1L), 1L, TRUE, 1L,
    as.double(qmc), points, 8L, tolerance, FALSE, 0L, quadrature, PACKAGE = "RoBMA")
  mcse <- native(.005)[["relative_mcse"]]
  z <- .selection_qmc_z()
  # fixture: the QMC route, accepted by one standard error at the default tolerance
  expect_gt(mcse, 0)
  expect_lte(mcse, .005)
  expect_gt(z * mcse, .005)
  expect_identical(native(z * mcse)[["relative_mcse"]], mcse)
  syntax <- paste0(
    "model { beta ~ dnorm(0,1)\n for(i in 1:2) { mu[i] <- beta + offset[i] }\n",
    "y[1:2] ~ dselnorm_mnorm_step(mu[],covariance[],sei[],omega[],",
    "z_lower[],z_upper[],bins[],1,1,1,qmc[,,],", points, ",8,tolerance,0,",
    "nodes[],log_weights[],orders[]) }")
  data <- function(tolerance) c(list(
    y = c(.8, 1.2), offset = c(-.1, .2), sei = c(.4, .6), omega = c(1, .1),
    z_lower = c(cutoff, -1e300), z_upper = c(1e300, cutoff), bins = c(1L, 1L),
    covariance = covariance[lower.tri(covariance, diag = TRUE)], qmc = qmc, tolerance = tolerance),
    quadrature[c("nodes", "log_weights", "orders")])
  expect_identical(.test_qmc_jags(syntax, data(.005)), paste0(
    "Selection normalizer was rejected by diagnostics: relative Monte Carlo standard error was ",
    sprintf("%g", mcse), " (accepted when 2.36 times the standard error is within 'relative_tolerance')",
    ". Increase 'points_per_scramble' or 'scrambles' in 'selection_control'."))
  expect_null(.test_qmc_jags(syntax, data(z * mcse * (1 + 1e-9))))
})

# Four rows loading on three factors and their sum: the conditioned-sampling normalizer
# is a randomized estimate (no deterministic rule applies). `initial` and `points` are
# the first and the largest per-scramble budget of its nested design.
.test_qmc_conditioned_arguments <- function(initial, points, tolerance) {

  loading <- rbind(diag(3), c(1, 1, 1))
  cluster_rules <- .selection_joint_cluster_quadrature_rules(SELNORM_CLUSTER_QUADRATURE_ORDERS)
  factor_rules  <- .selection_joint_factor_quadrature_rules()
  list(
    means = matrix(rep(0, 4), 1L), diagonal = matrix(rep(0, 4), 1L),
    loading = matrix(as.double(loading), 1L), rank = 3L, sei = rep(1, 4),
    omega = matrix(c(1, .2), 1L), z_lower = c(0, -Inf), z_upper = c(Inf, 0), sign = 1L,
    kernel_mode = 1L, telescope_probabilities = FALSE, vector_rule = 0L,
    group_index = rep(1L, 4),
    qmc = as.double(BayesTools::selection_qmc_design(
      dimensions = 8L, points = points, scrambles = 8L, seed = 341L)),
    initial_points = as.integer(initial), max_points = as.integer(points), scrambles = 8L,
    relative_tolerance = tolerance,
    cluster_nodes = as.double(cluster_rules$nodes), cluster_log_weights = as.double(cluster_rules$log_weights),
    cluster_orders = as.double(cluster_rules$orders), factor_nodes = as.double(factor_rules$nodes),
    factor_log_weights = as.double(factor_rules$log_weights), factor_orders = as.double(factor_rules$orders)
  )
}

.test_qmc_conditioned_native <- function(arguments) {

  do.call(.Call, c(list("RoBMA_selnorm_conditioned_normalizer_batch"),
                   unname(arguments), list(PACKAGE = "RoBMA")))
}

test_that("the conditioned-sampling JAGS wrapper and loop judge the scaled error", {

  skip_if_not_installed("rjags")
  z <- .selection_qmc_z()
  # a single stage of 64 points per scramble: the wrapper judges its diagnostics
  args <- .test_qmc_conditioned_arguments(64L, 64L, 1)
  diagnostic <- .test_qmc_conditioned_native(args)
  mcse   <- diagnostic[["relative_mcse"]]
  change <- diagnostic[["relative_change"]]
  low    <- max(mcse, change)
  # fixture: a band where one standard error and the change would both pass
  expect_gt(mcse, 0)
  expect_lt(low, z * mcse)
  inside <- sqrt(low * z * mcse)
  args[["relative_tolerance"]] <- inside
  expect_identical(.test_qmc_conditioned_native(args)[["relative_mcse"]], mcse)
  syntax <- paste0(
    "model {
 beta ~ dnorm(0,1)
 for (j in 1:4) {location[j] <- m[j]+beta}
",
    "y[1:4] ~ dselnorm_sampling_conditioned(",
    "location[1:4],C[1:10],d[1:4],L[1:12],3,e0[1:4],u0[1:4],",
    "se[1:4],w[1:2],lo[1:2],hi[1:2],bins[1:4],1,1,0,0,",
    "grp[1:4],qmc[1:Q],64,64,8,tol,cn,cl,co,fn,fl,fo) }")
  data <- function(tolerance) {
    V <- diag(4)
    list(y = rep(0, 4), m = rep(0, 4), C = V[lower.tri(V, diag = TRUE)],
         d = rep(0, 4), L = as.vector(args[["loading"]]), e0 = rep(0, 4), u0 = rep(0, 4),
         se = rep(1, 4), w = c(1, .2), lo = c(0, -1e300), hi = c(1e300, 0),
         bins = rep(1, 4), grp = rep(1, 4), qmc = args[["qmc"]], Q = length(args[["qmc"]]),
         tol = tolerance, cn = args[["cluster_nodes"]], cl = args[["cluster_log_weights"]],
         co = args[["cluster_orders"]], fn = args[["factor_nodes"]], fl = args[["factor_log_weights"]],
         fo = args[["factor_orders"]])
  }
  expect_identical(.test_qmc_jags(syntax, data(inside)), sprintf(paste0(
    "Selection normalizer was rejected by diagnostics: relative Monte Carlo standard error was %g",
    " (accepted when 2.36 times the standard error is within 'relative_tolerance'). ",
    "Increase 'max_points_per_scramble' or 'scrambles' in 'selection_control'."), mcse))
  expect_null(.test_qmc_jags(syntax, data(z * mcse * (1 + 1e-9))))

  # the loop itself, with room to refine from 64 to 256 points: the first stage is
  # accepted at the boundary and left behind inside the band
  long <- .test_qmc_conditioned_arguments(64L, 256L, 1)
  first <- .test_qmc_conditioned_native(long)
  low <- max(first[["relative_mcse"]], first[["relative_change"]])
  high <- z * first[["relative_mcse"]]
  expect_lt(low, high)
  long[["relative_tolerance"]] <- sqrt(low * high)
  continued <- .test_qmc_conditioned_native(long)
  expect_false(identical(continued[["relative_mcse"]], first[["relative_mcse"]]))
  long[["relative_tolerance"]] <- high
  stopped <- .test_qmc_conditioned_native(long)
  expect_identical(stopped[["relative_mcse"]], first[["relative_mcse"]])
  expect_identical(stopped[["log_normalizer"]], first[["log_normalizer"]])
})
