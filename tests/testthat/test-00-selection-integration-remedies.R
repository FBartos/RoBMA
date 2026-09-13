test_that("integration rejections name fitting and supported post-fit remedies", {

  selection <- .selection_spec(
    list(outcome = list(bias = BayesTools::prior_weightfunction(
      "one-sided", .025, BayesTools::wf_fixed(c(1, .5))
    ))),
    c(0, 0), c(1, 1), effect_direction = "positive", signed_data = FALSE
  )
  control <- set_selection_likelihood_control(
    points_per_scramble = 8L, max_points_per_scramble = 8L, scrambles = 2L
  )
  plan <- .selection_joint_execution_plan(
    row_blocks = list(1:2), block_methods = "factor", factor_ranks = 2L,
    selection_control = control, sampling = NULL,
    sampling_factor_blocks = NULL, random_covariance = NULL
  )
  common <- list(
    yi = c(0, 0), means = matrix(0, 1L, 2L), sei = c(1, 1),
    selection_context = selection, execution_plan = plan
  )
  native_result <- NULL
  native_environment <- new.env(parent = environment(.selection_joint_dense_loglik_block))
  # Inject diagnostics at the native boundary; no numerical calculation is replaced
  # in production, and these tests certify message formatting, not integration.
  native_environment[[".Call"]] <- function(...) native_result
  dense <- .selection_joint_dense_loglik_block
  factor <- .selection_joint_factor_loglik_block
  cluster <- .selection_joint_cluster_loglik_block
  environment(dense) <- native_environment
  environment(factor) <- native_environment
  environment(cluster) <- native_environment
  cluster_plan <- .selection_joint_execution_plan(
    row_blocks = list(1:2), block_methods = "rank_one", factor_ranks = 1L,
    selection_control = control, sampling = NULL,
    sampling_factor_blocks = NULL, random_covariance = NULL
  )
  remedies <- paste0(
    "'set_selection_likelihood_control()'. Pass the control as ",
    "'selection_control' when fitting, 'integration_control' in zplot(), or ",
    "'density_control$integration_control' for supported post-fit ",
    "densities and hypotheses."
  )

  for (mcse in c(.006, Inf)) {
    native_result <- list(
      log_density = 0, relative_mcse = mcse, log_normalizer = 0,
      integration_diagnostics = matrix(0, 1L, 4L)
    )
    condition <- tryCatch(
      do.call(dense, c(common, list(
        covariance_lower = matrix(c(1, 0, 1), 1L), block_size = 2L
      ))), error = identity
    )
    expect_identical(conditionMessage(condition), paste0(
      "Selection normalizer was rejected by diagnostics: relative ",
      "Monte Carlo standard error was ", if (is.finite(mcse)) "0.006" else "Inf",
      ". Increase 'points_per_scramble' in ", remedies
    ))
    expect_null(conditionCall(condition))

    native_result <- list(
      log_density = 0, relative_mcse = mcse, relative_change = .007,
      log_normalizer = 0
    )
    condition <- tryCatch(
      do.call(factor, c(common, list(
        residual_sd = matrix(1, 1L, 2L), loading = matrix(0, 1L, 4L),
        block_index = 1L
      ))), error = identity
    )
    expect_identical(conditionMessage(condition), paste0(
      "Selection factor normalizer was rejected by diagnostics: ",
      "relative Monte Carlo standard error was ",
      if (is.finite(mcse)) "0.006" else "Inf",
      " and nested-design relative change was 0.007. ",
      "Increase 'max_points_per_scramble' or 'scrambles' in ", remedies
    ))
    expect_null(conditionCall(condition))

    cluster_common <- common
    cluster_common$execution_plan <- cluster_plan
    condition <- tryCatch(
      do.call(cluster, c(cluster_common, list(
        residual_sd = matrix(1, 1L, 2L), loading = matrix(0, 1L, 2L)
      ))), error = identity
    )
    expect_identical(conditionMessage(condition), paste0(
      "Selection cluster normalizer was rejected by diagnostics: ",
      "relative Monte Carlo standard error was ",
      if (is.finite(mcse)) "0.006" else "Inf",
      " and nested-design relative change was 0.007. ",
      "Increase 'max_points_per_scramble' or 'scrambles' in ", remedies
    ))
    expect_null(conditionCall(condition))
  }
})


test_that("JAGS selection rejections report complete model-neutral diagnostics", {

  skip_if_not_installed("rjags")
  # These deliberately small integration designs exercise message formatting.
  # The raw R-native result supplies the reported metric, not a numerical oracle.
  quadrature <- .selection_joint_cluster_quadrature_rules(c(1L, 3L, 5L))
  qmc <- BayesTools::selection_qmc_design(
    dimensions = 4L, points = 8L, scrambles = 2L, seed = 1L
  )
  cluster_qmc <- BayesTools::selection_qmc_design(
    dimensions = 2L, points = 8L, scrambles = 2L, seed = 1L
  )
  cutoff <- stats::qnorm(.05, lower.tail = FALSE)
  data <- list(
    y = c(.8, 1.2), offset = c(-.1, .2), sei = c(.4, .6),
    omega = c(1, .1), z_lower = c(cutoff, -1e300),
    z_upper = c(1e300, cutoff), bins = c(1L, 1L)
  )
  residual_sd <- c(.2, .25)
  loading <- matrix(c(.8, .8, .3, -.3), 2L)
  covariance <- diag(residual_sd^2) + tcrossprod(loading)
  common <- list(data$sei, matrix(data$omega, 1L), c(cutoff, -Inf),
                   c(Inf, cutoff), data$bins, 1L, TRUE, SELKERNEL_STEP)
  rules <- list(quadrature$nodes, quadrature$log_weights,
                  as.double(quadrature$orders))
  inputs <- list(
    mnorm = c(list(data$y, matrix(data$offset, 1L),
                    matrix(covariance[lower.tri(covariance, diag = TRUE)], 1L)),
                common, list(as.double(qmc), 8L, 2L, .005, FALSE)),
    cluster = c(list(data$y, matrix(data$offset, 1L), matrix(residual_sd, 1L),
                      matrix(loading[, 1L], 1L)), common, rules,
                  list(as.double(cluster_qmc), 8L, 8L, 2L, .005, FALSE)),
    factor = c(list(data$y, matrix(data$offset, 1L), matrix(residual_sd, 1L),
                     matrix(as.double(loading), 1L)), common, rules,
                 list(as.double(length(quadrature$orders)), as.double(qmc),
                        8L, 8L, 2L, .005, FALSE))
  )
  for (method in names(inputs)) {
    tail_arguments <- if (method == "mnorm") list(0L, quadrature) else list(0L)
    raw <- do.call(.Call, c(list(
      paste0("RoBMA_selnorm_", method, "_step_loglik_batch")
    ), inputs[[method]], tail_arguments, list(PACKAGE = "RoBMA")))
    metrics <- unlist(raw[intersect(c("relative_mcse", "relative_change"), names(raw))])
    expect_true(all(is.finite(metrics)), info = method)
    expect_gt(max(metrics), .005)
    message <- switch(method,
      mnorm = paste0(
        "Selection normalizer was rejected by diagnostics: relative Monte Carlo ",
        "standard error was ", sprintf("%.6f", raw$relative_mcse),
        ". Increase 'points_per_scramble' or 'scrambles' in 'selection_control'."
      ),
      cluster = paste0(
        "Selection cluster normalizer was rejected by diagnostics: relative Monte Carlo ",
        "standard error was ", sprintf("%.6f", raw$relative_mcse),
        " and nested-design relative change was ", sprintf("%.6f", raw$relative_change),
        ". Increase 'max_points_per_scramble' or 'scrambles' in 'selection_control'."
      ),
      factor = paste0(
        "Selection factor normalizer was rejected by diagnostics: relative Monte Carlo ",
        "standard error was ", sprintf("%.6f", raw$relative_mcse),
        " and nested-design relative change was ", sprintf("%.6f", raw$relative_change),
        ". Increase 'max_points_per_scramble' or 'scrambles' in 'selection_control'."
      )
    )
    current_data <- data
    if (method == "mnorm") {
      current_data$covariance <- covariance[lower.tri(covariance, diag = TRUE)]
      current_data <- c(current_data, quadrature[c("nodes", "log_weights", "orders")])
    } else {
      current_data$residual_sd <- residual_sd
      current_data$loading <- if (method == "cluster") loading[, 1L] else loading
      current_data <- c(current_data, quadrature[c("nodes", "log_weights", "orders")])
    }
    current_data$qmc <- if (method == "cluster") cluster_qmc else qmc
    if (method == "factor") current_data$rule_counts <- length(quadrature$orders)
    covariance_arguments <- switch(method,
      mnorm = "mu[],covariance[],", cluster = "mu[],residual_sd[],loading[],",
      factor = "mu[],residual_sd[],loading[,],"
    )
    integration_arguments <- switch(method,
      mnorm = "qmc[,,],8,2,.005,0,nodes[],log_weights[],orders[]",
      cluster = "nodes[],log_weights[],orders[],qmc[,,],8,8,2,.005,0",
      factor = "nodes[],log_weights[],orders[],rule_counts,qmc[,,],8,8,2,.005,0"
    )
    syntax <- paste0(
      "model { beta ~ dnorm(0,1)\n",
      "for(i in 1:2) { mu[i] <- beta + offset[i] }\n",
      "y[1:2] ~ dselnorm_", method, "_step(", covariance_arguments,
      "sei[],omega[],z_lower[],z_upper[],bins[],1,1,1,", integration_arguments, ") }"
    )
    connection <- textConnection(syntax)
    condition <- tryCatch(rjags::jags.model(
      connection, data = current_data, n.chains = 1L, n.adapt = 0L, quiet = TRUE,
      inits = list(beta = 0, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 1L)
    ), error = identity, finally = close(connection))
    expect_s3_class(condition, "error")
    # JAGS owns its runtime-error header and trailing line breaks. Compare the
    # complete package-owned diagnostic, including its observed value and action.
    expect_true(startsWith(conditionMessage(condition), "RUNTIME ERROR:\n"), info = method)
    expect_identical(trimws(sub("^RUNTIME ERROR:\n", "", conditionMessage(condition))),
                       message, info = method)
  }
})


test_that("invalid zplot factor supports have a complete structural message", {

  for (support in list(matrix(1, 2L, 2L), matrix(NA, 2L, 2L), matrix(TRUE, 2L, 1L))) {
    condition <- tryCatch(.zplot_joint_block(
      z = 0, mean = matrix(0, 1L, 2L), covariance_lower = matrix(c(1, 0, 1), 1L),
      sei = c(1, 1), selection = NULL, probability = FALSE,
      control = set_selection_likelihood_control(),
      factors = list(residual_sd = matrix(1, 1L, 2L),
                       loading = matrix(0, 1L, 4L), loading_support = support)
    ), error = identity)
    expect_identical(conditionMessage(condition),
                       "Selection factor loading supports are invalid.")
    expect_null(conditionCall(condition))
  }
})
