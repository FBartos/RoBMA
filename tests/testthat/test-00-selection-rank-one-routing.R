context("Joint selection rank-one routing")
skip_on_cran()

test_that("transition resolution retains literal signs, nodes and exterior ownership", {

  prior <- BayesTools::prior_weightfunction("one-sided", .025, BayesTools::wf_fixed(c(1, .2)))
  nodes <- .gauss_hermite_nodes(31L)$nodes
  center <- nodes[17L]
  gap <- max(diff(nodes[16:18]))
  for (sign in c(1L, -1L)) for (loading in c(1, -1)) {
    spec <- .selection_spec(list(outcome = list(bias = prior)), 0, .7,
      effect_direction = if (sign == 1L) "positive" else "negative", signed_data = FALSE)
    spec$omega <- matrix(c(1, .2), 1L)
    mean <- sign * spec$z_lower[1L] * .7 - loading * center
    expect_false(.selection_transition_rule_resolves(nodes, mean, gap / 2, matrix(loading), .7, spec))
    expect_true(.selection_transition_rule_resolves(nodes, mean, gap * 1.01, matrix(loading), .7, spec))
    expect_true(.selection_transition_rule_resolves(nodes, mean, .001, matrix(0), .7, spec))
    exterior_mean <- sign * spec$z_lower[1L] * .7 - loading * (max(nodes) + .001)
    expect_true(.selection_transition_rule_resolves(nodes, exterior_mean, .001, matrix(loading), .7, spec))
    spec$omega[,] <- 1
    expect_true(.selection_transition_rule_resolves(nodes, mean, .001, matrix(loading), .7, spec))
  }
})

test_that("a 32-row equal-loading coincident ridge stays within the tolerance", {

  quadrature <- .selection_joint_factor_quadrature_rules()
  qmc <- BayesTools::selection_qmc_design(dimensions = 4L, points = 8192L, scrambles = 8L, seed = 1L)
  # 32 rows with equal loadings 1.2 on both axes and the step 2 SD below their mean.
  K <- 32L
  actual <- .Call("RoBMA_selnorm_factor_step_loglik_batch", rep(1, K), matrix(2, 1L, K),
    matrix(1, 1L, K), matrix(1.2, 1L, 2L * K), rep(1, K), matrix(c(1, 0), 1L),
    c(0, -Inf), c(Inf, 0), rep(1L, K), 1L, TRUE, 1L, quadrature$nodes, quadrature$log_weights,
    as.double(quadrature$orders), as.double(qmc), 512L, 8192L, 8L, .005, TRUE, 0L, PACKAGE = "RoBMA")
  # The rows depend on the latent axes eta only through t(a) %*% eta: a
  # one-dimensional integral over its standardized value is the reference.
  f <- function(u) exp(stats::dnorm(u, log = TRUE) + K * stats::pnorm((2 + sqrt(2) * 1.2 * u) / 1, log.p = TRUE))
  breaks <- c(-Inf, seq(-10, 10, by = .25), Inf)
  reference <- log(sum(vapply(seq_len(length(breaks) - 1L), function(i) stats::integrate(f, breaks[i],
    breaks[i + 1L], rel.tol = 1e-12, abs.tol = 0)$value, numeric(1L))))
  expect_identical(actual$relative_mcse, 0)
  expect_lt(abs(expm1(actual$log_normalizer - reference)), .005)
})

test_that("eight nearly parallel coincident rows keep the normalizer within the tolerance", {

  # Seven rows loading (.224, .975) and one (.224, .976), so the block has full
  # rank; mean -4.08, residual SD .344, selection SE 1, weights (1, 0) at the
  # cutoff 0. A per-row lattice exponent criterion accepted order 31 with relative
  # error 7.9e-3 (Sol, review 33).
  K <- 8L
  loading <- matrix(rep(c(.224, .975), each = K), K, 2L)
  loading[K, 2L] <- .976
  quadrature <- .selection_joint_factor_quadrature_rules()
  qmc <- BayesTools::selection_qmc_design(dimensions = 4L, points = 8192L, scrambles = 8L, seed = 1L)
  actual <- .Call("RoBMA_selnorm_factor_step_loglik_batch", rep(1, K), matrix(-4.08, 1L, K),
    matrix(.344, 1L, K), matrix(as.double(loading), 1L), rep(1, K), matrix(c(1, 0), 1L),
    c(0, -Inf), c(Inf, 0), rep(1L, K), 1L, TRUE, 1L, quadrature$nodes, quadrature$log_weights,
    as.double(quadrature$orders), as.double(qmc), 512L, 8192L, 8L, .005, TRUE, 0L, PACKAGE = "RoBMA")
  # Independent nested adaptive integration over the two latent axes.
  selected <- function(score) stats::pnorm((score - 4.08) / .344)
  inner <- function(u) vapply(u, function(u0) stats::integrate(function(v) {
    score <- .224 * u0 + .975 * v
    stats::dnorm(v) * selected(score)^(K - 1L) * selected(score + .001 * v)
  }, -Inf, Inf, rel.tol = 1e-10, abs.tol = 0, subdivisions = 500L)$value, numeric(1L))
  reference <- log(stats::integrate(function(u) stats::dnorm(u) * inner(u), -Inf, Inf,
    rel.tol = 1e-10, abs.tol = 0, subdivisions = 500L)$value)
  expect_equal(reference, -12.49313108908495, tolerance = 1e-8)
  expect_lt(abs(expm1(actual$log_normalizer - reference)), .005)
})


test_that("the native factor ladder accepts the rung the R certificate ladder predicts", {

  # Acceptance is the shift-amplitude certificate: the first rung that is finite,
  # has both successive changes, passes the tail gate and whose half- and
  # quarter-shift amplitudes and combined allowance fit the tolerance. The R
  # emulation takes its rungs and amplitudes from plain tensor Gauss-Hermite
  # rules, and the accepted value is within the tolerance of an independent
  # integral. A block with no certified rung reaches the randomized fallback.
  quadrature <- .selection_joint_factor_quadrature_rules()
  orders <- quadrature$orders
  qmc <- BayesTools::selection_qmc_design(dimensions = 4L, points = 512L, scrambles = 8L, seed = 1L)
  # an Assink-like pair, three coincident copies, a ridge, a support-one row and
  # an equal sharp pair
  blocks <- list(
    list(mean = c(.2, -.4), s = sqrt(.3) * c(1, 1.2), loading = rbind(c(sqrt(.5), sqrt(.2)), 1.2 * c(sqrt(.5), sqrt(.2))), w2 = .4),
    list(mean = rep(0, 3L), s = rep(sqrt(.3), 3L), loading = matrix(c(sqrt(.5), sqrt(.2)), 3L, 2L, byrow = TRUE), w2 = .2),
    list(mean = rep(1, 6L), s = rep(.6, 6L), loading = matrix(.7, 6L, 2L), w2 = 0),
    list(mean = c(0, -3), s = c(1, 2), loading = rbind(c(1, 0), c(1, .5)), w2 = .3),
    list(mean = c(0, 0), s = c(.1, .1), loading = matrix(1, 2L, 2L), w2 = .5))
  routes <- character()
  for (block in blocks) {
    K <- length(block$mean)
    omega <- matrix(c(1, block$w2), 1L)
    actual <- .Call("RoBMA_selnorm_factor_step_loglik_batch", rep(1, K), matrix(block$mean, 1L),
      matrix(block$s, 1L), matrix(as.double(block$loading), 1L), rep(1, K), omega,
      c(0, -Inf), c(Inf, 0), rep(1L, K), 1L, TRUE, 1L, quadrature$nodes, quadrature$log_weights,
      as.double(orders), as.double(qmc), 512L, 512L, 8L, .005, TRUE, 0L, PACKAGE = "RoBMA")
    certified <- .test_cert_block(lapply(seq_len(K), function(i) {
      list(count = 1L, loading = block$loading[i, ], sd = block$s[i], mean = block$mean[i], sei = 1)
    }), c(1, block$w2))
    expected <- .test_cert_emulate_ladder(certified, .005)
    if (is.null(expected)) {
      routes <- c(routes, "fallback")
      expect_gt(actual$relative_mcse, 0)
    } else {
      routes <- c(routes, "deterministic")
      expect_identical(actual$relative_mcse, 0)
      expect_equal(actual$log_normalizer, expected$value, tolerance = 1e-10)
      # The accepted rung keeps the change diagnostic of the ladder.
      expect_equal(actual$relative_change, expected$change, tolerance = 1e-12)
      expect_lt(.test_cert_relative_error(actual$log_normalizer, .test_cert_reference(certified)), .005)
    }
  }
  expect_true(all(c("deterministic", "fallback") %in% routes))
})


test_that("coupled factor blocks keep their independent integral within the tolerance", {

  weight <- function(mu, sd, cut, omega) omega[2L] + (omega[1L] - omega[2L]) * stats::pnorm((mu - cut) / sd)
  nested <- function(mean, sd, loading, cut, omega) {
    inner <- function(u) vapply(u, function(u0) stats::integrate(function(w) {
      density <- stats::dnorm(w)
      for (i in seq_along(mean)) {
        density <- density * weight(mean[i] + loading[i, 1L] * u0 + loading[i, 2L] * w, sd[i], cut[i], omega)
      }
      density
    }, -Inf, Inf, rel.tol = 1e-12, abs.tol = 0, subdivisions = 1000L)$value, numeric(1L))
    log(stats::integrate(function(u) stats::dnorm(u) * inner(u), -Inf, Inf,
      rel.tol = 1e-11, abs.tol = 0, subdivisions = 1000L)$value)
  }
  quadrature <- .selection_joint_factor_quadrature_rules()
  qmc <- BayesTools::selection_qmc_design(dimensions = 4L, points = 8192L, scrambles = 8L, seed = 1L)
  cutoff <- stats::qnorm(.975)
  evaluate <- function(y, mean, sd, loading, sei, omega) {
    .Call("RoBMA_selnorm_factor_step_loglik_batch", y, matrix(mean, 1L), matrix(sd, 1L),
      matrix(as.double(loading), 1L), sei, matrix(omega, 1L), c(cutoff, -Inf), c(Inf, cutoff),
      ifelse(y / sei >= cutoff, 1L, 2L), 1L, TRUE, 1L, quadrature$nodes, quadrature$log_weights,
      as.double(quadrature$orders), as.double(qmc), 512L, 8192L, 8L, .005, TRUE, 0L, PACKAGE = "RoBMA")
  }

  # An Assink block shape (common and study loadings sqrt(.5) and sqrt(.2), residual
  # sqrt(.3), per standard error): a certified deterministic rung.
  sei <- c(.2, .3, .4)
  mean <- c(.15, .25, .1)
  loading <- cbind(sqrt(.5) * sei, sqrt(.2) * sei)
  actual <- evaluate(c(.5, .1, .9), mean, sqrt(.3) * sei, loading, sei, c(1, .3))
  reference <- nested(mean, sqrt(.3) * sei, loading, cutoff * sei, c(1, .3))
  expect_identical(actual$relative_mcse, 0)
  expect_lt(abs(actual$log_normalizer - reference), 1e-4)

  # V7-01 rows: the geometric rule routed them to the randomized fallback; the
  # shift-amplitude certificate accepts a deterministic rung, within the tolerance.
  sei <- c(.1, .1183)
  for (tau in c(.38, .62)) {
    loading <- cbind(c(.0581, .0688), tau)
    actual <- evaluate(c(.2, .25), c(.08, .16), c(.0814, .0963), loading, sei, c(1, .5))
    reference <- nested(c(.08, .16), c(.0814, .0963), loading, cutoff * sei, c(1, .5))
    expect_identical(actual$relative_mcse, 0)
    expect_lt(abs(actual$log_normalizer - reference), .005)
  }
})

test_that("the vectorized scalar-axis transition rule equals the rule applied row by row", {

  # Dyadic rows put the transition at cut 2 with selection SE .5 exactly on a node
  # (center = node); the others are generic. Rows mix unchanged weights, normal
  # kernels, zero, missing and invalid scales, exterior centers and (through the
  # selection SE and a zero cut) missing centers.
  grid <- c(-2, -1, -.5, 0, .5, 1, 2)
  dyadic <- expand.grid(node = grid, coefficient = c(1, 2, .5, -1), residual_sd = c(.25, 1, 4))
  dyadic$mean <- 1 - dyadic$coefficient * dyadic$node
  set.seed(31)
  n_random <- 160L
  random <- data.frame(
    coefficient = stats::rnorm(n_random) * 2,
    residual_sd = exp(stats::runif(n_random, -3, 1)),
    mean = stats::rnorm(n_random, 0, 2))
  random$coefficient[c(3, 40, 77)] <- c(0, NA, Inf)
  random$residual_sd[c(5, 61, 99, 120)] <- c(0, NA, Inf, -1)
  random$mean[c(7, 88)] <- c(NA, -Inf)
  draws <- rbind(random, dyadic[, names(random)])
  n <- nrow(draws)
  omega <- cbind(1, matrix(sample(c(.2, .5, 1), 2L * n, replace = TRUE), n, 2L))
  omega[seq(1L, n, by = 9L), ] <- 1
  omega[seq(2L, n, by = 11L), 3L] <- omega[seq(2L, n, by = 11L), 2L]
  kernel_mode <- rep(SELKERNEL_STEP, n)
  kernel_mode[seq(4L, n, by = 13L)] <- SELKERNEL_NORMAL
  seen <- logical()
  for (z_lower in list(c(2, .5, 9), c(2, -Inf, 9), c(-Inf, -Inf, 9), c(0, .5, 9))) for (effect_sign in c(1L, -1L)) {
    selection <- list(kernel_mode = kernel_mode, omega = omega, z_lower = z_lower, sign = effect_sign)
    for (sei in c(.5, NA, Inf)) for (nodes in list(.gauss_hermite_nodes(15L)$nodes,
        .gauss_hermite_nodes(127L)$nodes, grid, rev(grid), c(grid, grid[1L]), c(grid, Inf))) {
      rows <- c(seq_len(n), 5L, 5L, n)
      reference <- vapply(rows, function(i) {
        .selection_transition_rule_resolves(nodes, draws$mean[i], draws$residual_sd[i],
          matrix(draws$coefficient[i], 1L), sei,
          list(kernel_mode = kernel_mode[i], omega = omega[i, , drop = FALSE],
               z_lower = z_lower, sign = effect_sign))
      }, logical(1))
      actual <- .selection_scalar_transition_rule_resolves(nodes, draws$mean[rows],
        draws$residual_sd[rows], draws$coefficient[rows], sei, selection, rows)
      expect_identical(actual, reference)
      seen <- c(seen, reference)
    }
  }
  # Both outcomes occur, so the comparison is not vacuous.
  expect_true(any(seen) && any(!seen))
})

test_that("rank-one acceptance returns the rung its certificate ladder accepts", {

  rho <- .99
  K <- 3L
  mean <- rep(-1.5, K)
  sd <- rep(sqrt(1 - rho), K)
  loading <- rep(sqrt(rho), K)
  rules <- .selection_joint_cluster_quadrature_rules(SELNORM_CLUSTER_QUADRATURE_ORDERS)
  spec <- .selection_spec(list(outcome = list(bias = BayesTools::prior_weightfunction(
    "one-sided", .5, BayesTools::wf_fixed(c(1, 0))))), rep(.3, K), rep(1, K),
    effect_direction = "positive", signed_data = FALSE)
  spec$omega <- matrix(c(1, 0), 1L)
  # The geometric transition rule, kept by the zplot mirrors, rejects order 511
  # and resolves order 1023 here; the kernel no longer consults it.
  expect_false(.selection_transition_rule_resolves(.gauss_hermite_nodes(511L)$nodes,
    mean, sd, matrix(loading, K), rep(1, K), spec))
  expect_true(.selection_transition_rule_resolves(.gauss_hermite_nodes(1023L)$nodes,
    mean, sd, matrix(loading, K), rep(1, K), spec))
  actual <- .Call("RoBMA_selnorm_cluster_step_loglik_batch", rep(.3, K), matrix(mean, 1L),
    matrix(sd, 1L), matrix(loading, 1L), rep(1, K), spec$omega, c(0, -Inf), c(Inf, 0),
    rep(1L, K), 1L, TRUE, 1L, rules$nodes, rules$log_weights, as.double(rules$orders),
    as.double(BayesTools::selection_qmc_design(2L, 8192L, 8L, seed = 173L)),
    512L, 8192L, 8L, .005, TRUE, 0L, PACKAGE = "RoBMA")
  reference <- stats::integrate(function(u) stats::dnorm(u) *
    stats::pnorm((-1.5 + sqrt(rho) * u) / sqrt(1 - rho))^K,
    -Inf, Inf, rel.tol = 1e-12, abs.tol = 1e-14)$value
  # The certificate ladder (one change test) accepts the first rung with a
  # certified amplitude; the kernel returns that rung.
  block <- .test_cert_block(list(list(count = K, loading = loading[1L], sd = sd[1L],
    mean = mean[1L], sei = 1)), c(1, 0))
  expected <- .test_cert_emulate_ladder(block, .005,
    orders = SELNORM_CLUSTER_QUADRATURE_ORDERS, changes = 1L)
  expect_false(is.null(expected))
  expect_identical(actual$relative_mcse, 0)
  expect_equal(actual$log_normalizer, expected$value, tolerance = 1e-12)
  expect_lt(abs(expm1(log(reference) - actual$log_normalizer)), .005)
})

test_that("nested intercept random effects use the exact rank-one kernel", {

  dat <- data.frame(
    yi    = c(0.10, 0.20, 0.05, 0.15),
    vi    = rep(0.01, 4L),
    study = factor(c("a", "a", "b", "b")),
    esid  = factor(seq_len(4L))
  )
  object <- bselmodel.mv(
    yi                        = yi,
    vi                        = vi,
    random                    = ~ 1 | study / esid,
    data                      = dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", steps = .025, weights = BayesTools::wf_cumulative(c(1, 1)),
      model = BayesTools::selection_model(
        other_random_effects = "integrate", known_sampling_variance = "integrate", group = "study"
      )
    ),
    only_priors               = TRUE,
    silent                    = TRUE
  )
  plan   <- .data_selection_execution_plan(object[["data"]])
  syntax <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_identical(plan[["exactness"]], "E1")
  expect_identical(
    plan[["random_covariance"]][["representation"]],
    "diagonal_factor"
  )
  expect_match(syntax, "dselnorm_cluster_step", fixed = TRUE)
  expect_false(grepl("dselnorm_mnorm_step", syntax, fixed = TRUE))
  expect_identical(plan[["schema_version"]], 5L)
  expect_identical(plan[["design_keys"]], rep("factor_1", 2L))
  expect_identical(dim(plan[["designs"]][["factor_1"]]), c(8L, 4096L, 2L))
  expect_match(syntax, "sel_joint_qmc_factor_1[1:8,1:4096,1:2],256,4096,8,",
               fixed = TRUE)
  expect_identical(.create_fit_data(object[["data"]], object[["priors"]])[[
    "sel_joint_qmc_factor_1"
  ]], plan[["designs"]][["factor_1"]])
})


test_that("rank-one quadrature rejection reuses controlled factor QMC", {

  # This valid five-row input triggered a GH1023 rejection during Assink
  # adaptation. Independent one-dimensional integration certifies the fallback.
  yi <- c(.2994, .2992, .2989, .291, .217)
  sei <- sqrt(c(.0041, .0042, .0041, .0042, .0041))
  means <- matrix(-.0024773474653312366, 1L, 5L)
  loading <- matrix(1.3700092963713908, 1L, 5L)
  omega <- c(1, .47948435347286433)
  selection <- .selection_spec(
    list(outcome = list(bias = BayesTools::prior_weightfunction(
      "one-sided", .025, BayesTools::wf_fixed(omega)
    ))), yi, sei, effect_direction = "positive", signed_data = FALSE
  )
  selection$omega       <- matrix(omega, 1L)
  selection$vector_rule <- 0L
  plan <- .selection_joint_execution_plan(
    row_blocks = list(1:5), block_methods = "rank_one", factor_ranks = 1L,
    selection_control = set_selection_likelihood_control(), sampling = NULL,
    sampling_factor_blocks = NULL, random_covariance = NULL
  )
  input <- list(yi = yi, means = means, residual_sd = matrix(sei, 1L),
    loading = loading, sei = sei, selection_context = selection,
    execution_plan = plan, return_normalizer = TRUE)
  actual <- do.call(.selection_joint_cluster_loglik_block, input)
  cutoff <- stats::qnorm(.025, lower.tail = FALSE)
  reference <- stats::integrate(function(z) {
    probability <- vapply(seq_along(sei), function(i) {
      p <- stats::pnorm(cutoff * sei[[i]], means[[i]] + loading[[i]] * z,
                         sei[[i]], lower.tail = FALSE)
      omega[[2L]] + (omega[[1L]] - omega[[2L]]) * p
    }, numeric(length(z)))
    exp(stats::dnorm(z, log = TRUE) + rowSums(log(probability)))
  }, -Inf, Inf, rel.tol = 1e-10)
  expect_lt(reference$abs.error / reference$value, 1e-8)
  expect_gt(actual$relative_mcse, 0)
  expect_lte(max(actual$relative_mcse, actual$relative_change),
               plan$relative_tolerance)
  expect_lt(abs(expm1(actual$log_normalizer - log(reference$value))),
              plan$relative_tolerance)
  expected <- mvtnorm::dmvnorm(yi, as.double(means),
    diag(sei^2) + tcrossprod(as.double(loading)), log = TRUE) - log(reference$value)
  expect_equal(actual$log_density, expected, tolerance = .005)

  malformed <- plan
  malformed$designs$factor_1 <- numeric()
  input$execution_plan <- malformed
  condition <- tryCatch(do.call(.selection_joint_cluster_loglik_block, input),
                        error = identity)
  expect_identical(conditionMessage(condition),
    "'qmc' dimensions do not match the cluster integration settings.")
})


test_that("two-coefficient random covariance uses the same rank-one route", {

  dat <- data.frame(
    yi      = c(-0.2, 0.1, 0.3, 0.5),
    vi      = rep(0.04, 4L),
    outcome = factor(c("sensitivity", "specificity",
                       "sensitivity", "specificity")),
    study   = factor(c("s1", "s1", "s2", "s2"))
  )
  object <- bselmodel.mv(
    yi                        = yi,
    vi                        = vi,
    random                    = ~ us(0 + outcome | study),
    prior_heterogeneity       = BayesTools::prior_random(
      study = BayesTools::random_block(
        contrasts = c(outcome = "independent")
      )
    ),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", steps = .025, weights = BayesTools::wf_cumulative(c(1, 1)),
      model = BayesTools::selection_model(
        other_random_effects = "integrate", known_sampling_variance = "integrate", group = "study"
      )
    ),
    only_priors               = TRUE,
    silent                    = TRUE
  )
  plan   <- .data_selection_execution_plan(object[["data"]])
  syntax <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_identical(plan[["exactness"]], "E1")
  expect_identical(
    plan[["random_covariance"]][["representation"]],
    "diagonal_factor"
  )
  expect_match(syntax, "_factor_basis", fixed = TRUE)
  expect_match(syntax, "dselnorm_cluster_step", fixed = TRUE)
  expect_false(grepl("dselnorm_mnorm_step", syntax, fixed = TRUE))
})


test_that("certified factors route by structural rank and dense V stays dense", {

  dat <- data.frame(
    yi    = c(0.10, 0.20, 0.05, 0.15, 0.12, 0.22),
    study = factor(rep(c("a", "b"), each = 3L)),
    esid  = factor(seq_len(6L)),
    time  = rep(1:3, 2L)
  )
  V <- diag(rep(0.01, 6L))
  V[1L, 2L] <- V[2L, 1L] <- 0.002
  V[2L, 3L] <- V[3L, 2L] <- 0.001
  V[4L, 5L] <- V[5L, 4L] <- 0.002
  V[5L, 6L] <- V[6L, 5L] <- 0.001

  correlated <- bselmodel.mv(
    yi                        = yi,
    V                         = V,
    random                    = ~ 1 | study / esid,
    data                      = dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", steps = .025, weights = BayesTools::wf_cumulative(c(1, 1)),
      model = BayesTools::selection_model(
        other_random_effects = "integrate", known_sampling_variance = "integrate", group = "study"
      )
    ),
    only_priors               = TRUE,
    silent                    = TRUE
  )
  correlated_plan   <- .data_selection_execution_plan(correlated[["data"]])
  correlated_syntax <- .create_model_syntax(
    correlated[["data"]],
    correlated[["priors"]]
  )
  expect_identical(correlated_plan[["exactness"]], "E2")
  expect_identical(
    correlated_plan[["random_covariance"]][["representation"]],
    "dense"
  )
  expect_match(correlated_syntax, "dselnorm_mnorm_step", fixed = TRUE)
  expect_false(grepl("dselnorm_cluster_step", correlated_syntax, fixed = TRUE))

  higher_rank <- bselmodel.mv(
    yi                        = yi,
    vi                        = rep(0.01, 6L),
    random                    = ~ har(time | study),
    data                      = dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", steps = .025, weights = BayesTools::wf_cumulative(c(1, 1)),
      model = BayesTools::selection_model(
        other_random_effects = "integrate", known_sampling_variance = "integrate", group = "study"
      )
    ),
    only_priors               = TRUE,
    silent                    = TRUE
  )
  higher_rank_plan   <- .data_selection_execution_plan(higher_rank[["data"]])
  higher_rank_syntax <- .create_model_syntax(
    higher_rank[["data"]],
    higher_rank[["priors"]]
  )
  expect_identical(higher_rank_plan[["exactness"]], "EF")
  expect_identical(
    higher_rank_plan[["random_covariance"]][["representation"]],
    "diagonal_factor"
  )
  expect_identical(
    higher_rank_plan[["factor_ranks"]],
    c(2L, 2L)
  )
  expect_match(higher_rank_syntax, "dselnorm_factor_step", fixed = TRUE)
  expect_false(grepl("dselnorm_mnorm_step", higher_rank_syntax, fixed = TRUE))
  expect_false(grepl("dselnorm_cluster_step", higher_rank_syntax, fixed = TRUE))
})


test_that("crossed random dependencies combine through the generic factor plan", {

  dat <- data.frame(
    yi       = c(-.2, .1, .3, .5),
    vi       = rep(.04, 4L),
    study    = factor(c("s1", "s1", "s2", "s2")),
    observer = factor(c("a", "b", "a", "b"))
  )
  object <- bselmodel.mv(
    yi                        = yi,
    vi                        = vi,
    random                    = list(
      study    = ~ 1 | study,
      observer = ~ 1 | observer
    ),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", steps = .025, weights = BayesTools::wf_cumulative(c(1, 1)),
      model = BayesTools::selection_model(
        other_random_effects = "integrate", known_sampling_variance = "integrate", group = "study"
      )
    ),
    only_priors               = TRUE,
    silent                    = TRUE
  )
  plan <- .data_selection_execution_plan(object[["data"]])
  syntax <- .create_model_syntax(object[["data"]], object[["priors"]])

  expect_identical(plan[["row_blocks"]], list(1:4))
  expect_identical(plan[["exactness"]], "EF")
  expect_identical(
    plan[["factor_ranks"]],
    4L
  )
  expect_match(syntax, "dselnorm_factor_step", fixed = TRUE)
  expect_false(grepl("dselnorm_mnorm_step", syntax, fixed = TRUE))
})


test_that("diagonal block-list V retains generic factor routing", {

  dat <- data.frame(
    yi    = c(-.2, .1, .3, .5),
    study = factor(c("s1", "s1", "s2", "s2"))
  )
  object <- bselmodel.mv(
    yi                        = yi,
    V                         = list(diag(c(.02, .03)), diag(c(.04, .05))),
    random                    = ~ 1 | study,
    data                      = dat,
    measure                   = "SMD",
    prior_unit_information_sd = 1,
    prior_bias = BayesTools::prior_weightfunction(
      "one-sided", steps = .025, weights = BayesTools::wf_cumulative(c(1, 1)),
      model = BayesTools::selection_model(
        other_random_effects = "integrate", known_sampling_variance = "integrate", group = "study"
      )
    ),
    only_priors               = TRUE,
    silent                    = TRUE
  )
  plan <- .data_selection_execution_plan(object[["data"]])

  expect_identical(plan[["exactness"]], "E1")
  expect_identical(
    plan[["block_methods"]],
    rep("rank_one", 2L)
  )
  expect_identical(
    plan[["random_covariance"]][["representation"]],
    "diagonal_factor"
  )
})


test_that("the rank-one Gaussian numerator subtracts no explained quadratic", {

  # Two rows share the unit latent loading 1 and the residual sd 2^-30, so the
  # covariance eigenvalues are d = 2^-60 and 2 + d; the residual (1, 1) lies
  # along the loading. In the normal-kernel mode the density is the numerator.
  # Subtracting the explained quadratic from the total returned 18.61.
  spec <- .selection_spec(
    list(outcome = list(bias = BayesTools::prior_weightfunction(
      "one-sided", .025, BayesTools::wf_fixed(c(1, .5))))),
    c(1, 1), c(1, 1), effect_direction = "positive", signed_data = FALSE
  )
  spec$omega       <- matrix(c(1, .5), 1L)
  spec$kernel_mode <- 0L
  spec$vector_rule <- 0L
  plan <- .selection_joint_execution_plan(
    row_blocks = list(1:2), block_methods = "rank_one", factor_ranks = 1L,
    selection_control = set_selection_likelihood_control(), sampling = NULL,
    sampling_factor_blocks = NULL, random_covariance = NULL
  )
  evaluate <- function(loading) {
    .selection_joint_cluster_loglik_block(
      yi = c(1, 1), means = matrix(0, 1L, 2L),
      residual_sd = matrix(2^-30, 1L, 2L), loading = loading, sei = c(1, 1),
      selection_context = spec, execution_plan = plan, return_normalizer = TRUE
    )
  }
  d <- 2^-60
  reference <- -log(2 * pi) - .5 * (log(d) + log(2 + d) + 2 / (2 + d))
  expect_equal(reference, 18.10996476010904, tolerance = 1e-15)
  observed <- evaluate(matrix(1, 1L, 2L))
  expect_equal(observed$log_density, reference, tolerance = 1e-14)
  expect_identical(as.numeric(observed$log_normalizer), 0)

  # Row-specific loadings against the dense Gaussian density at a moderate scale.
  residual_sd <- c(.3, .5)
  loading <- c(.7, -.4)
  dense <- mvtnorm::dmvnorm(c(1, 1), c(0, 0),
    diag(residual_sd^2) + tcrossprod(loading), log = TRUE)
  moderate <- .selection_joint_cluster_loglik_block(
    yi = c(1, 1), means = matrix(0, 1L, 2L), residual_sd = matrix(residual_sd, 1L),
    loading = matrix(loading, 1L), sei = c(1, 1), selection_context = spec,
    execution_plan = plan, return_normalizer = TRUE
  )
  expect_equal(moderate$log_density, dense, tolerance = 1e-12)
})
