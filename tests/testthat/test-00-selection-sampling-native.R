context("Whole sampling-error selection native kernels")

.sampling_native_arguments <- function(mean, diagonal, loading = matrix(0, length(mean), 0),
                                       omega = c(1, .2), sign = 1L, rule = 0L,
                                       groups = rep(1L, length(mean)), points = 1024L) {

  K <- length(mean)
  cluster_rules <- .selection_joint_cluster_quadrature_rules(SELNORM_CLUSTER_QUADRATURE_ORDERS)
  factor_rules  <- .selection_joint_factor_quadrature_rules()
  list(
    means                  = matrix(as.double(mean), 1L),
    diagonal               = matrix(as.double(diagonal), 1L),
    loading                = matrix(as.double(loading), 1L),
    rank                   = as.integer(ncol(loading)),
    sei                    = rep(1, K),
    omega                  = matrix(as.double(omega), 1L),
    z_lower                = c(0, -Inf),
    z_upper                = c(Inf, 0),
    sign                   = as.integer(sign),
    kernel_mode            = 1L,
    telescope_probabilities = FALSE,
    vector_rule            = as.integer(rule),
    group_index            = as.integer(groups),
    qmc                    = as.double(BayesTools::selection_qmc_design(
      dimensions = max(2L * K, ncol(loading)), points = points,
      scrambles = 8L, seed = 341L
    )),
    initial_points         = as.integer(points),
    max_points             = as.integer(points),
    scrambles              = 8L,
    relative_tolerance     = .005,
    cluster_nodes          = as.double(cluster_rules$nodes),
    cluster_log_weights    = as.double(cluster_rules$log_weights),
    cluster_orders         = as.double(cluster_rules$orders),
    factor_nodes           = as.double(factor_rules$nodes),
    factor_log_weights     = as.double(factor_rules$log_weights),
    factor_orders          = as.double(factor_rules$orders)
  )
}

.sampling_native_mass <- function(arguments) {

  do.call(.Call, c(list("RoBMA_selnorm_conditioned_normalizer_batch"),
                  unname(arguments), list(PACKAGE = "RoBMA")))
}

test_that("analytic sampling selection does not require a QMC grid", {

  args <- .sampling_native_arguments(c(.2, -.1), c(.25, .49))
  args$qmc <- 0.5
  expected <- sum(log(.2 + .8 * pnorm(c(.2, -.1) / c(.5, .7))))
  expect_equal(.sampling_native_mass(args)$log_normalizer, expected,
               tolerance = 1e-14)
  args$qmc <- 0
  expect_error(.sampling_native_mass(args),
    "'qmc' must contain finite values strictly between zero and one.", fixed = TRUE)
  args$qmc <- .5
  args$rank <- 1L
  args$loading <- matrix(c(.3, .1), 1L)
  expect_error(.sampling_native_mass(args),
    "Conditional sampling selection inputs have invalid dimensions.", fixed = TRUE)
})

test_that("selection standard errors scale each row's bin boundaries", {

  # A boundary at z = 0 sits at 0 in outcome units whatever the selection
  # standard error is, so an unequal-SE check at that cut cannot fail. Use a
  # nonzero cut, where the SE genuinely moves the boundary, and pin the value
  # against the analytic independent-row normalizer.
  cut       <- 1.2
  mean      <- c(.2, -.1)
  variance  <- c(.25, .49)
  omega     <- c(1, .2)
  equal     <- c(1, 1)
  unequal   <- c(.5, 1.4)

  analytic <- function(sei) {
    selected <- stats::pnorm(cut * sei, mean, sqrt(variance), lower.tail = FALSE)
    sum(log(omega[1] * selected + omega[2] * (1 - selected)))
  }
  normalizer <- function(sei, z_lower = c(cut, -Inf), z_upper = c(Inf, cut)) {
    args <- .sampling_native_arguments(mean, variance, omega = omega)
    args$sei     <- sei
    args$z_lower <- z_lower
    args$z_upper <- z_upper
    .sampling_native_mass(args)$log_normalizer
  }

  expect_equal(normalizer(equal),   analytic(equal),   tolerance = 1e-14)
  expect_equal(normalizer(unequal), analytic(unequal), tolerance = 1e-14)

  # The comparison only discriminates because the boundary actually moved.
  expect_gt(abs(normalizer(unequal) - normalizer(equal)), .1)

  # It does not move at a zero cut, which is why one is not used above.
  expect_equal(normalizer(unequal, c(0, -Inf), c(Inf, 0)),
               normalizer(equal,   c(0, -Inf), c(Inf, 0)),
               tolerance = 1e-14)
})

test_that("a shared factor keeps each row's own selection standard error", {

  cut        <- .8
  mean       <- c(.2, -.1, .4)
  variance   <- c(.25, .49, .16)
  loading    <- matrix(c(.5, .4, .3), 3L, 1L)
  omega      <- c(1, .3)
  sei        <- c(.4, 1, 1.8)
  covariance <- diag(variance) + tcrossprod(loading)

  bins <- expand.grid(rep(list(seq_along(omega)), length(mean)))
  reference <- log(sum(apply(bins, 1L, function(bin) {
    lower <- ifelse(bin == 1L, cut * sei, -Inf)
    upper <- ifelse(bin == 1L, Inf, cut * sei)
    prod(omega[bin]) * suppressWarnings(as.numeric(mvtnorm::pmvnorm(
      lower     = lower,
      upper     = upper,
      mean      = mean,
      sigma     = covariance,
      algorithm = mvtnorm::Miwa(steps = 4096L)
    )))
  })))

  args <- .sampling_native_arguments(mean, variance, loading, omega = omega)
  args$sei     <- sei
  args$z_lower <- c(cut, -Inf)
  args$z_upper <- c(Inf, cut)
  expect_equal(.sampling_native_mass(args)$log_normalizer, reference,
               tolerance = 1e-6)

  equal <- args
  equal$sei <- rep(1, length(mean))
  expect_gt(abs(.sampling_native_mass(equal)$log_normalizer - reference), .1)
})

test_that("Matheron likelihood preserves singular sampling covariance", {

  V    <- matrix(1, 2, 2)
  d    <- c(.25, .49)
  mean <- c(.1, -.2)
  yi   <- c(.9, -.3)
  e0   <- c(.3, .3)
  u0   <- c(-.1, .4)
  args <- .sampling_native_arguments(mean, d)
  args$qmc <- 0.5
  delta <- as.vector(solve(V + diag(d), yi - mean - e0 - u0))
  e     <- as.vector(e0 + V %*% delta)
  for (sign in c(-1L, 1L)) {
    args$sign <- sign
    bins <- ifelse(sign * yi >= 0, 1L, 2L)
    likelihood_args <- c(list(
      as.double(yi), args$means, as.double(V[lower.tri(V, diag = TRUE)]),
      args$diagonal, args$loading, args$rank,
      matrix(e0, 1L), matrix(u0, 1L), args$sei, args$omega,
      args$z_lower, args$z_upper, as.integer(bins)
    ), unname(args[9:length(args)]))
    out <- do.call(.Call, c(list("RoBMA_selnorm_sampling_conditioned_batch"),
      likelihood_args, list(PACKAGE = "RoBMA")))
    normalizer <- prod(.2 + .8 * pnorm(sign * (mean + e) / sqrt(d)))
    expected <- mvtnorm::dmvnorm(yi, mean, V + diag(d), log = TRUE) +
      sum(log(c(1, .2)[bins])) - log(normalizer)
    expect_equal(as.vector(out$delta), delta, tolerance = 1e-13)
    expect_equal(as.vector(out$e_star), e, tolerance = 1e-13)
    expect_equal(out$log_lik, expected, tolerance = 1e-13)
    expect_equal(out$log_normalizer, log(normalizer), tolerance = 1e-13)
  }
  likelihood_args[[4]]  <- matrix(c(0, .49), 1L)
  likelihood_args[[10]] <- matrix(c(1, 0), 1L)
  expect_error(do.call(.Call, c(list("RoBMA_selnorm_sampling_conditioned_batch"),
    likelihood_args, list(PACKAGE = "RoBMA"))),
    "Conditional sampling selection is unavailable because a fully retained outcome has zero acceptance probability. Use strictly positive selection weights or integrate an outcome-generating source.",
    fixed = TRUE)
})

test_that("whole-sampling likelihoods use relative nonmonotone selection weights", {

  args <- .sampling_native_arguments(c(.1, -.2), c(.25, .49), omega = c(1, .15, 2.5))
  args$z_lower <- c(1.5, -.4, -Inf)
  args$z_upper <- c(Inf, 1.5, -.4)
  call <- c(list(
    c(.4, -.6), args$means, c(1, 1, 1), args$diagonal, args$loading, args$rank,
    matrix(c(.3, .3), 1L), matrix(c(-.1, .4), 1L), args$sei, args$omega,
    args$z_lower, args$z_upper, c(2L, 3L)
  ), unname(args[9:length(args)]))
  evaluate <- function(values) {

    do.call(.Call, c(list("RoBMA_selnorm_sampling_conditioned_batch"), values,
                    list(PACKAGE = "RoBMA")))
  }
  for (rule in c(0L, 1L)) {
    call[[17]] <- rule
    original <- evaluate(call)
    scaled <- call
    scaled[[10]] <- 20 * call[[10]]
    rescaled <- evaluate(scaled)
    expect_equal(rescaled$log_lik, original$log_lik, tolerance = 1e-13)
    expect_equal(rescaled$log_normalizer,
                 original$log_normalizer + if (rule == 0L) 2 * log(20) else log(20),
                 tolerance = 1e-13)
    expect_gt(rescaled$log_normalizer, 0)
    expect_identical(rescaled$e_star, original$e_star)
    expect_identical(rescaled$delta, original$delta)
  }
})

test_that("deterministic and rank-one candidate laws retain exact event probabilities", {

  args <- .sampling_native_arguments(c(.3, -.4, .2), c(0, 0, 0))
  expect_equal(.sampling_native_mass(args)$log_normalizer, log(.2), tolerance = 1e-14)
  args$omega <- matrix(c(1, 0), 1L)
  expect_identical(.sampling_native_mass(args)$log_normalizer, -Inf)
  args$omega <- matrix(c(1, .2), 1L)
  loading <- c(.5, 1, -.2)
  args$loading <- matrix(loading, 1L)
  args$rank <- 1L
  # Independent scalar partition at the three threshold crossings.
  cuts <- sort(c(-Inf, -as.vector(args$means) / loading, Inf))
  mids <- c(cuts[2] - 1, (cuts[2] + cuts[3]) / 2,
            (cuts[3] + cuts[4]) / 2, cuts[4] + 1)
  weights <- vapply(mids, function(z) {

    prod(ifelse(as.vector(args$means) + loading * z >= 0, 1, .2))
  }, numeric(1))
  reference <- sum(diff(pnorm(cuts)) * weights)
  out <- .sampling_native_mass(args)
  expect_equal(exp(out$log_normalizer), reference, tolerance = 1e-13)
  expect_identical(out$relative_mcse, 0)
  expect_identical(out$relative_change, 0)
})

test_that("shared Gaussian events evaluate deterministic laws and event bounds exactly", {

  args <- .sampling_native_arguments(c(-1, .5), c(0, 0), points = 8L)
  event_args <- list(args$means, matrix(0, 1L, 3L), args$sei, args$omega,
    args$z_lower, args$z_upper, args$sign, args$kernel_mode, args$vector_rule,
    matrix(c(-2, .5), 1L), matrix(c(0, 1), 1L), args$qmc, 8L, 8L, NULL, NULL, .005)
  evaluate <- function(values) {

    do.call(.Call, c(list("RoBMA_selnorm_gaussian_event_mass_batch"), values,
                    list(PACKAGE = "RoBMA")))
  }
  inside <- evaluate(event_args)
  expect_equal(inside$log_mass, log(.2), tolerance = 1e-14)
  expect_identical(inside$relative_mcse, 0)
  event_args[[11]] <- matrix(c(0, .4), 1L)
  outside <- evaluate(event_args)
  expect_identical(outside$log_mass, -Inf)
  expect_identical(outside$relative_mcse, 0)
})

test_that("general singular factor normalizers agree with Gaussian orthant identities", {

  loading <- rbind(c(1, 0), c(0, 1), c(1, 1))
  args <- .sampling_native_arguments(rep(0, 3), rep(0, 3), loading,
                                    points = 4096L)
  out <- .sampling_native_mass(args)
  # P(Y_i > 0)=1/2; sum pair probabilities is 1; trivariate is 1/4.
  # The pair and trivariate formulas follow the centered Gaussian arcsine law.
  w <- .2
  reference <- w^3 + 1.5 * w^2 * (1 - w) + w * (1 - w)^2 + .25 * (1 - w)^3
  expect_equal(exp(out$log_normalizer), reference, tolerance = 5e-4)
  expect_lt(out$relative_mcse, .005)
  expect_lt(out$relative_change, .005)
})

test_that("rank-two factor quadrature partitions deterministic-row jumps", {

  mean <- c(.11, -.25, .33)
  loading <- .3 * rbind(c(1, -1), c(1, 0), c(1, 1))
  args <- .sampling_native_arguments(mean, rep(0, 3), loading, points = 8L)
  args$relative_tolerance <- 1e-6
  integrand <- function(t) {

    lo <- -(mean[3] + .3 * t) / .3
    hi <- (mean[1] + .3 * t) / .3
    joint <- ifelse(lo < hi, pnorm(hi) - pnorm(lo), 0)
    pair_weight <- .2^2 + .2 * .8 * (pnorm(hi) + pnorm(lo, lower.tail = FALSE)) + .8^2 * joint
    dnorm(t) * ifelse(mean[2] + .3 * t >= 0, 1, .2) * pair_weight
  }
  cuts <- c(-Inf, -(mean[1] + mean[3]) / .6, -mean[2] / .3, Inf)
  reference <- sum(vapply(seq_len(3), function(i) {

    integrate(integrand, cuts[i], cuts[i + 1], rel.tol = 1e-10)$value
  }, numeric(1)))
  out <- .sampling_native_mass(args)
  expect_equal(exp(out$log_normalizer), reference, tolerance = 1e-7)
  expect_identical(out$relative_mcse, 0)
  expect_lt(out$relative_change, args$relative_tolerance)
})

test_that("PSD projections preserve exact CDF and tail endpoints and empty selection events", {

  args <- .sampling_native_arguments(rep(0, 3), rep(0, 3),
    rbind(c(1, 0), c(0, 1), c(1, 1)), omega = c(1, 0), points = 64L)
  project <- function(z, kind) {

    do.call(.Call, c(list("RoBMA_selnorm_factor_projection_batch"),
      unname(args), list(as.double(z), as.integer(kind)), list(PACKAGE = "RoBMA")))
  }
  cdf <- project(c(-Inf, -1, 0, Inf), 1L)
  expect_equal(as.vector(cdf$density), c(0, 0, 0, 1), tolerance = 0)
  expect_lt(cdf$relative_mcse, args$relative_tolerance)
  tail <- project(0, 2L)
  expect_equal(as.vector(tail$density), 1, tolerance = 0)
  expect_lt(tail$relative_mcse, args$relative_tolerance)
})

test_that("positive-residual small factors reuse accurate deterministic quadrature", {

  mean    <- c(.1, -.2, .4)
  d       <- c(.2, .3, .4)
  loading <- c(1, .8, 1.2)
  args <- .sampling_native_arguments(mean, d, matrix(loading, ncol = 1), points = 8L)
  args$relative_tolerance <- 1e-8
  out <- .sampling_native_mass(args)
  reference <- integrate(function(z) {

    vapply(z, function(value) {

      dnorm(value) * prod(.2 + .8 * pnorm((mean + loading * value) / sqrt(d)))
    }, numeric(1))
  }, -Inf, Inf, rel.tol = 1e-11)$value
  expect_equal(exp(out$log_normalizer), reference, tolerance = 1e-8)
  expect_identical(out$relative_mcse, 0)
  expect_lt(out$relative_change, args$relative_tolerance)

  loading <- rbind(c(1, .2), c(.3, .8), c(.5, -.4))
  args <- .sampling_native_arguments(rep(0, 3), d, loading, points = 8L)
  # Preserve the existing two-successive-change factor criterion at the
  # integration fixture's 1e-6 tolerance; check the result more tightly below.
  args$relative_tolerance <- 1e-6
  out <- .sampling_native_mass(args)
  correlation <- cov2cor(diag(d) + tcrossprod(loading))
  arcsines <- sum(asin(correlation[lower.tri(correlation)]))
  sum_pairs <- 3 / 4 + arcsines / (2 * pi)
  triple    <- 1 / 8 + arcsines / (4 * pi)
  reference <- .2^3 + 1.5 * .2^2 * .8 + .2 * .8^2 * sum_pairs + .8^3 * triple
  expect_equal(exp(out$log_normalizer), reference, tolerance = 1e-8)
  expect_identical(out$relative_mcse, 0)
  expect_lt(out$relative_change, args$relative_tolerance)
})

test_that("best rules preserve publication groups when only sampling connects them", {

  args <- .sampling_native_arguments(c(.1, -.4, .6), c(.3, .7, 0),
                                    rule = 1L, groups = c(1, 1, 2))
  expected <- (1 - .8 * prod(pnorm(-c(.1, -.4) / sqrt(c(.3, .7))))) * 1
  expect_equal(exp(.sampling_native_mass(args)$log_normalizer), expected,
               tolerance = 1e-13)
})

test_that("fully conditioned fits cancel positive weights and reject impossible contexts", {

  args <- .sampling_native_arguments(.4, 0, omega = c(1, 1e-300))
  call <- c(list(.7, args$means, 1, args$diagonal, args$loading, 0L,
                 matrix(.2, 1L), matrix(0, 1L), args$sei, args$omega,
                 args$z_lower, args$z_upper, 1L), unname(args[9:length(args)]))
  out <- do.call(.Call, c(list("RoBMA_selnorm_sampling_conditioned_batch"),
                        call, list(PACKAGE = "RoBMA")))
  expect_equal(out$log_lik, dnorm(.7, .4, 1, log = TRUE), tolerance = 1e-14)
  expect_equal(as.vector(out$e_star), .3, tolerance = 1e-14)
  # Retained covariance can vary by posterior row when it includes random
  # effects. The source correction remains defined when weights cancel.
  varying <- call
  varying[[2]] <- matrix(.4, 2L, 1L)
  varying[[3]] <- matrix(c(1, 4), 2L, 1L)
  varying[[4]] <- matrix(0, 2L, 1L)
  varying[[5]] <- matrix(0, 2L, 0L)
  varying[[7]] <- matrix(.2, 2L, 1L)
  varying[[8]] <- matrix(0, 2L, 1L)
  varying[[10]] <- matrix(rep(c(1, 1e-300), each = 2L), 2L)
  varied <- do.call(.Call, c(list("RoBMA_selnorm_sampling_conditioned_batch"),
                           varying, list(PACKAGE = "RoBMA")))
  expect_equal(varied$log_lik, dnorm(.7, .4, c(1, 2), log = TRUE), tolerance = 1e-14)
  expect_equal(as.vector(varied$delta), .1 / c(1, 4), tolerance = 1e-14)
  call[[10]] <- matrix(c(1, 0), 1L)
  expect_error(do.call(.Call, c(list("RoBMA_selnorm_sampling_conditioned_batch"),
    call, list(PACKAGE = "RoBMA"))),
    "Conditional sampling selection is unavailable because a fully retained outcome has zero acceptance probability. Use strictly positive selection weights or integrate an outcome-generating source.",
    fixed = TRUE)
})

test_that("deterministic selected response generation does not retry positive weights", {

  inputs <- list(matrix(-.2, 1L), array(0, c(1L, 1L, 1L)), 1,
                 matrix(c(1, 1e-300), 1L), c(0, -Inf), c(Inf, 0),
                 1L, 1L, list(1L), 1L, 0L)
  out <- do.call(.Call, c(list("RoBMA_selnorm_mnorm_step_rng_batch"), inputs,
                        list(PACKAGE = "RoBMA")))
  expect_identical(out$draws, matrix(-.2, 1L))
  expect_identical(out$failure_code, 0L)
  inputs[[4]] <- matrix(c(1, 0), 1L)
  expect_error(do.call(.Call, c(list("RoBMA_selnorm_mnorm_step_rng_batch"), inputs,
    list(PACKAGE = "RoBMA"))),
    "Selected response simulation is unavailable because a fully retained outcome has zero acceptance probability. Use strictly positive selection weights or integrate an outcome-generating source.",
    fixed = TRUE)
})

test_that("JAGS reports the failed sampling-normalizer diagnostic and applicable remedy", {

  skip_on_cran()
  skip_if_not_installed("rjags")
  loading <- rbind(diag(3), c(1, 1, 1))
  args <- .sampling_native_arguments(rep(0, 4), rep(0, 4), loading, points = 8L)
  args$relative_tolerance <- 1e-12
  syntax <- paste0(
    "model {\n beta ~ dnorm(0,1)\n for (j in 1:4) {location[j] <- m[j]+beta}\n",
    "y[1:4] ~ dselnorm_sampling_conditioned(",
    "location[1:4],C[1:10],d[1:4],L[1:12],3,e0[1:4],u0[1:4],",
    "se[1:4],w[1:2],lo[1:2],hi[1:2],bins[1:4],1,1,0,0,",
    "grp[1:4],qmc[1:Q],8,8,8,tol,cn,cl,co,fn,fl,fo) }"
  )
  for (criterion in c("mcse", "change")) {
    if (criterion == "change") {
      design <- array(args$qmc, c(8L, 8L, 8L))
      for (scramble in 2:8) design[scramble, , ] <- design[1, , ]
      args$qmc <- as.double(design)
    }
    diagnostic <- .sampling_native_mass(args)
    if (criterion == "mcse") {
      expect_gt(diagnostic$relative_mcse, args$relative_tolerance)
      expected <- sprintf(paste0(
        "Selection normalizer was rejected by diagnostics: relative Monte Carlo standard error was %.6f. ",
        "Increase 'max_points_per_scramble' or 'scrambles' in 'selection_control'."
      ), diagnostic$relative_mcse)
    } else {
      expect_lte(diagnostic$relative_mcse, args$relative_tolerance)
      expect_gt(diagnostic$relative_change, args$relative_tolerance)
      expected <- sprintf(paste0(
        "Selection normalizer was rejected by diagnostics: relative nested-design change was %.6f. ",
        "Increase 'max_points_per_scramble' in 'selection_control'."
      ), diagnostic$relative_change)
    }
    V <- diag(4)
    data <- list(y = rep(0, 4), m = rep(0, 4), C = V[lower.tri(V, diag = TRUE)],
                 d = rep(0, 4), L = as.vector(loading), e0 = rep(0, 4), u0 = rep(0, 4),
                 se = rep(1, 4), w = c(1, .2), lo = c(0, -1e300), hi = c(1e300, 0),
                 bins = rep(1, 4), grp = rep(1, 4), qmc = args$qmc, Q = length(args$qmc),
                 tol = args$relative_tolerance,
                 cn = args$cluster_nodes, cl = args$cluster_log_weights, co = args$cluster_orders,
                 fn = args$factor_nodes, fl = args$factor_log_weights,
                 fo = args$factor_orders)
    connection <- textConnection(syntax)
    expect_error({
      model <- rjags::jags.model(connection, data = data, inits = list(beta = 0),
                                n.chains = 1L, n.adapt = 0L, quiet = TRUE)
      stats::update(model, 1L, progress.bar = "none")
    }, expected, fixed = TRUE)
    close(connection)
  }
})
