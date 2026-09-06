# Small, independently generated selection-model certification examples.
# No package selected-normal RNG or likelihood is used by the simulator/oracle.



.selection_recovery_design <- function(shape, studies = 40L) {

  shape <- match.arg(shape, c("diagonal", "structured", "dense"))
  sei   <- c(.16, .24, .32)
  dat   <- data.frame(
    study = rep(seq_len(studies), each = 3L),
    esid  = seq_len(3L * studies),
    obs   = rep(seq_len(3L), studies),
    type  = rep("outcome", 3L * studies),
    vi    = rep(sei^2, studies)
  )
  correlation <- switch(
    shape,
    diagonal   = diag(3L),
    structured = .35 + .65 * diag(3L),
    dense      = matrix(c(1, .25, -.10, .25, 1, .15, -.10, .15, 1), 3L)
  )
  block <- correlation * tcrossprod(sei)
  V <- switch(
    shape,
    diagonal = diag(dat$vi),
    structured = vcalc2(vi, cluster = study, type = type, obs = obs,
                        rho = c(.35, .20), data = dat),
    dense = kronecker(diag(studies), block)
  )
  # Independently specified conditioning decomposition: the structured case
  # has one common sampling factor; the ordinary dense matrix retains its
  # requested residual fraction. No package decomposition is used here.
  residual <- switch(shape, diagonal = diag(block),
                       structured = .65 * diag(block), dense = .1 * diag(block))
  loading <- switch(shape, diagonal = matrix(numeric(), 3L, 0L),
                      structured = matrix(sqrt(.35) * sei, ncol = 1L),
                      dense = t(chol(block - diag(residual))))
  latent <- tcrossprod(loading)
  stopifnot(min(eigen(correlation, symmetric = TRUE, only.values = TRUE)$values) > .1)
  stopifnot(isTRUE(all.equal(diag(residual) + latent, block, tolerance = 1e-14)))
  list(data = dat, V = V, block = block, residual = residual, latent = latent,
         loading = loading)
}


.selection_recovery_simulate <- function(design, target, truth, seed) {

  set.seed(seed)
  dat       <- design$data
  block     <- design$block
  threshold <- qnorm(.05, lower.tail = FALSE) * sqrt(diag(block))
  direction <- if (truth[["mu"]] < 0) -1 else 1
  weight    <- truth[["omega"]]
  studies   <- unique(dat$study)
  yi        <- numeric(nrow(dat))
  proposals <- 0L
  for (study in studies) {
    rows <- which(dat$study == study)
    if (target == "exact") {
      # Redraw ALL Gaussian effects when a complete candidate block is rejected.
      covariance <- block + truth[["study"]]^2 + diag(truth[["effect"]]^2, 3L)
      root       <- chol(covariance)
      repeat {
        candidate <- truth[["mu"]] + as.numeric(rnorm(3L) %*% root)
        accept    <- ifelse(direction * candidate > threshold, 1, weight)
        proposals <- proposals + 1L
        if (runif(1L) < prod(accept)) break
        if (proposals > 1000000L) stop("Recovery simulator exhausted its proposal budget.")
      }
      yi[rows] <- candidate
    } else {
      # Draw the conditioning effects once, then reject each row separately.
      # Estimate-specific heterogeneity is marginalized in the fitted target.
      sampling <- as.numeric(design$loading %*% rnorm(ncol(design$loading)))
      location <- truth[["mu"]] + rnorm(1L, sd = truth[["study"]]) + sampling
      sd       <- sqrt(design$residual + truth[["effect"]]^2)
      for (j in seq_len(3L)) {
        repeat {
          candidate <- rnorm(1L, location[[j]], sd[[j]])
          accept    <- if (direction * candidate > threshold[[j]]) 1 else weight
          proposals <- proposals + 1L
          if (runif(1L) < accept) break
          if (proposals > 1000000L) stop("Recovery simulator exhausted its proposal budget.")
        }
        yi[rows[[j]]] <- candidate
      }
    }
  }
  dat$yi <- yi
  list(data = dat, proposals = proposals)
}


.selection_recovery_fit <- function(design, dat, target, random, seed,
                                    sample = 6000L, fixed = NULL,
                                    effect_direction = "positive",
                                    selection_control = set_selection_likelihood_control()) {

  sd_prior <- BayesTools::prior("normal", list(0, .5), truncation = list(0, Inf))
  args <- list(
    yi                        = dat$yi,
    V                         = design$V,
    data                      = dat,
    measure                   = "GEN",
    prior_effect              = BayesTools::prior("normal", list(0, 1)),
    prior_unit_information_sd = 1,
    prior_bias                = BayesTools::prior_weightfunction(
      side    = "one-sided",
      steps   = .05,
      weights = if (is.null(fixed)) {
        BayesTools::wf_independent(BayesTools::prior("beta", list(1, 1)))
      } else BayesTools::wf_fixed(c(1, fixed[["omega"]]))
    ),
    effect_direction          = effect_direction,
    selection_likelihood      = target,
    selection_control         = selection_control,
    marginalize_estimate_level = TRUE,
    chains                    = 3L,
    parallel                  = TRUE,
    sample                    = sample,
    burnin                    = 1500L,
    adapt                     = 500L,
    seed                      = seed,
    silent                    = TRUE,
    # Diagnose the selected semantic parameters explicitly below, using rank
    # normalized R-hat and MCSE; never extend or retry until a check passes.
    convergence_checks        = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL)
  )
  if (random) {
    args$random <- list(study = ~ 1 | study, effect = ~ 1 | esid)
    args$prior_heterogeneity <- BayesTools::prior_random(
      study = BayesTools::random_block(sd = if (is.null(fixed)) sd_prior else
        BayesTools::prior("point", list(fixed[["study"]]))),
      effect = BayesTools::random_block(sd = if (is.null(fixed)) sd_prior else
        BayesTools::prior("point", list(fixed[["effect"]])))
    )
  } else {
    args$prior_heterogeneity <- BayesTools::prior("point", list(0))
  }
  do.call(bselmodel.mv, args)
}


.selection_recovery_draws <- function(fit, alias) {

  selection <- BayesTools::parameter_catalog_resolve(
    BayesTools::parameter_catalog(fit$fit), alias = alias, simplify_names = TRUE
  )
  chains <- BayesTools::parameter_draws(fit$fit, selection)
  do.call(cbind, lapply(chains, as.numeric))
}


.selection_recovery_summary <- function(draws, truth, case, parameter) {

  transform   <- switch(parameter, omega = qlogis, study = log, effect = log, identity)
  recovery    <- transform(draws)
  recovery_sd <- sd(as.numeric(recovery))
  data.frame(
    case          = case,
    parameter     = parameter,
    truth         = truth,
    mean          = mean(draws),
    sd            = sd(as.numeric(draws)),
    lower         = unname(quantile(draws, .025)),
    upper         = unname(quantile(draws, .975)),
    mcse          = posterior::mcse_mean(draws),
    rhat          = posterior::rhat(draws),
    ess_bulk      = posterior::ess_bulk(draws),
    ess_tail      = posterior::ess_tail(draws),
    recovery_z    = (mean(recovery) - transform(truth)) / recovery_sd,
    recovery_mcse = posterior::mcse_mean(recovery) / recovery_sd
  )
}


.selection_recovery_run <- function(target) {

  results <- list()
  index   <- 0L
  for (shape in c("diagonal", "structured")) {
    for (random in c(FALSE, TRUE)) {
      index <- index + 1L
      label <- paste(target, shape, if (random) "nested" else "fixed", sep = "/")
      truth <- c(mu = if (random) -.20 else .20, omega = .45,
                 study = if (random) .25 else 0, effect = if (random) .18 else 0)
      design    <- .selection_recovery_design(shape)
      expect_equal(as.numeric(design$V),
                   as.numeric(kronecker(diag(40L), design$block)), tolerance = 1e-14)
      if (shape == "structured") {
        expect_identical(attr(design$V, "RoBMA_vcalc_metadata")$factor_status, "certified")
      }
      simulated <- .selection_recovery_simulate(design, target, truth, 6100L + index)
      started   <- proc.time()[["elapsed"]]
      cat("Selection recovery:", label, "\n")
      fit <- .selection_recovery_fit(design, simulated$data, target, random,
                                     seed = 7100L + index,
                                     sample = if (target == "approximate" &&
                                                  shape == "structured") 24000L else 6000L,
                                     effect_direction = if (random) "negative" else "positive")
      elapsed <- proc.time()[["elapsed"]] - started
      expect_identical(fit$selection_likelihood$type, target, info = label)
      if (target == "exact" && shape == "structured") {
        expect_identical(fit$selection_likelihood$exactness,
                         if (random) "EF" else "E1", info = label)
      }
      aliases <- c(mu = if (random) "mu_intercept" else "mu", omega = "omega[2]")
      if (random) aliases <- c(aliases, study = "study: sd", effect = "effect: sd")
      rows <- lapply(names(aliases), function(parameter) {
        .selection_recovery_summary(
          .selection_recovery_draws(fit, aliases[[parameter]]),
          truth[[parameter]], label, parameter
        )
      })
      result <- do.call(rbind, rows)
      result$elapsed <- elapsed
      result$proposals <- simulated$proposals
      results[[index]] <- result
      print(result, row.names = FALSE, digits = 4)
      expect_true(all(is.finite(as.matrix(result[, 3:13]))), info = label)
      expect_true(all(result$rhat < 1.01), info = label)
      expect_true(all(result$ess_bulk > 1000 & result$ess_tail > 500), info = label)
      expect_true(all(result$mcse / result$sd < .04), info = label)
      # Finite-data recovery is not a comparison to truth at MCMC precision.
      # Use unconstrained coordinates: raw SDs are brittle for skewed positive
      # or bounded parameters. Four posterior SDs is a conservative sentinel
      # across 24 checks, not a claim of calibrated frequentist coverage.
      expect_true(all(abs(result$recovery_z) < 4 + 4 * result$recovery_mcse),
                  info = label)
      # A broad, nearly prior-only posterior is not informative recovery.
      expect_true(all(result$sd < c(mu = .15, omega = .20,
                                    study = .12, effect = .10)[result$parameter]),
                  info = label)
      if (target == "exact" && shape == "structured" && !random) {
        reference <- .selection_recovery_joint_reference(design, simulated$data$yi, 31L)
        refined   <- .selection_recovery_joint_reference(design, simulated$data$yi, 51L)
        expect_equal(reference, refined, tolerance = 1e-6)
        for (parameter in c("mu", "omega")) {
          draws <- .selection_recovery_draws(fit, aliases[[parameter]])
          for (power in 1:2) {
            values <- draws^power
            expect_lte(abs(mean(values) - refined[parameter, power]),
                       4 * posterior::mcse_mean(values) + 1e-6)
          }
        }
        print(refined)
      }
    }
  }
  do.call(rbind, results)
}


.selection_recovery_reference_density <- function(design, yi, truth, target, order) {

  # Independent 3-dimensional oracle, with nuisance parameters fixed. Exact:
  # expand prod[w + (1-w) I(Y>cut)] into eight Gaussian orthant probabilities.
  # Approximate: integrate reciprocal row normalizers over A | Y under the
  # *unselected* Gaussian model. This is Bayes' identity, not the package kernel.
  V          <- design$block
  n          <- length(yi)
  residual   <- design$residual + truth[["effect"]]^2
  latent     <- design$latent + truth[["study"]]^2
  covariance <- V + truth[["study"]]^2 + diag(truth[["effect"]]^2, n)
  threshold  <- qnorm(.05, lower.tail = FALSE) * sqrt(diag(V))
  omega      <- truth[["omega"]]
  log_weight <- sum(log(ifelse(yi > threshold, 1, omega)))
  subsets    <- as.matrix(expand.grid(rep(list(c(FALSE, TRUE)), n)))
  if (target == "approximate") {
    gain <- latent %*% solve(covariance)
    conditional_covariance <- latent - gain %*% latent
    # Algebraically symmetric; do not repair, jitter, or truncate eigenvalues.
    root <- chol(conditional_covariance)
    rule <- statmod::gauss.quad.prob(order, dist = "normal")
    grid <- as.matrix(expand.grid(rep(list(seq_len(order)), n)))
    nodes <- matrix(rule$nodes[grid], ncol = n) %*% root
    weights <- apply(matrix(rule$weights[grid], ncol = n), 1L, prod)
  }
  log_density <- function(mu) {

    location <- rep(mu, n)
    gaussian <- mvtnorm::dmvnorm(yi, location, covariance, log = TRUE)
    if (target == "exact") {
      normalizer <- 0
      for (i in seq_len(nrow(subsets))) {
        selected <- which(subsets[i, ])
        count    <- length(selected)
        probability <- if (count == 0L) 1 else if (count == 1L) {
          pnorm(threshold[selected], mu, sqrt(covariance[selected, selected]),
                lower.tail = FALSE)
        } else {
          as.numeric(mvtnorm::pmvnorm(
            lower = threshold[selected] - mu, upper = rep(Inf, count),
            sigma = covariance[selected, selected, drop = FALSE],
            algorithm = mvtnorm::TVPACK(abseps = 1e-10)
          ))
        }
        normalizer <- normalizer + omega^(n - count) * (1 - omega)^count * probability
      }
      return(gaussian + log_weight - log(normalizer))
    }
    conditional_mean <- location + as.numeric(gain %*% (yi - location))
    integrand <- weights
    for (j in seq_len(n)) {
      score <- (conditional_mean[[j]] + nodes[, j] - threshold[[j]]) / sqrt(residual[[j]])
      integrand <- integrand / (omega + (1 - omega) * pnorm(score))
    }
    gaussian + log_weight + log(sum(integrand))
  }
  log_density
}


.selection_recovery_reference <- function(design, yi, truth, target, order) {

  log_density <- .selection_recovery_reference_density(design, yi, truth, target, order)
  integrand <- function(mu, power = 0L) {

    vapply(mu, function(value) {
      value^power * exp(log_density(value) + dnorm(value, 0, 1, log = TRUE))
    }, numeric(1L))
  }
  integrals <- lapply(0:2, function(power) {
    integrate(integrand, -Inf, Inf, power = power, rel.tol = 1e-7, abs.tol = 1e-9)
  })
  values <- vapply(integrals, `[[`, numeric(1L), "value")
  list(moments = values[2:3] / values[[1L]],
       integration_error = vapply(integrals, `[[`, numeric(1L), "abs.error"))
}


.selection_recovery_joint_reference <- function(design, yi, order) {

  # Independent joint posterior for mu and omega, without random effects.
  # The Gaussian proposal only locates quadrature nodes: importance weights
  # retain the entire posterior, including the logit Jacobian for omega.
  y          <- matrix(yi, ncol = 3L, byrow = TRUE)
  covariance <- design$block
  precision  <- solve(covariance)
  threshold  <- qnorm(.05, lower.tail = FALSE) * sqrt(diag(covariance))
  nonsig     <- sum(sweep(y, 2L, threshold, "<"))
  first_nonsig <- sum(y[1L, ] < threshold)
  gaussian_constant <- -.5 * (length(yi) * log(2 * pi) +
    nrow(y) * as.numeric(determinant(covariance, logarithm = TRUE)$modulus) +
    sum((y %*% precision) * y))
  linear    <- sum(precision %*% colSums(y))
  quadratic <- nrow(y) * sum(precision)
  log_posterior <- function(parameters) {

    mu    <- parameters[[1L]]
    omega <- plogis(parameters[[2L]])
    truth <- c(mu = mu, omega = omega, study = 0, effect = 0)
    first <- .selection_recovery_reference_density(design, y[1L, ], truth, "exact", 25L)(mu)
    first_gaussian <- mvtnorm::dmvnorm(y[1L, ], rep(mu, 3L), covariance, log = TRUE)
    # All blocks share a normalizer. Replace the first block's Gaussian and
    # selection counts by the sufficient statistics for the whole dataset.
    gaussian_constant + linear * mu - .5 * quadratic * mu^2 +
      nrow(y) * (first - first_gaussian) +
      (nonsig - nrow(y) * first_nonsig) * log(omega) +
      dnorm(mu, 0, 1, log = TRUE) + dlogis(parameters[[2L]], log = TRUE)
  }
  mode <- optim(c(mean(yi), 0), function(x) -log_posterior(x), method = "BFGS",
                hessian = TRUE, control = list(reltol = 1e-12))
  stopifnot(mode$convergence == 0L)
  covariance_proposal <- solve(mode$hessian)
  rule <- statmod::gauss.quad.prob(order, dist = "normal")
  grid <- as.matrix(expand.grid(seq_len(order), seq_len(order)))
  nodes <- matrix(rule$nodes[grid], ncol = 2L) %*% chol(covariance_proposal)
  nodes <- sweep(nodes, 2L, mode$par, "+")
  log_weights <- rowSums(log(matrix(rule$weights[grid], ncol = 2L))) +
    apply(nodes, 1L, log_posterior) -
    mvtnorm::dmvnorm(nodes, mode$par, covariance_proposal, log = TRUE)
  weights <- exp(log_weights - max(log_weights))
  weights <- weights / sum(weights)
  values <- cbind(mu = nodes[, 1L], omega = plogis(nodes[, 2L]))
  cbind(first = colSums(values * weights), second = colSums(values^2 * weights))
}


.selection_recovery_oracle_run <- function(target, random) {

  cat("Selection posterior oracle:", target,
      if (random) "known random SDs" else "no random effects", "\n")
  # One truly fully connected dense 3 x 3 V, not a diagonal/factor proxy.
  design <- .selection_recovery_design("dense", studies = 1L)
  truth  <- c(mu = .20, omega = .45, study = if (random) .25 else 0,
              effect = if (random) .18 else 0)
  simulated <- .selection_recovery_simulate(design, target, truth, 6201L)
  # Certify the independent oracle's no-selection limit against conjugate
  # Gaussian updating before using it to judge the JAGS selected posterior.
  unselected <- truth
  unselected[["omega"]] <- 1
  gaussian_reference <- .selection_recovery_reference(
    design, simulated$data$yi, unselected, target, order = 25L
  )
  covariance <- design$block + truth[["study"]]^2 + diag(truth[["effect"]]^2, 3L)
  precision <- solve(covariance)
  variance  <- 1 / (1 + sum(precision))
  mean      <- variance * sum(precision %*% simulated$data$yi)
  expect_equal(gaussian_reference$moments, c(mean, variance + mean^2), tolerance = 1e-7)
  reference <- .selection_recovery_reference(design, simulated$data$yi, truth,
                                             target, order = 25L)
  refined <- .selection_recovery_reference(design, simulated$data$yi, truth,
                                           target, order = 41L)
  expect_equal(reference$moments, refined$moments, tolerance = 1e-6)
  expect_true(all(refined$integration_error < 1e-6))
  started <- proc.time()[["elapsed"]]
  fit <- .selection_recovery_fit(design, simulated$data, target, random = random,
                                 seed = 7201L,
                                 sample = if (target == "exact") 6000L else 24000L,
                                 fixed = truth)
  elapsed <- proc.time()[["elapsed"]] - started
  draws <- .selection_recovery_draws(fit, if (random) "mu_intercept" else "mu")
  expect_lt(abs(mean(draws) - truth[["mu"]]),
            4 * sd(as.numeric(draws)) + 4 * posterior::mcse_mean(draws))
  results <- lapply(1:2, function(power) {
    values <- draws^power
    mcse   <- posterior::mcse_mean(values)
    expect_lt(posterior::rhat(values), 1.01)
    expect_gt(posterior::ess_bulk(values), 1000)
    expect_gt(posterior::ess_tail(values), 500)
    expect_lt(mcse / sd(as.numeric(values)), .04)
    expect_lte(abs(mean(values) - refined$moments[[power]]), 4 * mcse + 1e-6)
    data.frame(target = target, random = random, power = power, jags = mean(values),
               reference = refined$moments[[power]], mcse = mcse, elapsed = elapsed)
  })
  result <- do.call(rbind, results)
  print(result, row.names = FALSE, digits = 6)
  result
}
