#!/usr/bin/env Rscript
# Independent posterior references for the compiled public joint samplers.
# Rscript tools/verify-joint-samplers.R --out=<directory> [--reps=8] [--draws=6000]
# Add --oracle-only for quadrature without packages or fits. --cases selects a
# regular expression, --baseline disables joint owners, and --seed-offset
# selects a distinct set of fits for checking discrepancies in a new panel.
# A fit replicate is one independent serial chain. The one-chain R-hat notice
# is expected. The comparison reports replicate-jackknife MCSEs and observed
# support, with no pass/fail threshold on realized Monte Carlo errors.

args <- commandArgs(trailingOnly = TRUE)
option <- function(name, default) {
  value <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(value) == 0L) return(default)
  sub(paste0("^--", name, "="), "", value[[1L]])
}
out <- option("out", NULL)
reps <- as.integer(option("reps", "8"))
draws <- as.integer(option("draws", "6000"))
seed_offset <- as.integer(option("seed-offset", "0"))
oracle_only <- "--oracle-only" %in% args
baseline <- "--baseline" %in% args
stopifnot(length(out) == 1L, reps >= 2L, draws >= if (oracle_only) 1L else 100L, seed_offset >= 0L)
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Each cutoff listed here has observations on both sides. The six selection
# curves are explicit fixed-height fixture priors. PET and PEESE slopes and
# positive heterogeneity are sampled, and enter the reference as integrals.
dat <- data.frame(
  sei = c(.075, .10, .14, .19, .24, .30, .36, .41, .55, .63),
  z = c(-1.04, -.18, 1.272, 1.292, 1.635, 1.655, 1.950, 1.970, 1.31, 2.19),
  x = c(-1.25, .32, -.83, .68, -1.46, 1.14, -.48, .87, -.37, 1.38)
)
dat$yi <- dat$sei * dat$z
bias <- list(
  none = list(weight = .83),
  sel_one_025 = list(weight = .25, side = "one-sided", steps = .025, omega = c(1, .31)),
  sel_one_05 = list(weight = .48, side = "one-sided", steps = .05, omega = c(1, .27)),
  sel_two_05 = list(weight = .71, side = "two-sided", steps = .05, omega = c(1, .39)),
  sel_one_three = list(weight = .32, side = "one-sided", steps = c(.025, .05),
                       omega = c(1, .68, .23)),
  sel_one_four = list(weight = .56, side = "one-sided", steps = c(.025, .05, .10),
                      omega = c(1, .77, .51, .19)),
  sel_two_three = list(weight = .61, side = "two-sided", steps = c(.05, .10),
                       omega = c(1, .61, .21)),
  PET = list(weight = .017, scale = 1),
  PEESE = list(weight = .013, scale = 5)
)
effect <- list(weight = c(.72, 1.31), mean = 0, sd = .65)
tau <- list(weight = c(1.17, .64), sd = .25)
moderator <- list(weight = c(.65, 1.47), mean = -.08, sd = .44)
moments <- c("mass", "mu", "mu2", "mu_negative", "tau", "tau2", "PET", "PET2", "PEESE", "PEESE2")
moderator_moments <- c("beta", "beta2", "beta_negative", "mu_beta")

# Golub-Welsch for Legendre nodes and weights. Split the unit prior quantile
# where density changes most, with a panel toward each tail. Agreement is
# reported for two independently evaluated quadrature resolutions.
legendre <- function(n) {
  i <- seq_len(n - 1L)
  off <- i / sqrt(4 * i^2 - 1)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(i, i + 1L)] <- jacobi[cbind(i + 1L, i)] <- off
  spectrum <- eigen(jacobi, symmetric = TRUE)
  list(x = spectrum$values, w = 2 * spectrum$vectors[1L, ]^2)
}

rescale_rule <- function(rule, lower, upper) {
  half <- (upper - lower) / 2
  list(x = lower + half * (1 + rule$x), w = half * rule$w)
}

prior_rule <- function(rule, quantile) {
  edges <- c(0, .03, .2, .5, .8, .95, .99, .999, 1)
  pieces <- lapply(seq_len(length(edges) - 1L), function(i) {
    rescale_rule(rule, edges[i], edges[i + 1L])
  })
  list(x = quantile(unlist(lapply(pieces, `[[`, "x"))),
       w = unlist(lapply(pieces, `[[`, "w")))
}

normal_location <- function(y, variance, slab, location) {
  precision <- 1 / variance
  if (!slab) {
    return(list(log_evidence = sum(dnorm(y, sd = sqrt(variance), log = TRUE)),
                mean = 0, sd = 0))
  }
  posterior_sd <- sqrt(1 / (sum(precision) + 1 / location$sd^2))
  posterior_mean <- posterior_sd^2 * (sum(y * precision) + location$mean / location$sd^2)
  log_evidence <- sum(dnorm(y, posterior_mean, sqrt(variance), log = TRUE)) +
    dnorm(posterior_mean, location$mean, location$sd, log = TRUE) -
    dnorm(posterior_mean, posterior_mean, posterior_sd, log = TRUE)
  list(log_evidence = log_evidence, mean = posterior_mean, sd = posterior_sd)
}

# Integrate every physical interval of the step function against the normal
# distribution. Two-sided curves use both mirrored intervals independently.
selected_factor <- function(mu, sd, se, spec, observed_weight) {
  cutoffs <- qnorm(spec$steps / if (spec$side == "two-sided") 2 else 1,
                  lower.tail = FALSE)
  upper <- c(Inf, cutoffs)
  lower <- c(cutoffs, if (spec$side == "two-sided") 0 else -Inf)
  log_factor <- numeric(length(mu))
  for (study in seq_along(se)) {
    normalizer <- numeric(length(mu))
    for (bin in seq_along(spec$omega)) {
      lo <- se[study] * lower[bin]
      hi <- se[study] * upper[bin]
      probability <- pnorm(hi, mu, sd[study]) - pnorm(lo, mu, sd[study])
      if (spec$side == "two-sided") {
        probability <- probability + pnorm(-lo, mu, sd[study]) - pnorm(-hi, mu, sd[study])
      }
      normalizer <- normalizer + spec$omega[bin] * probability
    }
    log_factor <- log_factor + log(observed_weight[study]) - log(normalizer)
  }
  exp(log_factor)
}

selected_moments <- function(normal, sd, se, spec, observed_weight, rule) {
  if (normal$sd == 0) {
    return(c(selected_factor(0, sd, se, spec, observed_weight), 0, 0, 0))
  }
  # Split at physical zero so the signed probability has a smooth integrand.
  # ±10 standard deviations leave a negligible tail even after multiplication
  # by the maximum possible correction (product of 1 / min(omega)).
  crossing <- max(-10, min(10, -normal$mean / normal$sd))
  integrate_part <- function(lower, upper) {
    nodes <- rescale_rule(rule, lower, upper)
    mu <- normal$mean + normal$sd * nodes$x
    mass <- nodes$w * dnorm(nodes$x) * selected_factor(mu, sd, se, spec, observed_weight)
    c(sum(mass), sum(mu * mass), sum(mu^2 * mass))
  }
  negative <- integrate_part(-10, crossing)
  positive <- integrate_part(crossing, 10)
  c(negative + positive, negative[[1L]])
}

bias_reference <- function(n) {
  rule <- legendre(n)
  location_rule <- legendre(2L * n)
  tau_rule <- prior_rule(rule, function(p) tau$sd * qnorm((1 + p) / 2))
  cells <- expand.grid(bias = names(bias), effect = 0:1, heterogeneity = 0:1,
                       stringsAsFactors = FALSE)
  values <- lapply(seq_len(nrow(cells)), function(i) {
    cell <- cells[i, ]
    label <- cell$bias
    spec <- bias[[label]]
    slab <- cell$effect == 1L
    selection <- !is.null(spec$omega)
    if (selection) {
      p <- if (spec$side == "two-sided") 2 * pnorm(-abs(dat$z)) else pnorm(-dat$z)
      observed_weight <- spec$omega[1L + findInterval(p, spec$steps)]
    }
    slopes <- if (is.null(spec$scale)) list(x = 0, w = 1) else
      prior_rule(rule, function(p) spec$scale * tan(pi * p / 2))
    heterogeneity <- if (cell$heterogeneity == 0L) list(x = 0, w = 1) else tau_rule
    covariate <- switch(label, PET = dat$sei, PEESE = dat$sei^2, rep(0, nrow(dat)))
    at_tau <- lapply(heterogeneity$x, function(t) {
      variance <- dat$sei^2 + t^2
      at_slope <- vapply(slopes$x, function(slope) {
        normal <- normal_location(dat$yi - slope * covariate, variance, slab, effect)
        m <- if (selection) {
          selected_moments(normal, sqrt(variance), dat$sei, spec, observed_weight, location_rule)
        } else {
          c(1, normal$mean, normal$mean^2 + normal$sd^2,
            if (slab) pnorm(-normal$mean / normal$sd) else 0)
        }
        pet <- if (label == "PET") slope else 0
        peese <- if (label == "PEESE") slope else 0
        exp(normal$log_evidence) * c(m, m[1L] * c(t, t^2, pet, pet^2, peese, peese^2))
      }, numeric(length(moments)))
      drop(at_slope %*% slopes$w)
    })
    mass <- spec$weight * effect$weight[cell$effect + 1L] * tau$weight[cell$heterogeneity + 1L]
    mass * drop(do.call(cbind, at_tau) %*% heterogeneity$w)
  })
  cbind(cells, setNames(as.data.frame(do.call(rbind, values)), moments))
}

# Independent dense Gaussian marginal likelihood for the 2×2 moderator and
# intercept gates, evaluated at each positive or zero heterogeneity value.
gaussian_reference <- function(n) {
  tau_rule <- prior_rule(legendre(n), function(p) tau$sd * qnorm((1 + p) / 2))
  cells <- expand.grid(effect = 0:1, moderator = 0:1, heterogeneity = 0:1)
  design <- cbind(1, dat$x)
  values <- lapply(seq_len(nrow(cells)), function(i) {
    cell <- cells[i, ]
    active <- as.logical(c(cell$effect, cell$moderator))
    prior_mean <- c(effect$mean, moderator$mean)[active]
    prior_variance <- c(effect$sd^2, moderator$sd^2)[active]
    x <- design[, active, drop = FALSE]
    at_tau <- function(t) {
      variance <- dat$sei^2 + t^2
      covariance <- diag(variance) + tcrossprod(sweep(x, 2L, sqrt(prior_variance), `*`))
      factor <- chol(covariance)
      whitened <- forwardsolve(t(factor), dat$yi - drop(x %*% prior_mean))
      evidence <- exp(-sum(whitened^2) / 2 - sum(log(diag(factor))) - nrow(dat) * log(2 * pi) / 2)
      means <- variances <- negative <- numeric(2L)
      product <- 0
      if (any(active)) {
        posterior_v <- solve(crossprod(x, x / variance) + diag(1 / prior_variance, ncol(x)))
        means[active] <- drop(posterior_v %*% (crossprod(x, dat$yi / variance) + prior_mean / prior_variance))
        variances[active] <- diag(posterior_v)
        negative[active] <- pnorm(-means[active] / sqrt(variances[active]))
        if (all(active)) product <- prod(means) + posterior_v[1L, 2L]
      }
      squares <- variances + means^2
      evidence * c(1, means[1L], squares[1L], negative[1L], t, t^2, 0, 0, 0, 0,
                   means[2L], squares[2L], negative[2L], product)
    }
    h <- if (cell$heterogeneity == 0L) list(x = 0, w = 1) else tau_rule
    mass <- effect$weight[cell$effect + 1L] * moderator$weight[cell$moderator + 1L] *
      tau$weight[cell$heterogeneity + 1L]
    mass * drop(vapply(h$x, at_tau, numeric(length(c(moments, moderator_moments)))) %*% h$w)
  })
  cbind(cells, setNames(as.data.frame(do.call(rbind, values)), c(moments, moderator_moments)))
}

# The same numerators and denominators define the reference and logged draw
# summaries. Grouping by the scientific branch avoids conflating model labels
# with a particular order of the JAGS categorical prior.
endpoints <- function(cells) {
  has_bias <- "bias" %in% names(cells)
  terms <- if (has_bias) moments[-1L] else c(moments[2:6], moderator_moments)
  all_rows <- rep(TRUE, nrow(cells))
  totals <- function(rows, name) sum(cells[[name]][rows])
  results <- list()
  add <- function(condition, quantity, numerator, denominator) {
    results[[length(results) + 1L]] <<- data.frame(condition = condition, quantity = quantity,
      numerator = numerator, denominator = denominator, stringsAsFactors = FALSE)
  }
  report <- function(condition, rows, arm_moments = TRUE) {
    den <- totals(rows, "mass")
    add(condition, "probability", den, totals(all_rows, "mass"))
    if (!arm_moments) return(invisible(NULL))
    for (name in terms) add(condition, name, totals(rows, name), den)
    add(condition, "effect_on", totals(rows & cells$effect == 1L, "mass"), den)
    add(condition, "heterogeneity_on", totals(rows & cells$heterogeneity == 1L, "mass"), den)
    if (!has_bias) add(condition, "moderator_on", totals(rows & cells$moderator == 1L, "mass"), den)
  }
  report("all", all_rows)
  report("effect_slab", cells$effect == 1L)
  report("heterogeneity_slab", cells$heterogeneity == 1L)
  if (has_bias) {
    report("any_selection", grepl("^sel_", cells$bias))
    for (branch in names(bias)) {
      branch_rows <- cells$bias == branch
      report(branch, branch_rows)
      for (state in 0:1) {
        rows <- branch_rows & cells$effect == state
        key <- paste0(branch, "/effect_", state)
        report(key, rows)
        for (h in 0:1) report(paste0(key, "/tau_", h), rows & cells$heterogeneity == h, FALSE)
      }
    }
  } else {
    report("moderator_slab", cells$moderator == 1L)
    for (effect_state in 0:1) for (mod_state in 0:1) {
      rows <- cells$effect == effect_state & cells$moderator == mod_state
      key <- paste0("effect_", effect_state, "/moderator_", mod_state)
      report(key, rows)
      for (h in 0:1) report(paste0(key, "/tau_", h), rows & cells$heterogeneity == h, FALSE)
    }
  }
  do.call(rbind, results)
}

cases <- c("bias", "bias_permuted", "bias_reflected", "gaussian", "gaussian_permuted")
cases <- grep(option("cases", "."), cases, value = TRUE)
stopifnot(length(cases) > 0L)

# Verifies that the synthetic observations straddle all effective cutoffs.
cutoffs <- sort(unique(unlist(lapply(bias, function(spec) {
  qnorm(spec$steps / if (identical(spec$side, "two-sided")) 2 else 1, lower.tail = FALSE)
}))))
stopifnot(all(vapply(cutoffs, function(z) {
  any(dat$z > z & dat$z < z + .02) && any(dat$z < z & dat$z > z - .02)
}, logical(1L))))

reference_cases <- intersect(c("bias", "bias_reflected", "gaussian"), sub("_permuted$", "", cases))
if (any(grepl("^bias", reference_cases))) {
  bias_low <- bias_reference(24L)
  bias_high <- bias_reference(48L)
}
reflect <- function(cells) {
  # Reflect both the outcomes and their selection direction. Only the signed
  # physical effect changes; there is no mass at zero inside a slab cell.
  cells$mu <- -cells$mu
  cells$mu_negative <- ifelse(cells$effect == 1L, cells$mass - cells$mu_negative, 0)
  cells
}
reference <- lapply(reference_cases, function(case) {
  low <- switch(case, bias = bias_low, bias_reflected = reflect(bias_low), gaussian = gaussian_reference(24L))
  high <- switch(case, bias = bias_high, bias_reflected = reflect(bias_high), gaussian = gaussian_reference(48L))
  targets <- endpoints(high)
  low_targets <- endpoints(low)
  targets$reference <- targets$numerator / targets$denominator
  targets$quadrature_change <- abs(targets$reference - low_targets$numerator / low_targets$denominator)
  # The numerical reference must resolve even the rare conditionals well below
  # the scale of the MCMC panel; the realized MCMC errors have no fixed cutoff.
  stopifnot(all(targets$quadrature_change < 1e-7))
  utils::write.csv(high, file.path(out, paste0(case, "-quadrature-cells.csv")), row.names = FALSE)
  utils::write.csv(targets, file.path(out, paste0(case, "-reference.csv")), row.names = FALSE)
  cat(case, "max quadrature change", format(max(targets$quadrature_change), digits = 4), "\n")
  targets
})
names(reference) <- reference_cases
if (oracle_only) quit(save = "no")

stopifnot(rjags::jags.version() >= "5.0.0", utils::packageVersion("rjags") >= "5.5")
cat("RoBMA", as.character(utils::packageVersion("RoBMA")),
    "BayesTools", as.character(utils::packageVersion("BayesTools")),
    "JAGS", as.character(rjags::jags.version()), "\n")
for (package in c("RoBMA", "BayesTools")) loadNamespace(package)
dlls <- getLoadedDLLs()[c("RoBMA", "BayesTools", "rjags")]
paths <- vapply(dlls, `[[`, character(1L), "path")
utils::write.csv(data.frame(package = names(paths), path = unname(paths), md5 = unname(tools::md5sum(paths))),
                 file.path(out, if (baseline) "runtime-baseline.csv" else "runtime-joint.csv"), row.names = FALSE)

point <- function(weight) BayesTools::prior("spike", list(0), prior_weights = weight)
normal <- function(spec) {
  BayesTools::prior("normal", list(spec$mean, spec$sd), prior_weights = spec$weight[2L])
}

public_fit <- function(case, seed) {
  has_bias <- grepl("^bias", case)
  permuted <- grepl("permuted$", case)
  reflected <- grepl("reflected$", case)
  joint <- if (baseline) "off" else if (has_bias) "auto" else "gaussian"
  if ("sampler.joint" %in% names(RoBMA::RoBMA.options())) {
    RoBMA::RoBMA.options(sampler.joint = joint, sampler.nested.weights = FALSE,
                        sampler.probabilities = character())
  } else if (!baseline) {
    stop("The installed RoBMA does not expose the joint sampler control.", call. = FALSE)
  }
  RoBMA::RoBMA.options(selection.sampler = "slice.coarse")

  input <- dat
  if (reflected) input$yi <- -input$yi
  effect_prior <- list(point(effect$weight[1L]), normal(effect))
  moderator_prior <- list(point(moderator$weight[1L]), normal(moderator))
  tau_prior <- list(point(tau$weight[1L]), BayesTools::prior("normal", list(0, tau$sd),
    truncation = list(0, Inf), prior_weights = tau$weight[2L]))
  order <- if (permuted) 2:1 else 1:2
  common <- list(yi = quote(yi), sei = quote(sei), data = input, measure = "GEN",
    prior_effect_null = effect_prior[[order[1L]]], prior_effect = effect_prior[[order[2L]]],
    prior_heterogeneity_null = tau_prior[[order[1L]]], prior_heterogeneity = tau_prior[[order[2L]]],
    prior_unit_information_sd = 1, chains = 1L, adapt = 1000L, burnin = 2000L,
    sample = draws, parallel = FALSE, seed = seed, silent = TRUE)

  if (has_bias) {
    labels <- if (permuted) c("sel_one_four", "PEESE", "sel_two_three", "none",
      "sel_one_025", "PET", "sel_two_05", "sel_one_three", "sel_one_05") else names(bias)
    bias_priors <- lapply(labels, function(label) {
      spec <- bias[[label]]
      switch(label,
        none = BayesTools::prior_none(prior_weights = spec$weight),
        PET = BayesTools::prior_PET("cauchy", list(0, spec$scale), prior_weights = spec$weight),
        PEESE = BayesTools::prior_PEESE("cauchy", list(0, spec$scale), prior_weights = spec$weight),
        BayesTools::prior_weightfunction(spec$side, spec$steps,
          BayesTools::wf_fixed(spec$omega), prior_weights = spec$weight))
    })
    fit <- do.call(RoBMA::RoBMA, c(common, list(prior_bias = bias_priors[-1L],
      prior_bias_null = bias_priors[[1L]], effect_direction = if (reflected) "negative" else "positive")))
  } else {
    labels <- NULL
    fit <- do.call(RoBMA::BMA.norm, c(common, list(mods = ~ x,
      standardize_continuous_predictors = FALSE,
      prior_mods = list(x = moderator_prior[[order[2L]]]),
      prior_mods_null = list(x = moderator_prior[[order[1L]]]))))
  }
  backend <- fit[["fit"]]
  if (!inherits(backend, "BayesTools_fit")) stop("Public fit failed in ", case, ".", call. = FALSE)
  samplers <- backend[["samplers"]]
  expected <- if (has_bias) "RoBMA::IndependentBiasJoint" else "RoBMA::GaussianJoint"
  sampler <- samplers[samplers$Sampler == expected, , drop = FALSE]
  if (baseline) {
    stopifnot(nrow(sampler) == 0L)
  } else {
    sources <- if (has_bias) "mu" else c("mu_intercept", "mu_x", "tau")
    expected_nodes <- c(paste0(sources, "_indicator"), paste0(sources, "_component_", match(2L, order)))
    if (has_bias) expected_nodes <- c("bias_indicator", expected_nodes)
    if (length(unique(sampler$Index.No)) != 1L || !setequal(expected_nodes, sampler$Node) ||
        anyDuplicated(samplers$Node)) {
      writeLines(capture.output(print(samplers)), con = stderr())
      stop("The expected eligible sampler did not exclusively own the declared gates in ", case,
           ".", call. = FALSE)
    }
  }

  raw <- as.matrix(coda::as.mcmc.list(backend)[[1L]])
  location <- if (has_bias) "mu" else "mu_intercept"
  mu <- raw[, location]
  t <- raw[, "tau"]
  slab_label <- match(2L, order)
  effect_slab <- raw[, paste0(location, "_indicator")] == slab_label
  tau_slab <- raw[, "tau_indicator"] == slab_label
  stopifnot(identical(mu == 0, !effect_slab), identical(t == 0, !tau_slab))
  # Generate the same cell order as the independent reference.
  cell_def <- if (has_bias) {
    label <- labels[raw[, "bias_indicator"]]
    stopifnot(all((raw[, "PET"] == 0) == (label != "PET")),
              all((raw[, "PEESE"] == 0) == (label != "PEESE")))
    expand.grid(bias = names(bias), effect = 0:1, heterogeneity = 0:1, stringsAsFactors = FALSE)
  } else {
    beta <- raw[, "mu_x"]
    moderator_slab <- raw[, "mu_x_indicator"] == slab_label
    stopifnot(identical(beta == 0, !moderator_slab))
    expand.grid(effect = 0:1, moderator = 0:1, heterogeneity = 0:1)
  }
  values <- if (has_bias) {
    cbind(1, mu, mu^2, mu < 0, t, t^2, raw[, "PET"], raw[, "PET"]^2,
          raw[, "PEESE"], raw[, "PEESE"]^2)
  } else {
    cbind(1, mu, mu^2, mu < 0, t, t^2, 0, 0, 0, 0, beta, beta^2, beta < 0, mu * beta)
  }
  counts <- lapply(seq_len(nrow(cell_def)), function(i) {
    cell <- cell_def[i, ]
    rows <- effect_slab == cell$effect & tau_slab == cell$heterogeneity
    rows <- rows & if (has_bias) label == cell$bias else moderator_slab == cell$moderator
    colSums(values[rows, , drop = FALSE])
  })
  quantity_names <- if (has_bias) moments else c(moments, moderator_moments)
  list(endpoints = endpoints(cbind(cell_def, setNames(as.data.frame(do.call(rbind, counts)), quantity_names))),
       samplers = samplers)
}

compare <- function(target, fits) {
  numerator <- do.call(cbind, lapply(fits, `[[`, "numerator"))
  denominator <- do.call(cbind, lapply(fits, `[[`, "denominator"))
  n <- rowSums(numerator)
  d <- rowSums(denominator)
  estimate <- ifelse(d > 0, n / d, NA_real_)
  leave_one <- (n - numerator) / (d - denominator)
  # Replicate jackknife on the pooled ratio preserves the Monte Carlo variation
  # in the rare-family denominator. Support counts are preserved in the table.
  center <- rowMeans(leave_one)
  mcse <- sqrt((reps - 1) / reps * rowSums((leave_one - center)^2))
  mcse[!is.finite(mcse)] <- NA_real_
  error <- estimate - target$reference
  support <- denominator
  events <- target$quantity == "probability"
  support[events, ] <- numerator[events, ]
  binary <- target$quantity %in% c("probability", "mu_negative", "beta_negative",
                                    "effect_on", "moderator_on", "heterogeneity_on")
  rarer_outcome <- ifelse(binary, pmin(n, d - n), NA_real_)
  data.frame(condition = target$condition, quantity = target$quantity, reference = target$reference,
    quadrature_change = target$quadrature_change, estimate = estimate, error = error, mcse = mcse,
    error_mcse = ifelse(mcse > 0, error / mcse, NA_real_), denominator_draws = d,
    support_draws = rowSums(support), min_support_per_rep = apply(support, 1L, min),
    reps_with_support = rowSums(support > 0), rarer_binary_outcome_draws = rarer_outcome)
}

for (case in cases) {
  target <- reference[[sub("_permuted$", "", case)]]
  prefix <- paste0(case, if (baseline) "-baseline" else "-joint")
  fits <- vector("list", reps)
  for (rep in seq_len(reps)) {
    cat(prefix, "replicate", rep, "of", reps, "\n")
    flush.console()
    seed <- seed_offset + 78301L + 1297L * rep + 70001L * match(case, c(
      "bias", "bias_permuted", "bias_reflected", "gaussian", "gaussian_permuted"))
    result <- public_fit(case, seed)
    fits[[rep]] <- transform(result$endpoints, replicate = rep, seed = seed)
    if (rep == 1L) utils::write.csv(result$samplers, file.path(out, paste0(prefix, "-owners.csv")),
                                  row.names = FALSE)
  }
  comparison <- compare(target, fits)
  utils::write.csv(do.call(rbind, fits), file.path(out, paste0(prefix, "-replicates.csv")), row.names = FALSE)
  utils::write.csv(comparison, file.path(out, paste0(prefix, "-comparison.csv")), row.names = FALSE)
  show <- comparison$condition %in% c("all", "PET", "PEESE") & comparison$quantity %in%
    c("probability", "effect_on", "moderator_on", "heterogeneity_on", "mu", "mu_negative", "tau", "beta")
  print(comparison[show, c("condition", "quantity", "reference", "estimate", "mcse", "error_mcse")],
        row.names = FALSE)
}
