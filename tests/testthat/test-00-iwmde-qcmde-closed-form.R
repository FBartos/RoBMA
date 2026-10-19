context("(00) qCMDE closed-form normalization of exact Normal rows")

# References: the log normalizers of the unit rows and the seed ordinate were
# computed with mpmath at 60 significant digits from the definitions (the
# integral of exp(l) over the support by the error function; the mixture of
# the conditional normals), never from the expressions under test.


# A row whose likelihood factor is a Gaussian kernel and whose prior factor is
# a (truncated) normal. As the density route evaluates it, the log density up
# to the normalizer is the kernel around `current` from its level `baseline`
# there, plus the prior factor: the normal log density, less the log mass of its
# truncation, plus a constant of the row (the other priors' baselines) that the
# route reports with the kernel as the log prior.
.cf_row <- function(center, scale, current, baseline, prior_mean, prior_sd,
                    constant, lower = -Inf, upper = Inf) {

  quadratic <- 1 / scale^2
  list(
    current    = current,
    linear     = (center - current) * quadratic,
    quadratic  = quadratic,
    baseline   = baseline,
    prior_mean = prior_mean,
    prior_sd   = prior_sd,
    constant   = constant,
    lower      = lower,
    upper      = upper,
    prior      = BayesTools::prior("normal",
      list(mean = prior_mean, sd = prior_sd),
      truncation = list(lower = lower, upper = upper))
  )
}


.cf_log_prior <- function(row, x) {

  mass <- stats::pnorm(row[["upper"]], row[["prior_mean"]], row[["prior_sd"]]) -
    stats::pnorm(row[["lower"]], row[["prior_mean"]], row[["prior_sd"]])
  out <- stats::dnorm(x, row[["prior_mean"]], row[["prior_sd"]], log = TRUE) -
    log(mass) + row[["constant"]]
  out[x < row[["lower"]] | x > row[["upper"]]] <- -Inf

  out
}


.cf_log_q <- function(row, x, log_prior = .cf_log_prior(row, x)) {

  delta <- x - row[["current"]]

  row[["baseline"]] + row[["linear"]] * delta -
    .5 * row[["quadratic"]] * delta^2 + log_prior
}


.cf_rows <- function() {

  list(
    wide      = .cf_row(.30, .08, .25, -3.7, 0, 10, -2.3),
    truncated = .cf_row(-.10, .30, .05, 1.25, .2, .5, .7, lower = 0),
    narrow    = .cf_row(.40, 1e-12, .4 + 3e-13, 2, 0, 1, 5.1),
    bounded   = .cf_row(.10, .20, .12, -.5, .1, .3, -1.1, lower = -.2,
                        upper = .35)
  )
}


# mpmath references of the rows above: log of the integral of exp(l) over the
# support.
.cf_log_normalizer <- c(
  wide      = -10.63348320748018798319839,
  truncated =   0.9414941626452637882767829,
  narrow    = -20.5660264373951018726759,
  bounded   =  -1.844332369857075677400523
)


# The log prior of all rows at `values`, as a matrix with a column per row.
.cf_log_prior_matrix <- function(rows, values) {

  matrix(unlist(lapply(rows, .cf_log_prior, x = values), use.names = FALSE),
         nrow = length(values))
}


# One normalization pass over `rows` with the joint density replaced by the
# rows' own evaluation. `log_prior` replaces the prior factor the route
# evaluates (and reports); `priors` the prior forms the row states declare;
# `report` the log prior the route attaches to its kernel; `kernel` NULL has
# the route report no kernel.
.cf_pass <- function(rows, log_prior = NULL, priors = NULL, report = TRUE,
                     kernel = "reported", replacement = list(type = "scalar"),
                     n_points = 41L, display_grid = c(-.2, .3),
                     .env = parent.frame()) {

  if (is.null(log_prior)) {
    log_prior <- function(values) .cf_log_prior_matrix(rows, values)
  }
  if (is.null(priors)) {
    priors <- lapply(rows, `[[`, "prior")
  }
  if (identical(kernel, "reported")) {
    kernel <- list(
      current     = vapply(rows, `[[`, numeric(1), "current"),
      linear      = vapply(rows, `[[`, numeric(1), "linear"),
      quadratic   = vapply(rows, `[[`, numeric(1), "quadratic"),
      prior_route = rep("focal", length(rows)),
      baseline    = vapply(rows, `[[`, numeric(1), "baseline")
    )
  }
  testthat::local_mocked_bindings(
    .iwmde_log_q_grid = function(context, parameter, values, row_states,
                                 replacement) {

      prior <- log_prior(values)
      out   <- matrix(unlist(lapply(seq_along(rows), function(i) {
        .cf_log_q(rows[[i]], values, prior[, i])
      }), use.names = FALSE), nrow = length(values))
      if (!is.null(kernel)) {
        if (report) {
          kernel[["log_prior"]] <- prior
        }
        attr(out, "gaussian_kernel") <- kernel
      }
      out
    },
    .iwmde_qcmde_estimate_only = function(context, row_states, rows) FALSE,
    .selection_normalizer_grid = function(context, queries, rows) NULL,
    .selection_normalizer_grid_diagnostics = function(shared) NULL,
    .package = "RoBMA",
    .env     = .env
  )
  values  <- seq(-1, 1, length.out = n_points)
  states  <- lapply(seq_along(rows), function(i) {
    list(row_index = i, focal_prior = priors[[i]], use_focal_prior_delta = TRUE)
  })
  pass <- .iwmde_qcmde_normalization_pass(
    context            = list(),
    parameter          = "mu",
    display_grid       = display_grid,
    normalization_grid = list(x = values, z = values,
                              log_jacobian = rep(0, length(values))),
    transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
    normalization_prob = .999,
    row_states         = states,
    replacement        = replacement,
    estimator_rows     = seq_along(rows)
  )

  list(pass = pass, states = states)
}


# The trapezoid normalizer of a pass's nested grid, which is what a row without
# a closed form gets.
.cf_lattice_normalizer <- function(pass, row) {

  grid <- pass[["nested"]]

  .iwmde_log_trapz_columns(grid[["z"]],
    grid[["log_q"]][, row, drop = FALSE] + grid[["log_jacobian"]])
}


# A normal prior with the parameters of `prior` moved.
.cf_misreport <- function(prior, mean = 0, sd = 1) {

  BayesTools::prior("normal", list(mean = prior$parameters$mean + mean,
                                   sd = prior$parameters$sd * sd),
                    truncation = prior$truncation)
}


test_that("exact Normal rows are normalized by their closed-form mass on every grid", {

  rows <- .cf_rows()
  run  <- .cf_pass(rows)
  pass <- run[["pass"]]

  expect_identical(pass[["closed_form"]][["rows"]], 1:4)
  for (name in c("initial", "nodes", "nested")) {
    expect_equal(unname(pass[[name]][["log_normalizer"]]),
                 unname(pass[["nested"]][["log_normalizer"]]), tolerance = 0)
  }
  expect_lt(max(abs(pass[["nested"]][["log_normalizer"]] -
                      .cf_log_normalizer)), 1e-10)
  # The narrow row (SD 1e-12) sits between lattice nodes: its lattice
  # normalizer is far off, its closed form is not.
  expect_gt(abs(.cf_lattice_normalizer(pass, 3L) - .cf_log_normalizer[["narrow"]]), 1)
  expect_lt(abs(pass[["nested"]][["log_normalizer"]][[3L]] -
                  .cf_log_normalizer[["narrow"]]), 1e-10)
})


test_that("a truncation narrower than the sd of its normal keeps its closed-form normalizer", {

  # The effect prior normal(0, 1) truncated to [0, 1e-15] and two unit
  # observations at 0: the row is a normal with the sd 1 / sqrt(3) on an
  # interval of 1.7e-15 sds, and its ordinate is 1e15 up to 1e-16. The log mass
  # of the interval from the difference of two distribution function values was
  # 3.7% off (the ordinate 1.0373e15). Reference: 60-digit mpmath quadrature of
  # exp(l) over [0, 1e-15] for the log density l of the row.
  prior   <- BayesTools::prior("normal", list(mean = 0, sd = 1),
                               truncation = list(lower = 0, upper = 1e-15))
  current <- 5e-16
  kernel  <- list(current = current, linear = -2 * current, quadratic = 2,
                  prior_route = "focal", baseline = 0)
  states  <- list(list(focal_prior = prior, use_focal_prior_delta = TRUE))
  laws    <- .iwmde_qcmde_row_laws(list(), states, list(type = "scalar"), kernel)
  expect_identical(laws[["kind"]], "exact")

  x         <- c(0, current, 1e-15)
  log_prior <- matrix(stats::dnorm(x, log = TRUE), ncol = 1L)
  delta     <- x - current
  grid      <- list(x = x, log_q = matrix(
    kernel[["linear"]] * delta - delta^2 + log_prior[, 1L], ncol = 1L))
  closed <- .iwmde_qcmde_closed_form_rows(laws, kernel, grid, log_prior, 1L,
    .iwmde_parameter_transform(c(-Inf, Inf)))
  expect_identical(closed[["rows"]], 1L)
  expect_lt(abs(closed[["log_normalizer"]] + 35.457714928115358), 1e-12)
  expect_equal(exp(grid[["log_q"]][2L, 1L] - closed[["log_normalizer"]]),
    999999999999999.92, tolerance = 1e-12)
})


test_that("closed-form rows leave no truncation and report their analytic mass in the range", {

  rows    <- .cf_rows()[c("wide", "truncated", "bounded")]
  pass    <- .cf_pass(rows)[["pass"]]
  summary <- .iwmde_qcmde_normalization_summary(pass, "mu", 1:3, 1, 3L)

  expect_identical(pass[["closed_form"]][["rows"]], 1:3)
  expect_identical(summary[["truncation"]][["max"]], 0)
  expect_identical(summary[["truncation"]][["status"]], "exact")
  expect_identical(summary[["truncation"]][["ordinate_bound"]], 0)
  expect_identical(summary[["pilot_log_normalizer"]], summary[["log_normalizer"]])
  expect_identical(summary[["validation_log_normalizer"]],
                   summary[["log_normalizer"]])

  # The normalization integral is the mean of the rows' masses in the final
  # range, from direct normal-distribution arithmetic on the rows' normals.
  fits <- lapply(rows, function(row) {
    precision <- row[["quadratic"]] + 1 / row[["prior_sd"]]^2
    mean      <- (row[["quadratic"]] * row[["current"]] + row[["linear"]] +
                    row[["prior_mean"]] / row[["prior_sd"]]^2) / precision
    list(mean = mean, sd = 1 / sqrt(precision), lower = row[["lower"]],
         upper = row[["upper"]])
  })
  range <- summary[["normalization_range"]]
  inside <- vapply(fits, function(fit) {
    mass <- function(lower, upper) {
      stats::pnorm(upper, fit$mean, fit$sd) - stats::pnorm(lower, fit$mean, fit$sd)
    }
    mass(max(fit$lower, range[[1L]]), min(fit$upper, range[[2L]])) /
      mass(fit$lower, fit$upper)
  }, numeric(1))
  expect_equal(summary[["final_normalization_integral"]], mean(inside),
               tolerance = 1e-12)
  expect_lt(summary[["final_normalization_integral"]], 1)
  expect_gt(summary[["final_normalization_integral"]], .99)
})


test_that("a row whose reported prior is not a constant from its declared form keeps the lattice normalizer", {

  rows <- .cf_rows()
  # A misreported prior: the route evaluates the row's prior, the row declares
  # another. The narrow row's lattice points all lie hundreds of widths from
  # its centre, where a check relative to the size of the log density passes
  # anything.
  for (misreport in list(list(sd = 2), list(mean = .5), list(mean = .01))) {
    priors <- lapply(rows, `[[`, "prior")
    priors[[1L]] <- do.call(.cf_misreport, c(list(priors[[1L]]), misreport))
    priors[[3L]] <- do.call(.cf_misreport, c(list(priors[[3L]]), misreport))
    pass <- .cf_pass(rows, priors = priors)[["pass"]]
    expect_identical(pass[["closed_form"]][["rows"]], c(2L, 4L))
    expect_equal(pass[["nested"]][["log_normalizer"]][[1L]],
                 .cf_lattice_normalizer(pass, 1L), tolerance = 0)
    expect_equal(pass[["nested"]][["log_normalizer"]][[3L]],
                 .cf_lattice_normalizer(pass, 3L), tolerance = 0)
    expect_equal(unname(pass[["nested"]][["log_normalizer"]][c(2L, 4L)]),
                 unname(.cf_log_normalizer[c("truncated", "bounded")]),
                 tolerance = 1e-10)
  }

  # A prior that differs from the declared form by something that is not a
  # constant, however small.
  tilted <- function(values) {
    out <- .cf_log_prior_matrix(rows, values)
    out[, 4L] <- out[, 4L] + 1e-6 * values
    out
  }
  pass <- .cf_pass(rows, log_prior = tilted)[["pass"]]
  expect_identical(pass[["closed_form"]][["rows"]], 1:3)
})


test_that("a row whose evaluated support differs from its law's keeps the lattice normalizer", {

  rows <- .cf_rows()[c("wide", "bounded")]
  # Narrower: -Inf where the declared form is finite.
  narrower <- function(values) {
    out <- .cf_log_prior_matrix(rows, values)
    out[values < -.15, 1L] <- -Inf
    out
  }
  pass <- .cf_pass(rows, log_prior = narrower)[["pass"]]
  expect_identical(pass[["closed_form"]][["rows"]], 2L)
  expect_equal(pass[["nested"]][["log_normalizer"]][[1L]],
               .cf_lattice_normalizer(pass, 1L), tolerance = 0)
  expect_equal(pass[["nested"]][["log_normalizer"]][[2L]],
               .cf_log_normalizer[["bounded"]], tolerance = 1e-10)

  # Wider: finite where the declared form is -Inf (the bounded row's support
  # ends at .35, inside the lattice of the wide row).
  wider <- function(values) {
    out <- .cf_log_prior_matrix(rows, values)
    out[values > .35, 2L] <- 0
    out
  }
  pass <- .cf_pass(rows, log_prior = wider)[["pass"]]
  expect_identical(pass[["closed_form"]][["rows"]], 1L)
})


test_that("rows without an exact law, a reported log prior or a kernel keep the lattice normalizer", {

  normal <- .cf_rows()[["bounded"]]
  rows   <- list(normal, normal)
  laplace <- function(values) {
    cbind(.cf_log_prior(normal, values), -abs(values - .1) / .5)
  }
  kernel <- list(
    current     = c(normal[["current"]], NA),
    linear      = c(normal[["linear"]], NA),
    quadratic   = c(normal[["quadratic"]], NA),
    prior_route = c("focal", NA),
    baseline    = c(normal[["baseline"]], NA)
  )
  # The second row is a Laplace row whose route reports no kernel for it.
  pass <- .cf_pass(rows, log_prior = laplace, kernel = kernel)[["pass"]]
  expect_identical(pass[["closed_form"]][["rows"]], 1L)
  expect_equal(pass[["nested"]][["log_normalizer"]][[1L]],
               .cf_log_normalizer[["bounded"]], tolerance = 1e-10)
  expect_equal(pass[["nested"]][["log_normalizer"]][[2L]],
               .cf_lattice_normalizer(pass, 2L), tolerance = 0)

  # A route that reports the kernel without its log prior has no eligible row,
  # and one that reports no kernel leaves every row to the lattice.
  for (case in list(list(report = FALSE), list(kernel = NULL))) {
    pass <- do.call(.cf_pass, c(list(list(normal)), case))[["pass"]]
    expect_identical(pass[["closed_form"]][["rows"]], integer())
    expect_equal(pass[["nested"]][["log_normalizer"]][[1L]],
                 .cf_lattice_normalizer(pass, 1L), tolerance = 0)
  }
})


test_that("the pass treats the ordered-source replacement like any other for rows whose route reports their log prior", {

  # The ordered-source producer reports no log prior (its exact rows are
  # normalized by its own source normalizer before a pass), so none of its
  # rows is eligible in practice; the pass itself makes no exception for the
  # replacement.
  rows <- .cf_rows()[c("wide", "bounded")]
  pass <- .cf_pass(rows, replacement = list(type = "ordered_source"))[["pass"]]

  expect_identical(pass[["closed_form"]][["rows"]], 1:2)
  expect_equal(unname(pass[["nested"]][["log_normalizer"]]),
               unname(.cf_log_normalizer[c("wide", "bounded")]), tolerance = 1e-10)
})


test_that("the row constant is checked at every finite lattice point", {

  constant <- .iwmde_qcmde_closed_form_constant
  declared <- c(-1000, -3, -Inf, 2.5, 0)
  shift    <- 4.25
  log_q    <- ifelse(is.finite(declared), declared - 7, -Inf)

  # The constant is recovered.
  expect_equal(constant(log_q, declared + shift, declared), shift,
               tolerance = 1e-12)
  # Constant to rounding in the size of the values, and no more.
  expect_equal(constant(log_q, (declared + shift) * (1 + 20 * .Machine$double.eps),
                        declared), shift, tolerance = 1e-9)
  expect_identical(constant(log_q, declared + shift + 1e-6 * seq_along(declared),
                            declared), NA_real_)
  expect_identical(constant(log_q, declared + shift + c(0, 0, 0, 1e-9, 0),
                            declared), NA_real_)
  # The finite / -Inf pattern must agree in both directions, in every vector.
  expect_identical(constant(log_q, c(-1000, -3, 0, 2.5, 0) + shift, declared),
                   NA_real_)
  expect_identical(constant(log_q, c(-Inf, -3, -Inf, 2.5, 0) + shift, declared),
                   NA_real_)
  expect_identical(constant(c(-Inf, -3, -Inf, 2.5, 0), declared + shift, declared),
                   NA_real_)
  # At least three finite points; missing values and +Inf are not log densities.
  few <- c(-Inf, -3, -Inf, 2.5, -Inf)
  expect_identical(constant(few, few + shift, few), NA_real_)
  expect_identical(constant(c(NA, -3, 0, 2.5, 1), declared + shift, declared),
                   NA_real_)
  expect_identical(constant(log_q, c(Inf, -3, -Inf, 2.5, 0), declared), NA_real_)
})


test_that("a cached normalization with closed-form rows serves a second ordinate", {

  .iwmde_qcmde_cache_clear()
  withr::defer(.iwmde_qcmde_cache_clear())
  rows <- .cf_rows()[c("wide", "truncated", "bounded")]
  run  <- .cf_pass(rows)
  cache <- list(cache = new.env(), key = "closed-form-entry",
                display_in_key = FALSE)
  assign(".order", character(), envir = cache[["cache"]])
  values <- seq(-1, 1, length.out = 41L)
  normalization <- function(display_grid) {
    .iwmde_qcmde_normalization(
      context            = list(),
      parameter          = "mu",
      display_grid       = display_grid,
      normalization_grid = list(x = values, z = values,
                                log_jacobian = rep(0, length(values))),
      transform          = .iwmde_parameter_transform(c(-Inf, Inf)),
      normalization_prob = .999,
      row_states         = run[["states"]],
      replacement        = list(type = "scalar"),
      estimator_rows     = 1:3,
      active_mass        = 1,
      n_candidate_rows   = 3L,
      cache              = cache
    )
  }
  first  <- normalization(c(-.2, .3))
  second <- normalization(c(.05, .15))
  expect_identical(second[["summary"]], first[["summary"]])
  expect_equal(unname(second[["summary"]][["log_normalizer"]]),
               unname(.cf_log_normalizer[c("wide", "truncated", "bounded")]),
               tolerance = 1e-10)

  # The second ordinate is the mixture of the rows' normals at its values.
  density <- .iwmde_qcmde_pilot_density(
    second[["log_q_display"]], second[["summary"]][["log_normalizer"]],
    is.finite(second[["summary"]][["log_normalizer"]]), 1, 3L)
  expected <- vapply(c(.05, .15), function(x) {
    mean(vapply(seq_along(rows), function(i) {
      exp(.cf_log_q(rows[[i]], x) - .cf_log_normalizer[[names(rows)[[i]]]])
    }, numeric(1)))
  }, numeric(1))
  expect_equal(density, expected, tolerance = 1e-12)
})


# ---------------------------------------------------------------------------- #
# The seed replay: 499 narrow rows and one broad row on the public route
# ---------------------------------------------------------------------------- #

# A Normal random-effects model with two equally precise observations whose
# rows are the conditional laws of mu given tau: 499 rows at tau = 0 (SD
# `step / divisor`) and one at tau = 1 (SD 0.58). The target value sits at a
# narrow row's mean, a quarter of a broad-lattice step from a node. The
# independent reference is the mixture of the rows' conditional normals
# (mpmath, 60 digits).
.cf_seed <- function(divisor = 32) {

  control <- .density_control_normalize("qCMDE")
  control <- .hypothesis_plan_density_control(control)
  N <- control[["normalization_points"]]
  S <- control[["samples"]]
  K <- 2L
  tau_narrow <- 0
  tau_broad  <- 1
  q <- stats::qnorm((1 - control[["normalization_prob"]]) / 2, lower.tail = FALSE)
  h <- 2 * q / (N - 1L) / divisor
  sei <- sqrt(K * (h^2 / (1 - h^2)) - tau_narrow^2)
  precision <- function(tau) K / (sei^2 + tau^2)
  A <- precision(c(tau_narrow, tau_broad))
  conditional_sds <- 1 / sqrt(1 + A)
  mean_factors    <- A / (1 + A)
  step   <- 2 * q * conditional_sds[[2L]] / (N - 1L)
  offset <- .25 * step
  ybar   <- offset / diff(rev(mean_factors))
  conditional_means <- mean_factors * ybar
  object <- BMA.norm(
    yi = rep(ybar, K), sei = rep(sei, K), measure = "GEN",
    prior_effect = BayesTools::prior("normal", list(mean = 0, sd = 1)),
    prior_effect_null = NULL,
    prior_heterogeneity = BayesTools::prior("cauchy",
      list(location = 0, scale = 1), truncation = list(lower = 0, upper = Inf)),
    prior_heterogeneity_null = BayesTools::prior("point", list(location = 0)),
    only_priors = TRUE, silent = TRUE
  )
  samples <- cbind(
    mu  = c(conditional_means[[1L]] + conditional_sds[[1L]] *
              stats::qnorm((seq_len(S - 1L) - .5) / (S - 1L)),
            conditional_means[[2L]]),
    tau = c(rep(tau_narrow, S - 1L), tau_broad)
  )
  priors     <- .create_fit_priors(object[["data"]], object[["priors"]])
  components <- attr(priors[["tau"]], "components")
  samples <- cbind(samples,
    tau_indicator = c(rep(which(components == "null"), S - 1L),
                      which(components == "alternative")),
    mu_indicator  = rep(1, S))
  # The broad row goes into a complete batch of the batch-means MCSE estimator.
  samples <- samples[append(seq_len(S - 1L), S, after = 249L), , drop = FALSE]
  fit <- coda::mcmc.list(coda::mcmc(samples))
  attr(fit, "formula_design") <- object[["formula_design"]]
  attr(fit, "prior_list")     <- priors
  object[["fit"]] <- as_bayestools_fit(fit)

  list(object = object, control = control, S = S, N = N, K = K,
       value = conditional_means[[1L]], means = conditional_means,
       sds = conditional_sds)
}


.cf_seed_mixture <- function(seed, x) {

  ((seed$S - 1L) * stats::dnorm(x, seed$means[[1L]], seed$sds[[1L]]) +
     stats::dnorm(x, seed$means[[2L]], seed$sds[[2L]])) / seed$S
}


# The ordinate at `value`, with the closed-form rows of the normalization pass
# it ran (none when the cached normalization served it). A failed estimate is
# returned as status "error".
.cf_seed_ordinate <- function(seed, value, .env = parent.frame()) {

  recorded <- new.env()
  original <- .iwmde_qcmde_normalization_pass
  testthat::local_mocked_bindings(
    .iwmde_qcmde_normalization_pass = function(...) {
      pass <- original(...)
      recorded$rows <- pass[["closed_form"]][["rows"]]
      pass
    },
    .package = "RoBMA",
    .env     = .env
  )
  context <- .iwmde_context(seed$object)
  plan <- .iwmde_plan(context, "mu", "qCMDE", seed$control,
                      outputs = "ordinate", values = value)
  expect_identical(plan[["status"]], "ok")
  ordinate <- tryCatch(
    .iwmde_estimate(context, "mu", "qCMDE", seed$control,
      outputs = "ordinate", values = value)[["posterior_ordinate"]],
    error = function(e) list(status = "error", message = conditionMessage(e))
  )
  ordinate[["closed_form_rows"]] <- recorded$rows

  ordinate
}


test_that("narrow exact rows beside a broad one are not lost between lattice nodes", {

  .iwmde_qcmde_cache_clear()
  withr::defer(.iwmde_qcmde_cache_clear())
  seed <- .cf_seed()
  expect_equal(c(seed$N, seed$S), c(100L, 500L))
  value <- seed$value

  # mpmath, 60 digits: ((S - 1) N(value; m1, s1) + N(value; m2, s2)) / S.
  expected <- 191.660886315168541816153238827189813875988646
  expect_equal(.cf_seed_mixture(seed, value), expected, tolerance = 1e-12)

  ordinate <- .cf_seed_ordinate(seed, value)
  expect_identical(ordinate[["status"]], "ok")
  expect_equal(ordinate[["ordinate"]], expected, tolerance = 1e-6)
  expect_length(ordinate[["closed_form_rows"]], seed$S)
  diagnostics <- ordinate[["diagnostics"]]
  expect_identical(diagnostics[["normalization_truncation"]], 0)
  expect_identical(diagnostics[["normalization_truncation_status"]], "exact")
  expect_equal(diagnostics[["ordinate_relative_change"]], 0)

  # A second ordinate on the same fit reuses the cached normalization (no pass
  # runs), and is exact as well.
  second_value <- value + .5 * seed$sds[[1L]]
  second <- .cf_seed_ordinate(seed, second_value)
  expect_equal(second[["ordinate"]], .cf_seed_mixture(seed, second_value),
               tolerance = 1e-6)
  expect_null(second[["closed_form_rows"]])

  # The public test of the null at the target value: mpmath, 60 digits,
  # BF10 = N(value; 0, 1) / ordinate.
  .iwmde_qcmde_cache_clear()
  out <- hypothesis(seed$object, paste0("mu = ", sprintf("%.17g", value)),
    standardized_coefficients = TRUE, density_method = "qCMDE", columns = "all")
  expect_equal(unname(attr(out, "raw_BF")), 0.00208063860833742663616681819896,
               tolerance = 1e-6)
})


test_that("rows a hundred million times narrower than the step are normalized by the declared prior, which is checked", {

  .iwmde_qcmde_cache_clear()
  withr::defer(.iwmde_qcmde_cache_clear())
  seed <- .cf_seed(divisor = 1e8)
  value <- seed$value

  # mpmath, 60 digits, as above with the rows' SD = step / 1e8. Every lattice
  # point lies at least 1e7 widths from a narrow row's centre.
  expected <- 598935951.668174094285044618791566085282258921
  expect_equal(.cf_seed_mixture(seed, value), expected, tolerance = 1e-12)
  ordinate <- .cf_seed_ordinate(seed, value)
  expect_equal(ordinate[["ordinate"]], expected, tolerance = 1e-6)
  expect_length(ordinate[["closed_form_rows"]], seed$S)

  # A prior form misreported to the closed form (the wrong scale or location)
  # is not a constant from what the route evaluated, however far from the rows'
  # centres the lattice lies: no row is eligible and the lattice keeps them (the
  # lattice cannot resolve them, so the ordinate is not the exact one).
  original <- .iwmde_qcmde_normal_prior_parameters
  for (misreport in list(function(form) {form$sd <- 2 * form$sd; form},
                         function(form) {form$mean <- form$mean + .5; form})) {
    .iwmde_qcmde_cache_clear()
    testthat::local_mocked_bindings(
      .iwmde_qcmde_normal_prior_parameters = function(prior) {
        form <- original(prior)
        if (is.null(form)) form else misreport(form)
      },
      .package = "RoBMA"
    )
    misreported <- .cf_seed_ordinate(seed, value)
    expect_identical(misreported[["closed_form_rows"]], integer())
    expect_false(isTRUE(all.equal(misreported[["ordinate"]], expected,
                                  tolerance = 1e-6)))
  }
})

test_that("closed Normal rows retain narrow priors and scaled support tails", {

  evaluate <- function(current, linear, quadratic, baseline, form, values, prior) {
    delta <- values - current
    log_q <- baseline + linear * delta - .5 * quadratic * delta^2 + prior
    result <- .iwmde_qcmde_closed_form_rows(
      laws = list(kind = "exact", form = list(form), lower = form$lower, upper = form$upper),
      kernel = list(current = current, linear = linear, quadratic = quadratic, baseline = baseline),
      grid = list(x = values, log_q = matrix(log_q, ncol = 1L)),
      prior = matrix(prior, ncol = 1L), columns = 1L,
      transform = list(lower = form$lower, upper = form$upper))
    expect_identical(result$rows, 1L)
    expect_identical(result$center, current)
    result$log_normalizer
  }
  focal <- function(row, values, prior = .cf_log_prior(row, values)) {
    form <- .iwmde_qcmde_normal_prior_parameters(row$prior)
    evaluate(row$current, row$linear, row$quadratic, row$baseline, form, values, prior)
  }
  # Original narrow-prior fixture: kernelSD.5/current.1/baseline-1.3 and
  # declared unbounded Normal(0,1e-8) prior. The independent integral scales
  # only the prior coordinate; it never uses the completed Gaussian mean.
  row <- .cf_row(.1, .5, .1, -1.3, 0, 1e-8, 0)
  reference <- stats::integrate(function(u) stats::dnorm(u) *
    exp(-1.3 - 2 * (1e-8 * u - .1)^2), -Inf, Inf, rel.tol = 1e-12, abs.tol = 0)
  expected <- log(reference$value)
  allowance <- reference$abs.error / reference$value +
    128 * .Machine$double.eps * max(1, abs(expected))
  expect_lte(abs(focal(row, c(-3, -1, 1, 3) * 1e-8) - expected), allowance)

  # Original half-normal T1e8 fixture, ordinaryT3, and its reflected support.
  # Their physical row is sqrt(2/pi)*exp(-T*x-x^2); u=T*x supplies the reference.
  for (T in c(1e8, 3)) for (direction in c(1, -1)) {
    current <- direction / T
    row <- .cf_row(-direction * T, 1, current,
      -T * abs(current) - .5 * current^2, 0, 1, 0,
      lower = if (direction == 1) 0 else -Inf,
      upper = if (direction == 1) Inf else 0)
    reference <- stats::integrate(function(u) exp(-u - (u / T)^2),
      0, Inf, rel.tol = 1e-12, abs.tol = 0)
    expected <- log(sqrt(2 / pi)) - log(T) + log(reference$value)
    allowance <- reference$abs.error / reference$value +
      128 * .Machine$double.eps * max(1, abs(expected))
    expect_lte(abs(focal(row, direction * c(0, 1, 2, 3) / T) - expected), allowance)
  }

  # A finite near-boundary interval exercises the actual16-node helper. The
  # evaluated prior uses its explicitly declared physical normalizing mass,
  # avoiding a subtractive-CDF artifact in the older general test helper.
  T <- 1e8
  upper <- 1e-10
  current <- upper / 2
  row <- .cf_row(-T, 1, current, -T * current - .5 * current^2,
    0, 1, 0, lower = 0, upper = upper)
  values <- seq(0, upper, length.out = 4L)
  prior_integral <- stats::integrate(function(u) stats::dnorm(upper * u),
    0, 1, rel.tol = 1e-12, abs.tol = 0)
  prior_mass <- upper * prior_integral$value
  prior <- stats::dnorm(values, log = TRUE) - log(prior_mass)
  reference <- stats::integrate(function(u) exp(-T * upper * u - (upper * u)^2),
    0, 1, rel.tol = 1e-12, abs.tol = 0)
  expected <- log(upper) - log(2 * pi) / 2 + log(reference$value) - log(prior_mass)
  allowance <- reference$abs.error / reference$value +
    prior_integral$abs.error / prior_integral$value +
    128 * .Machine$double.eps * max(1, abs(expected))
  expect_lte(abs(focal(row, values, prior) - expected), allowance)

  # Independent80-digit Gaussian-product constants retained from177. These
  # are explicit zero/two-factor declarations of the same row-level contract.
  forms <- list(
    list(mean = numeric(), sd = numeric(), lower = -Inf, upper = Inf),
    list(mean = c(0, .2), sd = c(.7, 1.1), lower = -Inf, upper = Inf))
  expected <- c(-1.0542086473552726098253773208085074,
    -2.9296430659478288107944095688822882)
  values <- c(-.3, -.1, .2, .4)
  for (index in seq_along(forms)) {
    form <- forms[[index]]
    prior <- numeric(length(values))
    for (factor in seq_along(form$mean)) {
      prior <- prior + stats::dnorm(values, form$mean[[factor]], form$sd[[factor]], log = TRUE)
    }
    actual <- evaluate(.1, .4, 4, -1.3, form, values, prior)
    expect_lte(abs(actual - expected[[index]]),
      128 * .Machine$double.eps * max(1, abs(expected[[index]])))
  }
})

test_that("rejected exact laws report estimated tails on their actual retained lattice", {

  rows <- list(.cf_row(.3, .8, .25, -3.7, 0, 10, -2.3))
  declared <- list(.cf_misreport(rows[[1L]]$prior, sd = 2))
  pass <- .cf_pass(rows, priors = declared)$pass
  expect_identical(pass$closed_form$rows, integer())
  expect_identical(pass$laws$kind, "estimate")
  expect_equal(pass$extension_steps, c(lower = 0L, upper = 0L))
  grid <- pass$nested
  # Endpoint exponential envelopes from the independently evaluated row.
  values <- .cf_log_q(rows[[1L]], grid$x)
  n <- length(values)
  edge <- c(1L, n)
  adjacent <- c(2L, n - 1L)
  slopes <- (values[adjacent] - values[edge]) /
    abs(grid$z[adjacent] - grid$z[edge])
  expect_true(all(slopes > 0))
  relative_tail <- sum(exp(values[edge] - grid$log_normalizer) / slopes)
  expected <- relative_tail / (1 + relative_tail)
  summary <- .iwmde_qcmde_normalization_summary(pass, "mu", 1L, 1, 1L)
  expect_identical(summary$truncation$status, "estimate")
  expect_equal(summary$truncation$max, expected,
    tolerance = 64 * .Machine$double.eps)
  expect_equal(summary$normalization_range, range(grid$x), tolerance = 0)
})
