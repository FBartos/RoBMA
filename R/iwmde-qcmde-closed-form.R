# ============================================================================ #
# qCMDE Closed-Form Normalization of Normal Rows
# ============================================================================ #
#
# A row whose likelihood factor is the Gaussian kernel its route reports and
# whose prior factor is a (truncated) normal in the target is a truncated
# normal in the target, which R/iwmde-qcmde-range.R classifies as "exact". Its
# normalizer is a normal-distribution mass,
#
#   log Z = l(m) + log(2 pi) / 2 + log(s) + log(Phi(b) - Phi(a)),
#
# where l is the row's log density up to its normalizer, m and s are the mean
# and SD of the untruncated normal, and (a, b) is the standardized support. The
# lattice does not decide it, so a row far narrower than a lattice step, or far
# from the other rows, is normalized as accurately as a wide one, and the mass
# the lattice range leaves out of the row is known analytically.
#
# The route evaluates l(x) = baseline + kernel(x) + log_prior(x), where the
# kernel is the Gaussian log likelihood relative to its level at the current
# value and log_prior the row's prior factor, which the route reports next to
# the kernel. The kernel is the same expression in the evaluation and in the
# closed form. The prior factor is checked: the row's law is read from the
# declared prior (its normal form), and the closed form is used only when the
# reported log prior differs from that form by a constant over the lattice,
# which is the row's constant (the prior's normalizer and the baselines of the
# other priors). Any other row, and any route that reports no log prior, keeps
# the lattice normalizer.


# The log density of the Gaussian factors of a prior form at `x`, -Inf outside
# its support when `support` is TRUE.
.iwmde_qcmde_closed_form_log_prior <- function(x, form, support = TRUE) {

  out <- numeric(length(x))
  for (j in seq_along(form[["mean"]])) {
    out <- out + stats::dnorm(x, form[["mean"]][[j]], form[["sd"]][[j]],
                              log = TRUE)
  }
  if (isTRUE(support)) {
    out[x < form[["lower"]] | x > form[["upper"]]] <- -Inf
  }

  return(out)
}


# The row's constant: the difference of the log prior a route reported
# (`evaluated`) and the declared form (`declared`) over a lattice, or NA when
# the two are not a constant apart. Required: at least three finite log
# densities of the row (`log_q`), the same finite / -Inf pattern in the three
# vectors, and a difference that is constant to a few units of rounding in the
# size of the values (and 1e-12) at every finite point. Missing values and +Inf
# are not log densities.
.iwmde_qcmde_closed_form_constant <- function(log_q, evaluated, declared) {

  if (anyNA(log_q) || anyNA(evaluated) || anyNA(declared) ||
      any(log_q == Inf) || any(evaluated == Inf) || any(declared == Inf)) {
    return(NA_real_)
  }
  finite <- is.finite(declared)
  if (sum(is.finite(log_q)) < 3L || !identical(is.finite(log_q), finite) ||
      !identical(is.finite(evaluated), finite)) {
    return(NA_real_)
  }
  difference <- evaluated[finite] - declared[finite]
  size       <- abs(evaluated[finite]) + abs(declared[finite])
  anchor     <- which.min(size)
  allowance  <- 64 * .Machine$double.eps * (size + size[[anchor]]) + 1e-12
  if (!all(abs(difference - difference[[anchor]]) <= allowance)) {
    return(NA_real_)
  }

  return(difference[[anchor]])
}


# The exact rows of a pass whose reported log prior is the declared form's up to
# a constant on the lattice `grid`, with their normalizers: the mass of the row
# over its whole support (the prior's support within the chart's). `prior` holds
# the log prior the route reported at the values of `grid` for the rows
# `columns` (NULL when it reported none); `laws` and `kernel` are those of the
# pass.
#
# With u = x - current, the row's log density is
#
#   l(x) = baseline + constant + sum_j log N(x; m_j, s_j) + L u - Q u^2 / 2
#        = const + B u - P u^2 / 2,
#
# where L and Q are the kernel's coefficients, P = Q + sum_j 1 / s_j^2 and
# B = L + sum_j (m_j - current) / s_j^2. The normal has centre u* = B / P and SD
# 1 / sqrt(P), and
#
#   log Z = baseline + constant + sum_j [-log s_j - log(2 pi) / 2]
#           - sum_j (m_j - current)^2 / (2 s_j^2) + B^2 / (2 P)
#           + log(2 pi) / 2 + log(1 / sqrt(P)) + log(mass of the support),
#
# evaluated from the coefficients and not at a rounded mean, which a row narrower
# than the spacing of the doubles at its centre could not place.
.iwmde_qcmde_closed_form_rows <- function(laws, kernel, grid, prior, columns,
                                          transform) {

  out <- list(
    rows           = integer(),
    log_normalizer = numeric(),
    center         = numeric(),
    shift          = numeric(),
    sd             = numeric(),
    lower          = numeric(),
    upper          = numeric()
  )
  exact <- which(laws[["kind"]] %in% "exact")
  if (length(exact) == 0L || !is.list(kernel) || is.null(prior)) {
    return(out)
  }
  support <- c(transform[["lower"]], transform[["upper"]])

  for (row in intersect(exact, columns)) {
    form      <- laws[["form"]][[row]]
    current   <- kernel[["current"]][[row]]
    linear    <- kernel[["linear"]][[row]]
    quadratic <- kernel[["quadratic"]][[row]]
    baseline  <- kernel[["baseline"]][[row]]
    if (is.null(form) || length(baseline) != 1L ||
        !all(is.finite(c(current, linear, quadratic, baseline)))) {
      next
    }
    constant <- .iwmde_qcmde_closed_form_constant(
      log_q     = grid[["log_q"]][, row],
      evaluated = prior[, match(row, columns)],
      declared  = .iwmde_qcmde_closed_form_log_prior(grid[["x"]], form)
    )
    if (is.na(constant)) {
      next
    }

    offset    <- form[["mean"]] - current
    precision <- quadratic + sum(1 / form[["sd"]]^2)
    slope     <- linear + sum(offset / form[["sd"]]^2)
    sd        <- 1 / sqrt(precision)
    shift     <- slope / precision
    lower     <- max(laws[["lower"]][[row]], support[[1L]])
    upper     <- min(laws[["upper"]][[row]], support[[2L]])
    log_mass  <- .iwmde_qcmde_log_normal_mass(lower - current, upper - current,
                                              shift, sd)
    log_normalizer <- baseline + constant +
      sum(-log(form[["sd"]]) - log(2 * pi) / 2) -
      sum(offset^2 / (2 * form[["sd"]]^2)) + slope^2 / (2 * precision) +
      log(2 * pi) / 2 + log(sd) + log_mass
    if (!is.finite(log_normalizer)) {
      next
    }
    out[["rows"]]           <- c(out[["rows"]], row)
    out[["log_normalizer"]] <- c(out[["log_normalizer"]], log_normalizer)
    out[["center"]]         <- c(out[["center"]], current)
    out[["shift"]]          <- c(out[["shift"]], shift)
    out[["sd"]]             <- c(out[["sd"]], sd)
    out[["lower"]]          <- c(out[["lower"]], lower)
    out[["upper"]]          <- c(out[["upper"]], upper)
  }

  return(out)
}


# The fraction of the mass of each closed-form row that lies in the interval
# `x_range` of a grid.
.iwmde_qcmde_closed_form_range_mass <- function(closed, x_range) {

  if (length(closed[["rows"]]) == 0L) {
    return(numeric())
  }
  total  <- .iwmde_qcmde_log_normal_mass(closed[["lower"]] - closed[["center"]],
    closed[["upper"]] - closed[["center"]], closed[["shift"]], closed[["sd"]])
  inside <- .iwmde_qcmde_log_normal_mass(
    pmax(closed[["lower"]], min(x_range)) - closed[["center"]],
    pmin(closed[["upper"]], max(x_range)) - closed[["center"]],
    closed[["shift"]], closed[["sd"]])

  return(exp(inside - total))
}
