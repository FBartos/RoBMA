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
# its support.
.iwmde_qcmde_closed_form_log_prior <- function(x, form) {

  out <- numeric(length(x))
  for (j in seq_along(form[["mean"]])) {
    out <- out + stats::dnorm(x, form[["mean"]][[j]], form[["sd"]][[j]],
                              log = TRUE)
  }
  out[x < form[["lower"]] | x > form[["upper"]]] <- -Inf

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


# Localized positive-tail integral for the closed-form row consumer. This is
# the positive-real Mills expansion with its first omitted term, capped at 64.
.iwmde_qcmde_closed_form_log_mills <- function(a) {

  if (length(a) != 1L || !is.finite(a) || a < 16) {
    return(NA_real_)
  }
  total <- 1
  term  <- 1
  for (index in seq_len(64L)) {
    next_term <- term * (-(2 * index - 1) / a / a)
    if (abs(next_term) <= .Machine$double.eps / 8 * total) {
      if (!is.finite(total) || total <= 0) {
        return(NA_real_)
      }
      return(-log(a) + log(total))
    }
    total <- total + next_term
    term  <- next_term
  }

  return(NA_real_)
}


# log integral_0^width exp(-a*t - t^2/2) dt. Narrow finite intervals use the
# same 16-node Gauss-Legendre machinery as the ordinary normal-mass helper.
.iwmde_qcmde_closed_form_log_tail <- function(a, width) {

  log_mills <- .iwmde_qcmde_closed_form_log_mills(a)
  if (!is.finite(log_mills) || length(width) != 1L || is.na(width) ||
      width <= 0) {
    return(NA_real_)
  }
  if (is.infinite(width)) {
    return(log_mills)
  }
  decay <- a * width + .5 * width * width
  if (is.na(decay) || decay <= 0) {
    return(NA_real_)
  }
  if (decay < .1) {
    rule <- .gauss_legendre_nodes(16L)
    node <- width * rule[["nodes"]]
    terms <- rule[["log_weights"]] - a * width * rule[["nodes"]] -
      .5 * node * node
    top <- max(terms)
    return(log(width) + top + log(sum(exp(terms - top))))
  }
  if (is.infinite(decay)) {
    return(log_mills)
  }
  exponent <- -decay + .iwmde_qcmde_closed_form_log_mills(a + width) -
    log_mills
  if (!is.finite(exponent) || exponent >= 0) {
    return(NA_real_)
  }

  return(log_mills + log(-expm1(exponent)))
}


# Successive Gaussian products retain the target's current coordinate and a
# separate shift. For support wholly in a far tail, evaluate the original row
# at its finite nearest endpoint instead of cancelling height and tail mass.
.iwmde_qcmde_closed_form_normalizer <- function(current, linear, quadratic,
                                                baseline, constant, form,
                                                lower, upper) {

  factor_mean <- form[["mean"]]
  factor_sd   <- form[["sd"]]
  if (!all(is.finite(c(current, linear, quadratic, baseline, constant))) ||
      quadratic <= 0 || length(factor_mean) != length(factor_sd) ||
      any(!is.finite(factor_mean)) || any(!is.finite(factor_sd) | factor_sd <= 0) ||
      is.na(lower) || is.na(upper) || lower >= upper) {
    return(NULL)
  }
  root_quadratic <- sqrt(quadratic)
  kernel_sd <- 1 / root_quadratic
  shift <- linear / quadratic
  sd <- kernel_sd
  if (!is.finite(shift) || !is.finite(sd) || sd <= 0) {
    return(NULL)
  }
  # Save only the overlap inputs here. Its completed height is not evaluated
  # on the boundary-scaled branch, where that unused term may overflow.
  overlap_mean <- overlap_location <- overlap_sd <- numeric(length(factor_sd))
  for (index in seq_along(factor_sd)) {
    offset <- factor_mean[[index]] - current
    scale  <- factor_sd[[index]]
    largest <- max(sd, scale)
    smallest <- min(sd, scale)
    combined_sd <- largest * sqrt(1 + (smallest / largest)^2)
    if (!is.finite(offset) || !is.finite(combined_sd) || combined_sd <= 0) {
      return(NULL)
    }
    overlap_mean[[index]] <- shift
    overlap_location[[index]] <- offset
    overlap_sd[[index]] <- combined_sd
    shift <- if (sd <= scale) {
      shift + (offset - shift) * (sd / combined_sd)^2
    } else {
      offset + (shift - offset) * (scale / combined_sd)^2
    }
    sd <- smallest * (largest / combined_sd)
    if (!is.finite(shift) || !is.finite(sd) || sd <= 0) {
      return(NULL)
    }
  }
  right_distance <- ((lower - current) - shift) / sd
  left_distance  <- (shift - (upper - current)) / sd
  right_tail <- isTRUE(right_distance >= 16)
  left_tail  <- isTRUE(left_distance >= 16)
  if (right_tail || left_tail) {
    boundary <- if (right_tail) lower else upper
    distance <- if (right_tail) right_distance else left_distance
    if (!is.finite(boundary)) {
      return(NULL)
    }
    width <- (upper - lower) / sd
    delta <- boundary - current
    log_boundary <- baseline + constant +
      sum(stats::dnorm(boundary, factor_mean, factor_sd, log = TRUE)) +
      linear * delta - .5 * quadratic * delta * delta
    log_normalizer <- log_boundary + log(sd) +
      .iwmde_qcmde_closed_form_log_tail(distance, width)
  } else {
    score <- linear / root_quadratic
    log_kernel <- baseline + .5 * score * score + log(kernel_sd) + log(2 * pi) / 2
    for (index in seq_along(overlap_sd)) {
      log_kernel <- log_kernel + stats::dnorm(overlap_mean[[index]],
        overlap_location[[index]], overlap_sd[[index]], log = TRUE)
    }
    log_mass <- .iwmde_qcmde_log_normal_mass(lower - current, upper - current,
      shift, sd)
    log_normalizer <- log_kernel + constant + log_mass
  }
  if (!is.finite(log_normalizer)) {
    return(NULL)
  }

  return(list(log_normalizer = log_normalizer, shift = shift, sd = sd))
}


# The exact rows of a pass whose reported log prior is the declared form's up to
# a constant on the lattice `grid`, with their normalizers: the mass of the row
# over its whole support (the prior's support within the chart's). `prior` holds
# the log prior the route reported at the values of `grid` for the rows
# `columns` (NULL when it reported none); `laws` and `kernel` are those of the
# pass.
#
# The Gaussian likelihood and every declared Gaussian prior factor are combined
# successively. The separate current coordinate and shift avoid rounding away a
# narrow row's location. A support wholly beyond sixteen combined SDs uses the
# original row at its nearest finite endpoint with a boundary-scaled tail mass,
# avoiding cancellation of a large completed height and a tiny support mass.
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

  for (row in exact) {
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

    lower <- max(laws[["lower"]][[row]], support[[1L]])
    upper <- min(laws[["upper"]][[row]], support[[2L]])
    normalizer <- .iwmde_qcmde_closed_form_normalizer(
      current, linear, quadratic, baseline, constant, form, lower, upper
    )
    if (is.null(normalizer)) {
      next
    }
    log_normalizer <- normalizer[["log_normalizer"]]
    shift <- normalizer[["shift"]]
    sd <- normalizer[["sd"]]
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
