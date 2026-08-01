# ============================================================================ #
# IWMDE Aggregation and Monte Carlo Error
# ============================================================================ #

.iwmde_density_aggregate <- function(log_terms, active_mass, denominator) {

  if (!is.matrix(log_terms)) {
    log_terms <- as.matrix(log_terms)
  }
  if (anyNA(log_terms) || any(log_terms == Inf)) {
    stop(
      "IWMDE density aggregation encountered positive-infinite or undefined log terms.",
      call. = FALSE
    )
  }
  denominator <- as.integer(denominator[[1L]])
  if (!is.finite(denominator) || denominator < ncol(log_terms)) {
    denominator <- ncol(log_terms)
  }

  n_grid           <- nrow(log_terms)
  y                <- numeric(n_grid)
  max_log_ratio    <- rep(Inf, n_grid)
  ess              <- numeric(n_grid)
  max_weight_share <- rep(1, n_grid)
  finite           <- is.finite(log_terms)
  finite_terms     <- rowSums(finite)
  contributions    <- matrix(0, nrow = n_grid, ncol = denominator)

  active_rows <- finite_terms > 0L
  if (!any(active_rows)) {
    return(list(
      y                = y,
      finite_terms     = finite_terms,
      max_log_ratio    = max_log_ratio,
      ess              = ess,
      max_weight_share = max_weight_share,
      contributions    = contributions
    ))
  }

  safe_terms <- log_terms
  safe_terms[!finite] <- -Inf
  max_term <- safe_terms[cbind(
    seq_len(n_grid),
    max.col(safe_terms, ties.method = "first")
  )]

  scaled_terms <- exp(safe_terms - max_term)
  scaled_terms[!finite] <- 0
  scaled_terms[!active_rows, ] <- 0

  sum_scaled_terms <- rowSums(scaled_terms)
  sum_scaled_sq    <- rowSums(scaled_terms^2)
  max_scaled       <- scaled_terms[cbind(
    seq_len(n_grid),
    max.col(scaled_terms, ties.method = "first")
  )]

  y[active_rows] <- active_mass * exp(max_term[active_rows]) *
    sum_scaled_terms[active_rows] / denominator
  ess[active_rows] <- sum_scaled_terms[active_rows]^2 /
    sum_scaled_sq[active_rows]
  max_weight_share[active_rows] <- max_scaled[active_rows] /
    sum_scaled_terms[active_rows]
  if (ncol(log_terms) > 0L) {
    contribution_terms          <- active_mass * exp(log_terms)
    contribution_terms[!finite] <- 0
    contributions[, seq_len(ncol(log_terms))] <- contribution_terms
  }

  median_terms <- vapply(which(active_rows), function(row) {
    stats::median(log_terms[row, finite[row, ]])
  }, numeric(1))
  max_log_ratio[active_rows] <- max_term[active_rows] - median_terms

  return(list(
    y                = y,
    finite_terms     = finite_terms,
    max_log_ratio    = max_log_ratio,
    ess              = ess,
    max_weight_share = max_weight_share,
    contributions    = contributions
  ))
}


.iwmde_batch_mcse <- function(contributions) {

  n <- ncol(contributions)
  if (n < 4L) {
    return(list(
      mcse          = rep(NA_real_, nrow(contributions)),
      relative_mcse = rep(NA_real_, nrow(contributions)),
      batch_size    = NA_integer_,
      n_batches     = 0L
    ))
  }

  batch_size <- max(2L, floor(sqrt(n)))
  n_batches  <- floor(n / batch_size)
  if (n_batches < 2L) {
    return(list(
      mcse          = rep(NA_real_, nrow(contributions)),
      relative_mcse = rep(NA_real_, nrow(contributions)),
      batch_size    = batch_size,
      n_batches     = n_batches
    ))
  }

  keep          <- seq_len(n_batches * batch_size)
  contributions <- contributions[, keep, drop = FALSE]
  batch_index   <- rep(seq_len(n_batches), each = batch_size)
  batch_means   <- matrix(NA_real_, nrow = nrow(contributions), ncol = n_batches)

  for (batch in seq_len(n_batches)) {
    batch_means[, batch] <- rowMeans(
      contributions[, batch_index == batch, drop = FALSE]
    )
  }

  mcse <- apply(batch_means, 1, stats::sd) / sqrt(n_batches)
  mean_contribution <- rowMeans(contributions)
  relative_mcse <- rep(NA_real_, length(mcse))
  positive <- mean_contribution > 0
  relative_mcse[positive] <- mcse[positive] / mean_contribution[positive]

  return(list(
    mcse          = mcse,
    relative_mcse = relative_mcse,
    batch_size    = batch_size,
    n_batches     = n_batches
  ))
}


.iwmde_integral_mcse <- function(contributions, x) {

  n <- ncol(contributions)
  if (n < 4L) {
    return(list(mcse = NA_real_, relative_mcse = NA_real_))
  }

  integrals <- .iwmde_trapz_columns(
    x = x,
    y = contributions
  )
  finite <- is.finite(integrals)
  if (sum(finite) < 4L) {
    return(list(mcse = NA_real_, relative_mcse = NA_real_))
  }

  integrals  <- integrals[finite]
  batch_size <- max(2L, floor(sqrt(length(integrals))))
  n_batches  <- floor(length(integrals) / batch_size)
  if (n_batches < 2L) {
    return(list(mcse = NA_real_, relative_mcse = NA_real_))
  }

  keep        <- seq_len(n_batches * batch_size)
  integrals   <- integrals[keep]
  batch_index <- rep(seq_len(n_batches), each = batch_size)
  batch_means <- vapply(seq_len(n_batches), function(batch) {
    mean(integrals[batch_index == batch])
  }, numeric(1))

  mcse          <- stats::sd(batch_means) / sqrt(n_batches)
  mean_integral <- mean(integrals)
  relative_mcse <- if (mean_integral > 0) mcse / mean_integral else NA_real_

  return(list(mcse = mcse, relative_mcse = relative_mcse))
}
