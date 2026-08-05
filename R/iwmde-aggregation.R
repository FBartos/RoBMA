# ============================================================================ #
# IWMDE Aggregation and Monte Carlo Error
# ============================================================================ #

.iwmde_density_aggregate <- function(log_terms, active_mass, denominator,
                                     contribution_rows = NULL,
                                     population_rows = NULL,
                                     chain_id = NULL) {

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

  use_population <- !is.null(contribution_rows) ||
    !is.null(population_rows) || !is.null(chain_id)
  if (use_population &&
      (is.null(contribution_rows) || is.null(population_rows) ||
       is.null(chain_id))) {
    stop(
      "Contribution rows, population rows, and chain IDs must be supplied together.",
      call. = FALSE
    )
  }
  if (use_population) {
    if (length(contribution_rows) != ncol(log_terms) ||
        length(unique(population_rows)) != length(population_rows) ||
        length(chain_id) != length(population_rows) || anyNA(chain_id)) {
      stop("Invalid IWMDE contribution-sequence metadata.", call. = FALSE)
    }
    contribution_index <- match(contribution_rows, population_rows)
    if (anyNA(contribution_index) ||
        length(unique(contribution_index)) != length(contribution_index)) {
      stop(
        "IWMDE contribution rows must identify unique population rows.",
        call. = FALSE
      )
    }
    n_contributions   <- length(population_rows)
    design_factor     <- active_mass * n_contributions / denominator
  } else {
    n_contributions   <- denominator
    contribution_index <- seq_len(ncol(log_terms))
    design_factor     <- active_mass
    chain_id          <- rep(1L, n_contributions)
  }

  n_grid           <- nrow(log_terms)
  y                <- numeric(n_grid)
  max_log_ratio    <- rep(Inf, n_grid)
  ess              <- numeric(n_grid)
  max_weight_share <- rep(1, n_grid)
  finite           <- is.finite(log_terms)
  finite_terms     <- rowSums(finite)
  contributions    <- matrix(0, nrow = n_grid, ncol = n_contributions)
  if (use_population) {
    attr(contributions, "chain_id") <- chain_id
  }

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
    contribution_terms          <- design_factor * exp(log_terms)
    contribution_terms[!finite] <- 0
    contributions[, contribution_index] <- contribution_terms
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
  contribution_sum <- rowSums(contributions)
  contribution_sq  <- rowSums(contributions^2)
  kish_ess <- numeric(nrow(contributions))
  positive_sq <- contribution_sq > 0
  kish_ess[positive_sq] <- contribution_sum[positive_sq]^2 /
    contribution_sq[positive_sq]
  if (n < 4L) {
    return(list(
      mcse          = rep(NA_real_, nrow(contributions)),
      relative_mcse = rep(NA_real_, nrow(contributions)),
      ess           = kish_ess,
      batch_size    = NA_integer_,
      n_batches     = 0L
    ))
  }

  chain_id <- attr(contributions, "chain_id", exact = TRUE)
  if (is.null(chain_id)) {
    chain_id <- rep(1L, n)
  }
  if (length(chain_id) != n || anyNA(chain_id)) {
    stop("Invalid chain IDs for IWMDE contributions.", call. = FALSE)
  }
  chain_rows   <- split(seq_len(n), factor(chain_id, levels = unique(chain_id)))
  chain_length <- vapply(chain_rows, length, integer(1))
  batch_size   <- max(2L, floor(sqrt(min(chain_length))))
  chain_batches <- floor(chain_length / batch_size)
  n_batches     <- sum(chain_batches)
  if (n_batches < 2L) {
    return(list(
      mcse          = rep(NA_real_, nrow(contributions)),
      relative_mcse = rep(NA_real_, nrow(contributions)),
      ess           = kish_ess,
      batch_size    = batch_size,
      n_batches     = n_batches
    ))
  }

  batch_means   <- matrix(NA_real_, nrow = nrow(contributions), ncol = n_batches)
  batch <- 0L
  for (chain in seq_along(chain_rows)) {
    rows <- chain_rows[[chain]]
    for (chain_batch in seq_len(chain_batches[[chain]])) {
      batch <- batch + 1L
      offset <- (chain_batch - 1L) * batch_size
      batch_rows <- rows[offset + seq_len(batch_size)]
      batch_means[, batch] <- rowMeans(
        contributions[, batch_rows, drop = FALSE]
      )
    }
  }

  mcse <- apply(batch_means, 1, stats::sd) / sqrt(n_batches)
  mean_contribution <- rowMeans(contributions)
  relative_mcse <- rep(NA_real_, length(mcse))
  positive <- mean_contribution > 0
  relative_mcse[positive] <- mcse[positive] / mean_contribution[positive]
  contribution_variance <- apply(contributions, 1L, stats::var)
  inefficiency <- rep(1, length(mcse))
  estimated <- contribution_variance > 0 & is.finite(mcse)
  inefficiency[estimated] <- pmax(
    1,
    n * mcse[estimated]^2 / contribution_variance[estimated]
  )
  ess <- kish_ess / inefficiency

  return(list(
    mcse          = mcse,
    relative_mcse = relative_mcse,
    ess           = ess,
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

  chain_id <- attr(contributions, "chain_id", exact = TRUE)
  integrals <- matrix(integrals[finite], nrow = 1L)
  if (!is.null(chain_id)) {
    attr(integrals, "chain_id") <- chain_id[finite]
  }
  mcse_data <- .iwmde_batch_mcse(integrals)

  return(list(
    mcse          = mcse_data[["mcse"]][[1L]],
    relative_mcse = mcse_data[["relative_mcse"]][[1L]]
  ))
}
