# ============================================================================ #
# IWMDE Aggregation and Monte Carlo Error
# ============================================================================ #

.iwmde_density_aggregate <- function(log_terms, active_mass, denominator,
                                     contribution_rows = NULL,
                                     sampling_population_rows = NULL,
                                     chain_id = NULL,
                                     expected_chain_ids = NULL) {

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

  use_sampling <- !is.null(contribution_rows) ||
    !is.null(sampling_population_rows) || !is.null(chain_id)
  if (use_sampling &&
      (is.null(contribution_rows) || is.null(sampling_population_rows) ||
       is.null(chain_id))) {
    stop(
      "Contribution rows, sampling-population rows, and chain IDs must be supplied together.",
      call. = FALSE
    )
  }
  if (use_sampling) {
    if (!is.finite(denominator) || denominator != ncol(log_terms)) {
      stop(
        "Self-weighting IWMDE row sampling requires the selected-row denominator.",
        call. = FALSE
      )
    }
    if (is.null(expected_chain_ids)) {
      expected_chain_ids <- unique(chain_id)
    }
    if (length(contribution_rows) != ncol(log_terms) ||
        length(unique(contribution_rows)) != length(contribution_rows) ||
        length(unique(sampling_population_rows)) !=
          length(sampling_population_rows) ||
        length(chain_id) != length(contribution_rows) || anyNA(chain_id) ||
        length(expected_chain_ids) == 0L || anyNA(expected_chain_ids) ||
        anyNA(match(chain_id, expected_chain_ids)) ||
        anyNA(match(contribution_rows, sampling_population_rows))) {
      stop("Invalid IWMDE contribution-sampling metadata.", call. = FALSE)
    }
  } else {
    if (!is.finite(denominator) || denominator < ncol(log_terms)) {
      denominator <- ncol(log_terms)
    }
    contribution_rows        <- seq_len(ncol(log_terms))
    sampling_population_rows <- contribution_rows
    chain_id                 <- rep(1L, ncol(log_terms))
    expected_chain_ids       <- 1L
  }

  n_grid           <- nrow(log_terms)
  y                <- numeric(n_grid)
  max_log_ratio    <- rep(Inf, n_grid)
  ess              <- numeric(n_grid)
  max_weight_share <- rep(1, n_grid)
  finite           <- is.finite(log_terms)
  finite_terms     <- rowSums(finite)
  contributions    <- matrix(0, nrow = n_grid, ncol = ncol(log_terms))

  active_rows <- finite_terms > 0L
  if (any(active_rows)) {
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
      design_factor <- active_mass * ncol(log_terms) / denominator
      contributions <- design_factor * exp(log_terms)
      contributions[!finite] <- 0
    }

    median_terms <- vapply(which(active_rows), function(row) {
      stats::median(log_terms[row, finite[row, ]])
    }, numeric(1))
    max_log_ratio[active_rows] <- max_term[active_rows] - median_terms
  }

  attr(contributions, "chain_id")                 <- chain_id
  attr(contributions, "expected_chain_ids")       <- expected_chain_ids
  attr(contributions, "contribution_rows")        <- contribution_rows
  attr(contributions, "sampling_population_size") <-
    length(sampling_population_rows)
  attr(contributions, "finite_terms") <- finite_terms
  attr(contributions, "target")       <- y

  sampling_error <- .iwmde_sampling_mcse(contributions)

  return(list(
    y                        = y,
    finite_terms             = finite_terms,
    max_log_ratio            = max_log_ratio,
    ess                      = ess,
    max_weight_share         = max_weight_share,
    contributions            = contributions,
    sampling_mcse            = sampling_error[["mcse"]],
    sampling_relative_mcse   = sampling_error[["relative_mcse"]],
    sampling_fraction        = sampling_error[["sampling_fraction"]],
    sampling_uncertainty_type = "finite_population_srswor"
  ))
}


# This is design uncertainty from selecting a finite simple random sample of
# eligible posterior rows. It is distinct from MCMC uncertainty in that sample.
.iwmde_sampling_mcse <- function(contributions) {

  contribution_rows <- attr(contributions, "contribution_rows", exact = TRUE)
  population_size <- attr(
    contributions,
    "sampling_population_size",
    exact = TRUE
  )
  target <- attr(contributions, "target", exact = TRUE)
  n <- length(contribution_rows)
  if (population_size < n) {
    return(list(
      mcse              = rep(NA_real_, nrow(contributions)),
      relative_mcse     = rep(NA_real_, nrow(contributions)),
      sampling_fraction = if (population_size > 0L) n / population_size else NA_real_
    ))
  }

  sampling_fraction <- if (population_size == 0L) 1 else n / population_size
  if (n == population_size) {
    mcse <- rep(0, nrow(contributions))
  } else if (n < 2L) {
    mcse <- rep(NA_real_, nrow(contributions))
  } else {
    sample_variance <- apply(contributions, 1L, stats::var)
    mcse <- sqrt((1 - sampling_fraction) * sample_variance / n)
  }
  relative_mcse <- rep(NA_real_, length(mcse))
  positive <- is.finite(target) & target > 0
  relative_mcse[positive] <- mcse[positive] / target[positive]

  return(list(
    mcse              = mcse,
    relative_mcse     = relative_mcse,
    sampling_fraction = sampling_fraction
  ))
}


# Batch diagnostics here describe only the selected posterior-row sequence.
# They are not a full-chain MCSE/ESS when a finite row budget is used.
.iwmde_batch_mcse <- function(contributions) {

  n <- ncol(contributions)
  contribution_sum <- rowSums(contributions)
  contribution_sq  <- rowSums(contributions^2)
  kish_ess <- numeric(nrow(contributions))
  positive_sq <- contribution_sq > 0
  kish_ess[positive_sq] <- contribution_sum[positive_sq]^2 /
    contribution_sq[positive_sq]
  chain_id <- attr(contributions, "chain_id", exact = TRUE)
  if (is.null(chain_id)) {
    chain_id <- rep(1L, n)
  }
  expected_chain_ids <- attr(
    contributions,
    "expected_chain_ids",
    exact = TRUE
  )
  if (is.null(expected_chain_ids)) {
    expected_chain_ids <- unique(chain_id)
  }
  if (length(chain_id) != n || anyNA(chain_id) ||
      length(expected_chain_ids) == 0L || anyNA(expected_chain_ids)) {
    stop("Invalid chain IDs for IWMDE contributions.", call. = FALSE)
  }
  missing_chain_ids <- expected_chain_ids[
    is.na(match(expected_chain_ids, unique(chain_id)))
  ]
  if (length(missing_chain_ids) > 0L) {
    return(list(
      mcse          = rep(NA_real_, nrow(contributions)),
      relative_mcse = rep(NA_real_, nrow(contributions)),
      ess           = rep(NA_real_, nrow(contributions)),
      batch_size    = NA_integer_,
      n_batches     = 0L,
      uncertainty_scope = "unavailable_missing_selected_chain",
      uncertainty_status = "unavailable",
      uncertainty_reason = paste0(
        "selected SRS contains no rows from fitted chain(s): ",
        paste(missing_chain_ids, collapse = ", ")
      )
    ))
  }
  if (n < 4L) {
    return(list(
      mcse          = rep(NA_real_, nrow(contributions)),
      relative_mcse = rep(NA_real_, nrow(contributions)),
      ess           = kish_ess,
      batch_size    = NA_integer_,
      n_batches     = 0L,
      uncertainty_scope = "selected_continuous_rows_only",
      uncertainty_status = "partial",
      uncertainty_reason = "fewer than four selected continuous rows"
    ))
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
      n_batches     = n_batches,
      uncertainty_scope = "selected_continuous_rows_only",
      uncertainty_status = "partial",
      uncertainty_reason = "fewer than two selected-row batches"
    ))
  }

  batch_means <- matrix(NA_real_, nrow = nrow(contributions), ncol = n_batches)
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
  target <- attr(contributions, "target", exact = TRUE)
  if (is.null(target)) {
    target <- rowMeans(contributions)
  }
  relative_mcse <- rep(NA_real_, length(mcse))
  positive <- target > 0
  relative_mcse[positive] <- mcse[positive] / target[positive]
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
    n_batches     = n_batches,
    uncertainty_scope = "selected_continuous_rows_only",
    uncertainty_status = "available",
    uncertainty_reason = NULL
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
  expected_chain_ids <- attr(
    contributions,
    "expected_chain_ids",
    exact = TRUE
  )
  if (!is.null(expected_chain_ids)) {
    attr(integrals, "expected_chain_ids") <- expected_chain_ids
  }
  attr(integrals, "target") <- .iwmde_trapz(
    x,
    attr(contributions, "target", exact = TRUE)
  )
  mcse_data <- .iwmde_batch_mcse(integrals)

  return(list(
    mcse          = mcse_data[["mcse"]][[1L]],
    relative_mcse = mcse_data[["relative_mcse"]][[1L]]
  ))
}
