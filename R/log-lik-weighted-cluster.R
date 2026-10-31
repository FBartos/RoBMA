# Exact estimate deletion for the fitted powered Gaussian cluster likelihood.
# Keep this eligibility separate from the known-V covariance backend: selection
# and arbitrary weighted covariance models do not share this Gaussian integral.
.estimate_weighted_normal_cluster <- function(data, priors) {

  identical(.data_outcome_type(data), "norm") && .is_data_multilevel(data) &&
    .is_data_weights(data) && !.is_data_known_v(data) &&
    !.is_data_random(data) && !.is_priors_weightfunction(priors)
}

.estimate_normal_target_integrates_local_effects <- function(data, priors) {

  .estimate_normal_target_uses_covariance_backend(data, priors) ||
    .estimate_weighted_normal_cluster(data, priors)
}

.weighted_cluster_numeric_unavailable <- function() {

  stop(structure(list(
    message = paste0("Weighted Gaussian cluster deletion is unavailable at the supplied numerical scales. ",
                     "The inputs exceed reliable floating-point precision. ",
                     "Inspect the outcome scales and likelihood weights."),
    call = NULL
  ), class = c("RoBMA_weighted_gaussian_unavailable", "RoBMA_numerical_unavailable", "error", "condition")))
}

# Prefix/suffix information sums never subtract the deleted contribution from
# a dominating total. Both sums are linear and vectorized over posterior draws.
.weighted_normal_cluster_deletion <- function(residual, sd, loading, weights, cluster) {

  residual <- as.matrix(residual)
  sd <- as.matrix(sd)
  loading <- as.matrix(loading)
  S <- nrow(residual)
  K <- ncol(residual)
  if(!identical(dim(sd), dim(residual)) || !identical(dim(loading), dim(residual)) ||
     length(weights) != K || any(!is.finite(weights) | weights <= 0) ||
     any(!is.finite(sd) | sd <= 0) || any(!is.finite(residual)) || any(!is.finite(loading))) {
    .weighted_cluster_numeric_unavailable()
  }
  rows <- unlist(cluster, use.names = FALSE)
  if(!identical(sort(as.integer(rows)), seq_len(K))) {
    stop("Cluster metadata must assign every estimate to one cluster.", call. = FALSE)
  }
  variance <- sd^2
  # Do not support representability-limit inputs with exact-arithmetic machinery.
  # Squared SDs must be normal doubles; powers above 1/eps and residuals above
  # 1/sqrt(eps) observation SDs cannot retain unit-scale information reliably.
  if(any(!is.finite(variance) | variance < .Machine$double.xmin) ||
     any(weights > 1 / .Machine$double.eps) ||
     any(abs(residual / sd) > 1 / sqrt(.Machine$double.eps))) {
    .weighted_cluster_numeric_unavailable()
  }
  scaled_a <- loading / sd
  information <- sweep(scaled_a^2, 2L, weights, `*`)
  numerator <- sweep(scaled_a * (residual / sd), 2L, weights, `*`)
  if(any(!is.finite(information)) || any(!is.finite(numerator))) {
    .weighted_cluster_numeric_unavailable()
  }
  latent_variance <- shift <- matrix(0, S, K)
  error <- residual
  for(ids in cluster){
    n <- length(ids)
    prefix_p <- matrix(1, S, n + 1L) # unit-normal prior included exactly once
    prefix_num <- matrix(0, S, n + 1L)
    for(j in seq_len(n)){
      i <- ids[j]
      prefix_p[, j + 1L] <- prefix_p[, j] + information[, i]
      prefix_num[, j + 1L] <- prefix_num[, j] + numerator[, i]
    }
    suffix_p <- suffix_num <- numeric(S)
    for(j in rev(seq_len(n))){
      i <- ids[j]
      p <- prefix_p[, j] + suffix_p
      num <- prefix_num[, j] + suffix_num
      if(any(!is.finite(p)) || any(!is.finite(num))){
        .weighted_cluster_numeric_unavailable()
      }
      current_shift <- loading[, i] * (num / p)
      latent_variance[, i] <- (loading[, i] / sqrt(p))^2
      shift[, i] <- current_shift
      error[, i] <- residual[, i] - current_shift
      suffix_p <- suffix_p + information[, i]
      suffix_num <- suffix_num + numerator[, i]
    }
  }
  if(any(!is.finite(error)) || any(!is.finite(shift))) .weighted_cluster_numeric_unavailable()
  # Materialize variances only in the ordinary domain. log1p preserves even
  # arbitrarily small positive powers in the determinant penalty.
  powered_variance <- variance + sweep(latent_variance, 2L, weights, `*`)
  observation_variance <- variance + latent_variance
  log_sd <- log(sd)
  log_penalty <- .5 * log1p(sweep(latent_variance / variance, 2L, weights, `*`))
  log_score_sd <- .5 * log(powered_variance)
  scaled_square <- sweep((error / sqrt(powered_variance))^2, 2L, weights, `*`)
  log_score <- -.5 * sweep(log(2 * pi) + 2 * log_sd, 2L, weights, `*`) -
    log_penalty - .5 * scaled_square
  if(any(!is.finite(log_score)) || any(!is.finite(powered_variance) | powered_variance <= 0) ||
     any(!is.finite(observation_variance) | observation_variance <= 0) ||
     any(latent_variance == 0 & loading != 0)) {
    .weighted_cluster_numeric_unavailable()
  }
  list(log_score = log_score, residual = error, shift = shift,
       log_score_sd = log_score_sd, log_sampling_sd = log_sd, log_penalty = log_penalty,
       log_sd = .5 * log(observation_variance))
}

.weighted_normal_cluster_summary_from_setup <- function(setup,
    components = c("log_density", "cdf", "log_lower", "log_upper", "mean", "variance")) {

  location <- .estimate_normal_covariance_target_location_from_setup(setup)
  sd <- .root_sum_squares(matrix(setup[["sei"]], setup[["S"]], setup[["K"]], byrow = TRUE),
                          setup[["tau_within"]])
  y <- matrix(location$y, setup[["S"]], setup[["K"]], byrow = TRUE)
  # Refuse large opposing locations before their subtraction loses observation-
  # scale information. No second deletion pass on rounding remainders is needed.
  if(any(!is.finite(y)) || any(!is.finite(location$means)) ||
     any(!is.finite(sd) | sd <= 0) ||
     any(pmax(abs(y), abs(location$means)) / sd > 1 / sqrt(.Machine$double.eps))) {
    .weighted_cluster_numeric_unavailable()
  }
  state <- .weighted_normal_cluster_deletion(y - location$means, sd,
    setup[["tau_between"]], setup[["weights"]], setup[["cluster"]])
  residual <- state$residual
  out <- list()
  if("log_density" %in% components) out[["log_density"]] <- state$log_score
  observation_sd <- if(any(components %in% c("cdf", "log_lower", "log_upper"))) exp(state$log_sd) else NULL
  if(any(!is.finite(observation_sd) | observation_sd <= 0)) .weighted_cluster_numeric_unavailable()
  for(component in intersect(c("cdf", "log_lower", "log_upper"), components)){
    out[[component]] <- stats::pnorm(residual, sd = observation_sd,
      lower.tail = if(component == "log_upper") !location$lower_tail else location$lower_tail,
      log.p = component != "cdf")
  }
  if("mean" %in% components){
    mean <- location$means + state$shift
    use_residual <- pmax(abs(y), abs(residual)) < pmax(abs(location$means), abs(state$shift))
    mean[use_residual] <- (y - residual)[use_residual]
    if(!location$lower_tail) mean <- -mean
    if(any(!is.finite(mean))) .weighted_cluster_numeric_unavailable()
    out[["mean"]] <- mean
  }
  if("variance" %in% components) out[["variance"]] <- exp(2 * state$log_sd)
  attr(out, "retained_residual") <- if(location$lower_tail) residual else -residual
  attr(out, "dependency_blocks") <- unname(setup[["cluster"]])
  out
}
