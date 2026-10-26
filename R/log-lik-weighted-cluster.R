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

.weighted_cluster_log_add <- function(x, y) {

  top <- pmax(x, y)
  out <- top
  finite <- is.finite(top)
  out[finite] <- top[finite] + log1p(exp(pmin(x[finite], y[finite]) - top[finite]))
  out
}

# Signed log subtraction, including exact zero. Its inputs are logarithms of
# positive/negative parts, so neither large precision nor a tiny latent variance
# needs to be materialized.
.weighted_cluster_log_difference <- function(positive, negative) {

  top <- pmax(positive, negative)
  sign <- sign(positive - negative)
  sign[positive == negative] <- 0
  value <- rep(-Inf, length(top))
  nonzero <- sign != 0
  value[nonzero] <- top[nonzero] +
    log(-expm1(pmin(positive[nonzero], negative[nonzero]) - top[nonzero]))
  list(log_abs = value, sign = sign)
}

.weighted_cluster_numeric_unavailable <- function() {

  stop(structure(list(
    message = paste0("Weighted Gaussian cluster deletion is unavailable at the supplied numerical scales. ",
                     "Inspect the outcome scales and likelihood weights."),
    call = NULL
  ), class = c("RoBMA_weighted_gaussian_unavailable", "RoBMA_numerical_unavailable", "error", "condition")))
}

# Compensated signed sum for rare cancellation cases. Normalization prevents
# overflow; a lost nonzero summand is refused rather than silently erased.
.weighted_cluster_compensated_sum <- function(log_abs, signs) {

  keep <- signs != 0 & is.finite(log_abs)
  if(!any(keep)) return(list(log_abs = -Inf, sign = 0))
  top <- max(log_abs[keep])
  values <- signs[keep] * exp(log_abs[keep] - top)
  if(any(values == 0)) .weighted_cluster_numeric_unavailable()
  total <- correction <- 0
  for(value in values){
    next_total <- total + value
    correction <- correction + if(abs(total) >= abs(value)) {
      (total - next_total) + value
    } else (value - next_total) + total
    total <- next_total
  }
  total <- total + correction
  list(log_abs = if(total == 0) -Inf else top + log(abs(total)), sign = sign(total))
}

# Compensated 2x2 determinant of bounded, power-of-two-scaled inputs. Retain
# the multiplication roundoff when two products nearly cancel. Subnormal
# products have no error-free split in double precision and are refused here.
.weighted_cluster_product <- function(a, b) {
    value <- a * b
    if(a != 0 && b != 0 && abs(value) < .Machine$double.xmin) {
      .weighted_cluster_numeric_unavailable()
    }
    split_a <- 134217729 * a
    high_a <- split_a - (split_a - a)
    low_a <- a - high_a
    split_b <- 134217729 * b
    high_b <- split_b - (split_b - b)
    low_b <- b - high_b
    error <- ((high_a * high_b - value) + high_a * low_b + low_a * high_b) + low_a * low_b
    c(value, error)
}

.weighted_cluster_pair_difference <- function(x, y, u, v) {

  first <- .weighted_cluster_product(x, y)
  second <- .weighted_cluster_product(u, v)
  difference <- first[1L] - second[1L]
  part <- difference - first[1L]
  error <- (first[1L] - (difference - part)) - (second[1L] + part)
  difference + (error + first[2L] - second[2L])
}

# Rare signed cancellation is evaluated from original operands, never by
# exponentiating rounded logarithms. Binary SD rescalings admit compensated
# products; other cases need a small forward-error bound on the focal shift.
.weighted_cluster_cancelled_numerator <- function(r, sd, a, weights,
                                                  log_precision, focal_a,
                                                  focal_sd, focal_weight) {

  use <- r != 0 & a != 0
  if(!any(use)) return(list(log_abs = -Inf, sign = 0))
  r <- r[use]; sd <- sd[use]; a <- a[use]; weights <- weights[use]
  powers <- round(log2(sd) - log2(sd[1L]))
  binary <- all(abs(powers) <= 511) &&
    (all(sd == sd[1L]) || (all(sd >= .Machine$double.xmin) && all(sd == sd[1L] * 2^powers)))
  if(binary){
    power_a <- min(.Machine$double.max.exp - 1, floor(log2(max(abs(a)))))
    power_r <- min(.Machine$double.max.exp - 1, floor(log2(max(abs(r)))))
    power_w <- min(.Machine$double.max.exp - 1, floor(log2(max(weights))))
    aa <- a / 2^power_a
    rr <- r / 2^power_r
    ww <- weights / 2^power_w
    if(any(aa == 0) || any(rr == 0) || any(ww == 0)) .weighted_cluster_numeric_unavailable()
    top <- max(-2 * powers)
    values <- unlist(lapply(seq_along(r), function(j){
      pair <- .weighted_cluster_product(aa[j], rr[j])
      product <- unlist(lapply(pair, .weighted_cluster_product, b = ww[j]), use.names = FALSE)
      value <- product * 2^(-2 * powers[j] - top)
      if(any(value == 0 & product != 0)) .weighted_cluster_numeric_unavailable()
      value
    }), use.names = FALSE)
    log_scale <- (power_a + power_r + power_w + top) * log(2) - 2 * log(sd[1L])
  }else{
    terms <- weights * (a / sd) * (r / sd)
    if(any(!is.finite(terms) | terms == 0)) .weighted_cluster_numeric_unavailable()
    scale <- max(abs(terms))
    values <- terms / scale
    log_scale <- log(scale)
    log_bound <- log(8 * .Machine$double.eps) + log_scale + log(sum(abs(values)))
    log_latent <- log(abs(focal_a)) - .5 * log_precision
    log_observed <- .5 * .weighted_cluster_log_add(2 * log(focal_sd), 2 * log_latent)
    log_powered <- .5 * .weighted_cluster_log_add(2 * log(focal_sd) - log(focal_weight), 2 * log_latent)
    if(log(abs(focal_a)) - log_precision + log_bound >
       log(64 * .Machine$double.eps) + min(log_observed, log_powered)) {
      .weighted_cluster_numeric_unavailable()
    }
  }
  total <- correction <- 0
  for(value in values){
    next_total <- total + value
    correction <- correction + if(abs(total) >= abs(value)) {
      (total - next_total) + value
    } else (value - next_total) + total
    total <- next_total
  }
  total <- total + correction
  list(log_abs = if(total == 0) -Inf else log(abs(total)) + log_scale, sign = sign(total))
}

# Prefix/suffix information sums: the deleted row is never subtracted from a
# potentially dominating total. Ordinary work is O(S*K), vectorized over draws.
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
  log_sd <- log(sd)
  log_a <- log(abs(loading))
  log_w <- log(weights)
  information <- sweep(2 * (log_a - log_sd), 2L, log_w, `+`)
  numerator <- sweep(log_a + log(abs(residual)) - 2 * log_sd, 2L, log_w, `+`)
  signs <- sign(loading) * sign(residual)
  positive <- negative <- numerator
  positive[signs <= 0] <- -Inf
  negative[signs >= 0] <- -Inf
  log_p <- shift <- matrix(0, S, K)
  error <- residual

  for(ids in cluster){
    n <- length(ids)
    prefix_p <- prefix_pos <- prefix_neg <- matrix(-Inf, S, n + 1L)
    prefix_p[, 1L] <- 0 # unit-normal cluster prior, included exactly once
    for(j in seq_len(n)){
      i <- ids[j]
      prefix_p[, j + 1L] <- .weighted_cluster_log_add(prefix_p[, j], information[, i])
      prefix_pos[, j + 1L] <- .weighted_cluster_log_add(prefix_pos[, j], positive[, i])
      prefix_neg[, j + 1L] <- .weighted_cluster_log_add(prefix_neg[, j], negative[, i])
    }
    suffix_p <- suffix_pos <- suffix_neg <- rep(-Inf, S)
    for(j in rev(seq_len(n))){
      i <- ids[j]
      lp <- .weighted_cluster_log_add(prefix_p[, j], suffix_p)
      pos <- .weighted_cluster_log_add(prefix_pos[, j], suffix_pos)
      neg <- .weighted_cluster_log_add(prefix_neg[, j], suffix_neg)
      num <- .weighted_cluster_log_difference(pos, neg)
      cancellation <- is.finite(pmax(pos, neg)) &
        num$log_abs <= pmax(pos, neg) + log(32 * .Machine$double.eps)
      for(s in which(cancellation)){
        other <- ids[-j]
        exact <- .weighted_cluster_cancelled_numerator(residual[s, other], sd[s, other],
          loading[s, other], weights[other], lp[s], loading[s, i], sd[s, i], weights[i])
        num$log_abs[s] <- exact$log_abs
        num$sign[s] <- exact$sign
      }
      log_shift <- log_a[, i] + num$log_abs - lp
      log_shift[loading[, i] == 0 | num$sign == 0] <- -Inf
      shift[, i] <- sign(loading[, i]) * num$sign * exp(log_shift)
      error[, i] <- residual[, i] - shift[, i]
      # When the data align with a large cluster effect, form the retained
      # pair differences directly and keep the small prior contribution.
      roundoff <- 16 * .Machine$double.eps * pmax(1, abs(lp),
        ifelse(is.finite(log_shift), abs(log_shift), 0), abs(log_sd[, i]))
      close <- is.finite(shift[, i]) & loading[, i] != 0 & num$sign != 0 &
        abs(error[, i]) <= roundoff * pmax(abs(residual[, i]), abs(shift[, i]))
      for(s in which(close)){
        other <- ids[-j]
        terms <- log(abs(residual[s, i])) - lp[s]
        term_sign <- sign(residual[s, i])
        for(k in other){
          scale_a <- max(abs(loading[s, c(i, k)]))
          scale_r <- max(abs(residual[s, c(i, k)]))
          if(scale_r == 0 || loading[s, k] == 0) next
          scale_a <- 2^min(.Machine$double.max.exp - 1, floor(log2(scale_a)))
          scale_r <- 2^min(.Machine$double.max.exp - 1, floor(log2(scale_r)))
          a <- loading[s, c(i, k)] / scale_a
          r <- residual[s, c(i, k)] / scale_r
          if(any(a == 0 & loading[s, c(i, k)] != 0) ||
             any(r == 0 & residual[s, c(i, k)] != 0)) .weighted_cluster_numeric_unavailable()
          difference <- .weighted_cluster_pair_difference(a[2L], r[1L], a[1L], r[2L])
          if(difference == 0) next
          terms <- c(terms, log_w[k] + log_a[s, k] - 2 * log_sd[s, k] - lp[s] +
                       log(scale_a) + log(scale_r) + log(abs(difference)))
          term_sign <- c(term_sign, sign(loading[s, k]) * sign(difference))
        }
        pair <- .weighted_cluster_compensated_sum(terms, term_sign)
        error[s, i] <- pair$sign * exp(pair$log_abs)
      }
      log_p[, i] <- lp
      suffix_p <- .weighted_cluster_log_add(suffix_p, information[, i])
      suffix_pos <- .weighted_cluster_log_add(suffix_pos, positive[, i])
      suffix_neg <- .weighted_cluster_log_add(suffix_neg, negative[, i])
    }
  }
  if(any(!is.finite(error)) || any(!is.finite(shift))) .weighted_cluster_numeric_unavailable()
  log_latent_sd <- log_a - .5 * log_p
  log_variance_ratio <- sweep(2 * (log_latent_sd - log_sd), 2L, log_w, `+`)
  log_penalty <- .5 * .weighted_cluster_log_add(matrix(0, S, K), log_variance_ratio)
  log_score_sd <- .5 * .weighted_cluster_log_add(2 * log_sd,
    sweep(2 * log_latent_sd, 2L, log_w, `+`))
  scaled_square <- exp(sweep(2 * (log(abs(error)) - log_score_sd), 2L, log_w, `+`))
  log_score <- -.5 * sweep(log(2 * pi) + 2 * log_sd, 2L, weights, `*`) -
    log_penalty - .5 * scaled_square
  if(anyNA(log_score)) .weighted_cluster_numeric_unavailable()
  list(log_score = log_score, residual = error, shift = shift,
       log_score_sd = log_score_sd, log_sampling_sd = log_sd, log_penalty = log_penalty,
       log_sd = .5 * .weighted_cluster_log_add(2 * log_sd, 2 * log_latent_sd))
}

.weighted_normal_cluster_summary_from_setup <- function(setup,
    components = c("log_density", "cdf", "log_lower", "log_upper", "mean", "variance")) {

  location <- .estimate_normal_covariance_target_location_from_setup(setup)
  sd <- .root_sum_squares(matrix(setup[["sei"]], setup[["S"]], setup[["K"]], byrow = TRUE),
                          setup[["tau_within"]])
  y <- matrix(location$y, setup[["S"]], setup[["K"]], byrow = TRUE)
  # Carry the rounding remainder of y-mu through the linear deletion operator.
  # Operating on two huge opposing inputs instead would lose small differences
  # that subtraction represents exactly (Sterbenz's lemma).
  residual_high <- y - location$means
  if(any(!is.finite(residual_high))) .weighted_cluster_numeric_unavailable()
  part <- residual_high - y
  residual_low <- (y - (residual_high - part)) + (-location$means - part)
  state <- .weighted_normal_cluster_deletion(residual_high, sd,
    setup[["tau_between"]], setup[["weights"]], setup[["cluster"]])
  low_state <- if(any(residual_low != 0)) .weighted_normal_cluster_deletion(residual_low, sd,
    setup[["tau_between"]], setup[["weights"]], setup[["cluster"]]) else
      list(residual = matrix(0, setup$S, setup$K), shift = matrix(0, setup$S, setup$K))
  residual <- state$residual + low_state$residual
  if(any(!is.finite(residual))) .weighted_cluster_numeric_unavailable()
  scaled_square <- exp(sweep(2 * (log(abs(residual)) - state$log_score_sd),
                            2L, log(setup$weights), `+`))
  log_score <- -.5 * sweep(log(2 * pi) + 2 * state$log_sampling_sd,
                           2L, setup$weights, `*`) -
    state$log_penalty - .5 * scaled_square
  if(anyNA(log_score)) .weighted_cluster_numeric_unavailable()
  out <- list()
  if("log_density" %in% components) out[["log_density"]] <- log_score
  observation_sd <- if(any(components %in% c("cdf", "log_lower", "log_upper"))) exp(state$log_sd) else NULL
  if(any(!is.finite(observation_sd) | observation_sd <= 0)) .weighted_cluster_numeric_unavailable()
  for(component in intersect(c("cdf", "log_lower", "log_upper"), components)){
    out[[component]] <- stats::pnorm(residual, sd = observation_sd,
      lower.tail = if(component == "log_upper") !location$lower_tail else location$lower_tail,
      log.p = component != "cdf")
  }
  if("mean" %in% components){
    mean <- (location$means + state$shift) + low_state$shift
    residual_scale <- pmax(abs(y), abs(residual))
    direct_scale <- pmax(abs(location$means), abs(state$shift), abs(low_state$shift))
    use_residual <- residual_scale < direct_scale
    mean[use_residual] <- (y - residual)[use_residual]
    if(!location$lower_tail) mean <- -mean
    if(any(!is.finite(mean))) .weighted_cluster_numeric_unavailable()
    out[["mean"]] <- mean
  }
  if("variance" %in% components){
    variance <- exp(2 * state$log_sd)
    if(any(!is.finite(variance) | variance <= 0)) .weighted_cluster_numeric_unavailable()
    out[["variance"]] <- variance
  }
  attr(out, "retained_residual") <- if(location$lower_tail) residual else -residual
  attr(out, "dependency_blocks") <- unname(setup[["cluster"]])
  out
}
