#' @title Weighted normal distribution
#'
#' @description Density, distribution function, quantile function
#' and random generation for the weighted normal distribution with
#' \code{mean}, standard deviation \code{sd}, steps \code{steps}
#' (or critical values) \code{crit_x}), and weights \code{omega}.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length
#' is taken to be the number required.
#' @param mean mean
#' @param sd standard deviation.
#' @param steps vector of steps for the weight function.
#' @param omega vector of weights defining the probability
#' of observing a t-statistics between each of the two steps.
#' @param crit_x vector of critical values defining steps
#' (if \code{steps} are not supplied).
#' @param type type of weight function (defaults to \code{"two.sided"}).
#' @param log,log.p logical; if \code{TRUE}, probabilities
#' \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities
#' are \eqn{P[X \le x]}, otherwise, \eqn{P[X \ge x]}.
#'
#' @details The \code{mean}, \code{sd}, \code{steps}, \code{omega} can be
#' supplied as a vectors (\code{mean}, \code{sd}) or matrices (\code{steps},
#' \code{omega}) with length / number of rows equal to \code{x}/\code{q}/
#' \code{p}. Otherwise, they are recycled to the length of the result.
#'
#'
#' @return \code{dwnorm} gives the density, \code{dwnorm} gives the
#' distribution function, \code{qwnorm} gives the quantile function,
#' and \code{rwnorm} generates random deviates.
#'
#' @export dwnorm
#' @export pwnorm
#' @export qwnorm
#' @export rwnorm
#' @name weighted_normal
#'
#' @seealso \link[stats]{Normal}
#' @examples
#' # generate random samples from weighted normal distribution
#' samples <- rwnorm(n = 10000, mean = 0.15, sd = 0.10,
#'                   steps = c(0.025, 0.5), omega = c(0.1, 0.5, 1),
#'                   type = "one.sided")
#' # hist(samples)
#'
#' # compute density at specific values
#' density_vals <- dwnorm(x = c(-2, 0, 2), mean = 0, sd = 1,
#'                        steps = c(0.05), omega = c(0.5, 1))
#'
#' # compute cumulative probabilities
#' prob_vals <- pwnorm(q = c(-1, 0, 1), mean = 0, sd = 1,
#'                     steps = c(0.05), omega = c(0.5, 1))
NULL

#' @rdname weighted_normal
dwnorm <- function(x, mean, sd, steps = if(!is.null(crit_x)) NULL, omega, crit_x = if(!is.null(steps)) NULL, type = "two.sided", log = FALSE){

  # set and check parameters
  BayesTools::check_real(x, "x", check_length = FALSE)
  BayesTools::check_bool(log, "log")
  .Xwnorm_check_input(mean, sd, omega, steps, crit_x, type)

  # repeat mean/sd to match the input length
  if(length(mean) != length(x)){
    mean  <- rep_len(mean,  length(x))
  }
  if(length(sd) != length(x)){
    sd    <- rep_len(sd,  length(x))
  }
  if(!is.matrix(omega)){
    omega <- matrix(omega, ncol = length(omega), nrow = length(x), byrow = TRUE)
  }

  # check that the omega and steps/crit_x dimensions match
  crit_x <- .Xwnorm_get_crit_x(steps, crit_x, mean, sd, type)
  if(ncol(crit_x) + 1 != ncol(omega))
    stop("'omega' argument must have one more weight than the number of defined steps with 'steps'/'crit_x' argument.")



  # do the actual computations
  if(type == "two.sided"){

    # assign appropriate weights to the current values
    w <- sapply(1:length(x), function(i) .get_log_weight_twosided(x[i], crit_x = crit_x[i, ], omega = omega[i, ]))

    # compute the nominator
    nom <- stats::dnorm(x, mean, sd, log = TRUE) + w

    # compute the denominator
    denoms  <- matrix(stats::pnorm(crit_x[,1], mean, sd) - stats::pnorm(-crit_x[,1], mean, sd), ncol = 1)
    # check and correct for possibly negative numbers due to numerical imprecision
    denoms[denoms < 0, 1]  <- 0
    if(ncol(omega) > 2){
      for(j in 2:(ncol(omega)-1)){
        denoms <- cbind(denoms, stats::pnorm(crit_x[,j], mean, sd) - stats::pnorm(-crit_x[,j], mean, sd) - apply(denoms, 1, sum))
        # check and correct for possibly negative numbers due to numerical imprecision
        denoms[denoms[,j] < 0, j] <- 0
      }
    }
    denoms  <- cbind(denoms, 1 - apply(denoms, 1, sum))
    denoms  <- log(denoms) + log(omega)
    denom   <- log(apply(exp(denoms), 1, sum))

    log_lik <- nom - denom

  }else if(type == "one.sided"){

    # assign appropriate weights to the current values
    w <- sapply(1:length(x), function(i) .get_log_weight_onesided(x[i], crit_x = crit_x[i, ], omega = omega[i, ]))

    # compute the nominator
    nom <- stats::dnorm(x, mean, sd, log = TRUE) + w

    denoms  <- matrix(stats::pnorm(crit_x[,1], mean, sd), ncol = 1)
    # check and correct for possibly negative numbers due to numerical imprecision
    denoms[denoms < 0, 1]  <- 0
    if(ncol(omega) > 2){
      for(j in 2:(ncol(omega)-1)){
        denoms <- cbind(denoms, stats::pnorm(crit_x[,j], mean, sd) - apply(denoms, 1, sum))
        # check and correct for possibly negative numbers due to numerical imprecision
        denoms[denoms[,j] < 0, j] <- 0
      }
    }
    denoms  <- cbind(denoms, 1 - apply(denoms, 1, sum))
    denoms  <- log(denoms) + log(omega)
    denom   <- log(apply(exp(denoms), 1, sum))

    log_lik <- nom - denom
  }


  if(log){
    return(log_lik)
  }else{
    return(exp(log_lik))
  }
}
#' @rdname weighted_normal
pwnorm <- function(q, mean, sd, steps = if(!is.null(crit_x)) NULL, omega, crit_x = if(!is.null(steps)) NULL, type = "two.sided", lower.tail = TRUE, log.p = FALSE){

  # set and check parameters
  BayesTools::check_real(q, "q", check_length = FALSE)
  BayesTools::check_bool(lower.tail, "lower.tail")
  BayesTools::check_bool(log.p, "log.p")
  .Xwnorm_check_input(mean, sd, omega, steps, crit_x, type)

  # repeat mean/sd to match the input length
  if(length(mean) != length(q)){
    mean  <- rep_len(mean,  length(q))
  }
  if(length(sd) != length(q)){
    sd    <- rep_len(sd,  length(q))
  }
  if(!is.matrix(omega)){
    omega <- matrix(omega, ncol = length(omega), nrow = length(q), byrow = TRUE)
  }

  # check that the omega and steps/crit_x dimensions match
  crit_x <- .Xwnorm_get_crit_x(steps, crit_x, mean, sd, type)
  if(ncol(crit_x) + 1 != ncol(omega))
    stop("'omega' argument must have one more weight than the number of defined steps with 'steps'/'crit_x' argument.")


  # work with two-sided as with one-sided
  if(type == "two.sided"){
    crit_x <- cbind(matrix(-crit_x[,ncol(crit_x):1], ncol = ncol(crit_x)), 0,  crit_x)
    omega  <- cbind(matrix( omega[,ncol(omega):1],   ncol = ncol(omega)),      omega)
  }

  p <- .pwnorm(q, mean, sd, omega, crit_x, type, lower.tail)

  if(log.p){
    return(p)
  }else{
    return(exp(p))
  }
}
#' @rdname weighted_normal
qwnorm <- function(p, mean, sd, steps = if(!is.null(crit_x)) NULL, omega, crit_x = if(!is.null(steps)) NULL, type = "two.sided", lower.tail = TRUE, log.p = FALSE){

  # set and check parameters
  BayesTools::check_bool(log.p, "log.p")
  BayesTools::check_real(p, "p", check_length = FALSE, lower = if(log.p) -Inf else 0, upper = if(log.p) 0 else 1)
  BayesTools::check_bool(lower.tail, "lower.tail")
  .Xwnorm_check_input(mean, sd, omega, steps, crit_x, type)

  # make p normal again
  if(log.p){
    p <- exp(p)
  }

  # repeat mean/sd to match the input length
  if(length(mean) != length(p)){
    mean  <- rep_len(mean,  length(p))
  }
  if(length(sd) != length(p)){
    sd    <- rep_len(sd,  length(p))
  }
  if(!is.matrix(omega)){
    omega <- matrix(omega, ncol = length(omega), nrow = length(p), byrow = TRUE)
  }

  # check that the omega and steps/crit_x dimensions match
  crit_x <- .Xwnorm_get_crit_x(steps, crit_x, mean, sd, type)
  if(ncol(crit_x) + 1 != ncol(omega))
    stop("'omega' argument must have one more weight than the number of defined steps with 'steps'/'crit_x' argument.")



  # work with two-sided as with one-sided
  if(type == "two.sided"){
    crit_x <- cbind(matrix(-crit_x[,ncol(crit_x):1], ncol = ncol(crit_x)), 0,  crit_x)
    omega  <- cbind(matrix( omega[,ncol(omega):1],   ncol = ncol(omega)),      omega)
  }


  # yeah, this is slow and inefficient, but I was too lazy to write it in a better way for now
  q <- rep(NA, length(p))
  for(i in 1:length(p)){

    # simple solution for 0 and 1
    if(p[i] == 0){
      q[i] <- ifelse(lower.tail, -Inf, Inf)
      next
    }else if(p[i] == 1){
      q[i] <- ifelse(lower.tail, Inf, -Inf)
      next
    }

    # much more difficult solution for everything else
    q[i] <- .qwnorm(p[i], mean[i], sd[i], omega[i,], crit_x[i,], type, lower.tail)
  }

  return(q)
}
#' @rdname weighted_normal
rwnorm <- function(n, mean, sd, steps = if(!is.null(crit_x)) NULL, omega, crit_x = if(!is.null(steps)) NULL, type = "two.sided"){

  # set and check parameters
  BayesTools::check_int(n, "n")
  .Xwnorm_check_input(mean, sd, omega, steps, crit_x, type)

  # for consistency with other functions
  if(n == 0)
    return(numeric())

  # repeat mean/sd to match the input length
  if(length(mean) != n){
    mean  <- rep_len(mean, n)
  }
  if(length(sd) != n){
    sd    <- rep_len(sd,  n)
  }
  if(!is.matrix(omega)){
    omega <- matrix(omega, ncol = length(omega), nrow = n, byrow = TRUE)
  }

  # check that the omega and steps/crit_x dimensions match
  crit_x <- .Xwnorm_get_crit_x(steps, crit_x, mean, sd, type)
  if(ncol(crit_x) + 1 != ncol(omega))
    stop("'omega' argument must have one more weight than the number of defined steps with 'steps'/'crit_x' argument.")


  # work with two-sided as with one-sided
  if(type == "two.sided"){
    crit_x <- cbind(matrix(-crit_x[,ncol(crit_x):1], ncol = ncol(crit_x)), 0,  crit_x)
    omega  <- cbind(matrix( omega[,ncol(omega):1],   ncol = ncol(omega)),      omega)
  }

  # yeah, this is slow and inefficient, but I was too lazy to write it in a better way for now
  x <- rep(NA, n)
  for(i in 1:n){
    while(is.na(x[i])){
      # sample a normal observation
      temp_x <- stats::rnorm(1, mean = mean[i], sd = sd[i])

      # find it's weight
      temp_omega <- exp(.get_log_weight_onesided(temp_x, crit_x = crit_x[i, ], omega = omega[i, ]))

      if(stats::rbinom(1, 1, prob = temp_omega) == 1){
        x[i] <- temp_x
      }
    }
  }


  return(x)
}

# helper functions
.pwnorm <- function(q, mean, sd, omega, crit_x, type, lower.tail){

  p  <- rep(NA, length(q))

  if(lower.tail){
    for(i in 1:length(q)){

      # simple solution for +- Inf
      if(is.infinite(q[i])){
        if(q[i] < 0){
          p[i] <- 0
        }else if(q[i] > 0){
          p[i] <- 1
        }
        next
      }

      # much more difficult solution for everything else
      s      <- 1
      temp_p <- NULL
      # do jumps by every cut-off
      while(q[i] > crit_x[i,s]){

        temp_p_add <- stats::pnorm(crit_x[i,s], mean[i], sd[i]) - sum(temp_p)
        # check and correct for possibly negative numbers due to numerical imprecision
        if(temp_p_add < 0){
          temp_p_add <- 0
        }
        temp_p <- c(temp_p, temp_p_add)

        s <- s + 1
        if(s > ncol(crit_x))
          break
      }

      # add the last piece
      temp_p_add <- stats::pnorm(q[i], mean[i], sd[i]) - sum(temp_p)
      # check and correct for possibly negative numbers due to numerical imprecision
      if(temp_p_add < 0){
        temp_p_add <- 0
      }
      temp_p <- c(temp_p, temp_p_add)

      # weight by the weights
      p[i] <- log(sum(exp(log(temp_p) + log(omega[i,1:s]))))
    }
  }else{
    for(i in 1:length(q)){

      # simple solution for +- Inf
      if(is.infinite(q[i])){
        if(q[i] < 0){
          p[i] <- 1
        }else if(q[i] > 0){
          p[i] <- 0
        }
        next
      }

      # much more difficult solution for everything else
      s      <- ncol(omega)
      temp_p <- NULL
      # do jumps by every cut-off
      while(q[i] < crit_x[i,s]){

        temp_p_add <- stats::pnorm(crit_x[i,s], mean[i], sd[i], lower.tail = FALSE) - sum(temp_p)
        # check and correct for possibly negative numbers due to numerical imprecision
        if(temp_p_add < 0){
          temp_p_add <- 0
        }
        temp_p <- c(temp_p, temp_p_add)

        s <- s + 1
        if(s > ncol(crit_x))
          break
      }

      # add the last piece
      temp_p_add <- stats::pnorm(q[i], mean[i], sd[i], lower.tail = FALSE) - sum(temp_p)
      # check and correct for possibly negative numbers due to numerical imprecision
      if(temp_p_add < 0){
        temp_p_add <- 0
      }
      temp_p <- c(temp_p, temp_p_add)

      # weight by the weights
      p[i] <- log(sum(exp(log(temp_p) + log(omega[i,s:ncol(omega)]))))
    }
  }

  # standardize
  denoms  <- matrix(stats::pnorm(crit_x[,1], mean, sd), ncol = 1)
  # check and correct for possibly negative numbers due to numerical imprecision
  denoms[denoms < 0, 1]  <- 0
  if(ncol(omega) > 2){
    for(j in 2:(ncol(omega)-1)){
      denoms <- cbind(denoms, stats::pnorm(crit_x[,j], mean, sd) - apply(denoms, 1, sum))
      # check and correct for possibly negative numbers due to numerical imprecision
      denoms[denoms[,j] < 0, j] <- 0
    }
  }
  denoms  <- cbind(denoms, 1 - apply(denoms, 1, sum))
  denoms  <- log(denoms) + log(omega)
  denom   <- log(apply(exp(denoms), 1, sum))

  p <- p - denom

  return(p)
}
.qwnorm <- function(p, mean, sd, omega, crit_x, type, lower.tail){
  return(stats::optim(
    par = 0,
    fn  = function(x, p, mean, sd, omega, crit_x, type, lower.tail){
      (exp(.pwnorm(x, mean, sd, omega, crit_x, type, lower.tail)) - p)^2
    },p = p, mean = mean, sd = sd, omega = matrix(omega, nrow = 1), crit_x = matrix(crit_x, nrow = 1), type = type, lower.tail = lower.tail,
    method = "L-BFGS-B",
    lower = -Inf, upper = Inf,
    control = list(factr = 1e-12))$par)
}

# fast computation - no input check, pre-formatted for bridge-sampling
.dwnorm_fast <- function(x, mean, sd, omega, crit_x, type = "two.sided", log = TRUE){


  if(type == "two.sided"){

    # assign appropriate weights to the current values
    w <- sapply(1:length(x), function(i) .get_log_weight_twosided(x[i], crit_x = crit_x[i, ], omega = omega))

    # compute the nominator
    nom <- stats::dnorm(x, mean, sd, log = TRUE) + w

    # compute the denominator
    denoms  <- matrix(stats::pnorm(crit_x[,1], mean, sd) - stats::pnorm(-crit_x[,1], mean, sd), ncol = 1)
    denoms[denoms < 0, 1]  <- 0 # check and correct for possibly negative numbers due to numerical imprecision
    if(length(omega) > 2){
      for(j in 2:(length(omega)-1)){
        denoms <- cbind(denoms, stats::pnorm(crit_x[,j], mean, sd) - stats::pnorm(-crit_x[,j], mean, sd) - apply(denoms, 1, sum))
        denoms[denoms[,j] < 0, j] <- 0 # check and correct for possibly negative numbers due to numerical imprecision
      }
    }
    denoms  <- cbind(denoms, 1 - apply(denoms, 1, sum))
    denoms  <- log(denoms) + log(matrix(omega, ncol = ncol(denoms), nrow = nrow(denoms), byrow = TRUE))
    denom   <- log(apply(exp(denoms), 1, sum))

    log_lik <- nom - denom

  }else if(type == "one.sided"){

    # assign appropriate weights to the current values
    w <- sapply(1:length(x), function(i) .get_log_weight_onesided(x[i], crit_x = crit_x[i, ], omega = omega))

    # compute the nominator
    nom <- stats::dnorm(x, mean, sd, log = TRUE) + w

    denoms  <- matrix(stats::pnorm(crit_x[,1], mean, sd), ncol = 1)
    denoms[denoms < 0, 1]  <- 0 # check and correct for possibly negative numbers due to numerical imprecision
    if(length(omega) > 2){
      for(j in 2:(length(omega)-1)){
        denoms <- cbind(denoms, stats::pnorm(crit_x[,j], mean, sd) - apply(denoms, 1, sum))
        denoms[denoms[,j] < 0, j] <- 0 # check and correct for possibly negative numbers due to numerical imprecision
      }
    }
    denoms  <- cbind(denoms, 1 - apply(denoms, 1, sum))
    denoms  <- log(denoms) + log(matrix(omega, ncol = ncol(denoms), nrow = nrow(denoms), byrow = TRUE))
    denom   <- log(apply(exp(denoms), 1, sum))

    log_lik <- nom - denom
  }


  if(log){
    return(log_lik)
  }else{
    return(exp(log_lik))
  }
}

# fast computation - no input check, pre-formatted for posterior predictive
.rwnorm_predict_fast      <- function(mean, sd, omega, crit_x, iter = 1){

  if(iter >= 50){
    # avoid getting stuck in an infinite recursion
    return(.rwnorm_predict_fast2(mean, sd, omega, crit_x))
  }

  # samples
  x <- stats::rnorm(length(mean), mean = mean, sd = sd)

  # find correct weight
  w <- rep(NA, length(mean))
  for(i in rev(seq_along(crit_x))){
    w[is.na(w) & x >= crit_x[i]] <- omega[is.na(w) & x >= crit_x[i], i + 1]
  }
  w[is.na(w) & x < crit_x[i]] <- omega[is.na(w) & x < crit_x[i], 1]

  # deal with computer precision errors from JAGS
  w[w > 1] <- 1
  w[w < 0] <- 0

  # assign publication status
  p <- stats::rbinom(length(mean), 1, prob = w) == 1

  # deal with possibility of sd = 0
  # (sampling never finishes, insert the mean value instead)
  x[!p & sd < sqrt(.Machine$double.ep)] <- mean
  p[!p & sd < sqrt(.Machine$double.ep)] <- TRUE

  # re-sample the missing estimates
  x[!p] <- NA

  if(any(!p)){
    x[!p] <- .rwnorm_predict_fast(mean[!p], sd[!p], omega[!p,,drop=FALSE], crit_x, iter = iter + 1)
  }

  return(x)
}
.rwnorm_predict_true_fast <- function(mean, tau, se, omega, crit_x){

  # samples
  xt <- stats::rnorm(length(mean), mean = mean, sd = tau)
  xo <- stats::rnorm(length(mean), mean = xt,   sd = se)

  # find correct weight
  w <- rep(NA, length(mean))
  for(i in rev(seq_along(crit_x))){
    w[is.na(w) & xo >= crit_x[i]] <- omega[is.na(w) & xo >= crit_x[i], i + 1]
  }
  w[is.na(w) & xo < crit_x[i]] <- omega[is.na(w) & xo < crit_x[i], 1]

  # deal with computer precision errors from JAGS
  w[w > 1] <- 1
  w[w < 0] <- 0

  # assign publication status
  p <- stats::rbinom(length(mean), 1, prob = w) == 1

  # re-sample the missing estimates
  xt[!p] <- NA

  if(any(!p)){
    xt[!p] <- .rwnorm_predict_true_fast(mean[!p], tau[!p], se, omega[!p,,drop=FALSE], crit_x)
  }

  return(xt)
}

# alternative version using inverse probability transform
# (helps when .rwnorm_predict_fast gets stuck)
.rwnorm_predict_fast2  <- function(mean, sd, omega, crit_x){

  # generate uniform random numbers
  u <- stats::runif(length(mean))

  # use fast quantile function
  x <- .qwnorm_fast(u, mean, sd, omega, crit_x)

  return(x)
}

# fast preformated p and q functions
.pwnorm_fast <- function(q, mean, sd, omega, crit_x){

  n     <- length(q)
  log_p <- rep(0, n)

  # compute standardizing constant once for all observations
  n_cuts    <- length(crit_x)
  n_weights <- ncol(omega)

  # pre-compute normal CDFs for all cutpoints and means/sds
  cdf_cuts <- matrix(0, nrow = n, ncol = n_cuts)
  for(j in seq_len(n_cuts)){
    cdf_cuts[, j] <- stats::pnorm(crit_x[j], mean, sd)
  }

  # compute denominators (standardizing constants)
  denoms <- matrix(0, nrow = n, ncol = n_weights)
  denoms[, 1] <- cdf_cuts[, 1]
  if(n_weights > 2){
    for(j in 2:(n_weights - 1)){
      denoms[, j] <- cdf_cuts[, j] - rowSums(denoms[, 1:(j-1), drop = FALSE])
    }
  }
  denoms[, n_weights] <- 1 - rowSums(denoms[, 1:(n_weights-1), drop = FALSE])

  # ensure non-negative (numerical precision)
  denoms[denoms < 0] <- 0

  log_denoms      <- log(denoms) + log(omega)
  log_denom_total <- log(rowSums(exp(log_denoms)))

  # compute numerators for each observation
  for(i in seq_len(n)){

    if(is.infinite(q[i])){
      log_p[i] <- if(q[i] < 0) -Inf else 0
      next
    }

    # find which interval q[i] falls into
    s <- 1
    temp_p <- numeric(0)

    # accumulate probabilities up to cutpoints
    while(s <= n_cuts && q[i] > crit_x[s]){
      temp_p_add <- cdf_cuts[i, s] - sum(temp_p)
      temp_p_add <- max(temp_p_add, 0)  # ensure non-negative
      temp_p <- c(temp_p, temp_p_add)
      s <- s + 1
    }

    # add final piece
    temp_p_add <- stats::pnorm(q[i], mean[i], sd[i]) - sum(temp_p)
    temp_p_add <- max(temp_p_add, 0)
    temp_p <- c(temp_p, temp_p_add)

    # weight by omega values
    log_p[i] <- log(sum(exp(log(temp_p) + log(omega[i, 1:length(temp_p)]))))
  }

  # subtract standardizing constant
  return(log_p - log_denom_total)
}
.qwnorm_fast <- function(p, mean, sd, omega, crit_x){

  n <- length(p)
  q <- rep(0, n)

  # vectorized bounds for optimization
  lower_bounds <- mean - 6 * sd  # approximately -6 sigma
  upper_bounds <- mean + 6 * sd  # approximately +6 sigma

  # for each observation, use Brent's method for root finding
  for(i in seq_len(n)){

    if(p[i] <= 0){
      q[i] <- -Inf
    } else if(p[i] >= 1){
      q[i] <- Inf
    } else {

      # objective function: F(x) - p = 0
      obj_fun <- function(x){
        exp(.pwnorm_fast(x, mean[i], sd[i], matrix(omega[i,], nrow = 1), crit_x)) - p[i]
      }

      # use uniroot for fast root finding
      tryCatch({
        result <- stats::uniroot(obj_fun,
                                 interval = c(lower_bounds[i], upper_bounds[i]),
                                 tol = 1e-8)
        q[i] <- result$root
      }, error = function(e){
        # fallback to wider search if initial bounds fail
        result <- stats::uniroot(obj_fun,
                                 interval = c(mean[i] - 10*sd[i], mean[i] + 10*sd[i]),
                                 tol = 1e-6)
        q[i] <- result$root
      })
    }
  }

  return(q)
}

# helper functions
.Xwnorm_check_input <- function(mean, sd, omega, steps, crit_x, type){

  BayesTools::check_real(mean, "mean", check_length = FALSE)
  BayesTools::check_real(sd, "sd", lower = 0, check_length = FALSE)
  BayesTools::check_real(as.vector(omega), "omega", lower = 0, upper = 1, check_length = FALSE)
  BayesTools::check_real(as.vector(steps), "steps", allow_NULL = !is.null(crit_x), lower = 0, upper = 1, check_length = FALSE, allow_bound = FALSE)
  BayesTools::check_char(type, "type", allow_values = c("two.sided", "one.sided"))
  BayesTools::check_real(as.vector(crit_x), "crit_x", allow_NULL = !is.null(steps), lower = if(type == "two.sided") 0 else -Inf, check_length = FALSE)

  return()
}
.Xwnorm_get_crit_x  <- function(steps, crit_x, mean, sd, type){

  if(!is.null(steps) & !is.null(crit_x))
    stop("Either 'steps' or 'crit_x' need to be specified.", call. = FALSE)


  if(!is.null(crit_x)){

    # use steps directly
    if(!is.matrix(crit_x)){
      crit_x <- matrix(crit_x, ncol = length(crit_x), nrow = length(mean), byrow = TRUE)
    }


  }else if(!is.null(steps)){

    # obtain critical values based on steps
    if(!is.matrix(steps)){
      steps <- matrix(steps, ncol = length(steps), nrow = length(mean), byrow = TRUE)
    }

    # reverse order of the steps
    steps <- steps[,ncol(steps):1, drop = FALSE]

    crit_z <- matrix(ncol = 0, nrow = length(mean))

    for(i in 1:ncol(steps)){
      if(type == "one.sided"){
        crit_z <- cbind(crit_z, stats::qnorm(steps[,i],   lower.tail = FALSE))
      }else if(type == "two.sided"){
        crit_z <- cbind(crit_z, stats::qnorm(steps[,i]/2, lower.tail = FALSE))
      }
    }

    crit_x <- crit_z * matrix(sd, ncol = ncol(crit_z), nrow = nrow(crit_z))
  }

  # check that the steps are increasing
  if(!all(sapply(1:nrow(crit_x), function(i) all(crit_x[i,] == cummax(crit_x[i,])))))
    stop("'steps'/'crit_x' argument must be inreasing.", call. = FALSE)

  return(crit_x)
}


#' @title Weighted multivariate normal distribution
#'
#' @description Density function for the weighted multivariate normal
#' distribution with \code{mean}, covariance matrix \code{sigma},
#' critical values \code{crit_x}, and weights \code{omega}.
#'
#' @param x quantiles.
#' @param p vector of probabilities.
#' @param mean mean
#' @param sigma covariance matrix.
#' @param crit_x vector of critical values defining steps.
#' @param omega vector of weights defining the probability
#' of observing a t-statistics between each of the two steps.
#' @param type type of weight function (defaults to \code{"two.sided"}).
#' @param log,log.p logical; if \code{TRUE}, probabilities
#' \code{p} are given as \code{log(p)}.
#'
#'
#' @return \code{.dwmnorm_fast} returns a density of the multivariate
#' weighted normal distribution.
#'
#' @name weighted_multivariate_normal
#'
#' @seealso \link[stats]{Normal}, [weighted_normal]
NULL

# fast computation - no input check, pre-formatted for bridge-sampling
.dwmnorm_fast   <- function(x, mean, sigma, omega, crit_x, type = "two.sided", log = TRUE){

  if(type == "two.sided"){

    log_w   <- sum(sapply(1:length(mean), function(k) .get_log_weight_twosided(x[k], crit_x[,k], omega)))
    log_lik <- mvtnorm::dmvnorm(x = x, mean = mean, sigma = sigma, log = TRUE) + log_w;

    log_std_constant <- .dwmnorm_log_std_constant_twosided(x, mean, sigma, crit_x, omega)

    log_lik <- log_lik - log_std_constant

  }else if(type == "one.sided"){

    log_w   <- sum(sapply(1:length(mean), function(k) .get_log_weight_onesided(x[k], crit_x[,k], omega)))
    log_lik <- mvtnorm::dmvnorm(x = x, mean = mean, sigma = sigma, log = TRUE) + log_w;

    log_std_constant <- .dwmnorm_log_std_constant_onesided(x, mean, sigma, crit_x, omega)

    log_lik <- log_lik - log_std_constant

  }


  if(log){
    return(log_lik)
  }else{
    return(exp(log_lik))
  }
}
.dwmnorm_v_fast <- function(x_v, mean_v, se2_v, tau2, rho, crit_x_v, omega, indx_v, type = "two.sided", log = TRUE){

  log_lik <- 0

  for(i in seq_along(indx_v)){

    # break into the individual observations
    if(i == 1){
      temp_K <- indx_v[i]
    }else{
      temp_K <- indx_v[i] - indx_v[i - 1]
    }
    indx_start <- indx_v[i] - temp_K + 1

    temp_index <- indx_start:indx_v[i]

    temp_x      <- x_v[temp_index]
    temp_mu     <- mean_v[temp_index]
    temp_sigma  <- diag(se2_v[temp_index] + tau2 * (1-rho), nrow = temp_K, ncol = temp_K) + matrix(tau2 * rho, nrow = temp_K, ncol = temp_K)
    if(type != "none"){
      temp_crit_x <- crit_x_v[,temp_index,drop = FALSE]
    }


    if(type == "two.sided"){

      temp_log_w   <- sum(sapply(1:length(temp_mu), function(k) .get_log_weight_twosided(temp_x[k], temp_crit_x[,k], omega)))
      temp_log_lik <- mvtnorm::dmvnorm(x = temp_x, mean = temp_mu, sigma = temp_sigma, log = TRUE) + temp_log_w;

      temp_log_std_constant <- .dwmnorm_log_std_constant_twosided(temp_x, temp_mu, temp_sigma, temp_crit_x, omega)

      log_lik <- log_lik + (temp_log_lik - temp_log_std_constant)

    }else if(type == "one.sided"){

      temp_log_w   <- sum(sapply(1:length(temp_mu), function(k) .get_log_weight_onesided(temp_x[k], temp_crit_x[,k], omega)))
      temp_log_lik <- mvtnorm::dmvnorm(x = temp_x, mean = temp_mu, sigma = temp_sigma, log = TRUE) + temp_log_w;

      temp_log_std_constant <- .dwmnorm_log_std_constant_onesided(temp_x, temp_mu, temp_sigma, temp_crit_x, omega)

      log_lik <- log_lik + (temp_log_lik - temp_log_std_constant)

    }else if(type == "none"){

      log_lik <- log_lik + mvtnorm::dmvnorm(x = temp_x, mean = temp_mu, sigma = temp_sigma, log = TRUE)

    }
  }


  if(log){
    return(log_lik)
  }else{
    return(exp(log_lik))
  }
}


# standardizing constant calculation
.dwmnorm_log_std_constant_onesided <- function(x, mu, sigma, crit_x, omega){

  std_constant <- 0

  # create a matrix indexing all sub-spaces created by the cut-points
  indexes <- as.matrix(do.call(expand.grid, args = lapply(1:length(mu), function(i) 1:length(omega))))

  for(i in 1:nrow(indexes)){

    # get the current lower and upper boundaries
    temp_lower <- .dwmnorm_lower_bound(crit_x, indexes[i,])
    temp_upper <- .dwmnorm_upper_bound(crit_x, indexes[i,])

    # get the probability and weights
    temp_log_prob   <- log(mvtnorm::pmvnorm(lower = temp_lower, upper = temp_upper, mean = mu, sigma = sigma))
    temp_log_weight <- sum(log(omega[unlist(indexes[i,])]))

    # add to the constant
    std_constant <- std_constant + exp(temp_log_prob + temp_log_weight)
  }

  return(log(std_constant))
}
.dwmnorm_log_std_constant_twosided <- function(x, mu, sigma, crit_x, omega){

  # turn the two-sided selection into one-sided selection
  crit_x <- rbind(-crit_x[nrow(crit_x):1,,drop=FALSE], crit_x)
  omega  <- c(rev(omega[-1]), omega)

  log_std_constant <- .dwmnorm_log_std_constant_onesided(x = x, mu = mu, sigma = sigma, crit_x = crit_x, omega = omega)

  return(log_std_constant)
}

# functions for assigning bounds
.dwmnorm_lower_bound <- function(crit_x, index){
  return(sapply(1:ncol(crit_x), function(k){
    if(index[k] == 1){
      return(-Inf)
    }else{
      return(crit_x[index[k] - 1, k])
    }
  }))
}
.dwmnorm_upper_bound <- function(crit_x, index){
  return(sapply(1:ncol(crit_x), function(k){
    if(index[k] == nrow(crit_x) + 1){
      return(Inf)
    }else{
      return(crit_x[index[k], k])
    }
  }))
}


#### general helper functionts
# functions for assigning weights & bounds
.get_log_weight_onesided <- function(x, crit_x, omega){

  # number of weights
  J <- length(omega)

  if(x >= crit_x[J - 1]){
    # return the last omega if x > last crit_x
    w <- omega[J]
  }else if(x < crit_x[1]){
    # return the first omega if x < first crit_x
    w <- omega[1]
  }else{
    # check the remaining cutpoints sequentially
    for(i in 2:J){
      if( (x >= crit_x[i - 1]) && (x < crit_x[i]) ){
        w = omega[i]
        break
      }
    }
  }

  return(log(w))
}
.get_log_weight_twosided <- function(x, crit_x, omega){

  # number of weights
  J <- length(omega)

  if(abs(x) >= crit_x[J - 1]){
    # return the last omega if x > last crit_x
    w <- omega[J]
  }else if(abs(x) < crit_x[1]){
    # return the first omega if x < first crit_x
    w <- omega[1]
  }else{
    # check the remaining cutpoints sequentially
    for(i in 2:J){
      if( (abs(x) >= crit_x[i - 1]) && (abs(x) < crit_x[i]) ){
        w = omega[i]
        break
      }
    }
  }

  return(log(w))
}

#### legacy code for weighted-t distribution ####
# #' @title Weighted t distribution
# #'
# #' @description Density, distribution function, quantile function
# #' and random generation for the weighted t distribution with
# #' \code{df} degrees of freedom, non-centrality parameter
# #' \code{ncp}, steps \code{steps} (or critical t-values
# #' \code{crit_t}), and weights \code{omega}.
# #'
# #' @param x,q vector of quantiles.
# #' @param p vector of probabilities.
# #' @param n number of observations. If length(n) > 1, the length
# #' is taken to be the number required.
# #' @param df degrees of freedom (> 0, maybe non-integer).
# #' \code{df = Inf} is allowed.
# #' @param ncp non-centrality parameter delta.
# #' @param steps vector of steps for the weight function.
# #' @param omega vector of weights defining the probability
# #' of observing a t-statistics between each of the two steps.
# #' @param crit_t vector of t-values defining steps
# #' (if \code{steps} are not supplied).
# #' @param type type of weight function (defaults to \code{"two.sided"}).
# #' @param log,log.p logical; if \code{TRUE}, probabilities
# #' \code{p} are given as \code{log(p)}.
# #' @param lower.tail logical; if \code{TRUE} (default), probabilities
# #' are \eqn{P[X \le x]}, otherwise, \eqn{P[X \ge x]}.
# #'
# #' @details The \code{df}, \code{ncp}, \code{steps}, \code{omega} can be
# #' supplied as a vectors (\code{df}, \code{ncp}) or matrices (\code{steps},
# #' \code{omega}) with length / number of rows equal to \code{x}/\code{q}/
# #' \code{p}. Otherwise, they are recycled to the length of the result.
# #'
# #' The functions quickly lose precision in the tails since they depend on
# #' sums of distribution functions of t distibution
# #' \code{\link[stats:TDist]{stats::pt}}. In cases where the density of
# #' t distribution cannot be computed by \code{\link[stats:TDist]{stats::dt}},
# #' the implementation switches to \code{\link[DPQ:dnt]{DPQ::dnt}}.
# #'
# #' @export dwt
# #' @export pwt
# #' @export qwt
# #' @export rwt
# #' @name weighted_t
# #'
# #' @seealso \link[stats]{TDist}, \link[DPQ]{dnt}
# NULL
#
# # @rdname weighted_t
# dwt <- function(x, df, ncp, steps = if(!is.null(crit_t)) NULL, omega, crit_t = if(!is.null(steps)) NULL, type = "two.sided", log = FALSE){
#
#   # set and check parameters
#   if(!is.numeric(x)  | !is.vector(x))stop("'x' must be a numeric vector.")
#   if(length(df)  != length(x))df  <- rep_len(df,  length(x))
#   if(length(ncp) != length(x))ncp <- rep_len(ncp, length(x))
#   if(!is.matrix(omega))omega <- matrix(omega, ncol = length(omega), nrow = length(x), byrow = TRUE)
#
#   crit_t <- .Xwt_get_crit_t(crit_t, steps, type, df)
#   # set a last weight to 1 if user forgets to do so
#   if(ncol(crit_t) == ncol(omega))omega <- cbind(omega, 1)
#
#   .Xwt_check_input(df, ncp, steps, omega, crit_t, type)
#
#
#   # do the actual computations
#   if(type == "two.sided"){
#
#     # assign appropriate weights to the current values
#     w <- sapply(1:length(x),function(i){
#       if(abs(x[i]) >= crit_t[i, ncol(omega)-1]){
#         return(log(1))
#       }else if(abs(x[i]) < crit_t[i,1]){
#         return(log(omega[i,1]))
#       }else{
#         for(j in 2:ncol(omega)){
#           if(abs(x[i]) < crit_t[i,j] & abs(x[i]) >= crit_t[i,j-1]){
#             return(log(omega[i,j]))
#           }
#         }
#       }
#     })
#
#
#     # compute the nominator
#     nom <- stats::dt(x, df, ncp, log = T)+w
#     # shift to different t-distribution computation if the classical one returns -Inf
#     if(any(is.infinite(nom))){
#       if(!try(requireNamespace("DPQ")))
#         stop("DPQ package needs to be installed. Run 'install.packages('DPQ')'")
#       nom[is.infinite(nom)] <- DPQ::dntJKBf(x[is.infinite(nom)], df[is.infinite(nom)], ncp[is.infinite(nom)], log = TRUE)+w[is.infinite(nom)]
#     }
#
#     # compute the denominator
#     denoms  <- matrix(stats::pt(crit_t[,1], df, ncp) - stats::pt(-crit_t[,1], df, ncp), ncol = 1)
#     denoms[denoms < 0, 1]  <- 0 # check and correct for possibly negative numbers due to numerical imprecision
#     if(ncol(omega) > 2){
#       for(j in 2:(ncol(omega)-1)){
#         denoms <- cbind(denoms, stats::pt(crit_t[,j], df, ncp) - stats::pt(-crit_t[,j], df, ncp) - apply(denoms, 1, sum))
#         denoms[denoms[,j] < 0, j] <- 0 # check and correct for possibly negative numbers due to numerical imprecision
#       }
#     }
#     denoms  <- cbind(denoms, 1 - apply(denoms, 1, sum))
#     denoms  <- log(denoms) + log(omega)
#     denom   <- log(apply(exp(denoms), 1, sum))
#
#     log_lik <- nom - denom
#
#   }else if(type == "one.sided"){
#
#     # assign appropriate weights to the current values
#     w <- sapply(1:length(x),function(i){
#       if(x[i] >= crit_t[i, ncol(omega)-1]){
#         return(log(1))
#       }else if(x[i] < crit_t[i,1]){
#         return(log(omega[i,1]))
#       }else{
#         for(j in 2:ncol(omega)){
#           if(x[i] < crit_t[i,j] & x[i] >= crit_t[i,j-1]){
#             return(log(omega[i,j]))
#           }
#         }
#       }
#     })
#
#
#     # compute the nominator
#     nom <- stats::dt(x, df, ncp, log = T)+w
#     # shift to different t-distribution computation if the classical one returns -Inf
#     if(any(is.infinite(nom))){
#       if(!try(requireNamespace("DPQ")))
#         stop("DPQ package needs to be installed. Run 'install.packages('DPQ')'")
#       nom[is.infinite(nom)] <- DPQ::dntJKBf(x[is.infinite(nom)], df[is.infinite(nom)], ncp[is.infinite(nom)], log = TRUE)+w[is.infinite(nom)]
#     }
#
#     denoms  <- matrix(stats::pt(crit_t[,1], df, ncp), ncol = 1)
#     denoms[denoms < 0, 1]  <- 0 # check and correct for possibly negative numbers due to numerical imprecision
#     if(ncol(omega) > 2){
#       for(j in 2:(ncol(omega)-1)){
#         denoms <- cbind(denoms, stats::pt(crit_t[,j], df, ncp) - apply(denoms, 1, sum))
#         denoms[denoms[,j] < 0, j] <- 0 # check and correct for possibly negative numbers due to numerical imprecision
#       }
#     }
#     denoms  <- cbind(denoms, 1 - apply(denoms, 1, sum))
#     denoms  <- log(denoms) + log(omega)
#     denom   <- log(apply(exp(denoms), 1, sum))
#
#     log_lik <- nom - denom
#   }
#
#
#   if(log){
#     return(log_lik)
#   }else{
#     return(exp(log_lik))
#   }
# }
# # @rdname weighted_t
# pwt <- function(q, df, ncp, steps = if(!is.null(crit_t)) NULL, omega, crit_t = if(!is.null(steps)) NULL, type = "two.sided", lower.tail = TRUE, log.p = FALSE){
#
#   # set and check parameters
#   if(!is.numeric(q) | !is.vector(q))stop("'q' must be a numeric vector.")
#   if(length(df)  != length(q))df  <- rep_len(df,  length(q))
#   if(length(ncp) != length(q))ncp <- rep_len(ncp, length(q))
#   if(!is.matrix(omega))omega <- matrix(omega, ncol = length(omega), nrow = length(q), byrow = TRUE)
#
#   crit_t <- .Xwt_get_crit_t(crit_t, steps, type, df)
#   # set a last weight to 1 if user forgets to do so
#   if(ncol(crit_t) == ncol(omega))omega <- cbind(omega, 1)
#
#   .Xwt_check_input(df, ncp, steps, omega, crit_t, type)
#
#   # work with two-sided as with one-sided
#   if(type == "two.sided"){
#     crit_t <- cbind(matrix(-crit_t[,rev(1:ncol(crit_t))], ncol = ncol(crit_t)), 0,  crit_t)
#     omega  <- cbind(matrix( omega[,rev(1:ncol(omega))],   ncol = ncol(omega)),      omega)
#   }
#
#   p <- .pwt(q, df, ncp, omega, crit_t, type, lower.tail)
#
#   if(log.p){
#     return(p)
#   }else{
#     return(exp(p))
#   }
#
# }
# # @rdname weighted_t
# qwt <- function(p, df, ncp, steps = if(!is.null(crit_t)) NULL, omega, crit_t = if(!is.null(steps)) NULL, type = "two.sided", lower.tail = TRUE, log.p = FALSE){
#
#   # set and check parameters
#   if(!is.numeric(p) | !is.vector(p))stop("'p' must be a numeric vector.")
#   if(length(df)  != length(p))df  <- rep_len(df,  length(p))
#   if(length(ncp) != length(p))ncp <- rep_len(ncp, length(p))
#   if(!is.matrix(omega))omega <- matrix(omega, ncol = length(omega), nrow = length(p), byrow = TRUE)
#   if(log.p)p <- exp(p)
#
#   crit_t <- .Xwt_get_crit_t(crit_t, steps, type, df)
#   # set a last weight to 1 if user forgets to do so
#   if(ncol(crit_t) == ncol(omega))omega <- cbind(omega, 1)
#
#   .Xwt_check_input(df, ncp, steps, omega, crit_t, type)
#
#   # work with two-sided as with one-sided
#   if(type == "two.sided"){
#     crit_t <- cbind(matrix(-crit_t[,rev(1:ncol(crit_t))], ncol = ncol(crit_t)), 0,  crit_t)
#     omega  <- cbind(matrix( omega[,rev(1:ncol(omega))],   ncol = ncol(omega)),      omega)
#   }
#
#
#   # yeah, this is slow and inefficient, but I was too lazy to write it in a better for now
#   q <- rep(NA, length(p))
#   for(i in 1:length(p)){
#     q[i] <- .qwt(p[i], df[i], ncp[i], omega[i,], crit_t[i,], type, lower.tail)
#   }
#
#   return(q)
# }
# # @rdname weighted_t
# rwt <- function(n, df, ncp, steps = if(!is.null(crit_t)) NULL, omega, crit_t = if(!is.null(steps)) NULL, type = "two.sided"){
#
#   # set and check parameters
#   if(length(n) > 1)n <- length(n)
#   if(length(df)  != n)df  <- rep_len(df,  n)
#   if(length(ncp) != n)ncp <- rep_len(ncp, n)
#   if(!is.matrix(omega))omega <- matrix(omega, ncol = length(omega), nrow = n, byrow = TRUE)
#
#   crit_t <- .Xwt_get_crit_t(crit_t, steps, type, df)
#   # set a last weight to 1 if user forgets to do so
#   if(ncol(crit_t) == ncol(omega))omega <- cbind(omega, 1)
#
#   .Xwt_check_input(df, ncp, steps, omega, crit_t, type)
#
#
#   # work with two-sided as with one-sided
#   if(type == "two.sided"){
#     crit_t <- cbind(matrix(-crit_t[,rev(1:ncol(crit_t))], ncol = ncol(crit_t)), 0,  crit_t)
#     omega  <- cbind(matrix( omega[,rev(1:ncol(omega))],   ncol = ncol(omega)),      omega)
#   }
#
#   # yeah, this is slow and inefficient, but I was too lazy to write it in a better for now
#   x <- rep(NA, n)
#   for(i in 1:n){
#     while(is.na(x[i])){
#       # sample a t-stat
#       temp_x <- stats::rt(1, df = df[i], ncp = ncp[i])
#
#       # find it's weight
#       if(temp_x >= crit_t[i, ncol(omega)-1]){
#         temp_omega <- omega[i, ncol(omega)]
#       }else if(temp_x < crit_t[i,1]){
#         temp_omega <- omega[i,1]
#       }else{
#         for(j in 2:ncol(omega)){
#           if(temp_x < crit_t[i,j] & temp_x >= crit_t[i,j-1]){
#             temp_omega <- omega[i,j]
#             break
#           }
#         }
#       }
#
#       if(stats::rbinom(1, 1, prob = temp_omega) == 1){
#         x[i] <- temp_x
#       }
#     }
#   }
#
#   # using the quantile function is also slow and the numerical imprecision kicks in
#   # p <- stats::runif(n, 0, 1)
#   # x <- rep(NA, n)
#   # for(i in 1:n){
#   #   x[i] <- .qwt(p[i], df[i], ncp[i], omega[i,], crit_t[i,], type, TRUE)
#   # }
#
#   return(x)
# }
#
# # helper functions
# .pwt <- function(q, df, ncp, omega, crit_t, type, lower.tail){
#   p  <- rep(NA, length(q))
#
#   if(lower.tail){
#     for(i in 1:length(q)){
#
#       s      <- 1
#       temp_p <- NULL
#       # do jumps by every cut-off
#       while(q[i] > crit_t[i,s]){
#
#         temp_p_add <- stats::pt(crit_t[i,s], df[i], ncp[i]) - sum(temp_p)
#         if(temp_p_add < 0)temp_p_add <- 0# check and correct for possibly negative numbers due to numerical imprecision
#         temp_p <- c(temp_p, temp_p_add)
#
#         s <- s + 1
#         if(s > ncol(crit_t))break
#       }
#
#       # add the last piece
#       temp_p_add <- stats::pt(q[i], df[i], ncp[i]) - sum(temp_p)
#       if(temp_p_add < 0)temp_p_add <- 0# check and correct for possibly negative numbers due to numerical imprecision
#       temp_p <- c(temp_p, temp_p_add)
#
#       #s <- s + 1
#
#       # weight by the weights
#       p[i] <- log(sum(exp(log(temp_p) + log(omega[i,1:s]))))
#     }
#   }else{
#     for(i in 1:length(q)){
#
#       s      <- ncol(omega)
#       temp_p <- NULL
#       # do jumps by every cut-off
#       while(q[i] < crit_t[i,s]){
#
#         temp_p_add <- stats::pt(crit_t[i,s], df[i], ncp[i], lower.tail = FALSE) - sum(temp_p)
#         if(temp_p_add < 0)temp_p_add <- 0# check and correct for possibly negative numbers due to numerical imprecision
#         temp_p <- c(temp_p, temp_p_add)
#
#         s <- s + 1
#         if(s > ncol(crit_t))break
#       }
#
#       # add the last piece
#       temp_p_add <- stats::pt(q[i], df[i], ncp[i], lower.tail = FALSE) - sum(temp_p)
#       if(temp_p_add < 0)temp_p_add <- 0# check and correct for possibly negative numbers due to numerical imprecision
#       temp_p <- c(temp_p, temp_p_add)
#
#       #s <- s - 1
#       # weight by the weights
#       p[i] <- log(sum(exp(log(temp_p) + log(omega[i,s:ncol(omega)]))))
#     }
#   }
#
#   # standardize
#   denoms  <- matrix(stats::pt(crit_t[,1], df, ncp), ncol = 1)
#   denoms[denoms < 0, 1]  <- 0 # check and correct for possibly negative numbers due to numerical imprecision
#   if(ncol(omega) > 2){
#     for(j in 2:(ncol(omega)-1)){
#       denoms <- cbind(denoms, stats::pt(crit_t[,j], df, ncp) - apply(denoms, 1, sum))
#       denoms[denoms[,j] < 0, j] <- 0 # check and correct for possibly negative numbers due to numerical imprecision
#     }
#   }
#   denoms  <- cbind(denoms, 1 - apply(denoms, 1, sum))
#   denoms  <- log(denoms) + log(omega)
#   denom   <- log(apply(exp(denoms), 1, sum))
#
#   p <- p - denom
#
#   return(p)
# }
# .qwt <- function(p, df, ncp, omega, crit_t, type, lower.tail){
#   return(stats::optim(
#     par = 0,
#     fn  = function(x, p, df, ncp, omega, crit_t, type, lower.tail){
#         (exp(.pwt(x, df, ncp, omega, crit_t, type, lower.tail)) - p)^2
#     },p = p, df = df, ncp = ncp, omega = matrix(omega, nrow = 1), crit_t = matrix(crit_t, nrow = 1), type = type, lower.tail = lower.tail,
#     method = "L-BFGS-B",
#     lower = -Inf, upper = Inf,
#     control = list(factr = 1e-12))$par)
# }
#
# .Xwt_check_input <- function(df, ncp, steps, omega, crit_t, type){
#
#   # check input
#   if(!is.numeric(df) | !is.vector(df))stop("'df' must be a numeric vector.")
#   if(any(df <= 0))stop("'df' must be positive.")
#   if(!is.numeric(omega)  | !is.matrix(omega))stop("'omega' must be numeric")
#   if(any(omega < 0) | any(omega > 1))stop("all 'omega' must be between 0 and 1")
#   if(ncol(crit_t) != ncol(omega) - 1)stop("'omega' is not specified properly - there must be N + 1 weights for N cutoffs")
#
# }
# .Xwt_get_crit_t  <- function(crit_t, steps, type, df){
#
#   if(!is.null(crit_t)){
#     if(!is.numeric(crit_t) | !(is.vector(crit_t) | is.matrix(crit_t)))stop("'crit_t' must be either a numeric vector or matrix.")
#   }else if(!is.null(steps)){
#
#     if(!all(steps > 0 & steps < 1))stop("'steps' must be higher than 0 and lower than 1.")
#
#     crit_t <- matrix(ncol = 0, nrow = length(df))
#     for(step in steps){
#       if(type == "one.sided"){
#         crit_t <- cbind(crit_t, stats::qt(step, df, 0, lower.tail = FALSE))
#       }else if(type == "two.sided"){
#         crit_t <- cbind(crit_t, stats::qt(step/2, df, 0, lower.tail = FALSE))
#       }
#     }
#
#   }else{
#     stop("Either 'crit_t' or 'steps' need to be supplied.")
#   }
#
#   # modify crit_t into matrix form
#   if(type == "two.sided"){
#     crit_t <- abs(crit_t)
#   }
#   if(length(crit_t) == 1){
#     crit_t <- matrix(crit_t, nrow = length(df))
#   }else{
#     if(!is.matrix(crit_t)){
#       crit_t <- matrix(crit_t, nrow = length(t), ncol = length(crit_t), byrow = TRUE)
#     }else if(nrow(crit_t) != length(df))stop("Dimensions of 't' and 'crit_t' don't match.")
#   }
#
#   return(crit_t)
#
# }
#
# ### for inner usage in the package
# # - no input checking and reformatting
# # - only one omega
# .dwt_fast    <- function(t, df, ncp, omega, crit_t, type = "two.sided", log = TRUE){
#
#
#   if(type == "two.sided"){
#
#     # assign appropriate weights to the current values
#     w <- sapply(1:length(t),function(i){
#       if(abs(t[i]) >= crit_t[i, length(omega)-1]){
#         return(log(1))
#       }else if(abs(t[i]) < crit_t[i,1]){
#         return(log(omega[1]))
#       }else{
#         for(j in 2:length(omega)){
#           if(abs(t[i]) < crit_t[i,j] & abs(t[i]) >= crit_t[i,j-1]){
#             return(log(omega[j]))
#           }
#         }
#       }
#     })
#
#
#     # compute the nominator
#     nom <- stats::dt(t, df, ncp, log = T)+w
#     # shift to different t-distribution computation if the classical one returns -Inf
#     if(any(is.infinite(nom))){
#       if(!try(requireNamespace("DPQ")))
#         stop("DPQ package needs to be installed. Run 'install.packages('DPQ')'")
#       nom[is.infinite(nom)] <- DPQ::dntJKBf(t[is.infinite(nom)], df[is.infinite(nom)], ncp[is.infinite(nom)], log = TRUE)+w[is.infinite(nom)]
#     }
#
#     # compute the denominator
#     denoms  <- matrix(stats::pt(crit_t[,1], df, ncp) - stats::pt(-crit_t[,1], df, ncp), ncol = 1)
#     denoms[denoms < 0, 1]  <- 0 # check and correct for possibly negative numbers due to numerical imprecision
#     if(length(omega) > 2){
#       for(j in 2:(length(omega)-1)){
#         denoms <- cbind(denoms, stats::pt(crit_t[,j], df, ncp) - stats::pt(-crit_t[,j], df, ncp) - apply(denoms, 1, sum))
#         denoms[denoms[,j] < 0, j] <- 0 # check and correct for possibly negative numbers due to numerical imprecision
#       }
#     }
#     denoms  <- cbind(denoms, 1 - apply(denoms, 1, sum))
#     denoms  <- log(denoms) + log(matrix(omega, ncol = ncol(denoms), nrow = nrow(denoms), byrow = TRUE))
#     denom   <- log(apply(exp(denoms), 1, sum))
#
#     log_lik <- nom - denom
#
#   }else if(type == "one.sided"){
#
#     # assign appropriate weights to the current values
#     w <- sapply(1:length(t),function(i){
#       if(t[i] >= crit_t[i, length(omega)-1]){
#         return(log(1))
#       }else if(t[i] < crit_t[i,1]){
#         return(log(omega[1]))
#       }else{
#         for(j in 2:length(omega)){
#           if(t[i] < crit_t[i,j] & t[i] >= crit_t[i,j-1]){
#             return(log(omega[j]))
#           }
#         }
#       }
#     })
#
#     # compute the nominator
#     nom <- stats::dt(t, df, ncp, log = T)+w
#     # shift to different t-distribution computation if the classical one returns -Inf
#     if(any(is.infinite(nom))){
#       if(!try(requireNamespace("DPQ")))
#         stop("DPQ package needs to be installed. Run 'install.packages('DPQ')'")
#       nom[is.infinite(nom)] <- DPQ::dntJKBf(t[is.infinite(nom)], df[is.infinite(nom)], ncp[is.infinite(nom)], log = TRUE)+w[is.infinite(nom)]
#     }
#
#     denoms  <- matrix(stats::pt(crit_t[,1], df, ncp), ncol = 1)
#     denoms[denoms < 0, 1]  <- 0 # check and correct for possibly negative numbers due to numerical imprecision
#     if(length(omega) > 2){
#       for(j in 2:(length(omega)-1)){
#         denoms <- cbind(denoms, stats::pt(crit_t[,j], df, ncp) - apply(denoms, 1, sum))
#         denoms[denoms[,j] < 0, j] <- 0 # check and correct for possibly negative numbers due to numerical imprecision
#       }
#     }
#     denoms  <- cbind(denoms, 1 - apply(denoms, 1, sum))
#     denoms  <- log(denoms) + log(matrix(omega, ncol = ncol(denoms), nrow = nrow(denoms), byrow = TRUE))
#     denom   <- log(apply(exp(denoms), 1, sum))
#
#     log_lik <- nom - denom
#   }
#
#
#   if(log){
#     return(log_lik)
#   }else{
#     return(exp(log_lik))
#   }
# }
