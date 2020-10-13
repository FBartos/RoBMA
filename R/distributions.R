#' @title Weighted t distribution
#'
#' @description Density, distribution function, quantile function
#' and random generation for the weighted t distribution with
#' \code{df} degrees of freedom, non-centrality parameter
#' \code{ncp}, steps \code{steps} (or critical t-values
#' \code{crit_t}), and weights \code{omega}.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If length(n) > 1, the length
#' is taken to be the number required.
#' @param df degrees of freedom (> 0, maybe non-integer).
#' \code{df = Inf} is allowed.
#' @param ncp non-centrality parameter delta.
#' @param steps vector of steps for the weight function.
#' @param omega vector of weights defining the probability
#' of observing a t-statistics between each of the two steps.
#' @param crit_t vector of t-values defining steps
#' (if \code{steps} are not supplied).
#' @param type type of weight function (defaults to \code{"two.sided"}).
#' @param log,log.p logical; if \code{TRUE}, probabilities
#' \code{p} are given as \code{log(p)}.
#' @param lower.tail logical; if \code{TRUE} (default), probabilities
#' are \eqn{P[X \le x]}, otherwise, \eqn{P[X \ge x]}.
#'
#' @details The \code{df}, \code{ncp}, \code{steps}, \code{omega} can be
#' supplied as a vectors (\code{df}, \code{ncp}) or matrices (\code{steps},
#' \code{omega}) with length / number of rows equal to \code{x}/\code{q}/
#' \code{p}. Otherwise, they are recycled to the length of the result.
#'
#' The functions quickly lose precision in the tails since they depend on
#' sums of distribution functions of t distibution
#' \code{\link[stats:TDist]{stats::pt}}. In cases where the density of
#' t distribution cannot be computed by \code{\link[stats:TDist]{stats::dt}},
#' the implementation switches to \code{\link[DPQ:dnt]{DPQ::dnt}}.
#'
#' @export dwt
#' @export pwt
#' @export qwt
#' @export rwt
#' @name weightedt
#'
#' @seealso \link[stats]{Normal}, \link[DPQ]{dnt}
NULL

#' @rdname weightedt
dwt <- function(x, df, ncp, steps = if(!is.null(crit_t)) NULL, omega, crit_t = if(!is.null(steps)) NULL, type = "two.sided", log = FALSE){

  # set and check parameters
  if(!is.numeric(x)  | !is.vector(x))stop("'x' must be a numeric vector.")
  if(length(df)  != length(x))df  <- rep_len(df,  length(x))
  if(length(ncp) != length(x))ncp <- rep_len(ncp, length(x))
  if(!is.matrix(omega))omega <- matrix(omega, ncol = length(omega), nrow = length(x), byrow = TRUE)

  crit_t <- .Xwt_get_crit_t(crit_t, steps, type, df)
  # set a last weight to 1 if user forgets to do so
  if(ncol(crit_t) == ncol(omega))omega <- cbind(omega, 1)

  .Xwt_check_input(df, ncp, steps, omega, crit_t, type)


  # do the actual computations
  if(type == "two.sided"){

    # assign appropriate weights to the current values
    w <- sapply(1:length(x),function(i){
      if(abs(x[i]) >= crit_t[i, ncol(omega)-1]){
        return(log(1))
      }else if(abs(x[i]) < crit_t[i,1]){
        return(log(omega[i,1]))
      }else{
        for(j in 2:ncol(omega)){
          if(abs(x[i]) < crit_t[i,j] & abs(x[i]) >= crit_t[i,j-1]){
            return(log(omega[i,j]))
          }
        }
      }
    })


    # compute the nominator
    nom <- stats::dt(x, df, ncp, log = T)+w
    # shift to different t-distribution computation if the classical one returns -Inf
    if(any(is.infinite(nom))){
       nom[is.infinite(nom)] <- DPQ::dntJKBf(x[is.infinite(nom)], df[is.infinite(nom)], ncp[is.infinite(nom)], log = TRUE)+w[is.infinite(nom)]
    }

    # compute the denominator
    denoms  <- matrix(stats::pt(crit_t[,1], df, ncp) - stats::pt(-crit_t[,1], df, ncp), ncol = 1)
    denoms[denoms < 0, 1]  <- 0 # check and correct for possibly negative numbers due to numerical imprecission
    if(ncol(omega) > 2){
      for(j in 2:(ncol(omega)-1)){
        denoms <- cbind(denoms, stats::pt(crit_t[,j], df, ncp) - stats::pt(-crit_t[,j], df, ncp) - apply(denoms, 1, sum))
        denoms[denoms[,j] < 0, j] <- 0 # check and correct for possibly negative numbers due to numerical imprecission
      }
    }
    denoms  <- cbind(denoms, 1 - apply(denoms, 1, sum))
    denoms  <- log(denoms) + log(omega)
    denom   <- log(apply(exp(denoms), 1, sum))

    log_lik <- nom - denom

  }else if(type == "one.sided"){

    # assign appropriate weights to the current values
    w <- sapply(1:length(x),function(i){
      if(x[i] >= crit_t[i, ncol(omega)-1]){
        return(log(1))
      }else if(x[i] < crit_t[i,1]){
        return(log(omega[i,1]))
      }else{
        for(j in 2:ncol(omega)){
          if(x[i] < crit_t[i,j] & x[i] >= crit_t[i,j-1]){
            return(log(omega[i,j]))
          }
        }
      }
    })


    # compute the nominator
    nom <- stats::dt(x, df, ncp, log = T)+w
    # shift to different t-distribution computation if the classical one returns -Inf
    if(any(is.infinite(nom))){
      nom[is.infinite(nom)] <- DPQ::dntJKBf(x[is.infinite(nom)], df[is.infinite(nom)], ncp[is.infinite(nom)], log = TRUE)+w[is.infinite(nom)]
    }

    denoms  <- matrix(stats::pt(crit_t[,1], df, ncp), ncol = 1)
    denoms[denoms < 0, 1]  <- 0 # check and correct for possibly negative numbers due to numerical imprecission
    if(ncol(omega) > 2){
      for(j in 2:(ncol(omega)-1)){
        denoms <- cbind(denoms, stats::pt(crit_t[,j], df, ncp) - apply(denoms, 1, sum))
        denoms[denoms[,j] < 0, j] <- 0 # check and correct for possibly negative numbers due to numerical imprecission
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
#' @rdname weightedt
pwt <- function(q, df, ncp, steps = if(!is.null(crit_t)) NULL, omega, crit_t = if(!is.null(steps)) NULL, type = "two.sided", lower.tail = TRUE, log.p = FALSE){

  # set and check parameters
  if(!is.numeric(q) | !is.vector(q))stop("'q' must be a numeric vector.")
  if(length(df)  != length(q))df  <- rep_len(df,  length(q))
  if(length(ncp) != length(q))ncp <- rep_len(ncp, length(q))
  if(!is.matrix(omega))omega <- matrix(omega, ncol = length(omega), nrow = length(q), byrow = TRUE)

  crit_t <- .Xwt_get_crit_t(crit_t, steps, type, df)
  # set a last weight to 1 if user forgets to do so
  if(ncol(crit_t) == ncol(omega))omega <- cbind(omega, 1)

  .Xwt_check_input(df, ncp, steps, omega, crit_t, type)

  # work with two-sided as with one-sided
  if(type == "two.sided"){
    crit_t <- cbind(matrix(-crit_t[,rev(1:ncol(crit_t))], ncol = ncol(crit_t)), 0,  crit_t)
    omega  <- cbind(matrix( omega[,rev(1:ncol(omega))],   ncol = ncol(omega)),      omega)
  }

  p <- .pwt(q, df, ncp, omega, crit_t, type, lower.tail)

  if(log.p){
    return(p)
  }else{
    return(exp(p))
  }

}
#' @rdname weightedt
qwt <- function(p, df, ncp, steps = if(!is.null(crit_t)) NULL, omega, crit_t = if(!is.null(steps)) NULL, type = "two.sided", lower.tail = TRUE, log.p = FALSE){

  # set and check parameters
  if(!is.numeric(p) | !is.vector(p))stop("'p' must be a numeric vector.")
  if(length(df)  != length(p))df  <- rep_len(df,  length(p))
  if(length(ncp) != length(p))ncp <- rep_len(ncp, length(p))
  if(!is.matrix(omega))omega <- matrix(omega, ncol = length(omega), nrow = length(p), byrow = TRUE)
  if(log.p)p <- exp(p)

  crit_t <- .Xwt_get_crit_t(crit_t, steps, type, df)
  # set a last weight to 1 if user forgets to do so
  if(ncol(crit_t) == ncol(omega))omega <- cbind(omega, 1)

  .Xwt_check_input(df, ncp, steps, omega, crit_t, type)

  # work with two-sided as with one-sided
  if(type == "two.sided"){
    crit_t <- cbind(matrix(-crit_t[,rev(1:ncol(crit_t))], ncol = ncol(crit_t)), 0,  crit_t)
    omega  <- cbind(matrix( omega[,rev(1:ncol(omega))],   ncol = ncol(omega)),      omega)
  }


  # yeah, this is slow and inefficient, but I was too lazy to write it in a better for now
  q <- rep(NA, length(p))
  for(i in 1:length(p)){
    q[i] <- .qwt(p[i], df[i], ncp[i], omega[i,], crit_t[i,], type, lower.tail)
  }

  return(q)
}
#' @rdname weightedt
rwt <- function(n, df, ncp, steps = if(!is.null(crit_t)) NULL, omega, crit_t = if(!is.null(steps)) NULL, type = "two.sided"){

  # set and check parameters
  if(length(n) > 1)n <- length(n)
  if(length(df)  != n)df  <- rep_len(df,  n)
  if(length(ncp) != n)ncp <- rep_len(ncp, n)
  if(!is.matrix(omega))omega <- matrix(omega, ncol = length(omega), nrow = n, byrow = TRUE)

  crit_t <- .Xwt_get_crit_t(crit_t, steps, type, df)
  # set a last weight to 1 if user forgets to do so
  if(ncol(crit_t) == ncol(omega))omega <- cbind(omega, 1)

  .Xwt_check_input(df, ncp, steps, omega, crit_t, type)


  # work with two-sided as with one-sided
  if(type == "two.sided"){
    crit_t <- cbind(matrix(-crit_t[,rev(1:ncol(crit_t))], ncol = ncol(crit_t)), 0,  crit_t)
    omega  <- cbind(matrix( omega[,rev(1:ncol(omega))],   ncol = ncol(omega)),      omega)
  }

  # yeah, this is slow and inefficient, but I was too lazy to write it in a better for now
  x <- rep(NA, n)
  for(i in 1:n){
    while(is.na(x[i])){
      # sample a t-stat
      temp_x <- stats::rt(1, df = df[i], ncp = ncp[i])

      # find it's weight
      if(temp_x >= crit_t[i, ncol(omega)-1]){
        temp_omega <- omega[i, ncol(omega)]
      }else if(temp_x < crit_t[i,1]){
        temp_omega <- omega[i,1]
      }else{
        for(j in 2:ncol(omega)){
          if(temp_x < crit_t[i,j] & temp_x >= crit_t[i,j-1]){
            temp_omega <- omega[i,j]
            break
          }
        }
      }

      if(stats::rbinom(1, 1, prob = temp_omega) == 1){
        x[i] <- temp_x
      }
    }
  }

  # using the quantile function is also slow and the numerical imprecission kicks in
  # p <- stats::runif(n, 0, 1)
  # x <- rep(NA, n)
  # for(i in 1:n){
  #   x[i] <- .qwt(p[i], df[i], ncp[i], omega[i,], crit_t[i,], type, TRUE)
  # }

  return(x)
}

# helper functions
.pwt <- function(q, df, ncp, omega, crit_t, type, lower.tail){
  p  <- rep(NA, length(q))

  if(lower.tail){
    for(i in 1:length(q)){

      s      <- 1
      temp_p <- NULL
      # do jumps by every cut-off
      while(q[i] > crit_t[i,s]){

        temp_p_add <- stats::pt(crit_t[i,s], df[i], ncp[i]) - sum(temp_p)
        if(temp_p_add < 0)temp_p_add <- 0# check and correct for possibly negative numbers due to numerical imprecission
        temp_p <- c(temp_p, temp_p_add)

        s <- s + 1
        if(s > ncol(crit_t))break
      }

      # add the last piece
      temp_p_add <- stats::pt(q[i], df[i], ncp[i]) - sum(temp_p)
      if(temp_p_add < 0)temp_p_add <- 0# check and correct for possibly negative numbers due to numerical imprecission
      temp_p <- c(temp_p, temp_p_add)

      #s <- s + 1

      # weight by the weights
      p[i] <- log(sum(exp(log(temp_p) + log(omega[i,1:s]))))
    }
  }else{
    for(i in 1:length(q)){

      s      <- ncol(omega)
      temp_p <- NULL
      # do jumps by every cut-off
      while(q[i] < crit_t[i,s]){

        temp_p_add <- stats::pt(crit_t[i,s], df[i], ncp[i], lower.tail = FALSE) - sum(temp_p)
        if(temp_p_add < 0)temp_p_add <- 0# check and correct for possibly negative numbers due to numerical imprecission
        temp_p <- c(temp_p, temp_p_add)

        s <- s + 1
        if(s > ncol(crit_t))break
      }

      # add the last piece
      temp_p_add <- stats::pt(q[i], df[i], ncp[i], lower.tail = FALSE) - sum(temp_p)
      if(temp_p_add < 0)temp_p_add <- 0# check and correct for possibly negative numbers due to numerical imprecission
      temp_p <- c(temp_p, temp_p_add)

      #s <- s - 1
      # weight by the weights
      p[i] <- log(sum(exp(log(temp_p) + log(omega[i,s:ncol(omega)]))))
    }
  }

  # standardize
  denoms  <- matrix(stats::pt(crit_t[,1], df, ncp), ncol = 1)
  denoms[denoms < 0, 1]  <- 0 # check and correct for possibly negative numbers due to numerical imprecission
  if(ncol(omega) > 2){
    for(j in 2:(ncol(omega)-1)){
      denoms <- cbind(denoms, stats::pt(crit_t[,j], df, ncp) - apply(denoms, 1, sum))
      denoms[denoms[,j] < 0, j] <- 0 # check and correct for possibly negative numbers due to numerical imprecission
    }
  }
  denoms  <- cbind(denoms, 1 - apply(denoms, 1, sum))
  denoms  <- log(denoms) + log(omega)
  denom   <- log(apply(exp(denoms), 1, sum))

  p <- p - denom

  return(p)
}
.qwt <- function(p, df, ncp, omega, crit_t, type, lower.tail){
  return(stats::optim(
    par = 0,
    fn  = function(x, p, df, ncp, omega, crit_t, type, lower.tail){
        (exp(.pwt(x, df, ncp, omega, crit_t, type, lower.tail)) - p)^2
    },p = p, df = df, ncp = ncp, omega = matrix(omega, nrow = 1), crit_t = matrix(crit_t, nrow = 1), type = type, lower.tail = lower.tail,
    method = "L-BFGS-B",
    lower = -Inf, upper = Inf,
    control = list(factr = 1e-12))$par)
}

.Xwt_check_input <- function(df, ncp, steps, omega, crit_t, type){

  # check input
  if(!is.numeric(df) | !is.vector(df))stop("'df' must be a numeric vector.")
  if(any(df <= 0))stop("'df' must be positive.")
  if(!is.numeric(omega)  | !is.matrix(omega))stop("'omega' must be numeric")
  if(any(omega < 0) | any(omega > 1))stop("all 'omega' must be between 0 and 1")
  if(ncol(crit_t) != ncol(omega) - 1)stop("'omega' is not specified properly - there must be N + 1 weights for N cutoffs")

}
.Xwt_get_crit_t  <- function(crit_t, steps, type, df){

  if(!is.null(crit_t)){
    if(!is.numeric(crit_t) | !(is.vector(crit_t) | is.matrix(crit_t)))stop("'crit_t' must be either a numeric vector or matrix.")
  }else if(!is.null(steps)){

    if(!all(steps > 0 & steps < 1))stop("'steps' must be higher than 0 and lower than 1.")

    crit_t <- matrix(ncol = 0, nrow = length(df))
    for(step in steps){
      if(type == "one.sided"){
        crit_t <- cbind(crit_t, stats::qt(step, df, 0, lower.tail = FALSE))
      }else if(type == "two.sided"){
        crit_t <- cbind(crit_t, stats::qt(step/2, df, 0, lower.tail = FALSE))
      }
    }

  }else{
    stop("Either 'crit_t' or 'steps' need to be supplied.")
  }

  # modify crit_t into matrix form
  if(type == "two.sided"){
    crit_t <- abs(crit_t)
  }
  if(length(crit_t) == 1){
    crit_t <- matrix(crit_t, nrow = length(df))
  }else{
    if(!is.matrix(crit_t)){
      crit_t <- matrix(crit_t, nrow = length(t), ncol = length(crit_t), byrow = TRUE)
    }else if(nrow(crit_t) != length(df))stop("Dimensions of 't' and 'crit_t' don't match.")
  }

  return(crit_t)

}

### for inner usage in the package
# - no input checking and reformating
# - only one omega
.dwt_fast   <- function(t, df, ncp, omega, crit_t, type = "two.sided", log = TRUE){


  if(type == "two.sided"){

    # assign appropriate weights to the current values
    w <- sapply(1:length(t),function(i){
      if(abs(t[i]) >= crit_t[i, length(omega)-1]){
        return(log(1))
      }else if(abs(t[i]) < crit_t[i,1]){
        return(log(omega[1]))
      }else{
        for(j in 2:length(omega)){
          if(abs(t[i]) < crit_t[i,j] & abs(t[i]) >= crit_t[i,j-1]){
            return(log(omega[j]))
          }
        }
      }
    })


    # compute the nominator
    nom <- stats::dt(t, df, ncp, log = T)+w
    # shift to different t-distribution computation if the classical one returns -Inf
    if(any(is.infinite(nom))){
      nom[is.infinite(nom)] <- DPQ::dntJKBf(t[is.infinite(nom)], df[is.infinite(nom)], ncp[is.infinite(nom)], log = TRUE)+w[is.infinite(nom)]
    }

    # compute the denominator
    denoms  <- matrix(stats::pt(crit_t[,1], df, ncp) - stats::pt(-crit_t[,1], df, ncp), ncol = 1)
    denoms[denoms < 0, 1]  <- 0 # check and correct for possibly negative numbers due to numerical imprecission
    if(length(omega) > 2){
      for(j in 2:(length(omega)-1)){
        denoms <- cbind(denoms, stats::pt(crit_t[,j], df, ncp) - stats::pt(-crit_t[,j], df, ncp) - apply(denoms, 1, sum))
        denoms[denoms[,j] < 0, j] <- 0 # check and correct for possibly negative numbers due to numerical imprecission
      }
    }
    denoms  <- cbind(denoms, 1 - apply(denoms, 1, sum))
    denoms  <- log(denoms) + log(matrix(omega, ncol = ncol(denoms), nrow = nrow(denoms), byrow = TRUE))
    denom   <- log(apply(exp(denoms), 1, sum))

    log_lik <- nom - denom

  }else if(type == "one.sided"){

    # assign appropriate weights to the current values
    w <- sapply(1:length(t),function(i){
      if(t[i] >= crit_t[i, length(omega)-1]){
        return(log(1))
      }else if(t[i] < crit_t[i,1]){
        return(log(omega[1]))
      }else{
        for(j in 2:length(omega)){
          if(t[i] < crit_t[i,j] & t[i] >= crit_t[i,j-1]){
            return(log(omega[j]))
          }
        }
      }
    })

    # compute the nominator
    nom <- stats::dt(t, df, ncp, log = T)+w
    # shift to different t-distribution computation if the classical one returns -Inf
    if(any(is.infinite(nom))){
      nom[is.infinite(nom)] <- DPQ::dntJKBf(t[is.infinite(nom)], df[is.infinite(nom)], ncp[is.infinite(nom)], log = TRUE)+w[is.infinite(nom)]
    }

    denoms  <- matrix(stats::pt(crit_t[,1], df, ncp), ncol = 1)
    denoms[denoms < 0, 1]  <- 0 # check and correct for possibly negative numbers due to numerical imprecission
    if(length(omega) > 2){
      for(j in 2:(length(omega)-1)){
        denoms <- cbind(denoms, stats::pt(crit_t[,j], df, ncp) - apply(denoms, 1, sum))
        denoms[denoms[,j] < 0, j] <- 0 # check and correct for possibly negative numbers due to numerical imprecission
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
