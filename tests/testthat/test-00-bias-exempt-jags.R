probe_state("test-00-bias-exempt-jags.R")
test_that("JAGS selection blocks preserve the full Gaussian law of exempt outcomes", {

  skip_if_not_installed("rjags")
  .load_RoBMA_module()
  fixtures <- .native_selection_jags_fixtures()
  cluster_rules <- .selection_joint_cluster_quadrature_rules(c(7L, 15L, 31L, 63L))
  factor_rules <- .selection_joint_factor_quadrature_rules()
  # Integrate the Gaussian over an eligible rectangle, conditioning the second
  # coordinate on the first. This reference uses no selection kernels.
  rectangle_mass <- function(lo, hi, mean, covariance) {
    sd <- sqrt(diag(covariance))
    if (length(lo) == 1L) return(pnorm(hi, mean, sd) - pnorm(lo, mean, sd))
    slope <- covariance[2, 1] / covariance[1, 1]
    conditional_sd <- sqrt(covariance[2, 2] - slope * covariance[1, 2])
    integrate(function(x) {
      center <- mean[2] + slope * (x - mean[1])
      dnorm(x, mean[1], sd[1]) *
        (pnorm(hi[2], center, conditional_sd) - pnorm(lo[2], center, conditional_sd))
    }, lo[1], hi[1], rel.tol = 1e-10)$value
  }
  # One-sided p-bins for the product (0) and best (1) rules.
  bounds <- c(Inf, 1, 0, -Inf)
  cell_mass <- function(bins, mean, covariance) {
    rectangle_mass(bounds[bins + 1L], bounds[bins], mean, covariance)
  }
  # Two-sided best (2): the largest eligible |z| falls below .3, between .3
  # and 1.4, or above 1.4, with weights .8, .2 and 1.
  two_sided_limits <- c(.3, 1.4)
  two_sided_omega <- c(.8, .2, 1)
  # Three correlated outcomes with the most significant result in the middle.
  # The eligible rows' covariance is large enough that dropping it moves every
  # dense and conditioned reference by more than the tolerance.
  sampling <- matrix(c(.8, .32, -.4, .32, 1.1, .16, -.4, .16, .93), 3L)
  for (kind in names(fixtures$cases)) {
    spec <- fixtures$cases[[kind]]
    data <- utils::modifyList(fixtures$base, spec$extra)
    data$mu <- numeric(3L)
    data$sei <- rep(1, 3L)
    data$mode <- 1L
    data$covariance <- sampling[lower.tri(sampling, diag = TRUE)]
    # Larger conditioned loadings make a dropped eligible covariance visible.
    data$loading <- if (kind == "conditioned") c(.55, -.35, .48) else c(.25, -.35, .18)
    data$sd <- c(.6, .9, .73)
    data$factor_loading <- matrix(c(data$loading, .4, .13, -.22), 3L)
    data$residual <- c(.4, .5, .6)
    data$groups <- rep(1L, 3L)
    data$points <- data$maximum <- 2048L
    data$scrambles <- 8L
    data$qmc <- BayesTools::selection_qmc_design(
      switch(kind, dense = 6L, cluster = 2L, factor = 4L, conditioned = 6L),
      data$points, data$scrambles, seed = 112L)
    quadrature <- if (kind == "factor") factor_rules else cluster_rules
    controls <- c("nodes", "log_weights", "orders")
    data[controls] <- quadrature[controls]
    total <- switch(kind,
      dense = sampling,
      cluster = diag(data$sd^2) + tcrossprod(data$loading),
      factor = diag(data$sd^2) + tcrossprod(data$factor_loading),
      conditioned = sampling + diag(data$residual) + tcrossprod(data$loading)
    )
    candidate <- if (kind == "conditioned") total - sampling else total
    for (rule in 0:2) {
      data$rule <- rule
      two_sided <- rule == 2L
      # Under the two-sided rule the eligible rows' signed and absolute
      # winners differ, so a signed comparison changes the observed weight.
      data$y <- c(if (two_sided) .1 else .4, 1.7, -.6)
      gaussian <- mvtnorm::dmvnorm(data$y, sigma = total, log = TRUE)
      # Retained sources condition on every observed outcome, including exempt rows.
      mean <- if (kind == "conditioned") drop(sampling %*% solve(total, data$y)) else data$mu
      data$omega <- if (two_sided) c(1, .2, .8, .2, 1) else c(1, .2, .8)
      data$lower <- if (two_sided) c(1.4, .3, -.3, -1.4, -1e300) else c(1, 0, -1e300)
      data$upper <- if (two_sided) c(1e300, 1.4, .3, -.3, -1.4) else c(1e300, 1, 0)
      # The exempt middle row is the most significant; under the two-sided
      # rule the negative third row is the most significant eligible one.
      bin_sets <- if (two_sided) {
        list(c(3L, 0L, 4L), c(0L, 0L, 4L), rep(0L, 3L))
      } else list(c(2L, 0L, 3L), c(0L, 0L, 3L), rep(0L, 3L))
      for (bins in bin_sets) {
        data$bins <- bins
        active <- which(bins > 0)
        expected <- gaussian
        if (length(active) && two_sided) {
          eligible_mean <- mean[active]
          eligible_covariance <- candidate[active, active, drop = FALSE]
          below <- vapply(two_sided_limits, function(limit) rectangle_mass(
            rep(-limit, length(active)), rep(limit, length(active)),
            eligible_mean, eligible_covariance), numeric(1L))
          normalizer <- sum(diff(c(0, below, 1)) * two_sided_omega)
          numerator <- two_sided_omega[findInterval(max(abs(data$y[active])), two_sided_limits) + 1L]
          expected <- expected + log(numerator / normalizer)
        } else if (length(active)) {
          weight <- function(cell) {
            if (rule == 1L) data$omega[min(cell)] else prod(data$omega[cell])
          }
          cells <- expand.grid(rep(list(seq_along(data$omega)), length(active)))
          normalizer <- sum(apply(cells, 1L, function(cell) {
            weight(cell) * cell_mass(cell, mean[active], candidate[active, active, drop = FALSE])
          }))
          numerator <- weight(bins[active])
          expected <- expected + log(numerator / normalizer)
        }
        syntax <- paste0("model { y[1:3] ~ ", spec$distribution, "(",
          paste(spec$arguments, collapse = ","), ") }")
        con <- textConnection(syntax)
        model <- tryCatch(rjags::jags.model(con,
          data = data[unique(c("y", spec$arguments))], n.chains = 1L,
          n.adapt = 0L, quiet = TRUE), finally = close(con))
        deviance <- rjags::jags.samples(model, "y", stat = "deviance",
          summary = "mean", n.iter = 1L, progress.bar = "none")
        ll <- -as.numeric(deviance[["y"]]) / 2
        expect_equal(ll, expected, tolerance = if (length(active) <= 1L) 1e-11 else 1e-4,
          label = paste(kind, rule, paste(bins, collapse = ",")))
      }
    }
  }
})

test_that("conditioning cancels selection at a fixed eligible outcome's threshold", {

  skip_if_not_installed("rjags")
  .load_RoBMA_module()
  fixtures <- .native_selection_jags_fixtures()
  spec <- fixtures$cases$conditioned
  data <- utils::modifyList(fixtures$base, spec$extra)
  data$y <- c(1.2, 0)
  data$mu <- c(-.2, .37)
  data$covariance <- c(.8, .32, 1.1)
  data$residual <- c(.4, 0)
  data$loading <- c(.25, 0)
  data$bins <- c(0L, 1L)
  data$groups <- c(1L, 1L)
  data$mode <- 1L
  # Only the exempt row has integrated variance. The eligible observation is
  # exactly on the bin boundary, irrespective of the retained auxiliary draw.
  data$sampling_aux <- c(.45, -.17)
  data$random_aux <- c(.21, 0)
  spec$arguments[6:7] <- c("sampling_aux", "random_aux")
  total <- matrix(c(.8 + .4 + .25^2, .32, .32, 1.1), 2L)
  expected <- mvtnorm::dmvnorm(data$y, data$mu, total, log = TRUE)
  syntax <- paste0("model { y[1:2] ~ ", spec$distribution, "(",
    paste(spec$arguments, collapse = ","), ") }")
  for (rule in 0:1) {
    data$rule <- rule
    con <- textConnection(syntax)
    model <- tryCatch(rjags::jags.model(con,
      data = data[unique(c("y", spec$arguments))], n.chains = 1L,
      n.adapt = 0L, quiet = TRUE), finally = close(con))
    deviance <- rjags::jags.samples(model, "y", stat = "deviance",
      summary = "mean", n.iter = 1L, progress.bar = "none")
    ll <- -as.numeric(deviance[["y"]]) / 2
    expect_equal(ll, expected, tolerance = 1e-12)
  }
})
