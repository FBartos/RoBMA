test_that("JAGS selection blocks preserve the full Gaussian law of exempt outcomes", {

  skip_if_not_installed("rjags")
  .load_RoBMA_module()
  rjags::load.module("dic", quiet = TRUE)
  fixtures <- .native_selection_jags_fixtures()
  cluster_rules <- .selection_joint_cluster_quadrature_rules(c(7L, 15L, 31L, 63L))
  factor_rules <- .selection_joint_factor_quadrature_rules()
  # Integrate the Gaussian in each eligible p-bin rectangle, conditioning the
  # second coordinate on the first. This reference uses no selection kernels.
  bounds <- c(Inf, 1, 0, -Inf)
  cell_mass <- function(bins, mean, covariance) {
    lo <- bounds[bins + 1L]
    hi <- bounds[bins]
    sd <- sqrt(diag(covariance))
    if (length(bins) == 1L) return(pnorm(hi, mean, sd) - pnorm(lo, mean, sd))
    slope <- covariance[2, 1] / covariance[1, 1]
    conditional_sd <- sqrt(covariance[2, 2] - slope * covariance[1, 2])
    integrate(function(x) {
      center <- mean[2] + slope * (x - mean[1])
      dnorm(x, mean[1], sd[1]) *
        (pnorm(hi[2], center, conditional_sd) - pnorm(lo[2], center, conditional_sd))
    }, lo[1], hi[1], rel.tol = 1e-10)$value
  }
  # Three correlated outcomes with the most significant result in the middle.
  sampling <- matrix(c(.8, .32, -.1, .32, 1.1, .16, -.1, .16, .93), 3L)
  for (kind in names(fixtures$cases)) {
    spec <- fixtures$cases[[kind]]
    data <- utils::modifyList(fixtures$base, spec$extra)
    data$y <- c(.4, 1.7, -.6)
    data$omega <- c(1, .2, .8)
    data$lower <- c(1, 0, -1e300)
    data$upper <- c(1e300, 1, 0)
    data$mu <- numeric(3L)
    data$sei <- rep(1, 3L)
    data$mode <- 1L
    data$covariance <- sampling[lower.tri(sampling, diag = TRUE)]
    data$loading <- c(.25, -.35, .18)
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
    gaussian <- mvtnorm::dmvnorm(data$y, sigma = total, log = TRUE)
    candidate <- if (kind == "conditioned") total - sampling else total
    # Retained sources condition on every observed outcome, including exempt rows.
    mean <- if (kind == "conditioned") drop(sampling %*% solve(total, data$y)) else data$mu
    for (rule in 0:1) {
      data$rule <- rule
      for (bins in list(c(2L, 0L, 3L), c(0L, 0L, 3L), rep(0L, 3L))) {
        data$bins <- bins
        active <- which(bins > 0)
        expected <- gaussian
        if (length(active)) {
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
        ll <- -as.numeric(rjags::jags.samples(model, "deviance", n.iter = 1L,
          progress.bar = "none")$deviance) / 2
        expect_equal(ll, expected, tolerance = if (length(active) <= 1L) 1e-11 else 1e-4,
          label = paste(kind, rule, paste(bins, collapse = ",")))
      }
    }
  }
})

test_that("conditioning cancels selection at a fixed eligible outcome's threshold", {

  skip_if_not_installed("rjags")
  .load_RoBMA_module()
  rjags::load.module("dic", quiet = TRUE)
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
    ll <- -as.numeric(rjags::jags.samples(model, "deviance", n.iter = 1L,
      progress.bar = "none")$deviance) / 2
    expect_equal(ll, expected, tolerance = 1e-12)
  }
})
