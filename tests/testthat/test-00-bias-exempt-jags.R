test_that("JAGS selection blocks preserve the full Gaussian law of exempt outcomes", {

  skip_if_not_installed("rjags")
  .load_RoBMA_module()
  rjags::load.module("dic", quiet = TRUE)
  fixtures <- .native_selection_jags_fixtures()
  for (kind in names(fixtures$cases)) {
    spec <- fixtures$cases[[kind]]
    data <- utils::modifyList(fixtures$base, spec$extra)
    data$y <- c(1.7, -.6)
    data$mode <- 1L
    data$covariance <- c(.8, .32, 1.1)
    sampling <- matrix(c(.8, .32, .32, 1.1), 2L)
    data$loading <- c(.25, -.35)
    data$sd <- c(.6, .9)
    data$factor_loading <- c(.25, -.35, .4, .3)
    data$residual <- c(.4, .5)
    data$groups <- c(1L, 1L)
    total <- switch(kind,
      dense = sampling,
      cluster = diag(data$sd^2) + tcrossprod(data$loading),
      factor = diag(data$sd^2) + tcrossprod(matrix(data$factor_loading, 2L)),
      conditioned = sampling + diag(data$residual) + tcrossprod(data$loading)
    )
    gaussian <- mvtnorm::dmvnorm(data$y, sigma = total, log = TRUE)
    acceptance <- .6 # P(Y_2 > 0) + .2 * P(Y_2 < 0) for a centered marginal.
    if (kind == "conditioned") {
      # The retained-source reconstruction conditions on BOTH observed outcomes,
      # including the exempt first row.
      retained <- drop(sampling %*% solve(total, data$y))
      candidate_sd <- sqrt(data$residual[2L] + data$loading[2L]^2)
      acceptance <- .2 + .8 * pnorm(retained[2L] / candidate_sd)
    }
    for (rule in 0:1) {
      data$rule <- rule
      for (all_exempt in c(FALSE, TRUE)) {
        data$bins <- if (all_exempt) c(0L, 0L) else c(0L, 2L)
        syntax <- paste0("model { y[1:2] ~ ", spec$distribution, "(",
          paste(spec$arguments, collapse = ","), ") }")
        con <- textConnection(syntax)
        model <- tryCatch(rjags::jags.model(con,
          data = data[unique(c("y", spec$arguments))], n.chains = 1L,
          n.adapt = 0L, quiet = TRUE), finally = close(con))
        ll <- -as.numeric(rjags::jags.samples(model, "deviance", n.iter = 1L,
          progress.bar = "none")$deviance) / 2
        expected <- gaussian + if (all_exempt) 0 else log(.2 / acceptance)
        expect_equal(ll, expected, tolerance = 1e-11,
          label = paste(kind, rule, all_exempt))
      }
    }
  }
})
