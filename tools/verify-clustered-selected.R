#!/usr/bin/env Rscript
# Compare a selected-normal grouped update with deterministic quadrature.
# Run with the installed development packages and JAGS 5.
args <- commandArgs(trailingOnly = TRUE)
n_draws <- if (length(args) >= 1L) as.integer(args[1L]) else 1500L
n_seeds <- if (length(args) >= 2L) as.integer(args[2L]) else 12L
stopifnot(n_draws >= 100L, n_seeds >= 8L)
stopifnot(requireNamespace("RoBMA"), requireNamespace("rjags"))

data <- list(
  y = c(.18, .42, -.12, .64, .22), se = c(.14, .18, .22, .27, .31),
  x = c(-1, 1, -1, 1, -1), study = c(1L, 1L, 2L, 2L, 2L),
  tau = .48, rho = .45, w = c(1, .42), cut = stats::qnorm(.95),
  zlo = c(stats::qnorm(.95), -1e300), zhi = c(1e300, stats::qnorm(.95)),
  location = c(.1, -.2), precision = 1 / c(.75, .65)^2, odds = c(.9, 1.1)
)
data$bin <- 1L + as.integer(data$y / data$se < data$cut)

hermite <- function(n) {
  k <- seq_len(n - 1L)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(k, k + 1L)] <- jacobi[cbind(k + 1L, k)] <- sqrt(k)
  e <- eigen(jacobi, symmetric = TRUE)
  list(z = e$values, w = e$vectors[1L, ]^2)
}

reference <- function(n) {
  grid <- hermite(n)
  variance <- data$se^2 + data$tau^2 * (1 - data$rho)
  between <- data$tau * sqrt(data$rho)
  studies <- outer(data$study, 1:2, "==")
  result <- vapply(1:2, function(label) {
    included <- label == 2L
    design <- cbind(1, if (included) data$x, studies * between)
    center <- c(data$location[1L], if (included) data$location[2L], 0, 0)
    precision <- c(data$precision[1L], if (included) data$precision[2L], 1, 1)

    # Use the exact joint Gaussian posterior as the quadrature measure.
    # The only remaining factor is the product of row selection corrections.
    covariance_y <- diag(variance) + tcrossprod(sweep(design, 2L, sqrt(precision), "/"))
    root_y <- chol(covariance_y)
    z <- forwardsolve(t(root_y), data$y - drop(design %*% center))
    gaussian_evidence <- exp(-.5 * sum(z^2) - sum(log(diag(root_y))))
    posterior_covariance <- chol2inv(chol(crossprod(design, design / variance) + diag(precision)))
    posterior_mean <- posterior_covariance %*% (crossprod(design, data$y / variance) + precision * center)
    index <- as.matrix(expand.grid(rep(list(seq_len(n)), ncol(design))))
    latent <- sweep(matrix(grid$z[index], nrow(index)) %*% chol(posterior_covariance),
                    2L, drop(posterior_mean), "+")
    weights <- apply(matrix(grid$w[index], nrow(index)), 1L, prod)
    means <- latent %*% t(design)
    for (i in seq_along(data$y)) {
      above <- stats::pnorm(data$cut * data$se[i], means[, i], sqrt(variance[i]), lower.tail = FALSE)
      weights <- weights * data$w[data$bin[i]] / (above + (1 - above) * data$w[2L])
    }
    interpolated <- if (included) latent[, 2L] else 0
    data$odds[label] * gaussian_evidence * c(sum(weights), sum(weights * latent[, 1L]),
                                           sum(weights * interpolated))
  }, numeric(3L))
  mass <- sum(result[1L, ])
  c(inclusion = result[1L, 2L] / mass, intercept = sum(result[2L, ]) / mass,
    slope = sum(result[3L, ]) / mass)
}

model <- function(tags) {
  declaration <- if (tags) "
    m1 <- robma_cluster_member(g, 1, 1)
    m2 <- robma_cluster_member(a, 2, 1)
    m3 <- robma_cluster_member(b, 3, 1)
    m4 <- robma_cluster_member(u[1], 4, 1)
    m5 <- robma_cluster_member(u[2], 4, 2)
    setup <- robma_cluster_setup(g, x, study, y, se, variance, bin, zlo, zhi,
                                 w, 1, 1, 1, 1, 2, tau, rho, 2, 2, 0, 0, 2, 5, 0)
  " else ""
  paste0("model {
    g ~ dcat(odds)
    a ~ dnorm(location[1], precision[1])
    b ~ dnorm(location[2], precision[2])
    for (j in 1:2) { u[j] ~ dnorm(0, 1) }
    for (i in 1:5) {
      mu[i] <- a + equals(g, 2) * b * x[i] + tau * sqrt(rho) * u[study[i]]
      total[i] <- sqrt(pow(se[i], 2) + pow(tau, 2) * (1 - rho))
      y[i] ~ dselnorm_step_switch(mu[i], total[i], se[i], 1, w, zlo, zhi, bin[i], 1, 1, 1)
    }
    ", declaration, "
  }")
}

run <- function(tags, seed) {
  con <- textConnection(model(tags))
  on.exit(close(con))
  inputs <- data[setdiff(names(data), "cut")]
  if (tags) inputs$variance <- data$se^2
  fit <- rjags::jags.model(con, inputs,
    inits = list(g = 2L, a = .1, b = -.2, u = c(.05, -.05),
                 .RNG.name = "base::Mersenne-Twister", .RNG.seed = seed),
    n.chains = 1L, n.threads = 1L, n.adapt = 0L, quiet = TRUE)
  owned <- rjags::list.samplers(fit)[["RoBMA::JointCluster"]]
  stopifnot(if (tags) setequal(owned, c("g", "a", "b", "u[1]", "u[2]")) else is.null(owned))
  rjags::adapt(fit, n.iter = 100L, end.adaptation = TRUE, progress.bar = "none")
  stats::update(fit, n.iter = 100L, progress.bar = "none")
  draws <- as.matrix(rjags::coda.samples(fit, c("a", "b", "g"), n.iter = n_draws,
                                          progress.bar = "none"))
  active <- draws[, "g"] == 2
  c(inclusion = mean(active), intercept = mean(draws[, "a"]), slope = mean(draws[, "b"] * active))
}

exact <- reference(20L)
ref_difference <- max(abs(exact - reference(14L)))
stopifnot(ref_difference < 5e-7)
cat("Maximum quadrature difference", ref_difference, "\n")
for (tags in c(FALSE, TRUE)) {
  estimates <- t(vapply(seq.int(801L, length.out = n_seeds), function(seed) {
    run(tags, seed)
  }, numeric(3L)))
  mcse <- apply(estimates, 2L, stats::sd) / sqrt(n_seeds)
  z <- (colMeans(estimates) - exact) / mcse
  print(data.frame(tagged = tags, target = names(exact), exact, observed = colMeans(estimates),
                   mcse, z), row.names = FALSE)
  stopifnot(all(abs(z) < 5))
}
