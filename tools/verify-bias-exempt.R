# Fits all eight publication-bias constructors. At these fixed slope and
# heterogeneity settings the posterior of mu has an independent 1D reference.
library(RoBMA)

requested <- commandArgs(trailingOnly = TRUE)
constructors <- c("bselmodel", "bPET", "bPEESE", "RoBMA")
stopifnot(all(requested %in% c(constructors, paste0(constructors, ".mv"))))

dat <- data.frame(
  yi = c(-.36, .52, .22, -.16, .45, -.32, .38, -.18),
  sei = c(.32, .43, .27, .38, .45, .37, .35, .42),
  paper = rep(letters[1:4], each = 2L),
  registered = c(TRUE, FALSE, FALSE, TRUE, TRUE, FALSE, TRUE, TRUE)
)
point <- function(value) BayesTools::prior("point", list(value))
slope_pet <- .8
slope_peese <- 2
lower_weight <- .18

posterior_reference <- function(V, methods) {

  # Each paper contains at most one eligible outcome. Its whole publication
  # normalizer is the univariate Gaussian mass of that outcome, even if the
  # protected sibling has correlated sampling error.
  active <- !dat$registered
  inverse <- solve(V)
  posterior_density <- function(mu, method) {
    offset <- switch(method, PET = slope_pet * dat$sei,
                     PEESE = slope_peese * dat$sei^2, rep(0, nrow(dat)))
    offset[!active] <- 0
    residual <- matrix(dat$yi - offset, nrow(dat), length(mu)) -
      rep(mu, each = nrow(dat))
    log_density <- dnorm(mu, log = TRUE) - colSums(residual * (inverse %*% residual)) / 2
    if (method == "selection") {
      log_weight <- sum(log(ifelse(dat$yi[active] >= 0, 1, lower_weight)))
      mass <- lower_weight + (1 - lower_weight) * pnorm(outer(1 / dat$sei[active], mu))
      log_density <- log_density + log_weight - colSums(log(mass))
    }
    exp(log_density)
  }
  integral <- function(method, power) {
    integrate(function(mu) mu^power * posterior_density(mu, method),
              -Inf, Inf, rel.tol = 1e-10)$value
  }
  evidence <- vapply(methods, integral, numeric(1L), power = 0L)
  first <- vapply(methods, integral, numeric(1L), power = 1L)
  list(mu = sum(first) / sum(evidence), probability = evidence / sum(evidence))
}

check_mcmc_mean <- function(draws, reference) {

  ess <- as.numeric(coda::effectiveSize(draws))
  mcse <- sd(draws) / sqrt(ess)
  stopifnot(is.finite(mcse), ess > 50,
            abs(mean(draws) - reference) <= 6 * mcse)
  cat(sprintf("  observed %.5f, reference %.5f, MCSE %.5f\n",
              mean(draws), reference, mcse))
}

for (mv in c(FALSE, TRUE)) {
  R <- if (mv) kronecker(diag(4L), matrix(c(1, .62, .62, 1), 2L)) else diag(8L)
  V <- R * outer(dat$sei, dat$sei)
  rule <- if (mv) "best" else "product"
  priors <- list(
    selection = BayesTools::prior_weightfunction("one-sided", .5,
      BayesTools::wf_fixed(c(1, lower_weight)),
      model = selection_model(weight_rule = rule, group = "paper")),
    PET = BayesTools::prior_PET("spike", list(location = slope_pet)),
    PEESE = BayesTools::prior_PEESE("spike", list(location = slope_peese))
  )
  for (method in c("selection", "PET", "PEESE", "mixture")) {
    name <- paste0(switch(method, selection = "bselmodel", PET = "bPET",
                          PEESE = "bPEESE", mixture = "RoBMA"), if (mv) ".mv")
    if (length(requested) && !name %in% requested) next
    cat(name, "\n")
    args <- list(yi = quote(yi), data = dat, bias_exempt = quote(registered),
      measure = "GEN", prior_unit_information_sd = 1, effect_direction = "positive",
      prior_effect = BayesTools::prior("normal", list(0, 1)),
      sample = 2400L, adapt = 500L, burnin = 250L, chains = 1L, seed = 91711L,
      silent = TRUE, convergence_checks = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL))
    if (mv) {
      args$V <- V
    } else {
      args$sei <- quote(sei)
      args$prior_heterogeneity <- point(0)
    }
    methods <- method
    if (method == "mixture") {
      methods <- c("none", names(priors))
      args$prior_bias <- unname(priors)
      args$prior_bias_null <- BayesTools::prior_none()
      args["prior_effect_null"] <- list(NULL)
      if (!mv) args["prior_heterogeneity_null"] <- list(NULL)
    } else {
      args$prior_bias <- priors[[method]]
    }
    fit <- do.call(getExportedValue("RoBMA", name), args)
    reference <- posterior_reference(V, methods)
    base <- as.matrix(predict(fit, type = "terms", bias_adjusted = TRUE, quiet = TRUE))
    check_mcmc_mean(base[, 1L], reference$mu)
    samples <- RoBMA:::.get_posterior_samples(fit$fit)
    branch <- if (method == "mixture") methods[samples[, "bias_indicator"]] else
      rep(method, nrow(base))
    if (method == "mixture") {
      for (j in seq_along(methods)) {
        check_mcmc_mean(as.numeric(branch == methods[j]), reference$probability[j])
      }
    }
    for (explicit in c(FALSE, TRUE)) {
      pred <- list(object = fit, type = "terms", quiet = TRUE)
      exempt <- dat$registered
      if (explicit) {
        exempt <- !exempt
        pred$newdata <- data.frame(sei = dat$sei, paper = dat$paper, bias_exempt = exempt)
      }
      observed <- as.matrix(do.call(predict, pred))
      offset <- outer((branch == "PET") * slope_pet, dat$sei) +
        outer((branch == "PEESE") * slope_peese, dat$sei^2)
      offset[, exempt] <- 0
      stopifnot(max(abs(observed - base - offset)) < 1e-11)
      pred$type <- "response"
      if (mv && explicit) pred$V_new <- V
      response <- as.matrix(do.call(predict, pred))
      stopifnot(identical(dim(response), dim(base)), all(is.finite(response)))
    }
    cat("  fitted and explicit-newdata predictions passed\n")
  }
}
