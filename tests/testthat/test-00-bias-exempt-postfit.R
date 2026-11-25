test_that("an all-FALSE bias flag preserves the complete univariate and multivariate JAGS graphs", {

  dat <- data.frame(yi = c(.1, .4, -.2, .3), sei = c(.4, .3, .5, .4),
                    registered = rep(FALSE, 4L))
  for (mv in c(FALSE, TRUE)) {
    constructor <- if (mv) RoBMA.mv else RoBMA
    args <- list(yi = quote(yi), data = dat, measure = "GEN",
                 prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
    if (mv) {
      args$V <- diag(dat$sei^2)
      args$V[1, 2] <- args$V[2, 1] <- .06
    } else {
      args$sei <- quote(sei)
    }
    before <- do.call(constructor, args)
    after <- do.call(constructor, c(args, list(bias_exempt = quote(registered))))
    expect_null(after$data$outcome$bias_exempt)
    expect_identical(.create_fit_data(after$data, after$priors),
                     .create_fit_data(before$data, before$priors))
    expect_identical(.create_model_syntax(after$data, after$priors),
                     .create_model_syntax(before$data, before$priors))
  }
})


test_that("PET and PEESE use the fitted and explicitly predicted row exemptions", {

  dat <- data.frame(yi = c(.1, -.2, .4), sei = c(.5, .7, .9),
                    registered = c(FALSE, TRUE, FALSE))
  newdata <- data.frame(sei = rev(dat$sei), bias_exempt = !dat$registered)
  slopes <- c(.4, .9, .6)
  for (mv in c(FALSE, TRUE)) for (parameter in c("PET", "PEESE")) {
    constructor <- get(paste0("b", parameter, if (mv) ".mv"))
    specification <- list(yi = quote(yi), data = dat, bias_exempt = quote(registered),
      measure = "GEN", prior_unit_information_sd = 1, effect_direction = "negative",
      only_priors = TRUE, silent = TRUE)
    if (mv) {
      specification$V <- diag(dat$sei^2)
    } else {
      specification$sei <- quote(sei)
    }
    object <- do.call(constructor, specification)
    samples <- cbind(mu = c(.2, -.1, .4), tau = .3)
    samples <- cbind(samples, slopes)
    colnames(samples)[ncol(samples)] <- parameter
    without_flags <- newdata["sei"]
    expect_error(predict(object, newdata = without_flags, type = "terms", quiet = TRUE,
                         .posterior_samples = samples), "bias_exempt")
    corrected <- predict(object, newdata = without_flags, type = "terms",
                         bias_adjusted = TRUE, quiet = TRUE, .posterior_samples = samples)
    expect_equal(as.numeric(corrected), rep(samples[, "mu"], nrow(without_flags)))
    for (explicit in c(FALSE, TRUE)) {
      args <- list(object = object, type = "terms", quiet = TRUE,
                   .posterior_samples = samples)
      if (explicit) args$newdata <- newdata
      adjusted <- do.call(predict, c(args, list(bias_adjusted = TRUE)))
      observed <- do.call(predict, args)
      sei <- if (explicit) newdata$sei else dat$sei
      exempt <- if (explicit) newdata$bias_exempt else dat$registered
      offset <- if (parameter == "PET") -sei else -sei^2
      offset[exempt] <- 0
      expect_equal(as.numeric(observed - adjusted), as.numeric(outer(slopes, offset)),
                   tolerance = 1e-12)
    }
  }
})

test_that("publication exemptions preserve univariate observation weights and source roles", {

  dat <- data.frame(yi = c(.14, .6, -.23), sei = c(.3, .55, .42), power = c(2, .5, 3))
  samples <- cbind(mu = c(.1, -.2), tau = c(.4, .3),
                   "theta[1]" = c(-.6, .3), "theta[2]" = c(.8, .4), "theta[3]" = c(.1, -.7))
  for (mode in c("integrate", "condition")) for (exempt in list(c(FALSE, TRUE, FALSE), rep(TRUE, 3L))) {
    object <- bselmodel(yi = yi, sei = sei, weights = power, data = dat, bias_exempt = exempt,
      measure = "GEN", prior_unit_information_sd = 1, effect_direction = "negative",
      prior_bias = BayesTools::prior_weightfunction("one-sided", .5, BayesTools::wf_fixed(c(1, .2)),
        model = selection_model(estimate_random_effects = mode)), only_priors = TRUE, silent = TRUE)
    actual <- .log_lik_from_posterior_samples(object$fit, samples, object$data, object$priors)
    expected <- matrix(NA_real_, nrow(samples), nrow(dat))
    for (j in seq_len(nrow(dat))) {
      mu <- samples[, "mu"]
      sd <- dat$sei[j]
      if (mode == "integrate") {
        sd <- sqrt(sd^2 + samples[, "tau"]^2)
      } else {
        mu <- mu + samples[, paste0("theta[", j, "]")] * samples[, "tau"]
      }
      score <- stats::dnorm(dat$yi[j], mu, sd, log = TRUE)
      if (!exempt[j]) {
        score <- score + log(if (dat$yi[j] <= 0) 1 else .2) - log(.2 + .8 * stats::pnorm(-mu / sd))
      }
      expected[, j] <- dat$power[j] * score
    }
    expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-12)
  }
})

test_that("deleting an outcome retains the eligible publication event", {

  # The protected outcome has the smallest p-value. Of the two eligible
  # outcomes only the second row is positive, separating product from best.
  dat <- data.frame(yi = c(1.2, .3, -.38, -.8), paper = c("a", "a", "a", "b"),
                    registered = c(TRUE, FALSE, FALSE, TRUE))
  V <- diag(c(.25, .36, .49, .25))
  V[1:3, 1:3] <- matrix(c(.25, .08, -.06, .08, .36, .12, -.06, .12, .49), 3L)
  samples <- cbind(mu = c(.13, -.08), "omega[1]" = 1, "omega[2]" = .2)
  for (rule in c("product", "best")) {
    object <- bselmodel.mv(yi = yi, V = V, data = dat, bias_exempt = registered,
      measure = "GEN", prior_unit_information_sd = 1, effect_direction = "positive",
      prior_bias = BayesTools::prior_weightfunction("one-sided", .5,
        BayesTools::wf_fixed(c(1, .2)), model = selection_model(weight_rule = rule,
          group = "paper")), only_priors = TRUE, silent = TRUE)
    actual <- .log_lik_from_posterior_samples(object$fit, samples, object$data, object$priors)
    expected <- matrix(NA_real_, nrow(samples), nrow(dat))
    for (j in seq_len(nrow(dat))) {
      # Conditional Gaussian before publication. Fixing a positive eligible
      # sibling makes the best-rule selection function constant for this row.
      regression <- drop(solve(V[-j, -j], V[-j, j]))
      sd <- sqrt(V[j, j] - sum(regression * V[-j, j]))
      for (s in seq_len(nrow(samples))) {
        mu <- samples[s, "mu"]
        mean <- mu + sum(regression * (dat$yi[-j] - mu))
        score <- stats::dnorm(dat$yi[j], mean, sd, log = TRUE)
        siblings <- dat$paper == dat$paper[j] & seq_len(nrow(dat)) != j & !dat$registered
        if (!dat$registered[j] && (rule == "product" || !any(dat$yi[siblings] >= 0))) {
          mass <- .2 + .8 * stats::pnorm(mean / sd)
          score <- score + log(if (dat$yi[j] >= 0) 1 else .2) - log(mass)
        }
        expected[s, j] <- score
      }
    }
    expect_equal(as.numeric(actual), as.numeric(expected), tolerance = 1e-8)
  }
})

test_that("a protected response retains its covariance with a selected sibling", {

  dat <- data.frame(yi = c(.3, -.2), paper = c("a", "a"), registered = c(FALSE, TRUE))
  V <- matrix(c(1, .8, .8, 1), 2L)
  object <- bselmodel.mv(yi = yi, V = V, data = dat, bias_exempt = registered,
    measure = "GEN", prior_unit_information_sd = 1, effect_direction = "positive",
    prior_bias = BayesTools::prior_weightfunction("one-sided", .5,
      BayesTools::wf_fixed(c(1, 0)), model = selection_model(weight_rule = "best",
        group = "paper")), only_priors = TRUE, silent = TRUE)
  n <- 1500L
  samples <- cbind(mu = rep(0, n), "omega[1]" = 1, "omega[2]" = 0)
  withr::local_seed(4911)
  expect_error(predict(object, newdata = data.frame(paper = dat$paper),
    V_new = V, type = "response", quiet = TRUE, .posterior_samples = samples), "bias_exempt")
  for (explicit in c(FALSE, TRUE)) {
    args <- list(object = object, type = "response", quiet = TRUE, .posterior_samples = samples)
    eligible <- 1L
    if (explicit) {
      eligible <- 2L
      args$newdata <- data.frame(paper = dat$paper, sei = sqrt(diag(V)),
                                 bias_exempt = c(TRUE, FALSE))
      args$V_new <- V
    }
    y <- as.matrix(do.call(predict, args))
    # The eligible response is half normal. Conditional on it, the other
    # response is N(.8*y, 1-.8^2), including when the fitted flags are reversed.
    expect_true(all(y[, eligible] >= 0))
    expect_lt(abs(mean(y[, eligible]) - sqrt(2 / pi)), 6 * sqrt((1 - 2 / pi) / n))
    residual <- y[, -eligible] - .8 * y[, eligible]
    expect_lt(abs(mean(residual)), 6 * .6 / sqrt(n))
    expect_lt(abs(stats::var(residual) - .36), 6 * .36 * sqrt(2 / (n - 1)))
  }
})
