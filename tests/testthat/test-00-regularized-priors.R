probe_state("test-00-regularized-priors.R")
test_that("moderator and scale groups survive both passes of formula binding", {
  dat <- data.frame(y = c(0.1, -0.3, 0.4, 0.2, 0.5, -0.2), se = 0.2,
                    x = 1:6, z = c(1, 3, 2, 6, 5, 4), w = c(1, 0, 1, 0, 0, 1))
  global <- prior("cauchy", list(0, 0.2), list(0, Inf))
  group <- prior_regularized("horseshoe", scale = global, slab_scale = 1)
  fit <- brma(yi = y, sei = se, data = dat, mods = ~ x + z + w, scale = ~ x + z,
              measure = "SMD", prior_mods = list(pool = group, w = prior("normal", list(0, 1))),
              prior_scale = list(first = prior_regularized("ridge", terms = "x", scale = global),
                                 second = prior_regularized("lasso", terms = "z", scale = global)),
              only_priors = TRUE, silent = TRUE)
  expect_s3_class(fit$priors$mods$x, "prior.regularized_term")
  expect_s3_class(fit$priors$mods$z, "prior.regularized_term")
  expect_false(BayesTools::is.prior.regularized_term(fit$priors$mods$w))
  expect_identical(fit$priors$mods$w$parameters$sd, 1)
  args <- .create_jags_formula_args(fit$data, fit$priors)
  for (component in c("mu", "log_tau")) {
    rebound <- BayesTools::JAGS_formula(
      args$formula_list[[component]], component, args$formula_data_list[[component]],
      prior_list = args$formula_prior_list[[component]],
      formula_scale = args$formula_scale_list[[component]])
    for (term in c("x", "z")) {
      expect_s3_class(rebound$prior_list[[paste(component, term, sep = "_")]],
                      "prior.regularized_term")
    }
  }
})

test_that("regularized moderator alternatives retain explicit null weights and omissions", {
  dat <- data.frame(y = c(0.1, 0.2, -0.2, 0.3, -0.1, 0.4), se = 0.1,
                    x = 1:6, z = c(3, 2, 5, 1, 6, 4))
  global <- prior("cauchy", list(0, 0.2), list(0, Inf))
  fit <- BMA(yi = y, sei = se, data = dat, mods = ~ x + z, measure = "SMD",
             prior_mods = list(pool = prior_regularized("horseshoe", scale = global)),
             prior_mods_null = list(x = prior("point", list(0), prior_weights = 3), z = NULL),
             only_priors = TRUE, silent = TRUE)
  mixture <- fit$priors$mods$x
  expect_true(BayesTools::is.prior.mixture(mixture))
  expect_equal(unname(attr(mixture, "prior_weights")), c(3, 1))
  expect_s3_class(mixture[[2]], "prior.regularized_term")
  expect_length(fit$priors$mods$z, 1)
  expect_s3_class(fit$priors$mods$z[[1]], "prior.regularized_term")
  args <- .create_jags_formula_args(fit$data, fit$priors)
  expect_silent(BayesTools::JAGS_formula(
    args$formula_list$mu, "mu", args$formula_data_list$mu,
    prior_list = args$formula_prior_list$mu, formula_scale = args$formula_scale_list$mu))
})

test_that("fixed lasso rescales its Laplace coefficient scale", {
  p <- prior_regularized("lasso", scale = 0.25)
  expect_identical(.rescale_prior_object(p, 2)$parameters$scale, 0.5)
  expect_equal(.rescale_prior_object(prior("laplace", list(1, 0.3)), 2)$parameters,
               list(location = 1, scale = 0.6))
})

test_that("density updates of regularized slopes condition on each row's shared scale", {
  global <- prior("uniform", list(0.1, 1))
  prior_list <- BayesTools::JAGS_formula(
    ~ 0 + x + z, "mu", data.frame(x = c(0, 1, 3), z = c(1, 2, 0)),
    list(pool = prior_regularized("ridge", scale = global)))$prior_list
  expect_false(.iwmde_can_use_focal_prior_delta(prior_list$mu_x))
  samples <- cbind(mu_intercept = 0, mu_x = 0.1, mu_z = -0.2, mu_pool_global = c(0.2, 0.8))
  proposed <- samples
  proposed[, "mu_x"] <- 0.4
  log_density <- .iwmde_numerical_prior_rows_evaluator(prior_list)
  expect_equal(log_density(proposed) - log_density(samples),
               dnorm(0.4, 0, c(0.2, 0.8), log = TRUE) - dnorm(0.1, 0, c(0.2, 0.8), log = TRUE))
})

test_that("factor contrasts and scale multipliers survive regularization group assignment", {
  dat <- data.frame(y = c(.2, .1, -.1, .4, .2, -.3), se = .2,
                    f = factor(rep(c("a", "b", "c"), 2)))
  global <- prior("uniform", list(.2, .6))
  build <- function(group) {
    BMA(yi = y, sei = se, data = dat, mods = ~ f, measure = "SMD",
        prior_mods = group, rescale_priors = 2, only_priors = TRUE)
  }
  default <- build(prior_regularized("horseshoe", scale = global))$priors$mods$f
  expect_true(BayesTools::is.prior.meandif(default[[1L]]))
  expect_true(BayesTools::is.prior.meandif(default[[2L]]))
  alternative <- build(prior_regularized("horseshoe", scale = global,
                           slab_scale = 1, contrast = "treatment"))$priors$mods$f[[2L]]
  expect_true(BayesTools::is.prior.treatment(alternative))
  expect_equal(alternative$regularization$multiplier, 2)
  expect_identical(alternative$regularization$scale, global)
  expect_equal(alternative$regularization$slab_scale, 1)
  alternative <- build(prior_regularized("lasso", scale = .25,
                                         contrast = "independent"))$priors$mods$f[[2L]]
  expect_true(BayesTools::is.prior.independent(alternative))
  expect_equal(alternative$parameters$scale, .5)
})

test_that("centered factor priors inherit or override the location and scale contrast", {
  dat <- data.frame(y = c(.2, .1, -.1, .4, .2, -.3), se = .2,
                    f = factor(rep(c("a", "b", "c"), 2L)), x = c(-1, 2, 0, 2, 1, 0))
  global <- prior("uniform", list(.2, .6))
  for (coding in c("meandif", "orthonormal")) {
    other <- if (coding == "meandif") "orthonormal" else "meandif"
    fit <- BMA(yi = y, sei = se, data = dat, mods = ~ f + x, scale = ~ f,
               prior_mods = list(pool = prior_regularized("horseshoe", scale = global)),
               prior_scale = list(f = prior_regularized("lasso", scale = .3, contrast = other)),
               prior_heterogeneity = prior("point", list(.2)), prior_heterogeneity_null = NULL,
               set_contrast_factor_predictors = coding, measure = "SMD", only_priors = TRUE)
    for (part in c("mods", "scale")) {
      components <- fit$priors[[part]]$f
      expected <- if (part == "mods") coding else other
      expect_identical(.prior_factor_contrast(components[[1L]]), expected)
      expect_identical(.prior_factor_contrast(components[[2L]]), expected)
      expect_identical(unname(attr(components, "components")), c("null", "alternative"))
    }
    expect_s3_class(fit$priors$mods$x[[2L]], "prior.regularized_term")
    fixed <- brma(yi = y, sei = se, data = dat, mods = ~ f, measure = "SMD",
                  prior_mods = list(f = prior_regularized("lasso", scale = .3)),
                  set_contrast_factor_predictors = coding, only_priors = TRUE)
    expect_identical(.prior_factor_contrast(fixed$priors$mods$f), coding)
  }
})

test_that("a pooled interaction preserves the fitted bases of its two factors", {
  dat <- expand.grid(f = factor(c("A", "B", "C")), h = factor(c("U", "V", "W")))
  dat$y <- seq_len(nrow(dat)) / nrow(dat)
  dat$se <- .2
  fit <- BMA(yi = y, sei = se, mods = ~ f * h, data = dat, measure = "SMD",
             prior_mods = list(f = prior_factor("mnormal", list(0, 1), contrast = "orthonormal"),
                               h = prior_factor("normal", list(0, 1), contrast = "treatment"),
                               pool = prior_regularized("lasso", terms = "f:h", scale = .4)),
             only_priors = TRUE)
  alternative <- fit$priors$mods[["f:h"]][[2L]]
  expect_identical(attr(alternative, "factor_contrasts"),
                   c(f = "contr.orthonormal", h = "contr.treatment"))
  expect_equal(alternative$regularization$unit_index, c(1, 1, 2, 2))
  expect_equal(alternative$regularization$unit_dimension, c(2, 2))
  expect_length(fit$priors$mods[["f:h"]], 2L)
  expect_identical(.prior_factor_contrast(fit$priors$mods[["f:h"]][[1L]]),
                   .prior_factor_contrast(alternative))
  expect_true(BayesTools::is.prior.orthonormal(alternative))

  for (run in list(brma, BMA)) {
    condition <- tryCatch(
      run(yi = y, sei = se, mods = ~ f * h, data = dat, measure = "SMD",
          prior_mods = list(pool = prior_regularized("lasso", terms = "f:h", scale = .4,
                                                    contrast = "treatment")),
          set_contrast_factor_predictors = "meandif", only_priors = TRUE),
      error = identity)
    expect_s3_class(condition, "BayesTools_formula_transform_unavailable")
    expect_identical(condition$reason, "regularized_contrast_conflict")
  }
  expect_error(BMA(yi = y, sei = se, mods = ~ f, data = dat, measure = "SMD",
                   prior_mods = prior("laplace", list(0, .4)), only_priors = TRUE),
               'prior_regularized("lasso", scale = ...)', fixed = TRUE)
})

test_that("a shrinkage group in the null list supplies its factor and slope nulls", {
  dat <- data.frame(y = c(.2, .1, -.1, .4, .2, -.3), se = .2,
                    f = factor(rep(c("a", "b", "c"), 2L)), x = c(-1, 2, 0, 2, 1, 0))
  wide   <- prior("uniform", list(.3, .9))
  narrow <- prior("uniform", list(.03, .09))
  build <- function(prior_mods, prior_mods_null) {
    BMA(yi = y, sei = se, data = dat, mods = ~ f + x, measure = "SMD",
        prior_mods = prior_mods, prior_mods_null = prior_mods_null, only_priors = TRUE)
  }
  null_group <- function(fit, term) fit$priors$mods[[term]][[1L]]$regularization$group

  # a group without 'terms' supplies every null not named in the list
  fit <- build(list(slab = prior_regularized("ridge", wide)),
               list(nullpool = prior_regularized("ridge", narrow)))
  for (term in c("f", "x")) {
    expect_identical(null_group(fit, term), "nullpool", info = term)
    expect_identical(fit$priors$mods[[term]][[2L]]$regularization$group, "slab", info = term)
  }
  # a group naming the factor supplies only its null; the slope keeps the spike
  fit <- build(list(slab = prior_regularized("ridge", wide)),
               list(nullpool = prior_regularized("ridge", narrow, terms = "f")))
  expect_identical(null_group(fit, "f"), "nullpool")
  expect_true(BayesTools::is.prior.point(fit$priors$mods$x[[1L]]))

  # bare groups on both sides, here horseshoes with local scales, stay apart
  fit <- build(prior_regularized("horseshoe", wide), prior_regularized("horseshoe", narrow))
  expect_identical(null_group(fit, "x"), "shrinkage_null")
  expect_identical(fit$priors$mods$x[[2L]]$regularization$group, "shrinkage")
  args <- .create_jags_formula_args(fit$data, fit$priors)
  rebound <- BayesTools::JAGS_formula(
    args$formula_list$mu, "mu", args$formula_data_list$mu,
    prior_list = args$formula_prior_list$mu, formula_scale = args$formula_scale_list$mu)
  locals <- unlist(lapply(rebound$prior_list$mu_x, function(p) p$regularization$local_names))
  expect_identical(locals, c("mu_x_component_1_local_1", "mu_x_component_2_local_1"))
  expect_true(all(c("mu_shrinkage_global", "mu_shrinkage_null_global") %in%
                    BayesTools::JAGS_to_monitor(rebound$prior_list)))

  # generated names avoid the names given on the other side
  groups_of <- function(fit) {
    args <- .create_jags_formula_args(fit$data, fit$priors)
    rebound <- BayesTools::JAGS_formula(
      args$formula_list$mu, "mu", args$formula_data_list$mu,
      prior_list = args$formula_prior_list$mu, formula_scale = args$formula_scale_list$mu)
    vapply(rebound$prior_list$mu_x, function(p) p$regularization$group, character(1))
  }
  expect_identical(groups_of(build(list(shrinkage_null = prior_regularized("ridge", wide)),
                                   prior_regularized("ridge", narrow))),
                   c("shrinkage_null_1", "shrinkage_null"))
  expect_identical(groups_of(build(prior_regularized("ridge", wide),
                                   list(shrinkage = prior_regularized("ridge", narrow)))),
                   c("shrinkage", "shrinkage_1"))
})
