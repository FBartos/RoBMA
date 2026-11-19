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
  expect_false(inherits(fit$priors$mods$w, "prior.regularized_term"))
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

test_that("factor contrasts and scale multipliers survive regularization group assignment", {
  dat <- data.frame(y = c(.2, .1, -.1, .4, .2, -.3), se = .2,
                    f = factor(rep(c("a", "b", "c"), 2)))
  global <- prior("uniform", list(.2, .6))
  build <- function(group) {
    BMA(yi = y, sei = se, data = dat, mods = ~ f, measure = "SMD",
        prior_mods = group, rescale_priors = 2, only_priors = TRUE)
  }
  expect_error(build(prior_regularized("horseshoe", scale = global)),
               "requires treatment or independent contrasts", fixed = TRUE)
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
