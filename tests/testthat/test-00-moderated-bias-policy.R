test_that("selection meta-regression binds the exact additive moderator offset on request", {

  skip_on_cran()
  previous <- RoBMA.options()
  on.exit(do.call(RoBMA.options, previous), add = TRUE)
  RoBMA.options(sampler.nested.weights = FALSE, sampler.probabilities = character(),
                selection.sampler = "slice.coarse")
  dat <- data.frame(
    yi = c(-.44, -.8, .1, -.3, .6, -.4, -.1, .2, -.25),
    vi = c(.1, .04, .2, .3, .1, .2, .05, .12, .08),
    group = factor(rep(c("a", "b", "c"), 3L)), age = c(31, 22, 36, 21, 33, 41, 29, 26, 30),
    exempt = c(TRUE, rep(FALSE, 8L))
  )
  fit <- function(policy, mods = ~ group + age) {
    RoBMA.options(sampler.joint = policy)
    RoBMA(yi = yi, vi = vi, data = dat, measure = "SMD", mods = mods,
      effect_direction = "negative", bias_exempt = exempt, seed = 707L,
      chains = 1L, adapt = 500L, burnin = 50L, sample = 100L,
      convergence_checks = set_convergence_checks(max_Rhat = NULL, min_ESS = NULL))[["fit"]]
  }
  claimed <- function(backend) {
    assignment <- backend[["samplers"]]
    assignment[["Node"]][assignment[["Sampler"]] == "RoBMA::IndependentBiasJoint"]
  }

  fitted <- fit("bias")
  declarations <- attr(fitted, "sampler_declarations", exact = TRUE)
  marker <- declarations[["joint"]][["markers"]][[1L]]
  expect_setequal(claimed(fitted), c("bias_indicator", "mu_intercept_indicator", "mu_intercept_component_2"))
  expect_identical(declarations[["robma"]][["joint"]][["status"]], "active")
  expect_identical(marker[["tag"]], "robma_joint_bias_moderated")
  expect_identical(marker[["arguments"]][[5L]], -1)
  expect_identical(marker[["arguments"]][[10L]], c("mu_group[1]", "mu_group[2]", "mu_age"))
  data <- fitted[["data"]]
  emitted <- cbind(data[["mu_data_group"]], data[["mu_data_age"]])
  expect_identical(data[["RoBMA_joint_moderator_design"]],
                   as.numeric(t(emitted[data[["sel_joint_singleton_row"]], , drop = FALSE])))
  expect_true(any(data[["sel_joint_singleton_obs_bin"]] == 0L))

  expect_length(claimed(fit("auto")), 0L)
  interaction <- fit("bias", mods = ~ group * age)
  expect_length(claimed(interaction), 0L)
  expect_null(attr(interaction, "sampler_declarations", exact = TRUE)[["joint"]])
})
