test_that("F2-1 joint selection evaluates the compiled known-R diagonal law", {

  dat <- data.frame(yi = c(.35, -.1, .6, .2), esid = factor(paste0("e", 1:4)))
  se <- c(.2, .3, .25, .22)
  diagonal <- c(.5, 1, 2, 4)
  R <- diag(diagonal); dimnames(R) <- list(levels(dat$esid), levels(dat$esid))
  for (scale in c("none", "cor", "no-R")) {
    args <- list(yi = quote(yi), V = diag(se^2), random = ~1 | esid, data = dat,
      measure = "GEN", prior_unit_information_sd = 1, effect_direction = "positive",
      prior_bias = BayesTools::prior_weightfunction("one-sided", .025,
        BayesTools::wf_fixed(c(1, .4))), only_priors = TRUE, silent = TRUE)
    if (scale != "no-R") args <- c(args, list(R = list(esid = R), Rscale = scale))
    compiled <- review_compile_object(do.call(bselmodel.mv, args))
    term <- compiled$object$formula_design$mu$random_effects[[1L]]
    samples <- cbind(mu_intercept = c(.1, -.05, .2), c(.15, .3, .25))
    colnames(samples)[2L] <- term$sd_parameter_names
    initial <- review_attach_replay(compiled, samples)
    generated <- BayesTools::JAGS_evaluate_deterministic(initial$fit, draws = samples)
    samples <- cbind(samples, generated[, setdiff(colnames(generated), colnames(samples)), drop = FALSE])
    object <- review_attach_replay(compiled, samples)
    setup <- .log_lik_posterior_setup(object$fit, samples, object$data, object$priors, "estimate", NULL)
    setup$S <- 3L; setup$K <- 4L
    plan <- .data_selection_execution_plan(object$data)
    multiplier <- if (scale == "none") diagonal else rep(1, 4)
    actual <- .selection_joint_random_factor_samples(setup)
    expect_lte(max(abs(actual$diagonal - outer(samples[, 2L]^2, multiplier))), 1e-12)
    if (scale != "no-R") expect_null(.selection_joint_random_diagonal_samples(setup, plan))
    observed <- .selection_joint_block_loglik_from_setup(setup)
    reference <- matrix(0, 3L, 4L)
    for (s in 1:3) for (k in 1:4) {
      sd <- sqrt(se[k]^2 + samples[s, 2L]^2 * multiplier[k])
      cut <- qnorm(.975) * se[k]
      p <- pnorm(cut, samples[s, 1L], sd, lower.tail = FALSE)
      reference[s, k] <- dnorm(dat$yi[k], samples[s, 1L], sd, log = TRUE) +
        log(if (dat$yi[k] >= cut) 1 else .4) - log(p + .4 * (1 - p))
    }
    expect_lte(max(abs(observed - reference)), 1e-10)
  }
})

test_that("F2-1 shortcut retains row-multiplier and source guards", {

  K <- 4L; S <- 3L
  term <- list(block_name = "study", sd_parameter_names = "tau_study",
    model_matrix = matrix(1, K, 1L), n_columns = 1L, structure = "id",
    group_map = seq_len(K), n_groups = K,
    group_covariance = list(type = "known", kernel = diag(c(.5, 1, 2, 4))))
  setup <- function(term) list(data = list(outcome = data.frame(yi = numeric(K), sei = 1),
    location = structure(matrix(1, K, 1L), marginalized_random_effects = list(term))),
    posterior_samples = matrix(c(.1, .2, .3), S, 1L, dimnames = list(NULL, "tau_study")), S = S, K = K)
  plan <- list(random_covariance = list(representation = "diagonal_factor", loading_ranks = rep(0L, K),
    term_names = "study"), row_blocks = as.list(seq_len(K)))
  expect_null(.selection_joint_random_diagonal_samples(setup(term), plan))
  term$row_multiplier <- c(.5, 1, 2, 4); term$row_multiplier_name <- "study_R_diagonal"
  current <- setup(term)
  multiplied <- .selection_joint_random_diagonal_samples(current, plan)
  expect_lte(max(abs(multiplied$diagonal - outer(c(.1, .2, .3)^2, term$row_multiplier))), 1e-15)
  expect_null(.selection_joint_random_diagonal_samples(current, plan, inputs = list()))
  plan$random_covariance$loading_ranks[1L] <- 1L
  expect_null(.selection_joint_random_diagonal_samples(current, plan))
})
