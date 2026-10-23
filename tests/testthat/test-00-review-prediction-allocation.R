test_that("singleton prediction pricing uses exact cells and ordered partitions", {

  # The budget that admits chunks of `size` draws beyond the reserve. It inverts
  # the production mapping: .predict_brma_singleton_covariance_chunks() divides
  # the budget left after the reserve by 16 (R/predict.R) and
  # .known_v_covariance_cells_chunk_indices() charges 2 * size + 2 matrices of
  # 8 * K^2 bytes per chunk (R/known-v-gls.R).
  chunk_charge <- function(size, K) 16 * (2 * as.double(size) + 2) * 8 * K

  for (S in c(1L, 5L, 17L)) for (K in c(1L, 8L, 31L)) {
    reserve <- .predict_brma_singleton_covariance_reserve(S, K)
    for (size in unique(c(1L, 2L, S))) {
      budget <- reserve + chunk_charge(size, K)
      chunks <- .predict_brma_singleton_covariance_chunks(S, K, budget)
      expected <- split(seq_len(S), ceiling(seq_len(S) / size))
      expect_identical(chunks, unname(expected))
      expect_identical(unlist(chunks, use.names = FALSE), seq_len(S))
    }
    expect_identical(.predict_brma_singleton_covariance_chunks(S, K, Inf), list(seq_len(S)))
    expect_error(.predict_brma_singleton_covariance_chunks(S, K, reserve), paste0(
      "before Gaussian-parts construction. Increase option ",
      "'RoBMA.known_v_covariance_max_bytes'."
    ), fixed = TRUE)
    expect_error(.predict_brma_singleton_covariance_chunks(S, K, reserve + 1), "single known-V")
  }
  # The reserve (64 K^2 + 16 S K + 2048 K bytes) is pinned by value, not by a
  # copy of its formula, at the grid of the loop and at a large known-V size.
  expect_equal(.predict_brma_singleton_covariance_reserve(1L, 1L), 2128, tolerance = 0)
  expect_equal(.predict_brma_singleton_covariance_reserve(5L, 8L), 21120, tolerance = 0)
  expect_equal(.predict_brma_singleton_covariance_reserve(17L, 31L), 133424, tolerance = 0)
  reserve <- .predict_brma_singleton_covariance_reserve(3L, 6500L)
  expect_equal(reserve, 2717624000, tolerance = 0)
  expect_equal(4 * 1024^3 - reserve, 1577343296, tolerance = 0)
  expect_identical(.predict_brma_singleton_covariance_chunks(3L, 6500L, 4 * 1024^3), list(1:3))
  expect_error(.predict_brma_singleton_covariance_chunks(3L, 1e200, 4 * 1024^3), "before Gaussian-parts")
  for (invalid in list(0, -1, NA_real_, NaN, -Inf, "4e9", c(1, 2), logical(), TRUE)) {
    dense <- tryCatch(.known_v_covariance_chunk_indices(3L, 8L, invalid), error = conditionMessage)
    singleton <- tryCatch(.predict_brma_singleton_covariance_chunks(3L, 8L, invalid), error = conditionMessage)
    expect_identical(singleton, dense)
  }
})


test_that("only declared ordinary singleton prediction receives cheaper pricing", {

  dat <- data.frame(yi = c(.8, -.1, .3), vi = c(.04, .09, .05),
    x = c(-1, 0, 1), paper = c("a", "a", "b"), weights = c(.25, 1.5, 4))
  bias <- function(model = BayesTools::selection_model()) BayesTools::prior_weightfunction(
    "one-sided", .05, BayesTools::wf_fixed(c(1, .4)), model = model)
  make <- function(...) bselmodel(yi = yi, vi = vi, data = dat, measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, ...)
  samples <- cbind(mu = c(.1, .2, -.1), tau = .3)
  context <- function(object, newdata = NULL) .predict_brma_context(object,
    newdata = newdata, V_new = NULL, type = "estimate", conditioning_depth = "marginal",
    conditioning_depth_specified = FALSE, as_measure = TRUE, output_measure = NULL,
    transform = NULL, probs = c(.025, .975), bias_adjusted = FALSE, quiet = TRUE,
    conditional = FALSE, dots = list(.posterior_samples = samples))
  ordinary <- make(prior_bias = bias())
  eligible <- list(ordinary = context(ordinary),
    weighted = context(bselmodel(yi = yi, vi = vi, weights = weights, data = dat,
      measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE, prior_bias = bias())),
    moderators = context(make(mods = ~ x, prior_bias = bias())))
  for (cx in eligible) {
    expect_true(.predict_brma_singleton_covariance_eligible(cx))
    cx$posterior_samples[,] <- NaN
    expect_true(.predict_brma_singleton_covariance_eligible(cx))
  }
  excluded <- list(newdata = context(ordinary, dat),
    scale = context(make(scale = ~ x, prior_bias = bias())),
    multilevel = context(bselmodel(yi = yi, vi = vi, cluster = paper, data = dat,
      measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE, prior_bias = bias())),
    retained = context(make(prior_bias = bias(BayesTools::selection_model(
      estimate_random_effects = "condition")))),
    retained_sampling = context(make(prior_bias = bias(BayesTools::selection_model(
      known_sampling_variance = "condition")))),
    best = context(make(prior_bias = bias(BayesTools::selection_model(weight_rule = "best")))),
    known_v = context(bselmodel.mv(yi = yi, V = diag(dat$vi), data = dat, measure = "GEN",
      prior_unit_information_sd = 1, only_priors = TRUE, prior_bias = bias())),
    compiled = context(bselmodel.mv(yi = yi, vi = vi, random = ~ 1 | paper,
      data = dat, measure = "GEN", prior_unit_information_sd = 1,
      only_priors = TRUE, prior_bias = bias())))
  for (cx in excluded) {
    expect_false(.predict_brma_singleton_covariance_eligible(cx))
    for (budget in c(4 * 1024^3, Inf, 1)) {
      withr::local_options(RoBMA.known_v_covariance_max_bytes = budget)
      expected <- tryCatch(.known_v_covariance_chunk_indices(3L, cx$K, budget / 4), error = conditionMessage)
      actual <- tryCatch(.predict_brma_covariance_chunks(cx), error = conditionMessage)
      expect_identical(actual, expected)
    }
  }
  cx <- eligible$ordinary
  cx$posterior_samples <- samples
  producer <- .predict_brma_from_context
  sentinel <- new.env(parent = environment(producer))
  sentinel$.predict_joint_selection_gaussian_parts <- function(...) stop("parts were reached")
  environment(producer) <- sentinel
  withr::local_options(RoBMA.known_v_covariance_max_bytes =
    .predict_brma_singleton_covariance_reserve(3L, 3L))
  expect_error(producer(cx), "before Gaussian-parts", fixed = TRUE)
})
