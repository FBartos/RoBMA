context("Multilevel random-effect components")

test_that("ordinary-normal multilevel ranef uses one coherent joint BLUP", {

  S       <- 4L
  K       <- 5L
  cluster <- c(1L, 1L, 2L, 2L, 2L)
  yi      <- c(0.25, -0.10, 0.30, 0.05, -0.20)
  vi      <- c(0.03, 0.02, 0.04, 0.01, 0.05)

  mu_samples <- matrix(
    c(
      0.05, 0.00, 0.10, 0.05, -0.05,
      0.10, 0.05, 0.15, 0.10,  0.00,
      0.00, 0.00, 0.05, 0.00, -0.10,
      0.15, 0.10, 0.20, 0.15,  0.05
    ),
    nrow = S,
    byrow = TRUE
  )
  tau_within <- matrix(
    rep(c(0.20, 0.15, 0.25, 0.10, 0.18), S),
    nrow = S,
    byrow = TRUE
  )
  tau_between <- matrix(
    c(
      rep(0.12, K),
      rep(0.18, K),
      rep(0.10, K),
      rep(0.15, K)
    ),
    nrow = S,
    byrow = TRUE
  )

  object <- structure(
    list(data = list(outcome = data.frame(
      yi      = yi,
      sei     = sqrt(vi),
      cluster = cluster
    ))),
    class = "brma"
  )
  context <- list(
    object              = object,
    type                = "estimate",
    same_data           = TRUE,
    new_data            = object[["data"]],
    known_V_new         = NULL,
    priors              = list(mods = list()),
    is_mods             = FALSE,
    is_multilevel       = TRUE,
    is_random           = FALSE,
    is_PET              = FALSE,
    is_PEESE            = FALSE,
    is_weightfunction   = FALSE,
    is_weights          = FALSE,
    is_known_v          = FALSE,
    outcome_type        = "norm",
    effect_direction    = "positive",
    posterior_samples   = matrix(0, nrow = S, ncol = 1L),
    outcome_data        = object[["data"]][["outcome"]],
    fit_data            = list(cluster = cluster),
    K_original          = K,
    K                   = K,
    n_chains            = 1L,
    n_iter              = S,
    random_mv           = FALSE,
    bias_adjusted       = FALSE,
    conditional         = FALSE,
    quiet               = TRUE
  )
  scale_state <- list(within = tau_within, between = tau_between)

  original_blup <- .evaluate.brma.multilevel_blup.norm
  blup_calls    <- 0L
  testthat::local_mocked_bindings(
    .is_multilevel = function(...) TRUE,
    .outcome_type = function(...) "norm",
    .is_weightfunction = function(...) FALSE,
    .predict_brma_context = function(...) context,
    .predict_brma_scale_state = function(...) scale_state,
    .evaluate.brma.mu = function(...) mu_samples,
    .evaluate.brma.multilevel_blup.norm = function(...) {
      blup_calls <<- blup_calls + 1L
      original_blup(...)
    },
    .get_estimate_labels = function(...) paste0("estimate-", seq_len(K)),
    .get_cluster_labels = function(...) c(`1` = "cluster-1", `2` = "cluster-2"),
    .package = "RoBMA"
  )

  observed <- ranef.brma(object, simplify = FALSE, expand = TRUE)
  unique   <- ranef.brma(object, simplify = FALSE)
  observed_table  <- as.data.frame(observed)
  observed_tables <- as.data.frame(observed, format = "list")

  expected_cluster  <- matrix(0, nrow = S, ncol = K)
  expected_estimate <- matrix(0, nrow = S, ncol = K)
  for (s in seq_len(S)) {
    for (idx in split(seq_len(K), cluster)) {
      covariance <- diag(vi[idx] + tau_within[s, idx]^2) +
        tcrossprod(tau_between[s, idx])
      weighted_residual <- as.vector(
        solve(covariance, yi[idx] - mu_samples[s, idx])
      )
      expected_estimate[s, idx] <-
        tau_within[s, idx]^2 * weighted_residual
      expected_cluster[s, idx] <- tau_between[s, idx] *
        sum(tau_between[s, idx] * weighted_residual)
    }
  }

  expect_equal(blup_calls, 2L)
  expect_s3_class(observed, "brma_samples_list")
  expect_s3_class(observed_table, "data.frame")
  expect_setequal(
    observed_table[["component"]],
    names(observed)
  )
  expect_named(observed_tables, names(observed))
  expect_true(all(vapply(observed_tables, is.data.frame, logical(1))))
  for (component in names(observed)) {
    expect_identical(
      observed_tables[[component]],
      as.data.frame(observed[[component]])
    )
  }
  expect_identical(
    data.frame(observed),
    data.frame(observed_table)
  )
  expect_false(any(grepl("attr\\(,\\\"class\\\"\\)", capture.output(
    print(observed)
  ))))
  expect_equal(
    unname(as.matrix(observed[["cluster"]])),
    expected_cluster,
    tolerance = 1e-12
  )
  expect_equal(
    unname(as.matrix(observed[["estimate"]])),
    expected_estimate,
    tolerance = 1e-12
  )
  expect_equal(
    unname(as.matrix(observed[["cluster"]]) +
      as.matrix(observed[["estimate"]])),
    expected_cluster + expected_estimate,
    tolerance = 1e-12
  )
  expect_equal(
    unname(as.matrix(unique[["cluster"]])),
    expected_cluster[, c(1L, 3L), drop = FALSE],
    tolerance = 1e-12
  )
  expect_equal(
    unname(as.matrix(unique[["estimate"]])),
    expected_estimate,
    tolerance = 1e-12
  )
  expect_error(
    ranef.brma(object, component = "total"),
    "expand = TRUE",
    fixed = TRUE
  )
})

test_that("selected-normal multilevel ranef retains sampled latent components", {

  S       <- 3L
  K       <- 2L
  fixed   <- matrix(0.2, nrow = S, ncol = K)
  cluster <- matrix(c(-0.1, 0.0, 0.1), nrow = S, ncol = K)
  estimate <- matrix(c(0.05, -0.02, 0.03), nrow = S, ncol = K)
  data <- list(outcome = data.frame(cluster = c(1L, 1L)))
  object <- structure(list(data = data), class = "brma")

  make_samples <- function(samples) {
    .new_brma_samples(
      samples  = samples,
      n_chains = 1L,
      n_iter   = S,
      title    = "Mock:",
      data     = data
    )
  }

  testthat::local_mocked_bindings(
    .is_multilevel = function(...) TRUE,
    .outcome_type = function(...) "norm",
    .is_weightfunction = function(...) TRUE,
    .ranef_brma_multilevel_normal = function(...) {
      stop("The ordinary-normal shortcut must not be used.", call. = FALSE)
    },
    predict.brma = function(object, type, ...) {
      switch(type,
        estimate = make_samples(fixed + cluster + estimate),
        terms    = make_samples(fixed),
        cluster  = make_samples(fixed + cluster)
      )
    },
    .get_estimate_labels = function(...) c("estimate-1", "estimate-2"),
    .get_cluster_labels = function(...) c(`1` = "cluster-1"),
    .package = "RoBMA"
  )

  observed <- ranef.brma(object, simplify = FALSE, expand = TRUE)

  expect_equal(unname(as.matrix(observed[["cluster"]])), cluster)
  expect_equal(unname(as.matrix(observed[["estimate"]])), estimate)
})


test_that("ranef unique-level output rejects row-varying contributions", {

  samples <- matrix(
    c(
      1, 2, 1, 2,
      3, 4, 3, 4
    ),
    nrow  = 2L,
    byrow = TRUE
  )
  observed <- .ranef_unique_level_samples(
    samples      = samples,
    group_map    = c(1L, 2L, 1L, 2L),
    group_levels = c("one", "two"),
    block        = "study"
  )

  expect_equal(unname(observed), samples[, 1:2, drop = FALSE])
  expect_identical(colnames(observed), c("u_study[one]", "u_study[two]"))

  samples[, 3L] <- samples[, 3L] + 1
  expect_error(
    .ranef_unique_level_samples(
      samples      = samples,
      group_map    = c(1L, 2L, 1L, 2L),
      group_levels = c("one", "two"),
      block        = "study"
    ),
    "row-specific contributions",
    fixed = TRUE
  )
})


test_that("ranef unique-level output follows first fitted-row order", {

  samples <- matrix(seq_len(8L), nrow = 2L)
  observed <- .ranef_unique_level_samples(
    samples      = samples,
    group_map    = c(1L, 3L, 2L, 4L),
    group_levels = c("one", "two", "three", "four"),
    block        = "esid_study"
  )

  expect_equal(unname(observed), samples)
  expect_identical(
    colnames(observed),
    paste0("u_esid_study[", c("one", "three", "two", "four"), "]")
  )
})


test_that("ranef keeps any number of random-effect blocks flat", {

  block_names <- c("site", "study", "outcome")
  components  <- lapply(seq_along(block_names), function(i) {
    .new_brma_samples(
      samples  = matrix(i, nrow = 2L, ncol = 2L),
      n_chains = 1L,
      n_iter   = 2L,
      title    = "Random Effects"
    )
  })
  names(components) <- block_names

  all_blocks <- .select_ranef_components(
    components = components,
    component  = "all",
    simplify   = FALSE,
    labels     = c("one", "two"),
    n_chains   = 1L,
    n_iter     = 2L,
    probs      = c(.025, .975),
    data       = NULL
  )
  site <- .select_ranef_components(
    components = components,
    component  = "site",
    labels     = c("one", "two"),
    n_chains   = 1L,
    n_iter     = 2L,
    probs      = c(.025, .975),
    data       = NULL
  )

  expect_named(all_blocks, block_names)
  expect_identical(site, components[["site"]])
  expect_error(
    .select_ranef_components(
      components = components,
      component  = "total",
      labels     = c("one", "two"),
      n_chains   = 1L,
      n_iter     = 2L,
      probs      = c(.025, .975),
      data       = NULL
    ),
    "expand = TRUE",
    fixed = TRUE
  )
})
