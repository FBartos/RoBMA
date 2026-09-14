test_that("outcome-mode funnels reject study-specific fitted quantities", {

  make_object <- function(mods = FALSE, scale = FALSE,
                          classes = "brma") {

    data <- list(outcome = data.frame(yi = c(.1, .2), sei = c(.2, .3)))
    attr(data, "mods")  <- mods
    attr(data, "scale") <- scale
    structure(list(data = data), class = classes)
  }

  expected <- paste0(
    "Outcome-mode funnel plots are not supported for models with location ",
    "or scale predictors.*Use 'residual = TRUE'"
  )

  expect_error(
    funnel(make_object(mods = TRUE), residual = FALSE, as_data = TRUE),
    expected
  )
  expect_error(
    funnel(make_object(scale = TRUE), residual = FALSE, as_data = TRUE),
    expected
  )
})


test_that("intercept-only outcome funnels and automatic residual routing remain", {

  make_object <- function(mods = FALSE, scale = FALSE,
                          classes = "brma") {

    data <- list(outcome = data.frame(yi = c(.1, .2), sei = c(.2, .3)))
    attr(data, "mods")  <- mods
    attr(data, "scale") <- scale
    structure(list(data = data), class = classes)
  }

  testthat::local_mocked_bindings(
    .funnel_common_heterogeneity = function(...) list(common = TRUE),
    .funnel_data_outcome = function(...) list(mode = "outcome"),
    .funnel_data_residual = function(...) list(mode = "residual"),
    .package = "RoBMA"
  )

  intercept_only <- list(
    normal  = make_object(),
    known_V = make_object(classes = c("brma.mv", "brma")),
    GLMM    = make_object(classes = c("brma.glmm", "brma"))
  )
  for (object in intercept_only) {
    expect_identical(
      funnel(object, residual = FALSE, as_data = TRUE),
      list(mode = "outcome")
    )
  }

  expect_identical(
    funnel(make_object(mods = TRUE), type = "outcome", as_data = TRUE),
    list(mode = "residual")
  )
  expect_identical(
    funnel(make_object(scale = TRUE), type = "outcome", as_data = TRUE),
    list(mode = "residual")
  )
})


test_that("funnel location supports canonical multivariate intercepts", {

  data <- data.frame(
    yi    = c(.1, .2),
    study = c("a", "b")
  )
  object <- brma.mv(
    yi          = yi,
    V           = diag(c(.04, .09)),
    random      = ~ 1 | study,
    data        = data,
    measure     = "SMD",
    only_priors = TRUE,
    silent      = TRUE
  )
  posterior_samples <- matrix(
    c(.15, .25),
    ncol     = 1L,
    dimnames = list(NULL, "mu_intercept")
  )

  expect_equal(
    .funnel_mu_samples(object, posterior_samples),
    c(.15, .25)
  )
})

# Small analytic regression for the existing plot-family test files.
# No fitting, constructor sampling, or visual-baseline regeneration.
test_that("fully conditioned selection uses Gaussian funnel and regression quantiles", {

  make_object <- function(conditioned) {

    mode <- if (conditioned) "condition" else "integrate"
    data <- list(outcome = data.frame(yi = c(.1, -.1), sei = c(.2, .4)))
    attr(data, "selection_model") <- structure(list(
      schema_version = 3L, estimate_random_effects = mode,
      other_random_effects = "condition", known_sampling_variance = "condition",
      applicability = list(estimate_random_effects = TRUE,
        other_random_effects = TRUE, known_sampling_variance = TRUE),
      sources = list(random = list(
        list(name = "estimate", role = "estimate", retained = conditioned),
        list(name = "study", role = "other", retained = TRUE)))),
      class = c("RoBMA_selection_model", "list"))
    bias <- BayesTools::prior_weightfunction("one-sided", .025,
      BayesTools::wf_fixed(c(1, .2)), model = selection_model(
        estimate_random_effects = mode, other_random_effects = "condition",
        known_sampling_variance = "condition", group = "study"))
    structure(list(data = data, priors = list(outcome = list(bias = bias)), fit = list()),
      class = c("brma.mv", "brma"))
  }
  selected_calls <- 0L
  testthat::local_mocked_bindings(
    .extract_bias_indicator = function(object, posterior_samples) rep(1L, nrow(posterior_samples)),
    .funnel_mu_samples = function(x, posterior_samples) posterior_samples[, "mu"],
    .selection_context = function(object, posterior_samples) {
      selected_calls <<- selected_calls + 1L
      list(use_normal = rep(FALSE, nrow(posterior_samples)), selected_marker = TRUE)
    },
    .package = "RoBMA"
  )
  object <- make_object(TRUE)
  expect_true(.selection_all_sources_conditioned(object$data))
  se <- c(.2, .4)
  sd <- sqrt(.7^2 + se^2)
  for (direction in c("positive", "negative")) {
    location <- if (direction == "positive") .3 else -.3
    for (S in c(1L, 2L)) {
      posterior <- matrix(rep(location, S), S, 1L, dimnames = list(NULL, "mu"))
      setup <- .funnel_setup_from_samples(object, posterior, rep(.7, S),
        sampling_heterogeneity = TRUE, sampling_bias = TRUE, weights = rep(1, S))
      expect_null(setup$selection)
      expect_false(any(setup$is_weightfunction))
      funnel <- .get_funnel_quantiles_from_setup(se, setup, direction)
      expected_lower <- stats::qnorm(.025, location, sd)
      expected_upper <- stats::qnorm(.975, location, sd)
      expect_equal(funnel$lower, expected_lower, tolerance = 1e-10)
      expect_equal(funnel$upper, expected_upper, tolerance = 1e-10)
      expect_equal(funnel$mid, rep(location, length(se)), tolerance = 1e-10)
      regression <- .regplot_selection_mixture_interval_quantiles(object,
        mean_samples = matrix(location, S, length(se)),
        sd_samples = matrix(rep(sd, each = S), S, length(se)),
        se = se, probs = c(.025, .975), posterior_samples = posterior)
      expect_equal(regression$lower, expected_lower, tolerance = 1e-10)
      expect_equal(regression$upper, expected_upper, tolerance = 1e-10)
    }
  }
  expect_identical(selected_calls, 0L)

})

test_that("scalar selected contour guards preserve supported and Gaussian cases", {

  make <- function(retained = FALSE, role = "estimate", sampling = "integrate",
                   dependent = FALSE, publication = FALSE, legacy = FALSE) {

    data <- list(outcome = data.frame(yi = c(.1, -.1), sei = c(.2, .4)))
    if (!legacy) {
      attr(data, "selection_model") <- structure(list(schema_version = 3L,
        estimate_random_effects = if (retained && role == "estimate") "condition" else "integrate",
        other_random_effects = if (retained && role == "other") "condition" else "integrate",
        known_sampling_variance = sampling,
        applicability = list(estimate_random_effects = role == "estimate",
          other_random_effects = role == "other", known_sampling_variance = TRUE),
        sources = list(random = list(list(name = role, role = role, retained = retained))),
        groups = list(row_blocks = if (publication) list(1:2) else list(1L, 2L))),
        class = c("RoBMA_selection_model", "list"))
      attr(data, "selection_execution_plan") <- structure(list(schema_version = 5L,
        row_blocks = if (dependent) list(1:2) else list(1L, 2L)),
        class = c("RoBMA_selection_execution_plan", "list"))
    }
    list(data = data, priors = list(outcome = list()), fit = list())
  }
  context <- list(use_normal = c(FALSE, FALSE), omega = rbind(c(1, .2), c(1, .3)),
    vector_rule = c(0L, 0L))
  funnel_message <- paste0(
    "Selected funnel contours are unavailable for this joint selection configuration. ",
    "Set 'sampling_bias = FALSE', or use 'zplot()' to view its marginal selected distribution.")
  regression_message <- paste0(
    "Selected regression-plot sampling intervals are unavailable for this joint selection configuration. ",
    "Set 'sampling_bias = FALSE' to draw bias-adjusted sampling intervals.")
  # A conditioned source no longer blocks a funnel contour: its band is the
  # mixture over that source's population law. Regression-plot intervals keep
  # the stricter scalar rule, so they still refuse the same inputs.
  conditioned <- list(
    make(retained = TRUE),
    make(sampling = "condition")
  )
  for (object in conditioned) {
    expect_null(.plot_check_scalar_selection_target(object, context, "funnel"))
    condition <- tryCatch(.plot_check_scalar_selection_target(object, context, "regplot"),
      error = function(e) e)
    expect_identical(conditionMessage(condition), regression_message)
    expect_null(conditionCall(condition))
  }

  # A joint publication event still has no scalar law, for either family.
  rejected <- list(
    make(dependent = TRUE),
    make(publication = TRUE)
  )
  for (i in seq_along(rejected)) {
    selected <- context
    if (i == 2L) selected$vector_rule[] <- 1L
    expect_error(.plot_check_scalar_selection_target(rejected[[i]], selected, "funnel"),
      funnel_message, fixed = TRUE)
    condition <- tryCatch(.plot_check_scalar_selection_target(rejected[[i]], selected, "regplot"),
      error = function(e) e)
    expect_identical(conditionMessage(condition), regression_message)
    expect_null(conditionCall(condition))
  }

  # Full-event product independence and no retained source are scalar-safe.
  expect_null(.plot_check_scalar_selection_target(make(), context, "funnel"))
  expect_null(.plot_check_scalar_selection_target(make(publication = TRUE), context, "regplot"))
  best <- context
  best$vector_rule[] <- 2L
  expect_null(.plot_check_scalar_selection_target(make(), best, "regplot"))
  # Preserve released scalar selection objects without a source specification.
  expect_null(.plot_check_scalar_selection_target(make(legacy = TRUE), best, "funnel"))
  expect_null(.plot_check_scalar_selection_target(make(legacy = TRUE), best, "regplot"))

  # Original disjoint supports cannot certify a changed regression design.
  other <- make(role = "other")
  expect_null(.plot_check_scalar_selection_target(other, context, "funnel"))
  expect_error(.plot_check_scalar_selection_target(other, context, "regplot"),
    regression_message, fixed = TRUE)
  # A fixed-zero scale is structural metadata, not a posterior-draw inference.
  other$priors$outcome$tau <- BayesTools::prior("point", list(location = 0))
  expect_null(.plot_check_scalar_selection_target(other, context, "regplot"))

  impossible_scalar <- make(retained = TRUE, sampling = "condition", dependent = TRUE,
    publication = TRUE)
  expect_null(.plot_check_scalar_selection_target(impossible_scalar, context, "funnel"))
  expect_null(.plot_check_scalar_selection_target(impossible_scalar, context, "regplot"))
  inactive <- context
  inactive$use_normal[] <- TRUE
  flat <- context
  flat$omega[,] <- 2
  for (selected in list(NULL, inactive, flat)) {
    expect_null(.plot_check_scalar_selection_target(make(dependent = TRUE), selected, "funnel"))
    expect_null(.plot_check_scalar_selection_target(make(dependent = TRUE), selected, "regplot"))
  }

  # Exercise the setup owners so direct helper coverage cannot hide a missing
  # guard call. No numerical CDF is invoked for these unavailable configurations.
  calls <- 0L
  testthat::local_mocked_bindings(
    .extract_bias_indicator = function(object, posterior_samples) rep(1L, nrow(posterior_samples)),
    .funnel_mu_samples = function(x, posterior_samples) posterior_samples[, "mu"],
    .selection_context = function(object, posterior_samples) {
      calls <<- calls + 1L
      context["omega"] <- list(context$omega[seq_len(nrow(posterior_samples)), , drop = FALSE])
      context$use_normal <- context$use_normal[seq_len(nrow(posterior_samples))]
      context$vector_rule <- context$vector_rule[seq_len(nrow(posterior_samples))]
      context
    }, .package = "RoBMA"
  )
  object <- make(dependent = TRUE)
  object$priors$outcome$bias <- BayesTools::prior_weightfunction("one-sided", .025,
    BayesTools::wf_fixed(c(1, .2)), model = selection_model(group = "paper"))
  posterior <- matrix(.3, 1L, 1L, dimnames = list(NULL, "mu"))
  expect_error(.funnel_setup_from_samples(object, posterior, .7, TRUE, TRUE, 1),
    funnel_message, fixed = TRUE)
  expect_error(.regplot_selection_setup(object, posterior), regression_message, fixed = TRUE)
  expect_identical(calls, 2L)
  expect_null(.funnel_setup_from_samples(object, posterior, .7, TRUE, FALSE, 1)$selection)
  expect_identical(calls, 2L)
})


test_that("conditioned sources give the funnel band their mixture, not their spread", {

  skip_if_not(.has_native_selnorm_kernel())

  spec <- .test_step_spec(c(.2, .1), c(.2, .2))
  selection <- spec
  selection[["omega"]]       <- matrix(spec[["fixed_omega"]], nrow = 1L)
  selection[["alpha"]]       <- 0
  selection[["phack_kind"]]  <- 0L
  selection[["kernel_mode"]] <- SELKERNEL_STEP
  selection[["use_normal"]]  <- FALSE
  selection[["vector_rule"]] <- 0L

  mu <- 0.30
  base <- list(
    mu = mu, tau = 0, PET = 0, PEESE = 0,
    is_weightfunction = TRUE, selection = selection, weights = 1
  )
  band <- function(integrated, conditioned, retains_sampling, se) {
    setup <- base
    setup[["sources"]] <- list(
      integrated = integrated, conditioned = conditioned,
      retains_sampling = retains_sampling, common = TRUE
    )
    .funnel_model_averaged_quantiles_native(se, setup, "positive")
  }

  # Reference: mix selected normals over the conditioned offset, each carrying
  # its own normalizer, which is what a retained source implies.
  static <- BayesTools::selection_native_static_args(selection)
  omega  <- as.numeric(selection[["omega"]])
  reference <- function(v_integrated, v_conditioned, se, probs = c(.025, .975)) {
    edges <- function(bin) {
      c(static[["z_lower"]][bin] * se, static[["z_upper"]][bin] * se)
    }
    weight_of <- function(y) {
      z <- y / se
      omega[vapply(z, function(value) {

        which(value >= static[["z_lower"]] & value < static[["z_upper"]])[1L]
      }, integer(1))]
    }
    rule  <- .gauss_hermite_nodes(301L)
    nodes <- rule[["nodes"]]
    node_weights <- rule[["weights"]] / sum(rule[["weights"]])
    sd_kernel <- sqrt(v_integrated)
    spread <- sqrt(v_integrated + v_conditioned)
    grid <- seq(mu - 14 * spread, mu + 14 * spread, length.out = 60001L)
    mass <- weight_of(grid)
    density <- numeric(length(grid))
    for (node in seq_along(nodes)) {
      centre <- mu + nodes[node] * sqrt(v_conditioned)
      normalizer <- sum(vapply(seq_along(omega), function(bin) {

        edge <- edges(bin)
        omega[bin] * (stats::pnorm(edge[2L], centre, sd_kernel) -
                        stats::pnorm(edge[1L], centre, sd_kernel))
      }, numeric(1)))
      density <- density + node_weights[node] * mass *
        stats::dnorm(grid, centre, sd_kernel) / normalizer
    }
    cumulative <- c(0, cumsum(
      (density[-1L] + density[-length(density)]) / 2 * diff(grid)
    ))
    cumulative <- cumulative / cumulative[length(cumulative)]
    # The far tails are flat to working precision; drop the ties so the
    # inversion is well defined rather than silently interpolated across them.
    keep <- !duplicated(cumulative)
    stats::approx(cumulative[keep], grid[keep], xout = probs)$y
  }

  se <- 0.2
  # Every source integrated: the band is the ordinary selected normal, and the
  # previous scalar construction from a total tau reproduces it exactly.
  integrated_only <- band(0.04, 0, FALSE, se)
  legacy <- base
  legacy[["tau"]] <- sqrt(0.04)
  scalar <- .funnel_model_averaged_quantiles_native(se, legacy, "positive")
  expect_equal(integrated_only[["lower"]], scalar[["lower"]], tolerance = 1e-12)
  expect_equal(integrated_only[["upper"]], scalar[["upper"]], tolerance = 1e-12)

  # A retained random source: the band is the mixture, not the same spread.
  conditioned <- band(0, 0.04, FALSE, se)
  target <- reference(se^2, 0.04, se)
  expect_equal(conditioned[["lower"]], target[1L], tolerance = 1e-4)
  expect_equal(conditioned[["upper"]], target[2L], tolerance = 1e-4)
  expect_gt(abs(conditioned[["lower"]] - integrated_only[["lower"]]), 0.01)

  # Retained sampling moves the standard error into the mixed bucket instead.
  retained_sampling <- band(0.04, 0, TRUE, se)
  expect_equal(retained_sampling[["lower"]],
               reference(0.04, se^2, se)[1L], tolerance = 1e-4)
  expect_gt(abs(retained_sampling[["lower"]] - integrated_only[["lower"]]), 0.01)

  # Both roles at once still sum to the same total spread as the legacy band.
  mixed <- band(0.02, 0.02, FALSE, se)
  expect_true(is.finite(mixed[["lower"]]) && is.finite(mixed[["upper"]]))
  expect_lt(mixed[["lower"]], mixed[["upper"]])
})
