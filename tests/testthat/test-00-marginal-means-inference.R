context("Marginal-means inference")


.fixed_marginal_posterior <- function(value = 0, n = 20L) {

  posterior <- structure(
    rep(value, n),
    class = c("marginal_posterior", "numeric")
  )
  BayesTools::posterior_metadata(posterior, "atoms") <-
    BayesTools::posterior_atom_attribute(data.frame(x = value, mass = 1))
  return(posterior)
}


test_that("BF-free marginal means retain structurally fixed cells", {

  posterior <- list(zero = .fixed_marginal_posterior())
  inference <- .marginal_means_inclusion_bf(
    posterior       = posterior,
    null_hypothesis = 0,
    compute         = FALSE
  )

  expect_named(inference, "zero")
  expect_true(is.na(inference[["zero"]]))
})


test_that("structurally fixed marginal cells have unavailable BFs", {

  posterior <- list(zero = .fixed_marginal_posterior())
  refusal <- NULL
  testthat::local_mocked_bindings(
    Savage_Dickey_BF = function(...) stop(refusal),
    .package = "BayesTools"
  )

  # BayesTools classes the refusal of a posterior point mass at the null; the
  # class, not the message, marks the cell as structurally fixed.
  refusal <- structure(
    class = c(
      "BayesTools_posterior_point_mass_at_null",
      "BayesTools_hypothesis_ordinate", "error", "condition"
    ),
    list(message = "Reworded refusal.", call = NULL)
  )
  inference <- .marginal_means_inclusion_bf(
    posterior       = posterior,
    null_hypothesis = 0,
    compute         = TRUE
  )

  expect_named(inference, "zero")
  expect_true(is.na(inference[["zero"]]))
  expect_match(
    attr(inference[["zero"]], "warnings", exact = TRUE),
    "structurally fixed"
  )

  refusal <- simpleError(paste0(
    "The posterior contains a declared point mass at the exact null ",
    "hypothesis value. The ordinary Savage-Dickey density ratio is invalid."
  ))
  expect_error(
    .marginal_means_inclusion_bf(
      posterior       = posterior,
      null_hypothesis = 0,
      compute         = TRUE
    ),
    "declared point mass at the exact null",
    fixed = TRUE
  )
})


test_that("marginal means whose Savage-Dickey Bayes factor is refused have NA rows", {

  # One row per refusal class of BayesTools' ordinate check, next to a row
  # with a Bayes factor: each refused row has an NA Bayes factor with the
  # reason of its class, and the other row is computed.
  refusals <- c(
    zero      = "BayesTools_zero_ordinate",
    infinite  = "BayesTools_infinite_ordinate",
    undefined = "BayesTools_undefined_ordinate",
    inexact   = "BayesTools_inexact_ordinate",
    mass      = "BayesTools_point_mass_at_null",
    fixed     = "BayesTools_posterior_point_mass_at_null"
  )
  posterior <- lapply(seq_len(length(refusals) + 1L), function(i) {
    structure(rep(i, 20L), class = c("marginal_posterior", "numeric"))
  })
  names(posterior) <- c(names(refusals), "computed")
  testthat::local_mocked_bindings(
    Savage_Dickey_BF = function(posterior, ...) {
      i <- posterior[[1L]]
      if (i > length(refusals)) {
        return(structure(2.5, warnings = "A computed row."))
      }
      stop(structure(
        class = c(refusals[[i]], "BayesTools_hypothesis_ordinate", "error", "condition"),
        list(message = "Reworded refusal.", call = NULL)
      ))
    },
    .package = "BayesTools"
  )
  inference <- .marginal_means_inclusion_bf(
    posterior       = posterior,
    null_hypothesis = 0,
    compute         = TRUE
  )

  expect_named(inference, names(posterior))
  expect_identical(inference[["computed"]], structure(2.5, warnings = "A computed row."))
  reasons <- c(
    zero      = paste0("The prior density of the marginal mean at the null ",
                       "hypothesis is zero; its inclusion Bayes factor is undefined."),
    infinite  = paste0("The prior density of the marginal mean at the null ",
                       "hypothesis is infinite; its inclusion Bayes factor is undefined."),
    undefined = paste0("The prior density of the marginal mean at the null ",
                       "hypothesis is undefined; its inclusion Bayes factor is undefined."),
    inexact   = paste0("The prior density of the marginal mean at the null ",
                       "hypothesis has no exact value; its inclusion Bayes factor ",
                       "is unavailable. Test a region hypothesis with 'hypothesis()' ",
                       "instead."),
    mass      = paste0("The prior of the marginal mean has a point mass at the ",
                       "null hypothesis; its inclusion Bayes factor is undefined."),
    fixed     = paste0("The marginal mean is structurally fixed at the null ",
                       "hypothesis; its inclusion Bayes factor is undefined.")
  )
  for (row in names(refusals)) {
    expect_true(is.na(inference[[row]]), info = row)
    expect_identical(attr(inference[[row]], "warnings", exact = TRUE), reasons[[row]],
                     info = row)
  }

  # qCMDE/IWMDE: each level with a precomputed ordinate is evaluated on its
  # own, so a refused level does not stop the others.
  ordinate <- BayesTools::posterior_ordinate_attribute(
    0, 1, "q_grid_cmde", "qCMDE",
    diagnostics = list(estimator = "q_grid_cmde", ordinate_relative_change = 0)
  )
  precomputed <- lapply(posterior[c("inexact", "computed")], function(level) {
    BayesTools::posterior_metadata(level, "posterior_ordinate") <- ordinate
    level
  })
  testthat::local_mocked_bindings(
    .iwmde_posterior_ordinate_matches_request = function(...) TRUE,
    .package = "RoBMA"
  )
  bf <- .marginal_means_iwmde_bf(precomputed, 0, density_method = "qCMDE")
  expect_identical(as.numeric(bf[["computed"]]), 2.5)
  expect_true(is.na(bf[["inexact"]]))
  expect_identical(attr(bf[["inexact"]], "warnings", exact = TRUE), reasons[["inexact"]])
  scalar <- .marginal_means_iwmde_bf(precomputed[["inexact"]], 0,
                                     density_method = "qCMDE")
  expect_true(is.na(scalar))
  expect_identical(attr(scalar, "warnings", exact = TRUE), reasons[["inexact"]])

  # Any other failure of a row's Bayes factor stops the table unchanged, on
  # the KDE and on the qCMDE/IWMDE paths: only the ordinate refusals above
  # give NA rows.
  failures <- list(
    simpleError("Another failure."),
    structure(
      class = c("BayesTools_refit_required", "error", "condition"),
      list(message = "Refit the model.", call = NULL)
    )
  )
  for (failure in failures) {
    local({
      testthat::local_mocked_bindings(
        Savage_Dickey_BF = function(...) stop(failure),
        .package = "BayesTools"
      )
      info <- class(failure)[[1L]]
      expect_identical(
        tryCatch(.marginal_means_inclusion_bf(posterior["computed"], 0, TRUE),
                 error = identity),
        failure,
        info = info
      )
      expect_identical(
        tryCatch(.marginal_means_iwmde_bf(precomputed, 0, density_method = "qCMDE"),
                 error = identity),
        failure,
        info = info
      )
      expect_identical(
        tryCatch(.marginal_means_iwmde_bf(precomputed[["computed"]], 0,
                                          density_method = "qCMDE"),
                 error = identity),
        failure,
        info = info
      )
    })
  }
})


test_that("marginal_means() reports rows with a zero prior ordinate at the null as NA", {

  skip_on_cran()
  set.seed(3)
  k    <- 30L
  data <- data.frame(x = stats::rnorm(k), sei = stats::runif(k, 0.1, 0.3))
  data[["yi"]] <- stats::rnorm(k, 0.3, data[["sei"]])
  # The effect prior is truncated at 0: the intercept (and the mean at the
  # mean of 'x') has a zero prior density at -0.5, the means at -1 and +1 SD
  # of 'x' do not.
  fit <- suppressWarnings(brma(
    yi = yi, sei = sei, mods = ~ x, data = data, measure = "SMD",
    prior_effect = prior("normal", list(0, 1), truncation = list(0, Inf)),
    chains = 1, sample = 1000, burnin = 200, adapt = 100, seed = 1,
    silent = TRUE
  ))
  means <- suppressWarnings(marginal_means(fit, null_hypothesis = -0.5, bf = TRUE))
  inference <- means[["inference"]]
  reason <- paste0("The prior density of the marginal mean at the null ",
                   "hypothesis is zero; its inclusion Bayes factor is undefined.")
  refused <- list(c("mu_intercept", "intercept"), c("mu_x", "0SD"))
  computed <- list(c("mu_x", "-1SD"), c("mu_x", "1SD"))
  bf_of <- function(cell) {
    value <- inference[["inference"]][[cell[[1L]]]]
    if (is.list(value)) value[[cell[[2L]]]] else value
  }
  for (cell in refused) {
    info <- paste(cell, collapse = " ")
    expect_true(is.na(bf_of(cell)), info = info)
    expect_identical(attr(bf_of(cell), "warnings", exact = TRUE), reason, info = info)
  }
  # The other rows are BayesTools' Savage-Dickey Bayes factors of their own
  # conditional marginal posteriors, as before.
  for (cell in computed) {
    posterior <- inference[["conditional"]][[cell[[1L]]]][[cell[[2L]]]]
    class(posterior) <- unique(c(class(posterior), "marginal_posterior"))
    expect_identical(
      bf_of(cell),
      BayesTools::Savage_Dickey_BF(posterior, null_hypothesis = -0.5,
                                   silent = TRUE, density_method = "KDE"),
      info = paste(cell, collapse = " ")
    )
  }

  # The summary prints the reason after the row label; the data frame has NA.
  table <- summary(means)
  expect_true(all(c(paste0("intercept: ", reason), paste0("x[0SD]: ", reason)) %in%
                    attr(table, "warnings")))
  frame <- as.data.frame(means)
  expect_identical(frame[["parameter"]], c("intercept", "x[-1SD]", "x[0SD]", "x[1SD]"))
  expect_identical(is.na(frame[["inclusion_BF"]]), c(TRUE, FALSE, TRUE, FALSE))
  expect_identical(as.numeric(frame[["inclusion_BF"]][c(2L, 4L)]),
                   vapply(computed, function(cell) as.numeric(bf_of(cell)), numeric(1)))

  # A point hypothesis on a refused mean stops with BayesTools' class.
  expect_error(hypothesis(means, "intercept = -0.5"), class = "BayesTools_zero_ordinate")
  expect_error(hypothesis(means, "x[0SD] = -0.5"), class = "BayesTools_zero_ordinate")
})


test_that("interaction marginals condition on every contributing coefficient", {

  terms <- c("intercept", "a", "b", "a:b")
  parameters <- c("mu_intercept", "mu_a", "mu_b", "mu_ab")
  spike_and_slab_priors <- stats::setNames(lapply(parameters, function(parameter) {

    BayesTools::prior_spike_and_slab(
      prior_parameter = BayesTools::prior("normal", list(mean = 0, sd = 1))
    )
  }), parameters)
  object <- list(fit = structure(list(), prior_list = spike_and_slab_priors))
  conditional_list <- .marginal_means_conditional_list(
    object     = object,
    terms      = terms,
    parameters = parameters
  )

  cell_weights <- rbind(
    base  = c(1, 0, 0, 0),
    a     = c(1, 1, 0, 0),
    b     = c(1, 0, 1, 0),
    `a:b` = c(1, 1, 1, 1)
  )
  colnames(cell_weights) <- parameters
  marginal <- lapply(seq_len(nrow(cell_weights)), function(i) {

    with_draw_metadata(
      numeric(20),
      linear_weights = cell_weights[i, ]
    )
  })
  names(marginal) <- rownames(cell_weights)
  prior_list <- stats::setNames(lapply(parameters, function(parameter) {

    BayesTools::prior("normal", parameters = list(mean = 0, sd = 1))
  }), parameters)

  effective <- BayesTools:::.marginal_inference_level_conditionals(
    marginal    = marginal,
    prior_list  = prior_list,
    conditional = conditional_list[["mu_ab"]]
  )
  expected <- lapply(seq_len(nrow(cell_weights)), function(i) {

    names(cell_weights[i, ])[cell_weights[i, ] != 0]
  })
  names(expected) <- rownames(cell_weights)
  expect_identical(effective, expected)

  # Independent finite-mixture reference for a point-null BF in the interaction
  # cell. Each active coefficient contributes an independent Normal variate.
  states <- as.matrix(expand.grid(
    mu_intercept = 0:1,
    mu_a         = 0:1,
    mu_b         = 0:1,
    mu_ab        = 0:1
  ))
  coefficient_mean <- c(
    mu_intercept = 0.4,
    mu_a         = -1.5,
    mu_b         = 1.8,
    mu_ab        = 0.8
  )
  coefficient_sd <- c(
    mu_intercept = 0.6,
    mu_a         = 0.2,
    mu_b         = 0.2,
    mu_ab        = 0.5
  )
  prior_probability <- rep(1 / nrow(states), nrow(states))
  posterior_log_weight <- as.numeric(
    states %*% c(mu_intercept = -1, mu_a = 2.5, mu_b = 2, mu_ab = -2)
  ) + 2 * states[, "mu_a"] * states[, "mu_b"]
  posterior_probability <- exp(
    posterior_log_weight - max(posterior_log_weight)
  )
  posterior_probability <- posterior_probability / sum(posterior_probability)

  conditional_density <- function(probability, conditional) {

    event <- rowSums(states[, conditional, drop = FALSE]) > 0
    state_mean <- as.numeric(states %*% coefficient_mean)
    state_sd <- sqrt(as.numeric((states^2) %*% (coefficient_sd^2)))
    event_probability <- probability[event] / sum(probability[event])
    sum(event_probability * stats::dnorm(
      0,
      mean = state_mean[event],
      sd   = state_sd[event]
    ))
  }
  mixture_bf <- function(conditional) {

    conditional_density(prior_probability, conditional) /
      conditional_density(posterior_probability, conditional)
  }

  contributing <- names(cell_weights["a:b", ])[cell_weights["a:b", ] != 0]
  actual_bf     <- mixture_bf(effective[["a:b"]])
  reference_bf  <- mixture_bf(contributing)
  incomplete_bf <- mixture_bf(c("mu_intercept", "mu_ab"))

  expect_equal(actual_bf, reference_bf, tolerance = 1e-14)
  expect_equal(actual_bf, 0.3565699, tolerance = 1e-7)
  expect_gt(abs(log(actual_bf / incomplete_bf)), log(2))
})


test_that("only coefficients with inclusion indicators are conditional candidates", {

  parameters <- c("mu_intercept", "mu_a", "mu_b")
  terms      <- c("intercept", "a", "b")
  make_object <- function(prior_list) {

    list(fit = structure(list(), prior_list = prior_list))
  }
  plain          <- BayesTools::prior("normal", list(mean = 0, sd = 1))
  spike_and_slab <- BayesTools::prior_spike_and_slab(prior_parameter = plain)
  null_mixture   <- BayesTools::prior_mixture(
    list(
      BayesTools::prior("spike",  list(location = 0)),
      BayesTools::prior("normal", list(mean = 0, sd = 1))
    ),
    components = c("null", "alternative")
  )

  # Plain priors are always included in every model, so BayesTools rejects them
  # as conditional labels; they must not be offered as candidates.
  plain_list <- .marginal_means_conditional_list(
    object     = make_object(stats::setNames(
      list(plain, plain, plain), parameters
    )),
    terms      = terms,
    parameters = parameters
  )
  expect_named(plain_list, parameters)
  expect_identical(plain_list[["mu_a"]], character(0))

  mixed_list <- .marginal_means_conditional_list(
    object     = make_object(stats::setNames(
      list(plain, spike_and_slab, null_mixture), parameters
    )),
    terms      = terms,
    parameters = parameters
  )
  expect_identical(mixed_list[["mu_intercept"]], c("mu_a", "mu_b"))
  expect_identical(mixed_list[["mu_b"]], c("mu_a", "mu_b"))

  # A parameter missing from the prior list is not a conditional candidate.
  partial_list <- .marginal_means_conditional_list(
    object     = make_object(list(mu_a = spike_and_slab)),
    terms      = terms,
    parameters = parameters
  )
  expect_identical(partial_list[["mu_a"]], "mu_a")
})
