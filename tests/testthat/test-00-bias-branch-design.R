context("Bias-member design classes")

# A RoBMA() ensemble mixes members with different bias terms: no member holds
# both the PET and the PEESE column. The design-based diagnostics use, for each
# posterior draw, the design of the draw's member. The references below are
# written in base R from the declared prior list, the posterior columns, and
# the data, not from the package's design helpers.
#
# The tests that need fitted ensembles are in test-02-bias-branch-design.R,
# which reads the fits prepared by test-01-bias-branch-design.R.

skip_on_cran()

.bbd_data <- function() {

  withr::local_seed(2)
  k <- 20
  d <- data.frame(sei = stats::runif(k, 0.08, 0.30), x = stats::rnorm(k))
  d[["yi"]] <- 0.25 + 0.1 * d[["x"]] + 0.5 * d[["sei"]] +
    stats::rnorm(k, 0, sqrt(0.02^2 + d[["sei"]]^2))
  d[["ni"]] <- 50

  return(d)
}

bbd_data <- .bbd_data()
bbd_args <- list(
  yi = quote(yi), sei = quote(sei), ni = quote(ni), mods = ~ x,
  data = bbd_data, measure = "GEN", standardize_continuous_predictors = FALSE
)

.bbd_only_priors <- function(mods = TRUE) {

  keep <- c("yi", "sei", "ni", if (mods) "mods", "data", "measure",
            "standardize_continuous_predictors")

  return(c(bbd_args[keep], list(only_priors = TRUE, silent = TRUE)))
}

test_that("class designs equal the designs of single-bias models", {

  for (mods in c(TRUE, FALSE)) {
    args     <- .bbd_only_priors(mods)
    ensemble <- suppressWarnings(do.call(RoBMA, c(args, list(model_type = "PP"))))
    single   <- list(
      suppressWarnings(do.call(brma, args)),
      suppressWarnings(do.call(bPET, args)),
      suppressWarnings(do.call(bPEESE, args))
    )

    union <- .get_model_matrix(ensemble)
    expect_identical(
      colnames(union),
      c("(Intercept)", if (mods) "x", "PET", "PEESE")
    )
    for (class in 1:3) {
      expect_equal(
        .bias_class_design(union, TRUE, TRUE, class),
        .get_model_matrix(single[[class]]),
        info = paste("moderators:", mods, "class:", class)
      )
    }
  }

  # an ensemble of the null and PET members holds one bias column
  args    <- .bbd_only_priors()
  pet_null <- suppressWarnings(do.call(RoBMA, c(
    args,
    list(prior_bias = BayesTools::prior_PET("cauchy", list(0, 1), list(0, Inf)))
  )))
  union <- .get_model_matrix(pet_null)
  expect_identical(colnames(union), c("(Intercept)", "x", "PET"))
  expect_identical(.bias_design_member_classes(pet_null[["priors"]]), c(1L, 2L))
  expect_equal(
    .bias_class_design(union, TRUE, FALSE, 2L),
    .get_model_matrix(suppressWarnings(do.call(bPET, args)))
  )
  expect_equal(
    .bias_class_design(union, TRUE, FALSE, 1L),
    .get_model_matrix(suppressWarnings(do.call(brma, args)))
  )
})

test_that("models with one design use the union design", {

  args <- .bbd_only_priors()
  for (constructor in list(brma, bPET, bPEESE)) {
    expect_null(.bias_design_classes(suppressWarnings(do.call(constructor, args))))
  }

  # members of one class form one design
  members <- list(
    BayesTools::prior_PET("cauchy", list(0, 1), list(0, Inf)),
    BayesTools::prior_PET("normal", list(0, 1), list(0, Inf))
  )
  two_pet <- suppressWarnings(do.call(RoBMA, c(
    args, list(prior_bias = members, prior_bias_null = NULL)
  )))
  expect_identical(.bias_design_member_classes(two_pet[["priors"]]), c(2L, 2L))
  expect_null(.bias_design_classes(two_pet))

  # the weight-function members share the location-only design with the null
  selection <- suppressWarnings(do.call(bselmodel, c(
    args[c("yi", "sei", "ni", "data", "measure", "only_priors", "silent")]
  )))
  expect_true(all(.bias_design_member_classes(selection[["priors"]]) == 1L))
  expect_null(.bias_design_classes(selection))

  expect_null(.bias_design_classes(
    structure(list(priors = list(outcome = list(bias = NULL))), class = "brma")
  ))
})

test_that("the ensemble table weights class GVIFs by their draw counts", {

  tables <- list(
    data.frame(term = "a", df = 1L, GVIF = 1, "GVIF^(1/(2*df))" = 1, check.names = FALSE),
    data.frame(term = c("a", "PET"), df = c(1L, 1L), GVIF = c(2, 3),
               "GVIF^(1/(2*df))" = sqrt(c(2, 3)), check.names = FALSE),
    NULL
  )
  designs <- list(matrix(0), matrix(0), matrix(0))

  full <- .vif_ensemble_table(tables, n_draws = c(2L, 6L, 0L), designs = designs)

  expect_identical(full[["term"]], c("a", "PET", "PEESE"))
  expect_equal(full[["GVIF"]], c(0.25 * 1 + 0.75 * 2, 3, NA_real_))
  expect_equal(full[["GVIF^(1/(2*df))"]], sqrt(full[["GVIF"]]))
  expect_identical(names(full), c("term", "df", "GVIF", "GVIF^(1/(2*df))"))
})

# ---------------------------------------------------------------------------- #
# known-V design classes: failures name the draws of the request
# ---------------------------------------------------------------------------- #

# A replay of a known-V ensemble (null, PET and PEESE members, a random
# intercept) with 19 draws whose design classes interleave. The known
# covariance is ill-conditioned only in the draw `ill_row`, where the random
# SD is 0 and V alone holds a pair of estimates with correlation 1 - 1e-12.
# The odd rows are the draws that 'max_samples = 10' keeps, with the classes
# 1, 2, 1, 2, 3, 1, 2, 3, 1, 2.
.bbd_ill_conditioned_replay <- function(ill_row) {

  dat <- data.frame(
    yi  = c(.1, .3, .2, .4, .15, .2),
    x   = c(-1, 0, 1, 2, .3, -.4),
    sei = c(.2, .3, .4, .25, .15, .22),
    study = rep(1:3, each = 2)
  )
  V <- diag(dat[["sei"]]^2)
  V[1L, 2L] <- V[2L, 1L] <- (1 - 1e-12) * dat[["sei"]][1L] * dat[["sei"]][2L]
  object <- RoBMA.mv(
    yi = yi, V = V, mods = ~ x, random = ~ 1 | study, data = dat,
    measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE,
    known_v_parameterization = "block_mvn",
    prior_heterogeneity = BayesTools::prior("normal", list(0, 1), list(0, Inf)),
    prior_heterogeneity_null = NULL, standardize_continuous_predictors = FALSE,
    silent = TRUE,
    prior_bias = list(
      BayesTools::prior_PET("normal", list(0, 1), list(0, Inf)),
      BayesTools::prior_PEESE("normal", list(0, 1), list(0, Inf))
    )
  )
  setup  <- review_compile_object(object)
  priors <- c(
    setup[["compiled"]][["mu"]][["prior_list"]],
    .create_fit_priors(object[["data"]], object[["priors"]])
  )
  gates   <- BayesTools:::.bt_random_variance_allocation_inclusion_indicator_names(
    setup[["object"]][["formula_design"]]
  )
  columns <- unique(c(
    unlist(lapply(names(priors), function(name) {
      BayesTools:::.prior_linear_prior_columns(name, priors[[name]])
    })),
    "bias_indicator", "PET", "PEESE", "tau", "tau_indicator", gates
  ))

  n_draws <- 19L
  samples <- matrix(.1, n_draws, length(columns), dimnames = list(NULL, columns))
  for (name in intersect(names(priors), columns)) {
    location <- priors[[name]][["parameters"]][["location"]]
    if (!is.null(location)) {
      samples[, name] <- location
    }
  }
  samples[, grepl("indicator", columns)] <- 2
  if (length(gates)) {
    samples[, gates] <- 1
  }
  samples[, "bias_indicator"] <- rep(c(3, 3, 2, 1, 1, 2, 3, 1, 2), length.out = n_draws)
  samples[seq(1L, n_draws, by = 2L), "bias_indicator"] <- c(1, 2, 1, 2, 3, 1, 2, 3, 1, 2)
  tau <- rep(.2, n_draws)
  tau[ill_row] <- 0
  samples[, "tau"] <- tau
  samples[, "tau_indicator"] <- ifelse(tau > 0, 2, 1)
  for (column in columns[grepl("xREx.*intercept$|__allocation_sd$", columns)]) {
    samples[, column] <- tau
  }

  replay <- review_attach_replay(setup, samples)
  class(replay) <- setdiff(class(replay), "only_priors.brma")

  return(replay)
}

# the hat computation without the monitored random coefficients that the
# likelihood setup of a replay lacks (they matter at the estimate depth only)
.bbd_replay_hat_setup <- function(object, bias_adjusted = FALSE,
                                  posterior_samples = NULL,
                                  condition_local_effects = TRUE) {

  samples <- .get_posterior_samples(object[["fit"]], posterior_samples)

  return(list(
    posterior_samples = samples,
    data              = object[["data"]],
    yi                = .outcome_data_yi(object),
    S                 = nrow(samples),
    K                 = length(.outcome_data_yi(object))
  ))
}

test_that("known-V design classes report precision failures at the draws of the request", {

  skip_on_cran()
  unavailable <- function(expr) {
    tryCatch(expr, RoBMA_known_v_precision_unavailable = identity)
  }
  # request position of each ill-conditioned row, with and without thinning,
  # and the class of the draw
  cases <- list(
    list(row = 11L, thinned = 6L, class = 1L),
    list(row = 13L, thinned = 7L, class = 2L),
    list(row = 15L, thinned = 8L, class = 3L)
  )
  for (case in cases) {
    replay  <- .bbd_ill_conditioned_replay(case[["row"]])
    samples <- .get_posterior_samples(replay[["fit"]])
    expect_identical(
      .bias_design_classes(replay)[["class"]][case[["row"]]], case[["class"]]
    )

    # the draw of the class that fails, as the class batch itself reports it
    class_rows <- which(.bias_design_classes(replay)[["class"]] == case[["class"]])
    design     <- .bias_design_classes(replay)[["designs"]][[case[["class"]]]]
    local      <- unavailable(.vif_vcov_from_factor_plan(
      replay, design, samples[class_rows, , drop = FALSE],
      .data_known_v_data(replay[["data"]])
    ))
    expect_s3_class(local, "RoBMA_known_v_precision_unavailable")
    expect_identical(local[["rows"]], match(case[["row"]], class_rows))
    expect_lt(local[["reciprocal_condition"]], 1e-9)

    testthat::local_mocked_bindings(
      .estimate_likelihood_setup.brma = .bbd_replay_hat_setup,
      .package = "RoBMA"
    )
    for (request in c("all", "thinned")) {
      expected <- if (request == "all") case[["row"]] else case[["thinned"]]
      max_samples <- if (request == "all") Inf else 10
      calls <- list(
        vif = function() vif(
          replay, posterior_correlation = FALSE, max_samples = max_samples
        ),
        hat = function() .compute_hat_matrix_samples(
          replay, conditioning_depth = "marginal", max_samples = max_samples
        ),
        resid = function() .compute_hat_matrix_samples(
          replay, conditioning_depth = "marginal", max_samples = max_samples,
          return_se = TRUE, return_resid = TRUE
        )
      )
      for (name in names(calls)) {
        condition <- unavailable(suppressWarnings(calls[[name]]()))
        info <- paste(name, "row", case[["row"]], request)

        expect_s3_class(condition, "RoBMA_known_v_precision_unavailable", exact = FALSE)
        expect_identical(
          class(condition),
          c("RoBMA_known_v_precision_unavailable", "RoBMA_numerical_unavailable",
            "error", "condition")
        )
        expect_identical(condition[["rows"]], expected, info = info)
        expect_identical(
          condition[["reciprocal_condition"]], local[["reciprocal_condition"]],
          info = info
        )
        expect_match(
          conditionMessage(condition),
          paste0("at draw ", expected, " among the posterior draws it was applied to"),
          fixed = TRUE, info = info
        )
      }
    }
  }
})
