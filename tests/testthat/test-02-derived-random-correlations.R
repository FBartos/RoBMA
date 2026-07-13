context("Derived scalar random-effect correlations")

source(testthat::test_path("common-functions.R"))

.derived_random_term <- function(structure = "ar1", n_columns = 3L,
                                 monitor_correlation = TRUE,
                                 sample_fixed = NULL,
                                 rho_scale = "rho",
                                 bounds = c(lower = -1, upper = 1),
                                 distance_matrix = NULL,
                                 time_values = NULL,
                                 block_name = "study",
                                 parameter_stem = paste0(
                                   "mu__xREx__", block_name
                                 )) {

  list(
    structure = structure,
    block_name = block_name,
    parameter_stem = parameter_stem,
    sd_parameter_names = paste0(parameter_stem, "_sd"),
    n_columns = n_columns,
    monitor = list(correlation = monitor_correlation),
    correlation = list(
      type = "rho",
      structure = structure,
      rho_name = paste0(parameter_stem, "_rho"),
      sample_name = paste0(parameter_stem, "_rho"),
      sample_fixed = sample_fixed,
      rho_scale = rho_scale,
      bounds = bounds
    ),
    car = if (is.null(distance_matrix) && is.null(time_values)) {
      NULL
    } else {
      list(
        distance_matrix = distance_matrix,
        time_values = time_values
      )
    }
  )
}

.derived_correlation_mcmc <- function(rho = c(-0.4, 0, 0.6)) {

  values <- cbind(
    mu = seq_along(rho),
    mu__xREx__study_rho = rho
  )
  coda::mcmc.list(coda::mcmc(values))
}

.derived_correlation_matrix <- function(draws, draw, n_columns) {

  base <- "mu__xREx__study_xRE_CORx_R"
  names <- unlist(lapply(seq_len(n_columns), function(column) {
    paste0(base, "[", seq_len(n_columns), ",", column, "]")
  }), use.names = FALSE)
  matrix(as.matrix(draws[[1L]])[draw, names], nrow = n_columns)
}

test_that("compact rho draws reconstruct scalar correlation matrices", {

  rho <- c(-0.4, 0, 0.6)
  for (structure in c("cs", "hcs", "ar1", "har")) {
    term <- .derived_random_term(structure = structure)
    draws <- RoBMA:::.brma_append_derived_random_correlation(
      mcmc_list = .derived_correlation_mcmc(rho),
      random_term = term
    )
    distance <- if (structure %in% c("cs", "hcs")) {
      1 - diag(3L)
    } else {
      abs(outer(seq_len(3L), seq_len(3L), "-"))
    }

    for (draw in seq_along(rho)) {
      expect_equal(
        .derived_correlation_matrix(draws, draw, 3L),
        rho[[draw]]^distance,
        tolerance = 0,
        info = structure
      )
    }
  }

  car_time     <- c(0, 0.5, 2)
  car_distance <- abs(outer(car_time, car_time, "-"))
  term <- .derived_random_term(
    structure = "car",
    time_values = car_time
  )
  draws <- RoBMA:::.brma_append_derived_random_correlation(
    mcmc_list = .derived_correlation_mcmc(rho = c(0.25, 0.81)),
    random_term = term
  )
  expect_equal(
    .derived_correlation_matrix(draws, 2L, 3L),
    0.81^car_distance,
    tolerance = 1e-15
  )

  legacy_term <- .derived_random_term(
    structure = "car",
    distance_matrix = car_distance
  )
  legacy_draws <- RoBMA:::.brma_append_derived_random_correlation(
    mcmc_list = .derived_correlation_mcmc(rho = 0.81),
    random_term = legacy_term
  )
  expect_equal(
    .derived_correlation_matrix(legacy_draws, 1L, 3L),
    0.81^car_distance,
    tolerance = 1e-15
  )
})

test_that("fixed rho coordinates and one-column structures are reconstructed", {

  raw <- coda::mcmc.list(coda::mcmc(cbind(mu = 1:2)))
  fisher <- .derived_random_term(
    structure = "cs",
    sample_fixed = 0.5,
    rho_scale = "fisher_z"
  )
  fisher_draws <- RoBMA:::.brma_append_derived_random_correlation(raw, fisher)
  expected <- matrix(tanh(0.5), nrow = 3L, ncol = 3L)
  diag(expected) <- 1
  expect_equal(
    .derived_correlation_matrix(fisher_draws, 1L, 3L),
    expected,
    tolerance = 1e-15
  )

  logit <- .derived_random_term(
    structure = "ar1",
    sample_fixed = 0,
    rho_scale = "logit",
    bounds = c(lower = 0, upper = 1)
  )
  logit_draws <- RoBMA:::.brma_append_derived_random_correlation(raw, logit)
  expect_equal(
    .derived_correlation_matrix(logit_draws, 1L, 3L),
    0.5^abs(outer(seq_len(3L), seq_len(3L), "-")),
    tolerance = 1e-15
  )

  scalar <- .derived_random_term(structure = "cs", n_columns = 1L)
  scalar_draws <- RoBMA:::.brma_append_derived_random_correlation(raw, scalar)
  expect_equal(.derived_correlation_matrix(scalar_draws, 1L, 1L), matrix(1))
})

test_that("monitoring and quadratic expansion guards preserve the schema", {

  raw <- .derived_correlation_mcmc()
  unmonitored <- .derived_random_term(monitor_correlation = FALSE)
  expect_identical(
    RoBMA:::.brma_append_derived_random_correlation(raw, unmonitored),
    raw
  )

  partial <- as.matrix(raw[[1L]])
  partial <- cbind(partial, mu__xREx__study_xRE_CORx_R.1.1. = 1)
  colnames(partial)[ncol(partial)] <-
    "mu__xREx__study_xRE_CORx_R[1,1]"
  partial <- coda::mcmc.list(coda::mcmc(partial))
  expect_error(
    RoBMA:::.brma_append_derived_random_correlation(
      partial,
      .derived_random_term()
    ),
    "incomplete random-effect correlation matrix"
  )

  existing <- RoBMA:::.brma_append_derived_random_correlation(
    raw,
    .derived_random_term()
  )
  expect_identical(
    RoBMA:::.brma_append_derived_random_correlation(
      existing,
      .derived_random_term()
    ),
    existing
  )

  old <- getOption("RoBMA.max_derived_correlation_cells")
  on.exit(options(RoBMA.max_derived_correlation_cells = old), add = TRUE)
  options(RoBMA.max_derived_correlation_cells = 4L)
  expect_error(
    RoBMA:::.brma_check_derived_random_correlation_budget(
      random_terms = list(.derived_random_term(n_columns = 3L)),
      mcmc_list    = raw
    ),
    "would derive 27 dense random-correlation cells"
  )

  options(RoBMA.max_derived_correlation_cells = 1e6)
  expect_error(
    RoBMA:::.brma_check_derived_random_correlation_budget(
      random_terms = list(.derived_random_term(n_columns = 4096L)),
      mcmc_list    = raw
    ),
    "would derive 50331648 dense random-correlation cells"
  )

  options(RoBMA.max_derived_correlation_cells = 100L)
  malformed <- .derived_random_term()
  malformed[["n_columns"]] <- 2.5
  expect_error(
    RoBMA:::.brma_append_derived_random_correlation(raw, malformed),
    "Invalid scalar random-correlation metadata"
  )
})

test_that("multiple compact blocks preserve block-local order across chains", {

  study <- .derived_random_term(
    structure = "ar1",
    n_columns = 3L,
    block_name = "study"
  )
  lab <- .derived_random_term(
    structure = "cs",
    n_columns = 2L,
    block_name = "lab"
  )
  make_chain <- function(study_rho, lab_rho) {

    values <- cbind(
      mu = seq_along(study_rho),
      study_sd = 0.3,
      study_rho = study_rho,
      lab_sd = 0.4,
      lab_rho = lab_rho,
      sigma = 1
    )
    colnames(values)[2:5] <- c(
      study[["sd_parameter_names"]],
      study[["correlation"]][["rho_name"]],
      lab[["sd_parameter_names"]],
      lab[["correlation"]][["rho_name"]]
    )
    return(coda::mcmc(values, start = 10, thin = 2))
  }
  raw <- coda::mcmc.list(
    make_chain(c(0.2, 0.4, 0.6), c(-0.1, 0.1, 0.3)),
    make_chain(c(0.1, 0.3, 0.5), c(-0.2, 0.2, 0.4))
  )
  draws <- RoBMA:::.brma_append_derived_random_correlation_terms(
    mcmc_list    = raw,
    random_terms = list(study, lab)
  )
  study_names <- RoBMA:::.brma_derived_random_correlation_names(
    RoBMA:::.brma_derived_random_correlation_spec(study)
  )
  lab_names <- RoBMA:::.brma_derived_random_correlation_names(
    RoBMA:::.brma_derived_random_correlation_spec(lab)
  )
  expect_identical(
    colnames(as.matrix(draws[[1L]])),
    c(
      "mu",
      study[["sd_parameter_names"]],
      study[["correlation"]][["rho_name"]],
      study_names,
      lab[["sd_parameter_names"]],
      lab[["correlation"]][["rho_name"]],
      lab_names,
      "sigma"
    )
  )
  expect_identical(
    colnames(as.matrix(draws[[2L]])),
    colnames(as.matrix(draws[[1L]]))
  )
  expect_equal(coda::mcpar(draws[[1L]]), coda::mcpar(raw[[1L]]))
  expect_equal(coda::mcpar(draws[[2L]]), coda::mcpar(raw[[2L]]))
  expect_equal(
    matrix(as.matrix(draws[[2L]])[3L, study_names], nrow = 3L),
    0.5^abs(outer(seq_len(3L), seq_len(3L), "-")),
    tolerance = 1e-15
  )
  expected_lab <- matrix(0.4, nrow = 2L, ncol = 2L)
  diag(expected_lab) <- 1
  expect_equal(
    matrix(as.matrix(draws[[2L]])[3L, lab_names], nrow = 2L),
    expected_lab,
    tolerance = 1e-15
  )

  old <- getOption("RoBMA.max_derived_correlation_cells")
  on.exit(options(RoBMA.max_derived_correlation_cells = old), add = TRUE)
  options(RoBMA.max_derived_correlation_cells = 77L)
  expect_error(
    RoBMA:::.brma_check_derived_random_correlation_budget(
      random_terms = list(study, lab),
      mcmc_list    = raw
    ),
    "would derive 78 dense random-correlation cells"
  )
})

test_that("default draws derive public matrices while raw draws stay compact", {

  fit <- .derived_correlation_mcmc()
  attr(fit, "formula_design") <- list(
    mu = list(random_effects = list(.derived_random_term()))
  )
  object <- structure(list(fit = fit), class = "brma")

  default <- RoBMA:::.brma_to_mcmc.list(object)
  raw     <- RoBMA:::.brma_to_mcmc.list(object, include_auxiliary = TRUE)
  expect_true(any(grepl("_xRE_CORx_R[", colnames(as.matrix(default[[1L]])),
                        fixed = TRUE)))
  expect_false(any(grepl("_xRE_CORx_R[", colnames(as.matrix(raw[[1L]])),
                         fixed = TRUE)))
})

test_that("all public draw converters retain scalar correlation semantics", {

  fit <- .derived_correlation_mcmc()
  attr(fit, "formula_design") <- list(
    mu = list(random_effects = list(.derived_random_term()))
  )
  object <- structure(list(fit = fit), class = "brma")
  correlation_names <- unlist(lapply(seq_len(3L), function(column) {
    paste0(
      "mu__xREx__study_xRE_CORx_R[",
      seq_len(3L), ",", column, "]"
    )
  }), use.names = FALSE)
  converters <- list(
    as_draws        = as_draws,
    as_draws_array  = as_draws_array,
    as_draws_df     = as_draws_df,
    as_draws_list   = as_draws_list,
    as_draws_matrix = as_draws_matrix
  )
  for (converter in names(converters)) {
    converted <- converters[[converter]](object)
    expect_true(
      all(correlation_names %in% posterior::variables(converted)),
      info = converter
    )
  }

  rvars <- as_draws_rvars(object)
  correlation_base <- "mu__xREx__study_xRE_CORx_R"
  expect_true(correlation_base %in% names(rvars))
  expect_identical(dim(rvars[[correlation_base]]), c(3L, 3L))
})

test_that("compiled CAR metadata drives public correlation reconstruction", {

  dat <- data.frame(
    yi    = seq(-0.2, 0.3, length.out = 6L),
    study = rep(c("s1", "s2"), each = 3L),
    time  = rep(c(0, 0.5, 2), 2L)
  )
  object <- brma.mv(
    yi                        = yi,
    V                         = diag(rep(0.04, nrow(dat))),
    random                    = ~ car(time | study),
    data                      = dat,
    measure                   = "GEN",
    prior_unit_information_sd = 1,
    only_priors               = TRUE
  )
  term <- .fitted_formula_design(object, "mu")[["random_effects"]][[1L]]
  rho  <- c(0.25, 0.81)
  values <- cbind(mu = seq_along(rho), rho)
  colnames(values)[[2L]] <- term[["correlation"]][["rho_name"]]
  fit <- coda::mcmc.list(coda::mcmc(values))
  attr(fit, "formula_design") <- list(
    mu = list(random_effects = list(term))
  )
  object[["fit"]] <- fit

  draws <- as_draws_matrix(object)
  distance <- abs(outer(c(0, 0.5, 2), c(0, 0.5, 2), "-"))
  names <- unlist(lapply(seq_len(3L), function(column) {
    paste0(
      term[["parameter_stem"]], "_xRE_CORx_R[",
      seq_len(3L), ",", column, "]"
    )
  }), use.names = FALSE)
  expect_equal(
    matrix(as.matrix(draws)[2L, names], nrow = 3L),
    0.81^distance,
    tolerance = 0
  )
})
