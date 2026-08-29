context("Multivariate PET and PEESE constructors")


.pet_peese_mv_input <- function() {

  data <- data.frame(
    yi    = c(0.05, 0.16, 0.22, 0.34, 0.08, 0.19),
    study = factor(rep(c("s1", "s2", "s3"), each = 2L)),
    obs   = factor(seq_len(6L)),
    x     = rep(c(0, 1), 3L)
  )
  V <- kronecker(
    diag(3L),
    matrix(c(0.010, 0.004, 0.004, 0.014), nrow = 2L)
  )

  list(data = data, V = V)
}


test_that("bPET.mv and bPEESE.mv expose the multivariate class hierarchy", {

  input <- .pet_peese_mv_input()
  constructors <- list(
    PET = list(
      fun       = bPET.mv,
      class     = c("bPET.mv", "bPET", "brma.mv", "brma.norm", "brma"),
      predicate = BayesTools::is.prior.PET
    ),
    PEESE = list(
      fun       = bPEESE.mv,
      class     = c("bPEESE.mv", "bPEESE", "brma.mv", "brma.norm", "brma"),
      predicate = BayesTools::is.prior.PEESE
    )
  )

  expect_true(all(
    c("bPET.mv", "bPEESE.mv") %in% getNamespaceExports("RoBMA")
  ))
  for (bias_type in names(constructors)) {
    constructor <- constructors[[bias_type]]
    object <- constructor[["fun"]](
      yi = yi, V = input[["V"]], mods = ~ x,
      random = list(study = ~ 1 | study, estimate = ~ 1 | obs),
      data = input[["data"]], measure = "GEN",
      prior_unit_information_sd = 1,
      only_priors = TRUE, silent = TRUE
    )

    expect_identical(
      class(object)[-1L],
      constructor[["class"]]
    )
    expect_true(constructor[["predicate"]](
      object[["priors"]][["outcome"]][["bias"]]
    ))
    expect_true(.is_data_known_v(object[["data"]]))
    expect_true(.is_data_random(object[["data"]]))
    expect_true(inherits(
      object[["formula_design"]][["mu"]],
      "BayesTools_formula_design"
    ))
  }
})


test_that("PET and PEESE predictors are available to every known-V backend", {

  input <- .pet_peese_mv_input()
  specifications <- list(
    PET   = list(fun = bPET.mv, parameter = "PET", power = 1),
    PEESE = list(fun = bPEESE.mv, parameter = "PEESE", power = 2)
  )

  for (bias_type in names(specifications)) {
    specification <- specifications[[bias_type]]
    for (backend in c("latent", "whitened", "block_mvn")) {
      object <- specification[["fun"]](
        yi = yi, V = input[["V"]], random = ~ 1 | study,
        data = input[["data"]], measure = "GEN",
        prior_unit_information_sd = 1,
        known_v_parameterization = backend,
        only_priors = TRUE, silent = TRUE
      )
      fit_data <- .create_fit_data(
        data   = object[["data"]],
        priors = object[["priors"]]
      )
      syntax <- .create_model_syntax(
        data   = object[["data"]],
        priors = object[["priors"]]
      )

      expect_equal(fit_data[["sei"]], sqrt(diag(input[["V"]])))
      expect_match(
        syntax,
        paste0(
          specification[["parameter"]],
          " \\* ",
          if (specification[["power"]] == 1) {
            "sei\\[i\\]"
          } else {
            "pow\\(sei\\[i\\],2\\)"
          }
        )
      )
      expect_identical(
        .data_known_v_effective_backend(object[["data"]]),
        backend
      )
    }
  }
})


test_that("multivariate PET and PEESE preserve brma.mv input semantics", {

  input <- .pet_peese_mv_input()
  vi <- diag(input[["V"]])

  pet_data <- bPET.mv(
    yi = input[["data"]][["yi"]], vi = vi,
    only_data = TRUE, silent = TRUE
  )
  peese_data <- bPEESE.mv(
    yi = input[["data"]][["yi"]], sei = sqrt(vi),
    only_data = TRUE, silent = TRUE
  )

  expect_equal(
    .known_v_covariance_matrix(.data_known_v_data(pet_data[["data"]])),
    diag(vi)
  )
  expect_equal(
    .known_v_covariance_matrix(.data_known_v_data(peese_data[["data"]])),
    diag(vi)
  )
  expect_error(
    bPET.mv(
      yi = yi, V = input[["V"]], cluster = study,
      data = input[["data"]], measure = "GEN", only_data = TRUE
    ),
    "'cluster' is not supported in bPET.mv\\(\\)"
  )
  expect_error(
    bPEESE.mv(
      yi = yi, V = input[["V"]], weights = rep(1, 6L),
      data = input[["data"]], measure = "GEN", only_data = TRUE
    ),
    "'weights' are not supported in bPEESE.mv\\(\\)"
  )
})


test_that("known-V PET and PEESE bridge targets match joint normal densities", {

  input <- .pet_peese_mv_input()
  specifications <- list(
    PET = list(
      fun       = bPET.mv,
      parameter = "PET",
      value     = 0.35,
      direction = "positive",
      predictor = sqrt(diag(input[["V"]]))
    ),
    PEESE = list(
      fun       = bPEESE.mv,
      parameter = "PEESE",
      value     = 0.60,
      direction = "negative",
      predictor = diag(input[["V"]])
    )
  )
  mu  <- 0.12
  tau <- 0.08

  for (bias_type in names(specifications)) {
    specification <- specifications[[bias_type]]
    for (backend in c("latent", "whitened", "block_mvn")) {
      object <- specification[["fun"]](
        yi = yi, V = input[["V"]], data = input[["data"]],
        measure = "GEN",
        effect_direction = specification[["direction"]],
        prior_unit_information_sd = 1,
        known_v_parameterization = backend,
        only_priors = TRUE, silent = TRUE
      )
      fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
      fit_data <- .marglik_add_random_covariance_bridge_data(
        fit_data                      = fit_data,
        model_data                    = object[["data"]],
        marginalizing                = TRUE,
        sampling_latent_marginalized = TRUE
      )
      parameters <- list(mu = mu, tau = tau)
      parameters[[specification[["parameter"]]]] <- specification[["value"]]
      direction <- if (specification[["direction"]] == "negative") -1 else 1
      expected <- .marglik_mvn_log_density(
        y    = direction * input[["data"]][["yi"]],
        mean = direction * mu +
          specification[["value"]] * specification[["predictor"]],
        covariance = input[["V"]] + diag(tau^2, nrow(input[["V"]]))
      )
      actual <- .log_posterior(
        parameters                  = parameters,
        data                        = fit_data,
        is_scale                    = FALSE,
        is_multilevel               = FALSE,
        is_weights                  = FALSE,
        is_known_v                  = TRUE,
        is_PET                      = bias_type == "PET",
        is_PEESE                    = bias_type == "PEESE",
        is_weightfunction           = FALSE,
        effect_direction            = specification[["direction"]],
        outcome_type                = "norm",
        model_data                  = object[["data"]],
        sampling_latent_marginalized = TRUE
      )

      expect_equal(actual, expected, tolerance = 1e-12)
    }
  }
})


test_that("multivariate PEESE composes with scale and random formulas", {

  input <- .pet_peese_mv_input()
  object <- bPEESE.mv(
    yi = yi, V = input[["V"]], scale = ~ x, random = ~ 1 | study,
    data = input[["data"]], measure = "GEN",
    effect_direction = "negative",
    prior_unit_information_sd = 1,
    only_priors = TRUE, silent = TRUE
  )

  expect_true(.is_data_scale(object[["data"]]))
  expect_true(.is_data_random(object[["data"]]))
  expect_identical(
    .data_known_v_effective_backend(object[["data"]]),
    "block_mvn"
  )
  fit_data <- .create_fit_data(object[["data"]], object[["priors"]])
  blocks   <- .known_v_backend_blocks(
    .data_known_v_data(object[["data"]]),
    "block_mvn"
  )
  for (b in seq_along(blocks)) {
    expect_equal(
      fit_data[[paste0("known_v_y_", b)]],
      -input[["data"]][["yi"]][blocks[[b]][["index"]]]
    )
  }
})
