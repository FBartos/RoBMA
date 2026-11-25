test_that("missing random scale-block metadata has a targeted error", {

  object <- brma.mv(yi = c(.1, .2, .3), vi = rep(.01, 3),
                     random = ~ 1 | id, scale = ~ x,
                     data = data.frame(id = 1:3, x = 1:3),
                     measure = "SMD", only_priors = TRUE, silent = TRUE)
  design <- .fitted_formula_design(object, "mu")
  design[["random_effects"]][[1L]][["block_name"]] <- "absent"
  expect_error(.validate_prior_random_scale_sources(object[["data"]], design),
               "^Random-effect scale metadata contain an unknown block\\.$")
})

test_that("a row-shaped scale source must be bound to its targeted random block", {

  dat <- data.frame(
    yi = c(.1, .2, .3, .4), sei = rep(.2, 4),
    study = factor(c("a", "a", "b", "b")), esid = factor(1:4), x = c(0, 0, 1, 1)
  )
  point <- BayesTools::prior("point", list(.3))
  source <- BayesTools::random_sd_source("tau_esid_study", shape = "row")
  construct <- function(...) {
    brma.mv(
      yi = yi, sei = sei, random = ~ 1 | study / esid, scale = list(esid = ~ x),
      data = dat, measure = "GEN", prior_unit_information_sd = 1,
      only_priors = TRUE, silent = TRUE, ...
    )
  }
  bound_source <- function(object) {
    design <- .fitted_formula_design(object, "mu", required = TRUE)
    blocks <- vapply(design[["random_effects"]], `[[`, character(1), "block_name")
    stats::setNames(lapply(design[["random_effects"]], function(term) {
      term[["sd_binding"]][["source"]][["name"]]
    }), blocks)
  }

  # the source of the 'esid' scale model belongs to its block, not to 'study'
  wrong <- BayesTools::prior_random(
    study      = BayesTools::random_block(sd_source = source),
    esid_study = BayesTools::random_block(sd = point)
  )
  error <- tryCatch(construct(prior_heterogeneity = wrong), error = identity)
  expect_s3_class(error, "RoBMA_prior_scale_source_target")
  expect_s3_class(error, "RoBMA_prior_error")
  if (inherits(error, "error")) {
    expect_match(conditionMessage(error), "'study'", fixed = TRUE)
    expect_match(conditionMessage(error), "'esid_study'", fixed = TRUE)
  }
  swapped <- BayesTools::prior_random(
    study      = BayesTools::random_block(sd_source = source),
    esid_study = BayesTools::random_block(
      sd_source = BayesTools::random_sd_source("tau_study", shape = "row")
    )
  )
  expect_error(construct(prior_heterogeneity = swapped))

  # the targeted block keeps its source and the untargeted block its SD prior
  correct <- BayesTools::prior_random(
    study      = BayesTools::random_block(sd = point),
    esid_study = BayesTools::random_block(sd_source = source)
  )
  bound <- bound_source(construct(prior_heterogeneity = correct))
  expect_null(bound[["study"]])
  expect_identical(bound[["esid_study"]], "tau_esid_study")

  # the default prior binds the same block
  bound <- bound_source(construct())
  expect_null(bound[["study"]])
  expect_identical(bound[["esid_study"]], "tau_esid_study")
})


test_that("component scale printing limits each component frame", {

  dat <- data.frame(yi = seq(.1, .6, length.out = 6), x = 1:6,
                    study = rep(c("a", "b", "c"), each = 2), effect = 1:6)
  object <- brma.mv(yi = yi, V = diag(rep(.04, 6)),
    random = list(study = ~ 1 | study, effect = ~ 1 | effect),
    scale = list(study = ~ x, effect = ~ x), data = dat, only_data = TRUE)
  output <- capture.output(print(object[["data"]], n = 1))
  expect_true(any(grepl("$study", output, fixed = TRUE)))
  expect_true(any(grepl("$effect", output, fixed = TRUE)))
  expect_false(any(grepl("^[2-6][[:space:]]", output)))
})


test_that("partial scale keeps an independent SD root for each component", {

  dat <- data.frame(yi = c(.1, .2, .3, .4), x = 0:3,
                    study = rep(c("a", "b"), each = 2), effect = 1:4)
  args <- list(yi = dat$yi, V = diag(rep(.04, 4)),
    random = list(study = ~ 1 | study, effect = ~ 1 | effect), data = dat,
    prior_heterogeneity = BayesTools::prior("normal", list(0, 1), list(0, Inf)),
    measure = "SMD", only_priors = TRUE)
  shared <- do.call(brma.mv, args)[["priors"]][["random"]]
  partial <- do.call(brma.mv, c(args, list(scale = list(study = ~ x))))[["priors"]]
  all_scaled <- do.call(brma.mv, c(args,
    list(scale = list(study = ~ x, effect = ~ x))))[["priors"]][["random"]]
  expect_length(shared[["allocation"]], 1L)
  expect_identical(shared[["allocation"]][[1L]][["scale"]], "total_variance")
  expect_null(partial[["random"]][["allocation"]])
  expect_identical(partial[["random"]][["blocks"]][["study"]][["sd_source"]][["name"]], "tau_study")
  expect_equal(partial[["random"]][["blocks"]][["effect"]][["sd"]][["parameters"]][["sd"]], 1)
  expect_equal(partial[["scale"]][["log_tau_study"]][["intercept"]][["parameters"]][["sd"]], 1)
  expect_null(all_scaled[["allocation"]])
  expect_identical(all_scaled[["blocks"]][["study"]][["sd_source"]][["name"]], "tau_study")
  expect_identical(all_scaled[["blocks"]][["effect"]][["sd_source"]][["name"]], "tau_effect")
})
