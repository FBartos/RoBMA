context("User-facing output data-frame contract")

.output_contract_table <- function(parameters = c("alpha", "beta")) {

  data.frame(
    Mean        = seq_along(parameters),
    Median      = seq_along(parameters) + 0.1,
    `0.025`     = seq_along(parameters) - 0.2,
    `0.975`     = seq_along(parameters) + 0.2,
    row.names   = parameters,
    check.names = FALSE
  )
}


.expect_output_data_frame <- function(x, components, intervals = FALSE) {

  output      <- as.data.frame(x)
  constructor <- data.frame(x)

  expect_s3_class(output, "data.frame")
  expect_identical(constructor, output)
  expect_identical(names(output)[1:2], c("component", "parameter"))
  expect_setequal(unique(output[["component"]]), components)
  expect_type(output[["parameter"]], "character")
  expect_false(any(grepl("^[0-9]", names(output))))
  expect_false(any(grepl("^X[0-9]", names(output))))
  if (intervals) {
    expect_true(all(c("CI_0.025", "CI_0.975") %in% names(output)))
  }

  invisible(output)
}


test_that("posterior sample results and summaries satisfy the output contract", {

  samples <- .new_brma_samples(
    samples   = matrix(
      1:40,
      nrow     = 20L,
      dimnames = list(NULL, c("alpha", "beta"))
    ),
    n_chains  = 1L,
    n_iter    = 20L,
    title     = "Location",
    component = "location"
  )
  components <- .new_brma_samples_list(list(
    location = samples,
    scale    = .new_brma_samples(
      samples   = matrix(41:60, ncol = 1L, dimnames = list(NULL, "tau")),
      n_chains  = 1L,
      n_iter    = 20L,
      title     = "Scale",
      component = "scale"
    )
  ))

  samples_frame <- .expect_output_data_frame(
    samples,
    components = "location",
    intervals  = TRUE
  )
  summary_frame <- .expect_output_data_frame(
    summary(samples),
    components = "location",
    intervals  = TRUE
  )
  list_frame <- .expect_output_data_frame(
    components,
    components = c("location", "scale"),
    intervals  = TRUE
  )

  expect_identical(samples_frame, summary_frame)
  expect_identical(
    list_frame[["parameter"]],
    c("alpha", "beta", "tau")
  )
})


test_that("internal interpretation keeps the raw summary-table schema", {

  draws <- seq_len(100L)
  samples <- .new_brma_samples(
    samples   = matrix(draws, ncol = 1L, dimnames = list(NULL, "mu")),
    n_chains  = 1L,
    n_iter    = length(draws),
    title     = "Location",
    component = "location"
  )

  record <- .interpret_samples_record(
    samples      = samples,
    parameter    = "pooled effect",
    conditioning = NULL,
    probs        = c(.10, .90),
    central      = "mode"
  )

  expect_identical(record[["parameter"]], "pooled effect")
  expect_true(is.finite(record[["central_value"]]))
  expect_equal(
    c(record[["lower_value"]], record[["upper_value"]]),
    unname(stats::quantile(draws, c(.10, .90)))
  )
})


test_that("summary result families satisfy the output contract", {

  table <- .output_contract_table()
  sections <- c(
    "inclusion_components", "inclusion_mods", "inclusion_scale",
    "estimates", "estimates_conditional",
    "estimates_mods", "estimates_mods_conditional",
    "estimates_scale", "estimates_scale_conditional",
    "estimates_random", "estimates_random_conditional",
    "estimates_bias", "estimates_bias_conditional"
  )
  brma_summary <- c(
    list(name = "Model", known_v_backend = list()),
    as.list(stats::setNames(vector("list", length(sections)), sections))
  )
  brma_summary[["inclusion_mods"]] <- data.frame(
    Prior      = 0.5,
    Posterior  = 0.7,
    row.names  = "alpha",
    check.names = FALSE
  )
  brma_summary[["estimates_mods"]] <- table
  class(brma_summary) <- "summary.brma"

  heterogeneity <- structure(
    list(estimates = table, component = "location"),
    class = "summary_heterogeneity.brma"
  )
  heterogeneity_list <- structure(
    list(location = heterogeneity),
    class = c("summary_heterogeneity.brma_list", "list")
  )

  marginal_summary <- table
  class(marginal_summary) <- c(
    "summary.marginal_means.brma",
    "data.frame"
  )

  zplot_summary <- structure(
    list(estimates = table),
    class = "summary.zplot_brma"
  )

  brma_frame <- .expect_output_data_frame(
    brma_summary,
    components = c("inclusion location", "location"),
    intervals  = TRUE
  )
  testthat::local_mocked_bindings(
    summary.brma = function(object, ...) brma_summary,
    .package     = "RoBMA"
  )
  fit_frame <- .expect_output_data_frame(
    structure(list(), class = "brma"),
    components = c("inclusion location", "location"),
    intervals  = TRUE
  )
  .expect_output_data_frame(
    heterogeneity,
    components = "location",
    intervals  = TRUE
  )
  .expect_output_data_frame(
    heterogeneity_list,
    components = "location",
    intervals  = TRUE
  )
  .expect_output_data_frame(
    marginal_summary,
    components = "marginal means",
    intervals  = TRUE
  )
  .expect_output_data_frame(
    zplot_summary,
    components = "zplot",
    intervals  = TRUE
  )
  expect_true(is.na(brma_frame[brma_frame[["component"]] ==
    "inclusion location", "Mean"]))
  expect_true(is.na(brma_frame[brma_frame[["component"]] ==
    "location", "Prior"][1L]))
  expect_identical(fit_frame, brma_frame)
})


test_that("print-delegating result objects match their summary data frames", {

  marginal_summary <- .output_contract_table()
  class(marginal_summary) <- c(
    "summary.marginal_means.brma",
    "data.frame"
  )
  zplot_summary <- structure(
    list(estimates = .output_contract_table()),
    class = "summary.zplot_brma"
  )
  testthat::local_mocked_bindings(
    summary.marginal_means.brma = function(object, ...) marginal_summary,
    summary.zplot_brma          = function(object, ...) zplot_summary,
    .package = "RoBMA"
  )

  marginal_means_result <- structure(list(), class = "marginal_means.brma")
  zplot_result          <- structure(list(), class = "zplot_brma")

  expect_identical(
    as.data.frame(marginal_means_result),
    as.data.frame(marginal_summary)
  )
  expect_identical(
    data.frame(marginal_means_result),
    data.frame(marginal_summary)
  )
  expect_identical(
    as.data.frame(zplot_result),
    as.data.frame(zplot_summary)
  )
  expect_identical(
    data.frame(zplot_result),
    data.frame(zplot_summary)
  )
})


test_that("printable multi-table results satisfy the output contract", {

  marginal_models <- structure(
    list(
      name     = "Model",
      type     = "marginal",
      marginal = list(
        Effect = data.frame(
          Hypothesis = c("Null", "Alternative"),
          prior_prob = c(0.5, 0.5),
          post_prob  = c(0.4, 0.6),
          row.names  = c("null", "alternative")
        ),
        Heterogeneity = data.frame(
          Hypothesis = c("Null", "Alternative"),
          prior_prob = c(0.6, 0.4),
          post_prob  = c(0.5, 0.5),
          row.names  = c("null", "alternative")
        )
      )
    ),
    class = "summary_models.RoBMA"
  )
  individual_models <- structure(
    list(
      name       = "Model",
      type       = "individual",
      individual = data.frame(
        Effect     = c("null", "alternative"),
        prior_prob = c(0.5, 0.5),
        post_prob  = c(0.4, 0.6)
      )
    ),
    class = "summary_models.RoBMA"
  )

  dfbetas <- structure(
    data.frame(mu = c(0.1, 0.2), row.names = c("study 1", "study 2")),
    class = c("dfbetas.brma", "data.frame")
  )
  influence <- structure(
    list(
      inf = data.frame(
        hat       = c(0.1, 0.2),
        row.names = c("study 1", "study 2")
      ),
      dfbs = dfbetas
    ),
    class = "infl.brma"
  )
  vif_result <- structure(
    list(
      vif = data.frame(
        term = c("x", "y"),
        df   = c(1, 1),
        GVIF = c(1.2, 1.3)
      ),
      posterior_correlation = matrix(
        c(1, 0.2, 0.2, 1),
        nrow     = 2L,
        dimnames = list(c("x", "y"), c("x", "y"))
      )
    ),
    class = "vif.brma"
  )
  covratio_result <- structure(
    c("study 1" = 0.9, "study 2" = 1.1),
    class = c("covratio.brma", "numeric")
  )

  .expect_output_data_frame(
    marginal_models,
    components = c("Effect", "Heterogeneity")
  )
  individual_frame <- .expect_output_data_frame(
    individual_models,
    components = "individual models"
  )
  .expect_output_data_frame(dfbetas, components = "dfbetas")
  .expect_output_data_frame(
    influence,
    components = c("influence", "dfbetas")
  )
  vif_frame <- .expect_output_data_frame(
    vif_result,
    components = c("vif", "posterior correlation")
  )
  .expect_output_data_frame(covratio_result, components = "covratio")

  expect_identical(
    individual_frame[["parameter"]],
    c("model 1", "model 2")
  )
  expect_identical(
    names(vif_frame)[1:3],
    c("component", "parameter", "VIF")
  )
})
