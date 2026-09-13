.product_group_data <- function() {

  data.frame(
    yi = c(.1, .2, -.1, .3), sei = sqrt(c(.04, .06, .04, .06)),
    study = c("a", "a", "b", "b"), esid = 1:4,
    paper = c("p1", "p1", "p2", "p2"),
    crossed = c("p1", "p2", "p1", "p2"),
    incomplete = c(NA, "p1", "p2", "p2")
  )
}


.product_group_V <- function() {

  kronecker(diag(2L), matrix(c(.04, .01, .01, .06), 2L))
}


.product_group_prior <- function(model) {

  BayesTools::prior_weightfunction(
    "one-sided", steps = .025, weights = BayesTools::wf_fixed(c(1, .5)),
    model = model
  )
}


.product_group_model <- function(constructor = bselmodel.mv, ...) {

  constructor(
    yi = yi, V = .product_group_V(), random = ~ 1 | study / esid,
    data = .product_group_data(), measure = "GEN", prior_unit_information_sd = 1,
    only_priors = TRUE, silent = TRUE, ...
  )
}


test_that("product defaults need only the random formula for three-level mv models", {

  for (constructor in list(bselmodel.mv, RoBMA.mv)) {
    implicit <- .product_group_model(constructor)
    explicit <- .product_group_model(constructor, selection = selection_model())
    expect_identical(implicit$data, explicit$data)
    expect_identical(implicit$priors, explicit$priors)
    expect_identical(.create_fit_data(implicit$data, implicit$priors),
                     .create_fit_data(explicit$data, explicit$priors))
    model <- .data_selection_model(implicit$data)
    expect_true(all(vapply(model$branches[model$active_branches], function(branch) {
      is.null(branch$group) && identical(branch$weight_rule, "product")
    }, logical(1))))
    # Sampling correlation still couples both estimates within each study.
    expect_identical(.data_selection_execution_plan(implicit$data)$row_blocks,
                     list(1:2, 3:4))
  }
})


test_that("product group columns do not bind rows or change fitting inputs", {

  cells <- expand.grid(estimate_random_effects = c("condition", "integrate"),
    other_random_effects = c("condition", "integrate"),
    known_sampling_variance = c("condition", "integrate"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(cells))) {
    cell <- as.list(cells[i, , drop = FALSE])
    baseline <- .product_group_model(selection = do.call(selection_model, cell))
    expected <- .create_fit_data(baseline$data, baseline$priors)
    expected_syntax <- .create_model_syntax(baseline$data, baseline$priors)
    for (group in c("paper", "crossed", "incomplete", "not_in_data")) {
      object <- .product_group_model(
        selection = do.call(selection_model, c(cell, list(group = group))))
      expect_equal(object$data$outcome, baseline$data$outcome)
      expect_identical(.create_fit_data(object$data, object$priors), expected)
      expect_identical(.create_model_syntax(object$data, object$priors), expected_syntax)
    }
  }
  rows <- 2:4
  object <- .product_group_model(subset = rows,
    selection = selection_model(group = incomplete))
  expect_equal(object$data$outcome$yi, .product_group_data()$yi[rows])
  expect_identical(.data_selection_model(object$data)$groups$row_index, rows)
  expect_equal(.known_v_covariance_matrix(.data_known_v_data(object$data)),
               .product_group_V()[rows, rows])
})


test_that("product ensemble branches impose no publication partition agreement", {

  models <- list(selection_model(), selection_model(group = paper),
    selection_model(group = crossed), selection_model(group = not_in_data))
  object <- .product_group_model(RoBMA.mv,
    prior_bias = lapply(models, .product_group_prior))
  model <- .data_selection_model(object$data)
  expect_identical(model$branches[model$active_branches], models)
  expect_identical(.data_selection_execution_plan(object$data)$row_blocks,
                   list(1:2, 3:4))
})


test_that("only active best branches establish ensemble publication groups", {

  best <- selection_model(weight_rule = "best", group = paper)
  equivalent <- selection_model(weight_rule = "best", group = study)
  for (product in list(selection_model(), selection_model(group = not_in_data),
                       selection_model(group = crossed))) {
    object <- .product_group_model(RoBMA.mv,
      prior_bias = lapply(list(product, best, equivalent), .product_group_prior))
    model <- .data_selection_model(object$data)
    expect_identical(model$groups$row_blocks, list(1:2, 3:4))
    expect_identical(model$groups$row_labels, .product_group_data()$paper)
  }
})


test_that("zero-mass internal best branches do not activate publication binding", {

  # Public prior constructors require positive mass; exercise the preparation
  # boundary's inactive-branch contract directly with its resolved prior list.
  object <- .product_group_model(RoBMA.mv)
  inactive_best <- .product_group_prior(selection_model(weight_rule = "best"))
  inactive_best$prior_weights <- 0
  testthat::local_mocked_bindings(.selection_bias_priors = function(priors) {
    list(.product_group_prior(selection_model()), inactive_best)
  }, .package = "RoBMA")
  attr(object$data, "selection_binding") <- list(row_index = 1:4,
    input_data = .product_group_data(), cluster = NULL)
  object <- .prepare_selection_model_object(object)
  model <- .data_selection_model(object$data)
  expect_true(all(vapply(model$branches[model$active_branches], function(branch) {
    identical(branch$weight_rule, "product")
  }, logical(1))))
  expect_identical(model$groups$provenance, "inactive")
})


test_that("best still requires groups compatible with integrated covariance", {

  expect_error(.product_group_model(selection = selection_model(weight_rule = "best")),
    "Publication groups are unavailable", fixed = TRUE)
  expect_error(.product_group_model(selection = selection_model(
    weight_rule = "best", group = esid)),
    "connect publication groups", fixed = TRUE)
  expect_error(.product_group_model(selection = selection_model(
    weight_rule = "best", group = incomplete)),
    "Publication group identifiers must not be missing among retained data rows.",
    fixed = TRUE)
})


test_that("product prediction needs no publication column in newdata", {

  new_rows <- data.frame(study = c("new_a", "new_a", "new_b", "new_b"),
                         esid = 11:14, sei = .product_group_data()$sei)
  known_new <- .known_v_newdata_prepare(.product_group_V(), nrow(new_rows))
  for (constructor in list(bselmodel.mv, RoBMA.mv)) {
    object <- .product_group_model(constructor,
      selection = selection_model(group = paper))
    new_data <- .prepare_newdata(object, new_rows, type = "estimate",
      bias_adjusted = FALSE, include_scale = TRUE, include_random = TRUE)
    context <- list(object = object, same_data = FALSE, known_V_new = known_new,
      K = nrow(new_rows), raw_newdata = new_rows, outcome_data = new_data$outcome)
    expect_identical(.predict_joint_selection_groups(context), as.list(1:4))
    expect_equal(new_data$outcome$sei, new_rows$sei)
  }
  object <- .product_group_model(RoBMA.mv, prior_bias = list(
    .product_group_prior(selection_model(group = absent)),
    .product_group_prior(selection_model(weight_rule = "best", group = paper))))
  context$object <- object
  expect_error(.predict_joint_selection_groups(context),
    "The 'group' column 'paper' was not found in 'data'.", fixed = TRUE)
  context$raw_newdata$paper <- c("x", "x", "y", "y")
  expect_identical(.predict_joint_selection_groups(context), list(1:2, 3:4))
})
