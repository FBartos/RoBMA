test_that("ordered Gaussian sums retain finite supplied log weights", {

  testthat::local_mocked_bindings(
    .gauss_hermite_nodes = function(order) list(nodes = c(0, 1),
      weights = c(1, 0), log_weights = c(0, -800)),
    .iwmde_ordered_normal_completion = function(block, shares) list(
      mean = rep(0, length(shares)), sd = rep(1, length(shares)),
      log_mass = rep(0, length(shares)), max_log_scale = 0),
    .iwmde_ordered_block_log_likelihood = function(context, block, total, share) 800 * total,
    .iwmde_ordered_block_gaussian_log = function(block, total, share) rep(0, length(total)),
    .package = "RoBMA"
  )
  counts <- list2env(list(inner_calls = 0L, max_log_scale = 0))
  block <- list(counts = counts, normal = NULL, selected = TRUE,
    normal_probability = NULL, log_correction_upper = 800)
  scalar <- .iwmde_ordered_block_inner(list(), block, .5, order = 2L)
  batched <- .iwmde_ordered_block_gh_many(list(), block, c(.25, .5), order = 2L)
  expect_equal(scalar[["log_value"]], log(2), tolerance = 1e-14)
  expect_equal(vapply(batched, `[[`, numeric(1L), "log_value"), rep(log(2), 2L),
    tolerance = 1e-14)
  expect_identical(scalar[["method"]], "completed_gaussian_gh")
})


test_that("ordered normalization conditions retain their compatibility ancestor", {

  error <- tryCatch(.iwmde_ordered_normalization_error("mu", 1L,
    list(ratio_relative_error = .03)), error = identity)
  expect_s3_class(error, "RoBMA_ordered_normalization_error")
  expect_s3_class(error, "iwmde_ordered_normalization_error")
  expect_identical(error[["proof"]][["ratio_relative_error"]], .03)
  expect_identical(conditionMessage(error), paste0(
    "Ordered density for 'mu' was rejected by diagnostics: source integration ",
    "error 0.03. Inspect the condition's 'proof' field or use ",
    "'density_method = \"KDE\"'."))
})


test_that("ordered panel refusals report the observed work and supported remedy", {

  error <- tryCatch(.iwmde_ordered_positive_integral(NULL, seq(0, 1, length.out = 52L),
    function(...) stop("must not evaluate")), error = identity)
  expect_s3_class(error, "RoBMA_ordered_panel_budget_exceeded")
  expect_identical(conditionMessage(error), paste0(
    "Ordered density was rejected by diagnostics: the source partition requires ",
    "102 fine panels. Inspect the condition's 'proof' field or use ",
    "'density_method = \"KDE\"'."))
  expect_identical(error[["proof"]][["max_fine_panels"]], 100L)
})


test_that("ordinary ordered integral uncertainty does not invent a conditioning event", {

  log_terms <- rbind(log(c(1, 2, 3)), log(c(2, 3, 4)))
  density <- .iwmde_density_aggregate(log_terms, 1, 3L,
    contribution_rows = 1:3, sampling_population_rows = 1:3,
    chain_id = rep(1L, 3L), expected_chain_ids = 1L)
  execution <- list(sampling_population_rows = 1:3, population_rows = 1:4,
    conditioned_chain_id = rep(1L, 4L))
  result <- tryCatch(.iwmde_ordered_integral_error(log_terms, c(0, 1),
    density, execution, 1), error = identity)
  expect_false(inherits(result, "error"))
  if (!inherits(result, "error")) {
    expect_named(result, c("mcse", "relative_mcse", "log_mcse", "log_scale",
      "evaluation_scale_underflow"))
  }
})
