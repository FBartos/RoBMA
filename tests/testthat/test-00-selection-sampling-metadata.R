test_that("ordinary vcalc covariance follows retained data rows and explicit publication groups", {

  skip_if_not_installed("metafor")
  dat <- data.frame(
    yi = seq(-.2, .3, length.out = 6L),
    vi = c(.02, .08, .03, .07, .04, .06),
    study = c("a", "a", "a", "a", "b", "b"),
    cohort = c(1, 1, 2, 2, 1, 1),
    obs = c(1, 1, 2, 3, 1, 2)
  )
  V <- metafor::vcalc(
    vi, cluster = study, subgroup = cohort, obs = obs,
    rho = .5, data = dat, checkpd = FALSE
  )
  rows <- c(1L, 3L, 4L, 5L, 6L)
  object <- bselmodel.mv(
    yi = yi, V = V, data = dat, subset = rows,
    selection = selection_model(group = study), measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE
  )
  known_V <- .data_known_v_data(object$data)
  metadata <- .known_v_selection_metadata(known_V)
  binding <- .data_selection_model(object$data)$groups
  expect_identical(metadata$origin, "matrix")
  expect_identical(metadata$row_index, rows)
  expect_identical(binding$row_index, rows)
  expect_identical(binding$row_labels, dat$study[rows])
  expect_identical(binding$provenance, "explicit")
  expect_equal(.known_v_covariance_matrix(known_V), unclass(V[rows, rows]))

  # Subsetting an ordinary matrix creates new local input-row identities;
  # the model's common data filter above retains the original data-row map.
  reordered <- c(5L, 3L, 1L)
  subset <- .known_v_subset_input(V, reordered)
  known_subset <- .known_v_canonicalize(subset, warn_singular = FALSE)
  expect_equal(.known_v_covariance_matrix(known_subset), unclass(V[reordered, reordered]))
  expect_identical(.known_v_selection_metadata(known_subset)$row_index, seq_along(reordered))
  expect_identical(.known_v_selection_metadata(.known_v_canonicalize_newdata(subset)),
                   .known_v_selection_metadata(known_subset))
  expect_error(bselmodel.mv(
    yi = yi, V = V, data = dat, subset = rows, measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE
  ), "Publication groups are unavailable for this input. Specify 'group' in 'selection_model()'.",
     fixed = TRUE)
})


test_that("publication binding follows declared precedence and the original row map", {

  dat <- data.frame(
    vi = c(.02, .03, .04, .05), study = c("a", "a", "b", "b"),
    paper = c(NA, "p", "q", "p")
  )
  rows <- c(4L, 2L, 3L)
  automatic <- BayesTools::selection_model()
  explicit <- BayesTools::selection_model(group = paper)

  bound <- .selection_bind_groups(explicit, rows, dat, cluster = dat$study)
  expect_identical(bound$requested, "paper")
  expect_identical(bound$row_index, rows)
  expect_identical(bound$row_labels, c("p", "p", "q"))
  expect_identical(bound$group_index, c(1L, 1L, 2L))
  expect_identical(bound$row_blocks, list(1:2, 3L))
  expect_identical(bound$provenance, "explicit")
  bound <- .selection_bind_groups(automatic, rows, dat, cluster = dat$study)
  expect_identical(bound$provenance, "cluster")
  expect_identical(bound$row_labels, dat$study[rows])
  expect_identical(bound$row_blocks, list(c(1L, 3L), 2L))
  expect_identical(
    .selection_bind_groups(automatic, rows, dat, allow_singletons = TRUE)$row_blocks,
    list(1L, 2L, 3L)
  )
  expect_error(.selection_bind_groups(automatic, rows, dat),
    "Publication groups are unavailable for this input. Specify 'group' in 'selection_model()'.",
    fixed = TRUE)
  expect_error(.selection_bind_groups(explicit, 1:4, dat),
    "Publication group identifiers must not be missing among retained data rows.", fixed = TRUE)
  paper <- rep("global", 4L)
  expect_error(.selection_bind_groups(explicit, rows, dat[c("vi", "study")]),
    "The 'group' column 'paper' was not found in 'data'.", fixed = TRUE)

  restored <- unserialize(serialize(explicit, NULL))
  dat$paper <- letters[1:4]
  expect_identical(.selection_bind_groups(restored, rows, dat)$row_labels, letters[rows])
})

test_that("whole sampling error includes diagonal and singleton factor variation", {

  diagonal <- c(.01, 0, .03)
  loading <- diag(c(.2, .3, 0))
  input <- known_v_factor(diagonal, loading)
  known_V <- .known_v_prepare(input, rep(TRUE, 3L), "auto")
  expect_identical(.known_v_effective_backend(known_V), "diagonal")
  expect_equal(known_V$residual_variance, diagonal + rowSums(loading^2))
  resolved <- .known_v_resolve_selection_structure(known_V)
  sampling <- .known_v_check_selection_structure(resolved)
  expect_identical(sampling$residual_variance, numeric(3L))
  expect_identical(sampling$source_ids, "sampling_error")
  expect_identical(sampling$provenance$policy, "whole_sampling_error")
  retained <- matrix(0, 3L, sampling$rank)
  for (block in sampling$latent_blocks) {
    retained[block$index, seq.int(block$z_start, block$z_end)] <- block$B
  }
  expect_equal(tcrossprod(retained),
    diag(diagonal) + tcrossprod(loading), tolerance = 1e-15)
  expect_identical(.known_v_resolve_selection_structure(resolved), resolved)

  rows <- c(3L, 1L)
  subset <- .known_v_subset_input(input, rows)
  subset <- .known_v_canonicalize(subset)
  restored <- .known_v_as_input(subset)
  expect_identical(.known_v_factor_metadata(restored)$row_index, rows)
  expect_identical(.known_v_factor_metadata(restored)$source_ids,
    .known_v_factor_metadata(input)$source_ids)
  expect_identical(restored$loading, loading[rows, , drop = FALSE])
  stale <- input
  stale$loading[1L, 1L] <- .25
  expect_error(.known_v_canonicalize(stale),
    "The 'V' factor metadata no longer match its declaration.", fixed = TRUE)
  resolved$selection_structure$residual_variance[[1L]] <- .02
  expect_error(.known_v_check_selection_structure(resolved),
    "The conditional sampling structure is missing or no longer matches the retained rows.",
    fixed = TRUE)
})


test_that("whole sampling covariance is independent of the Gaussian backend", {

  V <- matrix(c(1, .4, .4, 2), 2L)
  resolved <- lapply(c("whitened", "block_mvn", "latent"), function(backend) {

    .known_v_resolve_selection_structure(.known_v_prepare(V, c(TRUE, TRUE), backend))
  })
  structures <- lapply(resolved, .known_v_check_selection_structure)
  for (sampling in structures) {
    expect_identical(sampling$residual_variance, numeric(2L))
    expect_equal(tcrossprod(sampling$latent_blocks[[1L]]$B),
      V, tolerance = 1e-14)
    expect_identical(sampling$provenance$policy, "whole_sampling_error")
  }
  capped <- matrix(c(1, .95, .95, 1), 2L)
  expect_warning(capped <- .known_v_resolve_selection_structure(.known_v_canonicalize(capped)), NA)
  sampling <- .known_v_check_selection_structure(capped)
  expect_identical(sampling$residual_variance, numeric(2L))
  expect_equal(tcrossprod(sampling$latent_blocks[[1L]]$B),
    matrix(c(1, .95, .95, 1), 2L), tolerance = 1e-14)
  rank_one <- tcrossprod(c(1, 2))
  sampling <- .known_v_check_selection_structure(.known_v_resolve_selection_structure(
    .known_v_canonicalize(rank_one, warn_singular = FALSE)
  ))
  expect_identical(sampling$residual_variance, c(0, 0))
  expect_identical(sampling$rank, 1L)
  expect_equal(tcrossprod(sampling$latent_blocks[[1L]]$B), rank_one, tolerance = 0)

  singular <- tcrossprod(rbind(c(1, 0), c(0, 1), c(1, 1)))
  expect_warning(integrated <- .known_v_canonicalize(singular, warn_singular = FALSE), NA)
  expect_null(integrated$selection_structure)
  expect_warning(resolved <- .known_v_resolve_selection_structure(integrated), NA)
  sampling <- .known_v_check_selection_structure(resolved)
  expect_identical(sampling$residual_variance, numeric(3L))
  expect_equal(tcrossprod(sampling$latent_blocks[[1L]]$B), singular,
    tolerance = 1e-14)
})


test_that("whole sampling resolution preserves provenance without decomposition warnings", {

  expect_warning(plain <- .known_v_resolve_selection_structure(.known_v_canonicalize(
    matrix(c(.04, .02, .02, .09), 2L)
  )), NA)
  expect_identical(plain$selection_structure$provenance$input_origin, "matrix")
  diagonal <- .known_v_resolve_selection_structure(.known_v_canonicalize(c(.04, .09)))
  expect_identical(diagonal$selection_structure$residual_variance, numeric(2L))
  expect_identical(diagonal$selection_structure$rank, 2L)

  skip_if_not_installed("metafor")
  dat <- data.frame(vi = c(.04, .09), study = c("a", "a"), obs = 1:2, type = 1:2)
  inputs <- list(
    metafor::vcalc(vi, cluster = study, obs = obs, rho = .5, data = dat),
    metafor::vcalc(vi, cluster = study, type = type, obs = obs,
      rho = c(.6, .3), data = dat)
  )
  for (V in inputs) {
    expect_warning(resolved <- .known_v_resolve_selection_structure(
      .known_v_canonicalize(V)), NA)
    expect_identical(resolved$selection_structure$provenance$input_origin, "matrix")
    expect_identical(resolved$selection_structure$provenance$factor_status, "undeclared")
    expect_identical(resolved$selection_structure$residual_variance, numeric(2L))
    expect_equal(tcrossprod(resolved$selection_structure$latent_blocks[[1L]]$B),
      matrix(as.numeric(V), nrow(V)), tolerance = 1e-14)
  }
})

test_that("selection summaries print tables without configuration descriptions", {

  make_model <- function(selection = selection_model()) {

    object <- bselmodel(
      yi = c(.1, .2), sei = c(.2, .3), data = data.frame(paper = c("a", "a")),
      selection = selection, measure = "GEN", prior_unit_information_sd = 1,
      only_priors = TRUE, silent = TRUE
    )
    .data_selection_model(object$data)
  }
  models <- list(
    make_model(),
    make_model(selection_model(estimate_random_effects = "condition")),
    make_model(selection_model(estimate_random_effects = "condition",
      other_random_effects = "integrate", known_sampling_variance = "condition",
      weight_rule = "best", group = paper))
  )
  output <- structure(list(
    name = "Selection model",
    inclusion_components = matrix(.7, 1L, 1L, dimnames = list("Effect", "PIP")),
    estimates = matrix(c(.1, .2), 2L, 1L, dimnames = list(c("mu", "tau"), "Mean")),
    estimates_bias = matrix(.5, 1L, 1L, dimnames = list("omega[0.025,1]", "Mean"))
  ), class = "summary.brma")
  expected <- capture.output(print(output))
  expect_true(all(vapply(c("Effect", "PIP", "mu", "tau", "omega"), function(label) {
    any(grepl(label, expected, fixed = TRUE))
  }, logical(1L))))
  for (model in models) {
    output$selection_model <- model
    actual <- capture.output(returned <- print(output))
    expect_identical(actual, expected)
    expect_identical(returned$selection_model, model)
  }
})
