test_that("two-level total random effects keep the estimate names", {

  object <- brma(
    yi = c(-.2, .1, .4), sei = c(.2, .2, .2), measure = "OR",
    prior_unit_information_sd = 1, only_priors = TRUE
  )
  samples <- cbind(mu = c(-.2, 0, .2, .4), tau = c(.1, .2, .3, .4))

  estimate <- ranef(object, component = "estimate", .posterior_samples = samples)
  total    <- ranef(object, component = "total", .posterior_samples = samples)
  expect_identical(colnames(total), colnames(estimate))
  expect_false(anyDuplicated(colnames(total)) > 0L)
  expect_equal(unname(as.matrix(total)), unname(as.matrix(estimate)), tolerance = 0)

  # Expanded totals are unchanged.
  expanded <- ranef(object, component = "total", expand = TRUE, .posterior_samples = samples)
  expect_identical(colnames(expanded), colnames(estimate))
})

.ranef_routing_mv <- function(scale = NULL, drop_unused = FALSE) {

  dat <- data.frame(
    yi = c(.1, .2, .3, .4), vi = c(.04, .05, .04, .05),
    id = factor(c("a", "a", "b", "b"), levels = c("a", "b", "unused")),
    x  = c(0, 1, 0, 1)
  )
  if (drop_unused) dat[["id"]] <- droplevels(dat[["id"]])
  args <- list(
    yi = quote(yi), V = quote(vi), data = dat, random = ~1 | id,
    known_v_parameterization = "block_mvn", measure = "GEN",
    prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE
  )
  if (!is.null(scale)) args[["scale"]] <- scale
  object <- do.call(brma.mv, args)
  design <- .fitted_formula_design(object, "mu", required = TRUE)
  term   <- design[["random_effects"]][[1L]]
  list(object = object, term = term)
}

test_that("unique random effects skip declared grouping levels without rows", {

  fixture <- .ranef_routing_mv()
  expect_identical(fixture$term[["group_levels"]], c("a", "b", "unused"))
  sd_name <- fixture$term[["sd_parameter_names"]][[1L]]
  samples <- matrix(c(.05, .2, .06, .3), 2L, byrow = TRUE,
                    dimnames = list(NULL, c("mu", sd_name)))

  unique_effects <- ranef(fixture$object, .posterior_samples = samples)
  expanded       <- ranef(fixture$object, expand = TRUE, .posterior_samples = samples)
  expect_identical(colnames(unique_effects), c("u_id[a]", "u_id[b]"))
  expect_equal(unname(as.matrix(unique_effects)), unname(as.matrix(expanded))[, c(1L, 3L)],
               tolerance = 1e-12)
})

test_that("observed grouping levels follow declared order with and without metadata", {

  term <- list(group_map = c(3L, 3L, 1L), group_levels = c("a", "unused", "c"))
  expect_identical(
    .ranef_observed_group_levels(term),
    list(group_map = c(2L, 2L, 1L), group_levels = c("a", "c"))
  )
  term[["group_observed_levels"]] <- c("c", "a")
  expect_identical(
    .ranef_observed_group_levels(term),
    list(group_map = c(2L, 2L, 1L), group_levels = c("a", "c"))
  )
})
