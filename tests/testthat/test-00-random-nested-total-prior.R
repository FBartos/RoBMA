# The total SD of a nested random structure under BMA.mv (random = ~ 1 | study
# / esid, as in the Assink scenario) is T * gate: one inclusion gate on the
# whole random part (a root allocation) and a Dirichlet split of the SD over
# the study and esid levels (its child allocation). Its exact prior is an atom
# at 0 of the gate's exclusion probability plus the scale prior T scaled by the
# gate's inclusion probability. Small single-chain fit; the expectations are
# the exact prior of that structure, computed from the fit's own priors.
.nested_total_cache <- new.env(parent = emptyenv())

.nested_total_fit <- function() {

  if (!is.null(.nested_total_cache[["fit"]])) {
    return(.nested_total_cache[["fit"]])
  }

  set.seed(3)
  k   <- 12L
  dat <- data.frame(
    study = factor(rep(sprintf("s%02d", seq_len(6L)), each = 2L)),
    esid  = factor(seq_len(k))
  )
  dat[["yi"]] <- 0.2 + stats::rnorm(6L, 0, 0.12)[as.integer(dat[["study"]])] +
    stats::rnorm(k, 0, 0.1)
  fit <- suppressWarnings(BMA.mv(
    yi = yi, V = diag(rep(0.0225, k)), random = ~ 1 | study / esid,
    data = dat, measure = "GEN", prior_unit_information_sd = 1,
    chains = 1, sample = 400, burnin = 150, adapt = 100, seed = 1,
    silent = TRUE
  ))
  .nested_total_cache[["fit"]] <- fit

  return(fit)
}

# The fit's own priors: the root allocation's scale prior T and the prior
# probability of its inclusion gate.
.nested_total_priors <- function(fit) {

  design   <- .fitted_formula_design(fit, "mu", required = TRUE)
  root     <- design[["random_allocations"]][[1L]]
  gate     <- root[["inclusion"]][[1L]][["prior"]]
  scale    <- attr(fit[["fit"]], "prior_list")[[root[["source_node"]]]]
  included <- mean(gate)

  list(
    nested   = !is.null(design[["random_allocations"]][[2L]][["parent"]]),
    included = included,
    scale    = scale,
    # the continuous part of the total SD at y > 0: P(gate = 1) f_T(y)
    density  = function(y) included * exp(BayesTools::lpdf(scale, y))
  )
}

# The layers of 'with_prior' that 'without_prior' does not have: the prior
# curve and the prior atom, as built ggplot layer data and geoms.
.plot_prior_only_layers <- function(with_prior, without_prior) {

  built_with    <- ggplot2::ggplot_build(with_prior)[["data"]]
  built_without <- ggplot2::ggplot_build(without_prior)[["data"]]
  is_shared <- vapply(built_with, function(layer) {
    any(vapply(built_without, identical, logical(1L), layer))
  }, logical(1L))
  geoms <- vapply(with_prior[["layers"]], function(layer) {
    class(layer[["geom"]])[[1L]]
  }, character(1L))

  list(data = built_with[!is_shared], geom = geoms[!is_shared])
}


test_that("the total of a nested gated allocation has its exact prior curve and atom in plots", {

  skip_on_cran()
  fit    <- .nested_total_fit()
  priors <- .nested_total_priors(fit)
  # The nested structure and the gate probability the expectations below rely on.
  expect_true(priors[["nested"]])
  expect_equal(priors[["included"]], 0.5)

  # The mixed posterior of the total carries its exact prior density: the atom
  # at 0 is the gate's exclusion probability.
  total <- .brma_random_parameter_mixed_posterior(fit, "tau_total")[[1L]]
  prior_density <- BayesTools::posterior_metadata(total, "prior_density")
  expect_equal(prior_density[["points"]][["x"]], 0)
  expect_equal(prior_density[["points"]][["p"]], 1 - priors[["included"]],
               tolerance = 1e-12)

  # The plot draws that density: no omitted-curve warning, the prior curve at
  # the scaled scale prior, and the prior atom at 0.
  without_prior <- plot(fit, "tau_total", plot_type = "ggplot")
  with_prior <- NULL
  expect_no_warning(
    with_prior <- plot(fit, "tau_total", prior = TRUE, plot_type = "ggplot"),
    class = "BayesTools_prior_curve_unavailable"
  )
  added <- .plot_prior_only_layers(with_prior, without_prior)
  expect_identical(sort(unname(added[["geom"]])), c("GeomLine", "GeomSegment"))

  curve <- added[["data"]][[which(added[["geom"]] == "GeomLine")]]
  # y > 0 drops the zero-height edge points of the curve at its support bound.
  continuous <- curve[["x"]] > 0 & curve[["y"]] > 0
  expect_gt(sum(continuous), 500L)
  expect_equal(
    curve[["y"]][continuous], priors[["density"]](curve[["x"]][continuous]),
    tolerance = 1e-8
  )

  atom <- added[["data"]][[which(added[["geom"]] == "GeomSegment")]]
  expect_equal(atom[["x"]], 0)
  expect_equal(atom[["xend"]], 0)
  expect_equal(atom[["y"]], 0)
  expect_gt(atom[["yend"]], 0)
})


test_that("point hypotheses on the total of a nested gated allocation use its exact prior", {

  skip_on_cran()
  fit    <- .nested_total_fit()
  priors <- .nested_total_priors(fit)

  # Point tests of the total and its variance are available: their prior has
  # an exact ordinate off the atom at 0.
  quantities <- hypothesis_quantities(fit)
  totals <- quantities[quantities[["parameter"]] %in%
                         c("(mu) tau_total", "(mu) tau2_total"), , drop = FALSE]
  expect_gt(nrow(totals), 0L)
  expect_true(all(totals[["point_test"]]))
  expect_true(all(totals[["direction_test"]]))

  # The prior ordinate of a point hypothesis away from the atom is the scaled
  # scale prior (the variance's by the Jacobian 1 / (2 y)).
  value <- 0.2
  total <- hypothesis(fit, "`tau_total` = 0.2", density_method = "KDE",
                      seed = 1, columns = "all")
  expect_s3_class(total, "BayesTools_hypothesis_BF")
  expect_equal(as.numeric(total[["prior"]]), priors[["density"]](value),
               tolerance = 1e-8)
  variance <- hypothesis(fit, "`tau2_total` = 0.04", density_method = "KDE",
                         seed = 1, columns = "all")
  expect_equal(as.numeric(variance[["prior"]]),
               priors[["density"]](value) / (2 * value), tolerance = 1e-8)

  # At the atom the point test is refused: the prior has a point mass at 0.
  for (parameter in c("tau_total", "tau2_total")) {
    expect_error(
      hypothesis(fit, paste0("`", parameter, "` = 0"), density_method = "KDE"),
      class = "BayesTools_point_mass_at_null",
      info  = parameter
    )
  }
})
