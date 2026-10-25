test_that("a random inclusion gate belongs only to its declared block", {

  result <- BayesTools::JAGS_formula(
    formula = ~ 1 + random(1 | study, name = "study", covariance = "diag") +
      random(1 | esid, name = "esid", covariance = "diag"),
    parameter = "mu",
    data = data.frame(study = factor(c("a", "a", "b", "b")), esid = factor(1:4)),
    prior_list = list(intercept = BayesTools::prior("normal", list(0, 1))),
    prior_random = BayesTools::prior_random(
      sd = BayesTools::prior("gamma", list(2, 2)),
      allocation = list(BayesTools::random_variance_allocation(
        name = "gated", terms = c(study = "study"),
        sd = BayesTools::prior("gamma", list(2, 2)),
        inclusion = list(study = BayesTools::prior("spike", list(location = .5)))
      ))
    )
  )
  design <- result[["formula_design"]]
  allocation <- design[["random_allocations"]][["gated"]]
  source <- allocation[["source_node"]]
  gate <- allocation[["inclusion"]][["study"]][["indicator_name"]]
  sd_name <- design[["random_effects"]][[2L]][["sd_parameter_names"]]
  samples <- cbind(
    mu_intercept = c(-.1, .1, .2, .3), source = c(.2, .3, .4, .5),
    gate = c(0, 1, 0, 1), sd = c(.6, .7, .8, .9)
  )
  colnames(samples)[2:4] <- c(source, gate, sd_name)
  fit <- coda::mcmc.list(coda::mcmc(samples))
  class(fit) <- c("BayesTools_fit", class(fit))
  attr(fit, "prior_list") <- result[["prior_list"]]
  attr(fit, "formula_design") <- list(mu = design)
  fit <- BayesTools:::.bt_attach_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  fit <- BayesTools:::.bt_attach_fit_contract(fit)
  object <- structure(list(
    fit = fit, data = structure(list(), random = TRUE)
  ), class = c("RoBMA", "brma.mv", "brma"))

  for (quantity in c("tau", "tau2")) {
    study_parameter <- paste0("study: ", quantity)
    esid_parameter <- paste0("esid: ", quantity)

    study_posterior <- .brma_random_parameter_mixed_posterior(object, study_parameter)[[1L]]
    esid_posterior <- .brma_random_parameter_mixed_posterior(object, esid_parameter)[[1L]]
    study_atoms <- BayesTools::posterior_metadata(study_posterior, "atoms")
    esid_atoms <- BayesTools::posterior_metadata(esid_posterior, "atoms")
    expect_equal(as.numeric(study_atoms[["locations"]]), 0)
    expect_equal(study_atoms[["mass"]], .5)
    expect_length(esid_atoms[["mass"]], 0L)
    expect_equal(as.numeric(esid_posterior), if (quantity == "tau") {
      samples[, sd_name]
    } else {
      samples[, sd_name]^2
    })

    # Conditioning keeps the draws with the study gate on; the ungated esid
    # block has no inclusion event.
    study_conditional <- .brma_random_parameter_mixed_posterior(
      object, study_parameter, conditional = TRUE
    )[[1L]]
    expect_length(study_conditional, 2L)
    expect_true(BayesTools::posterior_atoms_free(study_conditional))
    expect_error(
      .brma_random_parameter_mixed_posterior(object, esid_parameter, conditional = TRUE),
      "the quantity has no inclusion gate"
    )
  }
})


test_that("variance proportions require a possible positive parent allocation", {

  # Parent (root) gate (0, 0, 1, 1) and child study gate (0, 1, 0, 1): the
  # study share is defined only where the parent allocation is included,
  # also where the child gate is on.
  object  <- shared_gate_random_object(
    root_gate  = c(0, 0, 1, 1),
    study_gate = c(0, 1, 0, 1)
  )
  selected <- .brma_random_parameter_select(object, "(mu) split: tau2_prop(study)")
  state <- .iwmde_gate_states(
    list(object = object, posterior_samples = as.matrix(object[["fit"]][[1L]])),
    list(gate_selection = .brma_random_parameter_gate_selection(object, selected))
  )
  expect_identical(state[["defined"]], c(FALSE, FALSE, TRUE, TRUE))
  expect_identical(state[["point_zero"]], c(FALSE, FALSE, TRUE, FALSE))
  expect_identical(state[["continuous"]], c(FALSE, FALSE, FALSE, TRUE))
})


test_that("conditional random-effect plots refuse qCMDE/IWMDE with the density-method classes", {

  # A known-V stand-in passes the qCMDE/IWMDE capability check of the model;
  # its conditional random-effect plots are KDE-only.
  object <- shared_gate_random_object()
  attr(object[["data"]], "known_V") <- TRUE
  attr(object[["data"]], "measure") <- "GEN"
  expect_true(.iwmde_capability(
    object = object, density_method = "IWMDE"
  )[["available"]])
  for (method in c("qCMDE", "IWMDE")) {
    for (class in c("RoBMA_density_method_conditional_random",
                    "RoBMA_density_method_unavailable")) {
      expect_error(
        plot(object, "study: tau", component = "random", conditional = TRUE,
             density_method = method, plot_type = "ggplot"),
        class = class,
        info  = method
      )
    }
    # hypothesis() refuses the same cause with its own classes followed by
    # the density-method classes.
    error <- tryCatch(
      hypothesis(object, "`(mu) study: tau(intercept)` = 0.3",
                 component = "random", conditional = TRUE,
                 density_method = method, n_samples = 1000L, seed = 1),
      error = identity
    )
    expect_identical(
      class(error),
      c("RoBMA_hypothesis_method", "RoBMA_hypothesis_unavailable",
        "RoBMA_density_method_conditional_random",
        "RoBMA_density_method_unavailable", "error", "condition"),
      info = method
    )
  }

  # Without known V the model does not support qCMDE/IWMDE: both entry
  # points report that cause first, not the conditional statement.
  unknown_v  <- shared_gate_random_object()
  capability <- .iwmde_capability(object = unknown_v, density_method = "qCMDE")
  for (method in c("qCMDE", "IWMDE")) {
    plotted <- tryCatch(
      plot(unknown_v, "study: tau", component = "random", conditional = TRUE,
           density_method = method, plot_type = "ggplot"),
      error = identity
    )
    error <- tryCatch(
      hypothesis(unknown_v, "`(mu) study: tau(intercept)` = 0.3",
                 component = "random", conditional = TRUE,
                 density_method = method, n_samples = 1000L, seed = 1),
      error = identity
    )
    expect_identical(
      class(plotted),
      c("RoBMA_density_method_random_unknown_v",
        "RoBMA_density_method_unavailable", "error", "condition"),
      info = method
    )
    expect_identical(
      class(error),
      c("RoBMA_hypothesis_method", "RoBMA_hypothesis_unavailable",
        class(plotted)),
      info = method
    )
    expect_identical(conditionMessage(error), capability[["reason"]], info = method)
  }
})


test_that("random-effect quantities without a scalar coordinate refuse qCMDE/IWMDE plots with the density-method classes", {

  # The variance of a random component has no scalar random-component
  # coordinate for a qCMDE/IWMDE density curve.
  object <- shared_gate_random_object()
  attr(object[["data"]], "known_V") <- TRUE
  attr(object[["data"]], "measure") <- "GEN"
  for (method in c("qCMDE", "IWMDE")) {
    error <- tryCatch(
      plot(object, "(mu) study: tau2(intercept)", component = "random",
           density_method = method, plot_type = "ggplot"),
      error = identity
    )
    expect_identical(
      class(error),
      c("RoBMA_density_method_random_target",
        "RoBMA_density_method_unavailable", "error", "condition"),
      info = method
    )
    expect_match(
      conditionMessage(error),
      "because it has no supported scalar random-component coordinate",
      fixed = TRUE,
      info  = method
    )
  }
})

# Genuine public constructors and current compiler metadata; deterministic source
# replay verifies ownership and algebra, rather than posterior calibration.
.allocation_inclusion_replay <- function(mode, unnamed = character()) {

  point <- function(x) BayesTools::prior("point", list(location = x))
  if (mode == "common") {
    dat <- data.frame(yi = c(.10, .25, .15, .30),
      study = factor(c("s1", "s1", "s2", "s2")), group = factor(c("a", "b", "a", "b")))
    object <- BMA.mv(yi = yi, V = diag(.04, 4), data = dat,
      random = ~ us(0 + group | study),
      prior_heterogeneity = BayesTools::prior_random(
        study = BayesTools::random_block(contrasts = c(group = "independent"))),
      measure = "GEN", prior_unit_information_sd = 1, only_priors = TRUE, silent = TRUE)
  } else {
    dat <- data.frame(yi = rep(c(.1, .2), 4), study = factor(rep(c("s1", "s2"), each = 4)),
      lab = factor(rep(c("l1", "l1", "l2", "l2"), 2)), site = factor(rep(c("t1", "t2"), 4)))
    random <- list(study = ~ 1 | study, lab = ~ 1 | lab, site = ~ 1 | site)
    if (mode == "gate_only") random$site <- NULL
    root <- BayesTools::random_variance_allocation(name = "root",
      terms = if (mode == "gate_only") c(pool = "pool") else if ("root" %in% unnamed) c("pool", "site") else c(pool = "pool", site = "site"),
      sd = BayesTools::prior("gamma", list(shape = 2, rate = 2)),
      weights = if (mode == "gate_only") NULL else BayesTools::prior("dirichlet", list(alpha = c(1, 1))),
      inclusion = if (mode %in% c("root", "gate_only")) list(pool = point(.5)) else NULL)
    split <- BayesTools::random_variance_allocation(name = "split",
      terms = if ("split" %in% unnamed) c("study", "lab") else c(study = "study", lab = "lab"),
      parent = BayesTools::allocation_ref("root", "pool"),
      weights = BayesTools::prior("dirichlet", list(alpha = c(1, 1))),
      inclusion = if (mode == "children") list(study = point(.5), lab = point(.5)) else if (mode == "lab") list(lab = point(.5)) else NULL)
    object <- brma.mv(yi = yi, V = diag(.01, 8), random = random, data = dat, mods = ~ 1,
      measure = "GEN", prior_unit_information_sd = 1,
      prior_effect = BayesTools::prior("normal", list(mean = 0, sd = 1)),
      prior_heterogeneity = BayesTools::prior_random(allocation = list(root, split)),
      standardize_continuous_predictors = FALSE, only_priors = TRUE, silent = TRUE)
  }
  setup <- review_compile_object(object)
  compiled <- setup$compiled$mu
  allocations <- compiled$formula_design$random_allocations
  gates <- unlist(lapply(allocations, function(a) vapply(a$inclusion, `[[`, character(1), "indicator_name")), use.names = FALSE)
  auxiliaries <- unlist(lapply(allocations, function(a) {
    if (length(a$weight_name)) .iwmde_simplex_auxiliary_columns(a$weight_name, if (a$target == "sd_component") length(a$weights$parameters$alpha) else a$n_targets) else character()
  }), use.names = FALSE)
  correlation <- unlist(lapply(compiled$formula_design$random_effects, function(term) term$correlation$primitive_names), use.names = FALSE)
  columns <- unique(c(review_prior_columns(compiled$prior_list), gates, auxiliaries, correlation))
  samples <- matrix(.5, 24L, length(columns), dimnames = list(NULL, columns))
  samples[, review_coefficient_name(setup, "intercept")] <- seq(-.1, .1, length.out = 24)
  root <- allocations[[1L]]
  samples[, root$source$name] <- seq(.8, 1.2, length.out = 24)
  weights <- cbind(rep(c(.25, .5, .75), 8), rep(c(.75, .5, .25), 8))
  for (a in allocations) if (length(a$weight_name)) {
    stopifnot((if (a$target == "sd_component") length(a$weights$parameters$alpha) else a$n_targets) == 2L)
    samples[, paste0(a$weight_name, "[", 1:2, "]")] <- weights
    samples[, .iwmde_simplex_auxiliary_columns(a$weight_name, 2L)] <- 2 * weights
  }
  if (mode == "children") {
    samples[, allocations$split$inclusion$study$indicator_name] <- rep(c(0, 1, 0, 1), 6)
    samples[, allocations$split$inclusion$lab$indicator_name] <- rep(c(0, 0, 1, 1), 6)
  } else samples[, gates] <- if (mode == "lab") rep(c(0, 0, 1, 1), 6) else rep(0:1, 12)
  list(object = review_attach_replay(setup, samples), sources = samples,
    allocations = allocations, gates = gates)
}

for (.unnamed_case in c("lab", "children", "root", "recursive")) {
  local({
    case <- .unnamed_case
    test_that(paste("allocation inclusion follows declared positions for unnamed terms:", case), {

      mode <- if (case == "recursive") "root" else case
      unnamed <- switch(case, lab = "split", children = "split", root = "root",
        recursive = c("root", "split"))
      named <- .allocation_inclusion_replay(mode)
      replay <- .allocation_inclusion_replay(mode, unnamed)
      original <- serialize(replay$object, NULL)
      expect_identical(replay$sources, named$sources)
      for (owner in unnamed) {
        expect_null(names(replay$allocations[[owner]]$terms))
        expect_identical(unname(replay$allocations[[owner]]$component_labels),
          unname(named$allocations[[owner]]$component_labels))
      }
      expected <- switch(mode, lab = "lab: tau", children = c("study: tau", "lab: tau"),
        root = "split: tau_total")
      quantities <- BayesTools::parameter_catalog(replay$object$fit)$quantities
      rows <- which(quantities$role == "random_inclusion")
      named_public <- as.matrix(posterior::as_draws_matrix(as_draws(named$object)))
      public <- tryCatch(as.matrix(posterior::as_draws_matrix(as_draws(replay$object))), error = identity)
      labels <- tryCatch(.summary_random_inclusion_labels(replay$object,
        quantities$canonical_name[rows], quantities$canonical_name[rows]), error = identity)
      expect_identical(labels, expected)
      expect_false(inherits(public, "error"))
      if (!inherits(public, "error")) {
        expect_identical(public, named_public)
        expect_identical(nrow(public), 24L)
        expect_identical(anyDuplicated(colnames(public)), 0L)
        for (i in seq_along(replay$gates)) {
          expect_identical(as.numeric(public[, paste0(expected[[i]], "_indicator")]),
            unname(replay$sources[, replay$gates[[i]]]))
        }
        expect_false(any(replay$gates %in% colnames(public)))
      }
      result <- summary(replay$object)
      expect_s3_class(result, "summary.brma")
      expect_length(result$inclusion_random, 0L)
      expect_identical(serialize(replay$object, NULL), original)
    })
  })
}

for (.inclusion_mode in c("children", "root", "gate_only", "lab", "common")) {
  local({
    mode <- .inclusion_mode
    test_that(paste("allocation inclusion uses exact declared scalar SD ownership:", mode), {

      replay <- .allocation_inclusion_replay(mode)
      object <- replay$object
      quantities <- BayesTools::parameter_catalog(object$fit)$quantities
      rows <- which(quantities$role == "random_inclusion")
      expected <- switch(mode, children = c("study: tau", "lab: tau"),
        root = "split: tau_total", gate_only = "split: tau_total", lab = "lab: tau", common = "component_1: tau_common")
      mapping <- .random_inclusion_sd_names(object)
      expect_identical(unname(mapping[replay$gates]), expected)
      expect_identical(.summary_random_inclusion_labels(object,
        quantities$canonical_name[rows], quantities$canonical_name[rows]), expected)
      internal <- as.matrix(BayesTools::JAGS_materialize_draws(object$fit, include_internal = TRUE))
      public <- as.matrix(posterior::as_draws_matrix(as_draws(object)))
      # the component SDs derived from an allocation are auxiliary draws
      public_all <- as.matrix(posterior::as_draws_matrix(
        as_draws(object, include_auxiliary = TRUE)
      ))
      for (i in seq_along(replay$gates)) {
        gate <- replay$gates[[i]]
        label <- paste0(expected[[i]], "_indicator")
        expect_identical(sum(colnames(public) == label), 1L)
        expect_identical(unname(internal[, gate]), unname(replay$sources[, gate]))
        expect_identical(as.numeric(public[, label]), unname(replay$sources[, gate]))
        expect_setequal(unname(public[, label]), c(0, 1))
      }
      expect_false(any(replay$gates %in% colnames(public)))
      if (mode == "lab") {
        key <- quantities$extraction_key[[rows]]
        expect_identical(key$index, 1L)
        expect_identical(replay$allocations$split$inclusion$lab$index, 2L)
        inactive <- replay$sources[, replay$gates] == 0
        expect_identical(sum(inactive), 12L)
        expect_true(all(public_all[inactive, "lab: tau"] == 0))
        expect_true(all(public_all[inactive, "study: tau"] > 0))
        expect_true(all(public[inactive, "split: tau_total"] > 0))
        expect_false("split: tau_total_indicator" %in% colnames(public))
      }
      if (mode != "common") {
        expect_equal(unname(public[, "split: tau_total"]^2),
          unname(public_all[, "study: tau"]^2 + public_all[, "lab: tau"]^2))
        result <- summary(object)
        expect_s3_class(result, "summary.brma")
        # These genuine brma carriers do not exercise RoBMA's Component Inclusion table.
        expect_length(result$inclusion_random, 0L)
        expect_true(all(.brma_random_parameter_io_labels(quantities[rows, , drop = FALSE], "label") %in% rownames(result$estimates_random)))
      } else {
        child <- replay$allocations$study_components
        expect_identical(child$target, "sd_component")
        expect_identical(child$scale, "mean_variance")
        expect_identical(child$parent$allocation, replay$allocations[[1L]]$label)
      }
    })
  })
}

