expect_finite_vector <- function(x, n = NULL, info = NULL) {

  testthat::expect_type(x, "double")
  if (!is.null(n)) {
    testthat::expect_equal(length(x), n, info = info)
  }
  testthat::expect_true(all(is.finite(x)), info = info)
}

expect_finite_table <- function(x, cols = NULL, n = NULL, min_cols = NULL,
                                info = NULL) {

  testthat::expect_s3_class(x, "data.frame")
  if (!is.null(cols)) {
    testthat::expect_true(all(cols %in% names(x)), info = info)
  }
  if (!is.null(n)) {
    testthat::expect_equal(nrow(x), n, info = info)
  }
  if (!is.null(min_cols)) {
    testthat::expect_true(ncol(x) >= min_cols, info = info)
  }
  numeric_cols <- vapply(x, is.numeric, TRUE)
  if (any(numeric_cols)) {
    testthat::expect_true(all(is.finite(as.matrix(x[, numeric_cols, drop = FALSE]))),
                          info = info)
  }
}

expect_positive <- function(x, strict = TRUE, info = NULL) {

  testthat::expect_true(all(is.finite(x)), info = info)
  if (strict) {
    testthat::expect_true(all(x > 0), info = info)
  } else {
    testthat::expect_true(all(x >= 0), info = info)
  }
}

expect_probability <- function(x, open = FALSE, info = NULL) {

  testthat::expect_true(all(is.finite(x)), info = info)
  if (open) {
    testthat::expect_true(all(x > 0 & x < 1), info = info)
  } else {
    testthat::expect_true(all(x >= 0 & x <= 1), info = info)
  }
}

expect_monotone <- function(x, direction = "increasing", strict = FALSE,
                            info = NULL) {

  delta <- diff(x)
  if (direction == "decreasing") {
    delta <- -delta
  }
  if (strict) {
    testthat::expect_true(all(delta > 0), info = info)
  } else {
    testthat::expect_true(all(delta >= 0), info = info)
  }
}

expect_valid_indicator <- function(x, values, info = NULL) {

  testthat::expect_true(all(is.finite(x)), info = info)
  testthat::expect_true(all(x == as.integer(x)), info = info)
  testthat::expect_true(all(x %in% values), info = info)
}

expect_error_cases <- function(cases, envir = parent.frame()) {

  for (case in cases) {
    testthat::expect_error(
      eval(case[["expr"]], envir = envir),
      regexp = case[["regexp"]],
      info   = case[["label"]]
    )
  }

  invisible(NULL)
}

expect_residual_vector <- function(x, n, info = NULL) {

  expect_finite_vector(x, n = n, info = info)
}

expect_residual_table <- function(x, n, check_se = TRUE, info = NULL) {

  expect_finite_table(x, cols = c("resid", "se", "z"), n = n, info = info)
  if (check_se) {
    expect_positive(x[["se"]], info = info)
  }
}

expect_hatvalues_vector <- function(x, n, info = NULL) {

  expect_finite_vector(x, n = n, info = info)
  testthat::expect_true(all(x >= 0 & x <= 1 + sqrt(.Machine$double.eps)),
                        info = info)
}

expect_dfbetas_table <- function(x, n, min_cols = 1, info = NULL) {

  testthat::expect_s3_class(x, "data.frame")
  testthat::expect_equal(nrow(x), n, info = info)
  testthat::expect_true(ncol(x) >= min_cols, info = info)

  numeric_cols <- vapply(x, is.numeric, TRUE)
  if (any(numeric_cols)) {
    values <- as.matrix(x[, numeric_cols, drop = FALSE])
    testthat::expect_true(all(is.finite(values) | is.nan(values)), info = info)
    if (any(is.nan(values))) {
      testthat::expect_true(!is.null(attr(x, "note")), info = info)
      testthat::expect_true(nzchar(attr(x, "note")), info = info)
    }
  }
}

expect_vif_table <- function(x, n_terms = NULL, info = NULL) {

  cols <- c("term", "df", "GVIF", "GVIF^(1/(2*df))")
  expect_finite_table(x, cols = cols, n = n_terms, info = info)
  testthat::expect_true(all(nzchar(x[["term"]])), info = info)
  testthat::expect_true(all(x[["df"]] >= 1), info = info)
  testthat::expect_true(all(x[["GVIF"]] >= 1 - sqrt(.Machine$double.eps)),
                        info = info)
  testthat::expect_true(
    all(x[["GVIF^(1/(2*df))"]] >= 1 - sqrt(.Machine$double.eps)),
    info = info
  )
}

expect_influence_object <- function(x, n, inf_cols, min_dfbs_cols = 1,
                                    info = NULL) {

  testthat::expect_s3_class(x, "infl.brma")
  testthat::expect_true(all(inf_cols %in% names(x[["inf"]])), info = info)
  testthat::expect_equal(nrow(x[["inf"]]), n, info = info)
  inf_values <- as.matrix(x[["inf"]][, inf_cols, drop = FALSE])
  testthat::expect_true(all(is.finite(inf_values) | is.nan(inf_values)),
                        info = info)
  if (any(is.nan(inf_values))) {
    testthat::expect_true(!is.null(attr(x, "note")), info = info)
    testthat::expect_true(nzchar(attr(x, "note")), info = info)
  }
  expect_dfbetas_table(x[["dfbs"]], n = n, min_cols = min_dfbs_cols,
                       info = info)
}

expect_brma_samples_matrix <- function(x, n_col, info = NULL) {

  testthat::expect_s3_class(x, "brma_samples")
  testthat::expect_true(is.matrix(x), info = info)
  testthat::expect_equal(ncol(x), n_col, info = info)
  testthat::expect_true(all(is.finite(unclass(x))), info = info)
}

expect_summary_heterogeneity_structure <- function(heterogeneity, expected_rows,
                                                   name) {

  columns <- c("Mean", "Median", "0.025", "0.975")

  testthat::expect_true(
    inherits(heterogeneity, "summary_heterogeneity.brma"),
    info = paste0("summary_heterogeneity class for '", name, "'")
  )
  testthat::expect_equal(
    sort(rownames(heterogeneity$estimates)),
    sort(expected_rows),
    info = paste0("summary_heterogeneity rows for '", name, "'")
  )

  estimates <- heterogeneity$estimates[expected_rows, columns, drop = FALSE]
  values    <- as.matrix(estimates)

  testthat::expect_true(
    all(is.finite(values)),
    info = paste0("summary_heterogeneity finite estimates for '", name, "'")
  )
  testthat::expect_true(
    all(values >= 0),
    info = paste0("summary_heterogeneity non-negative estimates for '", name, "'")
  )

  i2_rows <- grep("^I2", expected_rows, value = TRUE)
  if (length(i2_rows) > 0) {
    i2_values <- as.matrix(heterogeneity$estimates[i2_rows, columns, drop = FALSE])
    testthat::expect_true(
      all(i2_values >= 0 & i2_values <= 100),
      info = paste0("summary_heterogeneity I2 bounds for '", name, "'")
    )
  }

  if ("rho" %in% expected_rows) {
    rho_values <- as.matrix(heterogeneity$estimates["rho", columns, drop = FALSE])
    testthat::expect_true(
      all(rho_values >= 0 & rho_values <= 1),
      info = paste0("summary_heterogeneity rho bounds for '", name, "'")
    )
  }

  testthat::expect_true(
    all(heterogeneity$estimates["H2", columns] >= 1),
    info = paste0("summary_heterogeneity H2 bounds for '", name, "'")
  )
}

# Posterior draws carrying BayesTools draw metadata. 'x' gets the class
# 'class' (when given) and each field of '...' through
# BayesTools::posterior_metadata(); conditioning fields go into the 'condition'
# list.
with_draw_metadata <- function(x, ..., class = NULL) {

  if (!is.null(class)) {
    class(x) <- class
  }
  fields <- list(...)
  for (field in names(fields)) {
    BayesTools::posterior_metadata(x, field) <- fields[[field]]
  }
  x
}

# A hand-built fit ('fit' holds its draws and carries its 'prior_list'
# attribute) with the BayesTools fit contract: its parameter map, draw
# geometry and contract, as JAGS_fit() attaches them.
as_bayestools_fit <- function(fit) {

  class(fit) <- unique(c("BayesTools_fit", class(fit)))
  fit <- BayesTools:::.bt_attach_parameter_map(fit)
  fit <- BayesTools:::.bt_attach_draw_geometry(fit)
  BayesTools:::.bt_attach_fit_contract(fit)
}

# The pointwise log-likelihood that known-V fits with an estimate-level random
# intercept target: the log score of each estimate conditional on the other
# estimate of its 2 x 2 known-V dependency block,
#   y_i | y_j ~ N(mu + S_ij / S_jj (y_j - mu), S_ii - S_ij^2 / S_jj)
# with S = V + tau^2 I, evaluated analytically at every draw of 'fit' (draws
# x estimates; 'tau' is the estimate-level SD 'sd_name').
known_v_pair_conditional_log_lik <- function(fit, yi, V,
                                             sd_name = "mu__xREx__estimate_intercept") {

  draws <- as.matrix(coda::as.mcmc.list(fit[["fit"]]))
  vapply(seq_len(nrow(V)), function(i) {
    j <- setdiff(which(V[i, ] != 0), i)
    vapply(seq_len(nrow(draws)), function(s) {
      mu    <- draws[s, "mu_intercept"]
      Sigma <- V + diag(draws[s, sd_name]^2, nrow(V))
      stats::dnorm(
        yi[[i]],
        mean = mu + Sigma[i, j] / Sigma[j, j] * (yi[[j]] - mu),
        sd   = sqrt(Sigma[i, i] - Sigma[i, j]^2 / Sigma[j, j]),
        log  = TRUE
      )
    }, numeric(1))
  }, numeric(nrow(draws)))
}

# A fitted-object stand-in with the BayesTools fit contract: a root allocation
# of a gamma SD to one gated component (inclusion gate with prior probability
# 0.5), split by a child allocation over the blocks 'study' and 'esid' with
# Dirichlet(1, 1) shares, and synthetic draws: SD 1:4, root gate
# 'root_gate', and shares (.25, .75). With 'study_gate' the child allocation
# also has an inclusion gate on 'study' with these draws.
shared_gate_random_object <- function(root_gate = c(0, 1, 0, 1),
                                      study_gate = NULL) {

  mods <- data.frame(study = factor(c("a", "a", "b", "b")),
                     esid = factor(1:4))
  result <- BayesTools::JAGS_formula(
    formula = ~ 1 + random(1 | study, name = "study", covariance = "diag") +
      random(1 | esid, name = "esid", covariance = "diag"),
    parameter = "mu",
    data = mods,
    prior_list = list(intercept = BayesTools::prior("normal", list(0, 1))),
    prior_random = BayesTools::prior_random(allocation = list(
      BayesTools::random_variance_allocation(
        name = "root", terms = c(component = "component_gated"),
        sd = BayesTools::prior("gamma", list(2, 2)),
        inclusion = list(component = BayesTools::prior("spike", list(location = .5)))
      ),
      BayesTools::random_variance_allocation(
        name = "split", terms = c(study = "study", esid = "esid"),
        parent = BayesTools::allocation_ref("root", "component"),
        weights = BayesTools::prior("dirichlet", list(alpha = c(1, 1))),
        inclusion = if (!is.null(study_gate)) {
          list(study = BayesTools::prior("spike", list(location = .5)))
        }
      )
    ))
  )
  design <- result[["formula_design"]]
  root   <- design[["random_allocations"]][["root"]]
  split  <- design[["random_allocations"]][["split"]]
  weight <- split[["weight_name"]]
  samples <- cbind(
    mu_intercept = c(.1, .2, .3, .4),
    tau          = 1:4,
    gate         = root_gate,
    w1           = .25,
    w2           = .75,
    eta1         = 1,
    eta2         = 3
  )
  colnames(samples) <- c(
    "mu_intercept", root[["source_node"]],
    root[["inclusion"]][["component"]][["indicator_name"]],
    paste0(weight, "[", 1:2, "]"),
    paste0("prior_par_eta_", weight, "[", 1:2, "]")
  )
  if (!is.null(study_gate)) {
    samples <- cbind(samples, study_gate)
    colnames(samples)[ncol(samples)] <-
      split[["inclusion"]][["study"]][["indicator_name"]]
  }
  fit <- coda::mcmc.list(coda::mcmc(samples))
  attr(fit, "prior_list")     <- result[["prior_list"]]
  attr(fit, "formula_design") <- list(mu = design)

  structure(
    list(fit = as_bayestools_fit(fit), data = structure(list(), random = TRUE)),
    class = c("RoBMA", "brma.mv", "brma")
  )
}

# A fitted-object stand-in with the BayesTools fit contract whose 'mu'
# formula has a diag block 'g' (intercept and slope of 'x') with an
# SD-component allocation 'gc': a half-normal SD split by Dirichlet(1, 2)
# weights on the mean-variance scale, so each component SD is the source SD
# times sqrt(2 * weight). 'draws' maps the roles "source", "weight" (two
# columns) and "eta" (two columns; optional) to their values.
sd_component_allocation_object <- function(draws) {

  result <- BayesTools::JAGS_formula(
    formula      = ~ 1 + random(1 + x | g, name = "g", covariance = "diag"),
    parameter    = "mu",
    data         = data.frame(g = factor(c("a", "a", "b", "b")),
                              x = c(-1, 0, 1, 2)),
    prior_list   = list(intercept = BayesTools::prior("normal", list(0, 1))),
    prior_random = BayesTools::prior_random(allocation = list(
      BayesTools::random_variance_allocation(
        name = "gc", terms = "g", target = "sd_component",
        scale = "mean_variance",
        sd = BayesTools::prior("normal", list(0, 1), list(0, Inf)),
        weights = BayesTools::prior("dirichlet", list(alpha = c(1, 2)))
      )
    ))
  )
  design     <- result[["formula_design"]]
  allocation <- design[["random_effects"]][[1L]][["sd_binding"]][["allocations"]][[1L]]
  weight     <- allocation[["weight_name"]]
  samples    <- cbind(
    mu_intercept = 0,
    draws[["source"]],
    draws[["weight"]],
    draws[["eta"]]
  )
  colnames(samples) <- c(
    "mu_intercept", allocation[["source"]][["name"]],
    paste0(weight, "[", 1:2, "]"),
    if (!is.null(draws[["eta"]])) paste0("prior_par_eta_", weight, "[", 1:2, "]")
  )
  fit <- coda::mcmc.list(coda::mcmc(samples))
  attr(fit, "prior_list")     <- result[["prior_list"]]
  attr(fit, "formula_design") <- list(mu = design)

  structure(
    list(fit = as_bayestools_fit(fit), data = structure(list(), random = TRUE)),
    class = c("RoBMA", "brma.mv", "brma")
  )
}

# BayesTools posterior density and ordinate metadata built from the fields
# of a hand-written fixture list through the BayesTools constructors (the
# only accepted form); a 'status' field is the constructors' own.
as_posterior_density <- function(fields) {

  fields[["status"]] <- NULL
  do.call(BayesTools::posterior_density_attribute, fields)
}

as_posterior_ordinate <- function(fields) {

  fields[["status"]] <- NULL
  do.call(BayesTools::posterior_ordinate_attribute, fields)
}

# BayesTools label parts of a hand-built catalog quantity: the structured
# parts from which BayesTools renders the quantity's labels, with the
# quantity's canonical name as their selector.
catalog_label_parts <- function(selector, components, formula_parameter = "",
                                levels = character(), random = NULL) {

  BayesTools:::.bt_label_parts(
    components        = components,
    formula_parameter = formula_parameter,
    levels            = levels,
    random            = random,
    selector          = selector
  )
}

# The BayesTools::prior_ordinate_status() table of values whose prior
# ordinates are regular and exactly classified (eligible point hypotheses).
eligible_ordinate_status <- function(prior_density, values, labels = NULL) {

  data.frame(
    value               = values,
    eligible            = TRUE,
    condition           = NA_character_,
    reason              = NA_character_,
    continuous_behavior = "regular",
    stringsAsFactors    = FALSE
  )
}

# The 'targets' table (schema version 2) of a
# BayesTools_formula_coefficient_transform fixture: each target's map type,
# and the support of its map (positive for maps with an exp output).
formula_transform_targets <- function(map_types, output_transforms) {

  targets <- names(map_types)
  out <- data.frame(
    target            = targets,
    structural_status = "dependent",
    fixed_value       = NA_real_,
    reason            = "",
    map_type          = unname(map_types),
    stringsAsFactors  = FALSE
  )
  out[["support"]] <- lapply(targets, function(target) {
    if (identical(unname(output_transforms[target]), "exp")) {
      c(0, Inf)
    } else {
      c(-Inf, Inf)
    }
  })
  out
}
