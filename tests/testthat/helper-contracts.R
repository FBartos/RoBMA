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

# A fitted-object stand-in with a gated total-variance allocation of a
# half-normal SD over 'study' (inclusion gate with prior probability 0.5) and
# 'esid' (ungated), Dirichlet(1, 1) shares, and synthetic draws.
gated_random_object <- function(n = 200L) {

  result <- BayesTools::JAGS_formula(
    formula = ~ 1 + random(1 | study, name = "study", covariance = "diag") +
      random(1 | esid, name = "esid", covariance = "diag"),
    parameter = "mu",
    data = data.frame(study = factor(c("a", "a", "b", "b")), esid = factor(1:4)),
    prior_list = list(intercept = BayesTools::prior("normal", list(0, 1))),
    prior_random = BayesTools::prior_random(
      sd = BayesTools::prior("gamma", list(2, 2)),
      allocation = list(BayesTools::random_variance_allocation(
        name = "split", terms = c(study = "study", esid = "esid"),
        sd = BayesTools::prior("normal", list(0, 1), list(0, Inf)),
        inclusion = list(study = BayesTools::prior("spike", list(location = .5)))
      ))
    )
  )
  share <- seq(0.02, 0.98, length.out = n)
  samples <- cbind(
    mu_intercept = rep(c(-0.1, 0.1), length.out = n),
    mu__xRE_ALLOCx_split__allocation_sd = 0.3 + 0.4 * share,
    "mu__xRE_ALLOCx_split__weight[1]" = share,
    "mu__xRE_ALLOCx_split__weight[2]" = 1 - share,
    mu__xRE_ALLOCx_split__include_study_indicator = rep(c(0, 1), length.out = n)
  )
  fit <- coda::mcmc.list(coda::mcmc(samples))
  attr(fit, "prior_list") <- result[["prior_list"]]
  attr(fit, "formula_design") <- list(mu = result[["formula_design"]])

  structure(
    list(fit = as_bayestools_fit(fit), data = structure(list(), random = TRUE)),
    class = c("RoBMA", "brma.mv", "brma")
  )
}

# A fitted-object stand-in with one ungated random-intercept SD with prior
# 'prior' and constant draws 0.2 (a sampled quantity whose retained draws
# happen to be constant, or a structural one for a spike prior).
single_sd_random_object <- function(prior) {

  result <- BayesTools::JAGS_formula(
    formula = ~ 1 + random(1 | study, name = "study", covariance = "diag"),
    parameter = "mu",
    data = data.frame(study = factor(c("a", "a", "b", "b"))),
    prior_list = list(intercept = BayesTools::prior("normal", list(0, 1))),
    prior_random = BayesTools::prior_random(sd = prior)
  )
  sd_name <- result[["formula_design"]][["random_effects"]][[1L]][["sd_parameter_names"]]
  samples <- cbind(mu_intercept = seq(-.2, .2, length.out = 20L), sd = .2)
  colnames(samples)[2L] <- sd_name
  fit <- coda::mcmc.list(coda::mcmc(samples))
  attr(fit, "prior_list") <- result[["prior_list"]]
  attr(fit, "formula_design") <- list(mu = result[["formula_design"]])

  structure(list(
    fit  = as_bayestools_fit(fit),
    data = structure(list(), random = TRUE)
  ), class = c("brma.mv", "brma"))
}

# A fitted-object stand-in with a nested variance allocation: a gamma SD
# allocated over the blocks 'study' and 'drug', and the study share split
# over the intercept and slope SD components of 'study'. The study
# components multiply two Dirichlet shares: BayesTools gives them a
# plotting density whose ordinates are not structurally exact. All draws are
# 0.5.
nested_allocation_random_object <- function() {

  data <- data.frame(
    study = factor(c("s1", "s1", "s2", "s2")),
    drug  = factor(c("a", "b", "a", "b")),
    x     = c(-1, 0, 1, 2)
  )
  result <- BayesTools::JAGS_formula(
    formula = ~ 1 + random(1 + x | study, name = "study", covariance = "diag") +
      random(1 | drug, name = "drug", covariance = "diag"),
    parameter = "mu",
    data = data,
    prior_list = list(intercept = BayesTools::prior("normal", list(0, 1))),
    prior_random = BayesTools::prior_random(
      root = BayesTools::random_variance_allocation(
        name = "root", terms = c(study = "study", drug = "drug"),
        sd = BayesTools::prior("gamma", list(2, 2))
      ),
      split = BayesTools::random_variance_allocation(
        name = "split", terms = "study", target = "sd_component",
        parent = BayesTools::allocation_ref("root", "study")
      )
    )
  )
  prior_list <- result[["prior_list"]]
  columns <- unique(unlist(lapply(names(prior_list), function(parameter) {
    BayesTools:::.prior_linear_prior_columns(parameter, prior_list[[parameter]])
  })))
  samples <- matrix(.5, 4L, length(columns), dimnames = list(NULL, columns))
  # A runjags-like fit, so that prior draws can replace its draws.
  fit <- structure(
    list(mcmc = coda::mcmc.list(coda::mcmc(samples)), sample = 4L),
    class = c("runjags", "list")
  )
  attr(fit, "prior_list") <- prior_list
  attr(fit, "formula_design") <- list(mu = result[["formula_design"]])

  structure(
    list(fit = as_bayestools_fit(fit), data = structure(list(), random = TRUE)),
    class = c("brma.mv", "brma")
  )
}

# Independent Savage-Dickey ingredients of the original-scale heterogeneity
# intercept of a log-intercept scale regression on one standardized
# continuous predictor 'slope': tau_0 = b0 exp(-c b1), with the fitted
# intercept b0 ~ Normal(m0, s0) truncated to [0, Inf), the fitted slope
# b1 ~ Normal(m1, s1), and c = mean / sd of the predictor. The prior ordinate
# at 'value' is the scale-product integral over the slope,
#   f(y) = int f_b0(y exp(c t)) exp(c t) phi(t; m1, s1) dt,
# by numerical integration over m1 +- 40 s1 in log space (relative tolerance
# 1e-12); the posterior ordinate is the Gaussian kernel sum of the draws of
# tau_0 (from the fitted MCMC columns) reflected at the support bound 0, with
# bandwidth bw.nrd0.
exp_affine_scale_intercept_reference <- function(fit, slope, value) {

  priors  <- fit[["priors"]][["scale"]]
  scaling <- attr(fit[["fit"]], "formula_scale")[["log_tau"]][[paste0("log_tau_", slope)]]
  shift   <- scaling[["mean"]] / scaling[["sd"]]
  b0 <- priors[["intercept"]]
  b1 <- priors[[slope]]
  stopifnot(
    identical(b0[["distribution"]], "normal"),
    identical(b0[["truncation"]][["lower"]], 0),
    identical(b0[["truncation"]][["upper"]], Inf),
    identical(b1[["distribution"]], "normal")
  )
  log_f_b0 <- function(u) {
    stats::dnorm(u, b0[["parameters"]][["mean"]], b0[["parameters"]][["sd"]], log = TRUE) -
      stats::pnorm(0, b0[["parameters"]][["mean"]], b0[["parameters"]][["sd"]],
                   lower.tail = FALSE, log.p = TRUE)
  }
  limits <- b1[["parameters"]][["mean"]] + c(-40, 40) * b1[["parameters"]][["sd"]]
  prior  <- stats::integrate(function(t) {
    exp(log_f_b0(value * exp(shift * t)) + shift * t +
          stats::dnorm(t, b1[["parameters"]][["mean"]], b1[["parameters"]][["sd"]], log = TRUE))
  }, limits[[1L]], limits[[2L]], rel.tol = 1e-12, subdivisions = 1000L)[["value"]]
  mcmc  <- as.matrix(fit[["fit"]][["mcmc"]])
  draws <- unname(mcmc[, "log_tau_intercept"] * exp(-shift * mcmc[, paste0("log_tau_", slope)]))
  bandwidth <- stats::bw.nrd0(draws)

  list(
    shift     = shift,
    draws     = draws,
    prior     = prior,
    posterior = mean(stats::dnorm(value, draws, bandwidth) +
                       stats::dnorm(value, -draws, bandwidth))
  )
}

# Every alias hypothesis_quantities() lists names its row's quantity: the
# plan of '<alias> = 0' with the row's component, which hypothesis()
# executes, is the plan of that quantity. The alias of a factor term (a row
# with a 'bracket') also names each of its levels: the plan of
# '<alias>[<level>] = 0' with the row's component resolves the level of the
# row's quantity.
.expect_aliases_resolve <- function(object, quantities, info,
                                    metadata = .brma_parameter_catalog_metadata(object),
                                    cache = .hypothesis_plan_cache(object)) {

  plan_of <- function(reference, component) {
    tryCatch(
      .hypothesis_plans(
        object     = object,
        hypothesis = paste0("`", reference, "` = 0"),
        component  = component,
        metadata   = metadata,
        cache      = cache
      )[[1L]],
      error = function(error) error
    )
  }
  expect_plan <- function(plan, row, case) {
    expect_false(inherits(plan, "error"), info = paste(
      case, if (inherits(plan, "error")) conditionMessage(plan)
    ))
    if (!inherits(plan, "error")) {
      expect_identical(
        c(plan[["parameter"]], plan[["component"]]),
        c(row[["parameter"]], row[["component"]]),
        info = case
      )
    }
    !inherits(plan, "error")
  }
  entries <- metadata[["entries"]]

  for (i in seq_len(nrow(quantities))) {
    row  <- quantities[i, , drop = FALSE]
    case <- paste(info, row[["component"]], row[["alias"]])
    expect_plan(plan_of(row[["alias"]], row[["component"]]), row, case)
    if (is.na(row[["bracket"]])) {
      next
    }
    entry <- entries[entries[["parameter"]] == row[["parameter"]] &
                       entries[["component"]] == row[["component"]], ,
                     drop = FALSE]
    levels <- .hypothesis_plan_term_levels(
      metadata,
      as.list(entry[1L, setdiff(names(entry), "aliases"), drop = FALSE])
    )
    expect_gt(nrow(levels), 0L)
    for (j in seq_len(nrow(levels))) {
      level      <- levels[["level"]][[j]]
      level_case <- paste0(case, "[", level, "]")
      plan <- plan_of(paste0(row[["alias"]], "[", level, "]"), row[["component"]])
      if (expect_plan(plan, row, level_case)) {
        expect_identical(
          unique(plan[["selected"]][["resolution"]][["occurrences"]][["quantity_id"]]),
          levels[["quantity_id"]][[j]],
          info = level_case
        )
      }
    }
  }

  invisible(quantities)
}


# hypothesis() evaluates a statement exactly when its plan admits the method,
# and otherwise stops with the plan's refusal (its first class and message).
# hypothesis_quantities() renders the same plans, and lists aliases that name
# their rows' quantities. qCMDE/IWMDE statements the plans admit are
# evaluated for the first point and contrast statement of each object
# ('run_precomputed'), with a small density budget.
.expect_plans_consistent <- function(object, info, run_precomputed = TRUE) {

  metadata   <- .brma_parameter_catalog_metadata(object)
  quantities <- hypothesis_quantities(object)
  cache      <- .hypothesis_plan_cache(object)
  .expect_aliases_resolve(object, quantities, info, metadata, cache)
  control    <- list(n_points = 20, samples = 50)
  run <- function(statement, component, method) {
    tryCatch(
      suppressWarnings(hypothesis(
        object, statement, component = component, density_method = method,
        density_control = if (method %in% c("qCMDE", "IWMDE")) control,
        seed = 1, n_samples = 500
      )),
      error = function(error) error
    )
  }
  ran_precomputed <- character()
  entries <- metadata[["entries"]]
  entries <- entries[entries[["component"]] != "bias", , drop = FALSE]
  for (i in seq_len(nrow(entries))) {
    entry <- as.list(entries[i, setdiff(names(entries), "aliases"), drop = FALSE])
    rows  <- quantities[quantities[["parameter"]] == entry[["parameter"]] &
                          quantities[["component"]] == entry[["component"]], ,
                        drop = FALSE]
    plans <- .hypothesis_quantities_plans(object, entry, metadata, cache)
    rendered <- .hypothesis_quantities_render_plans(
      plans   = plans,
      bracket = identical(entry[["role"]], "formula_coefficient_group")
    )
    row_info <- paste(info, entry[["parameter"]])
    expect_gt(nrow(rows), 0L)
    for (column in c("point_test", "direction_test", "contrast_test",
                     "point_test_methods", "contrast_test_methods", "reason")) {
      expect_identical(
        unique(rows[[column]]), rendered[[column]],
        info = paste(row_info, column)
      )
    }
    for (type in c("point", "region", "contrast")) {
      for (plan in plans[[type]]) {
        statement <- BayesTools::hypothesis_render(plan[["statement"]])
        for (method in .hypothesis_plan_advertised_methods()) {
          refusal <- .hypothesis_plan_status(plan, method)
          case    <- paste(row_info, statement, method)
          precomputed <- method %in% c("qCMDE", "IWMDE")
          if (is.null(refusal) && precomputed &&
              (!run_precomputed || paste(type, method) %in% ran_precomputed ||
                 identical(type, "region"))) {
            next
          }
          out <- run(statement, entry[["component"]], method)
          if (is.null(refusal)) {
            # An admitted statement runs; qCMDE/IWMDE ordinates may still be
            # rejected by their numerical diagnostics.
            expect_false(
              inherits(out, "error") &&
                !(precomputed && inherits(out, "RoBMA_density_ordinate_error")),
              info = paste(case, if (inherits(out, "error")) conditionMessage(out))
            )
            if (precomputed) {
              ran_precomputed <- c(ran_precomputed, paste(type, method))
            }
          } else {
            expect_s3_class(out, refusal[["class"]][[1L]])
            expect_identical(conditionMessage(out), refusal[["reason"]], info = case)
          }
        }
      }
    }
  }

  invisible(quantities)
}
