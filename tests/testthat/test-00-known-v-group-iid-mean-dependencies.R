test_that("group-IID q shortcuts require declared fixed-mean invariance", {

  r <- function(name) get(name, envir = asNamespace("RoBMA"), inherits = FALSE)
  helpers <- list(review_compile_object = review_compile_object,
    review_prior_columns = review_prior_columns,
    review_coefficient_name = review_coefficient_name,
    review_attach_replay = review_attach_replay)
  make_fixture <- function(dependent) {
    dat <- data.frame(yi = rep(0, 4), x = c(1, 2, 1, 2), study = factor(c('a', 'a', 'b', 'b')))
    gamma <- BayesTools::prior('gamma', list(shape = 2, rate = 1))
    constructor <- function(x_prior) brma.mv(yi = yi, V = diag(4), mods = ~0 + x,
      random = ~id(1 | study), data = dat, measure = 'GEN', prior_unit_information_sd = 1,
      prior_mods = list(x = x_prior), prior_effect = BayesTools::prior('point', list(location = 0)),
      prior_heterogeneity = BayesTools::prior_random(sd = gamma),
      standardize_continuous_predictors = FALSE, only_priors = TRUE, silent = TRUE)
    discovery <- helpers$review_compile_object(constructor(gamma))
    source <- discovery$compiled$mu$formula_design$random_effects[[1L]]$sd_parameter_names[[1L]]
    x_prior <- gamma
    if (dependent) attr(x_prior, 'multiply_by') <- source
    setup <- helpers$review_compile_object(constructor(x_prior))
    priors <- c(setup$compiled$mu$prior_list, r('.create_fit_priors')(setup$object$data, setup$object$priors))
    term <- setup$compiled$mu$formula_design$random_effects[[1L]]
    name_map <- setup$compiled$mu$formula_design$name_map
    latent_base <- name_map$jags_name[name_map$kind == 'random' & name_map$role == '_xRE_Zx']
    stopifnot(length(latent_base) == 1L, identical(setup$compiled$mu$add_parameters, latent_base),
      term$n_groups == 2L, term$n_columns == 1L)
    latent_columns <- unlist(lapply(seq_len(term$n_columns), function(j) {
      paste0(latent_base, '[', seq_len(term$n_groups), ',', j, ']')
    }), use.names = FALSE)
    columns <- unique(c(helpers$review_prior_columns(priors), latent_columns))
    samples <- matrix(.6, 24L, length(columns), dimnames = list(NULL, columns))
    x_name <- helpers$review_coefficient_name(setup, 'x')
    intercept <- helpers$review_coefficient_name(setup, 'intercept')
    samples[, intercept] <- 0
    samples[, x_name] <- c(1, seq(.65, 1.75, length.out = 23))
    samples[, source] <- c(1, seq(.7, 1.8, length.out = 23))
    stopifnot(all(samples[, latent_columns, drop = FALSE] == .6))
    required_monitors <- unique(c(BayesTools::JAGS_to_monitor(priors), setup$compiled$mu$add_parameters))
    supplied_monitor_bases <- unique(sub('\\[.*$', '', colnames(samples)))
    stopifnot(all(required_monitors %in% supplied_monitor_bases), length(latent_columns) == 2L)
    object <- helpers$review_attach_replay(setup, samples)
    BayesTools::JAGS_validate_fit_contract(object$fit, requires = c('name_encoding', 'formula_name_map',
      'formula_design', 'parameter_map', 'draw_geometry'))
    stopifnot(identical(class(object), class(setup$object)), all(is.finite(samples)))
    quantities <- BayesTools::parameter_catalog(object$fit)$quantities
    selected <- quantities[quantities$role == 'random_sd', , drop = FALSE]
    selected <- selected[vapply(selected$extraction_key, function(key) identical(key$source_parameter, source), logical(1)), , drop = FALSE]
    stopifnot(nrow(selected) == 1L)
    selector <- r('.brma_random_parameter_io_labels')(selected, 'selector')
    target <- r('.brma_random_parameter_density_target')(object, selector)
    context <- r('.iwmde_context')(object)
    spec <- r('.iwmde_parameter_spec')(context, target$parameter, target$parameter_spec)
    states <- r('.iwmde_row_states')(context, 1L, parameter = target$parameter, parameter_spec = spec)
    replacement <- r('.iwmde_replacement_spec')(context, target$parameter, spec)
    list(dependent = dependent, dat = dat, setup = setup, source = source, x_name = x_name,
      intercept = intercept, declared_priors = priors, samples = samples, object = object,
      quantities = quantities, selected = selected, selector = selector, target = target, context = context,
      spec = spec, states = states, replacement = replacement, latent_base = latent_base,
      latent_columns = latent_columns, required_monitors = required_monitors)
  }

  for (dependent in c(TRUE, FALSE)) {
    fixture <- make_fixture(dependent)
    args <- list(context = fixture$context, parameter = fixture$target$parameter,
      values = c(1, 2), row_states = fixture$states,
      replacement = fixture$replacement)
    active <- fixture$states[[1L]]$active_setup
    plan <- do.call(r(".iwmde_known_v_random_group_iid_plan"), list(
      context = fixture$context, parameter = args$parameter,
      replacement = fixture$replacement, priors = active$priors))
    candidates <- do.call(r(".iwmde_build_replacement_samples"), args)
    expect_identical(plan$target_mode, "direct_sd")
    expect_identical(r(".iwmde_q_grid_fixed_mu_is_invariant")(
      fixture$context, candidates), !dependent)
    expect_identical(fixture$states[[1L]]$state_scope, "global")
    expect_identical(fixture$states[[1L]]$likelihood_mode, "marginal")
    expect_setequal(names(fixture$states[[1L]]$prior_list),
      c(fixture$source, fixture$x_name, fixture$intercept))

    # Local normal priors are integrated by the global marginal target.
    # The complete oracle includes both Gamma factors and structural point0.
    Z <- cbind(a = as.integer(fixture$dat$study == "a"),
      b = as.integer(fixture$dat$study == "b"))
    expected <- vapply(args$values, function(v) {
      covariance <- diag(4) + v^2 * tcrossprod(Z)
      mu <- fixture$dat$x * if (dependent) v else 1
      factor <- chol(covariance)
      residual <- backsolve(factor, fixture$dat$yi - mu, transpose = TRUE)
      log_likelihood <- -.5 * (4 * log(2 * pi) +
        2 * sum(log(diag(factor))) + sum(residual^2))
      log_likelihood + dgamma(v, 2, 1, log = TRUE) +
        dgamma(1, 2, 1, log = TRUE) + 0
    }, numeric(1))
    generic <- do.call(r(".iwmde_log_q_grid_batch"), args)
    scalar <- do.call(r(".iwmde_log_q_grid_scalar"), args)
    dispatcher <- do.call(r(".iwmde_log_q_grid"), args)
    direct <- do.call(r(".iwmde_log_q_grid_known_v_random_group_iid"),
      c(args, list(active_setup = active)))
    expect_equal(as.numeric(generic), expected, tolerance = 1e-12)
    expect_equal(as.numeric(scalar), expected, tolerance = 1e-12)
    expect_equal(as.numeric(dispatcher), expected, tolerance = 1e-12)
    if (dependent) {
      expect_null(direct)
      expect_equal(expected, c(-8.7743664214868, -12.179831529595),
        tolerance = 1e-12)
    } else {
      expect_true(is.matrix(direct))
      expect_equal(as.numeric(direct), expected, tolerance = 1e-12)
    }
  }
})
