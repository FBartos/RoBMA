CERTIFICATION_CASE_TIMEOUT_SECONDS <- 60 * 60


.required_tests <- function(file, test) {

  data.frame(
    file             = rep(file, length(test)),
    test             = test,
    stringsAsFactors = FALSE
  )
}


.iwmde_fast_required_tests <- function() {

  .required_tests(
    "test-02-iwmde-fast-paths.R",
    c(
      "IWMDE batched q evaluation matches scalar fallback",
      "IWMDE predictor fast path matches scalar fallback",
      "IWMDE normal quadratic fast path matches scalar fallback"
    )
  )
}


.v14_metafor_required_tests <- function(include_random = TRUE,
                                         include_new_effect = TRUE) {

  tests <- c(
    "v14 brma.mv fixed effects match metafor references",
    "v14 brma.mv heterogeneity components match metafor references",
    "v14 brma.mv heterogeneity component selectors expose expected names",
    "v14 brma.mv fixed fitted values and residuals match metafor references",
    "v14 brma.mv ranef components track metafor references",
    "v14 brma.mv ranef, blup, and true_effects use consistent targets",
    "v14 brma.mv hatvalues track metafor leverages"
  )
  if (include_random) {
    tests <- c(
      tests,
      "v14 brma.mv random-covariance parameters match metafor references"
    )
  }
  if (include_new_effect) {
    tests <- c(
      tests,
      "explicit brma.mv new-effect predictions track metafor targets"
    )
  }

  return(.required_tests("test-02-brma-mv-metafor.R", tests))
}


.random_parameter_required_tests <- function(include_overlay = FALSE) {

  tests <- c(
    "random catalog matches summaries across structures",
    "random plots and MCMC diagnostics use semantic draws",
    "random directional hypotheses use induced joint priors",
    "random point hypotheses follow quantity-specific policy"
  )
  if (include_overlay) {
    tests <- c(
      tests,
      "random prior overlays and diagnostic labels are semantic"
    )
  }

  return(.required_tests("test-02-random-parameters.R", tests))
}


required_native_test_symbols <- function() {

  c(
    paste0("RoBMA_selnorm_kernel_", c(
      "loglik_matrix", "log_norm_matrix", "cdf_matrix", "moments_matrix",
      "rng_matrix", "weighted_summary", "log_norm_delta_grid", "loglik_row_sum"
    )),
    paste0("RoBMA_known_v_covariance_plan_", c(
      "create", "loglik_batch", "conditional_loglik_batch",
      "conditional_summary_batch", "affine_grid_loglik", "factor_grid_loglik"
    )),
    "RoBMA_norm_predictor_grid_loglik",
    "RoBMA_norm_cluster_analytic_loglik",
    "RoBMA_selnorm_sampling_conditioned_batch",
    "RoBMA_selnorm_sampling_deletion_loglik_batch",
    "RoBMA_selnorm_gaussian_event_mass_batch",
    "RoBMA_glmm_binom_aghq", "RoBMA_glmm_pois_aghq",
    "RoBMA_plot_normal_mixture_quantiles",
    "RoBMA_plot_selnorm_mixture_quantiles",
    "RoBMA_selnorm_zcurve_threshold_summary",
    "RoBMA_selnorm_zcurve_density_matrix"
  )
}


validate_native_test_symbols <- function(
    registered = names(getDLLRegisteredRoutines("RoBMA")[[".Call"]]),
    registration_source = NULL) {

  required <- required_native_test_symbols()
  if (!is.null(registration_source)) {
    declarations <- readLines(registration_source, warn = FALSE)
    pattern <- '^\\s*\\{\\s*"(RoBMA_[^"]+)"\\s*,\\s*\\(DL_FUNC\\)'
    matched <- regmatches(declarations, regexec(pattern, declarations, perl = TRUE))
    declared <- vapply(Filter(function(x) length(x) > 1L, matched), `[[`, character(1), 2L)
    if (length(declared) == 0L) {
      stop("No native entry-point declarations found in the test source tree.", call. = FALSE)
    }
    required <- union(required, declared)
  }
  missing <- setdiff(required, registered)
  if (length(missing) > 0L) {
    stop(
      "Required native test entry points are unavailable: ",
      paste(missing, collapse = ", "),
      ". Rebuild and reload the tested RoBMA DLL before running this profile.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}


.hypothesis_plan_required_tests <- function() {

  .required_tests("test-00-hypothesis-plan.R", c(
    "hypothesis_quantities() renders the plans that hypothesis() executes",
    "hypothesis_quantities() names quantities with shared aliases by their selector",
    "qCMDE/IWMDE point hypotheses are refused for several scale formulas",
    "factor levels of every contrast have point tests and level contrasts",
    "levels of a factor alias shared by the location and scale formulas resolve with 'component'",
    "statements that do not match 'component' are component mismatches at every entry point",
    "BayesTools resolution errors that escape planning are statement errors",
    "publication-bias parameters are refused as hypothesis targets",
    "catalog quantities without a RoBMA parameter are refused as targets, not as stale fits",
    "the selection of plot() refuses catalog quantities without a RoBMA parameter as hypothesis() does",
    "display labels of parameters are evaluable with and without 'component'",
    "statement refusals of the plans are statement errors, not unavailable tests",
    "unresolved references have the same classes on fits and marginal means",
    "statements on every level of a term with a fixed level are refused as fixed",
    "mean-difference level point hypotheses follow the exact Savage-Dickey computation",
    "levels of mean-difference multivariate t factors have exact point tests",
    "a two-level ordered factor tests its level and the level contrast alike",
    "ordered levels with a Beta share have an infinite prior ordinate at zero",
    "a random-effect variance and its standard deviation give the same Bayes factor",
    "point hypotheses with an inexact prior ordinate are refused with its class",
    "region prior probabilities come from prior draws only for densities without provenance",
    "level contrasts of model-averaged factors condition on the included term",
    "levels of model-averaged mean-difference factors have exact point tests"
  ))
}


standard_required_tests <- function() {

  rbind(
    .hypothesis_plan_required_tests(),
    .required_tests("test-00-carrier-integration.R", c(
      "public no-intercept moderator means use the fitted contrast basis",
      "whole marginal-means regions retain structural fixed-cell refusal",
      "expanded public factor groups retain every statement-major result row",
      "point and region marginal groups select their own execution method",
      "derived public OR targets carry their effective declaration to scalar attachment",
      "ordered global totals participate in marginal inclusion conditioning",
      "declared stochastic contribution refusals are planned for levels and combinations"
    )),
    .required_tests("test-00-formula-contribution-state.R", c(
      "public contribution combinations fold constants once for both ordinate methods",
      "actual brma coefficients and static-method routes use the compiled owner",
      "scalar discovery retains declared refusal rows and usable sibling quantities"
    )),
    .required_tests("test-00-exp-affine-supported-routes.R", c(
      "exp-affine coefficients refuse qCMDE/IWMDE with the original-scale classes"
    )),
    .required_tests("test-00-factor-level-targets.R", c(
      "fixed reference-level statements keep one reason under qCMDE for model-averaged fits",
      "factor-level point hypotheses use the level's own fitted coordinate",
      "mean-difference level point hypotheses are linear targets of the contrast coordinates",
      "contrast-coefficient selectors stop naming the level-label form",
      "factor-cell plots draw the selected level for every contrast",
      "repeated hypothesis rows are numbered by statement",
      "hypothesis_quantities reports point and contrast tests for the levels of every contrast",
      "linear-target refusals of BayesTools keep their classes",
      "level refusals name the level by its selector",
      "point hypotheses on the fixed reference level give one reason for every method"
    )),
    .required_tests("test-00-hypothesis-scalar-eligibility.R", c(
      "scalar point-test eligibility matches the hypotheses that run"
    )),
    .required_tests("test-00-input-data-mv-04-likelihood.R", c(
      "fitted known-V row SD sources read the coordinates of factor scale terms"
    )),
    .required_tests("test-00-random-correlation-routes.R", c(
      "correlation hypotheses use the correlation draws",
      "correlation hypotheses use the defined draws of a declared quantity",
      "summary footnotes follow the renamed random-effect rows",
      "correlation plots draw the exact LKJ prior on the fitted scale",
      "correlation plots without an exact prior density draw the posterior alone",
      "qCMDE/IWMDE stops for the correlation name the requested operation",
      "original-scale correlations of allocated blocks are atom-free and plot"
    )),
    .required_tests("test-00-summary_heterogeneity.R", c(
      "undefined original-scale correlations are summarized over defined draws"
    )),
    .required_tests("test-00-random-nested-total-prior.R", c(
      "the total of a nested gated allocation has its exact prior curve and atom in plots",
      "point hypotheses on the total of a nested gated allocation use its exact prior"
    )),
    .required_tests("test-00-plot-prior-random.R", c(
      "plot_prior and print_prior show the declared prior of a random SD",
      "plot_prior draws the prior of a variance as the square of the SD prior",
      "random selectors of gated allocations show prior laws with their point masses",
      "random selectors of a nested gated allocation derive from its allocation chain",
      "random correlations without an original-scale prior retain target refusals",
      "only_priors objects point to the random-effect specification"
    )),
    .required_tests("test-02-random-prior-plots.R", c(
      "random SD, total and share selectors of an allocation show their priors",
      "the SD of a scale-modelled block keeps resolving to the scale intercept"
    )),
    .required_tests("test-02-bias-branch-design.R", c(
      "product-space VIF uses each member's design",
      "product-space heterogeneity summaries use each member's design",
      "hat-based diagnostics still refuse ensembles with weight-function members",
      "known-V product-space hat values run on the design of each draw's member",
      "known-V product-space VIF averages the draws of each design class",
      "the ensemble has draws of every design class",
      "hat values average the hat matrix of each draw's member design",
      "a fitted single-bias model is one design",
      "a single bias term keeps its own hat values",
      "internally standardized residuals use each draw's member design",
      "the hat column of influence() uses each draw's member design",
      "VIF reports each member design's GVIF",
      "VIF of a bias term without draws is unavailable",
      "heterogeneity summaries use each member design's typical variance",
      "every per-draw output equals the single-design output on its class rows",
      "known-V hat values use each draw's member design",
      "known-V hat values of thinned draws use the designs of the selected draws",
      "known-V design classes follow the rows of the thinned posterior",
      "known-V standardized residuals use each draw's member design",
      "known-V VIF reports each member design's GVIF",
      "every known-V per-draw output equals the single-design output on its class rows",
      "intercept-only ensembles average the member designs' hat values"
    )),
    .required_tests("test-00-selection-kernel-fit.R", c(
      "selection reference weights are structural convergence parameters"
    )),
    .required_tests("test-02-glmm-aghq.R", c(
      "a fitted point-baserate GLMM supports log-likelihood criteria",
      "boundary-concentrated beta priors support likelihood criteria",
      "fitted truncated-nuisance GLMMs support likelihood criteria"
    )),
    .required_tests("test-00-marginal-means-inference.R", c(
      "marginal_means() reports rows with a zero prior ordinate at the null as NA",
      "declared unavailable marginal measures retain valid sibling BFs and late KDE requests"
    )),
    .required_tests("test-00-selection-kernel-step.R", c(
      "step selected-normal kernel matches an independent p-bin reference",
      "step normalizers retain central and extreme one-sided tails",
      "step selected-normal RNG follows exact bin masses and moments"
    )),
    .required_tests("test-00-selection-kernel-native.R",
      "native funnel model-averaged quantiles match R fallback"),
    .required_tests("test-00-known-v-joint-loglik.R",
      "native covariance plan returns exact Schur conditional densities"),
    .required_tests("test-02-iwmde-oracles.R",
      "qCMDE matches the conditional-normal random-effects oracle"),
    .required_tests("test-02-pr82-estimator-integration.R", c(
      "public density hypotheses and plots match conjugate Gaussian inference",
      "marginal means retain genuine qCMDE curves and usable point ordinates"
    ))
  )
}


validate_test_results <- function(results) {

  results <- as.data.frame(results)
  required_columns <- c("failed", "error", "warning")
  missing_columns  <- setdiff(required_columns, names(results))
  if (length(missing_columns) > 0L) {
    stop(
      "Test results are missing required columns: ",
      paste(missing_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  problems <- vapply(required_columns, function(name) {
    sum(as.integer(results[[name]]), na.rm = TRUE)
  }, integer(1))
  if (any(problems > 0L)) {
    stop(
      "Test execution recorded problems: failed=", problems[["failed"]],
      ", errors=", problems[["error"]],
      ", warnings=", problems[["warning"]], ".",
      call. = FALSE
    )
  }

  return(invisible(TRUE))
}


validate_certification_evidence <- function(results, required_tests,
                                            case_name) {

  .validate_required_test_evidence(
    results, required_tests, paste0("Certification case '", case_name, "'")
  )
}


validate_standard_evidence <- function(results) {

  .validate_required_test_evidence(
    results, standard_required_tests(), "Standard profile"
  )
}


.validate_required_test_evidence <- function(results, required_tests, evidence_label) {

  if (is.null(required_tests) || nrow(required_tests) == 0L) {
    return(invisible(TRUE))
  }

  results <- as.data.frame(results)
  required_columns <- c("file", "test", "skipped", "passed")
  missing_columns  <- setdiff(required_columns, names(results))
  if (length(missing_columns) > 0L) {
    stop(
      evidence_label, " results are missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  for (i in seq_len(nrow(required_tests))) {
    file <- required_tests[["file"]][[i]]
    test <- required_tests[["test"]][[i]]
    matches <- which(results[["file"]] == file & results[["test"]] == test)
    label <- paste0(file, " :: ", test)

    if (length(matches) != 1L) {
      stop(
        evidence_label, " requires exactly one result for ",
        label, ", but observed ", length(matches), ".",
        call. = FALSE
      )
    }
    if (isTRUE(results[["skipped"]][[matches]])) {
      stop(
        evidence_label, " skipped required test ", label,
        ".",
        call. = FALSE
      )
    }
    if (is.na(results[["passed"]][[matches]]) ||
        results[["passed"]][[matches]] < 1L) {
      stop(
        evidence_label,
        " recorded no passing expectations for required test ", label, ".",
        call. = FALSE
      )
    }
  }

  return(invisible(TRUE))
}


certification_cases <- function() {

  # testthat matches a verification filter against the file stems; each
  # filter is anchored so that it selects whole stems only: '03-loo' does not
  # also select '03-loo-exact-refits', nor 'predict' '02-predict-mv'. File
  # families are selected by explicit wildcards ('plot.*').
  anchored <- function(...) paste0("^(", paste0(...), ")$")
  multivariate_tests <- paste0(
    "02-(brma-mv.*|forest|funnel|",
    "heterogeneity-mv|hypothesis|influence|marginal_means|plot.*|",
    "predict-mv|qqnorm|random-parameters|regplot|residuals|summary.*|vif)|",
    "03-(bridgesampling|loo|zplot)"
  )
  multivariate_filter <- anchored(multivariate_tests)

  list(
    "numerical-kernels" = list(
      description = paste(
        "Independent selected-normal, known-V, distribution, and GLMM",
        "quadrature oracles."
      ),
      fit_sources = "test-01-ordinary-lane-fixtures.R",
      fit_names = c("lane_glmm_point_baserate", "lane_glmm_boundary_beta", "lane_glmm_truncated_beta", "lane_glmm_truncated_lograte", "lane_selection_reference"),
      test_filter = anchored(
        "00-(covariance-factorization|known-v-joint-loglik|",
        "selection-kernel.*|selection-probability-numerics|",
        "selection-factor-routing|selection-estimate-cells|",
        "selection-coarse-exactness|jags-build-requirements)|",
        "02-(distributions|glmm-aghq)"
      ),
      required_tests = rbind(
        .required_tests("test-00-covariance-factorization.R",
          "covariance sampling factor preserves singular covariance"),
        .required_tests("test-00-selection-kernel-fit.R",
          "selection reference weights are structural convergence parameters"),
        .required_tests("test-00-selection-coarse-exactness.R",
          "the corrected chain on a rank-two factor block targets the exact posterior"),
        .required_tests("test-02-glmm-aghq.R", c(
          "a fitted point-baserate GLMM supports log-likelihood criteria",
          "boundary-concentrated beta priors support likelihood criteria",
          "fitted truncated-nuisance GLMMs support likelihood criteria"
        ))
      )
    ),
    "normal-models" = list(
      description = paste(
        "Normal, PET-PEESE, selection, and model-averaged metafor parity",
        "and visual diagnostics."
      ),
      fit_sources = c(
        "test-01-brma.norm.R",
        "test-01-bPET.R",
        "test-01-bPEESE.R",
        "test-01-bselmodel.R",
        "test-01-BMA.norm.R",
        "test-01-RoBMA.R",
        "test-01-vif-parity.R",
        "test-01-iwmde-multilevel-selection.R"
      ),
      test_filter = anchored(
        "02-(dfbetas|forest|funnel|hatvalues|influence|iwmde-fast-paths|",
        "iwmde-normal-predictor-grid|marginal_means|",
        "plot.*|predict|qqnorm|radial|regplot|residuals|summary.*|vif)|",
        "03-(bridgesampling|loo|zplot)"
      ),
      required_tests = rbind(
        .required_tests(
          "test-02-iwmde-fast-paths.R",
          c(
            "multilevel weightfunction formula path matches scalar fallback",
            "IWMDE dispatcher evaluates the joint density of singleton-cluster selection fits"
          )
        ),
        .iwmde_fast_required_tests()
      )
    ),
    "glmm-models" = list(
      description = paste(
        "Binomial and Poisson GLMM fitting, quadrature, diagnostics, and",
        "model averaging."
      ),
      fit_sources = c(
        "test-01-brma.glmm.R",
        "test-01-BMA.glmm.R"
      ),
      test_filter = anchored(
        "02-(dfbetas|distributions|forest|funnel|glmm-aghq|hatvalues|",
        "influence|iwmde-fast-paths|iwmde-glmm-local|marginal_means|predict|",
        "qqnorm|residuals|summary.*|vif)|",
        "03-(bridgesampling|loo)"
      ),
      required_tests = rbind(
        .required_tests(
          "test-02-iwmde-fast-paths.R",
          "IWMDE batched q evaluation matches scalar fallback"
        ),
        .required_tests(
          "test-02-iwmde-glmm-local.R",
          "IWMDE conditions multilevel GLMM rows on sampled local states"
        )
      )
    ),
    "multivariate-core" = list(
      description = paste(
        "Core known-V multivariate parameterizations, prediction,",
        "diagnostics, LOO, and marginal likelihoods."
      ),
      fit_sources = "test-01-brma.mv.R",
      fit_names = c(
        "brma.mv_latent",
        "brma.mv_whitened",
        "brma.mv_block_mvn",
        "brma.mv_block_mvn_fixed_random_null",
        "brma.mv_block_mvn_random"
      ),
      test_filter = multivariate_filter,
      required_tests = .required_tests(
        "test-03-loo.R",
        "brma.mv known-V fits expose conditional estimate-unit LOO and WAIC"
      )
    ),
    "loo-exact-refits" = list(
      description = paste(
        "PSIS estimate-unit LOO compared with exact observation-deletion",
        "refits of Kearon US/HCS and Ishak HAR scenario models."
      ),
      fit_sources = character(),
      test_filter = anchored("03-loo-exact-refits"),
      required_tests = .required_tests(
        "test-03-loo-exact-refits.R",
        "scenario-model PSIS LOO agrees with five exact deletion refits"
      )
    ),
    "selection-recovery-marg" = list(
      description = paste(
        "Independent marginal-selection simulations with diagonal, structured,",
        "and dense V; fixed/nested parameter recovery and a dense posterior oracle."
      ),
      fit_sources = character(),
      test_filter = anchored("03-selection-recovery-marg"),
      required_tests = .required_tests(
        "test-03-selection-recovery-marg.R",
        c(
          "marginal selection recovers simulated fixed and nested parameters",
          "marginal selection JAGS posterior matches independent dense integration"
        )
      )
    ),
    "selection-recovery-cond" = list(
      description = paste(
        "Independent conditional-selection simulations with diagonal, structured,",
        "and dense V; fixed/nested parameter recovery and a dense posterior oracle."
      ),
      fit_sources = character(),
      test_filter = anchored("03-selection-recovery-cond"),
      required_tests = .required_tests(
        "test-03-selection-recovery-cond.R",
        c(
          "conditional selection recovers simulated fixed and nested parameters",
          "conditional selection JAGS posterior matches independent dense integration"
        )
      )
    ),
    "multivariate-extended" = list(
      description = paste(
        "Extended known-R, scale, allocation, and moderator multivariate",
        "parameterizations."
      ),
      fit_sources = "test-01-brma.mv.R",
      fit_names = c(
        "brma.mv_latent",
        "brma.mv_whitened",
        "brma.mv_block_mvn",
        "brma.mv_block_mvn_random",
        "brma.mv_block_mvn_random_sampled",
        "brma.mv_block_mvn_known_R",
        "brma.mv_latent_estimate_scale",
        "brma.mv_block_mvn_estimate_scale",
        "brma.mv_block_mvn_random_scale",
        "brma.mv_block_mvn_3lvl_scale_total",
        "brma.mv_block_mvn_3lvl_scale_top",
        "brma.mv_block_mvn_3lvl_scale_bottom",
        "brma.mv_block_mvn_mods",
        "brma.mv_block_mvn_random_mods_scale"
      ),
      test_filter = anchored(multivariate_tests, "|02-iwmde-api"),
      required_tests = rbind(
        .required_tests(
          "test-03-loo.R",
          "brma.mv known-V fits expose conditional estimate-unit LOO and WAIC"
        ),
        .required_tests(
          "test-02-iwmde-api.R",
          "known-V marginalized allocation weights stay in global IWMDE state"
        ),
        .required_tests(
          "test-02-random-parameters.R",
          "IWMDE disables focal prior delta for sampled random SD rows"
        ),
        .random_parameter_required_tests()
      )
    ),
    "multivariate-singular" = list(
      description = paste(
        "Structurally regularized singular known-V multivariate",
        "parameterizations."
      ),
      fit_sources = "test-01-brma.mv.R",
      fit_names = c(
        "brma.mv_singular_regularized_whitened",
        "brma.mv_singular_regularized_block_mvn"
      ),
      test_filter = multivariate_filter,
      required_tests = .required_tests(
        "test-02-summary.R",
        "summary.brma returns a stable object contract"
      )
    ),
    "multivariate-parity-cs" = list(
      description = "Metafor parity for the Konstantopoulos CS model.",
      fit_sources = "test-01-brma.mv.R",
      fit_names = "brma.mv_v14_konstantopoulos2011_cs",
      test_filter = multivariate_filter,
      required_tests = rbind(
        .v14_metafor_required_tests(),
        .random_parameter_required_tests(include_overlay = TRUE)
      )
    ),
    "multivariate-parity-nested" = list(
      description = "Metafor parity for the Assink nested model.",
      fit_sources = "test-01-brma.mv.R",
      fit_names = "brma.mv_v14_assink2016_nested",
      test_filter = multivariate_filter,
      required_tests = .v14_metafor_required_tests(include_random = FALSE)
    ),
    "multivariate-parity-har" = list(
      description = "Metafor parity for the Ishak HAR model.",
      fit_sources = "test-01-brma.mv.R",
      fit_names = "brma.mv_v14_ishak2007_har",
      test_filter = multivariate_filter,
      required_tests = rbind(
        .v14_metafor_required_tests(),
        .random_parameter_required_tests()
      )
    ),
    "multivariate-parity-treatment" = list(
      description = "Metafor parity for the Begg treatment CS model.",
      fit_sources = "test-01-brma.mv.R",
      fit_names = "brma.mv_v14_begg1989_study_treatment",
      test_filter = multivariate_filter,
      required_tests = rbind(
        .v14_metafor_required_tests(include_new_effect = FALSE),
        .random_parameter_required_tests(include_overlay = TRUE)
      )
    ),
    "iwmde-qcmde" = list(
      description = paste(
        "Fitted and analytic qCMDE/IWMDE density, ordinate, and bridge",
        "oracles."
      ),
      fit_sources = c(
        "test-01-hypothesis-plan.R",
        "test-01-brma.glmm.R",
        "test-01-bselmodel.R",
        "test-01-brma.mv.R",
        "test-01-iwmde-oracle-nested.R"
      ),
      fit_names = c(
        "hypothesis_plan_treatment",
        "hypothesis_plan_meandif",
        "hypothesis_plan_orthonormal",
        "hypothesis_plan_independent",
        "hypothesis_plan_ordered_2",
        "hypothesis_plan_ordered_3",
        "hypothesis_plan_ordered_4",
        "hypothesis_plan_scaled_interaction",
        "hypothesis_plan_location_scale",
        "hypothesis_plan_model_averaged",
        "hypothesis_plan_robma_mixture",
        "hypothesis_plan_two_scale",
        "hypothesis_plan_meandif_mt",
        "brma.mv_latent",
        "brma.mv_whitened",
        "brma.mv_block_mvn",
        "brma.mv_block_mvn_random",
        "brma.mv_block_mvn_fixed_random_null",
        "nielweise2008_glmm",
        "nielweise2008_glmm_effect_null",
        "dat.lehmann2018-3PSM",
        "dat.lehmann2018-3PSM_effect_null",
        "iwmde_known_v_tau_full",
        "iwmde_known_v_tau_null"
      ),
      test_filter = anchored(
        "00-(hypothesis-plan|iwmde-ordered-targets)|",
        "02-(hypothesis|iwmde.*|marginal_means|random-parameters)|",
        "03-(bridgesampling|loo)"
      ),
      required_tests = rbind(
        .hypothesis_plan_required_tests(),
        .required_tests(
          "test-00-hypothesis-plan.R",
          "synthetic continuous controls execute every advertised KDE plan"
        ),
        .required_tests(
          "test-00-iwmde-ordered-targets.R",
          c(
            "regular signed fixed-total interior retains positive conditional contributions",
            "admitted nonpure share endpoints use the allocation group source",
            "physical ordered block oracles cover shifted numerators and selected denominators",
            "ordered partial uncertainty survives result assembly and public consumers",
            "near-fit ordered Gaussian and observed step factors match full original replay",
            "nonzero Bernoulli interior retains the regular conditional source",
            "fixed allocations with multiple total slices do not certify constant contributions"
          )
        ),
        .required_tests(
          "test-02-iwmde-oracles.R",
          c(
            "qCMDE matches GLMM and both estimators match selection bridge factors",
            "qCMDE and IWMDE match the known-V estimate SD boundary bridge factor"
          )
        ),
        .iwmde_fast_required_tests(),
        .required_tests(
          "test-03-loo.R",
          "brma.mv known-V fits expose conditional estimate-unit LOO and WAIC"
        )
      )
    ),
    "tooling" = list(
      description = paste(
        "Maintainer scenario tooling (helpers, run safety, benchmark driver) and",
        "fit-cache producer fixtures: runs that copy the package source tree,",
        "start child R processes or build temporary scenario trees."
      ),
      fit_sources = character(),
      test_filter = anchored("00-scenario-(benchmark|helpers|safety)|02-cache-catalog"),
      required_tests = rbind(
        .required_tests(
          "test-00-scenario-benchmark.R",
          c(
            "benchmark workload and helpers belong to the driver tree",
            "BCG marginal preparation agrees for full and selected synthetic cached inputs",
            "benchmark option defaults and missing call metadata remain explicit",
            "benchmark plots use the canonical scenario SVG renderer",
            "benchmark memory reads MB maxima in ordinary and finite-heap layouts",
            "Bem-shaped fit-list construction keeps its named outputs and loop binding",
            "benchmark cache versions match literals without forcing expressions",
            "benchmark cache loading validates genuine persisted fit contracts",
            "benchmark tools record an unavailable git HEAD as NA"
          )
        ),
        .required_tests(
          "test-00-scenario-helpers.R",
          c(
            "scenario_fit invalidates changed calls and supports regeneration",
            "scenario defaults distinguish direct and managed interactive runs",
            "scenario_time returns values and records successful computations",
            "scenario_fit supports optional targeted cache versions",
            "scenario refitting and output updates are independent",
            "scenario_fit stops at a cached RoBMA fit the current build refuses",
            "scenario_fit loads current RoBMA fits and other cached objects",
            "scenario_text adds, compares, and regenerates tracked output",
            "scenario_text uses portable quotes and restores printing options",
            "scenario_text seeds stochastic snapshots independently",
            "scenario review combines table and figure decisions",
            "scenario_text suppresses evaluation messages",
            "scenario_plot locks interactive plots and seeds each rendering",
            "scenario runs keep locked figures and still write candidates"
          )
        ),
        .required_tests(
          "test-00-scenario-safety.R",
          c(
            "scenario and text names reserve the candidate suffix before evaluation",
            "ordinary text artifacts preserve locked bytes and candidate behavior",
            "plot aliases fail before evaluating or writing a second artifact",
            "default scenario loading prepares the release DLL without test helpers"
          )
        ),
        .required_tests(
          "test-02-cache-catalog.R",
          c(
            "derived cache entries declare and validate consumed artifacts",
            "warmed dependency objects retain first consumption conflicts"
          )
        )
      )
    )
  )
}


certification_case_names <- function() {

  names(certification_cases())
}


certification_case <- function(name) {

  cases <- certification_cases()
  if (!is.character(name) || length(name) != 1L || is.na(name) ||
      !name %in% names(cases)) {
    stop(
      "Unknown certification case: ", paste(name, collapse = ", "),
      ". Available cases: ", paste(names(cases), collapse = ", "),
      call. = FALSE
    )
  }

  return(cases[[name]])
}


certification_case_fit_names <- function(name, catalog = fit_catalog()) {

  case <- certification_case(name)
  if (!is.null(case[["fit_names"]])) {
    unknown <- setdiff(case[["fit_names"]], catalog[["name"]])
    if (length(unknown) > 0L) {
      stop(
        "Unknown cached fit in certification case '", name, "': ",
        paste(unknown, collapse = ", "),
        call. = FALSE
      )
    }

    return(case[["fit_names"]])
  }

  if (length(case[["fit_sources"]]) == 0L) {
    return(character())
  }

  catalog[["name"]][catalog[["source_file"]] %in% case[["fit_sources"]]]
}


certification_case_fit_filter <- function(name) {

  sources <- certification_case(name)[["fit_sources"]]
  if (length(sources) == 0L) {
    return(NULL)
  }

  # testthat matches the filter against the file stems: each stem is matched
  # whole and literally, so that '01-bselmodel' does not also select
  # 'test-01-bselmodel.mv.R'.
  stems <- sub("^test-", "", sources)
  stems <- sub("\\.[Rr]$", "", stems)
  stems <- gsub("([][{}()|^$.*+?\\\\])", "\\\\\\1", stems)
  paste0("^(", paste(stems, collapse = "|"), ")$")
}


certification_case_phases <- function(name) {

  case <- certification_case(name)
  if (length(case[["fit_sources"]]) == 0L) {
    return("verify")
  }

  return(c("prepare", "verify"))
}
