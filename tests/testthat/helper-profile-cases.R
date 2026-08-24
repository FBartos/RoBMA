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

  if (is.null(required_tests) || nrow(required_tests) == 0L) {
    return(invisible(TRUE))
  }

  results <- as.data.frame(results)
  required_columns <- c("file", "test", "skipped", "passed")
  missing_columns  <- setdiff(required_columns, names(results))
  if (length(missing_columns) > 0L) {
    stop(
      "Certification results are missing required columns: ",
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
        "Certification case '", case_name, "' requires exactly one result for ",
        label, ", but observed ", length(matches), ".",
        call. = FALSE
      )
    }
    if (isTRUE(results[["skipped"]][[matches]])) {
      stop(
        "Certification case '", case_name, "' skipped required test ", label,
        ".",
        call. = FALSE
      )
    }
    if (is.na(results[["passed"]][[matches]]) ||
        results[["passed"]][[matches]] < 1L) {
      stop(
        "Certification case '", case_name,
        "' recorded no passing expectations for required test ", label, ".",
        call. = FALSE
      )
    }
  }

  return(invisible(TRUE))
}


certification_cases <- function() {

  multivariate_filter <- paste0(
    "02-(brma-mv.*|derived-random-correlations|forest|funnel|",
    "heterogeneity-mv|hypothesis|influence|marginal_means|plot.*|",
    "predict-mv|qqnorm|random-parameters|regplot|residuals|summary.*|vif)|",
    "03-(bridgesampling|loo|zplot)"
  )

  list(
    "numerical-kernels" = list(
      description = paste(
        "Independent selected-normal, known-V, distribution, and GLMM",
        "quadrature oracles."
      ),
      fit_sources = character(),
      test_filter = paste0(
        "00-(covariance-factorization|known-v-joint-loglik|",
        "selection-kernel.*|selection-probability-numerics)|",
        "02-(distributions|glmm-aghq)"
      ),
      required_tests = .required_tests(
        "test-00-covariance-factorization.R",
        "covariance sampling factor preserves singular covariance"
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
        "test-01-vif-parity.R"
      ),
      test_filter = paste0(
        "02-(dfbetas|forest|funnel|hatvalues|influence|iwmde-fast-paths|marginal_means|",
        "plot.*|predict|qqnorm|radial|regplot|residuals|summary.*|vif)|",
        "03-(bridgesampling|loo|zplot)"
      ),
      required_tests = rbind(
        .required_tests(
          "test-02-iwmde-fast-paths.R",
          "multilevel weightfunction formula path matches scalar fallback"
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
      test_filter = paste0(
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
      test_filter = paste(multivariate_filter, "02-iwmde-api", sep = "|"),
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
        "test-01-brma.glmm.R",
        "test-01-bselmodel.R",
        "test-01-brma.mv.R",
        "test-01-iwmde-oracle-nested.R"
      ),
      fit_names = c(
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
      test_filter = paste0(
        "02-(hypothesis|iwmde.*|marginal_means|random-parameters)|",
        "03-(bridgesampling|loo)"
      ),
      required_tests = rbind(
        .required_tests(
          "test-02-iwmde-oracles.R",
          c(
            "qCMDE matches GLMM and both estimators match selection bridge factors",
            "qCMDE and IWMDE match the known-V tau boundary bridge factor"
          )
        ),
        .iwmde_fast_required_tests(),
        .required_tests(
          "test-03-loo.R",
          "brma.mv known-V fits expose conditional estimate-unit LOO and WAIC"
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

  stems <- sub("^test-", "", sources)
  stems <- sub("\\.[Rr]$", "", stems)
  paste(stems, collapse = "|")
}
