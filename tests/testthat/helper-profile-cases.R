CERTIFICATION_CASE_TIMEOUT_SECONDS <- 60 * 60


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
        "02-(dfbetas|forest|funnel|hatvalues|influence|marginal_means|",
        "plot.*|predict|qqnorm|radial|regplot|residuals|summary.*|vif)|",
        "03-(bridgesampling|loo|zplot)"
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
        "influence|marginal_means|predict|qqnorm|residuals|summary.*|vif)|",
        "03-(bridgesampling|loo)"
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
      test_filter = multivariate_filter
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
      test_filter = multivariate_filter
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
      test_filter = multivariate_filter
    ),
    "multivariate-parity-cs" = list(
      description = "Metafor parity for the Konstantopoulos CS model.",
      fit_sources = "test-01-brma.mv.R",
      fit_names = "brma.mv_v14_konstantopoulos2011_cs",
      test_filter = multivariate_filter
    ),
    "multivariate-parity-nested" = list(
      description = "Metafor parity for the Assink nested model.",
      fit_sources = "test-01-brma.mv.R",
      fit_names = "brma.mv_v14_assink2016_nested",
      test_filter = multivariate_filter
    ),
    "multivariate-parity-har" = list(
      description = "Metafor parity for the Ishak HAR model.",
      fit_sources = "test-01-brma.mv.R",
      fit_names = "brma.mv_v14_ishak2007_har",
      test_filter = multivariate_filter
    ),
    "multivariate-parity-treatment" = list(
      description = "Metafor parity for the Begg treatment CS model.",
      fit_sources = "test-01-brma.mv.R",
      fit_names = "brma.mv_v14_begg1989_study_treatment",
      test_filter = multivariate_filter
    ),
    "iwmde-qcmde" = list(
      description = paste(
        "Fitted and analytic qCMDE/IWMDE density, ordinate, and bridge",
        "oracles."
      ),
      fit_sources = c(
        "test-01-brma.mv.R",
        "test-01-iwmde-oracle-nested.R"
      ),
      fit_names = c(
        "brma.mv_latent",
        "brma.mv_whitened",
        "brma.mv_block_mvn",
        "brma.mv_block_mvn_random",
        "brma.mv_block_mvn_fixed_random_null",
        "nielweise2008_glmm_effect_null",
        "dat.lehmann2018-3PSM_effect_null",
        "iwmde_known_v_tau_full",
        "iwmde_known_v_tau_null"
      ),
      test_filter = paste0(
        "02-(hypothesis|iwmde.*|marginal_means|random-parameters)|",
        "03-(bridgesampling|loo)"
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
