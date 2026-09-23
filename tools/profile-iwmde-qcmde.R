# Profile qCMDE/IWMDE ordinate evaluation against validated cached fits.

devtools::load_all(quiet = TRUE)
source(testthat::test_path("common-functions.R"))

.profile_iwmde_estimator <- function(name, likelihood_family, parameter, value,
                                     max_samples, density_method = "qCMDE",
                                     n_points = 40L,
                                     normalization_points = 101L) {

  fit     <- load_fit(name, validate = TRUE)
  context <- .iwmde_context(fit)
  memory  <- tempfile(fileext = ".out")
  cpu     <- tempfile(fileext = ".out")

  gc()
  Rprofmem(memory)
  utils::Rprof(cpu, interval = .005)
  on.exit({
    Rprofmem(NULL)
    utils::Rprof(NULL)
    unlink(c(memory, cpu))
  }, add = TRUE)

  elapsed <- system.time({
    estimate <- .iwmde_estimate(
      context         = context,
      parameter       = parameter,
      density_method  = density_method,
      density_control = list(
        n_points             = n_points,
        max_samples          = max_samples,
        normalization_points = normalization_points,
        normalization_prob   = .9999,
        display_grid         = "ordinate"
      ),
      outputs          = "ordinate",
      values           = value,
      parameter_spec   = list(
        type             = "primitive",
        conditional      = NULL,
        conditional_rule = "AND"
      ),
      metadata         = list(parameter = parameter),
      cache            = .iwmde_estimate_cache()
    )
  })[["elapsed"]]
  Rprofmem(NULL)
  utils::Rprof(NULL)

  allocations <- readLines(memory, warn = FALSE)
  sizes       <- suppressWarnings(as.numeric(sub(" .*", "", allocations)))
  finite      <- sizes[is.finite(sizes)]
  profile_summary <- utils::summaryRprof(cpu)
  profile         <- list(
    by_total = utils::head(profile_summary[["by.total"]], 20L),
    by_self  = utils::head(profile_summary[["by.self"]], 20L)
  )
  ordinate <- estimate[["posterior_ordinate"]]
  status   <- "accepted"
  if (is.null(ordinate)) {
    ordinate <- estimate[["rejected_posterior_ordinate"]]
    status   <- "rejected"
  }
  if (is.null(ordinate)) {
    stop("Profiled ordinate was unavailable for '", name, "'.", call. = FALSE)
  }
  diagnostics      <- ordinate[["diagnostics"]]
  evaluation_value <- diagnostics[["evaluation_value"]]

  result <- data.frame(
    fit                   = name,
    likelihood_family     = likelihood_family,
    parameter             = parameter,
    density_method        = density_method,
    status                = status,
    requested_value       = value,
    evaluation_value      = evaluation_value,
    max_samples           = max_samples,
    n_points              = n_points,
    normalization_points  = normalization_points,
    elapsed_seconds       = elapsed,
    allocated_megabytes   = sum(finite) / 1024^2,
    largest_allocation_mb = if (length(finite) > 0L) {
      max(finite) / 1024^2
    } else {
      NA_real_
    },
    posterior_ordinate    = ordinate[["ordinate"]],
    stringsAsFactors      = FALSE
  )
  attr(result, "profile") <- profile

  return(result)
}

cases <- list(
  list(
    name                 = "brma.mv_block_mvn_fixed_random_null",
    likelihood_family    = "known-V normal",
    parameter            = "mu",
    value                = 0,
    max_samples          = 40L,
    normalization_points = 51L
  ),
  list(
    name                 = "brma.mv_block_mvn_fixed_random_null",
    likelihood_family    = "known-V normal",
    parameter            = "mu",
    value                = 0,
    max_samples          = 240L,
    normalization_points = 51L
  ),
  list(
    name                 = "brma.mv_block_mvn_fixed_random_null",
    likelihood_family    = "known-V normal",
    parameter            = "mu",
    value                = 0,
    max_samples          = 240L,
    n_points             = 80L,
    normalization_points = 201L
  ),
  list(
    name                 = "brma.mv_block_mvn",
    likelihood_family    = "known-V normal",
    parameter            = "tau",
    value                = 0,
    max_samples          = 240L,
    n_points             = 60L,
    normalization_points = 120L
  ),
  list(
    name                 = "brma.mv_block_mvn",
    likelihood_family    = "known-V normal",
    parameter            = "tau",
    value                = 0,
    max_samples          = 240L,
    density_method       = "IWMDE",
    n_points             = 60L,
    normalization_points = 120L
  ),
  list(
    name              = "bcg_meta-analysis",
    likelihood_family = "ordinary normal",
    parameter         = "mu",
    value             = 0,
    max_samples       = 500L
  ),
  list(
    name              = "nielweise2008_glmm",
    likelihood_family = "Poisson GLMM",
    parameter         = "mu",
    value             = 0,
    max_samples       = 500L
  ),
  list(
    name              = "dat.lehmann2018-3PSM",
    likelihood_family = "selection normal",
    parameter         = "mu",
    value             = 0,
    max_samples       = 500L
  )
)

available <- list_fits(validate = TRUE)
cases     <- Filter(function(case) {

  case[["name"]] %in% available
}, cases)
if (length(cases) == 0L) {
  stop(
    "Required cached fits are unavailable; run the test-01 cache suite first.",
    call. = FALSE
  )
}

results <- lapply(cases, function(case) {

  do.call(.profile_iwmde_estimator, case)
})
for (result in results) {
  print(result, row.names = FALSE, digits = 6)
  cat("Top functions by total CPU time:\n")
  print(attr(result, "profile")[["by_total"]], digits = 4)
  cat("Top functions by self CPU time:\n")
  print(attr(result, "profile")[["by_self"]], digits = 4)
}
