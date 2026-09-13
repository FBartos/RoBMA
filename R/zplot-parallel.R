# Split only posterior rows. Every worker retains the complete publication
# events, their original block order, and the supplied integration design seed.
.zplot_selection_marginal_parallel <- function(object, posterior_samples,
    z_sequence, conditioning_depth, integration_control, cores) {

  # Optional continuation snapshots do not belong in each worker payload.
  attr(object[["fit"]], "runtime_state") <- NULL
  BayesTools::check_int(cores, "cores", lower = 2L, allow_NA = FALSE)
  S <- nrow(posterior_samples)
  if (S < cores) {
    stop("Parallel zplot computation requires at least one posterior draw per worker.", call. = FALSE)
  }
  control <- .check_selection_likelihood_control(integration_control)
  memory_options <- list(
    RoBMA.known_v_covariance_max_bytes = .known_v_covariance_max_bytes(),
    BayesTools.random_effects_memory_limit_bytes = BayesTools:::.bt_random_effect_memory_limit_bytes()
  )
  runtime_setup <- .selection_runtime_setup()
  runtime_started <- FALSE
  message("Computing zplot densities using ", cores, " parallel workers.")
  cluster <- parallel::makePSOCKcluster(cores, rscript_args = "--vanilla")
  on.exit(BayesTools:::.JAGS_finish_runtime_setup(
    if (runtime_started) runtime_setup else NULL, chains = cores, cl = cluster,
    operation = "Parallel zplot density computation"), add = TRUE)
  initialize <- function(paths, memory_options) {

    .libPaths(paths)
    options(memory_options)
    NULL
  }
  environment(initialize) <- baseenv()
  parallel::clusterCall(cluster, initialize, .libPaths(), memory_options)
  BayesTools:::.JAGS_require_packages(c("BayesTools", "RoBMA"), cluster,
    operation = "Parallel zplot density computation")
  runtime_started <- TRUE
  BayesTools:::.JAGS_run_runtime_setup(runtime_setup, chains = cores, cl = cluster)

  worker <- function(rows, object, posterior_samples, z_sequence,
                     conditioning_depth, integration_control) {

    RoBMA:::.zplot_selection_marginal(
      object = object, posterior_samples = posterior_samples[rows, , drop = FALSE],
      z_sequence = z_sequence, z_threshold = NULL,
      conditioning_depth = conditioning_depth,
      integration_control = integration_control, extrapolate_only = FALSE)
  }
  environment(worker) <- baseenv()
  group <- as.integer(cut(seq_len(S), breaks = cores, labels = FALSE))
  chunks <- unname(split(seq_len(S), group))
  pieces <- parallel::clusterApply(cluster, chunks, worker,
    object = object, posterior_samples = posterior_samples, z_sequence = z_sequence,
    conditioning_depth = conditioning_depth, integration_control = control)
  result <- list(fitted = matrix(NA_real_, S, length(z_sequence)),
    extrapolated = matrix(NA_real_, S, length(z_sequence)), weights = numeric(S), EDR = NULL)
  for (index in seq_along(chunks)) {
    rows <- chunks[[index]]
    value <- pieces[[index]]
    if (!is.list(value) ||
        !identical(dim(value[["fitted"]]), c(length(rows), length(z_sequence))) ||
        !identical(dim(value[["extrapolated"]]), c(length(rows), length(z_sequence))) ||
        length(value[["weights"]]) != length(rows) || !is.null(value[["EDR"]])) {
      stop("Parallel zplot density results are unavailable because a worker returned inconsistent dimensions.",
        call. = FALSE)
    }
    result$fitted[rows, ] <- value$fitted
    result$extrapolated[rows, ] <- value$extrapolated
    result$weights[rows] <- value$weights
  }
  result
}
