context("Exact observation-deletion refits for LOO")

source(testthat::test_path("common-functions.R"))
skip_on_cran()
skip_if_not_certification(
  "Exact LOO refits are reserved for their dedicated certification case."
)


.exact_loo_convergence_checks <- function() {

  # These scenario models intentionally disable generic parameter thresholds.
  # The certification below diagnoses the predictive integral directly through
  # its batch-means MCSE and the full-fit PSIS diagnostics.
  set_convergence_checks(max_Rhat = NULL, min_ESS = NULL)
}


.exact_loo_kearon_data <- function() {

  data("dat.kearon1998", package = "metadat", envir = environment())
  dat <- suppressWarnings(metafor::to.long(
    measure = "OR",
    ai      = tp,
    n1i     = np,
    ci      = tn,
    n2i     = nn,
    data    = dat.kearon1998,
    subset  = patients == "asymptomatic",
    append  = FALSE
  ))
  dat$group <- factor(
    dat$group,
    levels = c(1, 2),
    labels = c("sensitivity", "specificity")
  )
  dat <- metafor::escalc(
    measure = "PLO",
    xi      = out1,
    mi      = out2,
    data    = dat,
    add     = 1 / 2,
    to      = "all"
  )

  return(dat)
}


.exact_loo_fit_kearon <- function(dat, structure, omitted, seed) {

  structure <- match.arg(structure, c("us", "hcs"))
  keep      <- setdiff(seq_len(nrow(dat)), omitted)
  fit_data  <- dat[keep, , drop = FALSE]
  uisd      <- estimate_unit_information_sd(
    sei = sqrt(dat$vi),
    ni  = dat$out1 + dat$out2
  )
  random_formula <- switch(
    structure,
    us  = ~ us(0 + group | study),
    hcs = ~ hcs(group | study)
  )
  fit_args <- list(
    yi                             = fit_data$yi,
    V                              = fit_data$vi,
    ni                             = fit_data$out1 + fit_data$out2,
    mods                           = ~ 0 + group,
    random                         = random_formula,
    set_contrast_factor_predictors = "independent",
    prior_unit_information_sd      = uisd,
    data                           = fit_data,
    measure                        = "GEN",
    chains                         = 3,
    sample                         = 5000,
    burnin                         = 2000,
    adapt                          = 500,
    seed                           = seed,
    silent                         = TRUE,
    convergence_checks             = .exact_loo_convergence_checks()
  )
  if (structure == "us") {
    fit_args[["prior_heterogeneity"]] <- prior_random(
      study = random_block(contrasts = c(group = "independent"))
    )
  }

  return(do.call(brma.mv, fit_args))
}


.exact_loo_ishak_data <- function() {

  data("dat.ishak2007", package = "metadat", envir = environment())
  dat <- reshape(
    dat.ishak2007,
    direction = "long",
    idvar     = "study",
    v.names   = c("yi", "vi"),
    varying   = list(c(2, 4, 6, 8), c(3, 5, 7, 9))
  )
  dat <- dat[order(dat$study, dat$time), ]
  dat <- dat[!is.na(dat$yi), ]
  rownames(dat) <- NULL
  dat$time_factor <- factor(dat$time)
  V <- metafor::vcalc(
    vi,
    cluster = study,
    time1   = time,
    phi     = 0.97,
    data    = dat
  )

  return(list(data = dat, V = V))
}


.exact_loo_fit_ishak_har <- function(dat, V, omitted, seed) {

  keep     <- setdiff(seq_len(nrow(dat)), omitted)
  fit_data <- dat[keep, , drop = FALSE]
  fit_V    <- V[keep, keep, drop = FALSE]
  uisd     <- estimate_unit_information_sd(sei = sqrt(14.3), ni = 15)
  time_prior <- prior_factor(
    distribution = "normal",
    parameters   = list(mean = 0, sd = uisd),
    contrast     = "independent"
  )

  return(brma.mv(
    yi                        = fit_data$yi,
    V                         = fit_V,
    measure                   = "GEN",
    mods                      = ~ 0 + time_factor,
    random                    = ~ har(time | study),
    prior_effect              = NULL,
    prior_mods                = time_prior,
    prior_unit_information_sd = uisd,
    data                      = fit_data,
    chains                    = 3,
    sample                    = 10000,
    burnin                    = 4000,
    adapt                     = 2000,
    seed                      = seed,
    silent                    = TRUE,
    convergence_checks        = .exact_loo_convergence_checks()
  ))
}


.exact_loo_log_mean_density <- function(log_density, chain_id) {

  if (length(log_density) != length(chain_id) ||
      length(log_density) < 4L || any(!is.finite(log_density)) ||
      anyNA(chain_id)) {
    stop("Exact-LOO density draws and chain IDs are invalid.", call. = FALSE)
  }

  offset        <- max(log_density)
  contributions <- exp(log_density - offset)
  density_mean  <- mean(contributions)
  chain_rows    <- split(seq_along(chain_id), chain_id)
  chain_lengths <- lengths(chain_rows)
  batch_size    <- max(2L, floor(sqrt(min(chain_lengths))))
  batch_counts  <- floor(chain_lengths / batch_size)
  if (any(batch_counts < 2L)) {
    stop("Exact-LOO refit chains do not contain enough complete batches.",
         call. = FALSE)
  }

  batch_means <- numeric(sum(batch_counts))
  batch       <- 0L
  for (chain in seq_along(chain_rows)) {
    rows <- chain_rows[[chain]]
    for (chain_batch in seq_len(batch_counts[[chain]])) {
      batch      <- batch + 1L
      batch_rows <- rows[
        (chain_batch - 1L) * batch_size + seq_len(batch_size)
      ]
      batch_means[[batch]] <- mean(contributions[batch_rows])
    }
  }

  density_mcse <- stats::sd(batch_means) / sqrt(length(batch_means))
  return(list(
    estimate   = offset + log(density_mean),
    mcse       = density_mcse / density_mean,
    batch_size = batch_size,
    n_batches  = length(batch_means)
  ))
}


.exact_loo_refit_score <- function(full_fit, refit, heldout) {

  refit_draws <- BayesTools::JAGS_materialize_draws(
    refit[["fit"]],
    include_internal = TRUE
  )
  posterior_samples <- as.matrix(refit_draws)
  chain_id <- rep(
    seq_along(refit_draws),
    vapply(refit_draws, nrow, integer(1))
  )
  setup <- .estimate_likelihood_setup_from_parts(
    fit                     = full_fit[["fit"]],
    data                    = full_fit[["data"]],
    priors                  = full_fit[["priors"]],
    posterior_samples       = posterior_samples,
    unit                    = "estimate",
    condition_local_effects = FALSE,
    caller                  = ".exact_loo_refit_score()"
  )
  log_density <- .log_lik_estimate_from_setup(setup)[, heldout]

  return(.exact_loo_log_mean_density(log_density, chain_id))
}


.exact_loo_good_k_threshold <- function(n_draws) {

  # Sample-size-dependent reliability threshold used by loo 2.10 diagnostics.
  return(min(1 - 1 / log10(n_draws), 0.7))
}


test_that("scenario-model PSIS LOO agrees with five exact deletion refits", {

  skip_if_not_installed("metadat")
  skip_if_not_installed("metafor")
  skip_if_not_installed("loo", minimum_version = "2.10.0")

  kearon <- .exact_loo_kearon_data()
  ishak  <- .exact_loo_ishak_data()
  models <- list(
    kearon_us = list(
      deletions = c(11L, 12L),
      fit = function(omitted, seed) {
        .exact_loo_fit_kearon(kearon, "us", omitted, seed)
      }
    ),
    kearon_hcs = list(
      deletions = 8L,
      fit = function(omitted, seed) {
        .exact_loo_fit_kearon(kearon, "hcs", omitted, seed)
      }
    ),
    ishak_har = list(
      # Interior observations from two complete four-time-point studies.
      deletions = c(42L, 78L),
      fit = function(omitted, seed) {
        .exact_loo_fit_ishak_har(
          ishak[["data"]], ishak[["V"]], omitted, seed
        )
      }
    )
  )

  results <- vector("list", sum(vapply(
    models,
    function(model) length(model[["deletions"]]),
    integer(1)
  )))
  result_index <- 0L
  model_index  <- 0L
  for (model_name in names(models)) {
    model_index <- model_index + 1L
    model       <- models[[model_name]]
    full_fit    <- model[["fit"]](integer(), 8100L + 100L * model_index)
    expect_no_warning(
      full_fit <- add_loo(full_fit),
      info = model_name
    )
    loo_result  <- loo(full_fit)
    full_log_lik <- log_lik(full_fit)
    target <- attr(full_log_lik, "RoBMA_target", exact = TRUE)

    expect_identical(target[["unit"]], "estimate", info = model_name)
    expect_identical(
      target[["retained_context"]],
      "remaining_data",
      info = model_name
    )
    expect_identical(
      target[["target"]],
      "estimate_log_score",
      info = model_name
    )

    pointwise <- loo_result[["pointwise"]]
    pareto_k  <- loo::pareto_k_values(loo_result)
    threshold <- .exact_loo_good_k_threshold(nrow(full_log_lik))
    for (heldout in model[["deletions"]]) {
      result_index <- result_index + 1L
      label <- paste0(model_name, " row ", heldout)
      expect_lt(pareto_k[[heldout]], threshold, label = label)

      refit <- model[["fit"]](heldout, 8100L + 100L * model_index + heldout)
      expect_equal(nobs(refit), nobs(full_fit) - 1L, info = label)
      exact <- .exact_loo_refit_score(full_fit, refit, heldout)
      psis_estimate <- pointwise[heldout, "elpd_loo"]
      psis_mcse     <- pointwise[heldout, "mcse_elpd_loo"]
      combined_mcse <- sqrt(psis_mcse^2 + exact[["mcse"]]^2)
      standardized_difference <-
        abs(psis_estimate - exact[["estimate"]]) / combined_mcse

      expect_true(is.finite(psis_mcse) && psis_mcse > 0, info = label)
      expect_true(
        is.finite(exact[["mcse"]]) && exact[["mcse"]] > 0,
        info = label
      )
      # A combined MCSE below 0.05 keeps the four-SE certification band below
      # 0.2 log predictive-density units; a noisier oracle is not evidence.
      expect_lt(combined_mcse, 0.05, label = label)
      expect_lte(standardized_difference, 4, label = label)

      results[[result_index]] <- data.frame(
        model        = model_name,
        heldout      = heldout,
        psis         = psis_estimate,
        exact_refit  = exact[["estimate"]],
        difference   = psis_estimate - exact[["estimate"]],
        combined_mcse = combined_mcse,
        z            = standardized_difference,
        pareto_k     = pareto_k[[heldout]],
        stringsAsFactors = FALSE
      )
      rm(refit)
      invisible(gc())
    }
    rm(full_fit)
    invisible(gc())
  }

  results <- do.call(rbind, results)
  expect_equal(nrow(results), 5L)
  expect_setequal(results[["model"]], names(models))
})
