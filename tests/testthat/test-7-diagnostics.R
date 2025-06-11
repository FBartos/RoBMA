context("(7) Diagnostics plots")
skip_on_cran()

# test objects - assuming that the fit function worked properly
temp_fits_dir <- Sys.getenv("ROBMA_TEST_FITS_DIR")
if (temp_fits_dir == "" || !dir.exists(temp_fits_dir)) {
  stop("Temporary fits directory not found. Run test-4-fit.R first.")
}

saved_files <- paste0("fit_", c(1, 15), ".RDS")

# Read all fits at the beginning for efficiency
fits <- list()
fits[["fit_1"]] <- fits[["fit_1"]]
fits[["fit_15"]] <- fits[["fit_15"]]

test_that("Diagnostic plots work", {

  # RoBMA
  chains_mu    <- diagnostics(fits[["fit_1"]], "mu",    "chains", plot_type = "ggplot")
  chains_tau   <- diagnostics(fits[["fit_1"]], "tau",   "chains", plot_type = "ggplot")
  chains_omega <- diagnostics(fits[["fit_1"]], "omega", "chains", plot_type = "ggplot")
  chains_PET   <- diagnostics(fits[["fit_1"]], "PET",   "chains", plot_type = "ggplot")
  chains_PEESE <- diagnostics(fits[["fit_1"]], "PEESE", "chains", plot_type = "ggplot")

  autocorrelation_mu   <- diagnostics(fits[["fit_1"]], "mu",    "autocorrelation", plot_type = "ggplot")
  densities_mu         <- diagnostics(fits[["fit_1"]], "mu",    "densities", plot_type = "ggplot")


  expect_equal(all(sapply(chains_mu[1:18],  is.null)), TRUE)
  expect_equal(all(sapply(chains_mu[19:36], ggplot2::is.ggplot)), TRUE)
  expect_equal(all(sapply(chains_tau[c(1:9,19:27)],    is.null)), TRUE)
  expect_equal(all(sapply(chains_tau[c(10:18,28:36)],  ggplot2::is.ggplot)), TRUE)
  expect_equal(all(sapply(chains_omega[c(1,8:9,10,17:18,19,26:27,28,35:36)], is.null)), TRUE)
  expect_equal(all(sapply(chains_omega[c(2,4,11,13,20,22,29,31)], ggplot2::is.ggplot)), TRUE)
  expect_equal(sum(unlist(sapply(chains_omega[c(3,5:7,12,14:16,21,23:25,30,32:34)], function(l) sapply(l, ggplot2::is.ggplot))) ), 36)
  expect_equal(all(sapply(chains_PET[c(1:7,9:16,18:25,27:34,36)], is.null)), TRUE)
  expect_equal(all(sapply(chains_PET[c(8,17,26,35)],  ggplot2::is.ggplot)), TRUE)
  expect_equal(all(sapply(chains_PEESE[c(1:8,10:17,19:26,28:35)], is.null)), TRUE)
  expect_equal(all(sapply(chains_PEESE[c(9,18,27,36)],  ggplot2::is.ggplot)), TRUE)

  expect_equal(all(sapply(autocorrelation_mu[1:18],  is.null)), TRUE)
  expect_equal(all(sapply(autocorrelation_mu[19:36], ggplot2::is.ggplot)), TRUE)
  expect_equal(all(sapply(densities_mu[1:18],  is.null)), TRUE)
  expect_equal(all(sapply(densities_mu[19:36], ggplot2::is.ggplot)), TRUE)


  vdiffr::expect_doppelganger(paste0("ggplot_chains_mu"),          chains_mu[[36]])
  vdiffr::expect_doppelganger(paste0("ggplot_autocorrelation_mu"), autocorrelation_mu[[36]])
  vdiffr::expect_doppelganger(paste0("ggplot_densities_mu"),       densities_mu[[36]])
  vdiffr::expect_doppelganger(paste0("ggplot_chains_tau"),         chains_tau[[36]])
  for(j in 1:3){
    vdiffr::expect_doppelganger(paste0("ggplot_chains_omega_",j),  chains_omega[[34]][[j]])
  }
  vdiffr::expect_doppelganger(paste0("ggplot_chains_PET"),         chains_PET[[35]])
  vdiffr::expect_doppelganger(paste0("ggplot_chains_PEESE"),       chains_PEESE[[36]])

  # base plots
  vdiffr::expect_doppelganger(paste0("plot_chains_mu"),          function()diagnostics(fits[["fit_1"]], "mu", "chains",          show_models = 36))
  vdiffr::expect_doppelganger(paste0("plot_autocorrelation_mu"), function()diagnostics(fits[["fit_1"]], "mu", "autocorrelation", show_models = 36))
  vdiffr::expect_doppelganger(paste0("plot_densities_mu"),       function()diagnostics(fits[["fit_1"]], "mu", "densities",       show_models = 36))

  vdiffr::expect_doppelganger(paste0("plot_chains_tau"),          function()diagnostics(fits[["fit_1"]], "tau", "chains",          show_models = 10))
  vdiffr::expect_doppelganger(paste0("plot_autocorrelation_tau"), function()diagnostics(fits[["fit_1"]], "tau", "autocorrelation", show_models = 10))
  vdiffr::expect_doppelganger(paste0("plot_densities_tau"),       function()diagnostics(fits[["fit_1"]], "tau", "densities",       show_models = 10))

  vdiffr::expect_doppelganger(paste0("plot_chains_omega"),        function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(3,1))
    suppressMessages(diagnostics(fits[["fit_1"]], "omega", "chains", show_models = 7))
  })
  vdiffr::expect_doppelganger(paste0("plot_autocorrelation_omega"),        function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(3,1))
    suppressMessages(diagnostics(fits[["fit_1"]], "omega", "autocorrelation", show_models = 7))
  })
  vdiffr::expect_doppelganger(paste0("plot_densities_omega"),        function(){
    oldpar <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(mfrow = oldpar[["mfrow"]]))
    par(mfrow = c(3,1))
    suppressMessages(diagnostics(fits[["fit_1"]], "omega", "densities", show_models = 7))
  })

  vdiffr::expect_doppelganger(paste0("plot_chains_omega_2"),           function()suppressMessages(diagnostics(fits[["fit_1"]], "omega", "chains",          show_models = 2)))
  vdiffr::expect_doppelganger(paste0("plot_autocorrelation_omega_2"),  function()suppressMessages(diagnostics(fits[["fit_1"]], "omega", "autocorrelation", show_models = 2)))
  vdiffr::expect_doppelganger(paste0("plot_densities_omega_2"),        function()suppressMessages(diagnostics(fits[["fit_1"]], "omega", "densities",       show_models = 2)))

  vdiffr::expect_doppelganger(paste0("plot_chains_PET"),          function()diagnostics(fits[["fit_1"]], "PET", "chains",          show_models = 8))
  vdiffr::expect_doppelganger(paste0("plot_autocorrelation_PET"), function()diagnostics(fits[["fit_1"]], "PET", "autocorrelation", show_models = 8))
  vdiffr::expect_doppelganger(paste0("plot_densities_PET"),       function()diagnostics(fits[["fit_1"]], "PET", "densities",       show_models = 8))

  vdiffr::expect_doppelganger(paste0("plot_chains_PEESE"),          function()diagnostics(fits[["fit_1"]], "PEESE", "chains",          show_models = 36))
  vdiffr::expect_doppelganger(paste0("plot_autocorrelation_PEESE"), function()diagnostics(fits[["fit_1"]], "PEESE", "autocorrelation", show_models = 36))
  vdiffr::expect_doppelganger(paste0("plot_densities_PEESE"),       function()diagnostics(fits[["fit_1"]], "PEESE", "densities",       show_models = 36))

  ### RoBMA.reg
  chains_mu       <- diagnostics(fits[["fit_15"]], "mu",        "chains", plot_type = "ggplot")
  chains_omega    <- diagnostics(fits[["fit_15"]], "omega",     "chains", plot_type = "ggplot")
  chains_PET      <- diagnostics(fits[["fit_15"]], "PET",       "chains", plot_type = "ggplot")
  chains_mod_con  <- diagnostics(fits[["fit_15"]], "mod_con",   "chains", plot_type = "ggplot")

  vdiffr::expect_doppelganger("plot_chains.reg_mu",      chains_mu)
  vdiffr::expect_doppelganger("plot_chains.reg_omega",   chains_omega[[1]])
  vdiffr::expect_doppelganger("plot_chains.reg_PET",     chains_PET)
  vdiffr::expect_doppelganger("plot_chains.reg_mod_con", chains_mod_con)
})
