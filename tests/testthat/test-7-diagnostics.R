context("(7) diagnostics plots")
skip_on_cran()


# test objects - assuming that the fit function worked properly
saved_fits    <- readRDS(file = "saved_fits.RDS")
set.seed(1)


test_that("Parameter plots works", {

  chains_mu    <- diagnostics(saved_fits[[1]], "mu",    "chains", plot_type = "ggplot")
  chains_tau   <- diagnostics(saved_fits[[1]], "tau",   "chains", plot_type = "ggplot")
  chains_omega <- diagnostics(saved_fits[[1]], "omega", "chains", plot_type = "ggplot")

  autocorrelation_mu    <- diagnostics(saved_fits[[1]], "mu",    "autocorrelation", plot_type = "ggplot")
  autocorrelation_tau   <- diagnostics(saved_fits[[1]], "tau",   "autocorrelation", plot_type = "ggplot")
  autocorrelation_omega <- diagnostics(saved_fits[[1]], "omega", "autocorrelation", plot_type = "ggplot")

  densities_mu    <- diagnostics(saved_fits[[1]], "mu",    "densities", plot_type = "ggplot")
  densities_tau   <- diagnostics(saved_fits[[1]], "tau",   "densities", plot_type = "ggplot")
  densities_omega <- diagnostics(saved_fits[[1]], "omega", "densities", plot_type = "ggplot")


  expect_equal(all(sapply(chains_mu[1:6],  is.null)), TRUE)
  expect_equal(all(sapply(chains_mu[7:12], ggplot2::is.ggplot)), TRUE)
  expect_equal(all(sapply(chains_tau[c(1,2,3,7,8,9)],    is.null)), TRUE)
  expect_equal(all(sapply(chains_tau[c(4,5,6,10,11,12)], ggplot2::is.ggplot)), TRUE)
  expect_equal(all(sapply(chains_omega[c(1, 4, 7, 10)], is.null)), TRUE)
  expect_equal(all(sapply(chains_omega[c(2, 5, 8, 11)], ggplot2::is.ggplot)), TRUE)
  expect_equal(all(sapply(chains_omega[c(3, 6, 9, 12)], length) == 2), TRUE)

  expect_equal(all(sapply(autocorrelation_mu[1:6],  is.null)), TRUE)
  expect_equal(all(sapply(autocorrelation_mu[7:12], ggplot2::is.ggplot)), TRUE)
  expect_equal(all(sapply(autocorrelation_tau[c(1,2,3,7,8,9)],    is.null)), TRUE)
  expect_equal(all(sapply(autocorrelation_tau[c(4,5,6,10,11,12)], ggplot2::is.ggplot)), TRUE)
  expect_equal(all(sapply(autocorrelation_omega[c(1, 4, 7, 10)], is.null)), TRUE)
  expect_equal(all(sapply(autocorrelation_omega[c(2, 5, 8, 11)], ggplot2::is.ggplot)), TRUE)
  expect_equal(all(sapply(autocorrelation_omega[c(3, 6, 9, 12)], length) == 2), TRUE)

  expect_equal(all(sapply(densities_mu[1:6],  is.null)), TRUE)
  expect_equal(all(sapply(densities_mu[7:12], ggplot2::is.ggplot)), TRUE)
  expect_equal(all(sapply(densities_tau[c(1,2,3,7,8,9)],    is.null)), TRUE)
  expect_equal(all(sapply(densities_tau[c(4,5,6,10,11,12)], ggplot2::is.ggplot)), TRUE)
  expect_equal(all(sapply(densities_omega[c(1, 4, 7, 10)], is.null)), TRUE)
  expect_equal(all(sapply(densities_omega[c(2, 5, 8, 11)], ggplot2::is.ggplot)), TRUE)
  expect_equal(all(sapply(densities_omega[c(3, 6, 9, 12)], length) == 2), TRUE)


  for(i in 1:12){

    if(i %in% c(7:12)){
      expect_doppelganger(paste0("chains_mu",i),          chains_mu[[i]])
      expect_doppelganger(paste0("autocorrelation_mu",i), autocorrelation_mu[[i]])
      expect_doppelganger(paste0("densities_mu",i),       densities_mu[[i]])
    }

    if(i %in% c(4,5,6,10,11,12)){
      expect_doppelganger(paste0("chains_tau",i),          chains_tau[[i]])
      expect_doppelganger(paste0("autocorrelation_tau",i), autocorrelation_tau[[i]])
      expect_doppelganger(paste0("densities_tau",i),       densities_tau[[i]])
    }

    if(i %in% c(2, 5, 8, 11)){
      expect_doppelganger(paste0("chains_omega",i),          chains_omega[[i]])
      expect_doppelganger(paste0("autocorrelation_omega",i), autocorrelation_omega[[i]])
      expect_doppelganger(paste0("densities_omega",i),       densities_omega[[i]])
    }

    if(i %in% c(3, 6, 9, 12)){
      for(j in 1:2){
        expect_doppelganger(paste0("chains_omega",i,"_j_",j),          chains_omega[[i]][[j]])
        expect_doppelganger(paste0("autocorrelation_omega",i,"_j_",j), autocorrelation_omega[[i]][[j]])
        expect_doppelganger(paste0("densities_omega",i,"_j_",j),       densities_omega[[i]][[j]])
      }
    }

  }


  expect_error(diagnostics(saved_fits[[2]], "mu", "chains"),"Diagnostics cannot be produced because individual model posteriors were not save during the fitting process.")

})
