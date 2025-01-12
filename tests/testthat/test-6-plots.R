context("(6) Plot functions")
skip_on_cran()

# the plotting functions are imported from BayesTools and tested henceforth
# test objects - assuming that the fit function worked properly
saved_files <- paste0("fit_", 1:16, ".RDS")
saved_fits  <- list()
for(i in seq_along(saved_files)){
  saved_fits[[i]] <- readRDS(file = file.path("../results/fits", saved_files[i]))
}

# alternative components present in the models:
effect             <- c(1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15, 16)
heterogeneity      <- c(1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13, 14, 16)
weightfunctions    <- c(1, 2, 4, 5, 6, 7, 10, 11, 13, 15)
PETPEESE           <- c(1, 3, 4, 5, 6, 7, 11, 13, 15)
no_weightfunctions <- c(3, 8, 12, 14, 16)
no_PETPEESE        <- c(2, 8, 10, 12, 14, 16)
no_PET             <- c(2, 8, 10, 12, 14, 16)
no_PEESE           <- c(2, 7, 8, 9, 10, 11, 12, 14, 15, 16)
metaregression     <- c(14, 15)

test_that("Parameter plots work", {

  ### effect
  # default ggplot2
  for(i in 1:length(saved_fits)){

    vdiffr::expect_doppelganger(paste0("ggplot_mu1_",i), plot(saved_fits[[i]], "mu", plot_type = "ggplot"))
    vdiffr::expect_doppelganger(paste0("ggplot_mu2_",i), plot(saved_fits[[i]], "mu", prior = TRUE, plot_type = "ggplot"))

    if(i %in% effect){
      vdiffr::expect_doppelganger(paste0("ggplot_mu3_",i), plot(saved_fits[[i]], "mu", conditional = TRUE, plot_type = "ggplot"))
      vdiffr::expect_doppelganger(paste0("ggplot_mu4_",i), plot(saved_fits[[i]], "mu", conditional = TRUE, prior = TRUE, plot_type = "ggplot"))
    }else{
      expect_error(plot(saved_fits[[i]], "mu", conditional = TRUE, plot_type = "ggplot"),
                   "The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the effect. Please, verify that you specified at least one model assuming the presence of the effect.")
    }
  }

  # default base plot
  i <- 1
  vdiffr::expect_doppelganger(paste0("plot_mu1_",i), function()plot(saved_fits[[i]], "mu"))
  vdiffr::expect_doppelganger(paste0("plot_mu2_",i), function()plot(saved_fits[[i]], "mu", prior = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_mu3_",i), function()plot(saved_fits[[i]], "mu", conditional = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_mu4_",i), function()plot(saved_fits[[i]], "mu", conditional = TRUE, prior = TRUE))

  # additional settings
  vdiffr::expect_doppelganger(paste0("plot_mu5_",i), function()plot(saved_fits[[i]], "mu", prior = TRUE, dots_prior = list(col = "blue", lty = 2), col = "red", lty = 2, xlim = c(0, 1), main = "Title"))

  # transformation
  vdiffr::expect_doppelganger(paste0("plot_mu6_",i), function()plot(saved_fits[[i]], "mu", output_scale = "fishers_z"))
  vdiffr::expect_doppelganger(paste0("plot_mu7_",i), function()plot(saved_fits[[i]], "mu", output_scale = "r", prior = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_mu8_",i), plot(saved_fits[[i]], "mu", output_scale = "logOR", prior = TRUE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger(paste0("plot_mu10_",i), function()plot(saved_fits[[i]], "mu", output_scale = "OR"))

  ### heterogeneity
  # default ggplot2
  for(i in 1:length(saved_fits)){

    vdiffr::expect_doppelganger(paste0("ggplot_tau1_",i), plot(saved_fits[[i]], "tau", plot_type = "ggplot"))
    vdiffr::expect_doppelganger(paste0("ggplot_tau2_",i), plot(saved_fits[[i]], "tau", prior = TRUE, plot_type = "ggplot"))

    if(i %in% heterogeneity){
      vdiffr::expect_doppelganger(paste0("ggplot_tau3_",i), plot(saved_fits[[i]], "tau", conditional = TRUE, plot_type = "ggplot"))
      vdiffr::expect_doppelganger(paste0("ggplot_tau4_",i), plot(saved_fits[[i]], "tau", conditional = TRUE, prior = TRUE, plot_type = "ggplot"))
    }else{
      expect_error(plot(saved_fits[[i]], "tau", conditional = TRUE, plot_type = "ggplot"),
                   "The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the heterogeneity. Please, verify that you specified at least one model assuming the presence of the heterogeneity.")
    }
  }

  # default base plot
  i <- 1
  vdiffr::expect_doppelganger(paste0("plot_tau1_",i), function()plot(saved_fits[[i]], "tau"))
  vdiffr::expect_doppelganger(paste0("plot_tau2_",i), function()plot(saved_fits[[i]], "tau", prior = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_tau3_",i), function()plot(saved_fits[[i]], "tau", conditional = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_tau4_",i), function()plot(saved_fits[[i]], "tau", conditional = TRUE, prior = TRUE))

  # transformation
  vdiffr::expect_doppelganger(paste0("plot_tau5_",i), function()plot(saved_fits[[i]], "tau", output_scale = "r"))
  vdiffr::expect_doppelganger(paste0("plot_tau6_",i), function()plot(saved_fits[[i]], "tau", output_scale = "logOR", prior = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_tau7_",i), function()plot(saved_fits[[i]], "tau", output_scale = "OR", prior = TRUE))

  ### weightfunctions
  # default ggplot2
  for(i in 1:length(saved_fits)){
    set.seed(1)
    if(i %in% no_weightfunctions){
      expect_error(plot(saved_fits[[i]], "omega", plot_type = "ggplot"),
                   "The ensemble does not contain any posterior samples model-averaged across the selection models publication bias adjustment. Please, verify that you specified at least one selection models publication bias adjustment.")
      next
    }

    vdiffr::expect_doppelganger(paste0("ggplot_omega1_",i), plot(saved_fits[[i]], "omega", plot_type = "ggplot"))
    vdiffr::expect_doppelganger(paste0("ggplot_omega2_",i), plot(saved_fits[[i]], "omega", prior = TRUE, plot_type = "ggplot"))

    if(i %in% weightfunctions){
      vdiffr::expect_doppelganger(paste0("ggplot_omega3_",i), plot(saved_fits[[i]], "omega", conditional = TRUE, plot_type = "ggplot"))
      vdiffr::expect_doppelganger(paste0("ggplot_omega4_",i), plot(saved_fits[[i]], "omega", conditional = TRUE, prior = TRUE, plot_type = "ggplot"))
    }else{
      expect_error(plot(saved_fits[[i]], "omega", conditional = TRUE, plot_type = "ggplot"),
                   "The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of selection models publication bias adjustment. Please, verify that you specified at least one model assuming the presence of selection models publication bias adjustment.")
    }
  }

  # default base plot
  i <- 1
  set.seed(1)
  vdiffr::expect_doppelganger(paste0("plot_omega1_",i), function()plot(saved_fits[[i]], "omega"))
  vdiffr::expect_doppelganger(paste0("plot_omega2_",i), function()plot(saved_fits[[i]], "omega", prior = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_omega3_",i), function()plot(saved_fits[[i]], "omega", conditional = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_omega4_",i), function()plot(saved_fits[[i]], "omega", conditional = TRUE, prior = TRUE))

  # additional settings
  vdiffr::expect_doppelganger(paste0("plot_omega5_",i), function()plot(saved_fits[[i]], "omega", prior = TRUE, dots_prior = list(col = "blue", lty = 2, col.fill = "orange"), col = "red", lty = 2, col.fill = "red", rescale_x = TRUE))


  ### PET and PEESE
  for(i in 1:length(saved_fits)){
    set.seed(1)
    if(i %in% no_PET){
      expect_error(plot(saved_fits[[i]], "PET", plot_type = "ggplot"),
                   "The ensemble does not contain any posterior samples model-averaged across the PET. Please, verify that you specified at least one model for the PET")
      next
    }

    vdiffr::expect_doppelganger(paste0("ggplot_PET1_",i), plot(saved_fits[[i]], "PET", plot_type = "ggplot"))
    vdiffr::expect_doppelganger(paste0("ggplot_PET2_",i), plot(saved_fits[[i]], "PET", prior = TRUE, plot_type = "ggplot"))

    # PET is specified as null
    if(i %in% c(no_PET, 9)){
      expect_error(plot(saved_fits[[i]], "PET", conditional = TRUE, plot_type = "ggplot"),
                   "The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the PET models. Please, verify that you specified at least one model assuming the presence of the PET models.")
    }else{
      vdiffr::expect_doppelganger(paste0("ggplot_PET3_",i), plot(saved_fits[[i]], "PET", conditional = TRUE, plot_type = "ggplot"))
      vdiffr::expect_doppelganger(paste0("ggplot_PET4_",i), plot(saved_fits[[i]], "PET", conditional = TRUE, prior = TRUE, plot_type = "ggplot"))
    }
  }

  # default base plot
  i <- 1
  set.seed(1)
  vdiffr::expect_doppelganger(paste0("plot_PET1_",i), function()plot(saved_fits[[i]], "PET"))
  vdiffr::expect_doppelganger(paste0("plot_PET2_",i), function()plot(saved_fits[[i]], "PET", prior = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_PET3_",i), function()plot(saved_fits[[i]], "PET", conditional = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_PET4_",i), function()plot(saved_fits[[i]], "PET", conditional = TRUE, prior = TRUE))

  # additional settings
  vdiffr::expect_doppelganger(paste0("plot_PET5_",i), function()plot(saved_fits[[i]], "PET", prior = TRUE, dots_prior = list(col = "blue", lty = 2, col.fill = "orange"), col = "red", lty = 2, col.fill = "red", rescale_x = TRUE))

  for(i in 1:length(saved_fits)){
    set.seed(1)
    if(i %in% no_PEESE){
      expect_error(plot(saved_fits[[i]], "PEESE", plot_type = "ggplot"),
                   "The ensemble does not contain any posterior samples model-averaged across the PEESE. Please, verify that you specified at least one model for the PEESE.")
      next
    }

    vdiffr::expect_doppelganger(paste0("ggplot_PEESE1_",i), plot(saved_fits[[i]], "PEESE", plot_type = "ggplot"))
    vdiffr::expect_doppelganger(paste0("ggplot_PEESE2_",i), plot(saved_fits[[i]], "PEESE", prior = TRUE, plot_type = "ggplot"))

    if(i %in% no_PEESE){
      expect_error(plot(saved_fits[[i]], "PEESE", conditional = TRUE, plot_type = "ggplot"),
                   "The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the PEESE models. Please, verify that you specified at least one model assuming the presence of the PEESE models.")
    }else{
      vdiffr::expect_doppelganger(paste0("ggplot_PEESE3_",i), plot(saved_fits[[i]], "PEESE", conditional = TRUE, plot_type = "ggplot"))
      vdiffr::expect_doppelganger(paste0("ggplot_PEESE4_",i), plot(saved_fits[[i]], "PEESE", conditional = TRUE, prior = TRUE, plot_type = "ggplot"))
    }
  }

  # default base plot
  i <- 1
  set.seed(1)
  vdiffr::expect_doppelganger(paste0("plot_PEESE1_",i), function()plot(saved_fits[[i]], "PEESE"))
  vdiffr::expect_doppelganger(paste0("plot_PEESE2_",i), function()plot(saved_fits[[i]], "PEESE", prior = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_PEESE3_",i), function()plot(saved_fits[[i]], "PEESE", conditional = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_PEESE4_",i), function()plot(saved_fits[[i]], "PEESE", conditional = TRUE, prior = TRUE))

  # additional settings
  vdiffr::expect_doppelganger(paste0("plot_PEESE5_",i), function()plot(saved_fits[[i]], "PEESE", prior = TRUE, dots_prior = list(col = "blue", lty = 2, col.fill = "orange"), col = "red", lty = 2, col.fill = "red", rescale_x = TRUE))


  ### PET-PEESE (slow)
  # default ggplot2
  for(i in c(1, 2, 3, 11)){
    set.seed(1)
    if(i %in% no_PETPEESE){
      expect_error(plot(saved_fits[[i]], "PETPEESE", plot_type = "ggplot"),
                   "The ensemble does not contain any posterior samples model-averaged across the PET-PEESE publication bias adjustment. Please, verify that you specified at least one PET-PEESE publication bias adjustment.")
      next
    }

    vdiffr::expect_doppelganger(paste0("ggplot_PETPEESE1_",i), plot(saved_fits[[i]], "PETPEESE", plot_type = "ggplot"))
    vdiffr::expect_doppelganger(paste0("ggplot_PETPEESE2_",i), plot(saved_fits[[i]], "PETPEESE", prior = TRUE, plot_type = "ggplot"))

    if(i %in% PETPEESE){
      vdiffr::expect_doppelganger(paste0("ggplot_PETPEESE3_",i), plot(saved_fits[[i]], "PETPEESE", conditional = TRUE, plot_type = "ggplot"))
      vdiffr::expect_doppelganger(paste0("ggplot_PETPEESE4_",i), plot(saved_fits[[i]], "PETPEESE", conditional = TRUE, prior = TRUE, plot_type = "ggplot"))
    }else{
      expect_error(plot(saved_fits[[i]], "PETPEESE", conditional = TRUE, plot_type = "ggplot"),
                   "The ensemble does not contain any posterior samples model-averaged across the PET-PEESE publication bias adjustment. Please, verify that you specified at least one PET-PEESE publication bias adjustment.")
    }
  }

  # default base plot
  i <- 1
  set.seed(1)
  vdiffr::expect_doppelganger(paste0("plot_PETPEESE1_",i), function()plot(saved_fits[[i]], "PETPEESE"))
  vdiffr::expect_doppelganger(paste0("plot_PETPEESE2_",i), function()plot(saved_fits[[i]], "PETPEESE", prior = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_PETPEESE3_",i), function()plot(saved_fits[[i]], "PETPEESE", conditional = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_PETPEESE4_",i), function()plot(saved_fits[[i]], "PETPEESE", conditional = TRUE, prior = TRUE))

  # transformation
  vdiffr::expect_doppelganger(paste0("plot_PETPEESE5_",i), function()plot(saved_fits[[i]], "PETPEESE", output_scale = "logOR"))
  vdiffr::expect_doppelganger(paste0("plot_PETPEESE6_",i), function()plot(saved_fits[[i]], "PETPEESE", output_scale = "r", prior = TRUE))

  ### 3-level structure
  vdiffr::expect_doppelganger(paste0("plot_rho_",13),  function()plot(saved_fits[[13]], "rho"))
  vdiffr::expect_doppelganger(paste0("plot_rho2_",13), function()plot(saved_fits[[13]], "rho"))

  ### meta-regression parameter plots
  i <- 14
  set.seed(1)

  # factors
  vdiffr::expect_doppelganger(paste0("ggplot_reg-fac-1_",i), plot(saved_fits[[i]], "mod_cat", plot_type = "ggplot"))
  vdiffr::expect_doppelganger(paste0("ggplot_reg-fac-2_",i), plot(saved_fits[[i]], "mod_cat", prior = TRUE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger(paste0("ggplot_reg-fac-3_",i), plot(saved_fits[[i]], "mod_cat", conditional = TRUE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger(paste0("ggplot_reg-fac-4_",i), plot(saved_fits[[i]], "mod_cat", conditional = TRUE, prior = TRUE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger(paste0("plot_reg-fac-1_",i), function()plot(saved_fits[[i]], "mod_cat"))
  vdiffr::expect_doppelganger(paste0("plot_reg-fac-2_",i), function()plot(saved_fits[[i]], "mod_cat", prior = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_reg-fac-3_",i), function()plot(saved_fits[[i]], "mod_cat", conditional = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_reg-fac-4_",i), function()plot(saved_fits[[i]], "mod_cat", conditional = TRUE, prior = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_reg-fac-5_",i), function()plot(saved_fits[[i]], "mod_cat", output_scale = "logOR"))
  vdiffr::expect_doppelganger(paste0("plot_reg-fac-6_",i), function()plot(saved_fits[[i]], "mod_cat", output_scale = "r", prior = TRUE))

  # continuous
  vdiffr::expect_doppelganger(paste0("ggplot_reg-con-1_",i), plot(saved_fits[[i]], "mod_con", plot_type = "ggplot"))
  vdiffr::expect_doppelganger(paste0("ggplot_reg-con-2_",i), plot(saved_fits[[i]], "mod_con", prior = TRUE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger(paste0("ggplot_reg-con-3_",i), plot(saved_fits[[i]], "mod_con", conditional = TRUE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger(paste0("ggplot_reg-con-4_",i), plot(saved_fits[[i]], "mod_con", conditional = TRUE, prior = TRUE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger(paste0("plot_reg-con-1_",i), function()plot(saved_fits[[i]], "mod_con"))
  vdiffr::expect_doppelganger(paste0("plot_reg-con-2_",i), function()plot(saved_fits[[i]], "mod_con", prior = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_reg-con-3_",i), function()plot(saved_fits[[i]], "mod_con", conditional = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_reg-con-4_",i), function()plot(saved_fits[[i]], "mod_con", conditional = TRUE, prior = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_reg-con-5_",i), function()plot(saved_fits[[i]], "mod_con", output_scale = "logOR"))
  vdiffr::expect_doppelganger(paste0("plot_reg-con-6_",i), function()plot(saved_fits[[i]], "mod_con", output_scale = "r", prior = TRUE))

  # continuous, alternative only
  i <- 15
  vdiffr::expect_doppelganger(paste0("ggplot_reg-con-alt-1_",i), plot(saved_fits[[i]], "mod_con", plot_type = "ggplot"))
  vdiffr::expect_doppelganger(paste0("ggplot_reg-con-alt-2_",i), plot(saved_fits[[i]], "mod_con", prior = TRUE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger(paste0("ggplot_reg-con-alt-3_",i), plot(saved_fits[[i]], "mod_con", conditional = TRUE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger(paste0("ggplot_reg-con-alt-4_",i), plot(saved_fits[[i]], "mod_con", conditional = TRUE, prior = TRUE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger(paste0("plot_reg-con-alt-1_",i), function()plot(saved_fits[[i]], "mod_con"))
  vdiffr::expect_doppelganger(paste0("plot_reg-con-alt-2_",i), function()plot(saved_fits[[i]], "mod_con", prior = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_reg-con-alt-3_",i), function()plot(saved_fits[[i]], "mod_con", conditional = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_reg-con-alt-4_",i), function()plot(saved_fits[[i]], "mod_con", conditional = TRUE, prior = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_reg-con-alt-5_",i), function()plot(saved_fits[[i]], "mod_con", output_scale = "logOR"))
  vdiffr::expect_doppelganger(paste0("plot_reg-con-alt-6_",i), function()plot(saved_fits[[i]], "mod_con", output_scale = "r", prior = TRUE))

  plot(saved_fits[[14]], "mod_cat", conditional = TRUE)
})


test_that("Individual model plots work", {

  # default ggplot2
  for(i in 1:length(saved_fits)){

    if(saved_fits[[i]]$add_info[["algorithm"]] == "ss"){
      expect_error(plot_models(saved_fits[[i]], plot_type = "ggplot"),
                   "The estimated model using the spike and slab style model-averaging")
      next
    }

    vdiffr::expect_doppelganger(paste0("ggplot_models1_",i), plot_models(saved_fits[[i]], plot_type = "ggplot"))

    if(i %in% effect){
      vdiffr::expect_doppelganger(paste0("ggplot_models2_",i), plot_models(saved_fits[[i]], conditional = TRUE, plot_type = "ggplot"))
    }else{
      expect_error(plot_models(saved_fits[[i]], conditional = TRUE, plot_type = "ggplot"),
                   "The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the effect. Please, verify that you specified at least one model assuming the presence of the effect.")
    }
  }

  # default base plot
  i <- 1
  vdiffr::expect_doppelganger(paste0("plot_models1_",i), function()plot_models(saved_fits[[i]]))
  vdiffr::expect_doppelganger(paste0("plot_models2_",i), function()plot_models(saved_fits[[i]], conditional = TRUE))

  # different output scale
  vdiffr::expect_doppelganger(paste0("plot_models3_",i), function()plot_models(saved_fits[[i]], output_scale = "fishers_z"))
  vdiffr::expect_doppelganger(paste0("plot_models4_",i), function()plot_models(saved_fits[[i]], output_scale = "fishers_z", conditional = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_models4-1_",i), function()plot_models(saved_fits[[i]], output_scale = "logOR", conditional = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_models4-2_",i), function()plot_models(saved_fits[[i]], output_scale = "OR", conditional = TRUE))

  # different ordering
  vdiffr::expect_doppelganger(paste0("plot_models5_",i), function()plot_models(saved_fits[[i]], order = "increasing", order_by = "estimate"))
  vdiffr::expect_doppelganger(paste0("plot_models6_",i), function()plot_models(saved_fits[[i]], order = "decreasing", order_by = "BF"))
  vdiffr::expect_doppelganger(paste0("plot_models7_",i), function()plot_models(saved_fits[[i]], order = "increasing", order_by = "probability"))

  # check tau parameter
  vdiffr::expect_doppelganger(paste0("plot_models1_tau_",i), function()plot_models(saved_fits[[i]], parameter = "tau"))
  vdiffr::expect_doppelganger(paste0("plot_models2_tau_",i), function()plot_models(saved_fits[[i]], parameter = "tau", conditional = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_models3_tau_",i), function()plot_models(saved_fits[[i]], parameter = "tau", output_scale = "logOR", conditional = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_models4_tau_",i), function()plot_models(saved_fits[[i]], parameter = "tau", output_scale = "OR", conditional = TRUE))


})


test_that("Forest plots work", {

  # default ggplot2
  for(i in 1:length(saved_fits)){

    vdiffr::expect_doppelganger(paste0("ggplot_forest1_",i), forest(saved_fits[[i]], plot_type = "ggplot"))

    if(i %in% effect){
      vdiffr::expect_doppelganger(paste0("ggplot_forest2_",i), forest(saved_fits[[i]], conditional = TRUE, plot_type = "ggplot"))
    }else{
      expect_error(forest(saved_fits[[i]], conditional = TRUE, plot_type = "ggplot"),
                   "The ensemble does not contain any posterior samples model-averaged across the models assuming the presence of the effect. Please, verify that you specified at least one model assuming the presence of the effect.")
    }
  }

  # default base plot
  i <- 1
  vdiffr::expect_doppelganger(paste0("plot_forest1_",i), function()forest(saved_fits[[i]]))
  vdiffr::expect_doppelganger(paste0("plot_forest2_",i), function()forest(saved_fits[[i]], conditional = TRUE))

  # different output scale
  vdiffr::expect_doppelganger(paste0("plot_forest3_",i), function()forest(saved_fits[[i]], output_scale = "fishers_z"))
  vdiffr::expect_doppelganger(paste0("plot_forest4_",i), function()forest(saved_fits[[i]], output_scale = "fishers_z", conditional = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_forest4-1_",i), function()forest(saved_fits[[i]], output_scale = "logOR", conditional = TRUE))
  vdiffr::expect_doppelganger(paste0("plot_forest4-2_",i), function()forest(saved_fits[[i]], output_scale = "OR", conditional = TRUE))

  # different ordering
  vdiffr::expect_doppelganger(paste0("plot_forest5_",i), function()forest(saved_fits[[i]], order = "increasing"))
  vdiffr::expect_doppelganger(paste0("plot_forest6_",i), function()forest(saved_fits[[i]], order = "decreasing"))
  saved_fits[[i]] <- update(saved_fits[[i]], study_names = c("a", "c", "b"))
  vdiffr::expect_doppelganger(paste0("plot_forest7_",i), function()forest(saved_fits[[i]], order = "alphabetical"))

})


test_that("Marginal posterior plots work", {


  expect_error(marginal_plot(saved_fits[[1]]), "'marginal_plot' function is available only for RoBMA regression models")
  expect_error(marginal_plot(saved_fits[[14]], "mu"), "The 'mu' values are not recognized by the 'parameter' argument.")

  vdiffr::expect_doppelganger("mm_ggplot_mod_cat_1", marginal_plot(saved_fits[[14]], "mod_cat", plot_type = "ggplot"))
  vdiffr::expect_doppelganger("mm_ggplot_mod_cat_2", marginal_plot(saved_fits[[14]], "mod_cat", prior = TRUE, plot_type = "ggplot"))
  vdiffr::expect_doppelganger("mm_ggplot_mod_cat_3", marginal_plot(saved_fits[[14]], "mod_cat", prior = TRUE, plot_type = "ggplot", output_scale = "r"))
  vdiffr::expect_doppelganger("mm_ggplot_mod_con_1", marginal_plot(saved_fits[[14]], "mod_con", prior = TRUE, plot_type = "ggplot", xlim = c(-1, 1)))
  vdiffr::expect_doppelganger("mm_ggplot_mod_con_2", function()marginal_plot(saved_fits[[15]], "mod_con", conditional = TRUE))
})
