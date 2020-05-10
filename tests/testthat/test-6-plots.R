context("(6) plot functions")
skip_on_cran()


# test objects - assuming that the fit function worked properly
saved_fits    <- readRDS(file = "saved_fits.RDS")
set.seed(1)

test_that("Parameter plots works", {

  for(i in 1:length(saved_fits)){

    expect_doppelganger(paste0("plot_mu1_",i),plot(saved_fits[[i]], "mu", plot_type = "ggplot"))
    expect_doppelganger(paste0("plot_mu2_",i),plot(saved_fits[[i]], "mu", prior = TRUE, plot_type = "ggplot"))

    if(i %in% c(2, 4, 5, 6, 7, 8)){
      expect_error(plot(saved_fits[[i]], "mu", type = "conditional", plot_type = "ggplot"),
                   "The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameter or models with the parameter did not converge.")
      expect_error(plot(saved_fits[[i]], "mu", type = "conditional", prior = TRUE, plot_type = "ggplot"),
                   "The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameter or models with the parameter did not converge.")
    }else{
      expect_doppelganger(paste0("plot_mu3_",i),plot(saved_fits[[i]], "mu", type = "conditional", plot_type = "ggplot"))
      expect_doppelganger(paste0("plot_mu4_",i),plot(saved_fits[[i]], "mu", type = "conditional", prior = TRUE, plot_type = "ggplot"))
    }

  }

  for(i in 1:length(saved_fits)){

    expect_doppelganger(paste0("plot_tau1_",i),plot(saved_fits[[i]], "tau", plot_type = "ggplot"))
    expect_doppelganger(paste0("plot_tau2_",i),plot(saved_fits[[i]], "tau", prior = TRUE, plot_type = "ggplot"))

    if(i %in% c(2, 4, 5, 7, 8)){
      expect_error(plot(saved_fits[[i]], "tau", type = "conditional", plot_type = "ggplot"),
                   "The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameter or models with the parameter did not converge.")
      expect_error(plot(saved_fits[[i]], "tau", type = "conditional", prior = TRUE, plot_type = "ggplot"),
                   "The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameter or models with the parameter did not converge.")
    }else{
      expect_doppelganger(paste0("plot_tau3_",i),plot(saved_fits[[i]], "tau", type = "conditional", plot_type = "ggplot"))
      expect_doppelganger(paste0("plot_tau4_",i),plot(saved_fits[[i]], "tau", type = "conditional", prior = TRUE, plot_type = "ggplot"))
    }

  }

  for(i in 1:length(saved_fits)){

    if(i %in% c(2, 5, 6, 7, 8)){

      expect_error(plot(saved_fits[[i]], "omega", plot_type = "ggplot"),
                   "The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameter or models with the parameter did not converge.")
      expect_error(plot(saved_fits[[i]], "omega", prior = TRUE, plot_type = "ggplot"),
                   "The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameter or models with the parameter did not converge.")
      expect_error(plot(saved_fits[[i]], "omega", type = "conditional", plot_type = "ggplot"),
                   "The ensemble cointains no non-null model adjusting for publication bias.")
      expect_error(plot(saved_fits[[i]], "omega", type = "conditional", prior = TRUE, plot_type = "ggplot"),
                   "The ensemble cointains no non-null model adjusting for publication bias.")

    }else{

      expect_doppelganger(paste0("plot_weight_function1_",i),plot(saved_fits[[i]], "omega", plot_type = "ggplot"))
      expect_doppelganger(paste0("plot_weight_function2_",i),plot(saved_fits[[i]], "omega", prior = TRUE, plot_type = "ggplot"))
      expect_doppelganger(paste0("plot_weight_function3_",i),plot(saved_fits[[i]], "omega", type = "conditional", plot_type = "ggplot"))
      expect_doppelganger(paste0("plot_weight_function4_",i),plot(saved_fits[[i]], "omega", type = "conditional", prior = TRUE, plot_type = "ggplot"))

    }
  }

  for(i in 1:length(saved_fits)){

    if(i %in% c(2, 5, 6, 7, 8)){
      expect_error(plot(saved_fits[[i]], "omega", plot_type = "ggplot", weights = TRUE),
                   "The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameter or models with the parameter did not converge.")
      expect_error(plot(saved_fits[[i]], "omega", prior = TRUE, plot_type = "ggplot", weights = TRUE),
                   "The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameter or models with the parameter did not converge.")
      expect_error(plot(saved_fits[[i]], "omega", type = "conditional", plot_type = "ggplot", weights = TRUE),
                   "The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameter or models with the parameter did not converge.")
      expect_error(plot(saved_fits[[i]], "omega", type = "conditional", prior = TRUE, plot_type = "ggplot", weights = TRUE),
                   "The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameter or models with the parameter did not converge.")
    }else{
      temp_plots1 <- plot(saved_fits[[i]], "omega", plot_type = "ggplot", weights = TRUE)
      temp_plots2 <- plot(saved_fits[[i]], "omega", prior = TRUE, plot_type = "ggplot", weights = TRUE)
      temp_plots3 <- plot(saved_fits[[i]], "omega", type = "conditional", plot_type = "ggplot", weights = TRUE)
      temp_plots4 <- plot(saved_fits[[i]], "omega", type = "conditional", prior = TRUE, plot_type = "ggplot", weights = TRUE)
    }

    for(j in 1:length(temp_plots1)){

      expect_doppelganger(paste0("plot_omega1_",i,"_j_",j), temp_plots1[[j]])
      expect_doppelganger(paste0("plot_omega2_",i,"_j_",j), temp_plots2[[j]])

      if(i %in% c(2, 5, 6, 7, 8)){
        next
      }else{
        expect_doppelganger(paste0("plot_omega3_",i,"_j_",j), temp_plots3[[j]])
        expect_doppelganger(paste0("plot_omega4_",i,"_j_",j), temp_plots4[[j]])
      }
    }

  }

  for(i in 1:length(saved_fits)){

    expect_doppelganger(paste0("plot_theta1_",i),plot(saved_fits[[i]],    "theta",  plot_type = "ggplot"))
    expect_doppelganger(paste0("plot_forest1_",i),plot(saved_fits[[i]],   "forest", plot_type = "ggplot"))
    expect_doppelganger(paste0("plot_combined1_",i),plot(saved_fits[[i]], c("theta", "forest"), plot_type = "ggplot"))

    if(i %in% c(2, 4, 5, 6, 7, 8)){
      expect_error(plot(saved_fits[[i]],    "theta",  type = "conditional", plot_type = "ggplot"),"The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameters or all of the models did not converge.")
      expect_error(plot(saved_fits[[i]],   "forest", type = "conditional", plot_type = "ggplot"),"The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameters or all of the models did not converge.")
      expect_error(plot(saved_fits[[i]], c("theta", "forest"), type = "conditional", plot_type = "ggplot"),"The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameters or all of the models did not converge.")
    }else{
      expect_doppelganger(paste0("plot_theta2_",i),plot(saved_fits[[i]],    "theta",  type = "conditional", plot_type = "ggplot"))
      expect_doppelganger(paste0("plot_forest2_",i),plot(saved_fits[[i]],   "forest", type = "conditional", plot_type = "ggplot"))
      expect_doppelganger(paste0("plot_combined2_",i),plot(saved_fits[[i]], c("theta", "forest"), type = "conditional", plot_type = "ggplot"))
    }

  }

})


test_that("Individual model plots works", {

  for(i in 1:length(saved_fits)){

    expect_doppelganger(paste0("plot_individual_mu1_",i),plot(saved_fits[[i]], "mu", type = "individual", plot_type = "ggplot"))
    expect_doppelganger(paste0("plot_individual_mu2_",i),plot(saved_fits[[i]], "mu", type = "individual", order = "prob", plot_type = "ggplot"))

    if(i %in% c(2, 4, 5, 6)){
      expect_error(plot(saved_fits[[i]], "mu", type = c("individual", "conditional"), plot_type = "ggplot"),"The ensemble contains no non-null model with the specified parameter.")
    }else if(i %in% c(7, 8)){
      expect_error(plot(saved_fits[[i]], "mu", type = c("individual", "conditional"), plot_type = "ggplot"),"The parameter could not be plotted because it is not in the ensemble. Possible cause might be trying to plot a parameter from an ensemble where either no model has the parameters or all of the models did not converge.")
    }else{
      expect_doppelganger(paste0("plot_individual_mu3_",i),plot(saved_fits[[i]], "mu", type = c("individual", "conditional"), plot_type = "ggplot"))
    }
  }


  for(i in 1:length(saved_fits)){

    expect_doppelganger(paste0("plot_individual_tau1_",i),plot(saved_fits[[i]], "tau", type = "individual", plot_type = "ggplot"))
    expect_doppelganger(paste0("plot_individual_tau2_",i),plot(saved_fits[[i]], "tau", type = "individual", order = "prob", plot_type = "ggplot"))

    if(i %in% c(2, 4, 5, 7, 8)){
      expect_error(plot(saved_fits[[i]], "tau", type = c("individual", "conditional"), plot_type = "ggplot"),"The ensemble contains no non-null model with the specified parameter.")
    }else{
      expect_doppelganger(paste0("plot_individual_tau3_",i),plot(saved_fits[[i]], "tau", type = c("individual", "conditional"), plot_type = "ggplot"))
    }

  }


  for(i in 1:length(saved_fits)){


    if(i %in% c(2, 5, 6, 7, 8)){
        expect_error(plot(saved_fits[[i]], "omega", type = "individual", plot_type = "ggplot"),
                     "The ensemble cointains no non-null model adjusting for publication bias.")
        expect_error(plot(saved_fits[[i]], "omega", type = "individual", order = "prob", plot_type = "ggplot"),
                     "The ensemble cointains no non-null model adjusting for publication bias.")
        expect_error(plot(saved_fits[[i]], "omega", type = c("individual", "conditional"), plot_type = "ggplot"),
                     "The ensemble cointains no non-null model adjusting for publication bias.")
      }else{
        temp_plots1 <- plot(saved_fits[[i]], "omega", type = "individual", plot_type = "ggplot")
        temp_plots2 <- plot(saved_fits[[i]], "omega", type = "individual", order = "prob", plot_type = "ggplot")
        temp_plots3 <- plot(saved_fits[[i]], "omega", type = c("individual", "conditional"), plot_type = "ggplot")
      }

      for(j in 1:length(temp_plots1)){

        expect_doppelganger(paste0("plot_individual_omega1_",i,"_j_",j), temp_plots1[[j]])
        expect_doppelganger(paste0("plot_individual_omega2_",i,"_j_",j), temp_plots2[[j]])

        if(i %in% c(2, 5, 6, 7, 8)){
          next
        }else{
          expect_doppelganger(paste0("plot_individual_omega3_",i,"_j_",j), temp_plots3[[j]])
        }
      }

    }

})


