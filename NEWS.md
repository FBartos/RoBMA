## version 3.5.1
### Features
- `summary()` function now supports a `standardized_coefficients` argument to report either standardized (default) or raw meta-regression coefficients
- `extract()` function extracts the posterior samples of the model parameters
- `true_effects()` function to summarize the true effect size estimates
- `predict()` method for `RoBMA()` and `RoBMA.reg()` models when using the `algorithm = "ss"`

### Fixes
- fitting a meta-regression using predictors with missing values result in a clear error message

### Changes
- improving the speed of unit tests

## version 3.5
### Features
- approximate and computationally feasibly 3lvl selection models via the `RoBMA()` and `RoBMA.reg()` functions with the `study_ids` argument when using `algorithm = "ss"`
- 3lvl binomial-normal models for binary data via the `BiBMA` and `BiBMA.reg` functions with the `study_ids` argument when using `algorithm = "ss"`
- `pooled_effect()` function to compute the pooled effect size from the `RoBMA.reg`, `NoBMA.reg`, and `BiBMA.reg` models
- `adjusted_effect()` function to compute the adjusted effect size from the `RoBMA.reg`, `NoBMA.reg`, and `BiBMA.reg` models
- enables `summary_heterogeneity()` for BiBMA models

### Fixes
- passing and checks of the `study_ids` and `study_labels` arguments
- PEESE prior distribution now scale as 1/scale instead of 1/scale^2 with the `rescale_priors` argument  
- the conditional prediction interval based on `summary_heterogeneity()` is now conditional on the presence of the effect
- additional minor prior handling fixes (i.e., missing marginal estimates when only alternative prior distributions were specified etc)
- diagnostics with mixture baseline priors when using `algorithm = "ss"`
- `summary_heterogeneity()` with only a single study does not produce relative heterogeneity instead of crashing

## version 3.4
### Features
- adding binomial-normal meta-regression models for binary data via the `BiBMA.reg` function
- the spike and slab algorithm for faster model estimation via the `algorithm = "ss"` argument for BiBMA models
- default prior distributions for all parameters of BiBMA models are now set via the `set_default_binomial_priors()` function

## version 3.3
### Features
- the spike and slab algorithm for faster model estimation via the `algorithm = "ss"` argument (see a new vignette for more details)
- refactoring of the JAGS C++ code of weighted distributions and exporting of the lpdfs into JAGS (maintenance)
- weights_mix JAGS prior distribution to sample a mixture of weight functions directly

### Fixes
- incorrectly omitting models with more than one predictor when computing conditional marginal summary

## version 3.2.1
### Features
- default prior distributions for all parameters are now set via the `set_default_priors()` function
- `rescale_priors` argument allows to conveniently re-scale the prior distributions for the effect, heterogeneity, and bias simultaneously

## version 3.2
### Features
- `summary_heterogeneity()` function to summarize the heterogeneity of the RoBMA models (prediction interval, tau, tau^2, I^2, and H^2)
- `check_RoBMA_convergence()` function to check the convergence of the RoBMA models
- adds informed prior distributions for binary and time-to-event outcomes via BayesTools 0.2.17

### Fixes
- checking and fixing the number of available cores upon loading the package (hopefully fixes some parallelization issues)
- `update()` function re-evaluates convergence checks of individual models (https://github.com/FBartos/RoBMA/issues/34) 
- typos and minor issues in the vignettes


## version 3.1
### Features
- binomial-normal models for binary data via the `BiBMA` function
- `NoBMA` and `NoBMA.reg()` functions as wrappers around `RoBMA` `RoBMA.reg()` functions for simpler specification of publication bias unadjusted Bayesian model-averaged meta-analysis
- adding odds ratios output transformation` 
- extending (instead of a complete refitting) of models via the `update.RoBMA()` function (only non-converged models by default or all by setting `extend_all = TRUE`)

### Fixes
- handling of non-converged models

## version 3.0.1
### Fixes (thanks to Don & Rens)
- compilation issues with Clang (https://github.com/FBartos/RoBMA/issues/28)
- lapack path specifications (https://github.com/FBartos/RoBMA/issues/24)

## version 3.0
### Features
- meta-regression with `RoBMA.reg()` function
- posterior marginal summary and plots for the `RoBMA.reg` models with `summary_marginal()` and `plot_marginal()` functions
- new vignette on hierarchical Bayesian model-averaged meta-analysis
- new vignette on robust Bayesian model-averaged meta-regression
- adding vignette from AMPPS tutorial
- faster implementation of JAGS multivariate normal distribution (based on the BUGS JAGS module)
- incorporating `weight` argument in the `RoBMA` and `combine_data` functions in order to pass `custom` likelihood weights
- ability to use inverse square weights in the weighted meta-analysis by setting a `weighted_type = "inverse_sqrt"` argument 

### Changes
- reworked interface for the hierarchical models. Prior distributions are now specified via the `priors_hierarchical` and `priors_hierarchical_null` arguments instead of `priors_rho` and `priors_rho_null`. The model summary now shows `Hierarchical` component summary.

## version 2.3.2
### Fixes
- suppressing start-up message 
- cleaning up imports

## version 2.3.1
### Fixes
- fixing weighted meta-analysis parameterization 

## version 2.3
### Features
- weighted meta-analysis by specifying `study_ids` argument in `RoBMA()` and setting `weighted = TRUE`. The likelihood contribution of estimates from each study is down-weighted proportionally to the number of estimates in that study. Note that this experimental feature is supposed to provide a conservative alternative for estimating RoBMA in cases with multiple estimates from a study where the multivariate option is not computationally feasible.

## version 2.2.3
### Fixes
- updating the Makevars to install with R 4.2 and JAGS 4.3.1

## version 2.2.2
### Fixes
- updating the C++ to compile on M1 Mac

## version 2.2.1
### Changes
- message about the effect size scale of parameter estimates is always shown
- compatibility with BayesTools 0.2.0+

## version 2.2
### Features
- three-level meta-analysis by specifying `study_ids` argument in `RoBMA`. However, note that this is (1) an experimental feature and (2) the computational expense of fitting selection models with clustering is extreme. As of now, it is almost impossible to have more than 2-3 estimates clustered within a single study).

## version 2.1.2
### Fixes
- adding Windows ucrt patch (thanks to Tomas Kalibera)
- adding BayesTools version check

## version 2.1.1
### Fixes
- incorrectly formatted citations in vignettes and capitalization

### Features
- adding `informed_prior()` function (from the BayesTools package) that allows specification of various informed prior distributions from the field of medicine and psychology
- adding a vignette reproducing the example of dentine sensitivity with the informed Bayesian model-averaged meta-analysis from Bartoš et al., 2021 ([open-access](https://onlinelibrary.wiley.com/doi/10.1002/sim.9170)),
- further reductions of fitted object size when setting `save = "min"`

## version 2.1
### Fixes
- more informative error message when the JAGS module fails to load
- correcting wrong PEESE transformation for the individual models summaries (issue #12)
- fixing error message for missing conditional PET-PEESE
- fixing incorrect lower bound check for log(OR)

### Features
- adding `interpret()` function (issue #11)
- adding effect size transformation via `output_scale` argument to `plot()` and `plot_models()` functions
- better handling of effect size transformations and scaling - BayesTools style back-end functions with Jacobian transformations

## version 2.0
Please notice that this is a major release that breaks backwards compatibility.

### Changes
 - naming of the arguments specifying prior distributions for the different parameters/components of the models changed (`priors_mu` -> `priors_effect`, `priors_tau` -> `priors_heterogeneity`, and `priors_omega` -> `priors_bias`),
 - prior distributions for specifying weight functions now use a dedicated function (`prior(distribution = "two.sided", parameters = ...)` -> `prior_weightfunction(distribution = "two.sided", parameters = ...)`),
 - new dedicated function for specifying no publication bias adjustment component / no heterogeneity component (`prior_none()`),
 - new dedicated functions for specifying models with the PET and PEESE publication bias adjustments (`prior_PET(distribution = "Cauchy", parameters = ...)` and `prior_PEESE(distribution = "Cauchy", parameters = ...)`),
 - new default prior distribution specification for the publication bias adjustment part of the models (corresponding to the RoBMA-PSMA model from Bartoš et al., 2021 [manuscript](https://doi.org/10.1002/jrsm.1594)),
 - new `model_type` argument allowing to specify different "pre-canned" models (`"PSMA"` = RoBMA-PSMA, `"PP"` = RoBMA-PP, `"2w"` = corresponding to Maier et al., in press , [manuscript](https://doi.org/10.1037/met0000405)),
 - `combine_data` function allows combination of different effect sizes / variability measures into a common effect size measure (also used from within the `RoBMA` function),
 - better and improved automatic fitting procedure now enabled by default (can be turned of with `autofit = FALSE`)
 - prior distributions can be specified on the different scale than the supplied effect sizes (the package fits the model on Fisher's z scale and back transforms the results back to the scale that was used for prior distributions specification, Cohen's d by default, but both of them can be overwritten with the `prior_scale` and `transformation` arguments),
 - new prior distributions, e.g., beta or fixed weight functions,
 - estimates from individual models are now plotted with the `plot_models()` function and the forest plot can be obtained with the `forest()` function,
 - the posterior distribution plots for the individual weights are no able supported, however, the weightfunction and the PET-PEESE publication bias adjustments can be visualized with the `plot.RoBMA()` function and `parameter = "weightfunction"` and `parameter = "PET-PEESE"`.

## version 1.2.1
### Fixes
- check_setup function not working at all

## version 1.2.0
### Changes
- the studies's true effects are now marginalized out of the random effects models and are no longer estimated (see Appendix A of our [manuscript](https://doi.org/10.1037/met0000405) for more details). As a results, arguments referring to the true effects are now disabled.
- all models are now being estimated using the likelihood of effect sizes (instead of test-statistics as usually defined). We reproduced the simulation study that we used to evaluate the method performance and it achieved identical results (up to MCMC error, before marginalizing out the true effects). A big advantage of using the normal likelihood for effect sizes is a considerable speed up of the whole estimation process.
- as a results of these two changes, the results of the models will differ to those of pre 1.2.0 version

### Fixes
- autofit being turn on if any control argument was specified

## version 1.1.2
### Fixes
- vdiffr not being used conditionally in unit tests

## version 1.1.1
### Fixes
- inability to fit a model without specifying a seed
- inability to produce individual model plots due to incompatibility with the newer versions of ggplot2  

## version 1.1.0
### Features
- parallel within and between model fitting using the parallel package with 'parallel = TRUE' argument

## version 1.0.5
### Fixes:
- models being fitted automatically until reaching R-hat lower than 1.05 without setting max_rhat and autofit control parameters
- bug preventing to draw a bivariate plot of mu and tau
- range for parameter estimates from individual models no containing 0 (or 1 in case of OR measured effect sizes)
- inability to fit a model with only null mu distributions if correlation or OR measured effect sizes were specified
- ordering of the estimated and observed effects when both of them are requested simultaneously
- formatting of this file (NEWS.md)

### Improvements:
- priors plot: parameter specification, default plotting range, clearer x-axis labels in cases when the parameter is defined on transformed scale
- parameters plots: probability scale always ends at the same spot as is the last tick on the density scale
- adding warnings if any of the specified models has Rhat higher than 1.05 or the specified value
- grouping the same warnings messages together

## version 1.0.4
### Fixes:
- inability to run models without the silent = TRUE control

## version 1.0.3
### Features:
- x-axis rescaling for the weight function plot (by setting 'rescale_x = TRUE' in the 'plot.RoBMA' function)
- setting expected direction of the effect in for RoBMA function

### Fixes:
- marginal likelihood calculation for models with spike prior distribution on mean parameter which location was not set to 0
- some additional error messages 

### CRAM requested changes:
- changing information messages from 'cat' to 'message' from plot related functions
- saving and returning the 'par' settings to the user defined one in the base plot functions

## version 1.0.2
### Fixes:
- the summary and plot function now shows quantile based confidence intervals for individual models instead of the HPD provided before (this affects only 'summary'/'plot' with 'type = "individual"', all other confidence intervals were quantile based before)

## version 1.0.1
### Fixes:
- summary function returning median instead of mean

## version 1.0.0 (vs the osf version)
### Fixes:
- incorrectly weighted theta estimates
- models with non-zero point prior distribution incorrectly plotted using when "models" option in case that the mu parameter was transformed

### Additional features:
- analyzing OR
- distributions implemented using boost library (helps with convergence issues)
- ability to mute the non-suppressible "precision not achieved" warning messages by using "silent" = TRUE inside of the control argument
- vignettes

### Notable changes:
- the way how the seed is set before model fitting (the simulation study will not be reproducible with the new version of the package)
