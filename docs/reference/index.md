# Package index

## Package

- [`RoBMA-package`](https://fbartos.github.io/RoBMA/reference/RoBMA-package.md)
  [`RoBMA_package`](https://fbartos.github.io/RoBMA/reference/RoBMA-package.md)
  [`RoBMA.package`](https://fbartos.github.io/RoBMA/reference/RoBMA-package.md)
  : RoBMA: Robust Bayesian Meta-Analysis
- [`RoBMA.options()`](https://fbartos.github.io/RoBMA/reference/RoBMA_options.md)
  [`RoBMA.get_option()`](https://fbartos.github.io/RoBMA/reference/RoBMA_options.md)
  : Options for the RoBMA package
- [`set_autofit_control()`](https://fbartos.github.io/RoBMA/reference/RoBMA_control.md)
  [`set_convergence_checks()`](https://fbartos.github.io/RoBMA/reference/RoBMA_control.md)
  : Control MCMC fitting process

## Model fitting

- [`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md) :
  Bayesian Meta-Analysis
- [`brma.glmm()`](https://fbartos.github.io/RoBMA/reference/brma.glmm.md)
  : Bayesian Generalized Meta-Analysis
- [`RoBMA()`](https://fbartos.github.io/RoBMA/reference/RoBMA.md) :
  Robust Bayesian Model-Averaged Meta-Analysis
- [`BMA()`](https://fbartos.github.io/RoBMA/reference/BMA.md) : Bayesian
  Model-Averaged Meta-Analysis
- [`BMA.glmm()`](https://fbartos.github.io/RoBMA/reference/BMA.glmm.md)
  : Bayesian Model-Averaged Generalized Meta-Analysis
- [`bselmodel()`](https://fbartos.github.io/RoBMA/reference/bselmodel.md)
  : Bayesian Selection Model
- [`bPET()`](https://fbartos.github.io/RoBMA/reference/bPET.md) :
  Bayesian Precision-Effect Test (PET) Model
- [`bPEESE()`](https://fbartos.github.io/RoBMA/reference/bPEESE.md) :
  Bayesian Precision-Effect Estimate with Standard Errors (PEESE) Model
- [`update(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/update.brma.md)
  : Update a brma Fit

## Inputs, priors, and model specification

- [`data_input`](https://fbartos.github.io/RoBMA/reference/data_input.md)
  : Input Data Specification
- [`fitting_specification`](https://fbartos.github.io/RoBMA/reference/fitting_specification.md)
  : Fitting specification
- [`prior_specification`](https://fbartos.github.io/RoBMA/reference/prior_specification.md)
  : Prior specification
- [`RoBMA_prior_specification`](https://fbartos.github.io/RoBMA/reference/RoBMA_prior_specification.md)
  : Prior specification for model-averaging
- [`prior()`](https://fbartos.github.io/RoBMA/reference/prior.md) :
  Prior Distribution
- [`prior_none()`](https://fbartos.github.io/RoBMA/reference/prior_none.md)
  : Empty Prior
- [`prior_factor()`](https://fbartos.github.io/RoBMA/reference/prior_factor.md)
  : Factor Prior
- [`prior_informed()`](https://fbartos.github.io/RoBMA/reference/prior_informed.md)
  : Informed Prior
- [`prior_PET()`](https://fbartos.github.io/RoBMA/reference/prior_PET.md)
  : PET Prior
- [`prior_PEESE()`](https://fbartos.github.io/RoBMA/reference/prior_PEESE.md)
  : PEESE Prior
- [`prior_weightfunction()`](https://fbartos.github.io/RoBMA/reference/prior_weightfunction.md)
  [`wf_cumulative()`](https://fbartos.github.io/RoBMA/reference/prior_weightfunction.md)
  [`wf_fixed()`](https://fbartos.github.io/RoBMA/reference/prior_weightfunction.md)
  [`wf_independent()`](https://fbartos.github.io/RoBMA/reference/prior_weightfunction.md)
  : Weightfunction Prior
- [`estimate_unit_information_sd()`](https://fbartos.github.io/RoBMA/reference/estimate_unit_information_sd.md)
  : Estimate Unit Information Standard Deviation
- [`contr.orthonormal()`](https://fbartos.github.io/RoBMA/reference/contr.BayesTools.md)
  [`contr.meandif()`](https://fbartos.github.io/RoBMA/reference/contr.BayesTools.md)
  [`contr.independent()`](https://fbartos.github.io/RoBMA/reference/contr.BayesTools.md)
  : BayesTools Contrast Matrices

## Summaries and estimates

- [`summary(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/summary.brma.md)
  [`print(`*`<summary.brma>`*`)`](https://fbartos.github.io/RoBMA/reference/summary.brma.md)
  [`print(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/summary.brma.md)
  : Summarize brma Object
- [`interpret()`](https://fbartos.github.io/RoBMA/reference/interpret.md)
  [`print(`*`<interpret.brma>`*`)`](https://fbartos.github.io/RoBMA/reference/interpret.md)
  : Interpret brma Results
- [`summary_models()`](https://fbartos.github.io/RoBMA/reference/summary_models.md)
  [`print(`*`<summary_models.RoBMA>`*`)`](https://fbartos.github.io/RoBMA/reference/summary_models.md)
  : Summarize Model-Averaged Component Weights
- [`summary_heterogeneity()`](https://fbartos.github.io/RoBMA/reference/summary_heterogeneity.md)
  : Summary of Heterogeneity
- [`summary_heterogeneity(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/summary_heterogeneity.brma.md)
  : Summary of Heterogeneity for brma Objects
- [`pooled_effect()`](https://fbartos.github.io/RoBMA/reference/pooled_effect.md)
  : Pooled Effect Size
- [`pooled_effect(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/pooled_effect.brma.md)
  : Pooled Effect Size for brma Objects
- [`pooled_heterogeneity()`](https://fbartos.github.io/RoBMA/reference/pooled_heterogeneity.md)
  : Pooled Heterogeneity
- [`pooled_heterogeneity(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/pooled_heterogeneity.brma.md)
  : Pooled Heterogeneity for brma Objects
- [`marginal_means()`](https://fbartos.github.io/RoBMA/reference/marginal_means.md)
  : Estimated Marginal Means
- [`marginal_means(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/marginal_means.brma.md)
  : Estimated Marginal Means for brma Objects
- [`summary(`*`<marginal_means.brma>`*`)`](https://fbartos.github.io/RoBMA/reference/summary.marginal_means.brma.md)
  : Summarize Estimated Marginal Means
- [`true_effects()`](https://fbartos.github.io/RoBMA/reference/true_effects.md)
  : True Effects
- [`true_effects(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/true_effects.brma.md)
  : True Effects for brma Objects
- [`ranef()`](https://fbartos.github.io/RoBMA/reference/ranef.md) :
  Random Effects
- [`ranef(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/ranef.brma.md)
  : Random Effects for brma Objects
- [`blup()`](https://fbartos.github.io/RoBMA/reference/blup.md) : Best
  Linear Unbiased Predictions (BLUPs)
- [`blup(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/blup.brma.md)
  : Best Linear Unbiased Predictions for brma Objects
- [`coef(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/coef.brma.md)
  : Extract Model Coefficients for brma Objects
- [`fitted(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/fitted.brma.md)
  : Fitted Values for brma Objects
- [`nobs(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/nobs.brma.md)
  : Number of Observations for brma Objects
- [`predict(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/predict.brma.md)
  : Predict From brma Object

## Model comparison

- [`add_marglik(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md)
  : Add Marginal Likelihood to brma Objects
- [`bridge_sampler(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/bridge_sampler.brma.md)
  : Bridge Sampling for brma Objects
- [`logml(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/logml.brma.md)
  : Log Marginal Likelihood for brma Objects
- [`post_prob(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/post_prob.brma.md)
  : Posterior Model Probabilities for brma Objects
- [`bf(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/bf.brma.md)
  [`bayes_factor(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/bf.brma.md)
  : Bayes Factor for brma Objects
- [`add_loo(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md)
  : Add LOO-PSIS to brma Objects
- [`loo(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/loo.brma.md)
  : LOO-PSIS for brma Objects
- [`loo_compare(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/loo_compare.brma.md)
  : Compare brma Models Using LOO
- [`loo_compare(`*`<loo>`*`)`](https://fbartos.github.io/RoBMA/reference/loo_compare.loo.md)
  : Compare loo Objects Using LOO
- [`loo_weights(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/loo_weights.brma.md)
  : Extract Normalized PSIS Weights from brma Object
- [`check_loo(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/check_loo.brma.md)
  : Check LOO Diagnostics for brma Objects
- [`add_waic(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/add_waic.brma.md)
  : Add WAIC to brma Objects
- [`waic(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/waic.brma.md)
  : WAIC for brma Objects
- [`logLik(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/logLik.brma.md)
  : Extract Log-Likelihood Matrix from brma Object
- [`reexports`](https://fbartos.github.io/RoBMA/reference/reexports.md)
  [`bridge_sampler`](https://fbartos.github.io/RoBMA/reference/reexports.md)
  [`logml`](https://fbartos.github.io/RoBMA/reference/reexports.md)
  [`post_prob`](https://fbartos.github.io/RoBMA/reference/reexports.md)
  [`bf`](https://fbartos.github.io/RoBMA/reference/reexports.md)
  [`bayes_factor`](https://fbartos.github.io/RoBMA/reference/reexports.md)
  : Objects exported from other packages

## Diagnostics and influence

- [`influence(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/influence.brma.md)
  : Measure Influence for brma Objects
- [`cooks.distance(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/cooks.distance.brma.md)
  : Cook's Distance for brma Objects
- [`dfbetas(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/dfbetas.brma.md)
  : DFBETAS for brma Objects
- [`dffits(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/dffits.brma.md)
  : DFFITS for brma Objects
- [`covratio(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/covratio.brma.md)
  : COVRATIO for brma Objects
- [`hatvalues(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/hatvalues.brma.md)
  : Hat Values for brma Objects
- [`residuals(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/residuals.brma.md)
  : Residuals for brma Objects
- [`rstandard(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/rstandard.brma.md)
  : Internally Standardized Residuals for brma Objects
- [`rstudent(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/rstudent.brma.md)
  : Externally Standardized (Studentized) Residuals for brma Objects
- [`vif()`](https://fbartos.github.io/RoBMA/reference/vif.md) : Variance
  Inflation Factors
- [`vif(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/vif.brma.md)
  : Variance Inflation Factors for brma Objects

## Plots

- [`plot(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/plot.brma.md)
  : Plots brma Object
- [`funnel(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/funnel.md)
  : Funnel Plot for brma Object
- [`radial(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/radial.md)
  [`galbraith(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/radial.md)
  : Radial (Galbraith) Plot for brma Object
- [`regplot(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/regplot.md)
  : Regression Plot (Bubble Plot) for brma Object
- [`qqnorm(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/qqnorm.brma.md)
  : Normal QQ Plot for brma Object
- [`plot_prior()`](https://fbartos.github.io/RoBMA/reference/plot_prior.md)
  : Plot Prior Distributions
- [`print_prior()`](https://fbartos.github.io/RoBMA/reference/print_prior.md)
  : Print Prior Distributions
- [`plot_weightfunction()`](https://fbartos.github.io/RoBMA/reference/plot_weightfunction.md)
  : Plots Weight Function of brma Object
- [`plot_pet_peese()`](https://fbartos.github.io/RoBMA/reference/plot_pet_peese.md)
  : Plot PET-PEESE Fit of brma Object
- [`plot_diagnostic()`](https://fbartos.github.io/RoBMA/reference/plot_diagnostic.md)
  [`plot_diagnostic_autocorrelation()`](https://fbartos.github.io/RoBMA/reference/plot_diagnostic.md)
  [`plot_diagnostic_trace()`](https://fbartos.github.io/RoBMA/reference/plot_diagnostic.md)
  [`plot_diagnostic_density()`](https://fbartos.github.io/RoBMA/reference/plot_diagnostic.md)
  : Plot MCMC Diagnostics
- [`plot(`*`<marginal_means.brma>`*`)`](https://fbartos.github.io/RoBMA/reference/plot.marginal_means.brma.md)
  : Plot Estimated Marginal Means

## Posterior samples and zplot diagnostics

- [`as_draws()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md)
  [`as_draws_array()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md)
  [`as_draws_df()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md)
  [`as_draws_list()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md)
  [`as_draws_matrix()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md)
  [`as_draws_rvars()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md)
  : Convert brma Objects to posterior Draws Formats
- [`as_draws(`*`<brma_samples>`*`)`](https://fbartos.github.io/RoBMA/reference/as_draws.brma_samples.md)
  [`as_draws_array(`*`<brma_samples>`*`)`](https://fbartos.github.io/RoBMA/reference/as_draws.brma_samples.md)
  [`as_draws_df(`*`<brma_samples>`*`)`](https://fbartos.github.io/RoBMA/reference/as_draws.brma_samples.md)
  [`as_draws_list(`*`<brma_samples>`*`)`](https://fbartos.github.io/RoBMA/reference/as_draws.brma_samples.md)
  [`as_draws_matrix(`*`<brma_samples>`*`)`](https://fbartos.github.io/RoBMA/reference/as_draws.brma_samples.md)
  [`as_draws_rvars(`*`<brma_samples>`*`)`](https://fbartos.github.io/RoBMA/reference/as_draws.brma_samples.md)
  : Convert brma_samples to posterior Draws Formats
- [`as.matrix(`*`<brma_samples>`*`)`](https://fbartos.github.io/RoBMA/reference/as.matrix.brma_samples.md)
  : Convert brma_samples to Matrix
- [`print(`*`<brma_samples>`*`)`](https://fbartos.github.io/RoBMA/reference/print.brma_samples.md)
  : Print brma_samples Object
- [`summary(`*`<brma_samples>`*`)`](https://fbartos.github.io/RoBMA/reference/summary.brma_samples.md)
  : Summarize brma_samples Object
- [`print(`*`<summary.brma_samples>`*`)`](https://fbartos.github.io/RoBMA/reference/print.summary.brma_samples.md)
  : Print summary.brma_samples Object
- [`as_zplot(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/as_zplot.brma.md)
  : Transform brma Object to Zplot
- [`zplot(`*`<brma>`*`)`](https://fbartos.github.io/RoBMA/reference/zplot.brma.md)
  [`zplot(`*`<zplot_brma>`*`)`](https://fbartos.github.io/RoBMA/reference/zplot.brma.md)
  : Plot Zplot Diagnostics Directly
- [`summary(`*`<zplot_brma>`*`)`](https://fbartos.github.io/RoBMA/reference/summary.zplot_brma.md)
  : Summarize Zplot Results
- [`print(`*`<summary.zplot_brma>`*`)`](https://fbartos.github.io/RoBMA/reference/print.summary.zplot_brma.md)
  : Print Zplot Summary
- [`plot(`*`<zplot_brma>`*`)`](https://fbartos.github.io/RoBMA/reference/plot.zplot_brma.md)
  : Plot Zplot Results
- [`hist(`*`<zplot_brma>`*`)`](https://fbartos.github.io/RoBMA/reference/hist.zplot_brma.md)
  : Histogram of Z-Statistics
- [`lines(`*`<zplot_brma>`*`)`](https://fbartos.github.io/RoBMA/reference/lines.zplot_brma.md)
  : Add Zplot Density Lines

## Example datasets

- [`Anderson2010`](https://fbartos.github.io/RoBMA/reference/Anderson2010.md)
  : 23 experimental studies from Anderson et al. (2010) that meet the
  best practice criteria
- [`Andrews2021`](https://fbartos.github.io/RoBMA/reference/Andrews2021.md)
  : 39 study rows on household chaos and child executive functions from
  a meta-analysis by Andrews et al. (2021)
- [`Bem2011`](https://fbartos.github.io/RoBMA/reference/Bem2011.md) : 9
  experimental studies from Bem (2011) as described in Bem et al. (2011)
- [`Havrankova2025`](https://fbartos.github.io/RoBMA/reference/Havrankova2025.md)
  : 1159 effect sizes from a meta-analysis of beauty and professional
  success by Havránková et al. (2025)
- [`Hoppen2025`](https://fbartos.github.io/RoBMA/reference/Hoppen2025.md)
  : 37 studies from a meta-analysis of social comparison as a behavior
  change technique by Hoppen et al. (2025)
- [`Johnides2025`](https://fbartos.github.io/RoBMA/reference/Johnides2025.md)
  : 412 effect sizes from a meta-analysis of secondary benefits of
  family-based treatments by Johnides et al. (2025)
- [`Kroupova2021`](https://fbartos.github.io/RoBMA/reference/Kroupova2021.md)
  : 881 estimates from 69 studies of a relationship between employment
  and educational outcomes collected by Kroupova et al. (2021)
- [`Lui2015`](https://fbartos.github.io/RoBMA/reference/Lui2015.md) : 18
  studies of a relationship between acculturation mismatch and
  intergenerational cultural conflict collected by Lui (2015)
- [`ManyLabs16`](https://fbartos.github.io/RoBMA/reference/ManyLabs16.md)
  : 55 effect sizes from Many Labs 2 replication studies of Tversky and
  Kahneman (1981) framing effects
- [`Poulsen2006`](https://fbartos.github.io/RoBMA/reference/Poulsen2006.md)
  : 5 studies with a tactile outcome assessment from Poulsen et
  al. (2006) of the effect of potassium-containing toothpaste on dentine
  hypersensitivity
- [`Wang2025`](https://fbartos.github.io/RoBMA/reference/Wang2025.md) :
  70 effect sizes from a meta-analysis of ChatGPT's impact on student
  learning by Wang and Fan (2025)
- [`Weingarten2018`](https://fbartos.github.io/RoBMA/reference/Weingarten2018.md)
  : 582 effect sizes examining the ease-of-retrieval effect from a
  meta-analysis by Weingarten and Hutchinson (2018)

## Print helpers

- [`print(`*`<RoBMA_data>`*`)`](https://fbartos.github.io/RoBMA/reference/print.RoBMA_data.md)
  : Print method for RoBMA_data objects
- [`print(`*`<marginal_means.brma>`*`)`](https://fbartos.github.io/RoBMA/reference/print.marginal_means.brma.md)
  : Print Estimated Marginal Means
- [`print(`*`<summary.marginal_means.brma>`*`)`](https://fbartos.github.io/RoBMA/reference/print.summary.marginal_means.brma.md)
  : Print Summary of Estimated Marginal Means
- [`print(`*`<summary_heterogeneity.brma>`*`)`](https://fbartos.github.io/RoBMA/reference/print.summary_heterogeneity.brma.md)
  : Print Summary of Heterogeneity
- [`print(`*`<vif.brma>`*`)`](https://fbartos.github.io/RoBMA/reference/print.vif.brma.md)
  : Print VIF Results
