## version 4.1.0 (IN PROGRESS)
### Breaking changes
- requires R >= 4.3.0, BayesTools >= 0.4.0, loo >= 2.10.0, and JAGS 4.x (>= 4.3.1), and imports Matrix, ps, and statmod
- models fitted with RoBMA 4.0.0 must be refitted; functions that read the fitted parameters stop on them with a `BayesTools_refit_required` error
- `AIC()` and `BIC()` stop with an explanation instead of returning `numeric(0)`
- removes the ineffective `seed` argument of `update()`
- `radial()` stops for likelihood-weighted fits, and cluster-unit LOO and WAIC of binomial and Poisson GLMMs stop (use `unit = "estimate"`)
- fits with mean-difference or orthonormal contrasts stop when a factor interacts with a predictor whose main effect is missing (e.g., `mods = ~ g + g:x`); use `mods = ~ g * x`, an independent prior for the term, or treatment contrasts
- explicitly supplied `prior_mods`, `prior_scale`, and `prior_heterogeneity_allocation` stop when the model has no corresponding component instead of being ignored
- nonfinite `weights` are rejected
- the default Poisson GLMM `prior_lograte` is `Normal(log(pooled crude rate), 1)`, independent of the exposure-time unit (set its SD with `RoBMA.options(default_lograte.sd = )`)
- `BMA.glmm()` uses mean-difference contrasts by default, as the other model-averaging constructors do
- `plot()` of `zplot()` results no longer draws the extrapolated curve by default (`plot_extrapolation = FALSE`)
- selection-model p-value cut points are compared exactly instead of with a 1e-12 tolerance
- seeded fits give different draws than 4.0.0 (chain seeding of BayesTools 0.4.0), and seeded fitting no longer resets R's global random-number stream
- `predict()` gains `conditioning_depth` (`"marginal"` by default, `"cluster"`, or `"estimate"`), which changes the default same-data `type = "estimate"` predictions and the clustered and GLMM response predictions
- `zplot()` targets the marginal posterior predictive distribution by default (`conditioning_depth`)
- `funnel()` draws plug-in contours; the posterior-averaged contours of 4.0.0 are drawn by the new `bfunnel()`
- estimate-unit LOO and WAIC of Gaussian models integrate the local random effects, conditioning each estimate on the other estimates of its cluster
- marginal-means Bayes factors follow the Savage-Dickey rules of BayesTools 0.4.0: a null hypothesis outside the posterior draws gives a finite Bayes factor, and means fixed at the null or without an exact prior ordinate give `NA` with the reason
- `as_draws()` returns the parameters shown by `print()` and `summary()` by default; auxiliary parameters require `include_auxiliary = TRUE`
- labels are rendered from the BayesTools parameter catalog: factor levels are selected by their labels, and contrast coefficients that are not levels are labelled `g{j}`
- the scale intercept of `BMA()` and `RoBMA()` summaries is labelled `exp(intercept)`, and a location intercept fixed at zero is omitted from meta-regression output

### Features
- adds `brma.mv()` for meta-analysis with a known sampling covariance matrix `V` and formula random effects (`random`, e.g., `~ 1 | study / esid`) with diagonal, unstructured, compound-symmetry, AR(1), heterogeneous, and continuous-time AR(1) structures, and known random-effect group covariance (`R`)
- adds `BMA.mv()`, `RoBMA.mv()`, `bPET.mv()`, `bPEESE.mv()`, and `bselmodel.mv()` for model-averaged and publication-bias-adjusted multivariate and multilevel models
- re-exports BayesTools' `prior_random()` interface for structured random-effect priors (variance allocations, inclusion gates, and LKJ priors), completed with RoBMA's default priors
- random-effect quantities (`tau`, `tau2`, `rho`, `tau_total`, `tau2_prop`, ...) are named consistently across summaries, plots, hypotheses, and `as_draws()`, and random-effect inclusion is reported in the component inclusion table
- adds `component` to `summary_heterogeneity()`, `pooled_heterogeneity()`, `ranef()` (also `expand` and `simplify`), `plot()`, `plot_prior()`, `print_prior()`, and `plot_diagnostic()`
- re-exports BayesTools' `prior_ordered()` for ordered-factor moderators
- `mods`, `scale`, and `yi ~` formulas accept `^` crossing (e.g., `(x + z)^2`)
- adds conditional Gaussian vector selection models specified with `selection = selection_model(...)`, with product or best p-value weight rules, publication groups, and deterministic quadrature with a controlled QMC fallback (`set_selection_likelihood_control()`)
- adds caching of selection normalizers (`selection_cache_info()`, `remove_selection_cache()`) and a coarse-corrected selection sampler (`RoBMA.options(selection.sampler = "coarse_corrected")`)
- adds `selection_sensitivity_diagnostics()` and `add_selection_sensitivity_diagnostics()`
- adds `hypothesis()` (aliases `bf_hypothesis()` and `BF_hypothesis()`) for point, interval, and directional hypotheses on parameters, factor levels, and their linear combinations in fitted models and marginal means, and `hypothesis_quantities()`, which lists the testable quantities
- adds the likelihood-aware qCMDE and IWMDE posterior density and ordinate estimators (`density_method`, `density_control`) for plots, hypotheses, and marginal means, with `density_diagnostics()`
- adds `check_indicators`, `monitor`, and `allow_not_assessable` to convergence checks and autofit
- adds `log_lik()`, a `loo_model_weights()` method, and `parallel`, `cores`, `repetitions`, `method`, `maxiter`, and `silent` to `add_marglik()`
- supports fully fixed point-prior models
- adds prediction intervals to `pooled_effect()`
- adds `max_samples` to `hatvalues()`, `residuals()`, `rstandard()`, `qqnorm()`, and `vif()`, and `component` and `parameter` to `dfbetas()` and `covratio()`
- `vif()` computes design-based VIF/GVIF in each posterior draw, with optional credible intervals (`probs`) and draws
- adds `bfunnel()`, a `metafor::forest()` method with `as_metafor_forest()`, and `lines()` methods for fitted models and marginal means
- `transform = "EXP"` and `"LOG"` apply to individual coefficients and heterogeneity parameters in plots
- `regplot()` accepts moderators that appear only in the scale formula, and `radial()` supports `brma.mv()` models
- credible bands of `regplot()`, `zplot()`, and `plot_pet_peese()` use Harrell-Davis quantiles
- `as.data.frame()` converts printable results to long data frames
- `update()` accepts named `NULL` values to clear convergence and autofit settings
- `set_contrast_factor_predictors = "independent"` defines one coefficient per factor level
- adds the `Hoogeveen2023` data set
- adds the options `RoBMA.options(native_threads = )` and `RoBMA.options(jags.worker_output = )`

### Changes
- `cooks.distance()` uses the unscaled squared Mahalanobis distance, matching metafor (values are multiplied by the fixed-effect rank relative to 4.0.0)
- influence diagnostics that reuse the cached LOO result warn about high Pareto k values
- methods warn about unknown arguments in `...` instead of ignoring them
- row selection (`subset`, missing values) drops unused factor levels, except for the index factor of structured random-effect terms
- unspecified factor null priors use the contrast of the alternative prior; conflicting contrasts stop with an error

### Fixes
- marginal-means Bayes factors use the exact kernel density at the null; the grid estimate understated Bayes factors of long-tailed posteriors
- estimate-unit GLMM log-likelihoods use adaptive Gauss-Hermite quadrature or a convergence-checked grid instead of a fixed 30-node grid, which was inaccurate for some nuisance priors and failed for sparse data
- GLMM predictions accept the four-cell binomial count representation and return raw count samples
- `predict()` stops when `newdata` lacks a formula variable instead of filling outcome-named columns with zeros
- a response formula without an intercept (`yi ~ 0`) fixes the location intercept at zero
- prior-only objects evaluate formulas with standardized predictors, and their `summary()` no longer fails
- weighted residual diagnostics use the original outcome variances
- funnel contours and location-scale influence heterogeneity use the root mean square heterogeneity
- cluster-level selection likelihood integrals are certified by successive quadrature rules instead of heuristic fallbacks
- numerically stable cluster, known-covariance, and selected-normal likelihoods for precise or aligned estimates
- cached LOO/WAIC results and marginal-likelihood comparisons check that the data and targets match
- `plot_weightfunction()` marks two-sided p-values for two-sided weight functions
- transformed posterior plots take `xlim` on the displayed scale
- `ranef(component = "total")` keeps the estimate names
- `vif()`, `hatvalues()`, and related diagnostics of `RoBMA()` ensembles use the bias design of each member
- `dfbetas()`, `influence()`, and `vif()` omit a fixed-zero intercept, and `covratio()` is finite for mean-difference and orthonormal contrasts
- `as_zplot()` and `zplot()` validate `significance_level`
- fitted objects no longer store the calling environment of their formulas
- the RoBMA JAGS module is unregistered at process exit, preventing intermittent crashes on Windows
- clearer input and diagnostic messages

### Performance
- saved fits no longer contain the compiled rjags model
- post-fit native kernels run in parallel (`RoBMA.options(native_threads = )`) and use AVX2/AVX-512 for selection normalizers
- faster multilevel cluster likelihoods, selection normalizers, funnel, regplot, and zplot computations, and LOO

### Documentation
- adds the vignettes "Posterior Densities and Parameter Tests" and "Multivariate and Multilevel Meta-Analysis"
- documents the qCMDE/IWMDE estimators, their diagnostics, and the Savage-Dickey requirements

## version 4.0.0
### Breaking changes
- rewrites the package around the unified `brma` class hierarchy. Single-model fits now use `brma()`, `brma.glmm()`, `bselmodel()`, `bPET()`, and `bPEESE()`; model-averaged fits use `BMA()`, `BMA.glmm()`, and `RoBMA()`.
- removes the legacy `RoBMA.reg()`, `NoBMA()`, `NoBMA.reg()`, `BiBMA()`, and `BiBMA.reg()` constructors. Use `mods`, `scale`, and `cluster` in the new constructors, `BMA()` for no-bias normal-likelihood model averaging, and `BMA.glmm()` for GLMM model averaging.
- replaces old input aliases such as `d`, `r`, `logOR`, `OR`, `z`, `y`, `se`, `v`, `n`, `study_names`, `study_ids`, `weight`, and `transformation` with `yi`, `vi`/`sei`, `ni`, `slab`, `cluster`, `weights`, `measure`, `output_measure`, and `transform`.
- removes legacy helper APIs including `combine_data()`, `check_setup()`, `extract_posterior()`, `marginal_summary()`, `marginal_plot()`, `plot_models()`, `adjusted_effect()`, `as_zcurve()`, and the old z-curve plotting methods.
- normal-likelihood fitting functions now require an explicit `measure` for fitted models. Use `measure = "GEN"` for generic effect sizes without a known unit-information scale.
- `update()` for `brma` objects now focuses on extending MCMC samples, updating labels, and refreshing cached quantities, not changing model structure.
- `set_convergence_checks()` no longer accepts the old `remove_failed` and `balance_probability` arguments.

### Features
- adds `brma()` / `brma.norm()` for single normal-likelihood Bayesian meta-analysis, including random-effects, meta-regression, multilevel, and location-scale models.
- adds `brma.glmm()` for binomial-normal and Poisson-normal GLMM meta-analysis from raw two-arm counts (`measure = "OR"` and `"IRR"`).
- adds single-model publication-bias constructors `bselmodel()`, `bPET()`, and `bPEESE()`.
- adds `BMA()` / `BMA.norm()` for Bayesian model averaging without publication-bias adjustment.
- adds `BMA.glmm()` for Bayesian model averaging of GLMM meta-analyses without publication-bias adjustment.
- rewrites `RoBMA()` as a product-space model-averaged ensemble over effect, heterogeneity, moderator, scale, and publication-bias components.
- adds formula/data-frame input handling for effect sizes, moderators, scale predictors, clusters, labels, subsets, likelihood weights, and raw GLMM counts.
- adds default prior construction from standardized effect-size measures, estimated or manually supplied unit-information standard deviations, and informed empirical priors.
- adds `prior_weightfunction()`, `wf_cumulative()`, `wf_fixed()`, and `wf_independent()` for BayesTools-backed selection-weightfunction priors.
- adds `prior_PET()`, `prior_PEESE()`, `prior_none()`, `prior_factor()`, `prior_informed()`, and BayesTools contrast helpers as package-level prior utilities.
- adds `posterior` package interfaces via `as_draws()`, `as_draws_array()`, `as_draws_df()`, `as_draws_list()`, `as_draws_matrix()`, and `as_draws_rvars()` for fitted models and `brma_samples`.
- adds the `brma_samples` posterior-sample class with print, summary, matrix, and `posterior` conversion methods.
- adds `predict.brma()` for posterior predictions of fixed terms, cluster effects, latent true effects, observed responses, and scale terms, with `newdata`, `conditional`, `bias_adjusted`, `output_measure`, and `transform` support.
- adds convenience wrappers `fitted()`, `pooled_effect()`, `pooled_heterogeneity()`, `blup()`, `true_effects()`, and `ranef()` for `brma` objects.
- adds model-comparison helpers `add_loo()`, `loo()`, `loo_compare()`, `loo_weights()`, `check_loo()`, `add_waic()`, and `waic()` using the `loo` package.
- adds bridge-sampling marginal likelihood support for single-model `brma` fits via `add_marglik()`, `bridge_sampler()`, `logml()`, `bf()`, `bayes_factor()`, and `post_prob()`.
- adds residual and influence diagnostics: `residuals()`, `rstandard()`, `rstudent()` / `LOO-PIT`, `hatvalues()`, `influence()`, `dfbetas()`, `dffits()`, `cooks.distance()`, `covratio()`, and `vif()`.
- adds plotting methods for `brma` objects: posterior/prior plots, `funnel()`, `regplot()`, `qqnorm()`, `radial()` / `galbraith()`, MCMC diagnostic plots, weightfunction plots, and PET-PEESE plots.
- adds `marginal_means()` with summary and plotting methods for moderator models.
- adds `summary_models()` for marginal and individual model-weight summaries of product-space `RoBMA`, `BMA`, and `BMA.glmm` objects.
- adds `interpret()` for concise textual interpretation of fitted `brma` and model-averaged objects.
- renames the zplot diagnostic API to `as_zplot()` and adds the direct plotting wrapper `zplot()`, with `plot()`, `hist()`, `lines()`, `summary()`, and print methods for zplot objects.
- adds `RoBMA.options()` and `RoBMA.get_option()` package options for defaults such as core count, automatic LOO/WAIC/marginal-likelihood computation, prior scaling defaults, and selection-bias defaults.

### Changes
- renames the multilevel clustering argument to `cluster`.
- renames study labels to `slab`, matching `metafor` naming.
- renames likelihood weights to `weights` for supported constructors and applies them consistently to posterior fitting, log-likelihoods, LOO, WAIC, and diagnostics; `brma.mv()` currently rejects likelihood weights.
- uses `measure`, `output_measure`, and `transform` for effect-size scale handling. Supported conversions include `SMD`, `COR`, `ZCOR`, and `OR`; `transform = "EXP"` exponentiates log ratio measures for display.
- standardizes continuous predictors by default and transforms reported coefficients back to the original scale unless standardized coefficients are requested.
- uses treatment contrasts by default for single-model constructors and mean-difference contrasts by default for model-averaged constructors.
- changes `predict.brma()` default to `type = "terms"`. GLMM `type = "response"` predictions return continuity-corrected effect-size estimators by default via `as_measure = TRUE`.
- separates output `unit` from `conditioning_depth` for residuals, fitted values, LOO, WAIC, and related diagnostics.
- supports estimate-level LOO/WAIC targets for `brma.mv()` and estimate-/cluster-level LOO/WAIC targets for multilevel `brma()` models, with target metadata to prevent invalid comparisons.
- keeps bridge-sampling marginal likelihoods for single-model `brma` objects; product-space `RoBMA`, `BMA`, and `BMA.glmm` objects rely on product-space only.
- routes selection-weightfunction priors through the BayesTools selection backend and selected-normal kernel, removing legacy weighted-normal mapping paths.
- uses `bias_indicator` and branch-aware selected-normal contexts for RoBMA publication-bias mixtures instead of inferring selection branches from `omega`.
- increases zplot default posterior thinning controls to `10000` samples and accepts `Inf` where full posterior evaluation is requested.
- adds `max_samples` controls to expensive funnel, regplot, and zplot summaries.
- updates the package startup message to point users to `vignette("v00-introduction", package = "RoBMA")`.
- requires BayesTools 0.3.0 for forward API and selection-backend support.
- adds `bridgesampling`, `loo`, `MASS`, and `parallel` as imports and `posterior` as a suggested package.

### Fixes
- fixes loading and runtime checks for the RoBMA JAGS module and native R routines.

### Performance and internals
- moves fitting to JAGS product-space models with mixture-prior indicators for model averaging.
- replaces legacy weighted-normal and multivariate-normal native code with selected-normal kernels shared by JAGS and R-native calls.
- adds native selected-normal routines for log likelihoods, normalizers, CDFs, moments, RNG, weighted summaries, funnel contours, regplot intervals, and zplot densities/threshold summaries.
- adds native estimate-unit GLMM marginal log-likelihood helpers for binomial and Poisson models.
- caches selected-normal normalizers and uses telescoping selection probabilities with log-space fallbacks for better numerical stability.
- relocates selected-normal C++ code to `src/selnorm/` and updates `Makevars*`, native registration, cleanup rules, and JAGS distribution registration.
- removes unused native matrix/LAPACK helper sources and older source-level transformation helpers.

### Documentation and tests
- reorganizes vignettes into numbered workflows covering introduction, prior distributions, baseline Bayesian meta-analysis, feature coverage, metafor parity, model averaging, RoBMA, multilevel models, medicine examples, and zplot diagnostics.
- regenerates roxygen documentation for the new constructors, priors, predictions, summaries, diagnostics, plots, model-comparison methods, and datasets.
- refreshes the README and pkgdown site for the 4.0.0 API.
- adds cached model fits under numbered vignette/model directories.
- refactors tests into ordered input, fitting, prediction, plotting, diagnostics, model-comparison, selected-normal kernel, and vignette-cache coverage.
- adds regression tests for selected-normal telescope probabilities, native/R fallback parity, posterior-row alignment, GLMM response conversion, LOO/WAIC targets, bridge sampling, and visual outputs.

## version 3.6.1
### Features
- `Explanation` vignette that helps navigate users through the vignettes
- two vignettes demonstrating robust Bayesian meta-analysis and meta-regressions
- `summary()` function now provides publication bias model type summary (`type = "models"`) for models fitted using `algorithm = "ss"`
- improves control over zplot diagnostics (i.e., specifying col, border, etc for the individual elements)

## version 3.6
### Features
- `funnel()` plot to visualize residuals vs the expected sampling distribution for `RoBMA()` and `RoBMA.reg()` models when using the `algorithm = "ss"`
- `residuals()` method for `RoBMA()` and `RoBMA.reg()` models when using the `algorithm = "ss"`
- `as_zplot()` function to transform meta-analytic models into a zplot object, only available for `RoBMA()` and `RoBMA.reg()` fitted using the `algorithm = "ss"`
- `plot()`, `summary()`, and `print()` functions for the `as_zplot` objects

## version 3.5.1
### Features
- `summary()` function now supports a `standardized_coefficients` argument to report either standardized (default) or raw meta-regression coefficients
- `extract()` function to extract the posterior samples of the model parameters
- `true_effects()` function to summarize the true effect size estimates of `RoBMA()` and `RoBMA.reg()` models when using the `algorithm = "ss"`
- `predict()` method for `RoBMA()` and `RoBMA.reg()` models when using the `algorithm = "ss"`

### Fixes
- fitting a meta-regression using predictors with missing values result in a clear error message

### Changes
- improving the speed of unit tests

## version 3.5
### Features
- approximate and computationally feasibly 3lvl selection models via the `RoBMA()` and `RoBMA.reg()` functions with the `cluster` argument when using `algorithm = "ss"`
- 3lvl binomial-normal models for binary data via the `BiBMA` and `BiBMA.reg` functions with the `cluster` argument when using `algorithm = "ss"`
- `pooled_effect()` function to compute the pooled effect size from the `RoBMA.reg`, `NoBMA.reg`, and `BiBMA.reg` models
- `adjusted_effect()` function to compute the adjusted effect size from the `RoBMA.reg`, `NoBMA.reg`, and `BiBMA.reg` models
- enables `summary_heterogeneity()` for BiBMA models

### Fixes
- passing and checks of the `cluster` and `study_labels` arguments
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
- weighted meta-analysis by specifying `cluster` argument in `RoBMA()` and setting `weighted = TRUE`. The likelihood contribution of estimates from each study is down-weighted proportionally to the number of estimates in that study. Note that this experimental feature is supposed to provide a conservative alternative for estimating RoBMA in cases with multiple estimates from a study where the multivariate option is not computationally feasible.

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
- three-level meta-analysis by specifying `cluster` argument in `RoBMA`. However, note that this is (1) an experimental feature and (2) the computational expense of fitting selection models with clustering is extreme. As of now, it is almost impossible to have more than 2-3 estimates clustered within a single study).

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
