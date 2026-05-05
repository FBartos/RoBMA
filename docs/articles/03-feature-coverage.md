# Feature Coverage

This vignette overviews current features of the `RoBMA` R package
([Bartoš & Maier, 2020](#ref-RoBMA)). The table is organized by
substantive capability. Some rows point to additional features not (yet)
available in the `RoBMA` R package but featured in the `metafor` package
([Viechtbauer, 2010](#ref-metafor)) as reference points.

Green ticks mark available features. Cells marked `limited` mean the
feature exists but with narrower model coverage.

| Topic | Feature | [brma.norm](https://fbartos.github.io/RoBMA/reference/brma.md) | [brma.glmm](https://fbartos.github.io/RoBMA/reference/brma.glmm.md) | [bselmodel](https://fbartos.github.io/RoBMA/reference/bselmodel.md) | [bPET](https://fbartos.github.io/RoBMA/reference/bPET.md)/[bPEESE](https://fbartos.github.io/RoBMA/reference/bPEESE.md) | [BMA.norm](https://fbartos.github.io/RoBMA/reference/BMA.md) | [BMA.glmm](https://fbartos.github.io/RoBMA/reference/BMA.glmm.md) | [RoBMA](https://fbartos.github.io/RoBMA/reference/RoBMA.md) |
|----|----|----|----|----|----|----|----|----|
| Model / Structure | Fixed- and random-effects models | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Moderation (`mods`) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Location-scale models (`scale`) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Multilevel models (`cluster`) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Model averaging | No | No | No | No | ✓ | ✓ | ✓ |
|  | Inclusion Bayes factors | No | No | No | No | ✓ | ✓ | ✓ |
|  | General multivariate / covariance structures | No | No | No | No | No | No | No |
|  | Robust / sandwich inference | No | No | No | No | No | No | No |
|  | Refit / update ([`update()`](https://fbartos.github.io/RoBMA/reference/update.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Priors | Default prior distributions | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Empirical/informed prior distributions ([`prior_informed()`](https://fbartos.github.io/RoBMA/reference/prior_informed.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Custom prior distributions ([`prior()`](https://fbartos.github.io/RoBMA/reference/prior.md), [`prior_factor()`](https://fbartos.github.io/RoBMA/reference/prior_factor.md), [`prior_weightfunction()`](https://fbartos.github.io/RoBMA/reference/prior_weightfunction.md), [`prior_none()`](https://fbartos.github.io/RoBMA/reference/prior_none.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | PET / PEESE prior distributions ([`prior_PET()`](https://fbartos.github.io/RoBMA/reference/prior_PET.md), [`prior_PEESE()`](https://fbartos.github.io/RoBMA/reference/prior_PEESE.md)) | No | No | No | ✓ | No | No | ✓ |
|  | Prior-only inspection (`only_priors = TRUE`) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Weight-function shape ([`wf_cumulative()`](https://fbartos.github.io/RoBMA/reference/prior_weightfunction.md), [`wf_fixed()`](https://fbartos.github.io/RoBMA/reference/prior_weightfunction.md), [`wf_independent()`](https://fbartos.github.io/RoBMA/reference/prior_weightfunction.md)) | No | No | ✓ | No | No | No | ✓ |
|  | Factor contrasts ([`contr.treatment()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/contrast.html), [`contr.meandif()`](https://fbartos.github.io/RoBMA/reference/contr.BayesTools.md), [`contr.orthonormal()`](https://fbartos.github.io/RoBMA/reference/contr.BayesTools.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | UISD estimation ([`estimate_unit_information_sd()`](https://fbartos.github.io/RoBMA/reference/estimate_unit_information_sd.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Publication Bias | Selection models | No | No | ✓ | No | No | No | ✓ |
|  | PET / PEESE models | No | No | No | ✓ | No | No | ✓ |
|  | Model averaging over bias models | No | No | No | No | No | No | ✓ |
|  | Bias-adjusted summaries / predictions | No | No | ✓ | ✓ | No | No | ✓ |
|  | Trim-fill / fail-safe N | No | No | No | No | No | No | No |
| Prediction | Pooled effect / heterogeneity summaries ([`pooled_effect()`](https://fbartos.github.io/RoBMA/reference/pooled_effect.brma.md), [`pooled_heterogeneity()`](https://fbartos.github.io/RoBMA/reference/pooled_heterogeneity.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Heterogeneity decomposition ([`summary_heterogeneity()`](https://fbartos.github.io/RoBMA/reference/summary_heterogeneity.brma.md), tau², I², H²) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Fitted-value extraction ([`fitted()`](https://fbartos.github.io/RoBMA/reference/fitted.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Prediction for new covariate values ([`predict()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Posterior predictive response summaries ([`predict()`](https://fbartos.github.io/RoBMA/reference/predict.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | True effects / BLUPs / random effects ([`true_effects()`](https://fbartos.github.io/RoBMA/reference/true_effects.brma.md), [`blup()`](https://fbartos.github.io/RoBMA/reference/blup.brma.md), [`ranef()`](https://fbartos.github.io/RoBMA/reference/ranef.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Estimated marginal means ([`marginal_means()`](https://fbartos.github.io/RoBMA/reference/marginal_means.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Plots | Posterior plots ([`plot()`](https://fbartos.github.io/RoBMA/reference/plot.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Prior plots ([`plot_prior()`](https://fbartos.github.io/RoBMA/reference/plot_prior.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Marginal means plots ([`marginal_means()`](https://fbartos.github.io/RoBMA/reference/marginal_means.brma.md) + [`plot()`](https://fbartos.github.io/RoBMA/reference/plot.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Regression plots ([`regplot()`](https://fbartos.github.io/RoBMA/reference/regplot.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Funnel plots ([`funnel()`](https://fbartos.github.io/RoBMA/reference/funnel.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | zplot diagnostics ([`as_zplot()`](https://fbartos.github.io/RoBMA/reference/as_zplot.brma.md)/[`zplot()`](https://fbartos.github.io/RoBMA/reference/zplot.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Weight-function plots ([`plot_weightfunction()`](https://fbartos.github.io/RoBMA/reference/plot_weightfunction.md)) | No | No | ✓ | No | No | No | ✓ |
|  | PET-PEESE plots ([`plot_pet_peese()`](https://fbartos.github.io/RoBMA/reference/plot_pet_peese.md)) | No | No | No | ✓ | No | No | ✓ |
|  | Radial / Galbraith plots ([`radial()`](https://fbartos.github.io/RoBMA/reference/radial.md), [`galbraith()`](https://fbartos.github.io/RoBMA/reference/radial.md)) | ✓ | limited | ✓ | ✓ | ✓ | limited | ✓ |
|  | Forest plots ([`forest()`](https://wviechtb.github.io/metafor/reference/forest.rma.html)) | No | No | No | No | No | No | No |
|  | Baujat plots ([`baujat()`](https://wviechtb.github.io/metafor/reference/baujat.html)) | No | No | No | No | No | No | No |
|  | GOSH plots ([`gosh()`](https://wviechtb.github.io/metafor/reference/gosh.html)) | No | No | No | No | No | No | No |
|  | L’Abbe plots ([`labbe()`](https://wviechtb.github.io/metafor/reference/labbe.html)) | No | No | No | No | No | No | No |
| Residuals / Diagnostics | Raw residuals ([`residuals()`](https://fbartos.github.io/RoBMA/reference/residuals.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Pearson / standardized residuals ([`rstandard()`](https://fbartos.github.io/RoBMA/reference/rstandard.brma.md)) | ✓ | No | No | ✓ | ✓ | No | limited |
|  | Studentized residuals (LOO-PIT) ([`rstudent()`](https://fbartos.github.io/RoBMA/reference/rstudent.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Q-Q plots ([`qqnorm()`](https://fbartos.github.io/RoBMA/reference/qqnorm.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Cook’s distances ([`cooks.distance()`](https://fbartos.github.io/RoBMA/reference/cooks.distance.brma.md)) | ✓ | No | No | ✓ | ✓ | No | No |
|  | DFBETAS ([`dfbetas()`](https://fbartos.github.io/RoBMA/reference/dfbetas.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | DFFITS ([`dffits()`](https://fbartos.github.io/RoBMA/reference/dffits.brma.md)) | ✓ | No | No | ✓ | ✓ | No | No |
|  | Covariance ratios ([`covratio()`](https://fbartos.github.io/RoBMA/reference/covratio.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Hat values ([`hatvalues()`](https://fbartos.github.io/RoBMA/reference/hatvalues.brma.md)) | ✓ | No | No | ✓ | ✓ | No | No |
|  | Combined influence summary ([`influence()`](https://fbartos.github.io/RoBMA/reference/influence.brma.md)) | ✓ | limited | limited | ✓ | ✓ | limited | limited |
|  | LOO / WAIC diagnostics ([`check_loo()`](https://fbartos.github.io/RoBMA/reference/check_loo.brma.md), [`loo::pareto_k_ids()`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html), [`loo::pareto_k_table()`](https://mc-stan.org/loo/reference/pareto-k-diagnostic.html)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Moderator collinearity diagnostics ([`vif()`](https://fbartos.github.io/RoBMA/reference/vif.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Refit leave-one-out / permutation tests | No | No | No | No | No | No | No |
| MCMC Diagnostics | Posterior summaries (Rhat, ESS, MCMC error) ([`summary()`](https://fbartos.github.io/RoBMA/reference/summary.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Trace plots ([`plot_diagnostic_trace()`](https://fbartos.github.io/RoBMA/reference/plot_diagnostic.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Density plots ([`plot_diagnostic_density()`](https://fbartos.github.io/RoBMA/reference/plot_diagnostic.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Autocorrelation plots ([`plot_diagnostic_autocorrelation()`](https://fbartos.github.io/RoBMA/reference/plot_diagnostic.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Diagnostic-plot wrapper ([`plot_diagnostic()`](https://fbartos.github.io/RoBMA/reference/plot_diagnostic.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Model Comparison | Marginal likelihood ([`add_marglik()`](https://fbartos.github.io/RoBMA/reference/add_marglik.brma.md), [`logml()`](https://fbartos.github.io/RoBMA/reference/logml.brma.md)) | ✓ | ✓ | ✓ | ✓ | No | No | No |
|  | WAIC/LOO ([`add_waic()`](https://fbartos.github.io/RoBMA/reference/add_waic.brma.md), [`add_loo()`](https://fbartos.github.io/RoBMA/reference/add_loo.brma.md), [`waic()`](https://fbartos.github.io/RoBMA/reference/waic.brma.md), [`loo()`](https://fbartos.github.io/RoBMA/reference/loo.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Bayes factors ([`bf()`](https://fbartos.github.io/RoBMA/reference/bf.brma.md), [`bayes_factor()`](https://fbartos.github.io/RoBMA/reference/bf.brma.md), [`post_prob()`](https://fbartos.github.io/RoBMA/reference/post_prob.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | WAIC comparison ([`loo_compare()`](https://fbartos.github.io/RoBMA/reference/loo_compare.brma.md), [`loo_weights()`](https://fbartos.github.io/RoBMA/reference/loo_weights.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | LOO comparison ([`loo_compare()`](https://fbartos.github.io/RoBMA/reference/loo_compare.brma.md), [`loo_weights()`](https://fbartos.github.io/RoBMA/reference/loo_weights.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | AIC / BIC | No | No | No | No | No | No | No |
| Reporting | Plain-text interpretation ([`interpret()`](https://fbartos.github.io/RoBMA/reference/interpret.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Sub-model summary ([`summary_models()`](https://fbartos.github.io/RoBMA/reference/summary_models.md)) | No | No | No | No | ✓ | ✓ | ✓ |
|  | Prior inspection ([`print_prior()`](https://fbartos.github.io/RoBMA/reference/print_prior.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| Extraction | Posterior draw extraction ([`as_draws()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md), [`as_draws_array()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md), [`as_draws_df()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md), [`as_draws_list()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md), [`as_draws_matrix()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md), [`as_draws_rvars()`](https://fbartos.github.io/RoBMA/reference/as_draws.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Coefficients ([`coef()`](https://fbartos.github.io/RoBMA/reference/coef.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Point-wise log-likelihood ([`logLik()`](https://fbartos.github.io/RoBMA/reference/logLik.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Sample size ([`nobs()`](https://fbartos.github.io/RoBMA/reference/nobs.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Variance-covariance matrix ([`vcov()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/vcov.html)) | No | No | No | No | No | No | No |
|  | Credible intervals ([`summary()`](https://fbartos.github.io/RoBMA/reference/summary.brma.md), [`summary_heterogeneity()`](https://fbartos.github.io/RoBMA/reference/summary_heterogeneity.brma.md), [`pooled_effect()`](https://fbartos.github.io/RoBMA/reference/pooled_effect.brma.md), [`pooled_heterogeneity()`](https://fbartos.github.io/RoBMA/reference/pooled_heterogeneity.brma.md)) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
|  | Model weights ([`weights()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/weights.html)) | No | No | No | No | No | No | No |
|  | Simulated responses ([`simulate()`](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/simulate.html)) | No | No | No | No | No | No | No |

## References

Bartoš, F., & Maier, M. (2020). *RoBMA: An R package for robust Bayesian
meta-analyses*. <https://CRAN.R-project.org/package=RoBMA>

Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor
package. *Journal of Statistical Software*, *36*(3), 1–48.
<https://www.jstatsoft.org/v36/i03/>
