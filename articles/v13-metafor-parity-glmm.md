# Generalized Linear Mixed-Effects Meta-Analysis

This vignette illustrates how to use the `RoBMA` R package ([Bartoš &
Maier, 2020](#ref-RoBMA)) when the outcome is a count rather than an
effect-size estimate. We walk through the BCG vaccine dataset from the
`metafor` package ([Viechtbauer, 2010](#ref-metafor)), showing the
[`brma.glmm()`](https://fbartos.github.io/RoBMA/reference/brma.glmm.md)
call alongside its
[`metafor::rma.glmm()`](https://wviechtb.github.io/metafor/reference/rma.glmm.html)
counterpart for log odds ratios. The interface mirrors
[`metafor::rma.glmm()`](https://wviechtb.github.io/metafor/reference/rma.glmm.html):
counts (`ai`, `bi`, `ci`, `di`) are passed in directly, and no
`escalc()` step is needed.

All
[`brma.glmm()`](https://fbartos.github.io/RoBMA/reference/brma.glmm.md)
calls keep the default prior distributions; see the [*Prior
Distributions*](https://fbartos.github.io/RoBMA/articles/v01-prior-distributions.md)
vignette for the prior definitions and customization options. For the
non-GLMM workflow showing more features (also applicable to GLMM models)
see the [*Bayesian
Meta-Analysis*](https://fbartos.github.io/RoBMA/articles/v02-bayesian-meta-analysis.md)
vignette.

[`brma.glmm()`](https://fbartos.github.io/RoBMA/reference/brma.glmm.md)
currently supports binomial outcomes (log odds ratio, `measure = "OR"`)
and Poisson outcomes (log incidence rate ratio, `measure = "IRR"`).

## Binomial Random-Effects Model

The BCG vaccine dataset from the `metadat` package contains 13
randomized trials on tuberculosis prevention. Treated and control event
counts (`tpos`, `tneg`, `cpos`, `cneg`) are passed in directly.

``` r

data("dat.bcg", package = "metadat")
head(dat.bcg)
#>   trial               author year tpos  tneg cpos  cneg ablat     alloc
#> 1     1              Aronson 1948    4   119   11   128    44    random
#> 2     2     Ferguson & Simes 1949    6   300   29   274    55    random
#> 3     3      Rosenthal et al 1960    3   228   11   209    42    random
#> 4     4    Hart & Sutherland 1977   62 13536  248 12619    52    random
#> 5     5 Frimodt-Moller et al 1973   33  5036   47  5761    13 alternate
#> 6     6      Stein & Aronson 1953  180  1361  372  1079    44 alternate
```

[`metafor::rma.glmm()`](https://wviechtb.github.io/metafor/reference/rma.glmm.html)
fits the random-effects logistic model. We use `model = "UM.FS"` to
obtain the unconditional model with fixed study-specific intercepts, the
formulation that corresponds to the binomial GLMM in
[`brma.glmm()`](https://fbartos.github.io/RoBMA/reference/brma.glmm.md)
(Model 4 in Jackson et al. ([2018](#ref-jackson2018comparison))).

``` r

fit1_metafor <- metafor::rma.glmm(
  ai = tpos, bi = tneg, ci = cpos, di = cneg, measure = "OR", 
  model = "UM.FS", data = dat.bcg
)
fit1_metafor
#> 
#> Random-Effects Model (k = 13; tau^2 estimator: ML)
#> Model Type: Unconditional Model with Fixed Study Effects
#> 
#> tau^2 (estimated amount of total heterogeneity): 0.2949
#> tau (square root of estimated tau^2 value):      0.5430
#> I^2 (total heterogeneity / total variability):   91.02%
#> H^2 (total variability / sampling variability):  11.14
#> 
#> Tests for Heterogeneity:
#> Wld(df = 12) = 163.1649, p-val < .0001
#> LRT(df = 12) = 176.9544, p-val < .0001
#> 
#> Model Results:
#> 
#> estimate      se     zval    pval    ci.lb    ci.ub      
#>  -0.7450  0.1756  -4.2434  <.0001  -1.0891  -0.4009  *** 
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The matching Bayesian fit takes the same count columns and
`measure = "OR"`.

``` r

fit1_brma <- brma.glmm(
  ai = tpos, bi = tneg, ci = cpos, di = cneg, measure = "OR",
  data = dat.bcg, seed = 1
)
```

[`summary()`](https://rdrr.io/r/base/summary.html) reports posterior
means, credible intervals, and (suppressed here) MCMC convergence
diagnostics for the pooled log odds ratio `mu` and the heterogeneity
`tau`.

``` r

summary(fit1_brma, include_mcmc_diagnostics = FALSE)
#> 
#> Bayesian Random-Effects Model (k = 13)
#> 
#> Estimates
#>       Mean    SD  0.025    0.5  0.975
#> mu  -0.719 0.183 -1.089 -0.718 -0.353
#> tau  0.589 0.141  0.361  0.569  0.922
```

The posterior mean for `mu` tracks the maximum-likelihood estimate from
the `metafor` package, with comparable heterogeneity `tau`.
[`pooled_effect()`](https://fbartos.github.io/RoBMA/reference/pooled_effect.md)
returns the same summary backtransformed to the odds-ratio scale via
`transform = "EXP"`:

``` r

pooled_effect(fit1_brma, transform = "EXP")
#> 
#> Pooled Effect Size (odds ratio)
#>     Mean Median 0.025 0.975
#> mu 0.495  0.488 0.337 0.703
#> Effect estimates are summarized on the odds ratio scale using EXP on the log odds ratio measure.
```

## Default Prior Distributions and the Auxiliary Parameter

The default prior distributions on `mu` and `tau` follow the
unit-information set-up as elsewhere in the package, scaled by the known
unit information standard deviation for the chosen `measure` (`sqrt(4)`
for `OR` and `IRR`). Beyond `mu` and `tau`, GLMM models include
estimate-specific nuisance parameters for the per-study baseline level.
For binomial outcomes (`measure = "OR"`), `logit(pi[i])` is the midpoint
of the two arm logits before applying half of the study-specific log
odds ratio to each arm. The default `prior_baserate = Beta(1, 1)` is
placed independently on each `pi[i]`. For Poisson outcomes
(`measure = "IRR"`), `phi[i]` is the midpoint of the two arm log
incidence rates before applying half of the study-specific log incidence
rate ratio to each arm. Because this midpoint log-rate is unbounded,
there is no bounded flat default analogous to `Beta(1, 1)`, so the
default `prior_lograte` is a proper data-based unit-information normal
used independently for each `phi[i]`. Thus the default
nuisance-parameter prior distributions are exchangeable in the iid
sense, but they are not hierarchical prior distributions with a learned
common baseline-rate distribution. Both defaults can be overridden via
the corresponding `prior_*` argument; see the [*Prior
Distributions*](https://fbartos.github.io/RoBMA/articles/v01-prior-distributions.md)
vignette. The Poisson model corresponds to Bagos & Nikolopoulos
([2009](#ref-bagos2009mixed)) and uses event counts and exposure times
(`x1i`, `t1i`, `x2i`, `t2i`) in place of the four-way count layout.

## Other Inference Helpers

All inference helpers available for any other
[`brma()`](https://fbartos.github.io/RoBMA/reference/brma.md) fit work
the same way for a
[`brma.glmm()`](https://fbartos.github.io/RoBMA/reference/brma.glmm.md)
fit. Meta-regression via `mods`, multilevel structures via `cluster`,
location-scale models via `scale`, posterior summaries, plots,
predictions, residuals, influence, and LOO comparisons all carry over.
See the [*Bayesian
Meta-Analysis*](https://fbartos.github.io/RoBMA/articles/v02-bayesian-meta-analysis.md)
vignette for the full walkthrough and the [*Multilevel
Meta-Analysis*](https://fbartos.github.io/RoBMA/articles/v10-metafor-parity-multilevel.md)
vignette for the multilevel extension; the only thing that changes here
is the count-based input and the auxiliary baseline-rate prior.

## References

Bagos, P. G., & Nikolopoulos, G. K. (2009). Mixed-effects poisson
regression models for meta-analysis of follow-up studies with constant
or varying durations. *The International Journal of Biostatistics*,
*5*(1), 21. <https://doi.org/10.2202/1557-4679.1168>

Bartoš, F., & Maier, M. (2020). *RoBMA: An R package for robust Bayesian
meta-analyses*. <https://CRAN.R-project.org/package=RoBMA>

Jackson, D., Law, M., Stijnen, T., Viechtbauer, W., & White, I. R.
(2018). A comparison of seven random-effects models for meta-analyses
that estimate the summary odds ratio. *Statistics in Medicine*, *37*(7),
1059–1085. <https://doi.org/10.1002/sim.7588>

Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor
package. *Journal of Statistical Software*, *36*(3), 1–48.
<https://doi.org/10.18637/jss.v036.i03>
