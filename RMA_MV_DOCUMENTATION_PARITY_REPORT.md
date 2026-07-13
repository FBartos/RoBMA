# `metafor::rma.mv()` / `RoBMA::brma.mv()` documentation parity audit

Generated: 2026-07-13 07:12 CEST

> **Post-audit engineering update (2026-07-13).** The tables below preserve the
> original audit run. A subsequent root-cause investigation fixed the recurrent
> six-outcome UN/LKJ runtime error, added a fail-fast guard for unsafe dense
> random-effect graphs, and distinguished covariance-equivalent MCMC failures
> from translation defects. See the
> [cross-package engineering follow-up](RMA_MV_ENGINEERING_FOLLOWUP.md).

## Scope and interpretation

This audit treats every literal `rma.mv()` documentation occurrence as inventory, then deduplicates numerically identical statistical models for MCMC. It covers all 25 calls in current `metafor` help and all 69 stable `metadat` calls; the complete per-topic inventory and 73-call official website inventory are linked below. Optimizer-only repeats remain listed but share one Bayesian fit.

Versions: installed `metafor` 5.0-1 for numerical references; upstream `metafor` 5.1-12 at commit `71ae8f5`; stable `metadat` 1.6-0, with the 1.7-2 development delta recorded separately. Upstream `metafor` contains 129 test-only literal calls, summarized as a feature matrix rather than mislabelled as documentation examples.

Live sources spot-checked: [`rma.mv()` reference](https://wviechtb.github.io/metafor/reference/rma.mv.html), [`metadat` reference index](https://wviechtb.github.io/metadat/), and the official metafor analyses/tips pages listed in the website inventory.

Unique-model run status: close=60, error=8, provisional=7, not fitted=4, broadly close=3.

`metafor` values are ML/REML/GLS estimates. `brma.mv()` values are posterior means 
under a broad location prior and RoBMA's unit-information-scaled SD prior, so interval 
equality is neither expected nor claimed. 
`max |difference| / metafor SE` is the coefficient-level comparison. `close` means 
at most 0.25 SE and max R-hat at most 1.10.
`metafor` component values are variances; named `brma.mv()` `sd(...)` values are 
posterior mean standard deviations and are labelled explicitly.

Shared prior objects used in the displayed syntax:

```r
weak_location_prior <- BayesTools::prior("normal", list(mean = 0, sd = 20 * prior_scale))
unit_sd_prior <- BayesTools::prior("normal", list(mean = 0, sd = prior_scale),
                                    truncation = list(lower = 0, upper = Inf))
zero_sd_prior <- BayesTools::prior("spike", list(location = 0))
```

`fixed_correlation_prior` is a model-specific `BayesTools::prior_random()` object; 
`fixed_sigma_prior` similarly spikes selected SDs after nested grouping variables 
are expanded with `interaction()`. Executable construction is in 
`tools/rma-mv-parity/mapping.R`.

## Syntax crosswalk

| `metafor` | `brma.mv()` |
|---|---|
| `random = ~ x &#124; g, struct="CS"` | `random = ~ cs(x &#124; g)` |
| `HCS`, `UN`, `ID`, `DIAG` | `hcs()`, `us()`, `id()`, `diag()` |
| `AR`, `HAR`, `CAR` | `ar1()`, `har()`, `car()` |
| `GEN`, `GDIAG` | `us(1 + x | g)`, `diag(1 + x | g)` |
| no `random` | `prior_heterogeneity = spike(0)` |
| `R` known group covariance | same `R`; currently random intercepts only |
| spatial structures | unsupported |
| `W` / arbitrary weights | unsupported |
| fixed nested `sigma2` | expand nested group IDs; spike selected block SDs |
| fixed coefficient-specific `tau2` | partial public-API parity |

## How to read failures

`provisional` means the likelihood/syntax mapping exists but the retained MCMC 
run did not satisfy coefficient or component R-hat thresholds; it is not evidence 
that the two likelihoods differ. Runtime errors and resource bounds are preserved 
verbatim. Recurrent observed cases include a six-outcome UN JAGS `Invalid parent 
values` failure, poor mixing in fixed-correlation network models and a GEN random-
slope model, and bounded high-dimensional examples. Standardizing continuous 
design columns resolved the otherwise poorly mixing age/spline meta-regressions.

## Help-example comparisons

### `metafor/man/anova.rma.Rd`

#### `help-001`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; district/school, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~1 &#124; district/school, data = dat, measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1847 | 0.1907 |
| SE / posterior comparison | 0.08456 | max standardized difference: 0.07033 |
| variance/correlation | sigma2=0.06506, 0.03274 | sd_total(random_total)=0.3188; sd(intercept &#124; school:district)=0.1871; sd(intercept &#124; district)=0.2536 |
| status | `equivalent` | close |

#### `help-002`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; district/school, data = dat, sigma2 = c(0, NA))</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = list(sigma_001 = ~1 &#124; .metafor_sigma_group_001, sigma_002 = ~1 &#124; <br>    .metafor_sigma_group_002), data = dat, measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = fixed_sigma_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1279 | 0.1277 |
| SE / posterior comparison | 0.04387 | max standardized difference: 0.004991 |
| variance/correlation | sigma2=0, 0.08844 | sd(intercept &#124; .metafor_sigma_group_001)=0; sd(intercept &#124; .metafor_sigma_group_002)=0.3010 |
| status | `equivalent` | close |

#### `help-003`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; district/school, data = dat, sigma2 = c(NA, 0))</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = list(sigma_001 = ~1 &#124; .metafor_sigma_group_001, sigma_002 = ~1 &#124; <br>    .metafor_sigma_group_002), data = dat, measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = fixed_sigma_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1960 | 0.2030 |
| SE / posterior comparison | 0.08998 | max standardized difference: 0.07692 |
| variance/correlation | sigma2=0.08285, 0 | sd(intercept &#124; .metafor_sigma_group_001)=0.3090; sd(intercept &#124; .metafor_sigma_group_002)=0 |
| status | `equivalent` | close |

#### `help-004`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, data = dat, measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = zero_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.04641 | 0.04648 |
| SE / posterior comparison | 0.009190 | max standardized difference: 0.008336 |
| variance/correlation | none | none |
| status | `equivalent` | close |

#### `help-005`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, random = ~1 &#124; study/esid, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, random = ~1 &#124; study/esid, data = dat, measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.3678 | 0.3690 |
| SE / posterior comparison | 0.09653 | max standardized difference: 0.01274 |
| variance/correlation | sigma2=0.08073, 0.1545 | sd_total(random_total)=0.5060; sd(intercept &#124; esid:study)=0.3959; sd(intercept &#124; study)=0.3050 |
| status | `equivalent` | close |

#### `help-006`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, random = ~1 &#124; study/esid, data = dat, sigma2 = c(0, NA))</code></pre> | <pre><code>brma.mv(yi = yi, V = V, random = list(sigma_001 = ~1 &#124; .metafor_sigma_group_001, sigma_002 = ~1 &#124; <br>    .metafor_sigma_group_002), data = dat, measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = fixed_sigma_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.3452 | 0.3475 |
| SE / posterior comparison | 0.06240 | max standardized difference: 0.03683 |
| variance/correlation | sigma2=0, 0.1630 | sd(intercept &#124; .metafor_sigma_group_001)=0; sd(intercept &#124; .metafor_sigma_group_002)=0.4093 |
| status | `equivalent` | close |

#### `help-007`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, random = ~1 &#124; study/esid, data = dat, sigma2 = c(NA, 0))</code></pre> | <pre><code>brma.mv(yi = yi, V = V, random = list(sigma_001 = ~1 &#124; .metafor_sigma_group_001, sigma_002 = ~1 &#124; <br>    .metafor_sigma_group_002), data = dat, measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = fixed_sigma_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.2845 | 0.2841 |
| SE / posterior comparison | 0.08013 | max standardized difference: 0.004567 |
| variance/correlation | sigma2=0.08389, 0 | sd(intercept &#124; .metafor_sigma_group_001)=0.3202; sd(intercept &#124; .metafor_sigma_group_002)=0 |
| status | `equivalent` | close |

#### `help-008`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~deltype, random = ~1 &#124; study/esid, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~deltype, random = ~1 &#124; study/esid, data = dat, measure = "SMD", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_mods = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -0.2902, 0.7062, 0.4500 | -0.2986, 0.7156, 0.4595 |
| SE / posterior comparison | 0.2083, 0.1962, 0.2098 | max standardized difference: 0.04794 |
| variance/correlation | sigma2=0.08578, 0.1292 | sd_total(random_total)=0.4850; sd(intercept &#124; esid:study)=0.3624; sd(intercept &#124; study)=0.3134 |
| status | `equivalent` | close |

### `metafor/man/confint.rma.Rd`

#### `help-009`; reuses fit `help-001`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; district/school, data = dat.konstantopoulos2011)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~1 &#124; district/school, data = dat.konstantopoulos2011, measure = "SMD", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1847 | 0.1907 |
| SE / posterior comparison | 0.08456 | max standardized difference: 0.07033 |
| variance/correlation | sigma2=0.06506, 0.03274 | sd_total(random_total)=0.3188; sd(intercept &#124; school:district)=0.1871; sd(intercept &#124; district)=0.2536 |
| status | `equivalent` | close |

#### `help-010`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~school &#124; district, data = dat.konstantopoulos2011)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~cs(school &#124; district), data = dat.konstantopoulos2011, <br>    measure = "SMD", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1847 | 0.1858 |
| SE / posterior comparison | 0.08456 | max standardized difference: 0.01232 |
| variance/correlation | tau2=0.09780; rho=0.6653 | sd(shared &#124; district)=0.3188; rho(district)=0.6270 |
| status | `equivalent` | close |

### `metafor/man/deltamethod.Rd`

#### `help-011`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + var1.var2, random = ~var1.var2 &#124; study, struct = "UN", <br>    data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + var1.var2, random = ~us(var1.var2 &#124; study), struct = "UN", <br>    data = dat, measure = "GEN", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -0.06002, -0.1423, 0.3167, 0.5671, -0.4888, -0.4750 | -- |
| SE / posterior comparison | 0.1408, 0.09167, 0.08470, 0.03668, 0.05089, 0.05058 | max standardized difference: -- |
| variance/correlation | tau2=0.1611, 0.06045, 0.04678, 0.004663, 0.01252, 0.01107; rho=0.9497, -0.6178, 0.5491, 0.04323, 0.3532, -0.5969, 0.4604, -0.04949, 0.2688, -0.9345, 0.7023, -0.1311, -0.6961, -0.08908, 0.4193 | -- |
| status | `equivalent` | error The following error was encountered while attempting to run the JAGS model:  
   Error in update.jags(object, n.iter, ...) : 
  Error in node mu__xREx__study_xRE_CORx_L_flat
Invalid parent values


 |

#### `help-012`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + var1.var2, random = ~var1.var2 &#124; study, struct = "UN", <br>    data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + var1.var2, random = ~us(var1.var2 &#124; study), struct = "UN", <br>    data = dat, measure = "GEN", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -0.07457, -0.1597, 0.3460, 0.6203, -0.5135, -0.4872 | -- |
| SE / posterior comparison | 0.1493, 0.09272, 0.09654, 0.05262, 0.06234, 0.05803 | max standardized difference: -- |
| variance/correlation | tau2=0.1747, 0.05808, 0.05945, 0.009489, 0.01663, 0.01097; rho=0.9426, -0.5875, 0.6497, -0.06217, 0.2732, -0.5611, 0.4954, -0.1297, 0.2499, -0.9122, 0.8080, -0.1899, -0.6220, 0.1025, 0.1724 | -- |
| status | `equivalent` | error The following error was encountered while attempting to run the JAGS model:  
   Error in update.jags(object, n.iter, ...) : 
  Error in node mu__xREx__study_xRE_CORx_L_flat
Invalid parent values


 |

### `metafor/man/emmprep.Rd`

#### `help-013`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~factor(type), random = ~1 &#124; study, data = dat, test = "t", <br>    dfs = "contain")</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~factor(type), random = ~1 &#124; study, data = dat, measure = "ZCOR", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_mods = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.2474, -0.1228, 0.05732 | 0.2475, -0.1223, 0.05743 |
| SE / posterior comparison | 0.01873, 0.05817, 0.05979 | max standardized difference: 0.008586 |
| variance/correlation | sigma2=0.02824 | sd(intercept &#124; study)=0.1692 |
| status | `equivalent` | close |

### `metafor/man/forest.rma.Rd`

#### `help-014`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; district/school, data = dat, slab = paste0("District ", <br>    district, ", School: ", school))</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~1 &#124; district/school, data = dat, slab = paste0("District ", <br>    district, ", School: ", school), measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1847 | 0.1861 |
| SE / posterior comparison | 0.08456 | max standardized difference: 0.01645 |
| variance/correlation | sigma2=0.06506, 0.03274 | sd_total(random_total)=0.3204; sd(intercept &#124; school:district)=0.1864; sd(intercept &#124; district)=0.2561 |
| status | `equivalent` | close |

### `metafor/man/influence.rma.mv.Rd`

#### `help-015`; reuses fit `help-001`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; district/school, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~1 &#124; district/school, data = dat, measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1847 | 0.1907 |
| SE / posterior comparison | 0.08456 | max standardized difference: 0.07033 |
| variance/correlation | sigma2=0.06506, 0.03274 | sd_total(random_total)=0.3188; sd(intercept &#124; school:district)=0.1871; sd(intercept &#124; district)=0.2536 |
| status | `equivalent` | close |

### `metafor/man/matreg.Rd`

#### `help-016`; reuses fit `help-011`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + var1.var2, random = ~var1.var2 &#124; study, struct = "UN", <br>    data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + var1.var2, random = ~us(var1.var2 &#124; study), struct = "UN", <br>    data = dat, measure = "GEN", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -0.06002, -0.1423, 0.3167, 0.5671, -0.4888, -0.4750 | -- |
| SE / posterior comparison | 0.1408, 0.09167, 0.08470, 0.03668, 0.05089, 0.05058 | max standardized difference: -- |
| variance/correlation | tau2=0.1611, 0.06045, 0.04678, 0.004663, 0.01252, 0.01107; rho=0.9497, -0.6178, 0.5491, 0.04323, 0.3532, -0.5969, 0.4604, -0.04949, 0.2688, -0.9345, 0.7023, -0.1311, -0.6961, -0.08908, 0.4193 | -- |
| status | `equivalent` | error The following error was encountered while attempting to run the JAGS model:  
   Error in update.jags(object, n.iter, ...) : 
  Error in node mu__xREx__study_xRE_CORx_L_flat
Invalid parent values


 |

#### `help-017`; reuses fit `help-012`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + var1.var2, random = ~var1.var2 &#124; study, struct = "UN", <br>    data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + var1.var2, random = ~us(var1.var2 &#124; study), struct = "UN", <br>    data = dat, measure = "GEN", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -0.07457, -0.1597, 0.3460, 0.6203, -0.5135, -0.4872 | -- |
| SE / posterior comparison | 0.1493, 0.09272, 0.09654, 0.05262, 0.06234, 0.05803 | max standardized difference: -- |
| variance/correlation | tau2=0.1747, 0.05808, 0.05945, 0.009489, 0.01663, 0.01097; rho=0.9426, -0.5875, 0.6497, -0.06217, 0.2732, -0.5611, 0.4954, -0.1297, 0.2499, -0.9122, 0.8080, -0.1899, -0.6220, 0.1025, 0.1724 | -- |
| status | `equivalent` | error The following error was encountered while attempting to run the JAGS model:  
   Error in update.jags(object, n.iter, ...) : 
  Error in node mu__xREx__study_xRE_CORx_L_flat
Invalid parent values


 |

#### `help-018`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~0 + group, random = ~group &#124; study, struct = "UN", data = dat.long, <br>    method = "ML")</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~0 + group, random = ~us(group &#124; study), struct = "UN", data = dat.long, <br>    measure = "PLO", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -4.096, -4.834 | -4.083, -4.837 |
| SE / posterior comparison | 0.4347, 0.3396 | max standardized difference: 0.02955 |
| variance/correlation | tau2=2.407, 1.431; rho=0.9467 | sd_total(allocation)=1.281; sd(intercept &#124; study)=1.676; sd(group[exp] &#124; study)=0.6695; cor(intercept,group[exp] &#124; study)=-0.6170 |
| status | `equivalent` | close |

#### `help-019`; reuses fit `help-018`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~0 + group, random = ~group &#124; study, struct = "UN", data = dat.long, <br>    method = "ML", cvvc = "varcov", control = list(nearpd = TRUE))</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~0 + group, random = ~us(group &#124; study), struct = "UN", data = dat.long, <br>    measure = "PLO", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -4.096, -4.834 | -4.083, -4.837 |
| SE / posterior comparison | 0.4347, 0.3396 | max standardized difference: 0.02955 |
| variance/correlation | tau2=2.407, 1.431; rho=0.9467 | sd_total(allocation)=1.281; sd(intercept &#124; study)=1.676; sd(group[exp] &#124; study)=0.6695; cor(intercept,group[exp] &#124; study)=-0.6170 |
| status | `equivalent` | close |

### `metafor/man/methods.confint.rma.Rd`

#### `help-020`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + outcome, random = ~outcome &#124; trial, struct = "UN", data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + outcome, random = ~us(outcome &#124; trial), struct = "UN", <br>    data = dat, measure = "MD", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -0.3392, 0.3534 | -0.3503, 0.3541 |
| SE / posterior comparison | 0.08791, 0.05885 | max standardized difference: 0.1264 |
| variance/correlation | tau2=0.03265, 0.01173; rho=0.6088 | sd_total(allocation)=0.2153; sd(intercept &#124; trial)=0.2247; sd(outcome[PD] &#124; trial)=0.1937; cor(intercept,outcome[PD] &#124; trial)=-0.4426 |
| status | `equivalent` | close |

### `metafor/man/profile.rma.Rd`

#### `help-021`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~group, random = ~group &#124; study, struct = "UN", data = dat.long)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~group, random = ~us(group &#124; study), struct = "UN", data = dat.long, <br>    measure = "PLO", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -4.096, -0.7414 | -4.119, -0.7397 |
| SE / posterior comparison | 0.4529, 0.1880 | max standardized difference: 0.05123 |
| variance/correlation | tau2=2.617, 1.549; rho=0.9450 | sd_total(allocation)=1.290; sd(intercept &#124; study)=1.687; sd(group[exp] &#124; study)=0.6743; cor(intercept,group[exp] &#124; study)=-0.5916 |
| status | `equivalent` | close |

### `metafor/man/rma.mv.Rd`

#### `help-022`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; trial, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~1 &#124; trial, data = dat, measure = "OR", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -0.7452 | -0.7436 |
| SE / posterior comparison | 0.1860 | max standardized difference: 0.008603 |
| variance/correlation | sigma2=0.3378 | sd(intercept &#124; trial)=0.6210 |
| status | `equivalent` | close |

#### `help-023`; reuses fit `help-021`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~group, random = ~group &#124; study, struct = "UN", data = dat.long)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~group, random = ~us(group &#124; study), struct = "UN", data = dat.long, <br>    measure = "PLO", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -4.096, -0.7414 | -4.119, -0.7397 |
| SE / posterior comparison | 0.4529, 0.1880 | max standardized difference: 0.05123 |
| variance/correlation | tau2=2.617, 1.549; rho=0.9450 | sd_total(allocation)=1.290; sd(intercept &#124; study)=1.687; sd(group[exp] &#124; study)=0.6743; cor(intercept,group[exp] &#124; study)=-0.5916 |
| status | `equivalent` | close |

### `metafor/man/robust.Rd`

#### `help-024`; reuses fit `help-001`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; district/school, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~1 &#124; district/school, data = dat, measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1847 | 0.1907 |
| SE / posterior comparison | 0.08456 | max standardized difference: 0.07033 |
| variance/correlation | sigma2=0.06506, 0.03274 | sd_total(random_total)=0.3188; sd(intercept &#124; school:district)=0.1871; sd(intercept &#124; district)=0.2536 |
| status | `equivalent` | close |

#### `help-025`; reuses fit `help-020`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + outcome, random = ~outcome &#124; trial, struct = "UN", data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + outcome, random = ~us(outcome &#124; trial), struct = "UN", <br>    data = dat, measure = "MD", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -0.3392, 0.3534 | -0.3503, 0.3541 |
| SE / posterior comparison | 0.08791, 0.05885 | max standardized difference: 0.1264 |
| variance/correlation | tau2=0.03265, 0.01173; rho=0.6088 | sd_total(allocation)=0.2153; sd(intercept &#124; trial)=0.2247; sd(outcome[PD] &#124; trial)=0.1937; cor(intercept,outcome[PD] &#124; trial)=-0.4426 |
| status | `equivalent` | close |

### `metadat/man/dat.assink2016.Rd`

#### `help-026`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; study/esid, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~1 &#124; study/esid, data = dat, measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.4268 | 0.4320 |
| SE / posterior comparison | 0.1184 | max standardized difference: 0.04384 |
| variance/correlation | sigma2=0.1879, 0.1120 | sd_total(random_total)=0.5651; sd(intercept &#124; esid:study)=0.3398; sd(intercept &#124; study)=0.4468 |
| status | `equivalent` | close |

#### `help-027`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; study/esid, data = dat, sigma2 = c(0, NA))</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = list(sigma_001 = ~1 &#124; .metafor_sigma_group_001, sigma_002 = ~1 &#124; <br>    .metafor_sigma_group_002), data = dat, measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = fixed_sigma_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.5630 | 0.5621 |
| SE / posterior comparison | 0.06740 | max standardized difference: 0.01328 |
| variance/correlation | sigma2=0, 0.3859 | sd(intercept &#124; .metafor_sigma_group_001)=0; sd(intercept &#124; .metafor_sigma_group_002)=0.6249 |
| status | `equivalent` | close |

#### `help-028`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; study/esid, data = dat, sigma2 = c(NA, 0))</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = list(sigma_001 = ~1 &#124; .metafor_sigma_group_001, sigma_002 = ~1 &#124; <br>    .metafor_sigma_group_002), data = dat, measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = fixed_sigma_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.4032 | 0.4169 |
| SE / posterior comparison | 0.1110 | max standardized difference: 0.1231 |
| variance/correlation | sigma2=0.1925, 0 | sd(intercept &#124; .metafor_sigma_group_001)=0.4667; sd(intercept &#124; .metafor_sigma_group_002)=0 |
| status | `equivalent` | close |

#### `help-029`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~pubstatus, random = ~1 &#124; study/esid, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~pubstatus, random = ~1 &#124; study/esid, data = dat, measure = "SMD", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_mods = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.8117, -0.4474 | 0.8321, -0.4735 |
| SE / posterior comparison | 0.3056, 0.3294 | max standardized difference: 0.07925 |
| variance/correlation | sigma2=0.1712, 0.1128 | sd_total(random_total)=0.5466; sd(intercept &#124; esid:study)=0.3408; sd(intercept &#124; study)=0.4220 |
| status | `equivalent` | close |

#### `help-030`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~year, random = ~1 &#124; study/esid, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~year, random = ~1 &#124; study/esid, data = dat, measure = "SMD", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_mods = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.4257, -0.04207 | 0.4318, -0.04283 |
| SE / posterior comparison | 0.1040, 0.01800 | max standardized difference: 0.05887 |
| variance/correlation | sigma2=0.1351, 0.1128 | sd_total(random_total)=0.5174; sd(intercept &#124; esid:study)=0.3402; sd(intercept &#124; study)=0.3836 |
| status | `equivalent` | close |

#### `help-031`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~deltype, random = ~1 &#124; study/esid, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~deltype, random = ~1 &#124; study/esid, data = dat, measure = "SMD", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_mods = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.4702, -0.7297, -0.2219 | 0.4654, -0.7313, -0.2213 |
| SE / posterior comparison | 0.1180, 0.1923, 0.1392 | max standardized difference: 0.04070 |
| variance/correlation | sigma2=0.1899, 0.08469 | sd_total(random_total)=0.5374; sd(intercept &#124; esid:study)=0.2962; sd(intercept &#124; study)=0.4443 |
| status | `equivalent` | close |

#### `help-032`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~year + deltype, random = ~1 &#124; study/esid, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~year + deltype, random = ~1 &#124; study/esid, data = dat, measure = "SMD", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_mods = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.4656, -0.03796, -0.7094, -0.2040 | 0.4635, -0.03768, -0.7063, -0.2003 |
| SE / posterior comparison | 0.1071, 0.01827, 0.1914, 0.1385 | max standardized difference: 0.02634 |
| variance/correlation | sigma2=0.1493, 0.08526 | sd_total(random_total)=0.4972; sd(intercept &#124; esid:study)=0.2979; sd(intercept &#124; study)=0.3929 |
| status | `equivalent` | close |

#### `help-033`; reuses fit `help-005`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, random = ~1 &#124; study/esid, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, random = ~1 &#124; study/esid, data = dat, measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.3678 | 0.3690 |
| SE / posterior comparison | 0.09653 | max standardized difference: 0.01274 |
| variance/correlation | sigma2=0.08073, 0.1545 | sd_total(random_total)=0.5060; sd(intercept &#124; esid:study)=0.3959; sd(intercept &#124; study)=0.3050 |
| status | `equivalent` | close |

#### `help-034`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, random = ~1 &#124; study/esid, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, random = ~1 &#124; study/esid, data = dat, measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.3618 | 0.3622 |
| SE / posterior comparison | 0.09327 | max standardized difference: 0.004308 |
| variance/correlation | sigma2=0.07045, 0.1508 | sd_total(random_total)=0.4927; sd(intercept &#124; esid:study)=0.3914; sd(intercept &#124; study)=0.2886 |
| status | `equivalent` | close |

### `metadat/man/dat.bakdash2021.Rd`

#### `help-035`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(es.z, vi.z, mods = ~1, random = ~1 &#124; SampleID/Outcome, data = dat, test = "t")</code></pre> | <pre><code>brma.mv(yi = es.z, V = vi.z, mods = ~1, random = ~1 &#124; SampleID/Outcome, data = dat, measure = "GEN", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_mods = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.2692 | 0.2694 |
| SE / posterior comparison | 0.02397 | max standardized difference: 0.007646 |
| variance/correlation | sigma2=0.02976, 0.01468 | sd_total(random_total)=0.2124; sd(intercept &#124; Outcome:SampleID)=0.1215; sd(intercept &#124; SampleID)=0.1734 |
| status | `equivalent` | close |

#### `help-036`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(es.z, vi.z, mods = ~0 + SA.measure.type, random = ~1 &#124; SampleID/Outcome, <br>    data = dat, test = "t")</code></pre> | <pre><code>brma.mv(yi = es.z, V = vi.z, mods = ~0 + SA.measure.type, random = ~1 &#124; SampleID/Outcome, <br>    data = dat, measure = "GEN", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.3781, 0.2297, 0.2161, 0.1904, 0.4345, 0.3045, 0.2954, 0.3347, 0.1582, 0.2884 | 0.3824, 0.2323, 0.2185, 0.1933, 0.4352, 0.3031, 0.2939, 0.3341, 0.1569, 0.2868 |
| SE / posterior comparison | 0.1131, 0.07322, 0.05695, 0.09800, 0.05843, 0.09960, 0.03299, 0.1074, 0.03924, 0.03991 | max standardized difference: 0.04586 |
| variance/correlation | sigma2=0.03144, 0.01219 | sd_total(random_total)=0.2089; sd(intercept &#124; Outcome:SampleID)=0.1111; sd(intercept &#124; SampleID)=0.1762 |
| status | `equivalent` | close |

#### `help-037`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(es.z, vi.z, mods = ~1, random = ~factor(Outcome) &#124; factor(SampleID), data = dat, <br>    test = "t")</code></pre> | <pre><code>brma.mv(yi = es.z, V = vi.z, mods = ~1, random = ~cs(Outcome &#124; SampleID), data = dat, measure = "GEN", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_mods = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.2692 | -- |
| SE / posterior comparison | 0.02397 | max standardized difference: -- |
| variance/correlation | tau2=0.04445; rho=0.6697 | -- |
| status | `equivalent` | error Bounded execution: brma.mv exceeded 60 minutes for the 678-level Outcome &#124; SampleID CS model. |

### `metadat/man/dat.bartos2023.Rd`

#### `help-038`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + outcome, random = ~outcome &#124; person, struct = "UN", data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + outcome, random = ~us(outcome &#124; person), struct = "UN", <br>    data = dat, measure = "PR", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.5006, 0.5097 | 0.5008, 0.5096 |
| SE / posterior comparison | 0.001028, 0.002414 | max standardized difference: 0.1749 |
| variance/correlation | tau2=0.00001278, 0.0002295; rho=-0.3182 | sd_total(allocation)=0.01148; sd(intercept &#124; person)=0.003532; sd(outcome[same] &#124; person)=0.01579; cor(intercept,outcome[same] &#124; person)=-0.2163 |
| status | `equivalent` | close; components provisional (max R-hat 1.548) |

### `metadat/man/dat.begg1989.Rd`

#### `help-039`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~trt, random = ~1 &#124; study, data = dat, method = "ML", digits = 2)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~trt, random = ~1 &#124; study, data = dat, measure = "PR", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.3221, 0.1486 | 0.3229, 0.1525 |
| SE / posterior comparison | 0.01955, 0.03658 | max standardized difference: 0.1067 |
| variance/correlation | sigma2=0.001285 | sd(intercept &#124; study)=0.04965 |
| status | `equivalent` | close |

#### `help-040`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~trt + trt:arms, random = ~1 &#124; study, data = dat, method = "ML", <br>    digits = 2)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~trt + trt:arms, random = ~1 &#124; study, data = dat, measure = "PR", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_mods = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.3010, 0.1996, 0.02611, -0.04121 | 0.3095, 0.1979, 0.01943, -0.04312 |
| SE / posterior comparison | 0.04538, 0.07565, 0.04993, 0.07424 | max standardized difference: 0.1873 |
| variance/correlation | sigma2=0.001072 | sd(intercept &#124; study)=0.05488 |
| status | `equivalent` | close |

#### `help-041`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~trt, random = list(~1 &#124; study, ~trt &#124; study), struct = "UN", <br>    tau2 = c(0, NA), rho = 0, data = dat, method = "ML", digits = 2)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~trt, random = list(~1 &#124; study, ~us(trt &#124; study)), struct = "UN", <br>    data = dat, measure = "PR", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = fixed_correlation_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.3221, 0.1486 | -- |
| SE / posterior comparison | 0.01955, 0.03658 | max standardized difference: -- |
| variance/correlation | sigma2=0.001285; tau2=0, 0.00000000000000000000003482; rho=0 | -- |
| status | `partial` -- fixed coefficient-specific variance | not fitted: partial |

#### `help-042`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~trt, random = list(~1 &#124; study, ~trt &#124; study), struct = "CS", <br>    rho = 0, data = dat, method = "ML", digits = 2)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~trt, random = list(~1 &#124; study, ~cs(trt &#124; study)), struct = "CS", <br>    data = dat, measure = "PR", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = fixed_correlation_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.3225, 0.1463 | 0.3238, 0.1545 |
| SE / posterior comparison | 0.01956, 0.03704 | max standardized difference: 0.2225 |
| variance/correlation | sigma2=0.0000000001500; tau2=0.001278; rho=0 | sd(intercept &#124; study)=0.03980; sd(shared &#124; study)=0.04081; rho(study)=0 |
| status | `equivalent` | close |

### `metadat/man/dat.berkey1998.Rd`

#### `help-043`; reuses fit `help-020`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + outcome, random = ~outcome &#124; trial, struct = "UN", data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + outcome, random = ~us(outcome &#124; trial), struct = "UN", <br>    data = dat, measure = "MD", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -0.3392, 0.3534 | -0.3503, 0.3541 |
| SE / posterior comparison | 0.08791, 0.05885 | max standardized difference: 0.1264 |
| variance/correlation | tau2=0.03265, 0.01173; rho=0.6088 | sd_total(allocation)=0.2153; sd(intercept &#124; trial)=0.2247; sd(outcome[PD] &#124; trial)=0.1937; cor(intercept,outcome[PD] &#124; trial)=-0.4426 |
| status | `equivalent` | close |

#### `help-044`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + outcome + outcome:I(year - 1983), random = ~outcome &#124; <br>    trial, struct = "UN", data = dat, method = "ML")</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + outcome + outcome:I(year - 1983), random = ~us(outcome &#124; <br>    trial), struct = "UN", data = dat, measure = "MD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -0.3351, 0.3479, -0.01083, 0.0009747 | -0.3350, 0.3621, -0.01212, 0.002171 |
| SE / posterior comparison | 0.07865, 0.05197, 0.02433, 0.01544 | max standardized difference: 0.2729 |
| variance/correlation | tau2=0.02501, 0.008041; rho=0.6587 | sd_total(allocation)=0.2569; sd(intercept &#124; trial)=0.2616; sd(outcome[PD] &#124; trial)=0.2360; cor(intercept,outcome[PD] &#124; trial)=-0.3510 |
| status | `equivalent` | broadly close |

### `metadat/man/dat.besson2016.Rd`

#### `help-045`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi ~ 1, V = vi, random = list(~1 &#124; study_ID, ~1 &#124; comp_ID), data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = list(~1 &#124; study_ID, ~1 &#124; comp_ID), data = dat, measure = "SMD", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -0.1067 | -0.1080 |
| SE / posterior comparison | 0.07117 | max standardized difference: 0.01828 |
| variance/correlation | sigma2=0.1179, 0.3698 | sd_total(random_total)=0.7086; sd(intercept &#124; study_ID)=0.3591; sd(intercept &#124; comp_ID)=0.6073 |
| status | `equivalent` | close |

### `metadat/man/dat.bornmann2007.Rd`

#### `help-046`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; study/obs, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~1 &#124; study/obs, data = dat, measure = "OR", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -0.1010 | -0.09902 |
| SE / posterior comparison | 0.04174 | max standardized difference: 0.04713 |
| variance/correlation | sigma2=0.01609, 0.003752 | sd_total(random_total)=0.1464; sd(intercept &#124; obs:study)=0.07499; sd(intercept &#124; study)=0.1219 |
| status | `equivalent` | close |

#### `help-047`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~type, random = ~1 &#124; study/obs, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~type, random = ~1 &#124; study/obs, data = dat, measure = "OR", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_mods = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -0.2010, 0.1890 | -0.1975, 0.1874 |
| SE / posterior comparison | 0.04294, 0.05639 | max standardized difference: 0.08095 |
| variance/correlation | sigma2=0.004484, 0.003547 | sd_total(random_total)=0.1018; sd(intercept &#124; obs:study)=0.06618; sd(intercept &#124; study)=0.07278 |
| status | `equivalent` | close |

### `metadat/man/dat.bourassa1996.Rd`

#### `help-048`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; study/sample, data = dat, subset = sex == "combined")</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~1 &#124; study/sample, data = dat, subset = sex == "combined", <br>    measure = "OR", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 1.271 | 1.273 |
| SE / posterior comparison | 0.1133 | max standardized difference: 0.01684 |
| variance/correlation | sigma2=0.2097, 0.1985 | sd_total(random_total)=0.6568; sd(intercept &#124; sample:study)=0.4557; sd(intercept &#124; study)=0.4409 |
| status | `equivalent` | close |

#### `help-049`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; study/sample/sex, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~1 &#124; study/sample/sex, data = dat, measure = "OR", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 1.268 | 1.271 |
| SE / posterior comparison | 0.1130 | max standardized difference: 0.02645 |
| variance/correlation | sigma2=0.1828, 0.2259, 0.0000000009298 | sd_total(random_total)=0.6523; sd(intercept &#124; sex:sample:study)=0.2665; sd(intercept &#124; sample:study)=0.3697; sd(intercept &#124; study)=0.4151 |
| status | `equivalent` | close |

### `metadat/man/dat.craft2003.Rd`

#### `help-050`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + var1.var2, random = ~var1.var2 &#124; study, struct = "UN", <br>    data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + var1.var2, random = ~us(var1.var2 &#124; study), struct = "UN", <br>    data = dat, measure = "GEN", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.5671, -0.4888, -0.06002, -0.4750, -0.1423, 0.3167 | -- |
| SE / posterior comparison | 0.03668, 0.05089, 0.1408, 0.05058, 0.09167, 0.08470 | max standardized difference: -- |
| variance/correlation | tau2=0.004663, 0.01252, 0.1611, 0.01107, 0.06045, 0.04678; rho=-0.6961, 0.5491, -0.08908, 0.4604, -0.9345, 0.04323, 0.4193, -0.04949, 0.7023, 0.3532, 0.9497, -0.6178, 0.2688, -0.1311, -0.5969 | -- |
| status | `equivalent` | error The following error was encountered while attempting to run the JAGS model:  
   Error in update.jags(object, n.iter, ...) : 
  Error in node mu__xREx__study_xRE_CORx_L_flat
Invalid parent values


 |

### `metadat/man/dat.crede2010.Rd`

#### `help-051`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; studyid/sampleid, data = dat, subset = criterion == <br>    "grade")</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~1 &#124; studyid/sampleid, data = dat, subset = criterion == <br>    "grade", measure = "ZCOR", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.4798 | 0.4771 |
| SE / posterior comparison | 0.03305 | max standardized difference: 0.08223 |
| variance/correlation | sigma2=0.03760, 0.01586 | sd_total(random_total)=0.2330; sd(intercept &#124; sampleid:studyid)=0.1453; sd(intercept &#124; studyid)=0.1742 |
| status | `equivalent` | close |

#### `help-052`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; studyid/sampleid, data = dat, subset = criterion == <br>    "gpa")</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~1 &#124; studyid/sampleid, data = dat, subset = criterion == <br>    "gpa", measure = "ZCOR", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.3611 | 0.3666 |
| SE / posterior comparison | 0.05381 | max standardized difference: 0.1036 |
| variance/correlation | sigma2=0.05722, 0.006329 | sd_total(random_total)=0.2483; sd(intercept &#124; sampleid:studyid)=0.1400; sd(intercept &#124; studyid)=0.1911 |
| status | `equivalent` | close |

### `metadat/man/dat.dagostino1998.Rd`

#### `help-053`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + outcome, data = dat, subset = outcome %in% c("rnic1", <br>    "rnic2"))</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + outcome, data = dat, subset = outcome %in% c("rnic1", <br>    "rnic2"), measure = "OR", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = zero_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.2315, 0.2528 | 0.2313, 0.2520 |
| SE / posterior comparison | 0.06296, 0.06335 | max standardized difference: 0.01226 |
| variance/correlation | none | none |
| status | `equivalent` | close |

### `metadat/man/dat.fine1993.Rd`

#### `help-054`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, Vc, mods = ~0 + factor(time), method = "FE", data = dat.long)</code></pre> | <pre><code>brma.mv(yi = yi, V = Vc, mods = ~0 + factor(time), data = dat.long, measure = "GEN", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_mods = weak_location_prior, prior_heterogeneity = zero_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.2511, 0.4282, 0.3432, 0.2813 | 0.2493, 0.4257, 0.3398, 0.2771 |
| SE / posterior comparison | 0.1431, 0.1144, 0.1341, 0.1382 | max standardized difference: 0.03028 |
| variance/correlation | none | none |
| status | `equivalent` | close |

#### `help-055`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, Vc, mods = ~0 + factor(time), random = ~time &#124; study, struct = "HAR", <br>    data = dat.long, control = list(optimizer = "hjk"))</code></pre> | <pre><code>brma.mv(yi = yi, V = Vc, mods = ~0 + factor(time), random = ~har(time &#124; study), struct = "HAR", <br>    data = dat.long, measure = "GEN", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.2536, 0.4078, 0.3655, 0.3514 | 0.2589, 0.4264, 0.3551, 0.3059 |
| SE / posterior comparison | 0.1431, 0.1153, 0.1450, 0.1871 | max standardized difference: 0.2429 |
| variance/correlation | tau2=0, 0.0002293, 0.02943, 0.1967; rho=1.000 | sd_total(allocation)=0.1407; sd(time[1] &#124; study)=0.1202; sd(time[2] &#124; study)=0.1075; sd(time[3] &#124; study)=0.1253; sd(time[4] &#124; study)=0.1603; rho(study)=0.5174 |
| status | `equivalent` | close |

### `metadat/man/dat.hasselblad1998.Rd`

#### `help-056`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~0 + factor(study) + trt, random = ~trt &#124; study, rho = 1/2, <br>    data = dat, btt = "trt")</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~0 + factor(study) + trt, random = ~cs(trt &#124; study), data = dat, <br>    measure = "PLO", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = fixed_correlation_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -1.393, -2.730, -3.722, -2.271, -2.391, -0.9698, -2.086, -1.959, -2.349, -1.806, -2.262, -3.564, -2.345, -2.062, -2.458, -2.144, -2.681, -2.407, -1.444, -2.404, -2.436, -2.856, -2.327, -3.089, 0.3888, 0.6864, 0.8438 | -1.305, -2.728, -3.688, -2.331, -2.423, -0.9729, -2.113, -1.965, -2.333, -1.832, -2.230, -3.534, -2.387, -2.008, -2.477, -2.165, -2.686, -2.404, -1.479, -2.439, -2.478, -2.816, -2.313, -3.147, 0.4013, 0.6873, 0.8522 |
| SE / posterior comparison | 0.5820, 0.5842, 0.6681, 0.5829, 0.7374, 0.6421, 0.6369, 0.6674, 0.6032, 0.6296, 0.5978, 0.6139, 0.5836, 0.6492, 0.6868, 0.6707, 0.6448, 0.5964, 0.8108, 0.6331, 0.5817, 0.5991, 0.5838, 0.5847, 0.3221, 0.1904, 0.3641 | max standardized difference: 0.1502 |
| variance/correlation | tau2=0.4324; rho=0.5000 | sd(shared &#124; study)=0.7060; rho(study)=0.5000 |
| status | `equivalent` | close |

#### `help-057`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + self_help + ind_counseling + grp_counseling, random = ~comp &#124; <br>    study, rho = 1/2, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + self_help + ind_counseling + grp_counseling, random = ~cs(comp &#124; <br>    study), data = dat, measure = "OR", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = fixed_correlation_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.3888, 0.6864, 0.8438 | 0.4090, 0.6927, 0.8668 |
| SE / posterior comparison | 0.3221, 0.1904, 0.3641 | max standardized difference: 0.06325 |
| variance/correlation | tau2=0.4324; rho=0.5000 | sd(shared &#124; study)=0.6947; rho(study)=0.5000 |
| status | `equivalent` | close |

#### `help-058`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + self_help + ind_counseling + grp_counseling, random = list(~comp &#124; <br>    study, ~comp &#124; design), rho = 1/2, phi = 1/2, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + self_help + ind_counseling + grp_counseling, random = list(~cs(comp &#124; <br>    study), ~cs(comp &#124; design)), data = dat, measure = "OR", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_mods = weak_location_prior, prior_heterogeneity = fixed_correlation_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.3888, 0.6864, 0.8438 | 0.4720, 0.7643, 0.9391 |
| SE / posterior comparison | 0.3221, 0.1904, 0.3641 | max standardized difference: 0.4089 |
| variance/correlation | tau2=0.4324; rho=0.5000; gamma2=0.0000000004664; phi=0.5000 | sd(shared &#124; study)=0.7059; rho(study)=0.5000; sd(shared &#124; design)=0.3581; rho(design)=0.5000 |
| status | `equivalent` | broadly close |

### `metadat/man/dat.ishak2007.Rd`

#### `help-059`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + factor(time), random = ~time &#124; study, struct = "HAR", <br>    data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + factor(time), random = ~har(time &#124; study), struct = "HAR", <br>    data = dat, measure = "GEN", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -25.90, -27.46, -28.66, -26.49 | -25.83, -27.40, -28.70, -26.57 |
| SE / posterior comparison | 1.012, 1.141, 1.032, 1.382 | max standardized difference: 0.07338 |
| variance/correlation | tau2=22.73, 33.73, 26.14, 31.18; rho=0.8832 | sd_total(allocation)=5.147; sd(time[1] &#124; study)=4.553; sd(time[2] &#124; study)=5.526; sd(time[3] &#124; study)=4.958; sd(time[4] &#124; study)=5.309; rho(study)=0.8302 |
| status | `equivalent` | close |

### `metadat/man/dat.kalaian1996.Rd`

#### `help-060`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + outcome, random = ~outcome &#124; study, struct = "UN", data = dat, <br>    digits = 3)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + outcome, random = ~us(outcome &#124; study), struct = "UN", <br>    data = dat, measure = "GEN", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1379, 0.1168 | 0.1373, 0.1185 |
| SE / posterior comparison | 0.04338, 0.03374 | max standardized difference: 0.05221 |
| variance/correlation | tau2=0.01215, 0.002574; rho=-1.000 | sd_total(allocation)=0.1049; sd(intercept &#124; study)=0.07990; sd(outcome[verbal] &#124; study)=0.1219; cor(intercept,outcome[verbal] &#124; study)=-0.6734 |
| status | `equivalent` | close |

#### `help-061`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + outcome + outcome:loghrs, random = ~outcome &#124; study, struct = "UN", <br>    data = dat, digits = 3)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + outcome + outcome:loghrs, random = ~us(outcome &#124; study), <br>    struct = "UN", data = dat, measure = "GEN", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1017, 0.1098, 0.1694, 0.04900 | 0.1020, 0.1156, 0.1598, 0.05532 |
| SE / posterior comparison | 0.04870, 0.03467, 0.07254, 0.04587 | max standardized difference: 0.1674 |
| variance/correlation | tau2=0.01522, 0.001393; rho=-1.000 | sd_total(allocation)=0.1226; sd(intercept &#124; study)=0.1017; sd(outcome[verbal] &#124; study)=0.1374; cor(intercept,outcome[verbal] &#124; study)=-0.7520 |
| status | `equivalent` | close |

#### `help-062`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + outcome + outcome:loghrs, random = ~outcome &#124; study, struct = "UN", <br>    tau2 = c(NA, 0), data = dat, digits = 3)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + outcome + outcome:loghrs, random = ~us(outcome &#124; study), <br>    struct = "UN", data = dat, measure = "GEN", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_mods = weak_location_prior, prior_heterogeneity = fixed_correlation_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1070, 0.1140, 0.1874, 0.06022 | -- |
| SE / posterior comparison | 0.05088, 0.03358, 0.07610, 0.04438 | max standardized difference: -- |
| variance/correlation | tau2=0.02115, 0; rho=0 | -- |
| status | `partial` -- fixed coefficient-specific variance | not fitted: partial |

### `metadat/man/dat.kearon1998.Rd`

#### `help-063`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~0 + group, random = ~group &#124; study, struct = "UN", data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~0 + group, random = ~us(group &#124; study), struct = "UN", data = dat, <br>    measure = "PLO", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -0.04982, 3.009 | -0.05797, 3.037 |
| SE / posterior comparison | 0.2027, 0.2582 | max standardized difference: 0.1089 |
| variance/correlation | tau2=0.4231, 0.5408; rho=-0.4547 | sd_total(allocation)=0.9463; sd(intercept &#124; study)=0.6895; sd(group[specificity] &#124; study)=1.131; cor(intercept,group[specificity] &#124; study)=-0.6322 |
| status | `equivalent` | close |

### `metadat/man/dat.knapp2017.Rd`

#### `help-064`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; study/comp, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~1 &#124; study/comp, data = dat, measure = "GEN", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.6801 | 0.6793 |
| SE / posterior comparison | 0.05597 | max standardized difference: 0.01453 |
| variance/correlation | sigma2=0.03833, 0.02627 | sd_total(random_total)=0.2571; sd(intercept &#124; comp:study)=0.1632; sd(intercept &#124; study)=0.1862 |
| status | `equivalent` | close |

#### `help-065`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, random = ~1 &#124; study/comp, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, random = ~1 &#124; study/comp, data = dat, measure = "GEN", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.6696 | 0.6770 |
| SE / posterior comparison | 0.05455 | max standardized difference: 0.1358 |
| variance/correlation | sigma2=0.01018, 0.05568 | sd_total(random_total)=0.2799; sd(intercept &#124; comp:study)=0.2234; sd(intercept &#124; study)=0.1544 |
| status | `equivalent` | close |

#### `help-066`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~difficulty, random = ~1 &#124; study/comp, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~difficulty, random = ~1 &#124; study/comp, data = dat, measure = "GEN", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_mods = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1803, 0.1241 | 0.1968, 0.1186 |
| SE / posterior comparison | 0.1682, 0.03920 | max standardized difference: 0.1383 |
| variance/correlation | sigma2=0.04647, 0.02752 | sd_total(random_total)=0.2660; sd(intercept &#124; comp:study)=0.1723; sd(intercept &#124; study)=0.1893 |
| status | `equivalent` | close |

### `metadat/man/dat.konstantopoulos2011.Rd`

#### `help-067`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; study, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~1 &#124; study, data = dat, measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1279 | 0.1277 |
| SE / posterior comparison | 0.04387 | max standardized difference: 0.006100 |
| variance/correlation | sigma2=0.08844 | sd(intercept &#124; study)=0.3021 |
| status | `equivalent` | close |

#### `help-068`; reuses fit `help-001`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; district/school, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~1 &#124; district/school, data = dat, measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1847 | 0.1907 |
| SE / posterior comparison | 0.08456 | max standardized difference: 0.07033 |
| variance/correlation | sigma2=0.06506, 0.03274 | sd_total(random_total)=0.3188; sd(intercept &#124; school:district)=0.1871; sd(intercept &#124; district)=0.2536 |
| status | `equivalent` | close |

#### `help-069`; reuses fit `help-010`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~school &#124; district, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~cs(school &#124; district), data = dat, measure = "SMD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1847 | 0.1858 |
| SE / posterior comparison | 0.08456 | max standardized difference: 0.01232 |
| variance/correlation | tau2=0.09780; rho=0.6653 | sd(shared &#124; district)=0.3188; rho(district)=0.6270 |
| status | `equivalent` | close |

### `metadat/man/dat.lim2014.Rd`

#### `help-070`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = list(~1 &#124; article, ~1 &#124; esid, ~1 &#124; species, ~1 &#124; species.phy), <br>    R = list(species.phy = A), data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = list(~1 &#124; article, ~1 &#124; esid, ~1 &#124; species, ~1 &#124; species.phy), <br>    R = list(species.phy = A), data = dat, measure = "ZCOR", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -0.1564 | -0.1589 |
| SE / posterior comparison | 0.1277 | max standardized difference: 0.01999 |
| variance/correlation | sigma2=0.1387, 0.009256, 0.0000000001833, 0.05721 | sd_total(random_total)=0.4519; sd(intercept &#124; article)=0.3337; sd(intercept &#124; esid)=0.1185; sd(intercept &#124; species)=0.1242 |
| status | `equivalent` | close |

### `metadat/man/dat.lopez2019.Rd`

#### `help-071`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, data = dat, mods = ~0 + No.treatment + Wait.list + Placebo + F2F.CBT + <br>    Hybrid.CBT + Multimedia.CBT, random = ~comp &#124; study, rho = 1/2)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, data = dat, mods = ~0 + No.treatment + Wait.list + Placebo + F2F.CBT + <br>    Hybrid.CBT + Multimedia.CBT, random = ~cs(comp &#124; study), measure = "GEN", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_mods = weak_location_prior, prior_heterogeneity = fixed_correlation_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1975, 0.7030, -0.3699, -1.118, -1.078, -0.6017 | -1.912, -0.8193, -0.9516, -1.605, -1.443, -4.665 |
| SE / posterior comparison | 0.5905, 0.3471, 0.4623, 0.2730, 0.5190, 0.3399 | max standardized difference: 11.95 |
| variance/correlation | tau2=1.626; rho=0.5000 | sd(shared &#124; study)=1.237; rho(study)=0.5000 |
| status | `equivalent` | coefficients provisional (coefficient R-hat 10.042) |

#### `help-072`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, data = dat, mods = ~0 + No.treatment + Wait.list + Placebo + F2F.CBT + <br>    Hybrid.CBT + Multimedia.CBT, random = list(~comp &#124; study, ~comp &#124; design), rho = 1/2, <br>    phi = 1/2, control = list(optimizer = "BFGS"))</code></pre> | <pre><code>brma.mv(yi = yi, V = V, data = dat, mods = ~0 + No.treatment + Wait.list + Placebo + F2F.CBT + <br>    Hybrid.CBT + Multimedia.CBT, random = list(~cs(comp &#124; study), ~cs(comp &#124; design)), measure = "GEN", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_mods = weak_location_prior, <br>    prior_heterogeneity = fixed_correlation_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1975, 0.7030, -0.3699, -1.118, -1.078, -0.6017 | 5.564, -0.1596, 1.611, 0.4507, 2.982, -4.864 |
| SE / posterior comparison | 0.5905, 0.3471, 0.4623, 0.2730, 0.5190, 0.3399 | max standardized difference: 12.54 |
| variance/correlation | tau2=1.626; rho=0.5000; gamma2=0.00002913; phi=0.5000 | sd(shared &#124; study)=1.799; rho(study)=0.5000; sd(shared &#124; design)=0.9229; rho(design)=0.5000 |
| status | `equivalent` | coefficients provisional (coefficient R-hat 22.486) |

#### `help-073`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, data = dat, mods = ~0 + No.treatment + Wait.list + Placebo + CBT, <br>    random = ~comp &#124; study, rho = 1/2)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, data = dat, mods = ~0 + No.treatment + Wait.list + Placebo + CBT, <br>    random = ~cs(comp &#124; study), measure = "GEN", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_mods = weak_location_prior, prior_heterogeneity = fixed_correlation_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.3552, 0.7190, -0.2746, -0.9599 | 0.5996, 1.060, 0.1173, -0.6917 |
| SE / posterior comparison | 0.5792, 0.3447, 0.4584, 0.2454 | max standardized difference: 1.093 |
| variance/correlation | tau2=1.633; rho=0.5000 | sd(shared &#124; study)=1.330; rho(study)=0.5000 |
| status | `equivalent` | coefficients provisional (coefficient R-hat 6.100) |

#### `help-074`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, data = dat, mods = ~0 + No.treatment + Wait.list + Placebo + CBT + <br>    CBT:(multi + cog + ba + psed + home + prob + soc + relax + goal + final), random = ~comp &#124; <br>    study, rho = 1/2)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, data = dat, mods = ~0 + No.treatment + Wait.list + Placebo + CBT + <br>    CBT:(multi + cog + ba + psed + home + prob + soc + relax + goal + final), random = ~cs(comp &#124; <br>    study), measure = "GEN", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = fixed_correlation_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.2712, 0.6238, -0.5293, -1.491, 0.4097, -0.3086, 0.4766, 0.3588, -0.08963, 0.1888, 0.5758, 0.07352, 0.4226, -0.5714 | 9.103, 8.640, 8.929, -1.010, 11.27, 0.1408, 2.125, 8.526, 1.904, 2.729, 0.006335, -1.912, -6.749, 1.145 |
| SE / posterior comparison | 0.6131, 0.3774, 0.4815, 0.3407, 0.3685, 0.2460, 0.2590, 0.3486, 0.3705, 0.2408, 0.4485, 0.2422, 0.5676, 0.6091 | max standardized difference: 29.46 |
| variance/correlation | tau2=1.611; rho=0.5000 | sd(shared &#124; study)=1.702; rho(study)=0.5000 |
| status | `equivalent` | coefficients provisional (coefficient R-hat 21.242) |

### `metadat/man/dat.maire2019.Rd`

#### `help-075`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(s1, vars1, random = ~1 &#124; site_station, data = dat)</code></pre> | <pre><code>brma.mv(yi = s1, V = vars1, random = ~1 &#124; site_station, data = dat, measure = "GEN", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 187.6 | 186.9 |
| SE / posterior comparison | 20.04 | max standardized difference: 0.03365 |
| variance/correlation | sigma2=10522. | sd(intercept &#124; site_station)=104.7 |
| status | `equivalent` | close |

#### `help-076`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(s1, vars1, random = ~site_station &#124; const, struct = "SPGAU", data = dat, <br>    dist = list(dmat), control = list(rho.init = 10))</code></pre> | <pre><code>brma.mv(yi = s1, V = vars1, random = ~site_station &#124; const, struct = "SPGAU", data = dat, <br>    measure = "GEN", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 176.6 | -- |
| SE / posterior comparison | 26.40 | max standardized difference: -- |
| variance/correlation | tau2=    8946.; rho=15.06 | -- |
| status | `unsupported` -- spatial covariance structure | not fitted: unsupported |

#### `help-077`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(s1, vars1, random = list(~1 &#124; site/station, ~site_station &#124; const), struct = "SPGAU", <br>    data = dat, dist = list(dmat), control = list(rho.init = 10))</code></pre> | <pre><code>brma.mv(yi = s1, V = vars1, random = list(~1 &#124; site/station, ~site_station &#124; const), struct = "SPGAU", <br>    data = dat, measure = "GEN", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 182.1 | -- |
| SE / posterior comparison | 31.28 | max standardized difference: -- |
| variance/correlation | sigma2=7158., 0.00001511; tau2=    3426.; rho=12.88 | -- |
| status | `unsupported` -- spatial covariance structure | not fitted: unsupported |

### `metadat/man/dat.mccurdy2020.Rd`

#### `help-078`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~condition, random = list(~1 &#124; article/experiment/sample/id, <br>    ~1 &#124; pairing), data = dat, sparse = TRUE, digits = 3)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~condition, random = list(~1 &#124; article/experiment/sample/id, <br>    ~1 &#124; pairing), data = dat, measure = "GEN", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.4785, 0.1021 | -- |
| SE / posterior comparison | 0.01572, 0.004246 | max standardized difference: -- |
| variance/correlation | sigma2=0.02186, 0.006042, 0.00000000002129, 0.006375, 0.01651 | -- |
| status | `equivalent` | error Bounded execution: brma.mv exceeded 60 minutes for the 1,653-effect sparse nested/pairing model. |

### `metadat/man/dat.moura2021.Rd`

#### `help-079`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = list(~1 &#124; study.id, ~1 &#124; effect.size.id, ~1 &#124; species.id, <br>    ~1 &#124; species.id.phy), R = list(species.id.phy = A), data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = list(~1 &#124; study.id, ~1 &#124; effect.size.id, ~1 &#124; species.id, <br>    ~1 &#124; species.id.phy), R = list(species.id.phy = A), data = dat, measure = "ZCOR", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.3682 | -- |
| SE / posterior comparison | 0.1300 | max standardized difference: -- |
| variance/correlation | sigma2=0.01916, 0.01445, 0.05566, 0.05122 | -- |
| status | `equivalent` | error Bounded execution after corrected known-R mapping: brma.mv exceeded 20 minutes for the 1,828-effect four-level phylogenetic model. |

#### `help-080`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~spatially.pooled * temporally.pooled, random = list(~1 &#124; <br>    study.id, ~1 &#124; effect.size.id, ~1 &#124; species.id, ~1 &#124; species.id.phy), R = list(species.id.phy = A), <br>    data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~spatially.pooled * temporally.pooled, random = list(~1 &#124; <br>    study.id, ~1 &#124; effect.size.id, ~1 &#124; species.id, ~1 &#124; species.id.phy), R = list(species.id.phy = A), <br>    data = dat, measure = "ZCOR", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.3457, 0.08099, 0.05990, -0.07286 | -- |
| SE / posterior comparison | 0.1327, 0.03898, 0.02690, 0.04523 | max standardized difference: -- |
| variance/correlation | sigma2=0.01976, 0.01444, 0.05255, 0.05323 | -- |
| status | `equivalent` | error Bounded execution after corrected known-R mapping: brma.mv exceeded 25 minutes for the 1,828-effect four-level phylogenetic model. |

### `metadat/man/dat.nakagawa2007.Rd`

#### `help-081`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = ~1 &#124; StudyID, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = ~1 &#124; StudyID, data = dat, measure = "ZCOR", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.5205 | 0.5174 |
| SE / posterior comparison | 0.09636 | max standardized difference: 0.03290 |
| variance/correlation | sigma2=0.05123 | sd(intercept &#124; StudyID)=0.2164 |
| status | `equivalent` | close |

### `metadat/man/dat.obrien2003.Rd`

#### `help-082`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~bmicent, random = ~1 &#124; study/grp, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~bmicent, random = ~1 &#124; study/grp, data = dat, measure = "PR", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_mods = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 5.013, 0.5749 | 5.030, 0.5787 |
| SE / posterior comparison | 0.8068, 0.07612 | max standardized difference: 0.05000 |
| variance/correlation | sigma2=7.271, 2.519 | sd_total(random_total)=3.251; sd(intercept &#124; grp:study)=1.693; sd(intercept &#124; study)=2.745 |
| status | `equivalent` | close |

#### `help-083`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~bmicent, random = ~bmicent &#124; study, struct = "GEN", data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~bmicent, random = ~us(1 + bmicent &#124; study), struct = "GEN", <br>    data = dat, measure = "PR", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 4.734, 0.4887 | 5.092, -8.827 |
| SE / posterior comparison | 0.8085, 0.09133 | max standardized difference: 102.0 |
| variance/correlation | tau2=8.213, 0.09035; rho=0.8723 | sd_total(allocation)=7.386; sd(intercept &#124; study)=5.990; sd(bmicent &#124; study)=2.172; cor(intercept,bmicent &#124; study)=0.02169 |
| status | `equivalent` | coefficients provisional (coefficient R-hat 34.274) |

#### `help-084`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~rcs(bmicent, 4), random = ~1 &#124; study/grp, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~rcs(bmicent, 4), random = ~1 &#124; study/grp, data = dat, measure = "PR", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_mods = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 4.273, 0.3142, -0.05599, 1.110 | 4.228, 0.3060, -0.03849, 1.088 |
| SE / posterior comparison | 1.557, 0.3432, 0.8852, 2.203 | max standardized difference: 0.02871 |
| variance/correlation | sigma2=6.141, 2.479 | sd_total(random_total)=3.038; sd(intercept &#124; grp:study)=1.692; sd(intercept &#124; study)=2.489 |
| status | `equivalent` | close |

### `metadat/man/dat.pagliaro1992.Rd`

#### `help-085`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~0 + factor(study) + trt, random = ~trt &#124; study, rho = 1/2, <br>    data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~0 + factor(study) + trt, random = ~cs(trt &#124; study), data = dat, <br>    measure = "PLO", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = fixed_correlation_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -1.107, -1.095, -0.8816, -0.9237, -1.643, -1.317, -0.6125, -1.739, -0.9195, -1.396, -0.4776, -0.6253, -0.4554, -0.7814, 0.04556, -0.8312, -0.2181, -0.7510, -0.4790, -2.108, -0.1279, -1.150, -0.2472, 0.2882, -1.381, -1.218, -0.7163, -0.5839 | -1.125, -1.153, -0.8262, -0.9194, -1.607, -1.347, -0.5944, -1.694, -0.8968, -1.412, -0.4370, -0.6440, -0.4145, -0.7909, 0.04584, -0.8363, -0.2864, -0.8324, -0.5043, -2.076, -0.1920, -1.155, -0.2652, 0.2897, -1.420, -1.216, -0.7272, -0.5572 |
| SE / posterior comparison | 0.8684, 0.8483, 0.9712, 0.8980, 0.9724, 0.9221, 0.8997, 0.9479, 0.9816, 1.075, 0.9305, 0.9099, 0.9463, 0.9455, 0.9027, 0.9032, 0.8941, 0.8961, 0.9043, 1.076, 0.9071, 0.8861, 0.8957, 0.9746, 0.9557, 0.9737, 0.3871, 0.2684 | max standardized difference: 0.09973 |
| variance/correlation | tau2=0.9913; rho=0.5000 | sd(shared &#124; study)=0.9996; rho(study)=0.5000 |
| status | `equivalent` | close |

### `metadat/man/dat.senn2013.Rd`

#### `help-086`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~0 + study + treatment, data = dat, slab = paste0(study, treatment))</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~0 + study + treatment, data = dat, slab = paste0(study, <br>    treatment), measure = "MN", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = zero_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 1.092, 0.05602, -0.2775, 0.07089, -0.1932, 7.881, 0.07000, -0.5693, 8.620, 0.2700, 0.4319, 0.6071, 0.6919, 0.06423, 0.05122, 0.01585, 0.2892, 0.1220, 0.1352, -0.02407, 8.601, 0.4165, -1.388, 0.2502, 0.1396, -0.5819, -0.8273, -0.9051, -1.114, -0.9438, -1.066, -1.202, -0.5700, -0.4394, -0.7000 | 1.093, 0.05641, -0.2768, 0.07110, -0.1915, 7.883, 0.07048, -0.5637, 8.623, 0.2695, 0.4338, 0.6099, 0.6919, 0.06403, 0.05309, 0.01739, 0.2913, 0.1224, 0.1354, -0.02429, 8.602, 0.4168, -1.385, 0.2515, 0.1417, -0.5803, -0.8277, -0.9050, -1.116, -0.9458, -1.069, -1.203, -0.5683, -0.4402, -0.7004 |
| SE / posterior comparison | 0.08662, 0.05476, 0.1047, 0.07560, 0.07674, 0.1041, 0.09000, 0.2189, 0.07235, 0.09130, 0.1117, 0.1353, 0.1550, 0.05835, 0.1235, 0.06673, 0.06179, 0.1051, 0.1100, 0.1156, 0.1849, 0.09543, 0.1402, 0.06171, 0.1196, 0.08439, 0.1085, 0.1271, 0.05961, 0.1269, 0.07584, 0.04766, 0.1291, 0.09146, 0.1273 | max standardized difference: 0.04198 |
| variance/correlation | none | none |
| status | `equivalent` | close |

#### `help-087`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, mods = ~0 + study + treatment, random = ~treatment &#124; study, rho = 1/2, <br>    data = dat, slab = paste0(study, treatment))</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, mods = ~0 + study + treatment, random = ~cs(treatment &#124; study), <br>    data = dat, slab = paste0(study, treatment), measure = "MN", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_mods = weak_location_prior, prior_heterogeneity = fixed_correlation_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 1.088, 0.06719, -0.2770, 0.08686, -0.1858, 7.931, 0.07000, -0.5534, 8.659, 0.2700, 0.4349, 0.7039, 0.6949, 0.09056, 0.07822, 0.01801, 0.2400, 0.04707, 0.1287, -0.008130, 8.515, 0.3400, -1.377, 0.2668, 0.1614, -0.5392, -0.8414, -0.7369, -1.128, -0.9499, -1.129, -1.234, -0.5700, -0.4175, -0.7000 | 3.153, -5.054, 9.046, 0.6415, -0.7018, 8.769, -0.7590, -0.6700, 5.844, 0.8262, 0.5799, -0.4777, 0.7415, 0.2735, -0.04600, -1.779, -6.118, -1.830, 1.861, -0.05778, 7.385, 11.06, -0.7085, -4.009, 0.2810, -0.7004, -1.131, 0.8330, -1.361, -0.4525, -2.034, -1.481, -0.3784, -0.3115, 3.902 |
| SE / posterior comparison | 0.3221, 0.2850, 0.3118, 0.2896, 0.2924, 0.3197, 0.3287, 0.3567, 0.3139, 0.3290, 0.3100, 0.3228, 0.3280, 0.2857, 0.3180, 0.2988, 0.2887, 0.3172, 0.3461, 0.3026, 0.3528, 0.3210, 0.3168, 0.2863, 0.3161, 0.2926, 0.2384, 0.2776, 0.1494, 0.2253, 0.2119, 0.1235, 0.3414, 0.2326, 0.3408 | max standardized difference: 33.39 |
| variance/correlation | tau2=0.09991; rho=0.5000 | sd(shared &#124; study)=2.030; rho(study)=0.5000 |
| status | `equivalent` | coefficients provisional (coefficient R-hat 8.422) |

#### `help-088`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + acarbose + benfluorex + metformin + miglitol + pioglitazone + <br>    rosiglitazone + sitagliptin + sulfonylurea + vildagliptin, random = ~comp &#124; study, rho = 1/2, <br>    data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + acarbose + benfluorex + metformin + miglitol + pioglitazone + <br>    rosiglitazone + sitagliptin + sulfonylurea + vildagliptin, random = ~cs(comp &#124; study), <br>    data = dat, measure = "MD", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = fixed_correlation_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -0.8414, -0.7369, -1.128, -0.9499, -1.129, -1.234, -0.5700, -0.4175, -0.7000 | -0.8419, -0.7381, -1.135, -0.9441, -1.140, -1.238, -0.5511, -0.4214, -0.7098 |
| SE / posterior comparison | 0.2384, 0.2776, 0.1494, 0.2253, 0.2119, 0.1235, 0.3414, 0.2326, 0.3408 | max standardized difference: 0.05549 |
| variance/correlation | tau2=0.09991; rho=0.5000 | sd(shared &#124; study)=0.3326; rho(study)=0.5000 |
| status | `equivalent` | close |

#### `help-089`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~0 + acarbose + benfluorex + metformin + miglitol + pioglitazone + <br>    rosiglitazone + sitagliptin + sulfonylurea + vildagliptin, random = list(~comp &#124; study, <br>    ~comp &#124; design), rho = 1/2, phi = 1/2, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~0 + acarbose + benfluorex + metformin + miglitol + pioglitazone + <br>    rosiglitazone + sitagliptin + sulfonylurea + vildagliptin, random = list(~cs(comp &#124; study), <br>    ~cs(comp &#124; design)), data = dat, measure = "MD", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_mods = weak_location_prior, prior_heterogeneity = fixed_correlation_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -0.8414, -0.7369, -1.128, -0.9499, -1.129, -1.234, -0.5700, -0.4175, -0.7000 | -0.8524, -0.7224, -1.128, -0.9472, -1.148, -1.256, -0.5643, -0.4321, -0.7267 |
| SE / posterior comparison | 0.2384, 0.2776, 0.1494, 0.2253, 0.2119, 0.1235, 0.3414, 0.2326, 0.3408 | max standardized difference: 0.1763 |
| variance/correlation | tau2=0.09991; rho=0.5000; gamma2=0.00000000004172; phi=0.5000 | sd(shared &#124; study)=0.3389; rho(study)=0.5000; sd(shared &#124; design)=0.1589; rho(design)=0.5000 |
| status | `equivalent` | close |

### `metadat/man/dat.tannersmith2016.Rd`

#### `help-090`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, random = ~esid &#124; studyid, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, random = ~cs(esid &#124; studyid), data = dat, measure = "GEN", prior_unit_information_sd = prior_scale, <br>    prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1023 | -- |
| SE / posterior comparison | 0.02411 | max standardized difference: -- |
| variance/correlation | tau2=0.01113; rho=0.3424 | -- |
| status | `equivalent` | error Bounded execution: brma.mv exceeded 25 minutes and 16 GB working memory for the 113-level esid &#124; studyid CS model. |

#### `help-091`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~aget1 + aget2, random = ~1 &#124; studyid/esid, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~aget1 + aget2, random = ~1 &#124; studyid/esid, data = dat, measure = "GEN", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_mods = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.3960, 0.04512, -0.05621 | 0.4029, 0.04384, -0.05561 |
| SE / posterior comparison | 0.1754, 0.02148, 0.01929 | max standardized difference: 0.05973 |
| variance/correlation | sigma2=0.0007633, 0.007322 | sd_total(random_total)=0.09964; sd(intercept &#124; esid:studyid)=0.08504; sd(intercept &#124; studyid)=0.04756 |
| status | `equivalent` | close |

#### `help-092`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~aget1 * aget2, random = ~1 &#124; studyid/esid, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~aget1 * aget2, random = ~1 &#124; studyid/esid, data = dat, measure = "GEN", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_mods = weak_location_prior, <br>    prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -2.390, 0.2656, 0.1027, -0.01251 | -2.440, 0.2681, 0.1078, -0.01281 |
| SE / posterior comparison | 1.545, 0.1233, 0.08963, 0.006887 | max standardized difference: 0.05696 |
| variance/correlation | sigma2=0.0008025, 0.007265 | sd_total(random_total)=0.09869; sd(intercept &#124; esid:studyid)=0.08473; sd(intercept &#124; studyid)=0.04637 |
| status | `equivalent` | close |

#### `help-093`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, V, mods = ~aget1 * aget2 + sexmix, random = ~1 &#124; studyid/esid, data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = V, mods = ~aget1 * aget2 + sexmix, random = ~1 &#124; studyid/esid, data = dat, <br>    measure = "GEN", prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, <br>    prior_mods = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | -2.535, 0.2827, 0.1078, 0.03125, -0.007145, -0.01332 | -2.127, 0.2489, 0.08622, 0.02793, -0.007308, -0.01151 |
| SE / posterior comparison | 1.603, 0.1278, 0.09323, 0.05492, 0.05731, 0.007165 | max standardized difference: 0.2648 |
| variance/correlation | sigma2=0.001056, 0.007264 | sd_total(random_total)=0.1001; sd(intercept &#124; esid:studyid)=0.08479; sd(intercept &#124; studyid)=0.04888 |
| status | `equivalent` | broadly close |

### `metadat/man/dat.white2020.Rd`

#### `help-094`

| `metafor::rma.mv()` | `RoBMA::brma.mv()` |
|---|---|
| <pre><code>metafor::rma.mv(yi, vi, random = list(~1 &#124; study_id, ~1 &#124; obs), data = dat)</code></pre> | <pre><code>brma.mv(yi = yi, V = vi, random = list(~1 &#124; study_id, ~1 &#124; obs), data = dat, measure = "ZCOR", <br>    prior_unit_information_sd = prior_scale, prior_effect = weak_location_prior, prior_heterogeneity = unit_sd_prior)</code></pre> |

| comparison | `metafor` | `brma.mv()` |
|---|---|---|
| coefficients | 0.1573 | 0.1597 |
| SE / posterior comparison | 0.03294 | max standardized difference: 0.07530 |
| variance/correlation | sigma2=0.01531, 0.06493 | sd_total(random_total)=0.2897; sd(intercept &#124; study_id)=0.1344; sd(intercept &#124; obs)=0.2533 |
| status | `equivalent` | close |

## Stable `metadat` calls requiring optional dependencies

These four stable calls require `ape` or `rms` during example setup. They are included 
in the automated corpus; the table makes their additional setup explicit.

| example | `metafor::rma.mv()` | `RoBMA::brma.mv()` | reference answer |
|---|---|---|---|
| `dat.lim2014` | `rma.mv(yi, vi, random=list(~1&#124;article, ~1&#124;esid, ~1&#124;species, ~1&#124;species.phy), R=list(species.phy=A))` | same call with `V=vi`, `measure="ZCOR"`, audit priors | `b=-0.1564`, `SE=0.1277` |
| `dat.moura2021` base | `rma.mv(yi, vi, random=list(~1&#124;study.id, ~1&#124;effect.size.id, ~1&#124;species.id, ~1&#124;species.id.phy), R=...)` | same call with `V=vi`, `measure="ZCOR"`, audit priors | `b=0.3682`, `SE=0.1300` |
| `dat.moura2021` moderators | previous call plus `mods=~spatially.pooled*temporally.pooled` | same `mods`, known `R`, audit priors | `b=(0.3457,0.0810,0.0599,-0.0729)` |
| `dat.obrien2003` spline | `rma.mv(yi, vi, mods=~rcs(bmicent,4), random=~1&#124;study/grp)` | same with `V=vi`; use the identical spline design matrix | `b=(4.2726,0.3142,-0.0560,1.1102)` |

## Official website examples

Repeated optimizer or plotting calls are grouped by webpage; every one of the 73 
literal occurrences is enumerated in the linked website inventory. A separate 
execution corpus fitted 67 occurrences (65 unique fit objects); six explicit 
failures comprise three unavailable optional optimizers and three bounded large 
timing/plotting fits. See the linked web-corpus status for each call and reason.

| webpage example | representative `rma.mv()` | representative `brma.mv()` | parity |
|---|---|---|---|
| Berkey et al. | `random=~outcome&#124;trial, struct="UN"` | `random=~us(outcome&#124;trial)` | mapped; help corpus |
| Crede et al. | `random=~1&#124;studyid/sampleid` or `~sampleid&#124;studyid` | nested formula or `~cs(sampleid&#124;studyid)` | mapped |
| Gleser-Olkin | known `V`, no random term | known `V`, `zero_sd_prior` | mapped fixed GLS |
| Konstantopoulos | nested and CS formulations | nested and `cs()` formulations | mapped; equal likelihood models retained |
| van Houwelingen | `struct="UN"`, group-specific slopes | `us()` with identical design | mapped; `cvvc` is not a Bayesian fit option |
| Viechtbauer profile examples | ML/REML random intercept | Bayesian random intercept | same likelihood, different inference |
| independent subgroups | `DIAG` versus `ID` | `diag()` versus `id()` | mapped; LRT replaced by posterior/model comparison |
| convergence: `rma.mv()` | 17 optimizer variants of one nested model | one nested model; MCMC controls differ | one Bayesian model |
| convergence: `rma()` | two optimizer variants of one random intercept | one random-intercept call | mapped |
| subgroup heterogeneity | `struct="DIAG"` | `diag()` | mapped; `cvvc` omitted |
| aggregated forest | Assink known `V`, nested random intercept | same `V` and nesting | mapped |
| multilevel `I2` | nested, UN, and fixed models | nested, `us()`, and zero-SD models | fits map; scalar `I2` reporting differs |
| `rma.uni()` versus `rma.mv()` | fixed/random meta-regressions | zero-SD/random-intercept `brma.mv()` | mapped |
| speeding fitting | sparse plus known phylogenetic `R` | dense backend plus known `R` | model maps; sparse input is normalized |
| two-stage longitudinal | response formula plus `struct="UN"` | explicit `yi`/moderator design plus `us()` | mapped |
| weights | ordinary nested model; page derives weights | nested model maps | fit maps; frequentist weights API differs |
| grouped forest plot | Assink model repeated | same mapped model | plotting context only |
| orchard plot | two `DIAG` inner/outer blocks | two named `diag()` blocks | mapped; prediction selection syntax differs |

## Complete source inventories

- [`metadat` 1.6-0: all 29 topics / 69 calls](tools/rma-mv-parity/metadat-inventory.md)
- [Official `metafor` site: all 73 calls plus help/test audit](tools/rma-mv-parity/metafor-web-inventory.md)
- [Official website execution status: 67/73 fitted](tools/rma-mv-parity/web-reference-status.md)
- [Machine-readable upstream help-call inventory; includes two development-only `metadat` calls](tools/rma-mv-parity/inventory.csv)

## Conclusions

The directly mapped family includes fixed GLS, nested/crossed random intercepts, 
CS/HCS/UN/ID/DIAG/AR/HAR/CAR/GEN/GDIAG random structures, dense known `V`, 
moderators, subsets, and known `R` random-intercept models. Exact public-interface 
gaps are spatial covariance, arbitrary `W`, and selected fixed coefficient-specific 
variance/correlation patterns. Fixed nested `sigma2` values are reproduced by 
expanding nested group IDs and spiking the corresponding block SDs. ML/REML tests, 
Knapp-Hartung degrees of freedom, LRTs, and optimizer controls are frequentist 
inference machinery rather than alternate `brma.mv()` likelihoods.
