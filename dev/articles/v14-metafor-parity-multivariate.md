# Multivariate and Multilevel Meta-Analysis

This vignette is a starting point for fitting multivariate and
multilevel meta-analytic models with the `RoBMA` R package ([Bartoš &
Maier, 2020](#ref-RoBMA)). We walk through four datasets from the
`metadat` package ([Viechtbauer et al., 2026](#ref-metadat)), showing
the
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
call alongside its
[`metafor::rma.mv()`](https://wviechtb.github.io/metafor/reference/rma.mv.html)
counterpart ([Viechtbauer, 2010](#ref-metafor)). The examples move from
a simple three-level model to known sampling covariance matrices,
temporal correlation structures, and two simultaneous random structures.

The first two examples rely on default prior distributions and the last
two examples rely on unit-information-scaled prior distributions; see
the [*Prior
Distributions*](https://fbartos.github.io/RoBMA/dev/articles/v01-prior-distributions.md)
vignette for more details. For generic summaries, plots, and diagnostics
that also apply to
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
objects, see the [*Bayesian
Meta-Analysis*](https://fbartos.github.io/RoBMA/dev/articles/v02-bayesian-meta-analysis.md)
vignette.

## Three-Level Model

The `dat.konstantopoulos2011` dataset from the `metadat` package
contains 56 standardized mean differences from studies of modified
school calendars and student achievement. Positive values indicate
higher achievement for students on modified calendars, and schools are
clustered within districts. Effect sizes `yi` and sampling variances
`vi` are already provided.

``` r

data("dat.konstantopoulos2011", package = "metadat")
dat_konstantopoulos <- dat.konstantopoulos2011
head(dat_konstantopoulos)
#> 
#>   district school study year     yi    vi 
#> 1       11      1     1 1976 -0.180 0.118 
#> 2       11      2     2 1976 -0.220 0.118 
#> 3       11      3     3 1976  0.230 0.144 
#> 4       11      4     4 1976 -0.300 0.144 
#> 5       12      1     5 1989  0.130 0.014 
#> 6       12      2     6 1989 -0.260 0.014
```

The `metafor` package fit uses a compound-symmetric random effect for
schools within districts. We set `struct = "CS"` explicitly so that the
fitted random structure is visible in the call. The shorter
`brma(..., cluster = district)` syntax for this special case is shown in
the [*Multilevel
Meta-Analysis*](https://fbartos.github.io/RoBMA/dev/articles/v10-metafor-parity-multilevel.md)
vignette.

``` r

fit1_metafor <- metafor::rma.mv(
  yi, vi,
  random = ~ school | district,
  struct = "CS",
  data = dat_konstantopoulos
)
fit1_metafor
#> 
#> Multivariate Meta-Analysis Model (k = 56; method: REML)
#> 
#> Variance Components:
#> 
#> outer factor: district (nlvls = 11)
#> inner factor: school   (nlvls = 11)
#> 
#>             estim    sqrt  fixed 
#> tau^2      0.0978  0.3127     no 
#> rho        0.6653             no 
#> 
#> Test for Heterogeneity:
#> Q(df = 55) = 578.8640, p-val < .0001
#> 
#> Model Results:
#> 
#> estimate      se    zval    pval   ci.lb   ci.ub    
#>   0.1847  0.0846  2.1845  0.0289  0.0190  0.3504  * 
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The matching
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
call passes the sampling variances through `V` and uses the `cs()`
random-effect formula tag. The `cs(school | district)` term estimates
one shared heterogeneity standard deviation and one within-district
correlation for the school effects.

``` r

fit1_brma <- brma.mv(
  yi = yi, V = vi, measure = "SMD",
  random = ~ cs(school | district),
  data = dat_konstantopoulos, seed = 1,
  chains = 4, adapt = 2000, burnin = 5000, sample = 10000,
  parallel = TRUE, autofit = FALSE
)
```

[`summary()`](https://rdrr.io/r/base/summary.html) reports the pooled
effect and the random-effect parameters monitored by the fitted Bayesian
model. Here `tau` is the shared marginal heterogeneity SD corresponding
to the `metafor` package variance component, and `district: rho` is the
compound-symmetric correlation.

``` r

summary(fit1_brma, include_mcmc_diagnostics = FALSE)
#> 
#> Bayesian Multivariate Random-Effects Model (k = 56)
#> 
#> Estimates
#>      Mean    SD 0.025   0.5 0.975
#> mu  0.181 0.088 0.013 0.179 0.363
#> tau 0.321 0.054 0.237 0.314 0.449
#> rho 0.631 0.137 0.335 0.645 0.856
```

## Known Sampling Covariance

The `dat.assink2016` dataset from the `metadat` package contains
standardized mean differences from 17 studies on recidivism and mental
health among delinquent juveniles. Most studies contributed multiple
estimates for different delinquent behaviors, so the example needs both
a working sampling covariance matrix and random effect specification.
The working sampling covariance matrix and the random-effects formula
address different sources of dependence. The working sampling covariance
matrix encodes sampling-error dependence that is treated as known, while
the random-effects formula estimates residual heterogeneity across the
levels named in `random`.

The
[`metafor::vcalc()`](https://wviechtb.github.io/metafor/reference/vcalc.html)
call constructs a working sampling covariance matrix from the
effect-size variances, study labels, delinquency type, and effect-size
identifiers.

``` r

data("dat.assink2016", package = "metadat")
dat_assink <- dat.assink2016
V_assink   <- metafor::vcalc(
  vi, cluster = study, type = deltype, obs = esid,
  rho = c(0.7, 0.5), data = dat_assink
)
head(dat_assink)
#> 
#>   study esid id      yi     vi pubstatus year deltype 
#> 1     1    1  1  0.9066 0.0740         1  4.5 general 
#> 2     1    2  2  0.4295 0.0398         1  4.5 general 
#> 3     1    3  3  0.2679 0.0481         1  4.5 general 
#> 4     1    4  4  0.2078 0.0239         1  4.5 general 
#> 5     1    5  5  0.0526 0.0331         1  4.5 general 
#> 6     1    6  6 -0.0507 0.0886         1  4.5 general
```

The
[`metafor::rma.mv`](https://wviechtb.github.io/metafor/reference/rma.mv.html)
function then specifies the `V` matrix and nested random intercept for
studies and effect sizes within studies as the input.

``` r

fit2_metafor <- metafor::rma.mv(
  yi, V_assink,
  random = ~ 1 | study / esid,
  data = dat_assink
)
fit2_metafor
#> 
#> Multivariate Meta-Analysis Model (k = 100; method: REML)
#> 
#> Variance Components:
#> 
#>             estim    sqrt  nlvls  fixed      factor 
#> sigma^2.1  0.0704  0.2654     17     no       study 
#> sigma^2.2  0.1508  0.3883    100     no  study/esid 
#> 
#> Test for Heterogeneity:
#> Q(df = 99) = 840.9174, p-val < .0001
#> 
#> Model Results:
#> 
#> estimate      se    zval    pval   ci.lb   ci.ub      
#>   0.3618  0.0933  3.8794  0.0001  0.1790  0.5446  *** 
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
call mirrors the same inputs: `V` is the known sampling covariance
matrix and `random` describes the multilevel specification.

``` r

fit2_brma <- brma.mv(
  yi = yi, V = V_assink, measure = "SMD",
  random = ~ 1 | study / esid,
  data = dat_assink, seed = 2,
  chains = 4, adapt = 2000, burnin = 5000, sample = 10000
)
```

The `summary` function here deviates from the metafor input because it
summarizes parameters with directly specified prior distributions. Here,
it is the total between-study heterogeneity tau and variance allocation
of tau^2 into the individual levels.

``` r

summary(fit2_brma, include_mcmc_diagnostics = FALSE)
#> 
#> Bayesian Multivariate Random-Effects Model (k = 100)
#> 
#> Estimates
#>                        Mean    SD 0.025   0.5 0.975
#> mu                    0.355 0.097 0.169 0.353 0.551
#> tau_total             0.478 0.062 0.375 0.471 0.617
#> tau2_prop(esid_study) 0.676 0.138 0.398 0.682 0.926
#> tau2_prop(study)      0.324 0.138 0.074 0.318 0.602
```

A more detailed summary of the between-study heterogeneity can be
obtained via the `summary_heterogeneity` function. The first
heterogeneity table summarizes the shared `study/esid` heterogeneity
allocation, including the total heterogeneity and variance fractions;
the following tables summarize the block-specific standard deviations.

``` r

summary_heterogeneity(fit2_brma)
#> 
#> Heterogeneity Estimates (study/esid):
#>                        Mean Median 0.025 0.975
#> tau_total             0.478  0.471 0.375 0.617
#> tau2_total            0.232  0.222 0.141 0.381
#> tau2_prop(esid_study) 0.676  0.682 0.398 0.926
#> tau2_prop(study)      0.324  0.318 0.074 0.602
#> 
#> Heterogeneity Estimates (esid_study):
#>       Mean Median 0.025 0.975
#> tau  0.387  0.385 0.323 0.462
#> tau2 0.151  0.148 0.104 0.214
#> 
#> Heterogeneity Estimates (study):
#>       Mean Median 0.025 0.975
#> tau  0.270  0.264 0.111 0.467
#> tau2 0.081  0.070 0.012 0.218
```

## Heterogeneous AR(1) Structure

The `dat.ishak2007` dataset from the `metadat` package contains results
from 46 studies of deep-brain stimulation for Parkinson’s disease.
Negative mean differences indicate improvements in UPDRS motor-function
scores, and repeated follow-up measurements within studies motivate a
temporal correlation model. The data are first reshaped to long format,
and
[`metafor::vcalc()`](https://wviechtb.github.io/metafor/reference/vcalc.html)
uses a working autocorrelation of `0.97` for the sampling errors.

``` r

data("dat.ishak2007", package = "metadat")
dat_ishak <- dat.ishak2007
dat_ishak <- reshape(
  dat_ishak,
  direction = "long",
  idvar     = "study",
  v.names   = c("yi", "vi"),
  varying   = list(c(2, 4, 6, 8), c(3, 5, 7, 9))
)
dat_ishak <- dat_ishak[order(dat_ishak$study, dat_ishak$time), ]
dat_ishak <- dat_ishak[!is.na(dat_ishak$yi), ]
rownames(dat_ishak) <- NULL
dat_ishak$time_factor <- factor(dat_ishak$time)
V_ishak <- metafor::vcalc(
  vi, cluster = study, time1 = time,
  phi = 0.97, data = dat_ishak
)
head(dat_ishak)
#> 
#>               study mdur mbase time    yi   vi time_factor 
#> 1    Alegret (2001) 16.1  53.6    1 -33.4 14.3           1 
#> 2 Barichella (2003) 13.5  45.3    1 -20.0  7.3           1 
#> 3 Barichella (2003) 13.5  45.3    3 -30.0  5.7           3 
#> 4     Berney (2002) 13.6  45.6    1 -21.1  7.3           1 
#> 5   Burchiel (1999) 13.6  48.0    1 -20.0  8.0           1 
#> 6   Burchiel (1999) 13.6  48.0    2 -20.0  8.0           2
```

The `metafor` package uses `struct = "HAR"` for a heterogeneous AR(1)
true-effect structure over time within studies. The fixed effects
estimate one mean difference at each time point, with the intercept
suppressed.

``` r

fit3_metafor <- metafor::rma.mv(
  yi, V_ishak,
  mods = ~ 0 + time_factor,
  random = ~ time | study,
  struct = "HAR",
  data = dat_ishak
)
fit3_metafor
#> 
#> Multivariate Meta-Analysis Model (k = 82; method: REML)
#> 
#> Variance Components:
#> 
#> outer factor: study (nlvls = 46)
#> inner factor: time  (nlvls = 4)
#> 
#>              estim    sqrt  k.lvl  fixed  level 
#> tau^2.1    22.7258  4.7672     24     no      1 
#> tau^2.2    33.7294  5.8077     22     no      2 
#> tau^2.3    26.1425  5.1130     25     no      3 
#> tau^2.4    31.1804  5.5839     11     no      4 
#> rho         0.8832                    no        
#> 
#> Test for Residual Heterogeneity:
#> QE(df = 78) = 856.1635, p-val < .0001
#> 
#> Test of Moderators (coefficients 1:4):
#> QM(df = 4) = 873.1569, p-val < .0001
#> 
#> Model Results:
#> 
#>               estimate      se      zval    pval     ci.lb     ci.ub      
#> time_factor1  -25.9047  1.0118  -25.6021  <.0001  -27.8879  -23.9216  *** 
#> time_factor2  -27.4608  1.1406  -24.0768  <.0001  -29.6963  -25.2254  *** 
#> time_factor3  -28.6559  1.0321  -27.7634  <.0001  -30.6789  -26.6330  *** 
#> time_factor4  -26.4934  1.3817  -19.1743  <.0001  -29.2015  -23.7853  *** 
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
call uses the `har()` random-effects formula tag. To obtain separate
parameters for each time_factor level, we specify the independent factor
contrast via `prior_mods` argument (only suppressing the intercept via
-1 would fix the default’s level coefficient 0). Because this outcome is
a generic continuous outcome and sample sizes are not reported in the
source table, we use one of the primary-study sample sizes to estimate
the unit-information scale. The Alegret (2001) row reports a 3-month
sampling variance of `14.3`, and the primary study reports `15`
patients. We compute the corresponding unit-information standard
deviation with
[`estimate_unit_information_sd()`](https://fbartos.github.io/RoBMA/dev/reference/estimate_unit_information_sd.md)
and use it directly in the prior distribution specification.

``` r

ishak_unit_information_sd <- estimate_unit_information_sd(sei = sqrt(14.3), ni = 15)

fit3_brma <- brma.mv(
  yi = yi, V = V_ishak, measure = "GEN",
  mods = ~ 0 + time_factor,
  random = ~ har(time | study),
  prior_effect              = NULL,
  prior_mods                = prior_factor(
    distribution = "normal",
    parameters   = list(mean = 0, sd = ishak_unit_information_sd),
    contrast     = "independent"
  ),
  prior_unit_information_sd = ishak_unit_information_sd,
  data = dat_ishak, seed = 3,
  chains = 4, adapt = 2000, burnin = 5000, sample = 10000,
  parallel = TRUE, autofit = FALSE
)
```

``` r

summary(fit3_brma, include_mcmc_diagnostics = FALSE)
#> 
#> Bayesian Multivariate Mixed-Effect Model (k = 82)
#> 
#> Common Estimates
#>                     Mean    SD 0.025   0.5 0.975
#> tau_common         5.043 0.785 3.651 4.993 6.745
#> tau2_mult(time[1]) 0.834 0.245 0.428 0.811 1.382
#> tau2_mult(time[2]) 1.172 0.253 0.719 1.160 1.724
#> tau2_mult(time[3]) 0.949 0.213 0.600 0.926 1.427
#> tau2_mult(time[4]) 1.046 0.335 0.442 1.025 1.749
#> rho                0.827 0.073 0.647 0.840 0.929
#> 
#> Meta-Regression
#>                   Mean    SD   0.025     0.5   0.975
#> time_factor[1] -25.413 1.019 -27.399 -25.423 -23.362
#> time_factor[2] -26.924 1.123 -29.038 -26.941 -24.774
#> time_factor[3] -28.128 1.023 -30.113 -28.146 -26.130
#> time_factor[4] -25.782 1.388 -28.530 -25.795 -23.071
```

The `HAR` structure lets the heterogeneity differ by time point while
estimating one autocorrelation across neighboring time points. This
separates the temporal structure of the true effects from the working
temporal dependence already encoded in `V_ishak`.

## Study and Treatment Random Effects

The `dat.begg1989` dataset from the `metadat` package contains 20
estimates of two-year disease-free survival rates for acute
nonlymphocytic leukemia. The studies compare allogeneic bone-marrow
transplantation (BMT) and conventional chemotherapy (CMO), with both
controlled two-arm studies and uncontrolled one-arm studies. The example
uses the reported survival rates and their sampling variances directly
on the raw rate scale. Because the source table reports survival rates
and standard errors but not sample sizes, we recover effective sample
sizes and use them to scale the prior distributions.

``` r

data("dat.begg1989", package = "metadat")
dat_begg <- dat.begg1989
dat_begg$trt <- relevel(factor(dat_begg$trt), ref = "CMO")
dat_begg$ni  <- as.numeric(round(with(dat_begg, yi * (1 - yi) / sei^2)))
head(dat_begg)
#> 
#>   study trt arms     yi   sei     vi ni 
#> 1     1 BMT    2 0.4600 0.118 0.0139 18 
#> 2     1 CMO    2 0.2500 0.074 0.0055 34 
#> 3     2 BMT    2 0.5000 0.100 0.0100 25 
#> 4     2 CMO    2 0.2300 0.067 0.0045 39 
#> 5     3 BMT    2 0.4700 0.129 0.0166 15 
#> 6     3 CMO    2 0.4200 0.086 0.0074 33
```

The `metafor` package model combines a study intercept with a
treatment-specific random effect. The study intercept captures
heterogeneity in baseline survival. The treatment-specific component
allows the BMT and CMO estimates to vary within studies, with the
treatment-level correlation fixed at zero by `rho = 0`.

``` r

fit4_metafor <- metafor::rma.mv(
  yi, vi,
  mods = ~ trt,
  random = list(
    ~ 1 | study,
    ~ trt | study
  ),
  struct = "CS",
  rho    = 0,
  data   = dat_begg,
  method = "ML"
)
fit4_metafor
#> 
#> Multivariate Meta-Analysis Model (k = 20; method: ML)
#> 
#> Variance Components:
#> 
#>             estim    sqrt  nlvls  fixed  factor 
#> sigma^2    0.0000  0.0000     16     no   study 
#> 
#> outer factor: study (nlvls = 16)
#> inner factor: trt   (nlvls = 2)
#> 
#>             estim    sqrt  fixed 
#> tau^2      0.0013  0.0357     no 
#> rho        0.0000            yes 
#> 
#> Test for Residual Heterogeneity:
#> QE(df = 18) = 27.6534, p-val = 0.0675
#> 
#> Test of Moderators (coefficient 2):
#> QM(df = 1) = 15.5962, p-val < .0001
#> 
#> Model Results:
#> 
#>          estimate      se     zval    pval   ci.lb   ci.ub      
#> intrcpt    0.3225  0.0196  16.4880  <.0001  0.2841  0.3608  *** 
#> trtBMT     0.1463  0.0370   3.9492  <.0001  0.0737  0.2189  *** 
#> 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

The
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md)
sets names for the two random structures explicitly. The
[`prior_random()`](https://fbartos.github.io/RoBMA/dev/reference/prior_random.md)
call specifies the study and treatment blocks directly because it fixes
the treatment-level correlation at zero.

``` r

begg_unit_information_sd <- estimate_unit_information_sd(sei = dat_begg$sei, ni = dat_begg$ni)

sd_prior <- prior(
  distribution = "normal",
  parameters = list(mean = 0, sd = begg_unit_information_sd),
  truncation = list(lower = 0, upper = Inf)
)
cor_zero <- prior(
  distribution = "spike",
  parameters = list(location = 0)
)
cs_zero <- random_covariance(
  cor       = cor_zero,
  cor_scale = "cor"
)

fit4_brma <- brma.mv(
  yi = yi, V = vi, ni = ni, measure = "GEN",
  mods = ~ trt,
  random = list(
    study     = ~ 1 | study,
    treatment = ~ cs(trt | study)
  ),
  prior_heterogeneity = prior_random(
    study     = random_block(sd = sd_prior),
    treatment = random_block(
      sd         = sd_prior,
      covariance = cs_zero
    )
  ),
  prior_unit_information_sd = begg_unit_information_sd,
  data = dat_begg, seed = 4,
  chains = 4, adapt = 2000, burnin = 5000, sample = 10000
)
```

``` r

summary(fit4_brma, include_mcmc_diagnostics = FALSE)
#> 
#> Bayesian Multivariate Mixed-Effect Model (k = 20)
#> 
#> Common Estimates
#>                 Mean    SD 0.025   0.5 0.975
#> study: tau     0.040 0.028 0.002 0.036 0.105
#> treatment: tau 0.040 0.027 0.002 0.037 0.101
#> treatment: rho 0.000 0.000 0.000 0.000 0.000
#> 
#> Meta-Regression
#>            Mean    SD 0.025   0.5 0.975
#> intercept 0.327 0.026 0.276 0.326 0.380
#> trt[BMT]  0.136 0.043 0.052 0.136 0.220
```

The study component captures baseline heterogeneity across studies,
while the treatment component captures additional treatment-specific
variation within studies.

``` r

summary_heterogeneity(fit4_brma)
#> 
#> Heterogeneity Estimates (study):
#>       Mean Median 0.025 0.975
#> tau  0.040  0.036 0.002 0.105
#> tau2 0.002  0.001 0.000 0.011
#> 
#> Heterogeneity Estimates (treatment):
#>                 Mean Median 0.025 0.975
#> tau            0.040  0.037 0.002 0.101
#> tau2           0.002  0.001 0.000 0.010
#> treatment: rho 0.000  0.000 0.000 0.000
```

## References

Bartoš, F., & Maier, M. (2020). *RoBMA: An R package for robust Bayesian
meta-analyses*. <https://CRAN.R-project.org/package=RoBMA>

Viechtbauer, W. (2010). Conducting meta-analyses in R with the metafor
package. *Journal of Statistical Software*, *36*(3), 1–48.
<https://doi.org/10.18637/jss.v036.i03>

Viechtbauer, W., White, T., Noble, D., Senior, A., & Hamilton, W. K.
(2026). *Metadat: Meta-analysis datasets*.
<https://github.com/wviechtb/metadat>
