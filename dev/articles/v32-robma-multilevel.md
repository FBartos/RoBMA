# Multilevel Robust Bayesian Meta-Analysis

**This vignette accompanies the manuscript [Robust Bayesian Multilevel
Meta-Analysis: Adjusting for Publication Bias in the Presence of
Dependent Effect Sizes](https://doi.org/10.31234/osf.io/9tgp2_v1)
preprinted at *PsyArXiv* ([Bartoš et al.,
2026](#ref-bartos2025robust)).**

This vignette reproduces the first example from the manuscript. For
multilevel meta-regression with moderators, see the [*Multilevel Robust
Bayesian Model-Averaged
Meta-Regression*](https://fbartos.github.io/RoBMA/dev/articles/v33-robma-multilevel-metaregression.md)
vignette.

Multilevel Robust Bayesian Meta-Analysis (RoBMA) extends the standard
RoBMA framework to datasets containing multiple effect sizes from the
same studies. We demonstrate the method using data from Johnides et al.
([2025](#ref-johnides2025secondary)), who meta-analyzed 412 effect sizes
from 128 studies investigating secondary benefits of family-based
treatments for childhood disorders.

## When to Use Multilevel RoBMA

Multilevel RoBMA is appropriate for datasets that contain multiple
effect sizes from the same studies. The multilevel approach explicitly
models the clustering structure through the `cluster` argument, which:

- accounts for dependencies among effect sizes from the same study,
- partitions heterogeneity into within-study and between-study
  components,
- simultaneously adjusts for publication bias at the study level
  (corresponding to selective reporting).

## Loading the Data

We load the package and examine the dataset structure.

``` r

library(RoBMA)
data("Johnides2025", package = "RoBMA")
```

The dataset contains effect sizes (`d`), standard errors (`se`), and
study identifiers (`study`).

``` r

head(Johnides2025)
#>                 study       d        se
#> 1 Price et al. (2012)  0.0000 0.1337909
#> 2 Price et al. (2012) -0.1250 0.1337909
#> 3 Price et al. (2012)  0.0000 0.1337909
#> 4 Price et al. (2012)  0.2112 0.1337909
#> 5 Price et al. (2012)  0.0140 0.1337909
#> 6 Price et al. (2012) -0.0186 0.1337909
```

Multiple rows share the same study name, indicating that these effect
sizes come from the same study. For example, the first five effect sizes
are from “Price et al. (2012)” and might be more similar to each other
than to effect sizes from other studies.

## Prior Distributions

To reproduce the RoBMA 3.6.1 version results, which used a different
prior distribution formulation and scaling settings, we need to specify
the prior distributions manually (see the [*Prior
Distributions*](https://fbartos.github.io/RoBMA/dev/articles/v01-prior-distributions.md)
vignette for more details on specifying prior distributions).

``` r

prior_effect_361        <- prior("normal",   list(mean = 0, sd = 1))
prior_heterogeneity_361 <- prior("invgamma", list(shape = 1, scale = 0.15))
```

## Fitting the Multilevel Model

We specify the effect sizes as `yi` and standard errors as `sei`
argument. Setting `measure = "SMD"` specifies that we are using
standardized mean differences as the effect-size input, which is
important for the prior distribution specification (here the effect size
and heterogeneity prior distributions are specified explicitly since
they differ from the default specification in the 3.6.1 version; the
publication-bias specification remains the same). The `cluster` argument
specifies which effect sizes belong together.

``` r

fit <- RoBMA(
  yi = d, sei = se, measure = "SMD", cluster = study,
  prior_effect        = prior_effect_361,
  prior_heterogeneity = prior_heterogeneity_361,
  adapt = 5000, burnin = 5000, sample = 20000,
  parallel  = TRUE, seed = 1, data = Johnides2025
)
```

*Note: This model takes a few minutes with parallel processing enabled.
Without parallel processing, the chains run one after another and the
fit takes correspondingly longer.*

## Interpreting the Results

The [`summary()`](https://rdrr.io/r/base/summary.html) output contains
two main sections.

``` r

summary(fit)
#> 
#> Robust Bayesian Model-Averaged Multilevel Random-Effects Model (k = 412, clusters = 128)
#> 
#> Component Inclusion
#>                  Prior prob. Post. prob. Inclusion BF error%(Inclusion BF)
#> Effect                 0.500       0.623        1.656                8.286
#> Heterogeneity          0.500       1.000   >59999.000                   NA
#> Publication Bias       0.500       1.000   >59999.000                   NA
#> 
#> Estimates
#>      Mean    SD 0.025   0.5 0.975 error(MCMC) error(MCMC)/SD  ESS R-hat
#> mu  0.074 0.066 0.000 0.082 0.190     0.00241          0.036  759 1.001
#> tau 0.391 0.034 0.327 0.390 0.460     0.00063          0.019 2981 1.001
#> rho 0.492 0.080 0.328 0.494 0.640     0.00107          0.013 5519 1.000
#> 
#> Publication Bias
#>                    Mean    SD 0.025   0.5 0.975 error(MCMC) error(MCMC)/SD   ESS R-hat
#> omega[0,0.025]    1.000 0.000 1.000 1.000 1.000          NA             NA    NA    NA
#> omega[0.025,0.05] 0.996 0.016 0.946 1.000 1.000     0.00006          0.004 61524 1.000
#> omega[0.05,0.5]   0.942 0.050 0.815 0.955 0.998     0.00033          0.007 24268 1.000
#> omega[0.5,0.95]   0.213 0.058 0.124 0.205 0.345     0.00142          0.024  1683 1.000
#> omega[0.95,0.975] 0.213 0.058 0.124 0.205 0.345     0.00142          0.024  1683 1.000
#> omega[0.975,1]    0.213 0.058 0.124 0.205 0.345     0.00142          0.024  1683 1.000
#> PET               0.000 0.000 0.000 0.000 0.000     0.00000             NA     0    NA
#> PEESE             0.000 0.000 0.000 0.000 0.000     0.00000             NA     0    NA
#> P-value intervals for publication bias weights omega correspond to one-sided p-values.
```

### Components Summary

This table shows Bayes factors testing the presence of effect,
heterogeneity, and publication bias. Each row displays the prior
probability, posterior probability after observing the data, and the
inclusion Bayes factor comparing models with versus without that
component.

In our example, the inclusion BF for the effect is 1.66, indicating
undecided evidence for an effect. The inclusion BFs for heterogeneity
and publication bias are reported as larger than 59,999, indicating
extreme evidence for both components.

### Model-Averaged Estimates

This section provides meta-analytic estimates averaged across all
models, weighted by their posterior probabilities. The average effect
size is `mu`, the overall heterogeneity is `tau`, and the heterogeneity
allocation parameter is `rho`. The heterogeneity allocation parameter
`rho` indicates the proportion of heterogeneity allocated to the
cluster-level component: `rho = 0` means all heterogeneity is within
clusters, `rho = 1` means all heterogeneity is between clusters, and
`rho = 0.5` means equal split. The remaining parameters summarize the
weight function and PET/PEESE regression for publication bias.

In our example, we find a small effect size (*d* = 0.074, 95% CI
\[0.000, 0.190\]), substantial heterogeneity (τ = 0.391, 95% CI \[0.327,
0.460\]), and nearly balanced heterogeneity allocation (ρ = 0.492, 95%
CI \[0.328, 0.640\]). The large heterogeneity combined with the small
pooled effect suggests substantial variation in true effects. This
distribution of true effect sizes can be summarized with the
[`summary_heterogeneity()`](https://fbartos.github.io/RoBMA/dev/reference/summary_heterogeneity.md)
function. A prediction interval for a new true effect can be obtained
with
[`pooled_effect()`](https://fbartos.github.io/RoBMA/dev/reference/pooled_effect.md):

``` r

summary_heterogeneity(fit)
#> 
#> Heterogeneity Estimates:
#>                  Mean Median  0.025  0.975
#> tau             0.391  0.390  0.327  0.460
#> rho             0.492  0.494  0.328  0.640
#> tau [within]    0.277  0.277  0.222  0.335
#> tau [between]   0.273  0.272  0.207  0.348
#> tau2            0.154  0.152  0.107  0.211
#> tau2 [within]   0.078  0.077  0.049  0.112
#> tau2 [between]  0.076  0.074  0.043  0.121
#> I2             85.010 85.193 80.134 88.889
#> I2 [within]    43.176 42.937 30.751 56.902
#> I2 [between]   41.834 41.944 27.558 55.294
#> H2              6.821  6.754  5.034  9.000
pooled_effect(fit)
#> 
#> Pooled Effect Size
#>     Mean Median 0.025 0.975 PI 0.025 PI 0.975
#> mu 0.074  0.082 0.000 0.190   -0.717    0.847
```

The wide prediction interval includes both negative and positive
effects: some studies may show benefits while others show harm.

## Model Types Summary

To understand which publication bias mechanisms the data support, we
examine the posterior distribution across model types:

``` r

summary_models(fit)
#> 
#> Effect
#>              Hypothesis  Prior prob. Post. prob. Inclusion BF error%(Inclusion BF)
#> Spike(0)     Null              0.500       0.377        0.604                8.286
#> Normal(0, 1) Alternative       0.500       0.623        1.656                8.286
#> 
#> Heterogeneity
#>                   Hypothesis  Prior prob. Post. prob. Inclusion BF error%(Inclusion BF)
#> Spike(0)          Null              0.500       0.000        0.000                   NA
#> InvGamma(1, 0.15) Alternative       0.500       1.000          Inf                   NA
#> 
#> Publication Bias
#>                                 Hypothesis  Prior prob. Post. prob. Inclusion BF error%(Inclusion BF)
#> None                            Null              0.500       0.000        0.000                   NA
#> omega[two-sided: .05]           Alternative       0.042       0.000        0.000                   NA
#> omega[two-sided: .05, .1]       Alternative       0.042       0.000        0.000                   NA
#> omega[one-sided: .05]           Alternative       0.042       0.000        0.000                   NA
#> omega[one-sided: .025, .05]     Alternative       0.042       0.000        0.000                   NA
#> omega[one-sided: .05, .5]       Alternative       0.042       0.896      197.694                1.252
#> omega[one-sided: .025, .05, .5] Alternative       0.042       0.104        2.676                1.252
#> PET                             Alternative       0.125       0.000        0.000                   NA
#> PEESE                           Alternative       0.125       0.000        0.000                   NA
```

Essentially all of the posterior probability is allocated to selection
models (weight functions). The publication-bias adjustment therefore
reflects selective reporting rather than small-study effects (PET/PEESE
regression); specifically, one-sided selection operating on marginally
significant p-values and on the direction of the effect.

## Visualizing the Weight Function

The weight function shows how publication probability varies across
p-value ranges:

``` r

plot_weightfunction(fit)
```

![](v32-robma-multilevel_files/figure-html/unnamed-chunk-8-1.png)

The plot displays one-sided p-value cutoffs (x-axis) against relative
publication probability (y-axis), averaged across weight function
specifications. A flat line at 1.0 indicates no publication bias; values
below 1.0 indicate suppression.

In our example, we notice a small drop in the relative publication
probability at the cutoff corresponding to p = 0.05 (one-sided, i.e.,
0.10 two-sided: selection on marginally significant p-values) and a much
sharper drop corresponding to the direction of the effect (p = 0.50).

## Comparison with Single-Level RoBMA

To understand the importance of accounting for nested structure, we
compare our results with a standard single-level RoBMA that ignores
dependencies among effect sizes from the same study.

``` r

fit_simple <- RoBMA(
  yi = d, sei = se, measure = "SMD", 
  prior_effect        = prior_effect_361,
  prior_heterogeneity = prior_heterogeneity_361,
  adapt = 5000, burnin = 5000, sample = 20000,
  parallel  = TRUE, seed = 1, data = Johnides2025
)
```

``` r

summary(fit_simple)
#> 
#> Robust Bayesian Model-Averaged Random-Effects Model (k = 412)
#> 
#> Component Inclusion
#>                  Prior prob. Post. prob. Inclusion BF error%(Inclusion BF)
#> Effect                 0.500       0.054        0.057                3.414
#> Heterogeneity          0.500       1.000   >59999.000                   NA
#> Publication Bias       0.500       1.000   >59999.000                   NA
#> 
#> Estimates
#>      Mean    SD 0.025   0.5 0.975 error(MCMC) error(MCMC)/SD   ESS R-hat
#> mu  0.001 0.012 0.000 0.000 0.033     0.00012          0.010 10611 1.001
#> tau 0.371 0.022 0.329 0.371 0.414     0.00038          0.017 21305 1.006
#> 
#> Publication Bias
#>                    Mean    SD 0.025   0.5 0.975 error(MCMC) error(MCMC)/SD   ESS R-hat
#> omega[0,0.025]    1.000 0.000 1.000 1.000 1.000          NA             NA    NA    NA
#> omega[0.025,0.05] 0.997 0.013 0.958 1.000 1.000     0.00007          0.005 38172 1.001
#> omega[0.05,0.5]   0.950 0.043 0.840 0.962 0.999     0.00022          0.005 37700 1.000
#> omega[0.5,0.95]   0.225 0.058 0.168 0.220 0.293     0.00273          0.047 19830 1.192
#> omega[0.95,0.975] 0.225 0.058 0.168 0.220 0.293     0.00273          0.047 19830 1.192
#> omega[0.975,1]    0.225 0.058 0.168 0.220 0.293     0.00273          0.047 19830 1.192
#> PET               0.004 0.071 0.000 0.000 0.000     0.00367          0.052   123 1.298
#> PEESE             0.000 0.011 0.000 0.000 0.000     0.00005          0.004 20000 1.291
#> P-value intervals for publication bias weights omega correspond to one-sided p-values.
```

The single-level model differs by:

- omitting `cluster`: it treats all effect sizes as independent,
- estimating only one heterogeneity parameter: it cannot distinguish
  within-study from between-study variation,
- potentially biased inference: ignoring dependencies can lead to
  overconfident conclusions.

For the Johnides et al. data, the single-level RoBMA finds strong
evidence for the absence of an effect, while the multilevel model
provides only weak evidence against it. Properly accounting for data
structure leads to more conservative and appropriate inferences.

## References

Bartoš, F., Maier, M., & Wagenmakers, E.-J. (2026). Robust Bayesian
multilevel meta-analysis: Adjusting for publication bias in the presence
of dependent effect sizes. *Behavior Research Methods*.
<https://doi.org/10.31234/osf.io/9tgp2_v2>

Johnides, B. D., Borduin, C. M., Sheerin, K. M., & Kuppens, S. (2025).
Secondary benefits of family member participation in treatments for
childhood disorders: A multilevel meta-analytic review. *Psychological
Bulletin*, *151*(1), 1–32. <https://doi.org/10.1037/bul0000462>
