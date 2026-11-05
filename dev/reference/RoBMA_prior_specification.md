# Prior specification for model-averaging

The [`RoBMA`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md)
and
[`RoBMA.mv`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.mv.md)
functions extend the prior specification described in
[`prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/prior_specification.md)
by allowing mixture priors that define competing hypotheses for Bayesian
model-averaging. This enables multimodel inference via the product space
method (Carlin and Chib 1995; Lodewyckx et al. 2011) , where a single
MCMC chain explores a joint model space containing all competing models.

## Arguments

- prior_effect:

  prior distribution(s) for the alternative effect component(s).

- prior_heterogeneity:

  prior distribution(s) for the alternative heterogeneity component(s).

- prior_mods:

  prior distribution(s) for alternative moderator components. A single
  prior applies to all terms; a named list can specify term-specific
  components. This argument can be specified only when `mods` is
  specified. For a formula without an intercept (`~ 0 + x`), an
  `intercept` entry other than a point prior at zero stops with an error
  of class `RoBMA_prior_intercept_conflict` (also of class
  `RoBMA_prior_error`). A named alternative `intercept` must be a single
  simple or point prior; it is the complete intercept prior and replaces
  the global effect mixture. Specify it in one place: null-list
  intercepts and explicit global effect declarations left without a
  consumer are refused. Ordinary term entries replace only their
  alternative component.

- prior_scale:

  prior distribution(s) for alternative scale-regression components. A
  single prior applies to all terms; a named list can specify
  term-specific components. This argument can be specified only when
  `scale` is specified. A named alternative `intercept` must be a single
  positive simple or point prior; it is the complete SD prior used by
  `log(intercept)`. Global heterogeneity declarations remain meaningful
  when they supply another random root or its inclusion weights;
  declarations replaced without a consumer are refused. Null-list
  intercept declarations are unavailable.

- prior_heterogeneity_allocation:

  prior distribution(s) for the alternative cluster-level heterogeneity
  allocation component(s). This argument can be specified only when
  `cluster` is specified.

- prior_bias:

  prior distribution(s) for alternative publication-bias component(s),
  such as weight functions, PET, or PEESE.

  Alternative prior arguments can be:

  - A single prior distribution object (creates a mixture with one
    alternative)

  - A list of prior distributions (creates a mixture with multiple
    alternatives)

  - `NULL` or `FALSE` (omits the alternative hypothesis component)

  See
  [`publication_bias_prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/publication_bias_prior_specification.md)
  for details on specifying publication-bias priors and
  [`prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/prior_specification.md)
  for details on specifying meta-analytic parameter priors.

- prior_effect_null:

  prior distribution(s) for the null effect component(s).

- prior_heterogeneity_null:

  prior distribution(s) for the null heterogeneity component(s).

- prior_mods_null:

  prior distribution(s) for null moderator components. A named
  `intercept` entry, including `NULL` or `FALSE`, is unavailable: put
  the complete intercept prior in `prior_mods`, or use global prior
  arguments without a formula-list intercept. The same rule applies to
  `prior_scale_null`. A single prior applies to all terms; a named list
  can specify term-specific components. This argument can be specified
  only when `mods` is specified. The null and the alternative prior of a
  factor term need the same contrast; a conflict stops with an error of
  class `RoBMA_prior_contrast_conflict` (also of class
  `RoBMA_prior_error`) that names the term and the remedies. The same
  rule applies to `prior_scale_null`.

- prior_scale_null:

  prior distribution(s) for null scale-regression components. A single
  prior applies to all terms; a named list can specify term-specific
  components. This argument can be specified only when `scale` is
  specified.

- prior_heterogeneity_allocation_null:

  prior distribution(s) for the null cluster-level heterogeneity
  allocation component(s). This argument can be specified only when
  `cluster` is specified.

- prior_bias_null:

  prior distribution(s) for null publication-bias component(s), usually
  [`prior_none()`](https://fbartos.github.io/RoBMA/dev/reference/prior_none.md).
  See
  [`publication_bias_prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/publication_bias_prior_specification.md).

  Null prior arguments can be:

  - A single prior distribution object (creates a mixture with one null)

  - A list of prior distributions (creates a mixture with multiple
    nulls)

  - `NULL` or `FALSE` (omits the null hypothesis component)

  Defaults to a point mass (spike) at zero for effect and heterogeneity
  parameters.

- model_type:

  character string specifying predefined publication-bias model
  ensembles. One of:

  - `"PSMA"` (default): Full RoBMA-PSMA ensemble with 6 weight
    functions + PET + PEESE

  - `"6w"`: Six weight function models

  - `"2w"`: Two weight function models

  - `"PP"`: PET-PEESE models only

  Custom `prior_bias` replaces the preset alternative bias components.
  If `prior_bias` is omitted, `model_type` determines the default
  alternatives even when `prior_bias_null` is customized or omitted.

## Details

### The product space method for model-averaging

RoBMA implements Bayesian model-averaging using the product space method
(Carlin and Chib 1995; Lodewyckx et al. 2011) . Rather than fitting each
model separately and combining results after the fitting is complete,
this approach embeds all competing models within a single joint model
space. A discrete model indicator variable selects which model component
is "active" at each MCMC iteration, and the posterior probability of
each model is estimated from the proportion of iterations spent in that
model's subspace.

This is achieved by specifying **mixture priors** that combine:

- **Null hypothesis component(s)**: Typically point masses (spikes)
  representing the absence of an effect, heterogeneity, or publication
  bias

- **Alternative hypothesis component(s)**: Continuous distributions
  representing the presence of the parameter

The mixture structure enables direct computation of:

- **Inclusion Bayes factors**: Evidence for the presence vs. absence of
  each parameter (e.g., effect, heterogeneity, bias)

- **Model-averaged estimates**: Posterior estimates that account for
  model uncertainty by weighting across all models

- **Conditional estimates**: Posterior estimates given that a particular
  hypothesis is true

### Default mixture prior structure

By default,
[`RoBMA`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md)
creates the following mixture priors:

|  |  |  |
|----|----|----|
| **Parameter** | **Null component** | **Alternative component** |
| Effect (\\\mu\\) | Spike(0) | Normal(0, UISD-scaled) |
| Heterogeneity (\\\tau\\) | Spike(0) | Normal+(0, UISD-scaled) |
| Publication bias | No bias (prior_none) | Weight functions + PET + PEESE |
| Allocation (\\\rho\\) | (none by default) | Beta(1, 1) |

See
[`prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/prior_specification.md)
for details on how the UISD scaling is determined.

### Specifying custom mixture priors

#### Single prior vs. list of priors

Each `prior_*` and `prior_*_null` argument accepts either a single prior
or a list:

- **Single prior**: Creates one component in the mixture

- **List of priors**: Creates multiple components, enabling
  model-averaging across different prior specifications (e.g., different
  effect size scales)

When multiple priors are specified in a list, they are all treated as
alternative (or null) hypothesis components for the main inclusion Bayes
factor. The inclusion Bayes factor tests the combined evidence for all
alternative components against all null components. The detailed summary
output shows posterior probabilities and inclusion Bayes factors for
each individual component separately, allowing finer-grained inference
about which specific prior specification is most supported by the data.

#### Omitting hypothesis components

For top-level mixture components, setting a prior argument to `NULL` or
`FALSE` removes that hypothesis component:

- `prior_effect_null = NULL`: No null hypothesis for effect (assumes
  effect exists)

- `prior_effect = NULL`: No alternative hypothesis (assumes no effect)

- `prior_bias_null = NULL` or `FALSE`: No "no bias" model (assumes some
  bias mechanism)

- `prior_bias = NULL` or `FALSE`: No bias model (assumes no bias
  mechanism)

For moderator and scale formula priors, an omitted or whole-argument
`NULL` uses the corresponding default. Whole-argument `FALSE` is not a
supported shortcut for omitting every term. Omit a particular
alternative or null component with a named `NULL` or `FALSE` entry, for
example `prior_mods_null = list(x1 = NULL)`. Ordinary named terms
replace only that alternative or null component. A named alternative
`intercept` instead supplies the complete intercept prior, which must be
a single simple or point prior (for `prior_scale`, a single positive
prior). Null-list intercept declarations are refused, as are explicit
global priors discarded without a consumer. A global heterogeneity
declaration can still supply an unscaled root or component-inclusion
weights when scale intercepts have their own priors; retain it when it
has such a consumer.

An ordered term's mixture belongs on its scalar total. Supply the
complete `prior_ordered(total = ...)` in one term component and omit the
separate term component with a named `NULL` entry. For example:


    total <- BayesTools::prior_mixture(list(prior("point", list(0)),
                                             prior("normal", list(0, 1))),
                                        is_null = c(TRUE, FALSE))
    BMA(yi = yi, sei = sei, mods = ~ f, data = dat, measure = "SMD",
        prior_mods = list(f = prior_ordered(total = total)),
        prior_mods_null = list(f = NULL), only_priors = TRUE)

Combining an ordered container with a separate term null or alternative
is unavailable and raises `RoBMA_prior_error` with the corresponding
omission.

#### Parameter-specific customization for moderators

For moderator priors (`prior_mods`, `prior_mods_null`), you can specify
different priors for different predictors using named lists. A single
prior object applies to all moderator terms:


    RoBMA(...,
      mods = ~ x1 + x2,
      prior_mods_null = list(x1 = NULL)  # No null for x1, default null for x2
    )

This allows testing whether specific moderators have effects while
assuming others are always included or excluded.

### Mixture priors for different parameter types

#### Effect size (\\\mu\\) and heterogeneity (\\\tau\\)

Effect and heterogeneity priors combine spike-and-slab components:


    # Default: spike(0) null + normal alternative
    # Custom: multiple alternatives with different scales
    RoBMA(...,
      prior_effect = list(
        prior("normal", list(mean = 0, sd = 0.5)),
        prior("normal", list(mean = 0, sd = 1.0))
      )
    )

#### Publication bias

Bias priors combine different publication bias adjustment mechanisms:

- No publication bias (null hypothesis)

- Weight functions: Selection models with different step functions

- PET/PEESE: Regression-based bias adjustment

The `model_type` argument provides convenient presets:


    RoBMA(..., model_type = "PSMA")  # Full ensemble (default)
    RoBMA(..., model_type = "PP")    # Only PET-PEESE models

See
[`publication_bias_prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/publication_bias_prior_specification.md)
for the bias-prior constructors and preset model spaces.

#### Heterogeneity allocation (\\\rho\\) for multilevel models

When `cluster` is specified, \\\rho\\ controls the cluster-level
variance allocation. The decomposition is \\\tau\_{between} =
\tau\sqrt{\rho}\\ and \\\tau\_{within} = \tau\sqrt{1 - \rho}\\. By
default, only an alternative (Beta(1,1)) is specified, but null
hypotheses can be added:


    RoBMA(...,
      cluster = study,
      prior_heterogeneity_allocation_null = prior("spike", list(location = 0.5))
    )

#### Moderator and scale regression terms

When `mods` is specified, the effect prior becomes the intercept prior,
and each moderator receives a mixture prior. The **intercept** inherits
the effect mixture structure, while **regression coefficients** get
separate spike-and-slab mixtures.

Similarly, when `scale = ~ ...` is specified for location-scale models,
the heterogeneity prior becomes the scale intercept prior, and each
scale predictor receives a mixture prior following the same pattern as
moderator terms. Note that the scale intercept (i.e., baseline
heterogeneity) is always positive and the multiplicative scale
coefficients require non-zero prior distribution on the intercept.

### Model weights and prior odds

Each component in a mixture prior has an associated prior weight (prior
model probability). User-specified components receive equal weights
unless their prior objects set `prior_weights`. Publication-bias presets
use fixed weights: `"2w"` and `"6w"` split the alternative bias mass
equally across their weight functions, `"PP"` splits it equally across
PET and PEESE, and `"PSMA"` assigns half of the alternative bias mass to
the six weight functions combined and one quarter each to PET and PEESE.
Custom weights can be specified via the `prior_weights` argument in
individual prior objects.

The prior odds between null and alternative hypotheses affect the Bayes
factor interpretation. With equal prior weights on null and alternative,
the posterior odds equal the Bayes factor.

## References

Carlin BP, Chib S (1995). “Bayesian model choice via Markov chain Monte
Carlo methods.” *Journal of the Royal Statistical Society: Series B
(Methodological)*, **57**(3), 473–484.
[doi:10.1111/j.2517-6161.1995.tb02042.x](https://doi.org/10.1111/j.2517-6161.1995.tb02042.x)
.  
  
Lodewyckx T, Kim W, Lee MD, Tuerlinckx F, Kuppens P, Wagenmakers E
(2011). “A tutorial on Bayes factor estimation with the product space
method.” *Journal of Mathematical Psychology*, **55**(5), 331–347.
[doi:10.1016/j.jmp.2011.06.001](https://doi.org/10.1016/j.jmp.2011.06.001)
.

## See also

[`prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/prior_specification.md)
for base prior specification options,
[`publication_bias_prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/publication_bias_prior_specification.md)
for publication-bias priors,
[`RoBMA`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md) and
[`RoBMA.mv`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.mv.md)
for the main model-averaging functions,
[`prior`](https://fbartos.github.io/BayesTools/reference/prior.html) for
creating prior distribution objects

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("metadat", quietly = TRUE)) {
  data(dat.lehmann2018, package = "metadat")

  # Default mixture priors (spike-and-slab for effect and heterogeneity)
  fit1 <- RoBMA(
    yi      = yi,
    vi      = vi,
    data    = dat.lehmann2018,
    measure = "SMD",
    seed    = 1,
    silent  = TRUE
  )

  # Custom alternative priors (two effect size scales)
  fit2 <- RoBMA(
    yi           = yi,
    vi           = vi,
    data         = dat.lehmann2018,
    measure      = "SMD",
    prior_effect = list(
      prior("normal", list(mean = 0, sd = 0.35)),
      prior("normal", list(mean = 0, sd = 0.70))
    ),
    seed         = 1,
    silent       = TRUE
  )

  # No null hypothesis for effect (assume effect exists)
  fit3 <- RoBMA(
    yi                = yi,
    vi                = vi,
    data              = dat.lehmann2018,
    measure           = "SMD",
    prior_effect_null = NULL,
    seed              = 1,
    silent            = TRUE
  )

  # Only selection model bias adjustment (no PET-PEESE)
  fit4 <- RoBMA(
    yi         = yi,
    vi         = vi,
    data       = dat.lehmann2018,
    measure    = "SMD",
    model_type = "6w",
    seed       = 1,
    silent     = TRUE
  )

  # Meta-regression with different null specifications per moderator
  fit5 <- RoBMA(
    yi              = yi,
    vi              = vi,
    mods            = ~ Preregistered,
    data            = dat.lehmann2018,
    measure         = "SMD",
    prior_mods_null = list(Preregistered = NULL),
    seed            = 1,
    silent          = TRUE
  )
}
} # }
```
