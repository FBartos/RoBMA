# Selection Model

Declare the contextual sources retained by a selection model, how
estimate weights combine, and a publication-group column for
best-p-value selection when needed.

## Usage

``` r
selection_model(
  estimate_random_effects = "integrate",
  other_random_effects = "condition",
  known_sampling_variance = "integrate",
  weight_rule = "product",
  group = NULL
)
```

## Arguments

- estimate_random_effects:

  Either `"condition"` or `"integrate"` (default), for the declared
  estimate-level true-effect source. One-to-one grouping determines this
  role in multivariate models; it does not imply independent covariance.

- other_random_effects:

  Either `"condition"` (default) or `"integrate"`, for all remaining
  declared random-effect terms.

- known_sampling_variance:

  Either `"condition"` or `"integrate"` (default), for the complete
  sampling-error vector \\e \sim N(0,V)\\. Conditioning retains its
  realized value during selection normalization; integration uses the
  full sampling covariance. This choice also applies to univariate and
  independent sampling errors supplied through `vi` or `sei`. It does
  not change `V`.

- weight_rule:

  Either `"product"` (default), which multiplies estimate weights
  without publication grouping, or `"best"`, which uses the weight at
  the smallest p-value within each publication group.

- group:

  Optional unresolved data-column name, supplied as a bare name, a
  backtick-quoted name, or one character string. Active only for
  `weight_rule = "best"`, where the default `NULL` requests automatic
  group resolution when data are bound. Product models ignore this
  argument and require no group column.

## Value

A `selection_model` object containing fixed mode choices and a deferred
column name or `NULL`. The group is used only by the best rule.

## Details

This is RoBMA's re-export of
[`BayesTools::selection_model()`](https://fbartos.github.io/BayesTools/reference/selection_model.html).
The specification is stored in `prior_weightfunction(model = ...)`.
Constructor `selection = selection_model(...)` settings apply only to
generated default priors; explicit priors retain their own
specification. A group name is captured without evaluating data; RoBMA
binds it only for `weight_rule = "best"`. When `group = NULL`, best-rule
groups come from a supported constructor's `cluster` input or, for
independent models without a random structure, one group per estimate.
Other best-rule models require an explicit group column. Conditioned
sources remain unknown and estimated. Integrating estimate variation is
the default; conditioning on it defines a different selection model. The
specialized clustered interface retains its total heterogeneity,
allocation and within-/between-cluster I2 summaries.

With positive sampling standard errors, the sampling law can be written
as \\e = Sz\\, with \\z \sim N(0,R)\\ and \\V = SRS\\, where \\S\\
contains sampling standard errors and \\R\\ is their correlation matrix.
Conditioning retains all of \\e\\; no residual sampling error is
integrated separately. Covariance factorizations affect computation, not
the selection model. Selection thresholds continue to use the original
sampling standard errors.

When all outcome-generating sources are conditioned upon and weights are
positive almost surely, selection weights cancel and the observed law is
the ordinary Gaussian model. The likelihood then contains no information
about the selection weights. A retained context with zero acceptance
probability cannot complete the selection process. Integrating all
sources instead gives population-level selection among completed
studies. These choices are independent of prediction's
`conditioning_depth`. Non-unit observation `weights` require
`known_sampling_variance = "integrate"`; unit weights are equivalent to
omitting `weights`. Other observation-weight restrictions still apply.

`weight_rule = "best"` uses the weight of the smallest p-value's bin,
rather than the largest weight. `weight_rule = "product"` multiplies all
estimate weights. Its joint normalizer retains the dependencies in the
integrated covariance, regardless of grouping. In ensembles, only
best-rule branches require matching publication partitions; product-only
ensembles accept `group = NULL`.

See
[`bselmodel.mv()`](https://fbartos.github.io/RoBMA/dev/reference/bselmodel.mv.md)
for the supported conditioning cells, source and group requirements, and
covariance boundaries.

## See also

[`prior_weightfunction`](https://fbartos.github.io/RoBMA/dev/reference/prior_weightfunction.md),
[`bselmodel.mv()`](https://fbartos.github.io/RoBMA/dev/reference/bselmodel.mv.md),
[`RoBMA.mv()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.mv.md)

## Examples

``` r
prior_weightfunction(
  steps = .05,
  model = selection_model(known_sampling_variance = "condition")
)
#> omega[one-sided: .05] ~ CumDirichlet(1, 1)
#> Selection model:
#>   Estimate random effects: integrate (average effects before normalization).
#>   Other random effects: condition (retain unknown effects during normalization).
#>   Known sampling error: condition (retain full unknown error vector during normalization).
#>   Weight rule: product (Product of estimate weights).
#>   Group: unused for the product rule.
#>   Sources are resolved when model data are bound.
prior_weightfunction(
  steps = .05,
  model = selection_model(weight_rule = "best", group = paper_id)
)
#> omega[one-sided: .05] ~ CumDirichlet(1, 1)
#> Selection model:
#>   Estimate random effects: integrate (average effects before normalization).
#>   Other random effects: condition (retain unknown effects during normalization).
#>   Known sampling error: integrate (average full error vector before normalization).
#>   Weight rule: best (Weight of the best p-value).
#>   Group: 'paper_id' (unresolved data column).
#>   Sources are resolved when model data are bound.
```
