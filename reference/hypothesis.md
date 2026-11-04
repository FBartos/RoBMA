# Hypothesis Bayes Factors

Computes Bayes factors for scalar hypotheses written as expressions, for
example `"mu = 0"`, `"year > 0"`, or
`"mu_alloc[alternate] > mu_alloc[random]"`.

## Usage

``` r
hypothesis(object, ...)

# Default S3 method
hypothesis(object, hypothesis, ...)

# S3 method for class 'brma'
hypothesis(
  object,
  hypothesis,
  component = "auto",
  standardized_coefficients = FALSE,
  conditional = FALSE,
  logBF = FALSE,
  BF01 = FALSE,
  seed = NULL,
  density_method = "qCMDE",
  density_control = NULL,
  n_samples = 10000,
  columns = "default",
  ...
)

# S3 method for class 'marginal_means.brma'
hypothesis(
  object,
  hypothesis,
  parameter = NULL,
  logBF = FALSE,
  BF01 = FALSE,
  seed = NULL,
  density_method = NULL,
  density_control = NULL,
  columns = "default",
  ...
)

bf_hypothesis(object, ...)

BF_hypothesis(object, ...)
```

## Arguments

- object:

  a fitted `brma`, `BMA`, `RoBMA`, or a posterior object accepted by
  [`BayesTools::hypothesis_BF()`](https://fbartos.github.io/BayesTools/reference/hypothesis_BF.html).

- ...:

  additional arguments passed to methods. For fitted `brma` objects,
  unused arguments trigger a warning.

- hypothesis:

  character vector with scalar hypothesis statements. For
  `marginal_means.brma` objects, constants are specified on the fitted
  linear-predictor scale, even when the object stores a display
  transformation. Separate statements may reference different model
  parameters, also the marginal means of different parameters of
  `marginal_means.brma` objects (unless `parameter` is set); each
  statement itself must resolve to one scalar model target.

- component:

  parameter component. Defaults to `"auto"`, which infers the component
  when possible. Use `"mods"` (alias `"location"`), `"scale"`, or
  `"random"` to disambiguate terms used in multiple model components. A
  statement whose references still name several model parameters stops
  with an error of class `BayesTools_hypothesis_ambiguous`, with the
  parent class `BayesTools_hypothesis_statement`: the statement or
  `component` needs to change, the test itself is not unavailable. An
  ambiguous reference that BayesTools refuses with the class
  `BayesTools_parameter_ambiguous` (for example the contrast coefficient
  `"g{1}"` of a factor term in the location and the scale formula) has
  these classes too, before its BayesTools classes. A statement whose
  parameters belong to another component than `component`, including a
  name that is known only outside `component` next to references of
  `component` (with or without factor levels, e.g.
  `"log_tau_g[a] > mu_g[a]"` with `component = "mods"`), stops with the
  classes `BayesTools_hypothesis_statement` and
  `BayesTools_hypothesis_component_mismatch`
  ([`plot.brma()`](https://fbartos.github.io/RoBMA/reference/plot.brma.md)
  stops with `RoBMA_component_mismatch` for its own component
  mismatches). The exception is the alias of a factor term of both the
  location and the scale formula (for example `g` for
  `mods = ~ g, scale = ~ g`) next to a name that is unknown within
  `component`, such as a name of the other component or a display label
  (for example `"g > log_tau_intercept"` with `component = "mods"`): the
  statement is resolved without `component`, where the alias is
  ambiguous, and stops with `BayesTools_hypothesis_statement` followed
  by the classes of the BayesTools refusal of its first unresolved
  reference, the ambiguity of the alias
  (`BayesTools_parameter_ambiguous`, without
  `BayesTools_hypothesis_ambiguous`, as `component` is set already) or
  the unknown name (`BayesTools_parameter_not_found`, for example in
  `"log_tau_intercept < g"` or, with a level of the alias, in
  `"g[a] > log_tau_intercept"`). The display labels of the summary
  tables (for example `"(mu) intercept"`, `"exp(intercept)"`, or
  `"(mu) g[a]"`) name the parameters they label, with and without
  `component`. The random component supports point, interval, and
  directional hypotheses for semantic standard deviation, variance,
  correlation, and allocation quantities. Point-null hypotheses require
  a direct parameter reference, a factor level, or a linear combination
  of the levels of one factor term (for example `"g[a] = g[b]"` or
  `"2 * g[a] = 0.1"`). Point hypotheses on a random-effect variance are
  evaluated through its standard deviation, so that both give the same
  Bayes factor at positive variance values. At zero, the square map is
  singular: the variance is tested on its own scale and follows its own
  prior-ordinate eligibility rules. A statement comparing such a point
  with a region of the variance (for example
  `"tau2 = 0.09 vs tau2 > 0.09"`) evaluates the point through the
  standard deviation and the region on the variance draws. Certified
  `exp(affine)` fitted-scale hypotheses are available with KDE only for
  atom-free, unconditional scalar targets. Values at a prior point mass
  (for example, gated components at 0) have no Savage-Dickey Bayes
  factor; the Component Inclusion table compares the exclusion and
  inclusion of gated components. Point nulls at an exact support
  boundary (for example, a variance proportion at 0) use the one-sided
  prior ordinate when BayesTools classifies it as exact, finite, and
  positive. Publication-bias parameters are not supported (class
  `BayesTools_hypothesis_target`), including the quantities of the
  publication-bias prior such as the weight-function coordinates
  `"omega[0,0.025]"` and `"bias_indicator"`. Other quantities of the
  fitted model that are no tested model parameter, such as inclusion
  indicators (`"mu_x_indicator"`), the component inclusions of variance
  allocations (`"inclusion(study)"`), or latent cluster effects, are
  refused with the class `BayesTools_hypothesis_target` as well;
  [`hypothesis_quantities()`](https://fbartos.github.io/BayesTools/reference/hypothesis_quantities.html)
  lists the tested quantities.

- standardized_coefficients:

  whether moderator and scale coefficients are tested on the
  standardized predictor scale. Defaults to `FALSE`.

- conditional:

  whether to use the conditional posterior for product-space
  model-averaging objects. Defaults to `FALSE`. When this argument is
  omitted and the selected parameter has both null and alternative
  components, a warning notes that the full ensemble is used. Pass
  `FALSE` explicitly to retain that test without the warning, or `TRUE`
  to test only models where the parameter is active. Point hypotheses on
  linear combinations of factor levels of model-averaged objects require
  `conditional = TRUE`: the levels share the null component of the
  averaged prior.

- logBF:

  whether to display the Bayes factor on the log scale.

- BF01:

  whether to display the inverse Bayes factor.

- seed:

  optional seed used by BayesTools for sampled prior quantities.

- density_method:

  posterior density method. Defaults to `"qCMDE"` for fitted `brma`
  objects. `"KDE"` uses the standard BayesTools kernel density estimate.
  `"normal"` uses the BayesTools normal approximation for point-null
  hypotheses on fitted `brma` objects. `"qCMDE"` and `"IWMDE"` attach
  RoBMA likelihood-aware posterior ordinates for point-null hypotheses
  and propagate their `BF_error` estimates printed as `error%(BF)`.
  `BF_error` is a Monte Carlo diagnostic for the computed posterior
  ordinate, not a complete standard error. For a subsampled mixed
  posterior it is a worst-correlation delta upper bound combining the
  selected active-row and complete indicator-chain components. It
  excludes uncertainty in the prior ordinate and, for IWMDE, uncertainty
  from estimating the conditional weight function. At a location that is
  also a posterior atom, qCMDE/IWMDE reports the continuous ordinate;
  the atom probability is a separate discrete mass. Use
  [`density_diagnostics()`](https://fbartos.github.io/RoBMA/reference/density_diagnostics.md)
  to inspect the attached computation and reliability diagnostics. For
  `marginal_means.brma` objects, `"normal"` is not supported. Their
  default `NULL` reuses the density method stored by
  [`marginal_means()`](https://fbartos.github.io/RoBMA/reference/marginal_means.md);
  explicitly request `"KDE"` to override it. `"qCMDE"` and `"IWMDE"`
  compute missing point-null ordinates from the stored source model.
  Matching is case-insensitive. qCMDE/IWMDE ordinates use the exact
  fitted-to-original structural coefficient map and induced prior
  density stored by BayesTools. qCMDE/IWMDE supports exact linear maps:
  coefficients, factor levels of every contrast, and linear combinations
  of levels. For nonlinear joint maps, such as an exponentiated
  intercept that also depends on varying slopes, KDE point and
  directional hypotheses are available only when structural metadata
  certifies an atom-free, unconditional scalar target and the point is
  inside the open transformed support. Alternatively, use
  `standardized_coefficients = TRUE`.

  The statements are planned and evaluated by the hypothesis engine of
  BayesTools
  ([`BayesTools::hypothesis_BF_model()`](https://fbartos.github.io/BayesTools/reference/hypothesis_BF_model.html));
  RoBMA describes its models to it
  ([`BayesTools::hypothesis_provider()`](https://fbartos.github.io/BayesTools/reference/hypothesis_provider.html))
  and provides the density methods qCMDE and IWMDE
  ([`BayesTools::hypothesis_density_backend()`](https://fbartos.github.io/BayesTools/reference/hypothesis_density_backend.html)).
  Every statement is planned before it is evaluated: the plan records
  the tested target (its weights on the fitted coefficients, its
  BayesTools prior density and declared atoms), the classification of
  each point value by
  [`BayesTools::prior_ordinate_status()`](https://fbartos.github.io/BayesTools/reference/prior_ordinate_status.html),
  and whether each density method can evaluate it;
  [`hypothesis_quantities()`](https://fbartos.github.io/BayesTools/reference/hypothesis_quantities.html)
  renders the same plans. A point hypothesis whose prior ordinate is not
  exact, finite, and regular stops with the BayesTools condition class
  (for example `BayesTools_point_mass_at_null`,
  `BayesTools_infinite_ordinate`, or `BayesTools_inexact_ordinate`).
  Other refusals of a test that cannot be computed have the class
  `BayesTools_hypothesis_unavailable` and one of
  `BayesTools_hypothesis_fixed` (the quantity is fixed by the fitted
  model), `BayesTools_hypothesis_target` (no supported target or prior
  density), or `BayesTools_hypothesis_method` (the density method is
  unavailable for the target); linear combinations that BayesTools
  cannot certify keep the class `BayesTools_linear_target_unavailable`.
  A method refusal because the fitted model does not support qCMDE/IWMDE
  also has the classes of the refusals of these requests by
  [`plot.brma()`](https://fbartos.github.io/RoBMA/reference/plot.brma.md)
  and
  [`marginal_means()`](https://fbartos.github.io/RoBMA/reference/marginal_means.md):
  `RoBMA_density_method_unavailable` and one class naming the cause,
  `RoBMA_density_method_random_unknown_v`
  ([`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
  random-formula models without known `V`),
  `RoBMA_density_method_scale_components` (component-specific scale
  formulas), or `RoBMA_density_method_glmm` (IWMDE for binomial and
  Poisson GLMMs). qCMDE/IWMDE refusals for a cause for which
  [`plot.brma()`](https://fbartos.github.io/RoBMA/reference/plot.brma.md)
  refuses these methods have its classes too:
  `RoBMA_density_method_conditional_random` (conditional random-effect
  statements), `RoBMA_density_method_random_target` (random-effect
  quantities without a supported scalar random-component coordinate),
  and `RoBMA_density_method_original_scale` (original-scale coefficients
  whose fitted map is nonlinear, such as `exp(affine)` targets). As in
  [`plot.brma()`](https://fbartos.github.io/RoBMA/reference/plot.brma.md),
  allocated component SDs in a predictor-scaled block require
  `standardized_coefficients = TRUE` or `density_method = 'KDE'`; the
  block-wide restriction also covers an unchanged centered intercept. As
  in
  [`plot.brma()`](https://fbartos.github.io/RoBMA/reference/plot.brma.md),
  the capability of the fitted model is checked before these causes and
  before other qCMDE/IWMDE refusals of the tested quantity: a request to
  which both apply stops with the capability refusal. A qCMDE/IWMDE
  ordinate that cannot be constructed stops with
  `RoBMA_density_ordinate_error`; when a fitted chain has no eligible
  row of the target, the error also has the class
  `RoBMA_density_chain_coverage` and its message names
  `density_method = "KDE"` and longer chains as remedies. The refusal
  for fewer than 20 continuous active samples of the target takes
  precedence over this diagnosis: such a target is refused for that
  reason alone, without the class `RoBMA_density_chain_coverage`, even
  when a fitted chain has none of its rows.

  Statements that need to be restated, or an argument that needs to
  change, stop with the class `BayesTools_hypothesis_statement` without
  `BayesTools_hypothesis_unavailable`: point expressions that are no
  direct reference (for example `"2 * mu = 0"`), nonlinear expressions
  of factor levels, several point values of a random-effect variance in
  one statement, one quantity named by several names in one statement
  (for example `"mu > 0 & intercept < 1"` or `"g[a] > mu_g[a]"`; use one
  name), model-averaged marginal-means statements that mix point and
  region events or span several marginal means (levels or scalar
  marginal means), ambiguous references
  (`BayesTools_hypothesis_ambiguous`, see `component`), references to
  parameters of another component
  (`BayesTools_hypothesis_component_mismatch`), marginal-means
  statements that refer to another parameter than `parameter`
  (`BayesTools_hypothesis_parameter_mismatch`), factor contrast
  coefficients such as `"g{1}"` (state the hypothesis on factor levels;
  a selector of a coefficient that is a level, such as `"g{1}"` of a
  treatment factor, also has the classes of the BayesTools refusal
  `BayesTools_selector_unavailable`, which names the level), references
  to unknown names or factor levels, and statements that reference no
  parameter (followed by the classes of the BayesTools refusal of such a
  statement, `BayesTools_hypothesis_no_parameters` and
  `BayesTools_parameter_resolution_error`, on fitted objects and on
  marginal means). An unknown name or level also has the classes
  `BayesTools_parameter_not_found` and
  `BayesTools_parameter_resolution_error`, with which BayesTools refuses
  it (on fitted objects and in linear combinations of marginal-means
  levels, the BayesTools condition itself with its fields `alias` and
  `available`). Every other refusal of a statement's references by
  BayesTools (class `BayesTools_parameter_resolution_error`) stops
  `hypothesis()` too, on fitted objects and on marginal means, as the
  BayesTools condition with its classes and fields after
  `BayesTools_hypothesis_statement`, also when BayesTools raises it only
  when it evaluates the statement (for example the whole factor term
  next to one of its levels, `"g[a] > mu_g"`, as an unknown quantity).
  Missing or unsupported fitted metadata of the tested target (a
  parameter catalog without RoBMA parameters, the coefficient transform,
  the fitted coordinates, or the linear weights of a factor level; a fit
  of an older RoBMA/BayesTools build) stop with the class
  `RoBMA_refit_required`: refit the model. Its parent
  `BayesTools_refit_required` is the class of every error of RoBMA and
  BayesTools that asks for a refit; BayesTools' own refit errors stop
  `hypothesis()` unconverted.

- density_control:

  named list of qCMDE/IWMDE tuning settings. Supported entries are
  `n_points` (default `100`), `samples` (the fixed posterior-row sample
  size, default `500` for qCMDE and `1000` for IWMDE; use `Inf` for the
  eligible-row census), `target_relative_mcse` (default `0.05`),
  `normalization_points` (default `NULL`, resolved to
  `max(50, n_points)`), `normalization_prob` (default `0.999`; lower
  than 1 for qCMDE), and `display_grid` (default `"adaptive"`).
  `integration_control` (default `NULL`) retains the fitted
  selection-integration settings. Supply a control created by
  [`set_selection_likelihood_control()`](https://fbartos.github.io/RoBMA/reference/set_selection_likelihood_control.md)
  to change those settings for this post-fit calculation, for example
  `list(integration_control = set_selection_likelihood_control(max_points_per_scramble = 32768))`.
  This entry is available for Gaussian selection models with a fitted
  integration plan. The fitted object and posterior draws are unchanged.
  The maximum point budget controls factor QMC fallback; analytic and
  deterministic quadrature rules remain unchanged. Point ordinates use
  one state-independent simple random sample selected before ordinate
  contributions are evaluated. Multiple direct scalar point ordinates
  for one target share the same conditional-normalization pass while
  retaining separate diagnostics. A finite estimate is returned
  independently of its sample diagnostics. If the fixed sample does not
  meet the precision target, a warning recommends increasing `samples`
  or using the census; if the census does not meet the target, obtain
  more posterior draws. `display_grid` is immaterial for point-only
  requests. qCMDE normalizes each posterior row's conditional density
  over a range covering at least its central `normalization_prob` mass,
  on `normalization_points` nodes and the nested grid of those nodes and
  their midpoints. The mass of a row left outside the range, its
  truncation \\t\\, is therefore at most `1 - normalization_prob`: zero
  for rows whose Gaussian likelihood kernel meets a normal prior and
  whose evaluated values are those of the resulting truncated normal
  (their normalizer is its closed-form mass, so no lattice decides it),
  exactly computed for such rows that keep the lattice normalizer,
  bounded for a Gaussian kernel with another proper scalar prior, and
  estimated from the tails (each tail extended until its estimate is at
  most `(1 - normalization_prob) / 2`) otherwise. The tail extension has
  a limit, ends at values the joint density cannot evaluate, and cannot
  help a tail that does not decay, so the estimate for such a row can
  exceed that target and the bound below with it. A row normalizer
  missing the fraction \\t\\ overstates that row's density by the factor
  \\1 / (1 - t)\\, so the ordinate is overstated by at most \\t / (1 -
  t)\\ for the largest row truncation, which is at most
  `(1 - normalization_prob) / normalization_prob`; at the default
  `0.999` that is about 0.1 percent, and it is attained when all rows
  share one conditional law.
  [`density_diagnostics()`](https://fbartos.github.io/RoBMA/reference/density_diagnostics.md)
  report the largest row truncation (`normalization_truncation`), how it
  is known (`normalization_truncation_status`), and this bound
  (`truncation_ordinate_bound`), together with the change between the
  two nested grids (`ordinate_relative_change`), which estimates the
  discretization error of the lattice normalizers, and the
  contribution-share weighted lattice-resolution estimate
  (`lattice_resolution_relative_error`). The resolution estimate uses
  the Gaussian kernel width at the nested-grid maximum, or local log
  curvature when no kernel is available. Neither the resolution
  component nor the combined error estimate is a rigorous bound; a sharp
  prior core can be narrower than its likelihood kernel. Exact
  nonlattice normalizers have no lattice-resolution component. Set
  `normalization_prob` closer to 1 to reduce the truncation, and
  increase `normalization_points` to reduce the discretization. IWMDE
  checks its normalization mass on `normalization_points` points over
  the central `normalization_prob` quantile range of the target draws,
  widened by 10 percent on each side. A point ordinate warns at relative
  MCSE at least 5 percent, ESS below 100, maximum contribution share at
  least 20 percent, or fewer than 100 finite contributions. These
  same-sample diagnostics do not suppress a finite fixed-design
  estimate. The qCMDE ordinate error estimate (nested-grid change,
  lattice resolution and the tail-truncation bound) warns above 2.5
  percent and is rejected above 5 percent; IWMDE normalization error
  warns above 5 percent and is rejected above 10 percent. For
  retained-location conditional normalizers, adaptive-quadrature
  sensitivity warns above 2.5 percent and is rejected above 5 percent.

- n_samples:

  number of prior samples/grid points used by BayesTools for
  deterministic marginal prior densities.

- columns:

  output columns passed to
  [`BayesTools::hypothesis_BF()`](https://fbartos.github.io/BayesTools/reference/hypothesis_BF.html).
  `"default"` returns the compact BayesTools table with `Alternative`,
  `Null`, `BF`, and `error%(BF)` columns. `"all"` also returns `prior`,
  `posterior`, and `method` columns. The `prior` and `posterior` columns
  are diagnostics rather than probabilities for all hypothesis types.

- parameter:

  optional marginal-means parameter/term used to disambiguate the
  hypothesis expression. Without it, each statement is tested on the
  marginal means of the parameter it names, so that separate statements
  may refer to different parameters (for example
  `c("g[a] > 0", "h[b] > 0")`); the rows are those of the statements
  tested one at a time, in the order of the statements. One statement
  must refer to the marginal means of one parameter: a statement that
  names several (for example `"g[a] > h[b]"`), or an alias that names
  several parameters, stops with an error of class
  `BayesTools_hypothesis_ambiguous`, with the parent class
  `BayesTools_hypothesis_statement`; `parameter` selects the parameter
  of an ambiguous alias. With `parameter`, every statement is tested on
  that parameter, and a statement that refers to the marginal means of
  another parameter stops with an error of class
  `BayesTools_hypothesis_parameter_mismatch`, with the parent class
  `BayesTools_hypothesis_statement`; for a statement that names several
  parameters (`"g[a] > h[b]"` with `parameter = "g"`) its message is the
  one given without `parameter`. A statement that references no
  marginal-means parameter stops with the class
  `BayesTools_hypothesis_statement`; for a statement without parameter
  symbols (for example `"1 > 0"`, also with `parameter`) also with
  `BayesTools_hypothesis_no_parameters` and
  `BayesTools_parameter_resolution_error`, and for an unknown name, or
  an unknown level of a factor term, also with
  `BayesTools_parameter_not_found` and
  `BayesTools_parameter_resolution_error`, as on fitted objects.

## Value

A `BayesTools_hypothesis_BF` BayesTools table.

## Details

For a continuous focal parameter \\\theta\\, nuisance parameters
\\\xi\\, and an unnormalized joint posterior kernel \\q(\theta,\xi)\\,
qCMDE estimates the posterior ordinate as
\$\$\widehat{p}\_{Q}(\theta^{\*}\mid y)=\frac{1}{S}\sum\_{i=1}^{S}
\frac{q(\theta^{\*},\xi\_{i})}{\int\_{\Theta} q(t,\xi\_{i})\\dt}.\$\$
Each row is therefore normalized over the focal-parameter support before
the conditional densities are averaged. IWMDE instead uses
\$\$\widehat{p}\_{I}(\theta^{\*}\mid y)=\frac{1}{S}\sum\_{i=1}^{S}
w(\theta\_{i}\mid\xi\_{i})\frac{q(\theta^{\*},\xi\_{i})}
{q(\theta\_{i},\xi\_{i})},\$\$ where \\w(\theta\mid\xi)\\ is a
normalized, fitted conditional weight. For half-bounded targets, IWMDE
fits the distance from the declared boundary and normalizes a Normal
weight over its central 95 percent positive-tail interval, using the
actual computed endpoints. A draw outside that compact interval
contributes zero while retaining its original row and chain. The same
compact family is used when conditional fitting is unavailable. Compact
focal support does not certify a finite joint second moment or MCMC
long-run variance; reported uncertainty is conditional on the fitted
proposal and does not include uncertainty from learning it on the same
draws. Spike-mixture and product-space results additionally retain the
posterior mass of the continuous active branch. These estimators build
on conditional marginal density estimation and importance-weighted
marginal density estimation (Gelfand et al. 1992; Chen 1994) .

qCMDE is the preferred likelihood-aware method when it is supported and
its additional numerical normalization cost is acceptable. IWMDE can be
faster, but is more sensitive to the adequacy and concentration of its
fitted conditional weights. The normal approximation is useful as a
rough check for near-normal interior ordinates; it is not a validation
reference for skewed tails, bounded parameters, or one-sided support
boundaries. Selected weightfunction models near zero cumulative ordered
effects can take more than 10 minutes with qCMDE or IWMDE. qCMDE is the
default for fitted `brma` objects. Use `"KDE"`, `"normal"`, or `"IWMDE"`
explicitly to select another route. For binomial and Poisson GLMMs, the
nuisance state includes the sampled estimate-level heterogeneity effects
and baserates or log-rates. Averaging the normalized row-conditional
focal densities over the posterior draws of this nuisance state targets
the marginal posterior density by an exact identity; the estimate
approximates that expectation with the posterior draws and numerical
normalization, subject to the reported diagnostics. IWMDE is unavailable
for GLMM density estimation because its high-dimensional conditional
weights did not meet bridge-sampling certification. GLMM density curves
and point-null hypotheses require qCMDE. qCMDE and IWMDE evaluate one
scale formula: they are unavailable for
[`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md)
models with component-specific scale formulas (a named `scale` list,
e.g. one scale formula per random component), whose point-null
hypotheses use `density_method = "KDE"`.

For an unrestricted alternative against \\\theta=\theta_0\\, the
Savage-Dickey identity (Dickey 1971; Wagenmakers et al. 2010) gives
\\BF\_{10}=p(\theta_0\mid H_1)/p(\theta_0\mid y,H_1)\\ only when the
null model's nuisance-parameter prior is the conditional prior induced
by the unrestricted model at \\\theta_0\\. Otherwise, fit the
constrained model separately and compare marginal likelihoods, for
example by bridge sampling.

Marginal-means hypotheses are specified on the fitted linear-predictor
scale. Display transformations stored by
[`marginal_means()`](https://fbartos.github.io/RoBMA/reference/marginal_means.md)
do not transform hypothesis constants. Single-model hypotheses and
model-averaged region hypotheses use averaged marginal draws.
Model-averaged point-null hypotheses use alternative-conditioned draws
to avoid null atoms in averaged marginals. For model-averaged objects,
hypotheses that mix point and region events, and point hypotheses on
more than one marginal mean (several levels of a factor, or levels and
scalar marginal means of one parameter or of several), are not
supported. Separate region statements may refer to different parameters,
while each statement refers to the marginal means of one parameter.
Single-model hypotheses use one common averaged posterior and may
combine these events. Point hypotheses require the exact prior ordinate
that BayesTools attaches to the marginal mean, as for fitted objects.

## References

Chen M (1994). “Importance-weighted marginal Bayesian posterior density
estimation.” *Journal of the American Statistical Association*,
**89**(427), 818–824.
[doi:10.1080/01621459.1994.10476815](https://doi.org/10.1080/01621459.1994.10476815)
.  
  
Dickey JM (1971). “The weighted likelihood ratio, linear hypotheses on
normal location parameters.” *The Annals of Mathematical Statistics*,
204–223.
[doi:10.1214/aoms/1177693507](https://doi.org/10.1214/aoms/1177693507)
.  
  
Gelfand AE, Smith AF, Lee T (1992). “Bayesian analysis of constrained
parameter and truncated data problems using Gibbs sampling.” *Journal of
the American Statistical Association*, **87**(418), 523–532.
[doi:10.1080/01621459.1992.10475235](https://doi.org/10.1080/01621459.1992.10475235)
.  
  
Wagenmakers E, Lodewyckx T, Kuriyal H, Grasman R (2010). “Bayesian
hypothesis testing for psychologists: A tutorial on the Savage–Dickey
method.” *Cognitive Psychology*, **60**(3), 158–189.
[doi:10.1016/j.cogpsych.2009.12.001](https://doi.org/10.1016/j.cogpsych.2009.12.001)
.

## See also

[`density_diagnostics()`](https://fbartos.github.io/RoBMA/reference/density_diagnostics.md),
[`marginal_means()`](https://fbartos.github.io/RoBMA/reference/marginal_means.md),
[`plot.brma()`](https://fbartos.github.io/RoBMA/reference/plot.brma.md),
[`bridge_sampler.brma()`](https://fbartos.github.io/RoBMA/reference/bridge_sampler.brma.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- brma(
  yi      = c(0.12, 0.20, 0.05, 0.33, 0.18),
  sei     = c(0.10, 0.12, 0.09, 0.15, 0.11),
  measure = "GEN"
)

# Average-effect hypotheses
hypothesis(fit, "mu != 0 vs mu = 0")
hypothesis(fit, "mu > 0 vs mu = 0")
hypothesis(fit, "mu > 0 vs mu < 0")

# Heterogeneity hypothesis
hypothesis(fit, "tau = 0 vs tau > 0", component = "scale")

dat <- data.frame(
  yi    = c(0.12, 0.20, 0.05, 0.33, 0.18, 0.41),
  sei   = c(0.10, 0.12, 0.09, 0.15, 0.11, 0.14),
  group = factor(c("A", "A", "A", "B", "B", "B"))
)
fit_mod <- brma(
  yi      = yi,
  sei     = sei,
  mods    = ~ group,
  data    = dat,
  measure = "GEN"
)

# Marginal-mean hypothesis for a factor level
emm <- marginal_means(fit_mod, density_method = "qCMDE")
hypothesis(emm, "group[B] > 0 vs group[B] = 0")
} # }
```
