# Input Data Specification

Shared data-input arguments used by the RoBMA fitting functions.

Normal models use approximate effect-size estimates supplied through
`yi` with either `vi` or `sei`. GLMM models use the raw two-arm count
arguments for binomial (`measure = "OR"`) or Poisson (`measure = "IRR"`)
outcomes.

## Arguments

- yi:

  a vector of effect sizes, or a formula with the effect size on the
  left-hand side and location moderators on the right-hand side (for
  example `yi ~ x1 + x2`). If a formula is supplied, `mods` must not be
  specified.

- vi:

  a vector of sampling variances. Either `vi` or `sei` must be supplied
  for normal models.

- sei:

  a vector of standard errors. Either `vi` or `sei` must be supplied for
  normal models.

- V:

  a known working variance-covariance matrix, a list of block
  variance-covariance matrices, or an exact diagonal-plus-factor
  declaration created by
  [`known_v_factor()`](https://fbartos.github.io/RoBMA/reference/known_v_factor.md),
  used by
  [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md).
  Ordinary covariance matrices from
  [`metafor::vcalc()`](https://wviechtb.github.io/metafor/reference/vcalc.html)
  can be supplied directly. Product selection requires no publication
  groups. For best-p-value selection, declare groups separately with
  `selection_model(weight_rule = "best", group = ...)`. Matrix symmetry
  uses base R's numerical tolerance. Accepted matrices are stored
  symmetrically by averaging off-diagonal pairs; the supplied object and
  its diagonal are unchanged. Dependency blocks that are block constant
  on nested partitions are represented exactly as a diagonal-plus-factor
  structure without a
  [`known_v_factor()`](https://fbartos.github.io/RoBMA/reference/known_v_factor.md)
  declaration. Dimensions and row alignment are checked on the complete
  input; when rows are dropped (`subset` or missing values), the
  numerical covariance checks apply to the retained rows only, so
  acceptance says nothing about the entries of the dropped rows.

- known_v_parameterization:

  known-`V` backend used by
  [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md).
  `"auto"` chooses an exact backend when feasible; `"latent"` uses a
  latent `D + BB'` decomposition; `"whitened"` uses an eigen-rotated
  normal likelihood; `"block_mvn"` uses an exact native block
  multivariate-normal likelihood.

- weights:

  an optional vector of finite positive likelihood weights. For
  normal/effect-size models, each weight powers the estimate likelihood.
  For constructors with GLMM raw-count input, each weight powers the
  paired two-arm likelihood for one study. An integer weight \\w\\
  repeats that conditional likelihood contribution \\w\\ times, with the
  declared latent-effect and grouping structure held fixed. Positive
  fractional likelihood powers are also supported.

- ni:

  an optional vector of sample sizes. Used for `measure = "GEN"` or when
  estimating the unit information standard deviation.

- mods:

  an optional matrix, data.frame, or formula specifying location
  moderators (meta-regressors). Formula input is evaluated in `data` and
  supports explicit columns and standard formula operators, including
  `^` as R's crossing operator with a literal whole exponent of at least
  2: `(x + z)^2` expands to `x + z + x:z`, and `x^2` is not a quadratic
  term (`x + x^2` is `x`). Precompute transformed predictors, such as a
  squared predictor, as columns in `data`; inline calls such as
  [`I()`](https://rdrr.io/r/base/AsIs.html),
  [`poly()`](https://rdrr.io/r/stats/poly.html), spline constructors,
  and user functions are not supported.

- scale:

  an optional matrix, data.frame, or formula specifying scale predictors
  for location-scale models. Formula input is evaluated in `data` and
  follows the same explicit-column formula grammar as `mods`.
  Random-formula models also accept a named list targeting top-level
  random components or concrete random-effect blocks; see
  [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md).

- random:

  an optional formula or list of formulas specifying BayesTools
  random-effect terms for
  [`brma.mv()`](https://fbartos.github.io/RoBMA/reference/brma.mv.md).
  Use [random-effect formula structure
  tags](https://fbartos.github.io/RoBMA/reference/random_effect_formula_tags.md)
  such as [`diag()`](https://rdrr.io/r/base/diag.html), `us()`/`un()`,
  `cs()`, `hcs()`, `ar1()`/[`ar()`](https://rdrr.io/r/stats/ar.html),
  `har()`, or `car()` inside the formula. Plain `(expr | group)` syntax
  uses an unstructured random-coefficient block; `expr || group` uses a
  diagonal block. The coefficient-formula and structured-index tag
  families have different left-side semantics; see the linked
  documentation. A bare formula or unnamed one-entry list has no
  redundant component prefix. An explicitly named one-entry list retains
  its name; lists with two or more entries generate missing names as
  `component 1`, `component 2`, and so on.

- cluster:

  an optional vector of cluster identifiers for multilevel
  meta-analysis.

- data:

  an optional data frame containing the variables.

- slab:

  an optional vector of study labels.

- subset:

  an optional logical or numeric vector specifying a subset of data to
  be used. Subsetting precedes missing-row removal; retained-row outcome
  and predictor checks follow these steps. Type, shape, and existing
  extraction checks can still reject malformed inputs before subsetting.
  Row selection (by `subset` or missing-row removal) drops factor levels
  left without rows from `mods`, `scale`, and random-effect grouping
  factors. The index variable of a structured random-effect term
  (`cs()`, `hcs()`, `ar1()`, `har()` with one index variable) keeps its
  declared levels, which define the covariance structure; the variables
  of a multi-variable `cs(a + b | group)` index follow the ordinary
  rule, because that index is built from their observed combinations. If
  an index variable that keeps its levels is also used in `mods` or
  `scale` and loses one of its levels to the selection, the fit stops
  with an error of class `RoBMA_structured_index_level_removed` (also
  `RoBMA_data_error`); remove the level from the data with
  [`droplevels()`](https://rdrr.io/r/base/droplevels.html) or use a
  separate column for the moderator.

- measure:

  a character string specifying the effect size measure.
  Normal/effect-size constructors require an explicit value and support
  `"SMD"`, `"ZCOR"`, `"RR"`, `"OR"`, `"HR"`, `"RD"`, `"IRR"`, and
  `"GEN"`. Use `"GEN"` only for general effect sizes without a known
  unit information standard deviation. GLMM raw-count constructors
  support only `"OR"` and `"IRR"` and default to `"OR"`.

- effect_direction:

  direction used by publication-bias adjustments. `"positive"` assumes
  statistically significant positive estimates are more likely to be
  selected; `"negative"` mirrors the selection direction; `"detect"`
  infers the direction from the fitted data.

- ai:

  a vector of the number of events in the treatment or experimental
  group for binomial GLMM models.

- bi:

  a vector of the number of non-events in the treatment or experimental
  group for binomial GLMM models.

- ci:

  a vector of the number of events in the control group for binomial
  GLMM models.

- di:

  a vector of the number of non-events in the control group for binomial
  GLMM models.

- n1i:

  a vector of the sample size in the treatment or experimental group. If
  omitted for binomial GLMMs, it is computed as `ai + bi`.

- n2i:

  a vector of the sample size in the control group. If omitted for
  binomial GLMMs, it is computed as `ci + di`.

- x1i:

  a vector of the number of events in the treatment/experimental group
  (for Poisson data).

- x2i:

  a vector of the number of events in the control group (for Poisson
  data).

- t1i:

  a vector of the person-time in the treatment/experimental group.

- t2i:

  a vector of the person-time in the control group.
