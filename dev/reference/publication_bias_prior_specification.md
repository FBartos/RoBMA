# Publication-bias prior specification

Publication-bias prior distributions are specified separately from the
meta-analytic prior distributions documented in
[`prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/prior_specification.md).
They define selection-model weights, PET/PEESE regression coefficients,
or the publication-bias model space in
[`RoBMA`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md) and
[`RoBMA.mv`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.mv.md).

## Arguments

- prior_bias:

  prior distribution for publication-bias adjustment. For
  [`bselmodel`](https://fbartos.github.io/RoBMA/dev/reference/bselmodel.md),
  this is usually a weightfunction prior created by
  [`prior_weightfunction`](https://fbartos.github.io/RoBMA/dev/reference/prior_weightfunction.md).
  For [`bPET`](https://fbartos.github.io/RoBMA/dev/reference/bPET.md)
  and
  [`bPET.mv`](https://fbartos.github.io/RoBMA/dev/reference/bPET.mv.md),
  use
  [`prior_PET`](https://fbartos.github.io/RoBMA/dev/reference/prior_PET.md).
  For
  [`bPEESE`](https://fbartos.github.io/RoBMA/dev/reference/bPEESE.md)
  and
  [`bPEESE.mv`](https://fbartos.github.io/RoBMA/dev/reference/bPEESE.mv.md),
  use
  [`prior_PEESE`](https://fbartos.github.io/RoBMA/dev/reference/prior_PEESE.md).
  For [`RoBMA`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md)
  and
  [`RoBMA.mv`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.mv.md),
  this can be a single publication-bias prior distribution or a list of
  publication-bias prior distributions. In the single-model
  bias-adjustment constructors, omitted or `NULL` uses the corresponding
  default prior distribution.

- prior_bias_null:

  prior distribution(s) for null publication-bias component(s) in
  [`RoBMA`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md) or
  [`RoBMA.mv`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.mv.md),
  usually
  [`prior_none()`](https://fbartos.github.io/RoBMA/dev/reference/prior_none.md).
  A single prior distribution object creates one null component, a list
  creates multiple null components, and `NULL` or `FALSE` omits the null
  publication-bias component.

- model_type:

  character string specifying predefined publication-bias model
  ensembles for
  [`RoBMA`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md) and
  [`RoBMA.mv`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.mv.md).
  One of `"PSMA"`, `"6w"`, `"2w"`, or `"PP"`.

- steps:

  numeric vector of one-sided p-value cut points for the default
  [`bselmodel`](https://fbartos.github.io/RoBMA/dev/reference/bselmodel.md)
  weightfunction prior. If `prior_bias` is supplied, the prior
  distribution carries its own p-value cut points.

## Details

### Single-model bias-adjustment priors

[`bselmodel`](https://fbartos.github.io/RoBMA/dev/reference/bselmodel.md),
[`bselmodel.mv`](https://fbartos.github.io/RoBMA/dev/reference/bselmodel.mv.md),
[`bPET`](https://fbartos.github.io/RoBMA/dev/reference/bPET.md),
[`bPET.mv`](https://fbartos.github.io/RoBMA/dev/reference/bPET.mv.md),
[`bPEESE`](https://fbartos.github.io/RoBMA/dev/reference/bPEESE.md), and
[`bPEESE.mv`](https://fbartos.github.io/RoBMA/dev/reference/bPEESE.mv.md)
fit one publication-bias adjustment at a time. The `prior_bias` argument
must match the fitted bias-adjustment type.

|  |  |  |
|----|----|----|
| Constructor | Prior constructor | Default prior distribution |
| [`bselmodel()`](https://fbartos.github.io/RoBMA/dev/reference/bselmodel.md) | [`prior_weightfunction`](https://fbartos.github.io/RoBMA/dev/reference/prior_weightfunction.md) | one-sided cumulative weightfunction with `steps = 0.025` |
| [`bPET()`](https://fbartos.github.io/RoBMA/dev/reference/bPET.md) | [`prior_PET`](https://fbartos.github.io/RoBMA/dev/reference/prior_PET.md) | positive Cauchy centered at 0 |
| [`bPEESE()`](https://fbartos.github.io/RoBMA/dev/reference/bPEESE.md) | [`prior_PEESE`](https://fbartos.github.io/RoBMA/dev/reference/prior_PEESE.md) | positive Cauchy centered at 0, with UISD-adjusted scale |
| [`bPET.mv()`](https://fbartos.github.io/RoBMA/dev/reference/bPET.mv.md) | [`prior_PET`](https://fbartos.github.io/RoBMA/dev/reference/prior_PET.md) | same as [`bPET()`](https://fbartos.github.io/RoBMA/dev/reference/bPET.md) |
| [`bPEESE.mv()`](https://fbartos.github.io/RoBMA/dev/reference/bPEESE.mv.md) | [`prior_PEESE`](https://fbartos.github.io/RoBMA/dev/reference/prior_PEESE.md) | same as [`bPEESE()`](https://fbartos.github.io/RoBMA/dev/reference/bPEESE.md) |

The default PET prior distribution uses
`RoBMA.get_option("default_bias_PET.scale")`. The default PEESE prior
distribution uses `RoBMA.get_option("default_bias_PEESE.scale")` after
rescaling to the analyzed effect-size measure. The default
weightfunction prior distribution uses
`RoBMA.get_option("default_bias_weightfunction.alpha")` for each
cumulative-Dirichlet alpha parameter.

### Model-averaged publication-bias priors

[`RoBMA`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md) and
[`RoBMA.mv`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.mv.md)
average over publication-bias components. By default, `prior_bias_null`
is
[`prior_none()`](https://fbartos.github.io/RoBMA/dev/reference/prior_none.md)
and `model_type` determines the alternative components:

|          |                                      |
|----------|--------------------------------------|
| `"PSMA"` | six weight functions, PET, and PEESE |
| `"6w"`   | six weight functions                 |
| `"2w"`   | two two-sided weight functions       |
| `"PP"`   | PET and PEESE                        |

Custom `prior_bias` replaces the preset alternative components. If
`prior_bias` is omitted, `model_type` still supplies the default
alternative components even when `prior_bias_null` is customized or
omitted. Setting `prior_bias = NULL` or `prior_bias = FALSE` omits all
alternative bias-adjustment models. Setting `prior_bias_null = NULL` or
`prior_bias_null = FALSE` omits the no-bias component.

Each publication-bias component has a prior model weight. User-specified
components receive equal weights unless their prior objects set
`prior_weights`. The `"2w"` and `"6w"` presets split the alternative
bias mass equally across their weight functions, `"PP"` splits it
equally across PET and PEESE, and `"PSMA"` assigns half of the
alternative bias mass to the six weight functions combined and one
quarter each to PET and PEESE.

Publication-bias prior distributions are not rescaled by
`rescale_priors`. The exception is the default PEESE prior distribution,
whose scale is adjusted to the effect-size measure via the
unit-information standard deviation.

## See also

[`prior_weightfunction`](https://fbartos.github.io/RoBMA/dev/reference/prior_weightfunction.md),
[`prior_PET`](https://fbartos.github.io/RoBMA/dev/reference/prior_PET.md),
[`prior_PEESE`](https://fbartos.github.io/RoBMA/dev/reference/prior_PEESE.md),
[`prior_none`](https://fbartos.github.io/RoBMA/dev/reference/prior_none.md),
[`RoBMA_prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA_prior_specification.md),
[`prior_specification`](https://fbartos.github.io/RoBMA/dev/reference/prior_specification.md)
