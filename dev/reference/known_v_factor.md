# Declare a diagonal-plus-factor sampling covariance

Constructs an exact known sampling covariance representation \$\$V =
\mathrm{diag}(d) + U U^\mathsf{T}.\$\$ The supplied diagonal and loading
matrix define `V`; no numerical rank is inferred from a materialized
covariance matrix. This representation permits likelihood
implementations to use the declared factor dimension when that route is
supported. Gaussian and selection models retain the same target as a
conventional dense `V` input with the same source settings and, for
best-p-value selection, publication groups.

With `selection_model(known_sampling_variance = "condition")`, the
declared covariance defines the complete retained sampling-error vector,
including both diagonal and factor variation. With
`known_sampling_variance = "integrate"` (the default), the full sampling
error is integrated before selection normalization. The
diagonal-plus-factor split is a numerical representation; it does not
define which sampling sources are conditioned upon.

Ordinary covariance matrices, including results from
[`metafor::vcalc()`](https://wviechtb.github.io/metafor/reference/vcalc.html),
can also be supplied directly. Covariance inputs do not specify
publication groups. Only `weight_rule = "best"` uses publication groups;
declare them separately with
`selection_model(weight_rule = "best", group = ...)`. The default
product rule requires no publication grouping.

A declaration is not required to reach the factor routes. Dependency
blocks whose correlation is constant on nested partitions of their rows,
which is what compound symmetry and nested
[`metafor::vcalc()`](https://wviechtb.github.io/metafor/reference/vcalc.html)
`type`/`obs` structures produce, are recognized from an ordinary matrix
and represented exactly as \\\mathrm{diag}(d) + U U^\mathsf{T}\\. Other
blocks that are exactly a diagonal plus a few loading columns, such as
covariance from overlapping samples, shared control arms, or partially
shared raters, are recovered by a rank search that stops at the first
successful candidate, the kernel rank cap, or the Ledermann bound. It
finds the smallest rank reached by that search, without certifying the
minimum possible rank. The recovered representation must reproduce the
supplied entries to working precision; blocks that do not, such as
autoregressive `phi` structures, keep the general route with the
supplied entries. `V` remains the stored covariance either way. Use
`known_v_factor()` to declare a factor that cannot be recovered from the
matrix, or to name its sources.

Matrix inputs use base R's numerical symmetry tolerance. Accepted
matrices are stored symmetrically by averaging corresponding
off-diagonal entries, leaving the supplied R object and its diagonal
unchanged. Positive semidefiniteness is judged from the standardized
eigenvalues, with values inside the eigensolver's roundoff envelope
treated as zero. Exactly singular designs are therefore accepted and
factorized to their true rank, including
[`metafor::vcalc()`](https://wviechtb.github.io/metafor/reference/vcalc.html)
output with `rho = 1`, whose diagonal and off-diagonal reach working
precision along different rounding paths. A singular covariance still
requires model structure covering its null space; see
[`brma.mv()`](https://fbartos.github.io/RoBMA/dev/reference/brma.mv.md).
The stored covariance uses this symmetric representation and retains
input type and row mapping; it does not include a separate copy of
unequal input triangles.

## Usage

``` r
known_v_factor(diagonal, loading)
```

## Arguments

- diagonal:

  finite non-negative numeric vector `d`.

- loading:

  finite numeric matrix `U` with one row per element of `diagonal`. A
  zero-column matrix is allowed and represents a diagonal covariance.

## Value

A `RoBMA_known_v_factor` object accepted as the `V` or `V_new` argument
of multivariate-model interfaces.
