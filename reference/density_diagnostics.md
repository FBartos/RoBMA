# Extract qCMDE/IWMDE Density Diagnostics

Extracts the compact reliability diagnostics attached to a qCMDE/IWMDE
hypothesis result. The diagnostics describe the posterior ordinate
computation. They do not include uncertainty in the prior ordinate or,
for IWMDE, uncertainty from estimating its conditional weight function.

## Usage

``` r
density_diagnostics(object, ...)
```

## Arguments

- object:

  an object returned by
  [`hypothesis()`](https://fbartos.github.io/RoBMA/reference/hypothesis.md)
  with `density_method = "qCMDE"` or `density_method = "IWMDE"`, or an
  `RoBMA_density_ordinate_error` caught from a failed or rejected
  point-ordinate computation.

- ...:

  unused.

## Value

A data frame of class `RoBMA_density_diagnostics`, with one row per
computed or failed requested point ordinate. Columns identify the
estimator, parameter, requested and evaluated values, schema/source
provenance, row counts, active mass, relative MCSE, ESS, largest
contribution share, finite terms, normalization, truncation, and
quadrature checks, the qCMDE pilot-gate columns (the nested
normalization grids have no refinement sequence to stop, so
`pilot_gate_stopped` is `FALSE` and `pilot_bulk_ess` is `NA`),
fixed-sampling state, policy thresholds, weight fallbacks, status, and
warnings.

The exact columns, in order, are `schema_version`, `algorithm_version`,
`source_fingerprint`, `estimator`, `internal_method`,
`working_distribution`, `normalization_basis`, `ordered_error_status`,
`block_numerator_relative_error`, `block_denominator_relative_error`,
`block_max_order`, `block_max_panels`, `block_native_nodes`,
`density_method`, `parameter`, `level`, `requested_value`,
`evaluation_value`, `requested_samples`, `achieved_row_budget`,
`eligible_rows`, `evaluated_rows`, `finite_terms`, `active_mass`,
`relative_mcse`, `active_branch_relative_mcse`,
`active_mass_relative_mcse`, `mixture_mcse_type`,
`sampling_relative_mcse`, `sampling_fraction`, `mcmc_uncertainty_scope`,
`sampling_uncertainty_type`, `ess`, `max_weight_share`,
`normalization_relative_error`, `stability_metric`,
`stability_relative_error`, `ordinate_relative_change`,
`lattice_resolution_relative_error`, `normalization_truncation`,
`normalization_truncation_status`, `truncation_ordinate_bound`,
`quadrature_relative_change`, `pilot_gate_stopped`, `pilot_bulk_ess`,
`target_relative_mcse`, `stability_warning_threshold`,
`stability_rejection_threshold`, `quadrature_warning_threshold`,
`quadrature_rejection_threshold`, `warning_relative_mcse`,
`warning_min_finite_terms`, `warning_min_ess`,
`warning_max_weight_share`, `all_rows_used`, `target_met`,
`precision_target_met`, `sampling_target_met`, `bf_grade_met`,
`n_weight_fallbacks`, `weight_fallback_reasons`, `status`, `ordinate`,
`log_ordinate`, `numerical_status`, `failure_reason`, and `warnings`.

## Details

Point ordinates evaluate one fixed posterior-row sample whose size is
chosen by `samples` before any ordinate contributions are computed.
Precision, effective-sample-size, and contribution-concentration
diagnostics never determine whether that finite estimate is returned. An
unmet `target_relative_mcse` produces a warning; increase `samples` or
set it to `Inf` to use every eligible row. The table records the
requested and achieved sample sizes, `all_rows_used`, whether the
precision target was met, and whether the ordinate passed the
independent numerical-stability gates. Finite row budgets use a
reproducible simple random sample without replacement. For a purely
continuous ordinate, `relative_mcse` and `ess` describe the selected
continuous posterior-row sequence. For a mixed point/continuous ordinate
evaluated on every active continuous row, they use the full
conditioned-chain sequence, with zero contributions for inactive rows,
and therefore include active-mass uncertainty and its covariance with
continuous contributions. With a subset of active rows,
`active_branch_relative_mcse` describes the selected active-row
contribution sequence and `active_mass_relative_mcse` comes from the
complete binary indicator chain. `relative_mcse` is their conservative
worst-correlation delta upper bound on the unconditional density scale.
It is suitable as a descriptive precision diagnostic but does not
certify the unobserved full-mixture covariance, so `bf_grade_met`
remains false. `sampling_relative_mcse` separately reports the
finite-population simple-random-sampling component conditional on the
fitted rows (`sampling_uncertainty_type = "finite_population_srswor"`).
The sampling design is exact, but the standard error is an estimate: the
usual without-replacement formula applied to the sample variance of the
selected contributions, and zero when every eligible row is used. It is
not added to the selected-row sequence MCSE, which already describes the
actual selected estimator sequence.

These fields concern the continuous density ordinate. At a location that
is also a posterior atom, the atom probability remains a separate
discrete mass and is never added to the density ordinate.

`ordinate` and `log_ordinate` describe the posterior density at the
requested point, not the model's marginal likelihood. `log_ordinate` is
retained only when actually computed; it is not reconstructed for older
results. `numerical_status` distinguishes `finite`, `underflow`,
`overflow`, `arithmetic_failure`, `invalid_density`,
`zero_without_log_evidence`, `nonfinite`, and `not_computed`.
Underflow/overflow requires finite log evidence and confirmation that
exponentiating that log loses range. A zero without this evidence is not
classified as a mathematical zero. A small density alone does not
establish unreliability.

Failed requested points retain separate diagnostic records with a
`failure_reason`, including failures occurring before a final ordinate
was obtained. `status = "unavailable"` distinguishes these failures from
a finite positive ordinate rejected by accuracy diagnostics.
Ordinary-scale precision diagnostics lost to underflow/overflow are
unavailable, not zero error or proof of convergence. Retaining a log
estimate does not authorize a Bayes factor: the existing finite-positive
ordinate requirement is unchanged. Legacy tables remain inspectable with
unknown new fields set to `NA`. A hypothesis result without attached
diagnostics, such as one whose point ordinates were estimated with
`density_method = "KDE"` or one without a point hypothesis, stops with
the class `RoBMA_density_diagnostics_unavailable`. Density ordinate,
diagnostic-unavailable, plot, and construction
(`RoBMA_density_construction_error`) conditions share the
`RoBMA_density_error` parent while retaining their specific condition
classes.

The reliability policy warns when relative MCSE is at least 5 percent,
ESS is below 100, the largest contribution share is at least 20 percent,
or fewer than 100 finite contributions remain. These same-sample
diagnostics do not reject a finite fixed-design estimate. Their warning
thresholds are returned as columns rather than being implicit. For a
descriptive mixed density curve using the conservative partial-row
bound, a finite ESS below the usual curve rejection threshold is
retained as a warning when the quantitative mixture-scale MCSE checks
pass.

qCMDE's full method diagnostics distinguish
`pilot_normalization_integral` from `final_normalization_integral`;
IWMDE uses `support_grid_normalization_integral`. This compact accessor
reports the relative difference of the final integral from the active
mass as `normalization_relative_error`. For qCMDE it is a bookkeeping
check, not a discretization error: rows normalized on the lattice are
normalized by the same trapezoid rule, on the same grid, that integrates
them, so they sum to the active mass up to roundoff by construction, and
rows normalized in closed form contribute their analytic mass within the
range. The discretization and truncation errors of the lattice
normalizers are reported by `ordinate_relative_change`,
`lattice_resolution_relative_error`, and `truncation_ordinate_bound`.
The active policy check is identified by `stability_metric`: qCMDE
checks its posterior ordinate error estimate
(`ordinate_error_estimate`), whereas IWMDE checks normalization mass.
The qCMDE estimate adds three parts. `ordinate_relative_change` is the
ordinate's movement between the selected normalization grid and the
coarser nested grid over the same range. It estimates the discretization
error of the lattice normalizers; it is an estimate, not a bound,
because a row narrower than the lattice step can escape both grids
alike. Rows whose lattice normalization is under-resolved also
contribute an ordinate-share weighted resolution estimate
(`lattice_resolution_relative_error`), using the kernel width or
nested-grid log curvature at each row's maximum. This component and
their sum are estimates, not rigorous bounds: a sharp prior core can be
narrower than its likelihood kernel. Rows whose Gaussian likelihood
kernel meets a normal prior are normalized by the closed-form mass of
their support when their evaluated values are those of that closed form:
their normalizer is the same on every grid, so they add no
discretization change, lattice-resolution component or truncation.
`truncation_ordinate_bound` is \\t / (1 - t)\\ for the largest fraction
\\t\\ of any row's conditional mass outside the normalization range
(`normalization_truncation`): a row normalizer missing that fraction
overstates the row's density by at most this relative amount.
`normalization_truncation_status` states how \\t\\ is known: `"exact"`
for rows whose Gaussian likelihood kernel meets a normal prior,
`"bound"` for a Gaussian kernel with another proper prior, and
`"estimate"` (from the exponential envelope of the tails at the range
ends) otherwise; it reports the weakest status among the rows. The range
covers every row's central `normalization_prob` interval, so \\t\\ is at
most `1 - normalization_prob` (for estimated rows, as estimated once
their extension reaches its target) and `truncation_ordinate_bound` is
at most `(1 - normalization_prob) / normalization_prob`, about 0.1
percent at the default `0.999`. The extension of an estimated row stops
at its limit (four initial range widths per side and 40 passes), at a
value the joint density cannot evaluate, or for a tail that does not
decay; the reported `normalization_truncation` and
`truncation_ordinate_bound` can then exceed these values. Setting
`normalization_prob` closer to 1 widens the range and lowers the
truncation; increasing `normalization_points` refines the grid. The
table also reports the stability and adaptive-quadrature
warning/rejection thresholds.

Ordered effects use normalized conditional-source contributions for both
requested methods. Regular exact offsets can use a positive total/share
block. The complete numerical ratio error must be at most 0.025 for
every evaluated row, separately from the MCMC precision diagnostics. The
columns `internal_method`, `working_distribution`, and
`normalization_basis` identify this computation while `estimator`
retains the requested method. `ordered_error_status` preserves its
estimated error status. Block numerator and denominator errors, maximum
quadrature order/panel count, and actual native normalizer node count
are reported in the `block_*` columns. Exact offsets classified as
nonfinite or unsupported are excluded only from density curves; explicit
point requests retain their structural or precision refusal. The ordered
per-row limit of 0.025 is separate from the generic stability
warning/rejection thresholds reported in this table.
Source-normalization refusals have class
`RoBMA_ordered_normalization_error` and retain
`iwmde_ordered_normalization_error` as a compatibility ancestor. Their
condition's `proof` field contains the numerical integration evidence.
When only some ordered evaluations lack certified contribution
variation, their MCMC error and ESS stay unavailable while other
evaluations retain their metrics. A printed table footnote reports the
`partial` uncertainty status and its reason; the original selected-row
scope stays in its column. Such curves fail the MCMC availability gate,
and an affected point cannot meet its precision target. If every
evaluation is affected, the full unavailable scope is retained. A finite
log MCSE whose absolute value underflows remains a separate
representation limitation.

## See also

[`hypothesis()`](https://fbartos.github.io/RoBMA/reference/hypothesis.md)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- hypothesis(
  fit,
  hypothesis = "mu != 0 vs mu = 0"
)
density_diagnostics(result)

failed <- tryCatch(
  hypothesis(
    fit,
    hypothesis     = "mu != 0 vs mu = 0",
    density_method = "IWMDE"
  ),
  RoBMA_density_ordinate_error = identity
)
density_diagnostics(failed)
} # }
```
