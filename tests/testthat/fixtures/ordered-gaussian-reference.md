The binary fixture preserves all 248 original GH31 nodes (three and ninety
observations, ordinary and fixed-step selected, positive and negative effects)
and all nine original share-boundary controls, including all 38 conditioning
witnesses. Totals, shares, model inputs, primitive Gamma rows, priors, selection
context, and captured replay mean/SD values retain their original IEEE doubles.
Nodes are stored rather than regenerated or recovered from a rounded CSV.

Declared-target constants were independently computed on 2026-10-04 with
mpmath 1.3.0 at 100 decimal digits and python-flint 0.9.0 / Arb at 384 bits;
both use the exact binary input values. The Gaussian target is
sum_i p_i * [-log(2*pi*v_i)/2 - (y_i-m_i)^2/(2*v_i)], where
m_i = mu0_i + S*mu_c_i + (mu_a_i + S*mu_b_i)*T.
GH cases include the standard-Normal total prior once. Selected cases add the
original powered observed log weights and subtract powered log normalizers
computed from independent erf/erfc interval probabilities at the original signed
step boundaries. The persisted Gamma/normalization/design recipe independently
gives the same real means. No package likelihood or Gaussian helper supplied
these constants.

Separate constants use each actual rounded replay mean and the exact square of
its stored SD. These captured-input comparisons additionally run when local
replay mean/SD inputs match those captured values. Declared-target comparisons
always run with the exact same stored allowances.

Each declared mean additionally has a nearest-double midpoint and certified
radius obtained from its exact rational bilinear expression and checked with
Arb at 448 bits. Every midpoint agrees with the independent 100-digit cell
reference. No existing fixture field or node was replaced to add these fields.
Tests use the current replay mean/SD and these enclosures to recompute the
Gaussian loss bound
sum p*(abs(y-m)*abs(delta)/v + delta^2/(2*v)), then add the variance change at
m+delta. Selected cases add the global step-law sensitivities
TV(omega)/(sqrt(2*pi*v)*min(omega)) for the mean and
TV(omega)/(2*v*sqrt(2*pi*e)*min(omega)) for variance, using the smaller endpoint
variance. The actual mean difference includes subtraction rounding and the
certified radius. The exact square of each live SD is enclosed before bounding
its variance difference; the smaller variance is rounded outward. A finite
operation-count margin covers the bound's arithmetic and row reduction. This
input perturbation margin adds no broader numerical-evaluation floor and
leaves the original 64*eps*Nobs likelihood allowance unchanged.

Twelve standalone controls change captured nonzero means or SDs by one IEEE
ulp in either direction, for ordinary, selected and share-boundary witnesses.
Independent 120-digit physical likelihood oracles and the same analytic
perturbation formula verify their shifts. They exercise mismatched replay
inputs without fitting a model or manufacturing a certificate. No math-tool
dependency is needed at runtime or during tests.

The companion `ordered-normal-source-reference.rds` retains the exact original
ordinary and (.999,.9995) truncated literal Normal source cases, plus the
independently checked (-Inf,.9992) far-tail case, data and
primitive Gamma inputs, with independent MP100/Arb384 or MP120
source/prior/posterior-mass and nuisance constants. The source comparison
accounts for the selected-tail representation, endpoint standardization, the finite-range relative
log-CDF libm allowance, exp/subtraction amplification and physical-reference
coordinate rounding. Producer and physical reference are independently checked
against the source constant. These arithmetic allowances are development
verification, not rigorous libm enclosures; the tight complete bound remains
below 1e-9. They do not alter any of the original 257 node allowances.

Actual replay inputs additionally obey fixed affine-mean/variance operation
caps before the adaptive relation is calculated. The finite fixtures use three
Gamma coordinates and short predictor dot products; 64 rounded operations is
conservative for this scope, and deliberate larger perturbations fail the caps.
The independent `normal-tail-reference.rds` preserves all 32 exact interval
controls and 14 complete edge/error/shape records from the MP100/Arb384 probe.
Original complementary-CDF failures remain recorded; the scalar source helper
now reuses the existing tail-selecting block helper. Source diagnostics include
interval conditioning separately from residual arithmetic, preserving the
existing numerical error budget and making the far-tail rounding floor explicit.
