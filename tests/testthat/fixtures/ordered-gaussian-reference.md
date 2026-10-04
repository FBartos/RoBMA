The binary fixture preserves all 248 original GH31 nodes (three and ninety
observations, ordinary and fixed-step selected, positive and negative effects)
and all nine original share-boundary controls, including all 38 conditioning
witnesses. Totals, shares, model inputs, primitive Gamma rows, priors, selection
context, and replay mean/SD values retain their original IEEE doubles. Nodes are
stored rather than regenerated or recovered from a rounded CSV.

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
its stored SD. Stored signed mean and variance perturbations and declared
residuals come from the same independent high precision calculation. Tests
recompute the Gaussian loss bound
sum p*(abs(y-m)*abs(delta)/v + delta^2/(2*v)), then add the variance change at
m+delta. Selected cases add the global step-law sensitivities
TV(omega)/(sqrt(2*pi*v)*min(omega)) for the mean and
TV(omega)/(2*v*sqrt(2*pi*e)*min(omega)) for variance, using the smaller endpoint
variance. These conditioning checks retain the original 64*eps*Nobs allowance;
they do not redefine the declared statistical target. No math-tool dependency
is needed at runtime or during tests.
