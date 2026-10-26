# Maintainer decisions and rationale

These are the final decisions from the 2026-10-25/26 behavior review, including
the maintainer's later clarifications. They supersede conflicting intermediate
audit recommendations. Fix bugs while preserving these statistical targets;
changing a target or restoring a superseded convention requires a new maintainer
decision. Record a later replacement here with its reason. Implementation details
may improve; the known numerical limitations below are not declared resolved.

## Cook's distance follows the unscaled meta-analytic convention

**Decision.** Report `Delta' V_mu^+ Delta`, without dividing by fixed-effect rank.
Keep the Bayesian full-posterior covariance and PSIS deletion calculation.

**Why.** The maintainer explicitly wanted the convention used by metafor for
meta-analytic Cook's distance. Dividing by rank changes the reported scale
systematically; multiplying the reference values in tests would conceal the
disagreement. This is an intentional departure from RoBMA 4.0.0, not a release
compatibility regression to repair. Matching the convention does not make PSIS
posterior deletion identical to a classical refit.

`0daab96c` supersedes the rank-divisor restoration in `d5b0f0e0`. See
[model interfaces](model-interface.md), `R/cooks.distance.R`, and
`tests/testthat/test-02-influence.R`.

## One posterior-based VIF interface

**Decision.** `vif()` evaluates the design-based GLS VIF/GVIF for each posterior
draw, then reports posterior means, with optional `probs` and `include_samples`.
Use this rule across simple, location-scale, multilevel and known-V models.
Average GVIF and its adjusted root separately. Keep empirical posterior
coefficient correlation as a separate diagnostic.

**Why.** Users expect a Bayesian summary under the familiar diagnostic name,
consistent with posterior-averaged `hatvalues()`. A plug-in at mean heterogeneity,
VIF of the mean covariance, and mean draw-wise VIF are different nonlinear
functionals. Switching between them by model class makes equivalent modeling
interfaces behave unexpectedly. Empirical coefficient covariance is another
target: it also reflects prior shrinkage and must not replace the design
diagnostic merely because it is readily available.

The initial suggestion of separate `vif()`/`bvif()` was withdrawn. The maintainer
clarified that `funnel()`/`bfunnel()` exists for a substantial selection-model
speed trade-off. Measured VIF evaluation for 15,000 draws took about 1-1.3 seconds
in the ordinary examples, insufficient justification for a parallel API. This
is not a universal runtime guarantee: retain explicit deterministic thinning
and its warning for expensive designs. No blanket b-prefixed copies of influence
or residual diagnostics were approved, nor were new hatvalue intervals implemented.

Location terms retain their fitted design; PET/PEESE diagnostics condition on
their respective design class, with absent-class draws unavailable rather than
zero. See [model interfaces](model-interface.md), `R/vif.R`, and
`tests/testthat/test-02-vif.R`. Implemented in `f0a361eb`; test-reference repairs
are `11028f60` and `d37f4cdf`. The older simple-model plug-in restoration is
superseded, even if it happens to be closer to a metafor value in one example.

## Integrate the powered likelihood for weighted Gaussian cluster LOO

**Decision.** For supported non-selection Gaussian `brma(cluster=..., weights=...)`,
integrate the focal powered likelihood over the cluster effect conditional on
retained weighted rows. Preserve the power-dependent normalizing constant.
Estimate deletion removes the row's whole power, not one replicated copy.

**Why.** Conditional and integrated importance sampling can estimate the same
deleted-row target, but their finite-sample PSIS stability differs. A weight
slightly different from one should not force the Gaussian model back to sampled
local effects when an exact integral is available. Multiplying the old unweighted
marginal log score by the weight is incorrect: `E[f^w]` generally differs from
`E[f]^w`. Likelihood powers are also distinct from sampling inverse variances.

If the original within variance is `d`, then
`N(y; mu, d)^w = (2*pi*d)^((1-w)/2) / sqrt(w) * N(y; mu, d/w)`,
where the Normal's last argument denotes variance. The constant matters for
the score. The observation CDF and predictive residual moments still use `d`,
plus retained cluster uncertainty, not `d/w`. Combining integrated PSIS weights
with full-data sampled local effects would retain information from the deleted
row; update those consumers together. Global influence moments may use the
integrated weights directly.

The exact Gaussian helper avoids general quadrature and dense covariance
inversion; a representative 6,000-draw/12-row call cost 0.16 seconds CPU versus
0.02 seconds conditionally. It preserves unweighted routes. This does not promise
universal Pareto-k improvement or authorize approximate selection/GLMM integration.
Public `brma.mv()` still rejects weights. See
[predictive targets](multivariate-targets.md), `R/log-lik-weighted-cluster.R`,
and `tests/testthat/test-00-weighted-cluster-loo.R` and
`test-02-weighted-cluster-loo.R` (including exact Gaussian deletion references).
Implemented in `2559d084`.

Old stored weighted estimate-unit LOO/WAIC results require `add_loo()`/`add_waic()`
recomputation. Their integration marker prevents mixing old weights with new
consumers; it is not a new model-comparison target. This post-fit change does not
require a model refit. Distinguish it from BayesTools' parameter-map change,
which did require refitting older development fits and prompted the separately
authorized vignette-cache refresh.

## Pooled heterogeneity and RMS summaries intentionally differ

**Decision.** Retain the arithmetic mean of study-specific SDs for
`pooled_heterogeneity()`, including arithmetic averaging of random-formula scale
sources. Pooled prediction intervals use this pooled quantity. Parameter/design
summaries and the defined `influence()` tau.del calculations retain their RMS
conventions and availability limits.

**Why.** These functions answer different questions. An average SD describes
the pooled quantity; RMS summarizes variance over the specified design. The
maintainer explicitly reconfirmed that distinction. Pooled prediction intervals
are not mixtures over studies; composite random designs make such a replacement
confusing. Do not harmonize the outputs simply because their names are similar.
Random-formula pooling retains its average expanded-design calculation, not an
invented arithmetic average of every final row-level random SD.

The full definition and earlier history are in the pooled-functions section of
[model interfaces](model-interface.md). Relevant checks include
`tests/testthat/test-00-pooled-prediction.R` and `test-00-influence-tau-deletion.R`.

## Loose selection tolerances do not invalidate deterministic calibration

**Decision.** Product deterministic ladder/refinement acceptance checks use
`min(relative_tolerance, 0.005)` where specified in
[selected-normal kernels](selected-normal.md). Preserve default/tighter behavior
and other routes' policies; use the existing diagnostic refusal if the available
budget cannot produce an accepted result.

**Why.** Adversarial examples accepted errors of about 5.1% at a requested 2%
tolerance and nearly 100% at 5%. The deterministic acceptance checks had only
been calibrated at the tighter level; a looser user request did not validate
their coarser quadrature. Capping these checks fixed the reproduced cases while
22 default-result probes stayed bit-identical and matched timing showed no
material default slowdown. Implemented in `208a3af5`; see
`tests/testthat/test-00-selection-certificate-acceptance.R` and the related
selection-certificate tests.

**Rank-one follow-up.** Fourteen saved accepted errors all came from rank-one
product fallbacks, including ten still outside tolerance at actual default
budgets. After the existing fast rules decline, these integrals now use bounds
from exact row groups and nonnegative excess-weight bin extrema over the full
Gaussian line. The actual returned midpoint must meet tolerance against both
endpoints. The existing maximum point count and scramble count set a scalar
probability work ceiling; numerical defaults are unchanged. Exhausted work or
arithmetic ambiguity refuses, with no QMC rescue. Outward arithmetic and explicit
working-precision allowances are tested engineering safeguards, not a universal
libm theorem. First matched four-round release-DLL measurements gave median
CPU seconds before/after of .655/.205 for the fourteen outliers, 1.215/.380 for
the mixed 22 controls and 1.625/.325 for 21 heterogeneous controls. The one
explicitly refused mixed control already failed the baseline diagnostic.
These synthetic measurements do not promise universal speed or availability.

**Remaining limitation.** Higher-rank and best-rule randomized QMC retains an
empirical error diagnostic, not a guaranteed bound. Do not describe accepted
estimates as certified solely because the empirical MCSE passes.
