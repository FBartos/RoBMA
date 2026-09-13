# Regression evidence inventory

This inventory identifies the high-information evidence expected for public
post-fit behavior. Structural checks supplement these oracles; merely returning
an object or rendering without error is not sufficient representative evidence.

| Area | Standard evidence | Certification extension |
|---|---|---|
| Residuals and standardized residuals | Numerical metafor comparisons for normal, regression, multilevel, PET-PEESE, selection, and GLMM cases where the estimand exists | Extended interaction, negative-direction, multilevel-regression, and high-draw known-V cases |
| Predictions | Numerical metafor comparisons for fitted and new-data predictions, target-specific known-V identities, complete sampling-error source reconstruction, and singular Gaussian latent projections | Extended model families and high-draw known-V targets |
| Hat values, DFBETAS, VIF, and influence | Numerical metafor or independent GLS/moment oracles for representative supported families | Extended parameterizations and high-draw known-V diagnostics |
| Heterogeneity summaries | Numerical metafor comparisons for ordinary, scale, multilevel, selection, and GLMM models; deterministic AR/CS/CAR correlation moments, intervals, component ownership, and data-frame coercion | Negative-direction and secondary model variants |
| Funnel, Q-Q, radial, regression, forest, and z plots | Analytic funnel PIT-projection, CDF-mixture, and posterior-integration identities plus human-reviewed vdiffr snapshots, including side-by-side metafor views where displays share an estimand | Secondary model families and customization galleries |
| Prior plots | Base and ggplot snapshots for outcome, moderator, and publication-bias priors | Structural component-selection tests cover additional prior families |
| Marginal z densities | Independent one-dimensional integrations for jointly selected correlated Gaussian marginals and conditionally normalized models; nonunit-SE Jacobians, normal/selection branch routing, reflected directions, threshold probabilities, and extrapolated area | Cached Assink marginal known-V regression comparison with independent joint-selection quadrature |
| Multivariate diagnostics | Analytic known-V/known-R data checks and representative forest coverage | Human-reviewed funnel, Q-Q, forest, and regression snapshots from high-draw metafor fixtures |
| Estimate-unit LOO | Analytic Gaussian conditional-density identities and cached-fit target contracts; independent conditioned-sampling deletion integrals that integrate the deleted error conditional on retained remaining errors, including singular integrated random covariance | Five exact observation-deletion refits of Kearon US/HCS and Ishak HAR models in the `loo-exact-refits` case |
| Joint selection likelihood | Independent `mvtnorm` rectangle-probability identities for the dependent-block normalizer, analytic diagonal reduction, selected-response region masses and singular-covariance support, specialized cluster-posterior moments, dense-versus-compiled conditional random-effect algebra, a cached random-formula `bselmodel.mv` fit, and a cached `RoBMA.mv` product-space fit combining unadjusted, selection, PET, PEESE, moderator, and independently gated random branches. The two cached multivariate fixtures explicitly integrate both contextual sources with product weights; their fitted checks cover manual posterior averaging, model enumeration, summary, LOO, bridge availability contracts, prediction, diagnostics, plots, and LOO-PIT residuals. | Larger dependency blocks and high-draw bridge stability checks |
| Selection target binding | All eight source cells with product and best rules, explicit child-prior precedence, publication-group agreement, estimate-level random slopes and correlated known-R ownership, unchanged weight priors and model odds, and full-error covariance/likelihood invariance across dense and factor representations. Completion tests reject impossible hard-selection contexts and non-unit observation weights under conditioned sampling while retaining unit-weight identity. | Target-matched recovery and prior-invariance cases below |
| Whole sampling-error likelihood | Independent Gaussian Matheron source reconstruction, retained and integrated source covariance identities, ordinary-Gaussian cancellation with all sources conditioned, primitive evidence replay, analytic rank-one/deterministic normalizers, independent singular-factor orthants, PSD CDF/tail endpoints, and selected-response support. Fitted source semantics are tested independently of auxiliary backend coordinates. | Conditional recovery and dense posterior oracles below |
| GLMM quadrature and likelihoods | Independent R integrations for representative binomial and Poisson cases | The `numerical-kernels` and `glmm-models` certification cases |
| qCMDE/IWMDE | Analytic identities, scalar/batch parity, provenance, and failure-contract tests | Fitted bridge/ordinate evidence in the `iwmde-qcmde` case |

Removed RoBMA 4.0 GLMM Q-Q snapshots are intentionally not restored because
the current API rejects discrete outcomes until a discrete PIT convention is
defined. Removed radial snapshots for multilevel, GLMM, selection, PET, and
model-averaged objects likewise correspond to explicit unsupported-input errors,
not superficial replacement tests.

Ranks or finite-value checks may replace direct equality only when RoBMA and
metafor target different predictive quantities. Such cases must retain an
independent structural or moment identity and an explanatory case label.

## Simulated selection-model certification

Run `Rscript tools/test-profile.R certification selection-recovery-marg` or
`Rscript tools/test-profile.R certification selection-recovery-cond`.
These are independent, one-hour-bounded cases, excluded from the standard
profile. They fit fresh small models without scenario caches or snapshots.
BayesTools checks the actual fitting workers for matching RoBMA and BayesTools
versions and native DLL contents at JAGS startup. `load_all()` alone
does not update worker packages. Install the tested builds into an isolated
library and set `R_LIBS_USER` to that library before starting the runner.
Each fit uses three parallel chains; running both cases concurrently uses six
fitting workers. The simulator and independent posterior reference are in
`helper-selection-recovery.R`.

Each profile has six fits: four 40-study examples and two single-study dense
posterior-oracle examples. Together the two profiles cover twelve fitted
cases. The recovery result table retains twelve parameter rows per profile;
in the conditional profile, two fixed-effect selection-weight rows test prior
invariance instead of claiming recovery.

| Sampling covariance | Data | Estimated parameters | Reference |
|---|---|---|---|
| Diagonal matrix | 40 independent studies, three estimates each | Marginal: mean and selection weight; conditional fixed-effect: mean with selection-weight prior invariance; nested: mean, selection weight, study and estimate SDs | Known-parameter generation; conditional fixed-effect Gaussian conjugacy and exact Uniform prior moments |
| Structured matrix from `metafor::vcalc()` | Same size, within-study correlation 0.35; publication groups supplied separately | Same fixed and nested variants | Independently constructed full covariance and target-matched generation; conditional fixed-effect Gaussian/prior invariance |
| Fully connected dense matrix with unequal variances and mixed-sign correlations | Three estimates | Mean, with known selection weight; random SDs either zero or fixed positive values | Independent integrated posterior, including its Gaussian no-selection limit |

Fixed recovery examples use a positive mean; nested recovery examples use a
negative mean and negative effect direction. All non-reference priors are
explicit: Normal(0, 1) for the mean, independent Uniform(0, 1) selection weight,
and half-Normal(0, 0.5) component SDs. The two-level random formula has one study
effect and one estimate-specific effect, i.e., a three-level meta-analysis.
Dense examples fix nuisance parameters because a tiny fully connected dataset
cannot identify all variance components and selection strength reliably. Their
posterior-oracle comparison is stronger than merely checking whether a broad
credible interval includes the truth.

The marginal generator accepts/rejects a complete Gaussian block with probability
equal to the product of its row selection weights, redrawing all effects after
rejection. The conditional generator draws the complete sampling error
`e ~ Normal(0,V)` and study effect once, then accepts/rejects each estimate's
independent true-effect deviation. Only estimate-specific heterogeneity is
redrawn; no residual sampling innovation remains. This rule applies equally to
diagonal, structured, and dense `V`. When the estimate SD is zero, positive
weights cancel and the generator returns the ordinary Gaussian outcome
directly. In fixed-effect conditional fits, selection weights remain distributed
as their Uniform(0,1) prior. No package decomposition, selection RNG, or
selected-normal density generates the certification data.

Recovery is checked against finite-data posterior uncertainty, not just MCMC
error: on unconstrained coordinates (mean, log-odds weight, and log-SDs), the
difference from the generating parameter must be below four posterior SDs
plus four MCSEs. This avoids treating a skewed positive/bounded posterior as
Gaussian on its raw scale. This is a conservative regression sentinel,
not a repeated-simulation coverage or simulation-based calibration claim.
Posterior-width checks reject uninformative, nearly prior-only recovery for
identified parameters. Fixed-effect conditional weights instead have no
generating recovery target: their first and second posterior moments must agree
with `1/2` and `1/3` within four MCSEs plus `1e-6`. The corresponding mean's first
two moments must agree with its independent conjugate Gaussian posterior.
Every free parameter must also have rank-normalized R-hat below 1.01, bulk ESS
above 1,000, tail ESS above 500, and mean MCSE below 4% of posterior SD.
Each chain has 500 adaptation and 1,500 burn-in iterations, followed by 6,000
draws (24,000 for the correlated conditional examples, whose latent effects
mix more slowly). Seeds are fixed independently of observed recovery; the
test runner never retries or extends chains until expectations pass.

For the dense reference, marginal normalizers are an expansion into eight Gaussian
orthant probabilities evaluated by `mvtnorm::TVPACK()`. The conditional
reference uses retained covariance `C = V + study_sd^2 * 11'` and integrated
covariance `estimate_sd^2 * I`. With positive estimate SD, Gaussian conditioning
gives the unselected `C | Y` law, over which the oracle independently integrates
reciprocal row normalizers. With zero estimate SD and positive weights, it
returns the ordinary Gaussian density exactly. For the nondegenerate conditional
integral, Gaussian quadrature orders 25 and 41 must agree. Base R integration
then integrates the mean parameter over the entire real line with its prior.
JAGS posterior first and second moments must agree with this reference within
four MCSEs plus the independently checked integration allowance. Default
selection integration settings and all diagnostic tolerances are unchanged.
The structured marginal fixed-effect example also checks the joint posterior
first and second moments of the mean and selection weight against independent
two-dimensional quadrature, with orders 31 and 51 checked for agreement.

Use these cases before larger Assink scenario fits when changing selection
likelihoods or JAGS kernels. They do not certify arbitrary covariance dimension,
selection functions, singular boundaries, or repeated-simulation calibration;
the existing numerical-kernel and scenario evidence remains necessary.

Focused endpoint evidence is in `test-00-selection-sampling-native.R` (native
Gaussian/orthant references, deterministic and singular support, and PSD
projections), `test-00-selection-estimate-cells.R` (all eight source cells,
both weight rules, and covariance-representation invariance), and
`test-00-selection-completion.R` (positive-acceptance guarantees and observation
weights). `test-00-selection-sampling-evidence.R` checks marginal likelihood,
deleted-error integration, CDFs and moments against independent integrals;
`test-00-selection-sampling-postfit.R` checks covariance-owned source
reconstruction and singular conditional-mean allocation. These small references
supplement the expensive fitted certification; updating its generator does not
constitute a passing certification run.
