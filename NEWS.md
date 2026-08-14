## version 4.1.5 (IN PROGRESS)
### Features
- evaluates eligible AR1, HAR, and CAR random-effect covariance blocks with an
  exact linear-time Markov likelihood and tridiagonal conditional solver.
  The route retains every covariance parameter and the complete coefficient
  factor, supports missing and repeated states plus row-specific scales, and
  falls back to the existing full covariance paths for correlated sampling
  blocks or general random-effect designs.
- uses the exact diagonal coefficient-factor contract to construct transformed
  random-effect designs in linear rather than quadratic time in the number of
  coefficients. The native evaluator verifies every structural zero before
  taking this route and retains the general dense multiplication otherwise.
- caches the exact fixed random-effect basis and sampling-precision
  cross-products for eligible known-group covariance blocks with diagonal
  coefficient factors. Marginal likelihood and conditional LOO reuse this
  cache while nonzero state-specific sampling variances and general covariance
  structures retain the existing complete evaluators.
- caches exact Cholesky factors, log determinants, and precision diagonals for
  fixed sampling-covariance blocks, reusing them in zero-extra-variance bridge
  and LOO states while nonzero states retain their per-state factorization.
- reuses exact Matrix/CHOLMOD symbolic factorizations for structurally sparse
  latent Gaussian systems in known-`V` marginal likelihoods. Numeric values
  are assembled from every active random-effect coordinate at each state,
  with no diagonal bounds, covariance repair, or altered likelihood; dense
  latent systems retain the LAPACK route.
- extends the exact sparse latent solver through independently factored
  correlated sampling-covariance blocks. Cholesky whitening propagates the
  complete structural active set, so nested and crossed random effects retain
  every covariance contribution without materializing a dense latent system.
- applies the same exact diagonal- and block-base Woodbury plans to known-`V`
  LOO conditional densities. The native evaluator computes only `Q(y - mu)`
  and `diag(Q)` from dense or sparse latent solves, with full-covariance
  inversion retained as the numerical and unsupported-structure fallback.
- reuses the exactly centered varying log-likelihood matrix in `add_loo()` for
  relative effective sample sizes and PSIS instead of constructing the same
  matrix twice, reducing transient memory without changing either calculation.
- consumes BayesTools' compact exact bridge factor-state contract so invariant
  random-effect designs, group maps, and known group kernels are validated once
  while every draw-varying coefficient factor and row scale retains full
  support validation. Full covariance-factor inputs are retained and normalized
  once to the same internal plan/state representation; the native evaluator
  therefore consumes one immutable-geometry and dynamic-state contract while
  the generic full-factor and dense interfaces remain available.
- shares exact latent-system, fixed-known-group, and correlated sampling-block
  preparation between marginal-likelihood and conditional evaluators, and
  routes bridge block-MVN and post-fit joint known-`V` densities through the
  same reusable native covariance plan. Independent R Cholesky and Schur
  implementations remain as numerical references and invalid-covariance
  diagnostics.
- compiles the active latent coordinates of each observation and assembles
  exact diagonal-base Woodbury systems only over their structural nonzeros,
  substantially reducing repeated work for nested and crossed random effects
  without dropping latent coordinates or covariance parameters.
- evaluates exact low-rank random-effect updates that connect otherwise
  independent known sampling-covariance blocks by factoring those base blocks
  separately and applying block-whitened Woodbury algebra. The reusable plan
  retains exact dense-Cholesky fallback for singular base blocks and avoids a
  redundant global eigendecomposition.
- routes covariance factors whose latent width is not below the observation
  dimension through exact BLAS root products and observation-space Cholesky,
  avoiding entry-by-entry `ZGZ'` reconstruction while retaining the dense
  fallback for covariance sources without a root contract.
- stabilizes exact Woodbury quadratic forms when cancellation can dominate the
  usual subtraction, switching to the algebraically identical
  residual-plus-penalty evaluation without clamping or covariance repair.
- speeds exact selected-normal plotting without reducing posterior draws or
  plot grids: `zplot()` reuses invariant affine-density terms and evaluates
  fitted and extrapolated curves together when their predictive target is
  shared; `funnel()` and selection-model `regplot()` cache invariant step-bin
  CDF plans; and `regplot()` batches prediction-specific standard errors and
  reuses exact neighboring quantile roots with the original global fallback.
- converts numerically symmetric known sampling covariance matrices once to an
  exactly symmetric representation at input validation, so fitting,
  prediction, LOO, and bridge sampling use the same covariance values.
- speeds exact marginal-likelihood and LOO post-processing for multivariate
  models: `add_marglik()` inherits or overrides the fitted parallel/core
  settings and uses the nodes-only BayesTools bridge context without changing
  the marginalized parameterization, while diagonal known-`V` estimate-level
  log likelihoods are vectorized and invariant block factors are cached.
- exactly marginalizes sampled Gaussian location random effects in
  `brma.mv()` bridge targets as draw-specific `ZGZ'`, retaining every SD,
  allocation, correlation parameter, and prior. Validated covariance factors
  and reusable native covariance plans compile invariant geometry checks,
  retain per-state covariance and scale validation, avoid materializing
  structural zeros, pass only draw-varying factor state after setup, and use
  exact low-rank likelihood algebra where applicable. Eligible correlated
  known-`V` blocks additionally reuse an exact fixed-covariance eigensystem and
  evaluate low-rank updates with determinant/Woodbury identities; all other
  states retain the dense Cholesky path. Non-Gaussian selection likelihoods
  retain the joint latent target. General correlated known-`V` LOO blocks reuse
  the same plan without changing their Schur target.
- exactly integrates the standardized effects fitted by the legacy `cluster`
  argument for Gaussian marginal likelihoods, retaining heterogeneity, its
  allocation, row-specific scale regression, weights, and all associated
  priors through diagonal-plus-cluster-rank-one likelihood algebra.
- validates and resolves prepared known-`V` metadata once during bridge setup
  instead of repeating invariant representation checks for every bridge state.
- compiles scalar marginalized random-effect variance plans used by known-`V`
  bridge likelihoods, with the complete generic evaluator retained for
  row-indexed and other non-scalar SD sources.
- exposes bridge `repetitions`, `method`, `maxiter`, and `silent` controls from
  `add_marglik()`. Parallel settings default to the fitted model settings, and
  effective sample sizes are always obtained from the fitted chains.
- lets `hypothesis()` use unquoted formula interaction references such as
  `alloc:ablat[random]`, evaluates factor-level qCMDE/IWMDE point hypotheses
  on the exact displayed coefficient scale, and prints public hypothesis
  labels without redundant row names after internal parameter resolution.
- adds `as.data.frame()` and `data.frame()` summary coercion for individual
  `brma_samples` results. Multi-component results return a long data frame with
  component and parameter identifiers; `as.data.frame(format = "list")`
  retains the named, possibly nested table structure. Credible- and prediction-
  interval columns use syntactic names such as `CI_0.025` and `PI_0.025` while
  printed summaries retain their compact quantile labels.
- keeps known-`V` backend routing in summary metadata without printing this
  implementation detail in ordinary `brma.mv()` summaries.
- requests BayesTools' structured component labels for random-effect rows, so
  both single- and multiple-component models print labels such as
  `study: sd(intercept)` without downstream label rewriting.
- extends posterior `plot()` and `lines()` transformations: `EXP` now applies
  to individual log-scale ratio meta-regression coefficients and to
  scale-regression coefficients as multiplicative changes in heterogeneity,
  while `LOG` applies to positive heterogeneity intercepts; KDE, qCMDE, and
  IWMDE curves use the exact change-of-variable Jacobian.
- warns when `hypothesis()` uses the full product-space ensemble for a
  parameter with both null and alternative components because `conditional`
  was omitted; explicitly set `conditional = FALSE` to retain the full-ensemble
  test without the warning, or `conditional = TRUE` to test only active models.
- supports atom-free pairwise point contrasts across levels of one fitted model,
  including symbolic equalities and KDE, qCMDE, or IWMDE inference for fitted
  and marginal-means hypotheses, and accepts equivalent constant-left scalar
  relations such as `0 > theta` and `0 = theta`
- re-exports `loo::loo_model_weights()` and adds a `brma` method that computes
  stacking or pseudo-BMA weights directly from compatible stored LOO results,
  including models fitted with different numbers of posterior draws.
- adds `metafor::forest()` support for `brma` objects and `as_metafor_forest()` for preparing RoBMA forest-plot data.
- adds `brma.mv()` for normal-likelihood meta-analysis with known sampling
  covariance matrices, including latent, whitened, and block-MVN known-`V`
  backends and BayesTools random-effect formula integration.
- adds metafor-style `R`/`Rscale` support to `brma.mv()` for known
  random-effect group covariance in random-intercept blocks, including
  supported one-to-one marginalized known-`R` blocks under known-`V` normal
  models.
- adds opt-in likelihood-aware qCMDE/IWMDE posterior density and point-ordinate
  estimation for supported `plot()`, `hypothesis()`/`bf_hypothesis()`,
  `marginal_means()`, and `hypothesis.marginal_means()` workflows via
  `density_method`.
- adds fixed, user-controlled simple-random posterior-row samples for
  qCMDE/IWMDE density curves and point ordinates through the shared
  `density_control$samples` setting. Density-control lists reject unknown
  settings. Compact public diagnostics are available via
  `density_diagnostics()` and method-neutral `RoBMA_density_diagnostics` /
  `RoBMA_density_ordinate_error` classes. Diagnostics quantify uncertainty but
  never adaptively change the estimator's posterior-row sample.
- adds RoBMA-owned schema/version provenance to precomputed qCMDE/IWMDE density
  and ordinate attributes.
- adds `log_lik()` for pointwise posterior log-likelihood draws used by LOO and
  WAIC.
- adds semantic `component = "random"` parameter summaries, plots, MCMC
  diagnostics, and coherent interval/directional hypotheses for `brma.mv()`.
- lets semantic multivariate random-effect plots use qCMDE for direct sampled
  quantities and two-component allocation fractions, and draws exact
  Beta marginal priors for Dirichlet allocation fractions.
- adds `include_auxiliary` to all fitted-model `as_draws*()` methods; the
  default hides backend-only variables and `TRUE` exposes raw backend draws.

### Testing and development
- separates cached-fit refresh from the standard regression profile, keeps the
  cached standard suite under a hard 15-minute budget, and restores caches in
  CI using fitting-source and backend fingerprints.
- splits numerical certification into independently runnable one-hour cases
  without a total certification limit.
- adds maintainer analysis scenarios with automatic fit caching and tracked
  text and visual baselines, plus representative prior and multivariate visual
  regression coverage.
- lets `scenario_text()` automatically print visible returned values, removing
  the need for explicit `print()` calls around summaries and tables.
- allows interactive `scenario_plot()` calls outside a scenario runner to draw
  directly without requiring a vdiffr snapshot context.

### Breaking changes
- requires BayesTools 0.3.1.33 and R 4.3.0 for the multivariate random-effect
  backend, point-prior monitoring, exact zero-dimensional marginal likelihoods,
  scalable diagonal marginal variances, versioned fitted-formula identities,
  deterministic draw geometry, metadata-only parameter catalogs, hypothesis
  ASTs, structural prior-ordinate classification, and exact induced formula-
  coefficient prior densities, structured random-effect component labels, and
  consistent base posterior-overlay spike scaling.
- requires loo 2.10.0 internally while preserving RoBMA's released numeric
  `compare.loo` matrix and printing contract.
- removes transitional pre-release known-`V`, dense random-correlation, and
  BayesTools capability fallback paths; multivariate fits now use only the
  current canonical backend metadata and typed conditions.
- removes the released `marginal_means.brma(normal_approximation)` argument;
  marginal-means Bayes factors now use posterior ordinates from KDE/qCMDE/IWMDE
  routes instead of the removed normal-approximation switch.
- simplifies the marginal-means density contract: plot and hypothesis
  methods inherit the stored density method unless explicitly overridden,
  `bf = FALSE` skips qCMDE/IWMDE ordinate precomputation, and single-model
  objects reject explicit conditional density precomputation.
- removes the released `logLik.brma()` method. Use `log_lik()` for pointwise
  posterior log-likelihood draws; `AIC()` and `BIC()` now fail explicitly.
- removes the unreleased logical aggregate prediction mode. Explicit prediction
  rows now each represent one new true effect; use `pooled_effect()` for direct
  fitted-design aggregation.
- corrects the `BMA.glmm()` factor-contrast default from treatment coding to
  mean-difference coding, matching the other BMA/RoBMA constructors. This can
  change moderator coefficient interpretation for calls that did not specify
  `set_contrast_factor_predictors` explicitly.
- reports `cooks.distance()` as the unscaled squared posterior Mahalanobis
  distance, following metafor's chi-square-based meta-analytic convention.
  Values were previously divided by the fixed-effect model rank and therefore
  increase by that rank when it exceeds one.

### Fixes
- treats `xlim` as a displayed-scale range for transformed posterior and
  marginal-means plots, so `EXP` prior curves extend toward zero instead of
  starting at an exponentiated fitted-scale limit.
- lets `ylim2` and `ylab2` control the secondary point-mass probability axis,
  keeps its mapping fixed across prior and posterior `lines()` overlays from
  multiple objects, and warns when an overlay falls outside the active limits
- accelerates native funnel-contour and regression-plot interval inversion for
  continuous posterior mixtures with a shared machine-relative Brent solve,
  while retaining exact generalized-inverse bisection whenever the mixture
  contains an atomic component.
- makes `regplot(mod = "vi")` and `regplot(mod = "sei")` preserve the
  prediction-grid precision, derive the matching standard error or variance,
  and use pointwise precision in sampling intervals
- lets `regplot()` use moderators that appear only in a scale formula, showing
  a constant location effect with heterogeneity-dependent prediction and
  sampling intervals
- restores named factor-level resolution in fitted-model hypotheses through
  BayesTools' semantic factor catalog while retaining a single RoBMA-owned
  whole-term target for product-space inference, and recognizes native sampled
  or structural random-effect quantities. It also permits mixed
  point/region marginal-means hypotheses when all events use the same
  single-model averaged posterior and avoids density-ordinate work for pure
  region hypotheses.
- adds marginal-means `lines()` overlays and routes qCMDE/IWMDE coefficient
  plots through the exact fitted-to-original structural transform.
- makes fixed-budget qCMDE/IWMDE precision warnings explain how to request more
  evaluated posterior samples without presenting the diagnostic target as an
  adaptive stopping rule.
- permits mixed-case scenario artifact names such as `mu_BF_comparison`,
  matching the names already used by the maintainer scenarios.
- resets base-graphics overlay mode around interactive `scenario_plot()`
  evaluations so consecutive scenario figures do not draw over one another.
- extends existing fits through the supported BayesTools continuation contract,
  preserving each JAGS chain's stored random-number generator state;
  `update.brma(seed = ...)` now rejects attempts to reseed a continuation.
- lets named `NULL` values explicitly clear nullable convergence/autofit
  thresholds and targeted monitors during `update()`, while omitted controls
  continue to inherit their fitted values.
- rejects explicitly supplied moderator, scale-regression, and multilevel
  heterogeneity-allocation priors when the corresponding model component is
  absent, instead of silently discarding them.
- checks product-space model indicators by default through BayesTools'
  label-invariant state-occupancy diagnostics; `set_autofit_control()` and
  `set_convergence_checks()` retain an explicit opt-out and parameter routing,
  and targeted monitor selections no longer suppress eligible indicators.
- reports an actionable error when additional models are passed to
  `loo_weights()`, clarifying the distinction between single-model PSIS
  importance weights and the new `loo_model_weights()` method across models.
- evaluates qCMDE/IWMDE density-curve relative MCSE, effective sample size,
  and importance-weight concentration over the empirical 5%--95% bulk, records
  the 5% and 95% tail checkpoints, and retains a whole-curve absolute MCSE
  safeguard relative to the density peak. Extreme display endpoints remain
  available in raw diagnostics but no longer reject an otherwise reliable
  curve; requested Bayes-factor ordinates keep their strict local gates.
- accepts symbol-valued primitive hypothesis AST nodes when attaching
  qCMDE/IWMDE warning records from the current BayesTools contract.
- shortens qCMDE/IWMDE relative-MCSE and effective-sample-size warnings to
  report only the observed diagnostic value instead of repeating policy
  thresholds.
- supports fully fixed point-prior models across `brma`, publication-bias,
  GLMM, multilevel, and `brma.mv()` classes; fixed parameters remain available
  in posterior summaries, and zero-dimensional marginal likelihoods are exact.
- treats constant log-likelihood draws from fully fixed models as exact LOO
  importance ratios instead of reporting spurious infinite Pareto diagnostics.
- delegates fitted formula identities, parameter resolution, hypothesis parsing,
  and original-scale coefficient prior ordinates to their versioned BayesTools
  contracts, including grouped factor terms and overlapping component aliases.
- materializes deterministic and latent-only fitted draws from persisted chain
  geometry and structural values while keeping private backend anchors out of
  public posterior output.
- resolves formula-prior intercept precedence and supports omitting a complete
  formula-prior branch while continuing to reject ambiguous partial omission.
- binds cached marginal likelihoods to the fitted target, canonicalizes outcome
  identities, and preserves structured interpretation and diagnostic failures.
- aligns weighted heterogeneity, known-`R` scale prediction, rank-deficient VIF,
  influence, deleted-tau, and leverage diagnostics with their fitted posterior
  and model-averaged estimands.
- defines exact forest, funnel, and regplot prediction targets and intervals,
  including new-group random-effect semantics, pointwise predictive contours,
  current `metafor` compatibility, and consistent base/ggplot rendering.
- honors explicit zplot axes and cutpoints while preserving histogram-method
  delegation.
- uses one global nested simple-random sample for qCMDE/IWMDE and
  model-averaged funnel calculations, preserving empirical product-state
  weights without enumerating or forcibly retaining rare states.
- fits IWMDE proposals to the full eligible posterior population and samples
  only the expensive contribution evaluations, preserving the sampling target
  while avoiding repeated proposal fits. Boundary-risk Chen fallbacks are
  preflighted before constructing the full nuisance matrix.
- evaluates nested qCMDE integration grids incrementally and stops after the
  first certified refinement pair, reusing all previously evaluated nodes.
- aligns same-data `brma.mv()` response prediction with `brma()` semantics by
  using fixed means and marginal known-`V` plus random-effect covariance instead
  of conditioning on fitted random effects.
- requires `V_new` for explicit known-`V` response predictions, rejects unseen
  known-`R` levels without an `R_new` interface, and makes forest prediction
  intervals target the pooled average design when `newdata` is omitted or
  exactly one explicit new true effect when it is supplied.
- supports singular positive-semidefinite `V` when compiled priors or
  marginalized random effects structurally regularize every null direction;
  invalid negative eigenvalues now fail input validation instead of being
  projected to the PSD boundary. Exact rank-one blocks use one latent sampling
  factor so correlations of one and separate tiny diagonal variance remain
  structurally distinct.
- retains diagonal and block-list known-`V` inputs in compact canonical form
  while preserving numerical parity across fitting, prediction, likelihood,
  diagnostics, hashing, and bridge-sampling consumers.
- averages posterior coefficient covariance when computing VIF/GVIF instead of
  plugging posterior mean heterogeneity into a single covariance matrix.
- validates cached LOO/WAIC objects against the current data and target
  fingerprint before reuse.
- makes LOO/WAIC outcome fingerprints explicitly versioned and independent of
  the R serialization writer version, and refreshes affected precomputed fits.
- evaluates `pooled_heterogeneity()` at the average expanded scale and random
  design and adds the corresponding prediction-interval columns to
  `pooled_effect()`; observed-design summaries, funnel plots, radial plots, and
  influence diagnostics retain their within-draw RMS heterogeneity targets.
- hardens known-`V` estimate log-likelihood sums, diagnostics, funnel variance
  extraction, qCMDE/IWMDE random-SD handling, and `brma.mv()` target metadata
  regressions.
- separates qCMDE display, evaluation, and integration grids so plot ranges and
  boundary nulls no longer distort row-normalizer grids.
- makes qCMDE BF diagnostics use final-vs-validation ordinate movement instead
  of same-grid normalization mass.
- fails qCMDE/IWMDE estimation when any selected row lacks a valid normalizer,
  proposal density, or finite contribution instead of silently dropping the
  row or replacing its contribution by zero.
- freezes finite baseline row eligibility in qCMDE/IWMDE plans, cache keys, and
  diagnostics so denominator rows and estimator rows are reported separately.
- conditions GLMM qCMDE/IWMDE ordinates on sampled estimate-level effects and
  baserates or log-rates, then averages over posterior rows instead of using a
  prior-scale nuisance grid inside the density estimator.
- replaces ordinary estimate-unit GLMM fixed-grid integration with native joint
  adaptive quadrature and a convergence-checked prior-CDF fallback for
  supported truncated nuisance priors.
- rejects cluster-unit GLMM log-likelihood, LOO, and WAIC until certified nested
  adaptive quadrature is available; estimate-unit GLMM diagnostics remain
  supported.
- rejects GLMM IWMDE density estimation after all-row bridge and density-policy
  probes missed the approved certification tolerances; qCMDE remains supported
  for descriptive density curves and point-null Bayes factors.
- evaluates qCMDE/IWMDE boundary point nulls with explicit requested/evaluation
  value provenance for one-sided Savage-Dickey Bayes factors.
- classifies primitive and induced qCMDE/IWMDE prior ordinates before estimation
  and warns for zero, infinite, discrete, undefined, or structurally unknown
  target densities without reconstructing them from posterior draws.
- reuses marginal-means qCMDE/IWMDE ordinates only when provenance proves the
  method, settings, target, source object, and value are compatible.
- increases the `metafor::rma.glmm()` optimizer budget in the GLMM parity
  vignette so full vignette builds are reproducible without changing results.
- rejects unused arguments to multivariate pooled-heterogeneity methods instead
  of silently ignoring misspelled controls.
- accepts the same four-cell binomial count representation in explicit
  predictions as in model fitting and rejects inconsistent duplicate totals.
- returns raw binomial and Poisson response predictions as count samples without
  effect-size transformation metadata.
- uses within-draw RMS heterogeneity consistently for model-averaged funnel
  contours.
- normalizes selected-normal funnel and regression-plot CDF tails jointly,
  preventing valid low-standard-error contours from failing due to separately
  rounded normalization masses.
- certifies selected-normal cluster integrals by successive 7-, 15-, and
  31-node quadrature rules and fails explicitly when the requested accuracy is
  not established, replacing heuristic curvature fallbacks.
- preserves moderator-specific forest labels unless `mlab` is explicitly
  overridden.
- supplies `metafor` forest shade/dist styles with a deterministic density from
  the actual posterior predictive draws instead of a Normal reconstruction from
  interval endpoints.
- makes certification workers fail when required numerical evidence is missing,
  skipped, duplicated, or executes without a passing expectation.
- applies known-R row-variance multipliers consistently to scale predictions,
  random-effect draws, BLUPs, and pooled heterogeneity.
- evaluates random-slope `terms.scale` predictions from each explicit new row's
  fitted-scale random design instead of reusing the fitted design.
- computes deleted-row location-scale influence heterogeneity with the package's
  within-draw RMS target.
- recovers sparse valid binomial and Poisson GLMM likelihood columns with a
  convergence-checked prior-CDF quadrature fallback when AGHQ cannot certify
  them.
- evaluates boundary-concentrated and endpoint-truncated beta baserate priors
  with boundary-matched quadrature, preventing valid sparse GLMM fits from
  failing in `log_lik()`, LOO, and WAIC.
- samples qCMDE/IWMDE rows uniformly without replacement from all eligible
  continuous active posterior rows, preserving the fitted product-state mass
  without enumerating nuisance model states or forcing rare-state inclusion.
- separates finite-population row-sampling uncertainty from selected-row MCMC
  MCSE and effective sample size and requires per-chain batch coverage. Mixed
  point/continuous estimates may subsample the active rows: they report the
  selected active-branch MCSE and the complete indicator-chain active-mass
  MCSE separately and use their worst-correlation delta upper bound on the
  unconditional density scale. A full active-row census still batches the
  complete conditioned chain, including active-mass covariance, and is
  required for strict Bayes-factor precision grading.
- removes structurally inactive random-effect coordinates from marginal-
  likelihood replay when every corresponding standard deviation is fixed at
  zero, restoring exact zero-dimensional results.
- prunes individually fixed-zero random-formula blocks and independent latent
  coordinates from bridge proposals while reconstructing their structural zero
  values in the target.
- handles every finite constant log-likelihood column with exact uniform LOO
  importance ratios and zero Pareto-k diagnostics while retaining ordinary
  PSIS diagnostics for varying columns.
- restricts binomial GLMM baserate priors to point or beta families and Poisson
  GLMM log-rate priors to point or normal families, including their supported
  truncations, so every accepted nuisance prior has a certified post-fit
  likelihood route.
- rejects ambiguous marginal-means hypothesis aliases instead of silently
  resolving a canonical name to a different moderator coefficient.
- evaluates certified `exp(affine)` formula-coefficient hypotheses under KDE
  only for structurally atom-free, unconditional scalar targets. Point
  equalities are evaluated on the inverse log/affine scale, where the Jacobian
  cancels; nonpositive nulls, compound point expressions, and nonlinear
  qCMDE/IWMDE routes fail clearly.
- rejects random-parameter point hypotheses when the declared transformed
  prior contains an atom, while retaining coherent region and directional
  hypotheses.
- rewrites vector hypothesis aliases independently while limiting point-null
  syntax to direct parameters and levels so product-space atoms and
  conditioning metadata cannot be discarded.
- limits hypothesis discovery to supported estimators and syntax, and retains
  only qCMDE/IWMDE ordinates and diagnostics matching the current request.
- rejects unnamed, duplicate, nonnumeric, or non-finite qCMDE/IWMDE linear
  target weights before removing exact zero coefficients.
- fails qCMDE/IWMDE construction when a selected row lacks a finite baseline,
  certified normalizer, or normalized proposal density; mathematically exact
  zero ordinates remain valid `-Inf` log-density contributions.
- computes marginal means through the public BayesTools BF-free route when
  inclusion Bayes factors are disabled, retaining supported no-intercept and
  structurally fixed marginal cells.
- preserves the released secondary `marginal_means` result class and the
  released positional `type` argument of `predict.brma()` while adding the
  named `V_new` covariance argument.
- returns genuine upstream `bridge` objects from `bridge_sampler()` for
  bridge-backed fits; exact zero-dimensional fits retain exact `logml()`, Bayes
  factors, and posterior probabilities and report clearly that no bridge object
  exists.
- reconstructs sampled CS/HAR random-effect correlations in every public
  `as_draws*()` format while keeping compact backend coordinates private.
- detaches persisted formula and terms environments after their variables have
  been materialized, preventing fitted objects and vignette caches from
  retaining unrelated caller workspaces.
- uses nested qCMDE normalization grids and removes superseded row-drop recovery
  machinery, reducing repeated work and simplifying the fail-closed estimator.
- treats eigensolver null-space artifacts inside the certified covariance
  backward-error envelope as zero only in spectral factor representations,
  preserving valid low-rank known-`V` inputs without modifying submitted
  covariances.
- evaluates diagonal-plus-rank-one Gaussian quadratic forms through a stable
  positive-sum identity, avoiding cancellation near singular covariance
  boundaries.
- keeps parser-only prediction placeholders out of formula data and retains
  coherent `vi`/`sei` views from `V_new`, so missing formula predictors fail
  rather than being silently evaluated at zero.
- derives ordinary multilevel `ranef()` cluster and estimate components from
  one joint Gaussian BLUP solve instead of mixing analytic totals with sampled
  latent components.
- rejects moderator-dependent PET/PEESE curves and outcome-mode funnels with
  location or scale predictors, where a single unconditional curve is not a
  fitted estimand.
- marks GLMM outcome-mode funnel contours as descriptive normal effect-size
  approximations in warnings, plot labels, and returned metadata instead of
  implying coverage under the fitted discrete likelihood.

### Documentation
- documents the qCMDE/IWMDE estimating equations, Savage-Dickey nesting
  requirement, method selection, tuning controls, reliability diagnostics,
  unsupported targets, sensitivity workflow, and normal-approximation limits.
- clarifies the separate posterior-row sampling and selected-row MCMC errors;
  prior-ordinate and IWMDE proposal-weight estimation uncertainty remain
  outside `BF_error`.

## version 4.0.0
### Breaking changes
- rewrites the package around the unified `brma` class hierarchy. Single-model fits now use `brma()`, `brma.glmm()`, `bselmodel()`, `bPET()`, and `bPEESE()`; model-averaged fits use `BMA()`, `BMA.glmm()`, and `RoBMA()`.
- removes the legacy `RoBMA.reg()`, `NoBMA()`, `NoBMA.reg()`, `BiBMA()`, and `BiBMA.reg()` constructors. Use `mods`, `scale`, and `cluster` in the new constructors, `BMA()` for no-bias normal-likelihood model averaging, and `BMA.glmm()` for GLMM model averaging.
- replaces old input aliases such as `d`, `r`, `logOR`, `OR`, `z`, `y`, `se`, `v`, `n`, `study_names`, `study_ids`, `weight`, and `transformation` with `yi`, `vi`/`sei`, `ni`, `slab`, `cluster`, `weights`, `measure`, `output_measure`, and `transform`.
- removes legacy helper APIs including `combine_data()`, `check_setup()`, `extract_posterior()`, `marginal_summary()`, `marginal_plot()`, `plot_models()`, `adjusted_effect()`, `as_zcurve()`, and the old z-curve plotting methods.
- normal-likelihood fitting functions now require an explicit `measure` for fitted models. Use `measure = "GEN"` for generic effect sizes without a known unit-information scale.
- `update()` for `brma` objects now focuses on extending MCMC samples, updating labels, and refreshing cached quantities, not changing model structure.
- `set_convergence_checks()` no longer accepts the old `remove_failed` and `balance_probability` arguments.

### Features
- adds `brma()` / `brma.norm()` for single normal-likelihood Bayesian meta-analysis, including random-effects, meta-regression, multilevel, and location-scale models.
- adds `brma.glmm()` for binomial-normal and Poisson-normal GLMM meta-analysis from raw two-arm counts (`measure = "OR"` and `"IRR"`).
- adds single-model publication-bias constructors `bselmodel()`, `bPET()`, and `bPEESE()`.
- adds `BMA()` / `BMA.norm()` for Bayesian model averaging without publication-bias adjustment.
- adds `BMA.glmm()` for Bayesian model averaging of GLMM meta-analyses without publication-bias adjustment.
- rewrites `RoBMA()` as a product-space model-averaged ensemble over effect, heterogeneity, moderator, scale, and publication-bias components.
- adds formula/data-frame input handling for effect sizes, moderators, scale predictors, clusters, labels, subsets, likelihood weights, and raw GLMM counts.
- adds default prior construction from standardized effect-size measures, estimated or manually supplied unit-information standard deviations, and informed empirical priors.
- adds `prior_weightfunction()`, `wf_cumulative()`, `wf_fixed()`, and `wf_independent()` for BayesTools-backed selection-weightfunction priors.
- adds `prior_PET()`, `prior_PEESE()`, `prior_none()`, `prior_factor()`, `prior_informed()`, and BayesTools contrast helpers as package-level prior utilities.
- adds `posterior` package interfaces via `as_draws()`, `as_draws_array()`, `as_draws_df()`, `as_draws_list()`, `as_draws_matrix()`, and `as_draws_rvars()` for fitted models and `brma_samples`.
- adds the `brma_samples` posterior-sample class with print, summary, matrix, and `posterior` conversion methods.
- adds `predict.brma()` for posterior predictions of fixed terms, cluster effects, latent true effects, observed responses, and scale terms, with `newdata`, `conditional`, `bias_adjusted`, `output_measure`, and `transform` support.
- adds convenience wrappers `fitted()`, `pooled_effect()`, `pooled_heterogeneity()`, `blup()`, `true_effects()`, and `ranef()` for `brma` objects.
- adds model-comparison helpers `add_loo()`, `loo()`, `loo_compare()`, `loo_weights()`, `check_loo()`, `add_waic()`, and `waic()` using the `loo` package.
- adds bridge-sampling marginal likelihood support for single-model `brma` fits via `add_marglik()`, `bridge_sampler()`, `logml()`, `bf()`, `bayes_factor()`, and `post_prob()`.
- adds residual and influence diagnostics: `residuals()`, `rstandard()`, `rstudent()` / `LOO-PIT`, `hatvalues()`, `influence()`, `dfbetas()`, `dffits()`, `cooks.distance()`, `covratio()`, and `vif()`.
- adds plotting methods for `brma` objects: posterior/prior plots, `funnel()`, `regplot()`, `qqnorm()`, `radial()` / `galbraith()`, MCMC diagnostic plots, weightfunction plots, and PET-PEESE plots.
- adds `marginal_means()` with summary and plotting methods for moderator models.
- adds `summary_models()` for marginal and individual model-weight summaries of product-space `RoBMA`, `BMA`, and `BMA.glmm` objects.
- adds `interpret()` for concise textual interpretation of fitted `brma` and model-averaged objects.
- renames the zplot diagnostic API to `as_zplot()` and adds the direct plotting wrapper `zplot()`, with `plot()`, `hist()`, `lines()`, `summary()`, and print methods for zplot objects.
- adds `RoBMA.options()` and `RoBMA.get_option()` package options for defaults such as core count, automatic LOO/WAIC/marginal-likelihood computation, prior scaling defaults, and selection-bias defaults.

### Changes
- renames the multilevel clustering argument to `cluster`.
- renames study labels to `slab`, matching `metafor` naming.
- renames likelihood weights to `weights` for supported constructors and applies them consistently to posterior fitting, log-likelihoods, LOO, WAIC, and diagnostics; `brma.mv()` currently rejects likelihood weights.
- uses `measure`, `output_measure`, and `transform` for effect-size scale handling. Supported conversions include `SMD`, `COR`, `ZCOR`, and `OR`; `transform = "EXP"` exponentiates log ratio measures for display.
- standardizes continuous predictors by default and transforms reported coefficients back to the original scale unless standardized coefficients are requested.
- uses treatment contrasts by default for single-model constructors and mean-difference contrasts by default for model-averaged constructors.
- changes `predict.brma()` default to `type = "terms"`. GLMM `type = "response"` predictions return continuity-corrected effect-size estimators by default via `as_measure = TRUE`.
- separates output `unit` from `conditioning_depth` for residuals, fitted values, LOO, WAIC, and related diagnostics.
- supports estimate-level LOO/WAIC targets for `brma.mv()` and estimate-/cluster-level LOO/WAIC targets for multilevel `brma()` models, with target metadata to prevent invalid comparisons.
- keeps bridge-sampling marginal likelihoods for single-model `brma` objects; product-space `RoBMA`, `BMA`, and `BMA.glmm` objects rely on product-space only.
- routes selection-weightfunction priors through the BayesTools selection backend and selected-normal kernel, removing legacy weighted-normal mapping paths.
- uses `bias_indicator` and branch-aware selected-normal contexts for RoBMA publication-bias mixtures instead of inferring selection branches from `omega`.
- increases zplot default posterior thinning controls to `10000` samples and accepts `Inf` where full posterior evaluation is requested.
- adds `max_samples` controls to expensive funnel, regplot, and zplot summaries.
- updates the package startup message to point users to `vignette("v00-introduction", package = "RoBMA")`.
- requires BayesTools 0.3.0 for forward API and selection-backend support.
- adds `bridgesampling`, `loo`, `MASS`, and `parallel` as imports and `posterior` as a suggested package.

### Fixes
- fixes loading and runtime checks for the RoBMA JAGS module and native R routines.

### Performance and internals
- moves fitting to JAGS product-space models with mixture-prior indicators for model averaging.
- replaces legacy weighted-normal and multivariate-normal native code with selected-normal kernels shared by JAGS and R-native calls.
- adds native selected-normal routines for log likelihoods, normalizers, CDFs, moments, RNG, weighted summaries, funnel contours, regplot intervals, and zplot densities/threshold summaries.
- adds native estimate-unit GLMM marginal log-likelihood helpers for binomial and Poisson models.
- caches selected-normal normalizers and uses telescoping selection probabilities with log-space fallbacks for better numerical stability.
- relocates selected-normal C++ code to `src/selnorm/` and updates `Makevars*`, native registration, cleanup rules, and JAGS distribution registration.
- removes unused native matrix/LAPACK helper sources and older source-level transformation helpers.

### Documentation and tests
- reorganizes vignettes into numbered workflows covering introduction, prior distributions, baseline Bayesian meta-analysis, feature coverage, metafor parity, model averaging, RoBMA, multilevel models, medicine examples, and zplot diagnostics.
- regenerates roxygen documentation for the new constructors, priors, predictions, summaries, diagnostics, plots, model-comparison methods, and datasets.
- refreshes the README and pkgdown site for the 4.0.0 API.
- adds cached model fits under numbered vignette/model directories.
- refactors tests into ordered input, fitting, prediction, plotting, diagnostics, model-comparison, selected-normal kernel, and vignette-cache coverage.
- adds regression tests for selected-normal telescope probabilities, native/R fallback parity, posterior-row alignment, GLMM response conversion, LOO/WAIC targets, bridge sampling, and visual outputs.

## version 3.6.1
### Features
- `Explanation` vignette that helps navigate users through the vignettes
- two vignettes demonstrating robust Bayesian meta-analysis and meta-regressions
- `summary()` function now provides publication bias model type summary (`type = "models"`) for models fitted using `algorithm = "ss"`
- improves control over zplot diagnostics (i.e., specifying col, border, etc for the individual elements)

## version 3.6
### Features
- `funnel()` plot to visualize residuals vs the expected sampling distribution for `RoBMA()` and `RoBMA.reg()` models when using the `algorithm = "ss"`
- `residuals()` method for `RoBMA()` and `RoBMA.reg()` models when using the `algorithm = "ss"`
- `as_zplot()` function to transform meta-analytic models into a zplot object, only available for `RoBMA()` and `RoBMA.reg()` fitted using the `algorithm = "ss"`
- `plot()`, `summary()`, and `print()` functions for the `as_zplot` objects

## version 3.5.1
### Features
- `summary()` function now supports a `standardized_coefficients` argument to report either standardized (default) or raw meta-regression coefficients
- `extract()` function to extract the posterior samples of the model parameters
- `true_effects()` function to summarize the true effect size estimates of `RoBMA()` and `RoBMA.reg()` models when using the `algorithm = "ss"`
- `predict()` method for `RoBMA()` and `RoBMA.reg()` models when using the `algorithm = "ss"`

### Fixes
- fitting a meta-regression using predictors with missing values result in a clear error message

### Changes
- improving the speed of unit tests

## version 3.5
### Features
- approximate and computationally feasibly 3lvl selection models via the `RoBMA()` and `RoBMA.reg()` functions with the `cluster` argument when using `algorithm = "ss"`
- 3lvl binomial-normal models for binary data via the `BiBMA` and `BiBMA.reg` functions with the `cluster` argument when using `algorithm = "ss"`
- `pooled_effect()` function to compute the pooled effect size from the `RoBMA.reg`, `NoBMA.reg`, and `BiBMA.reg` models
- `adjusted_effect()` function to compute the adjusted effect size from the `RoBMA.reg`, `NoBMA.reg`, and `BiBMA.reg` models
- enables `summary_heterogeneity()` for BiBMA models

### Fixes
- passing and checks of the `cluster` and `study_labels` arguments
- PEESE prior distribution now scale as 1/scale instead of 1/scale^2 with the `rescale_priors` argument  
- the conditional prediction interval based on `summary_heterogeneity()` is now conditional on the presence of the effect
- additional minor prior handling fixes (i.e., missing marginal estimates when only alternative prior distributions were specified etc)
- diagnostics with mixture baseline priors when using `algorithm = "ss"`
- `summary_heterogeneity()` with only a single study does not produce relative heterogeneity instead of crashing

## version 3.4
### Features
- adding binomial-normal meta-regression models for binary data via the `BiBMA.reg` function
- the spike and slab algorithm for faster model estimation via the `algorithm = "ss"` argument for BiBMA models
- default prior distributions for all parameters of BiBMA models are now set via the `set_default_binomial_priors()` function

## version 3.3
### Features
- the spike and slab algorithm for faster model estimation via the `algorithm = "ss"` argument (see a new vignette for more details)
- refactoring of the JAGS C++ code of weighted distributions and exporting of the lpdfs into JAGS (maintenance)
- weights_mix JAGS prior distribution to sample a mixture of weight functions directly

### Fixes
- incorrectly omitting models with more than one predictor when computing conditional marginal summary

## version 3.2.1
### Features
- default prior distributions for all parameters are now set via the `set_default_priors()` function
- `rescale_priors` argument allows to conveniently re-scale the prior distributions for the effect, heterogeneity, and bias simultaneously

## version 3.2
### Features
- `summary_heterogeneity()` function to summarize the heterogeneity of the RoBMA models (prediction interval, tau, tau^2, I^2, and H^2)
- `check_RoBMA_convergence()` function to check the convergence of the RoBMA models
- adds informed prior distributions for binary and time-to-event outcomes via BayesTools 0.2.17

### Fixes
- checking and fixing the number of available cores upon loading the package (hopefully fixes some parallelization issues)
- `update()` function re-evaluates convergence checks of individual models (https://github.com/FBartos/RoBMA/issues/34) 
- typos and minor issues in the vignettes


## version 3.1
### Features
- binomial-normal models for binary data via the `BiBMA` function
- `NoBMA` and `NoBMA.reg()` functions as wrappers around `RoBMA` `RoBMA.reg()` functions for simpler specification of publication bias unadjusted Bayesian model-averaged meta-analysis
- adding odds ratios output transformation` 
- extending (instead of a complete refitting) of models via the `update.RoBMA()` function (only non-converged models by default or all by setting `extend_all = TRUE`)

### Fixes
- handling of non-converged models

## version 3.0.1
### Fixes (thanks to Don & Rens)
- compilation issues with Clang (https://github.com/FBartos/RoBMA/issues/28)
- lapack path specifications (https://github.com/FBartos/RoBMA/issues/24)

## version 3.0
### Features
- meta-regression with `RoBMA.reg()` function
- posterior marginal summary and plots for the `RoBMA.reg` models with `summary_marginal()` and `plot_marginal()` functions
- new vignette on hierarchical Bayesian model-averaged meta-analysis
- new vignette on robust Bayesian model-averaged meta-regression
- adding vignette from AMPPS tutorial
- faster implementation of JAGS multivariate normal distribution (based on the BUGS JAGS module)
- incorporating `weight` argument in the `RoBMA` and `combine_data` functions in order to pass `custom` likelihood weights
- ability to use inverse square weights in the weighted meta-analysis by setting a `weighted_type = "inverse_sqrt"` argument 

### Changes
- reworked interface for the hierarchical models. Prior distributions are now specified via the `priors_hierarchical` and `priors_hierarchical_null` arguments instead of `priors_rho` and `priors_rho_null`. The model summary now shows `Hierarchical` component summary.

## version 2.3.2
### Fixes
- suppressing start-up message 
- cleaning up imports

## version 2.3.1
### Fixes
- fixing weighted meta-analysis parameterization 

## version 2.3
### Features
- weighted meta-analysis by specifying `cluster` argument in `RoBMA()` and setting `weighted = TRUE`. The likelihood contribution of estimates from each study is down-weighted proportionally to the number of estimates in that study. Note that this experimental feature is supposed to provide a conservative alternative for estimating RoBMA in cases with multiple estimates from a study where the multivariate option is not computationally feasible.

## version 2.2.3
### Fixes
- updating the Makevars to install with R 4.2 and JAGS 4.3.1

## version 2.2.2
### Fixes
- updating the C++ to compile on M1 Mac

## version 2.2.1
### Changes
- message about the effect size scale of parameter estimates is always shown
- compatibility with BayesTools 0.2.0+

## version 2.2
### Features
- three-level meta-analysis by specifying `cluster` argument in `RoBMA`. However, note that this is (1) an experimental feature and (2) the computational expense of fitting selection models with clustering is extreme. As of now, it is almost impossible to have more than 2-3 estimates clustered within a single study).

## version 2.1.2
### Fixes
- adding Windows ucrt patch (thanks to Tomas Kalibera)
- adding BayesTools version check

## version 2.1.1
### Fixes
- incorrectly formatted citations in vignettes and capitalization

### Features
- adding `informed_prior()` function (from the BayesTools package) that allows specification of various informed prior distributions from the field of medicine and psychology
- adding a vignette reproducing the example of dentine sensitivity with the informed Bayesian model-averaged meta-analysis from Bartoš et al., 2021 ([open-access](https://onlinelibrary.wiley.com/doi/10.1002/sim.9170)),
- further reductions of fitted object size when setting `save = "min"`

## version 2.1
### Fixes
- more informative error message when the JAGS module fails to load
- correcting wrong PEESE transformation for the individual models summaries (issue #12)
- fixing error message for missing conditional PET-PEESE
- fixing incorrect lower bound check for log(OR)

### Features
- adding `interpret()` function (issue #11)
- adding effect size transformation via `output_scale` argument to `plot()` and `plot_models()` functions
- better handling of effect size transformations and scaling - BayesTools style back-end functions with Jacobian transformations

## version 2.0
Please notice that this is a major release that breaks backwards compatibility.

### Changes
 - naming of the arguments specifying prior distributions for the different parameters/components of the models changed (`priors_mu` -> `priors_effect`, `priors_tau` -> `priors_heterogeneity`, and `priors_omega` -> `priors_bias`),
 - prior distributions for specifying weight functions now use a dedicated function (`prior(distribution = "two.sided", parameters = ...)` -> `prior_weightfunction(distribution = "two.sided", parameters = ...)`),
 - new dedicated function for specifying no publication bias adjustment component / no heterogeneity component (`prior_none()`),
 - new dedicated functions for specifying models with the PET and PEESE publication bias adjustments (`prior_PET(distribution = "Cauchy", parameters = ...)` and `prior_PEESE(distribution = "Cauchy", parameters = ...)`),
 - new default prior distribution specification for the publication bias adjustment part of the models (corresponding to the RoBMA-PSMA model from Bartoš et al., 2021 [manuscript](https://doi.org/10.1002/jrsm.1594)),
 - new `model_type` argument allowing to specify different "pre-canned" models (`"PSMA"` = RoBMA-PSMA, `"PP"` = RoBMA-PP, `"2w"` = corresponding to Maier et al., in press , [manuscript](https://doi.org/10.1037/met0000405)),
 - `combine_data` function allows combination of different effect sizes / variability measures into a common effect size measure (also used from within the `RoBMA` function),
 - better and improved automatic fitting procedure now enabled by default (can be turned of with `autofit = FALSE`)
 - prior distributions can be specified on the different scale than the supplied effect sizes (the package fits the model on Fisher's z scale and back transforms the results back to the scale that was used for prior distributions specification, Cohen's d by default, but both of them can be overwritten with the `prior_scale` and `transformation` arguments),
 - new prior distributions, e.g., beta or fixed weight functions,
 - estimates from individual models are now plotted with the `plot_models()` function and the forest plot can be obtained with the `forest()` function,
 - the posterior distribution plots for the individual weights are no able supported, however, the weightfunction and the PET-PEESE publication bias adjustments can be visualized with the `plot.RoBMA()` function and `parameter = "weightfunction"` and `parameter = "PET-PEESE"`.

## version 1.2.1
### Fixes
- check_setup function not working at all

## version 1.2.0
### Changes
- the studies's true effects are now marginalized out of the random effects models and are no longer estimated (see Appendix A of our [manuscript](https://doi.org/10.1037/met0000405) for more details). As a results, arguments referring to the true effects are now disabled.
- all models are now being estimated using the likelihood of effect sizes (instead of test-statistics as usually defined). We reproduced the simulation study that we used to evaluate the method performance and it achieved identical results (up to MCMC error, before marginalizing out the true effects). A big advantage of using the normal likelihood for effect sizes is a considerable speed up of the whole estimation process.
- as a results of these two changes, the results of the models will differ to those of pre 1.2.0 version

### Fixes
- autofit being turn on if any control argument was specified

## version 1.1.2
### Fixes
- vdiffr not being used conditionally in unit tests

## version 1.1.1
### Fixes
- inability to fit a model without specifying a seed
- inability to produce individual model plots due to incompatibility with the newer versions of ggplot2  

## version 1.1.0
### Features
- parallel within and between model fitting using the parallel package with 'parallel = TRUE' argument

## version 1.0.5
### Fixes:
- models being fitted automatically until reaching R-hat lower than 1.05 without setting max_rhat and autofit control parameters
- bug preventing to draw a bivariate plot of mu and tau
- range for parameter estimates from individual models no containing 0 (or 1 in case of OR measured effect sizes)
- inability to fit a model with only null mu distributions if correlation or OR measured effect sizes were specified
- ordering of the estimated and observed effects when both of them are requested simultaneously
- formatting of this file (NEWS.md)

### Improvements:
- priors plot: parameter specification, default plotting range, clearer x-axis labels in cases when the parameter is defined on transformed scale
- parameters plots: probability scale always ends at the same spot as is the last tick on the density scale
- adding warnings if any of the specified models has Rhat higher than 1.05 or the specified value
- grouping the same warnings messages together

## version 1.0.4
### Fixes:
- inability to run models without the silent = TRUE control

## version 1.0.3
### Features:
- x-axis rescaling for the weight function plot (by setting 'rescale_x = TRUE' in the 'plot.RoBMA' function)
- setting expected direction of the effect in for RoBMA function

### Fixes:
- marginal likelihood calculation for models with spike prior distribution on mean parameter which location was not set to 0
- some additional error messages 

### CRAM requested changes:
- changing information messages from 'cat' to 'message' from plot related functions
- saving and returning the 'par' settings to the user defined one in the base plot functions

## version 1.0.2
### Fixes:
- the summary and plot function now shows quantile based confidence intervals for individual models instead of the HPD provided before (this affects only 'summary'/'plot' with 'type = "individual"', all other confidence intervals were quantile based before)

## version 1.0.1
### Fixes:
- summary function returning median instead of mean

## version 1.0.0 (vs the osf version)
### Fixes:
- incorrectly weighted theta estimates
- models with non-zero point prior distribution incorrectly plotted using when "models" option in case that the mu parameter was transformed

### Additional features:
- analyzing OR
- distributions implemented using boost library (helps with convergence issues)
- ability to mute the non-suppressible "precision not achieved" warning messages by using "silent" = TRUE inside of the control argument
- vignettes

### Notable changes:
- the way how the seed is set before model fitting (the simulation study will not be reproducible with the new version of the package)
