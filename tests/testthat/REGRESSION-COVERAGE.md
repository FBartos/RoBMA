# Regression evidence inventory

This inventory identifies the high-information evidence expected for public
post-fit behavior. Structural checks supplement these oracles; merely returning
an object or rendering without error is not sufficient representative evidence.

| Area | Standard evidence | Certification extension |
|---|---|---|
| Residuals and standardized residuals | Numerical metafor comparisons for normal, regression, multilevel, PET-PEESE, selection, and GLMM cases where the estimand exists | Extended interaction, negative-direction, multilevel-regression, and high-draw known-V cases |
| Predictions | Numerical metafor comparisons for fitted and new-data predictions, plus target-specific known-V identities | Extended model families and high-draw known-V targets |
| Hat values, DFBETAS, VIF, and influence | Numerical metafor or independent GLS/moment oracles for representative supported families | Extended parameterizations and high-draw known-V diagnostics |
| Heterogeneity summaries | Numerical metafor comparisons for ordinary, scale, multilevel, selection, and GLMM models | Negative-direction and secondary model variants |
| Funnel, Q-Q, radial, regression, forest, and z plots | Analytic funnel PIT-projection, CDF-mixture, and posterior-integration identities plus human-reviewed vdiffr snapshots, including side-by-side metafor views where displays share an estimand | Secondary model families and customization galleries |
| Prior plots | Base and ggplot snapshots for outcome, moderator, and publication-bias priors | Structural component-selection tests cover additional prior families |
| Multivariate diagnostics | Analytic known-V/known-R data checks and representative forest coverage | Human-reviewed funnel, Q-Q, forest, and regression snapshots from high-draw metafor fixtures |
| Estimate-unit LOO | Analytic Gaussian conditional-density identities and cached-fit target contracts | Five exact observation-deletion refits of Kearon US/HCS and Ishak HAR models in the `loo-exact-refits` case |
| Exact selection likelihood | Independent `mvtnorm` rectangle-probability identities for the dependent-block normalizer, analytic diagonal reduction, selected-response region masses and singular-covariance support, specialized cluster-posterior moments, dense-versus-compiled conditional random-effect algebra, a cached random-formula `bselmodel.mv` fit, and a cached `RoBMA.mv` product-space fit combining unadjusted, selection, PET, PEESE, moderator, and independently gated random branches. The fitted checks cover manual posterior averaging, model enumeration, summary, LOO, bridge availability contracts, prediction, diagnostics, plots, and LOO-PIT residuals. | Larger dependency blocks and high-draw bridge stability checks |
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
