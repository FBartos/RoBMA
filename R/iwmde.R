# ============================================================================ #
# IWMDE File Map
# ============================================================================ #

# The old monolithic iwmde.R implementation is split by concern:
# - iwmde-context.R: fitted-object context and availability checks
# - iwmde-density-attributes.R: method and result-attribute adapters
# - iwmde-parameter-diagnostics.R: display-grid normalization and diagnostic
#   relabeling (the diagnostic orchestration is in iwmde-estimate.R)
# - iwmde-cache.R: estimator caches and target keys
# - iwmde-control.R: controls and reliability thresholds
# - iwmde-estimate.R: estimator facade and fixed-budget ordinate evaluation
# - iwmde-plan.R: estimate plans, row selection, and baseline contracts
# - iwmde-result.R: result schema versions and provenance hashes
# - iwmde-schema.R: plan, row-state, density, and diagnostic schemas
# - iwmde-diagnostics.R: reliability decisions and presentation helpers
# - iwmde-ordinate-failures.R: numerical status and diagnostic rows of failed
#   point ordinates
# - iwmde-parameters.R: parameter specs, active components, and focal priors
# - iwmde-grid.R: support transforms, normalization grids, and display grids
# - iwmde-density.R: row states and density estimator assembly
# - iwmde-aggregation.R: density aggregation and Monte Carlo error
# - iwmde-weights.R: IWMDE conditional-weight kernels and typed fallbacks
# - iwmde-qcmde-normalization.R: qCMDE normalization and refinement
# - iwmde-qcmde-range.R: qCMDE row laws, normalization range, and tail extension
# - iwmde-qcmde-closed-form.R: closed-form normalization of exact Normal rows
# - iwmde-retained-location.R: exact Gaussian retained-intercept conditioning
# - iwmde-q-grid.R: q-grid batch dispatch and GLMM conditional batches
# - iwmde-predictor.R: predictor/location fast paths and candidate construction
# - iwmde-predictor-known-v.R: known-V predictor kernels
# - iwmde-random-affine.R: exact affine random-covariance q-grid
# - iwmde-random-factor.R: exact coefficient-state random-covariance q-grid
# - iwmde-replacement.R: scalar, linear, and indexed replacement maps
# - iwmde-likelihood.R: likelihood and prior ordinates for replacement rows
# - iwmde-active.R: active prior branches, row parameters, selection, and formula inputs
# - iwmde-ordered.R: ordered scalar targets, conditional charts, and kernels
# - iwmde-ordered-block.R: positive conditional block integration for regular
#   ordered offsets
# - iwmde-plot.R: KDE, histogram, and trapezoid helpers
# - marginal-means-iwmde.R: marginal-means density integration and provenance
#
# This file is intentionally kept as an index for maintainers and old mental maps.
# ============================================================================ #
