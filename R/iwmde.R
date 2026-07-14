# ============================================================================ #
# IWMDE File Map
# ============================================================================ #

# The old monolithic iwmde.R implementation is split by concern:
# - iwmde-api-diagnostics.R: exported developer diagnostic surfaces
# - iwmde-context.R: fitted-object context and availability checks
# - iwmde-density-attributes.R: method and result-attribute adapters
# - iwmde-parameter-diagnostics.R: parameter diagnostic orchestration
# - iwmde-cache.R: estimator caches, target keys, and row selection
# - iwmde-control.R: controls and reliability thresholds
# - iwmde-estimate.R: estimator facade and adaptive ordinate evaluation
# - iwmde-schema.R: plan, row-state, density, and diagnostic schemas
# - iwmde-diagnostics.R: reliability decisions and presentation helpers
# - iwmde-parameters.R: parameter specs, active components, and focal priors
# - iwmde-grid.R: support transforms, normalization grids, and display grids
# - iwmde-density.R: row states and density estimator assembly
# - iwmde-aggregation.R: density aggregation and Monte Carlo error
# - iwmde-weights.R: IWMDE conditional-weight kernels and typed fallbacks
# - iwmde-qcmde-normalization.R: qCMDE normalization and refinement
# - iwmde-q-grid.R: q-grid batch dispatch and GLMM conditional batches
# - iwmde-predictor.R: predictor/location fast paths and candidate construction
# - iwmde-predictor-known-v.R: known-V predictor kernels
# - iwmde-replacement.R: scalar, linear, and indexed replacement maps
# - iwmde-likelihood.R: likelihood and prior ordinates for replacement rows
# - iwmde-active.R: active prior branches, row parameters, selection, and formula inputs
# - iwmde-plot.R: KDE/histogram data and base plotting
# - marginal-means-iwmde.R: marginal-means density integration and provenance
#
# This file is intentionally kept as an index for maintainers and old mental maps.
# ============================================================================ #
