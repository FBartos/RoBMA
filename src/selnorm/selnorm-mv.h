#ifndef ROBMA_SELNORM_MV_H
#define ROBMA_SELNORM_MV_H

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <vector>

// Immutable quadrature weights shared by R and JAGS. A null node pointer
// requests exp(log weight); otherwise preserve the typed shifted expression.
// Exact owned keys (including node values and shift) decide reuse. Cache
// entries are limited to 32 and their accounted payload to 1 MiB per
// precision, plus bounded container metadata; caller owners
// remain valid after eviction. Double and long-double arithmetic stay separate.
std::shared_ptr<const std::vector<double> > cpp_selnorm_rule_weights_double(
  const double *nodes, const double *log_weights, int count, double shift);
std::shared_ptr<const std::vector<long double> > cpp_selnorm_rule_weights_long_double(
  const double *log_weights, int count);

// One process-wide exact/coarse budget. Global bytes include all directly
// allocated shared blocks; component counters do not duplicate that storage.
struct SelNormCacheStats {
  std::size_t entries = 0;
  std::uint64_t hits = 0, misses = 0, evictions = 0, resets = 0, allocation_failures = 0;
};
struct SelNormCacheInfo {
  std::size_t capacity_bytes = 0, allocated_bytes = 0, peak_bytes = 0;
  SelNormCacheStats exact, coarse;
};

// Clear mask: 0 neither, 1 exact, 2 coarse, 3 both. A changed capacity clears
// both components and releases all storage regardless of the mask.
SelNormCacheInfo cpp_selnorm_cache_control(
  bool set_capacity, std::size_t capacity_bytes, unsigned int clear_mask
);

// Snapshot buffers are call-owned views. R allocates export memory outside
// the cache mutex; restore validates every shard before replacing cache storage.
struct SelNormCacheBlob { const unsigned char *data; std::size_t size; };
struct SelNormCacheRestoreInfo {
  std::uint64_t snapshots = 0, available_exact = 0, available_coarse = 0;
  std::uint64_t restored_exact = 0, restored_coarse = 0, skipped = 0;
};
std::size_t cpp_selnorm_cache_snapshot_size();
bool cpp_selnorm_cache_snapshot_write(unsigned char *data, std::size_t size);
SelNormCacheRestoreInfo cpp_selnorm_cache_restore(const SelNormCacheBlob *snapshots, std::size_t count);

// Corrected factor-box rejection for a fixed retained context. Status:
// 1 accepted, 0 attempt cap exhausted, -1 optional geometry unavailable before
// RNG consumption, -2 invalid proposal arithmetic. The partition is call-owned.
int cpp_selnorm_factor_box_rng(
  const double *mean, const double *covariance, int dimension,
  const double *selection_se, const double *omega, int n_bins,
  const double *z_lower, const double *z_upper, int effect_sign,
  int max_attempts, double (*uniform)(), double *output
);

// Rao-Blackwellized scalar marginals of a finite-vector selected Gaussian.
// The error is the largest absolute density MCSE, relative to the curve peak.
struct SelNormZProjection {
  const double *z;
  int size;
  double *density;
  double relative_error;
  bool probability;
  const double *factor_loading;
  const int *loading_support;
  const double *residual_sd;
  int rank;
  const double *nodes;
  const double *log_weights;
  int order;
};

// Optional integration rules for an ordinary dense covariance. Covariance
// envelopes bound monotone product normalizers. Nonmonotone products require
// a covariance-perturbation bound: strictly positive weights permit a relative
// bound, while zero weights retain the absolute bound. Covariance error is
// separate from estimated interior quadrature error and bounded exterior mass.
struct SelNormDenseIntegration {
  const double *nodes = nullptr;
  const double *log_weights = nullptr;
  const double *orders = nullptr;
  int rule_count = 0;
  double tolerance = 0.0;
  bool used_envelope = false;
  double quadrature_change = 0.0;
  double covariance_width = 0.0;
  double tail_bound = 0.0;
  // Explicit R post-fit request; no sampler or RNG sets this flag.
  bool bounded_cdf = false;
  double cdf_log_error = 0.0;
  double cdf_relative_error = 0.0;
};

// Private corrected-sampler policy. Fixed before retained sampling; never an
// ambient approximation mode. Steps are relative to the largest original SE
// (mean) and its square (diagonal); weight steps are on the natural-log scale.
struct SelNormCoarseSettings {
  double mean_step = 0.05;
  double diagonal_step = 0.05;
  double log_weight_step = 0.01;
  unsigned int max_rules = 3;
};

double cpp_selnorm_mnorm_step_surrogate_lpdf(
  const double *x, const double *mean, const double *covariance_lower,
  int dimension, const double *selection_se, const double *omega, int n_bins,
  const double *z_lower, const double *z_upper, const int *obs_bin,
  int effect_sign, bool telescope_probabilities, int kernel_mode, int vector_rule,
  const SelNormDenseIntegration &quadrature, const SelNormCoarseSettings &settings
);

// Call-owned grid topology and scratch space for quadrature mixtures. The grid
// must remain unchanged throughout this workspace's lifetime. False leaves
// output untouched and requests direct evaluation of the same mixture. Only a
// serial R owner supplies a polling function; shared/JAGS owners leave it null.
class SelNormNormalProjectionWorkspace {
  struct Impl;
  std::unique_ptr<Impl> impl;
public:
  SelNormNormalProjectionWorkspace(const double *z, int size,
                                  void (*poll)() = nullptr,
                                  std::size_t poll_interval = 1000000);
  ~SelNormNormalProjectionWorkspace();
  bool density(const double *mean, int count, double sd, double sei,
               const double *log_coefficients, long double *output);
};

// Test hook: one rung of the product rule with the shift-amplitude certificate
// of the factor (or rank-one cluster) kernel, on the same routines the kernels
// use. Groups are the latent axes; `column` is the loading column of each.
struct SelNormCertificateProbe {
  int path = -1;  // 0 declined, 1 nested forest, 2 tensor direct, 3 tensor log, 4 cluster
  bool transitions = false;
  // The sharpness gate at this rung: the transition rule with every step width
  // multiplied by the scale of the ladder and of the refinement tier (the tier
  // scale narrows with the number of rows loading on an axis, see
  // selnorm-mv.cc).
  bool resolves_ladder = false;
  bool resolves_tier = false;
  double tier_width_scale = std::numeric_limits<double>::quiet_NaN();
  double log_normalizer = 0.0;
  double reference_log_normalizer = 0.0;
  double tail_bound = 0.0;
  double amplitude_sum = 0.0;
  bool finite = false;
  bool early_reject = false;
  bool accepted = false;
  bool accepted_without_early_exit = false;
  std::vector<int> column;
  std::vector<double> mean, gap, half, quarter, amplitude;
  std::vector<double> reuse_value, full_value;
  // Whether the rung (and its shifted copies) came from the double precision
  // node products, and the shifted values of the same requests computed with
  // the long double products.
  bool lean = false;
  std::vector<double> reference_value;
  // Averages over the shift grids M = 2 and 4 of the refinement tier, at the
  // gaps of the rung's tilted means (forests of depth <= 1, tensors of rank <= 2).
  bool refined = false;
  double refined_2 = std::numeric_limits<double>::quiet_NaN();
  double refined_4 = std::numeric_limits<double>::quiet_NaN();
  // Tilted moments of the rung in loading-column order (tensors and forests of
  // depth at most one, and the rank-one cluster rule).
  bool has_moments = false;
  std::vector<double> moment_mean, moment_covariance;
  // Rank-one cluster rule: the tilted mean and variance of the latent
  // coordinate under the shifted rule of each shift request (the nodes are
  // x_j + shift), from the evaluator of the ladder and from the long double
  // evaluator.
  std::vector<double> shifted_moment_mean, shifted_moment_variance;
  std::vector<double> reference_shifted_moment_mean, reference_shifted_moment_variance;
};

bool cpp_selnorm_certificate_probe(
  bool cluster, const double *mean, const double *residual_sd,
  const double *loading, int dimension, int rank, const double *selection_se,
  const double *omega, int n_bins, const double *z_lower, const double *z_upper,
  int effect_sign, bool telescope_probabilities, const double *quadrature_nodes,
  const double *quadrature_log_weights, const double *quadrature_orders,
  int quadrature_rule_count, int rule, double relative_tolerance,
  const int *shift_axis, const double *shift_delta, int shift_count,
  bool precomputed_weights, bool refinement, SelNormCertificateProbe *out
);

double cpp_selnorm_certificate_node_gap(const double *nodes, int order, double mean);

double cpp_selnorm_mnorm_step_lpdf(
  const double *x,
  const double *mean,
  const double *covariance_lower,
  int dimension,
  const double *selection_se,
  const double *omega,
  int n_bins,
  const double *z_lower,
  const double *z_upper,
  const int *obs_bin,
  int effect_sign,
  bool telescope_probabilities,
  int kernel_mode,
  const double *qmc,
  int points,
  int scrambles,
  double *relative_mcse,
  SelNormZProjection *projection = nullptr,
  double *log_normalizer = nullptr,
  int vector_rule = 0,
  SelNormDenseIntegration *integration = nullptr,
  const double *log_normalizer_override = nullptr
);

// Stable log of the powered Gaussian one-factor integral
//   log int prod_i N(residual_i; loading_i * g, sd_i)^power_i N(g; 0, 1) dg,
// the cluster likelihood with one shared unit-normal latent. Rows with power 0
// are omitted, a null `power` means unit powers, and the loadings may be
// signed. The residual after the weighted latent projection is formed
// directly, so nothing is subtracted from an explained quadratic, and no
// quotient by a loading is formed: the form is evaluated in standardized units
// when no loading exceeds the smallest sd and in raw units scaled by the
// largest loading otherwise (see the definition). `work` is optional scratch
// of 2 n values. Returns NaN for an invalid input (non-finite value, non-positive
// sd, negative power) and -Inf when the result is not representable.
double cpp_gaussian_latent_loglik(
  int n,
  const double *residual,
  const double *sd,
  const double *loading,
  const double *power,
  double *work = nullptr
);

double cpp_selnorm_cluster_step_lpdf(
  const double *x,
  const double *mean,
  const double *residual_sd,
  const double *loading,
  int dimension,
  const double *selection_se,
  const double *omega,
  int n_bins,
  const double *z_lower,
  const double *z_upper,
  const int *obs_bin,
  int effect_sign,
  bool telescope_probabilities,
  int kernel_mode,
  const double *quadrature_nodes,
  const double *quadrature_log_weights,
  const double *quadrature_orders,
  int quadrature_rule_count,
  const double *qmc,
  int initial_points,
  int max_points,
  int scrambles,
  double relative_tolerance,
  double *relative_mcse,
  double *relative_change,
  const long double *quadrature_weights = nullptr,
  double *log_normalizer = nullptr,
  int vector_rule = 0,
  bool *deterministic = nullptr,
  const double *log_normalizer_override = nullptr
);

// Mean and covariance (column-major) of the latent coordinates under the tilted
// integrand of the finest deterministic product rung, in the order of the
// loading columns, and the number of finite rungs the ladder evaluated. The
// randomized fallback of the factor kernel proposes from them (the prior, a
// Gaussian with these moments and a copy with 1.5 times their standard
// deviations) when at least three rungs were finite, the finest covers the
// exterior mass (the tail gate of the ladder) and the covariance is positive
// definite. The rank-one cluster kernel carries its own to the factor kernel it
// delegates to, empty when its finest rung does not cover the exterior mass.
struct SelNormTiltedMoments {
  int finite_rungs = 0;
  std::vector<double> mean;
  std::vector<double> covariance;
};

// Zero loading columns are compacted out before integration: they only add a
// latent axis that integrates to one, so `rank` is the declared rank and the
// evaluated one follows the columns that carry dependence.
//
// The optional `deterministic` flag of the factor and cluster kernels is set
// exactly when the returned density is a function of the call's inputs alone: it
// came from a deterministic ladder (of the product or the best weight rule), the
// refinement tier or the sparse grid.
// The randomized fallback never sets it, whatever its Monte Carlo standard
// error (which can be exactly zero), and neither does a call that returns
// before a normalizer is computed. Only such results may be remembered.
//
// A non-null `log_normalizer_override` (finite, used by the corrected sampler's
// surrogate and its test hook) replaces the selection normalizer: the call
// returns the Gaussian log density plus the observed log weights minus the
// override, after the NORMAL early return and before any normalizer work. The
// cluster kernel takes the same argument at the same place.
double cpp_selnorm_factor_step_lpdf(
  const double *x,
  const double *mean,
  const double *residual_sd,
  const double *declared_loading,
  int dimension,
  int declared_rank,
  const double *selection_se,
  const double *omega,
  int n_bins,
  const double *z_lower,
  const double *z_upper,
  const int *obs_bin,
  int effect_sign,
  bool telescope_probabilities,
  int kernel_mode,
  const double *quadrature_nodes,
  const double *quadrature_log_weights,
  const double *quadrature_orders,
  int quadrature_rule_count,
  const double *qmc,
  int initial_points,
  int max_points,
  int scrambles,
  double relative_tolerance,
  double *relative_mcse,
  double *relative_change,
  double *log_normalizer = nullptr,
  int vector_rule = 0,
  const SelNormTiltedMoments *proposal_moments = nullptr,
  bool *deterministic = nullptr,
  const double *log_normalizer_override = nullptr
);

// Surrogates of the rank-one cluster and factor blocks of the corrected sampler:
// the kernel's own numerator (Gaussian log density plus observed log weights)
// minus the log normalizer of ONE rule of the block's own quadrature schedule,
// rule min(settings.max_rules, rule_count) at the actual parameters, with no
// change test, certificate, refinement, sparse grid, QMC, memo or rounding. The
// result is a deterministic function of the arguments. `*exact_fallback` is set
// (and NaN returned) when that rule is unavailable: a non-product vector rule,
// a rule the node budget does not afford for the block, or a non-finite log
// normalizer; the caller then evaluates the exact density. A NORMAL kernel
// returns its numerator, and an invalid input -Inf, as the exact kernel does.
// `log_normalizer_out` receives the rule's log normalizer (NaN when unavailable,
// 0 for the NORMAL kernel).
double cpp_selnorm_cluster_step_surrogate_lpdf(
  const double *x, const double *mean, const double *residual_sd,
  const double *loading, int dimension, const double *selection_se,
  const double *omega, int n_bins, const double *z_lower,
  const double *z_upper, const int *obs_bin, int effect_sign,
  bool telescope_probabilities, int kernel_mode,
  const double *quadrature_nodes, const double *quadrature_log_weights,
  const double *quadrature_orders, int quadrature_rule_count, int vector_rule,
  const SelNormCoarseSettings &settings, bool *exact_fallback,
  double *log_normalizer_out = nullptr
);

double cpp_selnorm_factor_step_surrogate_lpdf(
  const double *x, const double *mean, const double *residual_sd,
  const double *declared_loading, int dimension, int declared_rank,
  const double *selection_se, const double *omega, int n_bins,
  const double *z_lower, const double *z_upper, const int *obs_bin,
  int effect_sign, bool telescope_probabilities, int kernel_mode,
  const double *quadrature_nodes, const double *quadrature_log_weights,
  const double *quadrature_orders, int quadrature_rule_count, int vector_rule,
  const SelNormCoarseSettings &settings, bool *exact_fallback,
  double *log_normalizer_out = nullptr
);

double cpp_selnorm_gaussian_event_log_mass(
  const double *mean, const double *covariance_lower, int dimension,
  const double *selection_se, const double *omega, int n_bins,
  const double *z_lower, const double *z_upper, int effect_sign,
  int kernel_mode, int vector_rule, const double *lower,
  const double *upper, const double *qmc, int points, int scrambles,
  double *relative_mcse, const double *rank_one_loading = nullptr,
  SelNormDenseIntegration *integration = nullptr
);

// Whole sampling-error conditioning, evaluated by Gaussian conditional
// simulation. The retained covariance and auxiliary include all retained
// Gaussian sources; callers recover individual sources from delta.
// The declared candidate covariance is diag(diagonal) + L L'.
// No density with respect to a singular candidate Gaussian is evaluated.
struct SelNormConditionedQuadrature {
  const double *cluster_nodes;
  const double *cluster_log_weights;
  const double *cluster_orders;
  int cluster_rule_count;
  const double *factor_nodes;
  const double *factor_log_weights;
  const double *factor_orders;
  int factor_rule_count;
};

double cpp_selnorm_sampling_conditioned_lpdf(
  const double *x, const double *mean, const double *sampling_lower,
  const double *diagonal, const double *loading, int dimension, int rank,
  const double *sampling_auxiliary, const double *random_auxiliary,
  const double *selection_se, const double *omega, int n_bins,
  const double *z_lower, const double *z_upper, const int *obs_bin,
  int effect_sign, int kernel_mode, bool telescope_probabilities,
  int vector_rule, const int *group_index, const double *qmc,
  int initial_points, int max_points, int scrambles, double tolerance,
  double *relative_mcse, double *relative_change, double *sampling_effect,
  double *delta, double *log_normalizer,
  const SelNormConditionedQuadrature &quadrature
);

double cpp_selnorm_conditioned_log_normalizer(
  const double *mean, const double *diagonal, const double *loading,
  int dimension, int rank, const double *selection_se, const double *omega,
  int n_bins, const double *z_lower, const double *z_upper, int effect_sign,
  int kernel_mode, bool telescope_probabilities, int vector_rule,
  const int *group_index, const double *qmc, int initial_points, int max_points,
  int scrambles, double tolerance, double *relative_mcse, double *relative_change,
  const SelNormConditionedQuadrature &quadrature,
  const double *lower = nullptr, const double *upper = nullptr
);

void cpp_selnorm_factor_projection(
  const double *mean, const double *diagonal, const double *loading,
  int dimension, int rank, const double *selection_se, const double *omega,
  int n_bins, const double *z_lower, const double *z_upper, int effect_sign,
  int kernel_mode, bool telescope_probabilities, int vector_rule,
  const int *group_index, const double *qmc, int initial_points, int max_points,
  int scrambles, double tolerance, const SelNormConditionedQuadrature &quadrature,
  const double *z, int size, int kind, double *density,
  double *log_normalizer, double *relative_error
);

bool cpp_selnorm_covariance_envelope_components(
  const double *covariance_lower, int dimension, const double *selection_se,
  bool upper, double *residual_sd, double *loading, int *groups, int *type_index
);

bool cpp_selnorm_context_star_projection(
  const double *means, int contexts, const double *covariance_lower, int dimension,
  const double *selection_se, const double *omega, int n_bins,
  const double *lower, const double *upper, int effect_sign, bool telescope,
  const double *log_normalizers, const double *context_log_weights,
  const double *nodes, const double *log_weights, int order,
  const double *z, int size, double allowance, double *density,
  double *mass, double *compression_error, std::size_t *components,
  double *compact_log_normalizers = nullptr,
  bool bounded_cdf = false, double *cdf_factor_log_errors = nullptr,
  double *input_omitted_mass_error = nullptr,
  int threads = 1, const double *cached_rules = nullptr, double *saved_rules = nullptr
);

#endif
