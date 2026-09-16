#ifndef ROBMA_SELNORM_MV_H
#define ROBMA_SELNORM_MV_H

#include <cstddef>
#include <cstdint>
#include <memory>

// One process-wide exact/coarse budget. Global bytes include all pool/index
// blocks requested upstream; component counters do not duplicate that storage.
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
// the cache mutex; restore validates every shard before replacing pool state.
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
// an absolute covariance-perturbation bound. That error is separate from the
// estimated interior quadrature error and bounded exterior Gaussian mass.
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
// output untouched and requests direct evaluation of the same mixture.
class SelNormNormalProjectionWorkspace {
  struct Impl;
  std::unique_ptr<Impl> impl;
public:
  SelNormNormalProjectionWorkspace(const double *z, int size);
  ~SelNormNormalProjectionWorkspace();
  bool density(const double *mean, int count, double sd, double sei,
               const double *log_coefficients, long double *output);
};

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
  int vector_rule = 0
);

double cpp_selnorm_factor_step_lpdf(
  const double *x,
  const double *mean,
  const double *residual_sd,
  const double *loading,
  int dimension,
  int rank,
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
  int vector_rule = 0
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
  int threads = 1
);

#endif
