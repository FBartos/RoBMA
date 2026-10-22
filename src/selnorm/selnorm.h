#ifndef SELNORM_H_
#define SELNORM_H_

#include <cstddef>
#include <vector>
#include "selnorm-fma.h"

// Longest step partition the batched tail kernel hoists cutoffs for; longer
// partitions use the unchanged per-element evaluator.
#define SELNORM_TAIL_MAX_BINS 64

// Largest certified factor rank the deterministic factor routes accept. The
// affordable quadrature at a given rank is decided by the node budget and the
// runtime support shape, not by this cap; the cap bounds the workspaces the
// kernels and the JAGS distributions allocate per block.
#define SELNORM_FACTOR_MAX_RANK 8

// Node budget one deterministic factor rule may spend. A tensor rule costs
// `order^rank` and the nested forest rule `order^(depth + 1)`; both are capped
// here. The value reproduces the per-rank rule counts the ladder used to carry
// as constants: `63^3` and `21^4` are the largest rules the rank three and
// rank four sequences reached, and the next rung of either exceeds it.
#define SELNORM_FACTOR_NODE_BUDGET 262144.0

// Whether one rule of this order fits the budget at the given cost exponent.
inline bool selnorm_factor_rule_affordable(int order, int exponent)
{
  double nodes = 1.0;
  for (int axis = 0; axis < exponent; ++axis) {
    nodes *= static_cast<double>(order);
    if (nodes > SELNORM_FACTOR_NODE_BUDGET) return false;
  }
  return true;
}

enum SelKernelMode {
  SELKERNEL_NORMAL           = 0,
  SELKERNEL_STEP             = 1,
  SELKERNEL_PHACK_POWER      = 2,
  SELKERNEL_STEP_PHACK_POWER = 3
};

enum SelVectorRule {
  SELVECTOR_PRODUCT = 0,
  SELVECTOR_BEST_ONE_SIDED = 1,
  SELVECTOR_BEST_TWO_SIDED = 2
};

double cpp_selnorm_normal_interval_log_prob(
  double lower, double upper, double mean, double sd
);

double cpp_selnorm_normal_interval_quantile(
  double lower, double upper, double mean, double sd, double u
);

struct SelNormKernelData {
  int n_bins;
  int n_segments;
  int effect_sign;
  int q;
  const double *z_lower;
  const double *z_upper;
  const double *phack_z_source;
  const double *phack_z_dest;
  const double *segment_bounds;
  const int *segment_step_bin;
  const int *segment_phack_region;
  const double *segment_step_bin_real;
  const double *segment_phack_region_real;
  bool trusted_step_partition;
  bool telescope_probabilities;
  // Set only on a call-owned R post-fit integration copy.
  bool bounded_cdf = false;
};

// Whether a finite cutoff separates two different adjacent step weights. A
// block without such a transition has a constant selection weight over the
// actual axes, so its quadrature has no aliasing to certify.
bool cpp_selnorm_step_has_transitions(
  const double *omega, const SelNormKernelData &data, int omega_stride = 1
);

// Additional finite-rule selection-transition heuristic. The actual axes are
// standard Gaussian; loading is column-major dimension by rank. Existing error
// and exterior-mass gates remain separate. This is not a global error bound.
// `width_scale` multiplies every step width (residual SD over |loading|) the
// rule compares with a node gap; the default 1 is the rule itself, a larger
// scale relaxes it (the product-rule certificate asks for edges no narrower
// than a fixed fraction of the gap, see selnorm-mv.cc).
bool cpp_selnorm_step_rule_resolves_transitions(
  const double *nodes, int offset, int order,
  const double *mean, const double *residual_sd, const double *loading,
  int dimension, int rank, const double *selection_se,
  const double *omega, const SelNormKernelData &data, int omega_stride = 1,
  double width_scale = 1.0
);

double cpp_selnorm_bounded_cdf_absolute_error();

bool selnorm_is_descending_step_partition(const double *z_lower,
                                          const double *z_upper,
                                          int n_bins);

double cpp_selnorm_normal_lpdf(double x, double mean, double sd);

double cpp_selnorm_affine_normal_lpdf_log_scale(
  double z, double sei, double mean, double sd, double log_scale
);

double cpp_selnorm_affine_normal_pdf_log_scale(
  double z, double sei, double mean, double sd, double log_scale
);

double cpp_selnorm_kernel_lpdf(
  double y,
  double mu_num,
  double sigma_num,
  double mu_norm,
  double sigma_norm,
  double sei,
  double weight,
  const double *omega,
  int obs_bin,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  int omega_stride = 1,
  bool validate_omega = true
);

double cpp_selnorm_kernel_log_norm(
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  int omega_stride = 1,
  bool validate_omega = true
);

double cpp_selnorm_kernel_cdf(
  double q,
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  int omega_stride = 1,
  bool lower_tail = true,
  bool validate_omega = true
);

bool cpp_selnorm_step_cdf_plan(
  double mean,
  double sd,
  double sei,
  const double *omega,
  const SelNormKernelData &data,
  double *lower_score,
  double *upper_score,
  double *log_weight,
  int *n_groups,
  int omega_stride = 1,
  bool validate_omega = true
);

double cpp_selnorm_step_log_norm(
  double mean,
  double sd,
  double sei,
  const double *omega,
  const SelNormKernelData &data,
  int omega_stride = 1,
  bool validate_omega = true
);

bool cpp_selnorm_step_cdf_telescope_plan(
  double mean,
  double sd,
  double sei,
  const double *omega,
  const SelNormKernelData &data,
  double *boundary_tail,
  double *omega_diff,
  double *omega_last,
  double *normalizer,
  int omega_stride = 1,
  bool validate_omega = true
);

// Multiply one row's normalizers over a quadrature grid. The caller validates
// omega once per likelihood state; false retains the scalar/log-scale fallback.
bool cpp_selnorm_step_normalizer_product(
  const double *mean, std::size_t count, double sd, double sei,
  const double *omega, const SelNormKernelData &data, long double *product
);

// Double-precision node products for the product ladders: the arithmetic of the
// batch routine behind cpp_selnorm_step_normalizer_product with the node
// products in double. The constants of the rows that do not depend on the node
// are built once per kernel call and shared by every rung and shifted copy of
// that call. A plan belongs to one call (the batch entries are threaded) and is
// never shared. It is not usable when `rows * log2(smallest positive omega)` is
// below -960, where the node products leave the double range for any row, or
// when a row's geometry is unsupported (the caller then keeps the long double
// routine above).
class SelNormDoubleProducts {
 public:
  void build(const double *residual_sd, const double *selection_se, int rows,
             const double *omega, const SelNormKernelData &data);
  bool usable() const { return usable_; }
  // Multiply the normalizers of `row` at the `count` locations in `mean` into
  // `product`. The normalizers of the whole row are computed into the
  // caller-owned `normalizer` (at least `count` doubles) and checked as the
  // batch routine does before the first product is multiplied; every updated
  // node product must lie in [2^-960, 2^960]. False means "out of range": the
  // products are untouched after a failed check and partly updated after a
  // product out of range, and the caller discards every product of the rung
  // and evaluates it with the long double routine.
  bool multiply(int row, const double *mean, std::size_t count, double *product,
                double *normalizer) const;

 private:
  bool usable_ = false;
  // The changed boundaries and their weight differences belong to the selection
  // weights alone. Per row, `row_stride_` doubles hold the inverse residual SD
  // and then one scaled cutoff per changed boundary.
  int boundaries_ = 0;
  int row_stride_ = 0;
  double last_ = 0.0;
  double difference_[SELNORM_TAIL_MAX_BINS];
  std::vector<double> rows_;
};

double cpp_selnorm_step_cdf_from_telescope_plan(
  double q,
  double mean,
  double sd,
  double sei,
  const SelNormKernelData &data,
  const double *boundary_tail,
  const double *omega_diff,
  double omega_last,
  double normalizer,
  bool lower_tail = true
);

double cpp_selnorm_step_cdf_from_plan(
  double q,
  double mean,
  double sd,
  double sei,
  const SelNormKernelData &data,
  const double *lower_score,
  const double *upper_score,
  const double *log_weight,
  int n_groups,
  bool lower_tail = true
);

void cpp_selnorm_kernel_log_tail_pair(
  double q,
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  double *log_lower,
  double *log_upper,
  int omega_stride = 1,
  bool validate_omega = true
);

void cpp_selnorm_kernel_moments(
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  double *moment_mean,
  double *moment_second,
  int omega_stride = 1,
  bool validate_omega = true,
  double *moment_variance = nullptr
);

double cpp_selnorm_kernel_threshold(
  double z_threshold,
  double mean,
  double sd,
  double sei,
  const double *omega,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  double *inverse_weight,
  int omega_stride = 1,
  bool validate_omega = true
);

double cpp_selnorm_kernel_rng(
  double mean,
  double sd,
  double sei,
  const double *omega,
  double u_bin,
  double u_interval,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  int omega_stride = 1,
  bool validate_omega = true
);

double cpp_selnorm_kernel_rng_workspace(
  double mean,
  double sd,
  double sei,
  const double *omega,
  double u_bin,
  double u_interval,
  double alpha,
  int phack_kind,
  int kernel_mode,
  const SelNormKernelData &data,
  double *mass,
  double *lower,
  double *upper,
  int omega_stride = 1,
  bool validate_omega = true
);

bool cpp_selnorm_step_rng_prepare_partition(
  double mean, double sd, double sei, const double *omega,
  const SelNormKernelData &data, double *mass, double *lower, double *upper,
  int *n_groups, double *normalizer, double *log_normalizer
);

double cpp_selnorm_step_rng_from_partition(
  double mean, double sd, int effect_sign, double u_bin, double u_interval,
  const double *mass, const double *lower, const double *upper,
  int n_groups, double normalizer
);

double cpp_selnorm_step_log_norm_rng_workspace(
  double mean,
  double sd,
  double sei,
  const double *omega,
  double u_bin,
  double u_interval,
  const SelNormKernelData &data,
  double *mass,
  double *lower,
  double *upper,
  double *log_normalizer,
  int omega_stride = 1,
  bool validate_omega = true
);

#endif
