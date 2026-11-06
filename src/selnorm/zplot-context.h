// Native zplot context projection.
//
// The certified point/continuous-context projection controller of
// .zplot_context_projection() (R/zplot-context.R), ported statement by
// statement so that a whole block of posterior draws can be projected in one
// call and in parallel over draws. The R controller remains the reference: the
// port reproduces its arithmetic exactly (R's long-double reductions, the BLAS
// and LAPACK routines behind chol(), forwardsolve(), %*%, tcrossprod() and
// norm(), and Rmath), its memo tables and their effect on control flow, the
// refinement routing and every acceptance gate. Change both together; the
// equivalence tests in test-00-zplot-context-native.R pin them.
//
// Nothing here calls the R API. Rmath functions are called only with arguments
// for which they cannot warn (see the controller source).
#ifndef ROBMA_SELNORM_ZPLOT_CONTEXT_H
#define ROBMA_SELNORM_ZPLOT_CONTEXT_H

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

// A failure the R route reports as an error, carrying its exact message.
class ZplotNativeFailure : public std::runtime_error {
 public:
  explicit ZplotNativeFailure(const std::string &message) : std::runtime_error(message) {}
};

// One call of the context star kernel (the body of
// RoBMA_selnorm_context_star_zplot without its R-owned raw-rule cache objects).
struct ZplotStarCall {
  const double *means = nullptr;          // contexts x k, column-major
  int contexts = 0;
  int k = 0;
  const double *covariance_lower = nullptr;  // packed lower triangle by columns
  const double *selection_se = nullptr;      // k
  const double *omega = nullptr;             // bins
  int bins = 0;
  const double *lower = nullptr;             // bins
  const double *upper = nullptr;             // bins
  int effect_sign = 1;
  bool telescope = false;
  const double *log_normalizers = nullptr;   // contexts, or nullptr for own-A
  const double *context_log_weights = nullptr;  // contexts
  const double *nodes = nullptr;             // order
  const double *log_weights = nullptr;       // order
  int order = 0;
  const double *z = nullptr;                 // size
  int size = 0;
  double allowance = 0.0;
  bool bounded_cdf = false;
  const double *cached_rules = nullptr;
  double *saved_rules = nullptr;
};

struct ZplotStarResult {
  bool available = false;
  std::vector<double> density;
  double mass = 0.0;
  double compression_error = 0.0;
  double omitted_mass_error = 0.0;
  double components = 0.0;
  std::vector<double> compact_log_normalizers;
  std::vector<double> cdf_factor_log_error;
  std::vector<double> cdf_log_error;
  std::vector<double> cdf_numerator_log_error;
};

// Input checks of the R entry, in its order; throws ZplotNativeFailure with
// the entry's messages.
void zplot_star_validate(const ZplotStarCall &call);
// The thread count the R entry gives the kernel.
int zplot_star_threads(const ZplotStarCall &call);
// Validated call -> result. A kernel exception becomes a ZplotNativeFailure
// with the entry's "Deterministic star projection failed: ..." message.
void zplot_star_evaluate(const ZplotStarCall &call, int threads, ZplotStarResult *result);

// RoBMA_selnorm_covariance_envelope_components without R objects.
struct ZplotEnvelopeComponents {
  bool available = false;
  int groups = 0;
  int rank = 0;
  std::vector<double> residual_sd;        // k
  std::vector<double> residual_variance;  // k
  std::vector<double> loading;            // k x rank, column-major
  std::vector<int> type_index;            // k, one-based
  std::vector<double> covariance_residual;         // k x k
  std::vector<double> covariance_residual_radius;  // k x k
};
// Throws ZplotNativeFailure("Compact envelope allocation failed.") when the
// R entry would report that.
void zplot_covariance_envelope_components(const double *covariance_lower, int k,
  const double *selection_se, bool upper, ZplotEnvelopeComponents *result);

// One Gauss-Hermite rule as .gauss_hermite_nodes() and
// .gauss_hermite_nodes_only() return it.
struct ZplotGaussHermiteRule {
  int order = 0;
  std::vector<double> nodes;
  std::vector<double> log_weights;
  std::vector<double> nodes_only;
};

// Inputs shared by every draw of a block.
struct ZplotProjectionSetup {
  const double *z = nullptr;
  int size = 0;
  const double *sei = nullptr;
  int k = 0;
  // BayesTools::selection_native_static_args() of the selection context.
  const double *z_lower = nullptr;
  const double *z_upper = nullptr;
  int bins = 0;
  int effect_sign = 1;
  bool telescope = false;
  // selection$z_lower and selection$sign as the transition rules read them.
  const double *rule_z_lower = nullptr;
  int rule_z_lower_length = 0;
  double rule_sign = 1.0;
  double tolerance = 0.0;
  bool has_absolute_tolerance = false;
  double absolute_tolerance = 0.0;
  const std::vector<ZplotGaussHermiteRule> *rules = nullptr;
  // The controller's ladders: c(3L, 7L, SELNORM_CLUSTER_QUADRATURE_ORDERS)
  // and c(7L, SELNORM_CLUSTER_QUADRATURE_ORDERS).
  std::vector<int> context_orders;
  std::vector<int> inner_orders;
};

// One draw: matrix arguments are column-major R matrices.
struct ZplotProjectionDraw {
  const double *mean = nullptr;            // k
  const double *covariance = nullptr;      // k x k
  const double *context_factor = nullptr;  // factor_rows x k
  int factor_rows = 1;
  const double *point_means = nullptr;     // point_rows x k, or nullptr
  int point_rows = 0;
  const double *omega = nullptr;           // bins
};

struct ZplotProjectionOutcome {
  enum Status { unavailable = 0, accepted = 1, failed = 2 };
  int status = unavailable;
  std::string message;
  // TRUE for the constant-weight Gaussian shortcut, whose list carries only
  // the absolute integration error and no quadrature orders.
  bool uniform = false;
  std::vector<double> density;
  double relative_error = 0.0;
  double mass_error = 0.0;
  double mass_omission_error = 0.0;
  // absolute, numerator, denominator, compression, cdf_curve, covariance,
  // pruning, context
  double integration_error[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  int quadrature[3] = {0, 0, 0};  // context, numerator, denominator
  bool bounded_cdf = true;
  // Refinement steps in order (zplot_trace_* codes), whatever the status.
  std::vector<int> trace;
};

enum ZplotTraceCode {
  trace_resolution_inner = 1, trace_resolution_context, trace_selection_inner,
  trace_selection_denominator, trace_selection_outer, trace_omission_allowance,
  trace_route_normalizer, trace_route_inner, trace_route_outer, trace_route_cdf,
  trace_route_compression, trace_route_pruning, trace_route_covariance, trace_own_normalizer
};

// .zplot_context_projection() for every draw of a block, after the R-side
// argument checks that do not depend on the draw. Never throws for a draw:
// failures are returned with their R message (only an interrupt unwinds).
// With 'threads' <= 1 each draw runs to completion on the main thread, which
// polls for interrupts between kernel calls; otherwise draws run in parallel
// rounds (kernel threads = 1) and stop only where they need a Gaussian
// resolution bound, which the main thread computes between rounds. Either way
// every draw performs the serial computation. Call from an R entry's native
// main-thread context.
void zplot_context_projection_batch(const ZplotProjectionSetup &setup,
  const std::vector<ZplotProjectionDraw> &draws, int threads, int chunk_rows,
  std::vector<ZplotProjectionOutcome> *outcomes);

#endif
