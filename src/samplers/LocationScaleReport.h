#ifndef ROBMA_LOCATION_SCALE_REPORT_H
#define ROBMA_LOCATION_SCALE_REPORT_H

// The per-fit report and the test hooks of the location block's heterogeneity
// scale ownership (LocationScale.h), without JAGS headers, for the R-native
// entry points (r-location-scale.cc).

#include <string>
#include <vector>

namespace jags {
namespace RoBMA {

// The per-fit counters of this process (relaxed atomics), read and cleared
// when clear is true, in the order of location_scale_report_names(): scale
// owners built and declarations declined; scale updates in a state where the
// declared sources are active (moves, slice calls that exhausted their
// interval included), the identity updates of the states where they are not,
// the identity updates of states outside the move's set (a weight or share
// beyond [2^-20, 1 - 2^-20], a source outside its support), the updates
// skipped because log f(0) was not finite and the slices exhausted in
// floating point; the density evaluations of the moves and the sum of the
// squared steps u = log c; with the timing switch the wall time of the
// embedded heterogeneity updates, of the scale moves and, within them, of the
// per-move check that no child's parent value changed.
std::vector<char const *> const &location_scale_report_names();
std::vector<double> location_scale_report(bool clear);

// One factory call that found a scale declaration: the layout ("allocation"
// for gamma weight sources, "clustered" for a scalar share, "single" without
// weights), L and K as declared, the outcome ("owned" or "declined"), the
// reason of a decline, the outcome of the invariance check ("passed",
// "failed: ..." or "not run") and the class names of the embedded samplers.
struct LocationScaleDecision {
  std::string layout, status, reason, check, embedded;
  int L = 0, K = 0;
};
// The decisions of this process since the last clearing read.
std::vector<LocationScaleDecision> location_scale_decisions(bool clear);

// Test hooks, not package interface. Switches frozen into the scale owners
// built from now on (a negative argument keeps a setting): mutant 0 none,
// 1 the exponent c^(1 - K) instead of c^(2 - K), 2 the factor a'^-(2L - 1)
// omitted, 3 the last latent left out of the transformation, 4 the gamma
// total not preserved (only eta_s scaled, by c^2); decline 1 makes the factory
// decline every scale declaration after reading it; timing 1 accumulates the
// wall times. Returns the previous three.
std::vector<int> location_scale_test_config(int mutant, int decline, int timing);
// The serial number of the most recently built live scale owner (owners are
// numbered from 0 in construction order), or -1 without one. The functions
// below take this serial number as `index`.
int location_scale_test_last();

// The state of one chain of a live scale owner: the layout, L, the sampled
// component s (1-based) and K; whether the declared sources are active; the
// move's coordinates in the order latents, a, weight sources; the weights w;
// the value of the sigma_s node; the log densities (PDF_FULL) of the
// stochastic children of the sources and latents (markers excluded); the sum
// of the sources' prior log densities from the nodes' own logDensity
// (PDF_PRIOR) and from the move's evaluation at u = 0.
struct LocationScaleState {
  std::string layout;
  int L = 0, s = 0, K = 0;
  bool active = false;
  std::vector<double> values, weights, children;
  double sigma = 0, prior_nodes = 0, prior_move = 0;
};
LocationScaleState location_scale_test_state(int index, unsigned int chain);
// Writes the move's coordinates (the order of the state) and updates their
// deterministic descendants.
void location_scale_test_set(int index, unsigned int chain, std::vector<double> const &values);
// Applies T_c with c = exp(u) (the owner's mutant included) with the move's
// check of every child's parent values (the classed error, entry restored);
// false when the current or the transformed state is outside the move's set
// (nothing is written).
bool location_scale_test_transform(int index, unsigned int chain, double u);
// log f(u) of the move at the chain's current state (the owner's mutant
// included); -Inf where the transformed state is outside the move's set.
std::vector<double> location_scale_test_density(int index, unsigned int chain, std::vector<double> const &u);
// One slice update of the move's kernel (width 0.5, budget 40) from each
// start x0 of a uniform target on (lower, upper), shifted so that the update
// starts at 0 as in the move; with fixed_cap the budget is split 20/20
// instead of at random (an invalid kernel). Returns the new points.
std::vector<double> location_scale_test_slice(double lower, double upper, std::vector<double> const &starts,
                                              bool fixed_cap, double (*uniform)(), double (*exponential)());

}
}
#endif
