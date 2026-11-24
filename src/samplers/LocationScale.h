#ifndef ROBMA_LOCATION_SCALE_H
#define ROBMA_LOCATION_SCALE_H

// Heterogeneity scale interweaving in the location block sampler
// (LocationBlock.h); the design is .work/tmp/fitopt/design/wp6-spec.md
// (revision 2) of the workspace.
//
// The location block owns the standardized latents Z_1..Z_K (dnorm(0, 1)) of
// the random-effect term it samples. With a the total heterogeneity scale, w
// the weights of the L variance components (v_j = a^2 w_j) and s the sampled
// component (sigma_s = a sqrt(w_s) multiplies the latents), the scale group
// T_c, c = exp(u) > 0, sets v_s to c^2 v_s, keeps the other v_j and sets Z to
// Z / c, so sigma_s Z and every other component SD are unchanged and no child
// density changes. With d = 1 - w_s + c^2 w_s:
//   a' = a sqrt(d), w_s' = c^2 w_s / d, w_j' = w_j / d (j != s), Z' = Z / c;
// gamma weight sources eta_j (w_j = eta_j / T) keep their total T and become
// eta_j' = T w_j'. u is drawn by one slice update from u = 0 (width 0.5, a
// total budget of 40 expansions split at random between the sides, Neal 2003,
// then shrinkage) of
//   log f(u) = log p_a(a') + log p_w(w') - (2L - 1) log a'
//              - sum(Z^2) / (2 c^2) + (2 - K) u,
// p the sources' own prior densities (PDF_PRIOR, fixed parameters), which by
// Liu and Sabatti (2000) leaves the posterior invariant. The move needs no
// likelihood evaluation, tuning or acceptance step. It runs after the location
// sweep and after the ordinary update of the heterogeneity sources, in states
// where every declared gate indicator is at its declared level; otherwise it
// is the identity.
//
// R declares the structure with one terminal tag anchored at the source of a:
//   robma_location_scale(a source, a, sigma_s, w nodes, share source,
//                        indicators, levels, bias indicator, c(L, s, K))
// w nodes: the L ratio nodes eta_j / sum(eta) (gamma sources, resolved from
// the ratios so that the tag adds no child to the gamma group), or the chart
// of the share w_s (with its scalar source on (0, 1), L = 2, s = 1), or 0
// (L = 1); indicators and levels: the gate and inclusion indicators and the
// level at which each makes the declared sources active, or 0; the bias
// indicator or 0. The factory verifies the declaration on the compiled graph,
// runs the invariance check once (chain 0, indicators forced to their levels,
// every bias label, c in {0.7, 1.3}: every parent value of every stochastic
// child of the sources and latents unchanged within 1e-9 (1 + |value|) and
// sigma_s' = c sigma_s; state restored) and builds the sources' ordinary
// samplers by asking the factories after the location block's, which the
// owner embeds. Anything else declines the scale ownership: the location block
// and every other node keep exactly what they had.

#include "JointSamplerTag.h"

#include <graph/Graph.h>
#include <graph/StochasticNode.h>
#include <sampler/Sampler.h>
#include <sampler/SamplerFactory.h>

#include <list>
#include <memory>
#include <vector>

namespace jags {
namespace RoBMA {

class LocationScaleTag : public JointSamplerTag {
public:
  LocationScaleTag();
  bool checkParameterLength(std::vector<unsigned long> const &lengths) const override;
};

// Called by the location block factory with the sampler it built: returns it
// unchanged without a scale declaration or when the declaration is declined,
// otherwise one sampler that owns the location members and the heterogeneity
// sources (its update: the location sweep, the embedded updates, the move).
std::unique_ptr<Sampler> location_scale_owner(std::unique_ptr<Sampler> location,
                                              std::list<StochasticNode *> const &free, Graph const &graph,
                                              std::vector<StochasticNode *> const &latents,
                                              SamplerFactory const *location_factory);

}
}
#endif
