#ifndef ROBMA_LOCATION_BLOCK_H
#define ROBMA_LOCATION_BLOCK_H

// The location block sampler (RoBMA::LocationBlock) of model-averaged `.mv`
// fits: one owner of the location gates, their Gaussian slabs (and ungated
// Gaussian location coefficients) and the standardized cluster latents of one
// random-effect term. One update is a sweep of component moves, each
// reversible for the surrogate target, in a fixed order or fully reversed,
// followed by one exact correction against the JAGS densities of the graph.
// The design is .work/tmp/fitopt/design/glb-spec-v2.md of the workspace.
//
// R declares the members through terminal tags (JointSamplerTag.h):
//   robma_location_group(gate, slab, label, physical, slot)   one per group
//   robma_location_latent(Z, slot)                            one per latent
//   robma_location_setup(first gate, groups, latents)         one per model
// An ungated scalar group passes its slab for the gate, slab and physical, and
// label 0. A slab is one dnorm or dmnorm node or an array of scalar dnorm nodes.

#include "JointSamplerTag.h"
#include "LocationBlockReport.h"

#include <sampler/SamplerFactory.h>

#include <vector>

namespace jags {
namespace RoBMA {

class LocationTag : public JointSamplerTag {
  const int _kind;
public:
  // kind 0 group, 1 latent, 2 setup
  explicit LocationTag(int kind);
  bool checkParameterLength(std::vector<unsigned long> const &lengths) const override;
};

// Borrowed process-lifetime factory. Module registration must not delete it.
SamplerFactory *location_block_factory();


}
}
#endif
