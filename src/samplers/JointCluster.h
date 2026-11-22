#ifndef ROBMA_JOINT_CLUSTER_H_
#define ROBMA_JOINT_CLUSTER_H_

#include "JointSamplerTag.h"

namespace jags {

class SamplerFactory;

namespace RoBMA {

// The member and setup declarations are unobserved deterministic terminals.
class ClusterTag : public JointSamplerTag {
  const bool _setup;
public:
  explicit ClusterTag(bool setup);
  bool checkParameterLength(std::vector<unsigned long> const &) const override;
};

SamplerFactory *joint_cluster_factory();

} // namespace RoBMA
} // namespace jags

#endif
