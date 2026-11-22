#ifndef ROBMA_JOINT_CLUSTER_H_
#define ROBMA_JOINT_CLUSTER_H_

#include <distribution/VectorDist.h>

namespace jags {

class SamplerFactory;

namespace RoBMA {

// The member and setup declarations are observed point masses at zero.
class ClusterMarker : public VectorDist {
  const bool _setup;
public:
  explicit ClusterMarker(bool setup);
  double logDensity(double const *, PDFType, std::vector<double const *> const &,
                    std::vector<unsigned long> const &) const override;
  void randomSample(double *, std::vector<double const *> const &,
                    std::vector<unsigned long> const &, RNG *) const override;
  void support(double *, double *, std::vector<double const *> const &,
               std::vector<unsigned long> const &) const override;
  bool checkParameterLength(std::vector<unsigned long> const &) const override;
  bool checkParameterValue(std::vector<double const *> const &,
                           std::vector<unsigned long> const &) const override;
  unsigned long length(std::vector<unsigned long> const &) const override;
  bool isSupportFixed(std::vector<bool> const &) const override;
  bool isDiscreteValued(std::vector<bool> const &) const override;
};

SamplerFactory *joint_cluster_factory();

} // namespace RoBMA
} // namespace jags

#endif
