#ifndef ROBMA_INDEPENDENT_BIAS_JOINT_H
#define ROBMA_INDEPENDENT_BIAS_JOINT_H

#include <distribution/VectorDist.h>
#include <sampler/SamplerFactory.h>

namespace jags::RoBMA {

// Observed zero with constant likelihood. The parents are the two indicators,
// raw location, slab label, location sign, branch families, branch-specific
// raw slopes, and branch-major weights on the common selection partition.
class DIndependentBiasJointMarker : public VectorDist {
public:
  DIndependentBiasJointMarker();
  double logDensity(double const *, PDFType, std::vector<double const *> const &,
                    std::vector<unsigned long> const &) const override;
  void randomSample(double *, std::vector<double const *> const &,
                    std::vector<unsigned long> const &, RNG *) const override;
  void support(double *, double *, std::vector<double const *> const &,
               std::vector<unsigned long> const &) const override;
  bool checkParameterValue(std::vector<double const *> const &,
                           std::vector<unsigned long> const &) const override;
  bool checkParameterLength(std::vector<unsigned long> const &) const override;
  bool checkParameterDiscrete(std::vector<bool> const &) const override;
  bool isDiscreteValued(std::vector<bool> const &) const override;
  bool isSupportFixed(std::vector<bool> const &) const override;
  unsigned long length(std::vector<unsigned long> const &) const override;
};

class IndependentBiasJointFactory : public SamplerFactory {
public:
  std::vector<Sampler *> makeSamplers(std::list<StochasticNode *> const &nodes,
                                     Graph const &graph) const override;
  std::string name() const override;
};

} // namespace jags::RoBMA

#endif
