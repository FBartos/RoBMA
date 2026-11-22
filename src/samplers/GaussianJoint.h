#ifndef ROBMA_GAUSSIAN_JOINT_H_
#define ROBMA_GAUSSIAN_JOINT_H_

#include <distribution/VectorDist.h>

namespace jags {

class SamplerFactory;

namespace RoBMA {

// Observed point-mass marker. The six scalar stochastic parents are paired
// indicator/slab leaves in intercept, moderator, heterogeneity order. The
// remaining parents are the three active labels, moderator column, observed
// outcomes, sampling SEs, and the effect direction (+1 or -1).
class GaussianJointMarker : public VectorDist {
public:
  GaussianJointMarker();
  double logDensity(const double *x, PDFType,
                    const std::vector<double const *> &,
                    const std::vector<unsigned long> &) const override;
  void randomSample(double *x, const std::vector<double const *> &,
                    const std::vector<unsigned long> &, RNG *) const override;
  void support(double *lower, double *upper, const std::vector<double const *> &,
               const std::vector<unsigned long> &) const override;
  bool checkParameterLength(const std::vector<unsigned long> &lengths) const override;
  bool checkParameterValue(const std::vector<double const *> &,
                           const std::vector<unsigned long> &) const override;
  unsigned long length(const std::vector<unsigned long> &) const override;
  bool isSupportFixed(const std::vector<bool> &) const override;
  bool isDiscreteValued(const std::vector<bool> &) const override;
};

SamplerFactory *gaussian_joint_factory();

} // namespace RoBMA
} // namespace jags

#endif
