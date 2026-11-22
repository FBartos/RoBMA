#ifndef DSELNORMSAMPLINGCONDITIONED_H_
#define DSELNORMSAMPLINGCONDITIONED_H_

#include <distribution/VectorDist.h>

namespace jags {
namespace RoBMA {

class DSELNORMSAMPLINGCONDITIONED : public VectorDist {
public:
  DSELNORMSAMPLINGCONDITIONED();
  double logDensity(double const *x, PDFType type,
    std::vector<double const *> const &par, std::vector<unsigned long> const &len) const override;
  void randomSample(double *x,
    std::vector<double const *> const &par, std::vector<unsigned long> const &len, RNG *rng) const override;
  bool checkParameterValue(std::vector<double const *> const &par,
    std::vector<unsigned long> const &len) const override;
  bool checkParameterLength(std::vector<unsigned long> const &len) const override;
  unsigned long length(std::vector<unsigned long> const &len) const override;
  void support(double *lower, double *upper,
    std::vector<double const *> const &par, std::vector<unsigned long> const &len) const override;
  bool isSupportFixed(std::vector<bool> const &fixmask) const override;
};

}
}
#endif
