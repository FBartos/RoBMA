#ifndef DSELNORMSAMPLINGCONDITIONED_H_
#define DSELNORMSAMPLINGCONDITIONED_H_

#include <distribution/VectorDist.h>

namespace jags {
namespace RoBMA {

class DSELNORMSAMPLINGCONDITIONED : public VectorDist {
public:
  DSELNORMSAMPLINGCONDITIONED();
  double logDensity(double const *x, unsigned int length, PDFType type,
    std::vector<double const *> const &par, std::vector<unsigned int> const &len,
    double const *lower, double const *upper) const;
  void randomSample(double *x, unsigned int length,
    std::vector<double const *> const &par, std::vector<unsigned int> const &len,
    double const *lower, double const *upper, RNG *rng) const;
  void typicalValue(double *x, unsigned int length,
    std::vector<double const *> const &par, std::vector<unsigned int> const &len,
    double const *lower, double const *upper) const;
  bool checkParameterValue(std::vector<double const *> const &par,
    std::vector<unsigned int> const &len) const;
  bool checkParameterLength(std::vector<unsigned int> const &len) const;
  unsigned int length(std::vector<unsigned int> const &len) const;
  void support(double *lower, double *upper, unsigned int length,
    std::vector<double const *> const &par, std::vector<unsigned int> const &len) const;
  bool isSupportFixed(std::vector<bool> const &fixmask) const;
};

}
}
#endif
