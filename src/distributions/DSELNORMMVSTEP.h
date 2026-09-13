#ifndef DSELNORMMVSTEP_H_
#define DSELNORMMVSTEP_H_

#include <distribution/VectorDist.h>

struct SelNormCoarseSettings;

namespace jags {
namespace RoBMA {

class DSELNORMMVSTEP : public VectorDist
  {
  public:
    DSELNORMMVSTEP();

  double logDensity(double const *x, unsigned int length, PDFType type,
                    std::vector<double const *> const &parameters,
                    std::vector<unsigned int> const &lengths,
                    double const *lower, double const *upper) const;
  double surrogateLogDensity(double const *x, unsigned int length, PDFType type,
                    std::vector<double const *> const &parameters,
                    std::vector<unsigned int> const &lengths,
                    double const *lower, double const *upper,
                    SelNormCoarseSettings const &settings) const;
  void randomSample(double *x, unsigned int length,
                    std::vector<double const *> const &parameters,
                    std::vector<unsigned int> const &lengths,
                    double const *lower, double const *upper, RNG *rng) const;
  void typicalValue(double *x, unsigned int length,
                    std::vector<double const *> const &parameters,
                    std::vector<unsigned int> const &lengths,
                    double const *lower, double const *upper) const;
  bool checkParameterValue(std::vector<double const *> const &parameters,
                           std::vector<unsigned int> const &lengths) const;
  // Lengths are fixed by the compiled node. The full validator always calls
  // both helpers; only a caller with a fixed-control certificate may split it.
  bool checkDynamicParameterValue(std::vector<double const *> const &parameters,
                                  std::vector<unsigned int> const &lengths) const;
  bool checkControlParameterValue(std::vector<double const *> const &parameters,
                                  std::vector<unsigned int> const &lengths) const;
  bool checkParameterLength(std::vector<unsigned int> const &lengths) const;
  unsigned int length(std::vector<unsigned int> const &lengths) const;
  void support(double *lower, double *upper, unsigned int length,
               std::vector<double const *> const &parameters,
               std::vector<unsigned int> const &lengths) const;
  bool isSupportFixed(std::vector<bool> const &fixmask) const;
  };

}
}

#endif
