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

  double logDensity(double const *x, PDFType type,
                    std::vector<double const *> const &parameters,
                    std::vector<unsigned long> const &lengths) const override;
  double surrogateLogDensity(double const *x, PDFType type,
                    std::vector<double const *> const &parameters,
                    std::vector<unsigned long> const &lengths,
                    SelNormCoarseSettings const &settings,
                    bool *surrogate_not_exact = nullptr) const;
  void randomSample(double *x,
                    std::vector<double const *> const &parameters,
                    std::vector<unsigned long> const &lengths, RNG *rng) const override;
  bool checkParameterValue(std::vector<double const *> const &parameters,
                           std::vector<unsigned long> const &lengths) const override;
  // Lengths are fixed by the compiled node. The full validator always calls
  // both helpers; only a caller with a fixed-control certificate may split it.
  bool checkDynamicParameterValue(std::vector<double const *> const &parameters,
                                  std::vector<unsigned long> const &lengths) const;
  bool checkControlParameterValue(std::vector<double const *> const &parameters,
                                  std::vector<unsigned long> const &lengths) const;
  bool checkParameterLength(std::vector<unsigned long> const &lengths) const override;
  unsigned long length(std::vector<unsigned long> const &lengths) const override;
  void support(double *lower, double *upper,
               std::vector<double const *> const &parameters,
               std::vector<unsigned long> const &lengths) const override;
  bool isSupportFixed(std::vector<bool> const &fixmask) const override;
  };

}
}

#endif
