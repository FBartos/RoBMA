#ifndef DSELNORMCLUSTERSTEP_H_
#define DSELNORMCLUSTERSTEP_H_

#include <distribution/VectorDist.h>

struct SelNormCoarseSettings;

namespace jags {
namespace RoBMA {

class DSELNORMCLUSTERSTEP : public VectorDist
  {
  public:
    DSELNORMCLUSTERSTEP();

  double logDensity(double const *x, PDFType type,
                    std::vector<double const *> const &parameters,
                    std::vector<unsigned long> const &lengths) const override;
  // The corrected sampler's surrogate density: the exact numerator minus the
  // log normalizer of one rule of the block's own schedule (rule number
  // settings.max_rules or the named settings.rule, or the last) at the actual
  // parameters, and the exact density where that rule is unavailable. Not
  // memoized (see logDensity).
  // surrogate_not_exact, when given, is set to true (never cleared) when the
  // value returned was computed with a normalizer that is not the exact one.
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
  // Dynamic: kernel mode, vector rule, mean, residual SD, loadings and weights.
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
