#ifndef DWMN1_H_
#define DWMN1_H_
#include <distribution/ArrayDist.h>

namespace jags {
namespace RoBMA { // module namespace

class DWMN1 : public ArrayDist {
public:
  DWMN1();

  double logDensity(double const *x,  PDFType type,
		    std::vector<double const *> const &parameters,
		    std::vector<std::vector<unsigned long> > const &dims) const override;
  void randomSample(double *x, 
		    std::vector<double const *> const &parameters,
		    std::vector<std::vector<unsigned long> > const &dims, RNG *rng) const override;
  bool checkParameterDim(std::vector<std::vector<unsigned long> > const &dims) const override;
  bool checkParameterValue(std::vector<double const *> const &parameters,
			std::vector<std::vector<unsigned long> > const &dims) const override;
  std::vector<unsigned long> dim(std::vector<std::vector<unsigned long> > const &dims) const override;
  void support(double *lower, double *upper, 
	       std::vector<double const *> const &parameters,
	       std::vector<std::vector<unsigned long> > const &dims) const override;
  bool isSupportFixed(std::vector<bool> const &fixmask) const override;
};

}}
#endif /* DWMN1_H_ */


