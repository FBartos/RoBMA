#ifndef DWMN2_H_
#define DWMN2_H_
#include <distribution/ArrayDist.h>

namespace jags {
namespace RoBMA { // module namespace

class DWMN2 : public ArrayDist {
public:
  DWMN2();

  double logDensity(double const *x, unsigned int length, PDFType type,
		    std::vector<double const *> const &parameters,
		    std::vector<std::vector<unsigned int> > const &dims,
		    double const *lower, double const *upper) const;
  void randomSample(double *x, unsigned int length,
		    std::vector<double const *> const &parameters,
		    std::vector<std::vector<unsigned int> > const &dims,
		    double const *lower, double const *upper, RNG *rng) const;
  void typicalValue(double *x, unsigned int length,
		    std::vector<double const *> const &parameters,
		    std::vector<std::vector<unsigned int> > const &dims,
		    double const *lower, double const *upper) const;
  bool checkParameterDim(std::vector<std::vector<unsigned int> > const &dims) const;
  bool checkParameterValue(std::vector<double const *> const &parameters,
			std::vector<std::vector<unsigned int> > const &dims) const;
  std::vector<unsigned int> dim(std::vector<std::vector<unsigned int> > const &dims) const;
  void support(double *lower, double *upper, unsigned int length,
	       std::vector<double const *> const &parameters,
	       std::vector<std::vector<unsigned int> > const &dims) const;
  bool isSupportFixed(std::vector<bool> const &fixmask) const;
};

}}
#endif /* DWMN2_H_ */


