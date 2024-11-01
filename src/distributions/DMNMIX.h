#ifndef DMNMIX_H_
#define DMNMIX_H_
#include <distribution/ArrayDist.h>

namespace jags {
namespace RoBMA { // module namespace

class DMNMIX : public ArrayDist {
public:
  DMNMIX();

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
  void simulate_mnorm(double *x, const double *mu, const double *sigma, int K, RNG *rng) const;
};

}}
#endif /* DMNMIX_H_ */


