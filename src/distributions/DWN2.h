#ifndef DWN2_H_
#define DWN2_H_
#include <distribution/VectorDist.h>

namespace jags {
namespace RoBMA { // module namespace

class DWN2 : public VectorDist
  {
  public:
    DWN2();
    
  double logDensity(double const *x,  PDFType type, 
		    std::vector<double const *> const &parameters,
		    std::vector<unsigned long> const &lengths) const override;
  void randomSample(double *x, 
		    std::vector<double const *> const &parameters,
		    std::vector<unsigned long> const &lengths, RNG *rng) const override;
  bool checkParameterValue(std::vector<double const *> const &parameters,
                           std::vector<unsigned long> const &lengths) const override;
  bool checkParameterLength(std::vector<unsigned long> const &lengths) const override;
  unsigned long length(std::vector<unsigned long> const &dim) const override;
  void support(double *lower, double *upper, 
	       std::vector<double const *> const &parameters,
	       std::vector<unsigned long> const &lengths) const override;
  bool isSupportFixed(std::vector<bool> const &fixmask) const override;
  };

}
}
#endif /* DWN2_H_ */


