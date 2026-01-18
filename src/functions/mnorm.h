#ifndef mnorm_H_
#define mnorm_H_

#include <function/ArrayFunction.h>

namespace jags {
  namespace RoBMA {

    class mnorm_lpdf : public ArrayFunction
    {
      public:
        mnorm_lpdf();

		void evaluate(double *value, std::vector<double const *> const &args, std::vector<std::vector<unsigned long> > const &dims) const override;
		bool checkParameterDim (std::vector<std::vector<unsigned long> > const &dims) const override;
		bool checkParameterValue(std::vector<double const *> const &par, std::vector<std::vector<unsigned long> > const &dims) const override;
		std::vector<unsigned long> dim(std::vector<std::vector<unsigned long> > const &dims, std::vector<double const *> const &values) const override;
    };

    class mnorm_v_lpdf : public ArrayFunction
    {
      public:
        mnorm_v_lpdf();

		void evaluate(double *value, std::vector<double const *> const &args, std::vector<std::vector<unsigned long> > const &dims) const override;
		bool checkParameterDim (std::vector<std::vector<unsigned long> > const &dims) const override;
		bool checkParameterValue(std::vector<double const *> const &par, std::vector<std::vector<unsigned long> > const &dims) const override;
		std::vector<unsigned long> dim(std::vector<std::vector<unsigned long> > const &dims, std::vector<double const *> const &values) const override;
    };
  }
}

#endif /* mnorm_H_ */
