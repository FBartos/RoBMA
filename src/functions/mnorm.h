#ifndef mnorm_H_
#define mnorm_H_

#include <function/ArrayFunction.h>

namespace jags {
  namespace RoBMA {

    class mnorm_lpdf : public ArrayFunction
    {
      public:
        mnorm_lpdf();

		void evaluate(double *value, std::vector<double const *> const &args, std::vector<std::vector<unsigned int> > const &dims) const;
		bool checkParameterDim (std::vector<std::vector<unsigned int> > const &dims) const;
		bool checkParameterValue(std::vector<double const *> const &par, std::vector<std::vector<unsigned int> > const &dims) const;
		std::vector<unsigned int> dim(std::vector<std::vector<unsigned int> > const &dims, std::vector<double const *> const &values) const;
    };

    class mnorm_v_lpdf : public ArrayFunction
    {
      public:
        mnorm_v_lpdf();

		void evaluate(double *value, std::vector<double const *> const &args, std::vector<std::vector<unsigned int> > const &dims) const;
		bool checkParameterDim (std::vector<std::vector<unsigned int> > const &dims) const;
		bool checkParameterValue(std::vector<double const *> const &par, std::vector<std::vector<unsigned int> > const &dims) const;
		std::vector<unsigned int> dim(std::vector<std::vector<unsigned int> > const &dims, std::vector<double const *> const &values) const;
    };
  }
}

#endif /* mnorm_H_ */
