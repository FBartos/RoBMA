#ifndef wmnorm_H_
#define wmnorm_H_

#include <function/ArrayFunction.h>

double cpp_wmnorm_1s_lpdf(double const *x, double const *mu, double const *sigma, double const *omega, double const *crit_x, const int K, const int J);
double cpp_wmnorm_2s_lpdf(double const *x, double const *mu, double const *sigma, double const *omega, double const *crit_x, const int K, const int J);


namespace jags {
  namespace RoBMA {

    class wmnorm_1s_lpdf : public ArrayFunction
    {
      public:
        wmnorm_1s_lpdf();

    void evaluate(double *value, std::vector<double const *> const &args, std::vector<std::vector<unsigned int> > const &dims) const;
		bool checkParameterDim (std::vector<std::vector<unsigned int> > const &dims) const;
		bool checkParameterValue(std::vector<double const *> const &par, std::vector<std::vector<unsigned int> > const &dims) const;
		std::vector<unsigned int> dim(std::vector<std::vector<unsigned int> > const &dims, std::vector<double const *> const &values) const;
    };

    class wmnorm_2s_lpdf : public ArrayFunction
    {
      public:
        wmnorm_2s_lpdf();

    void evaluate(double *value, std::vector<double const *> const &args, std::vector<std::vector<unsigned int> > const &dims) const;
		bool checkParameterDim (std::vector<std::vector<unsigned int> > const &dims) const;
		bool checkParameterValue(std::vector<double const *> const &par, std::vector<std::vector<unsigned int> > const &dims) const;
		std::vector<unsigned int> dim(std::vector<std::vector<unsigned int> > const &dims, std::vector<double const *> const &values) const;
    };

  }
}

#endif /* wmnorm_H_ */
