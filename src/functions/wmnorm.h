#ifndef wmnorm_H_
#define wmnorm_H_

#include <function/ArrayFunction.h>

namespace jags {
  namespace RoBMA {

    class wnorm_1s_lpdf : public ArrayFunction
    {
    public:
      wnorm_1s_lpdf();

      void evaluate(double *value, std::vector<double const *> const &args, std::vector<std::vector<unsigned int> > const &dims) const;
      bool checkParameterDim (std::vector<std::vector<unsigned int> > const &dims) const;
      bool checkParameterValue(std::vector<double const *> const &par, std::vector<std::vector<unsigned int> > const &dims) const;
      std::vector<unsigned int> dim(std::vector<std::vector<unsigned int> > const &dims, std::vector<double const *> const &values) const;
    };

    class wnorm_2s_lpdf : public ArrayFunction
    {
    public:
      wnorm_2s_lpdf();

      void evaluate(double *value, std::vector<double const *> const &args, std::vector<std::vector<unsigned int> > const &dims) const;
      bool checkParameterDim (std::vector<std::vector<unsigned int> > const &dims) const;
      bool checkParameterValue(std::vector<double const *> const &par, std::vector<std::vector<unsigned int> > const &dims) const;
      std::vector<unsigned int> dim(std::vector<std::vector<unsigned int> > const &dims, std::vector<double const *> const &values) const;
    };

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

    class wmnorm_1s_v_lpdf : public ArrayFunction
    {
      public:
        wmnorm_1s_v_lpdf();

      void evaluate(double *value, std::vector<double const *> const &args, std::vector<std::vector<unsigned int> > const &dims) const;
      bool checkParameterDim (std::vector<std::vector<unsigned int> > const &dims) const;
      bool checkParameterValue(std::vector<double const *> const &par, std::vector<std::vector<unsigned int> > const &dims) const;
      std::vector<unsigned int> dim(std::vector<std::vector<unsigned int> > const &dims, std::vector<double const *> const &values) const;
    };

    class wmnorm_2s_v_lpdf : public ArrayFunction
    {
      public:
        wmnorm_2s_v_lpdf();

      void evaluate(double *value, std::vector<double const *> const &args, std::vector<std::vector<unsigned int> > const &dims) const;
      bool checkParameterDim (std::vector<std::vector<unsigned int> > const &dims) const;
      bool checkParameterValue(std::vector<double const *> const &par, std::vector<std::vector<unsigned int> > const &dims) const;
      std::vector<unsigned int> dim(std::vector<std::vector<unsigned int> > const &dims, std::vector<double const *> const &values) const;
    };
  }
}

#endif /* wmnorm_H_ */
