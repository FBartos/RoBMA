#ifndef wmnorm_H_
#define wmnorm_H_

#include <function/ArrayFunction.h>

namespace jags {
  namespace RoBMA {

    class wnorm_1s_lpdf : public ArrayFunction
    {
    public:
      wnorm_1s_lpdf();

      void evaluate(double *value, std::vector<double const *> const &args, std::vector<std::vector<unsigned long> > const &dims) const override;
      bool checkParameterDim (std::vector<std::vector<unsigned long> > const &dims) const override;
      bool checkParameterValue(std::vector<double const *> const &par, std::vector<std::vector<unsigned long> > const &dims) const override;
      std::vector<unsigned long> dim(std::vector<std::vector<unsigned long> > const &dims, std::vector<double const *> const &values) const override;
    };

    class wnorm_2s_lpdf : public ArrayFunction
    {
    public:
      wnorm_2s_lpdf();

      void evaluate(double *value, std::vector<double const *> const &args, std::vector<std::vector<unsigned long> > const &dims) const override;
      bool checkParameterDim (std::vector<std::vector<unsigned long> > const &dims) const override;
      bool checkParameterValue(std::vector<double const *> const &par, std::vector<std::vector<unsigned long> > const &dims) const override;
      std::vector<unsigned long> dim(std::vector<std::vector<unsigned long> > const &dims, std::vector<double const *> const &values) const override;
    };

    class wmnorm_1s_lpdf : public ArrayFunction
    {
      public:
        wmnorm_1s_lpdf();

      void evaluate(double *value, std::vector<double const *> const &args, std::vector<std::vector<unsigned long> > const &dims) const override;
      bool checkParameterDim (std::vector<std::vector<unsigned long> > const &dims) const override;
      bool checkParameterValue(std::vector<double const *> const &par, std::vector<std::vector<unsigned long> > const &dims) const override;
      std::vector<unsigned long> dim(std::vector<std::vector<unsigned long> > const &dims, std::vector<double const *> const &values) const override;
    };

    class wmnorm_2s_lpdf : public ArrayFunction
    {
      public:
        wmnorm_2s_lpdf();

      void evaluate(double *value, std::vector<double const *> const &args, std::vector<std::vector<unsigned long> > const &dims) const override;
      bool checkParameterDim (std::vector<std::vector<unsigned long> > const &dims) const override;
      bool checkParameterValue(std::vector<double const *> const &par, std::vector<std::vector<unsigned long> > const &dims) const override;
      std::vector<unsigned long> dim(std::vector<std::vector<unsigned long> > const &dims, std::vector<double const *> const &values) const override;
    };

    class wmnorm_1s_v_lpdf : public ArrayFunction
    {
      public:
        wmnorm_1s_v_lpdf();

      void evaluate(double *value, std::vector<double const *> const &args, std::vector<std::vector<unsigned long> > const &dims) const override;
      bool checkParameterDim (std::vector<std::vector<unsigned long> > const &dims) const override;
      bool checkParameterValue(std::vector<double const *> const &par, std::vector<std::vector<unsigned long> > const &dims) const override;
      std::vector<unsigned long> dim(std::vector<std::vector<unsigned long> > const &dims, std::vector<double const *> const &values) const override;
    };

    class wmnorm_2s_v_lpdf : public ArrayFunction
    {
      public:
        wmnorm_2s_v_lpdf();

      void evaluate(double *value, std::vector<double const *> const &args, std::vector<std::vector<unsigned long> > const &dims) const override;
      bool checkParameterDim (std::vector<std::vector<unsigned long> > const &dims) const override;
      bool checkParameterValue(std::vector<double const *> const &par, std::vector<std::vector<unsigned long> > const &dims) const override;
      std::vector<unsigned long> dim(std::vector<std::vector<unsigned long> > const &dims, std::vector<double const *> const &values) const override;
    };
  }
}

#endif /* wmnorm_H_ */
