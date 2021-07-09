#ifndef r_H_
#define r_H_

#include <function/ScalarFunction.h>

namespace jags {
  namespace RoBMA {

    // effect sizes transformations
    class r2d : public ScalarFunction
    {
      public:
        r2d();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class r2z : public ScalarFunction
    {
      public:
        r2z();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class r2logOR : public ScalarFunction
    {
      public:
        r2logOR();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    // standard errors transformations
    class se_r2se_d : public ScalarFunction
    {
      public:
        se_r2se_d();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class se_r2se_z : public ScalarFunction
    {
      public:
        se_r2se_z();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class se_r2se_logOR : public ScalarFunction
    {
      public:
        se_r2se_logOR();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

  }
}

#endif /* r_H_ */
