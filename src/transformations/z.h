#ifndef z_H_
#define z_H_

#include <function/ScalarFunction.h>

namespace jags {
  namespace RoBMA {

    // effect sizes transformations
    class z2d : public ScalarFunction 
    {
      public:
        z2d();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class z2r : public ScalarFunction 
    {
      public:
        z2r();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class z2logOR : public ScalarFunction 
    {
      public:
        z2logOR();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    // standard errors transformations
    class se_z2se_d : public ScalarFunction 
    {
      public:
        se_z2se_d();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class se_z2se_r : public ScalarFunction 
    {
      public:
        se_z2se_r();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class se_z2se_logOR : public ScalarFunction 
    {
      public:
        se_z2se_logOR();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    // linear scaling function
    class scale_z2d : public ScalarFunction 
    {
      public:
        scale_z2d();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class scale_z2logOR : public ScalarFunction 
    {
      public:
        scale_z2logOR();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

  }
}

#endif /* z_H_ */
