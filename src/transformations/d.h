#ifndef d_H_
#define d_H_

#include <function/ScalarFunction.h>

namespace jags {
  namespace RoBMA {

    // effect sizes transformations
    class d2r : public ScalarFunction 
    {
      public:
        d2r();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class d2z : public ScalarFunction 
    {
      public:
        d2z();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class d2logOR : public ScalarFunction 
    {
      public:
        d2logOR();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    // standard errors transformations
    class se_d2se_r : public ScalarFunction 
    {
      public:
        se_d2se_r();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class se_d2se_z : public ScalarFunction 
    {
      public:
        se_d2se_z();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class se_d2se_logOR : public ScalarFunction 
    {
      public:
        se_d2se_logOR();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    // linear scaling function
    class scale_d2logOR : public ScalarFunction 
    {
      public:
        scale_d2logOR();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class scale_d2z : public ScalarFunction 
    {
      public:
        scale_d2z();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };
  }
}

#endif /* d_H_ */
