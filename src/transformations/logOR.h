#ifndef logOR_H_
#define logOR_H_

#include <function/ScalarFunction.h>

namespace jags {
  namespace RoBMA {

    // effect sizes transformations
    class logOR2d : public ScalarFunction 
    {
      public:
        logOR2d();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class logOR2r : public ScalarFunction 
    {
      public:
        logOR2r();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class logOR2z : public ScalarFunction 
    {
      public:
        logOR2z();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    // standard errors transformations
    class se_logOR2se_d : public ScalarFunction 
    {
      public:
        se_logOR2se_d();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class se_logOR2se_r : public ScalarFunction 
    {
      public:
        se_logOR2se_r();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class se_logOR2se_z : public ScalarFunction 
    {
      public:
        se_logOR2se_z();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    // linear scaling function
    class scale_logOR2d : public ScalarFunction 
    {
      public:
        scale_logOR2d();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class scale_logOR2z : public ScalarFunction 
    {
      public:
        scale_logOR2z();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

    class scale_logOR2r : public ScalarFunction 
    {
      public:
        scale_logOR2r();

        bool checkParameterValue(std::vector<double const *> const &args) const;
        double evaluate(std::vector<double const *> const &args) const;
    };

  }
}

#endif /* logOR_H_ */
