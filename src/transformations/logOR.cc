#include "logOR.h"
#include "../source/transformations.h"
#include <util/nainf.h>
#include <JRmath.h>
#include <cmath>



namespace jags {
  namespace RoBMA {

    // effect sizes transformations
    logOR2d::logOR2d() :ScalarFunction("logOR2d", 1)
    {}
    bool logOR2d::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double logOR2d::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_logOR2d(*args[0]);
    }

    logOR2z::logOR2z() :ScalarFunction("logOR2z", 1)
    {}
    bool logOR2z::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double logOR2z::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_logOR2z(*args[0]);
    }

    logOR2r::logOR2r() :ScalarFunction("logOR2r", 1)
    {}
    bool logOR2r::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double logOR2r::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_logOR2r(*args[0]);
    }

    // standard errors transformations
    se_logOR2se_d::se_logOR2se_d() :ScalarFunction("se_logOR2se_d", 2)
    {}
    bool se_logOR2se_d::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(*args[0] >= 0);
    }
    double se_logOR2se_d::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_se_logOR2se_d(*args[0]);
    }

    se_logOR2se_z::se_logOR2se_z() :ScalarFunction("se_logOR2se_z", 2)
    {}
    bool se_logOR2se_z::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(*args[0] >= 0);
    }
    double se_logOR2se_z::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_se_logOR2se_z(*args[0], *args[1]);
    }

    se_logOR2se_r::se_logOR2se_r() :ScalarFunction("se_logOR2se_r", 2)
    {}
    bool se_logOR2se_r::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(*args[0] >= 0);
    }
    double se_logOR2se_r::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_se_logOR2se_r(*args[0], *args[1]);
    }

    // linear scaling function
    scale_logOR2d::scale_logOR2d() :ScalarFunction("scale_logOR2d", 1)
    {}
    bool scale_logOR2d::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double scale_logOR2d::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_scale_logOR2d(*args[0]);
    }

    scale_logOR2z::scale_logOR2z() :ScalarFunction("scale_logOR2z", 1)
    {}
    bool scale_logOR2z::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double scale_logOR2z::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_scale_logOR2z(*args[0]);
    }

    scale_logOR2r::scale_logOR2r() :ScalarFunction("scale_logOR2r", 1)
    {}
    bool scale_logOR2r::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double scale_logOR2r::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_scale_logOR2r(*args[0]);
    }
  }
}
