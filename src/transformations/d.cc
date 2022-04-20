#include "d.h"
#include "../source/transformations.h"
#include <util/nainf.h>
#include <JRmath.h>
#include <cmath>



namespace jags {
  namespace RoBMA {

    // effect sizes transformations
    d2r::d2r() :ScalarFunction("d2r", 1)
    {}
    bool d2r::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double d2r::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_d2r(*args[0]);
    }

    d2z::d2z() :ScalarFunction("d2z", 1)
    {}
    bool d2z::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double d2z::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_d2z(*args[0]);
    }

    d2logOR::d2logOR() :ScalarFunction("d2logOR", 1)
    {}
    bool d2logOR::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double d2logOR::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_d2logOR(*args[0]);
    }

    // standard errors transformations
    se_d2se_r::se_d2se_r() :ScalarFunction("se_d2se_r", 2)
    {}
    bool se_d2se_r::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(*args[0] >= 0);
    }
    double se_d2se_r::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_se_d2se_r(*args[0], *args[1]);
    }

    se_d2se_z::se_d2se_z() :ScalarFunction("se_d2se_z", 2)
    {}
    bool se_d2se_z::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(*args[0] >= 0 && cpp_n_d(*args[1], *args[0]) > 3);
    }
    double se_d2se_z::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_se_d2se_z(*args[0], *args[1]);
    }

    se_d2se_logOR::se_d2se_logOR() :ScalarFunction("se_d2se_logOR", 2)
    {}
    bool se_d2se_logOR::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(*args[0] >= 0);
    }
    double se_d2se_logOR::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_se_d2se_logOR(*args[0]);
    }

    // linear scaling function
    scale_d2z::scale_d2z() :ScalarFunction("scale_d2z", 1)
    {}
    bool scale_d2z::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double scale_d2z::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_scale_d2z(*args[0]);
    }

    scale_d2logOR::scale_d2logOR() :ScalarFunction("scale_d2logOR", 1)
    {}
    bool scale_d2logOR::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double scale_d2logOR::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_scale_d2logOR(*args[0]);
    }

    scale_d2r::scale_d2r() :ScalarFunction("scale_d2r", 1)
    {}
    bool scale_d2r::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double scale_d2r::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_scale_d2r(*args[0]);
    }
  }
}
