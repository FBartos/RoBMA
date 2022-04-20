#include "z.h"
#include "../source/transformations.h"
#include <util/nainf.h>
#include <JRmath.h>
#include <cmath>



namespace jags {
  namespace RoBMA {

    // effect sizes transformations
    z2d::z2d() :ScalarFunction("z2d", 1)
    {}
    bool z2d::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double z2d::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_z2d(*args[0]);
    }

    z2r::z2r() :ScalarFunction("z2r", 1)
    {}
    bool z2r::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double z2r::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_z2r(*args[0]);
    }

    z2logOR::z2logOR() :ScalarFunction("z2logOR", 1)
    {}
    bool z2logOR::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double z2logOR::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_z2logOR(*args[0]);
    }

    // standard errors transformations
    se_z2se_d::se_z2se_d() :ScalarFunction("se_z2se_d", 2)
    {}
    bool se_z2se_d::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(*args[0] >= 0);
    }
    double se_z2se_d::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_se_z2se_d(*args[0], *args[1]);
    }

    se_z2se_r::se_z2se_r() :ScalarFunction("se_z2se_r", 2)
    {}
    bool se_z2se_r::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(*args[0] >= 0);
    }
    double se_z2se_r::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_se_z2se_r(*args[0], *args[1]);
    }

    se_z2se_logOR::se_z2se_logOR() :ScalarFunction("se_z2se_logOR", 2)
    {}
    bool se_z2se_logOR::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(*args[0] >= 0);
    }
    double se_z2se_logOR::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_se_z2se_logOR(*args[0], *args[1]);
    }

    // linear scaling function
    scale_z2d::scale_z2d() :ScalarFunction("scale_z2d", 1)
    {}
    bool scale_z2d::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double scale_z2d::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_scale_z2d(*args[0]);
    }

    scale_z2logOR::scale_z2logOR() :ScalarFunction("scale_z2logOR", 1)
    {}
    bool scale_z2logOR::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double scale_z2logOR::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_scale_z2logOR(*args[0]);
    }

    scale_z2r::scale_z2r() :ScalarFunction("scale_z2r", 1)
    {}
    bool scale_z2r::checkParameterValue(std::vector<double const *> const &args) const
    {
      return(true);
    }
    double scale_z2r::evaluate(std::vector<double const *> const &args) const
    {
      return cpp_scale_z2r(*args[0]);
    }

  }
}
