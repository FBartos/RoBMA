#include "r.h"
#include "../source/transformations.h"
#include <util/nainf.h>
#include <JRmath.h>
#include <cmath>

using std::vector;
using std::fabs;

namespace jags {
  namespace RoBMA {

    // effect sizes transformations
    r2d::r2d() :ScalarFunction("r2d", 1)
    {}
    bool r2d::checkParameterValue(vector<double const *> const &args) const
    {
      return(fabs(*args[0]) < 1);
    }
    double r2d::evaluate(vector<double const *> const &args) const
    {
      return cpp_r2d(*args[0]);
    }

    r2z::r2z() :ScalarFunction("r2z", 1)
    {}
    bool r2z::checkParameterValue(vector<double const *> const &args) const
    {
      return(fabs(*args[0]) < 1);
    }
    double r2z::evaluate(vector<double const *> const &args) const
    {
      return cpp_r2z(*args[0]);
    }

    r2logOR::r2logOR() :ScalarFunction("r2logOR", 1)
    {}
    bool r2logOR::checkParameterValue(vector<double const *> const &args) const
    {
      return(fabs(*args[0]) < 1);
    }
    double r2logOR::evaluate(vector<double const *> const &args) const
    {
      return cpp_r2logOR(*args[0]);
    }

    // standard errors transformations
    se_r2se_d::se_r2se_d() :ScalarFunction("se_r2se_d", 2)
    {}
    bool se_r2se_d::checkParameterValue(vector<double const *> const &args) const
    {
      return(*args[0] >= 0 && fabs(*args[1]) < 1);
    }
    double se_r2se_d::evaluate(vector<double const *> const &args) const
    {
      return cpp_se_r2se_d(*args[0], *args[1]);
    }

    se_r2se_z::se_r2se_z() :ScalarFunction("se_r2se_z", 2)
    {}
    bool se_r2se_z::checkParameterValue(vector<double const *> const &args) const
    {
      return(*args[0] >= 0 && fabs(*args[1]) < 1 &&  cpp_n_r(*args[1], *args[0]) > 3);
    }
    double se_r2se_z::evaluate(vector<double const *> const &args) const
    {
      return cpp_se_r2se_z(*args[0], *args[1]);
    }

    se_r2se_logOR::se_r2se_logOR() :ScalarFunction("se_r2se_logOR", 2)
    {}
    bool se_r2se_logOR::checkParameterValue(vector<double const *> const &args) const
    {
      return(*args[0] >= 0 && fabs(*args[1]) < 1);
    }
    double se_r2se_logOR::evaluate(vector<double const *> const &args) const
    {
      return cpp_se_r2se_logOR(*args[0], *args[1]);
    }

    // linear scaling function (not used, for completeness)
    scale_r2d::scale_r2d() :ScalarFunction("scale_r2d", 1)
    {}
    bool scale_r2d::checkParameterValue(vector<double const *> const &args) const
    {
      return(true);
    }
    double scale_r2d::evaluate(vector<double const *> const &args) const
    {
      return cpp_scale_r2d(*args[0]);
    }

    scale_r2z::scale_r2z() :ScalarFunction("scale_r2z", 1)
    {}
    bool scale_r2z::checkParameterValue(vector<double const *> const &args) const
    {
      return(true);
    }
    double scale_r2z::evaluate(vector<double const *> const &args) const
    {
      return cpp_scale_r2z(*args[0]);
    }

    scale_r2logOR::scale_r2logOR() :ScalarFunction("scale_r2logOR", 1)
    {}
    bool scale_r2logOR::checkParameterValue(vector<double const *> const &args) const
    {
      return(true);
    }
    double scale_r2logOR::evaluate(vector<double const *> const &args) const
    {
      return cpp_scale_r2logOR(*args[0]);
    }
  }
} 
