#include "LogWeightedtFun.h"
#include <util/nainf.h>
#include <JRmath.h>
#include <cmath>

using std::vector;
using std::string;

using std::vector;
using std::log;
using std::exp;


namespace jags {
namespace weightedt {

LogWeightedtFun::LogWeightedtFun() :ScalarFunction("logweightedt2", 5)
{}

bool LogWeightedtFun::checkParameterValue(vector<double const *> const &args) const
{
  return(*args[1] > 2 && *args[4] >= 0.0 && *args[4] <= 1.0);
}

double LogWeightedtFun::evaluate(vector<double const *> const &args) const
{
  double x      = *args[0];
  double df     = *args[1];
  double ncp    = *args[2];
  double crit_t = *args[3];
  double beta   = *args[4];

  double nom;
  double denom;
  double denom1;
  double denom2;
  double denom3;
  double log_lik;
  double w;
  
  if(x < crit_t){
    w = log(beta);
  }else{
    w = log(1);
  }
  
  if(ncp > 0){
  
    nom    = dnt(-x, df, -ncp, true)+w;
  
    denom1 = pnt(-crit_t, df, -ncp, true, false);
    denom2 = log(
      pnt( crit_t, df, -ncp, true, false) - 
      pnt(-crit_t, df, -ncp, true, false)
      ) + log(beta);
    denom3 = pnt( crit_t, df, -ncp, false, false);

  }else{

    nom    = dnt(x, df, ncp, true)+w;
  
    denom1 = pnt(-crit_t, df, ncp, true, false);
    denom2 = log(
      pnt( crit_t, df, ncp, true, false) - 
      pnt(-crit_t, df, ncp, true, false)
      ) + log(beta);
    denom3 = pnt( crit_t, df, ncp, false, false);
  }

  denom   = log(denom1 + exp(denom2) + denom3);
  log_lik = nom-denom;

  return log_lik;
}

}
} // namespace jags
