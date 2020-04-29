#ifndef WEIGHTEDT_FUNC_H_
#define WEIGHTEDT_FUNC_H_

#include <function/ScalarFunction.h>

namespace jags {
namespace weightedt {

class LogWeightedtFun : public ScalarFunction 
{
  public:
    LogWeightedtFun();

    bool checkParameterValue(std::vector<double const *> const &args) const;
    double evaluate(std::vector<double const *> const &args) const;
};

}
}

#endif /* WEIGHTEDT_FUNC_H_ */
