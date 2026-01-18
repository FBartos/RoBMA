#ifndef OMEGA_H_
#define OMEGA_H_

#include <function/VectorFunction.h>

namespace jags {
    namespace RoBMA {

    class eta2omega : public VectorFunction
    {
    public:
        eta2omega();
        void evaluate(double *value, 
                      std::vector<double const *> const &values,
                      std::vector<unsigned long> const &lengths) const;
        unsigned long length(std::vector<unsigned long> const &lengths,
                            std::vector<double const *> const &values) const;
        bool checkParameterLength(std::vector<unsigned long> const &args) const;
        bool isDiscreteValued(std::vector<bool> const &mask) const;
        bool checkParameterDiscrete(std::vector<bool> const &mask) const;
        bool checkParameterFixed(std::vector<bool> const &mask) const;
        bool checkParameterValue(std::vector<double const *> const &args,
                                 std::vector<unsigned long> const &lens) const;
    };

    }
}

#endif /* OMEGA_H_ */
