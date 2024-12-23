#ifndef DWEIGHTSMIX_H_
#define DWEIGHTSMIX_H_
#include <distribution/ArrayDist.h>

namespace jags {
namespace RoBMA {

class DWEIGHTSMIX : public ArrayDist {
public:
    DWEIGHTSMIX();

    double logDensity(double const *x, unsigned int length, PDFType type,
                      std::vector<double const *> const &par,
                      std::vector<std::vector<unsigned int>> const &dims,
                      double const *lower, double const *upper) const;

    void randomSample(double *x, unsigned int length,
                      std::vector<double const *> const &par,
                      std::vector<std::vector<unsigned int>> const &dims,
                      double const *lower, double const *upper,
                      RNG *rng) const;

    std::vector<unsigned int> dim(std::vector<std::vector<unsigned int>> const &dims) const;

    bool checkParameterDim(std::vector<std::vector<unsigned int>> const &dims) const;

    bool checkParameterValue(std::vector<double const *> const &par,
                             std::vector<std::vector<unsigned int>> const &dims) const;

    void support(double *lower, double *upper, unsigned int length,
                 std::vector<double const *> const &par,
                 std::vector<std::vector<unsigned int>> const &dims) const;

    void typicalValue(double *x, unsigned int length,
                      std::vector<double const *> const &par,
                      std::vector<std::vector<unsigned int>> const &dims,
                      double const *lower, double const *upper) const;

    bool isSupportFixed(std::vector<bool> const &fixmask) const;
};

}} // namespace jags::RoBMA

#endif /* DWEIGHTSMIX_H_ */
