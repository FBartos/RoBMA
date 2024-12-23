#include "DWEIGHTSMIX.h"

#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <JRmath.h> // For rgamma, etc.
#include <numeric> // For std::accumulate

namespace jags {
namespace RoBMA {

DWEIGHTSMIX::DWEIGHTSMIX() : ArrayDist("weights_mix", 4) {}

// Dimension of the output
std::vector<unsigned int> DWEIGHTSMIX::dim(std::vector<std::vector<unsigned int>> const &dims) const {
    return std::vector<unsigned int>(1,dims[1][0]);
}

// Checking parameter dimensions
bool DWEIGHTSMIX::checkParameterDim(std::vector<std::vector<unsigned int>> const &dims) const {

    bool alpha_and_index_OK = dims[1][1] == dims[0][1];
    bool row_indicator_OK   = dims[3][0] == 1;

    return alpha_and_index_OK && row_indicator_OK;
}

// Checking parameter values
bool DWEIGHTSMIX::checkParameterValue(std::vector<double const *> const &par,
                                      std::vector<std::vector<unsigned int>> const &dims) const {
    return true;
}

// Log Density
double DWEIGHTSMIX::logDensity(double const *x, unsigned int length, PDFType type,
                               std::vector<double const *> const &par,
                               std::vector<std::vector<unsigned int>> const &dims,
                               double const *lower, double const *upper) const {
    return 0.0; // Placeholder: Not implemented
}

// Random Sample
void DWEIGHTSMIX::randomSample(double *x, unsigned int length,
                               std::vector<double const *> const &par,
                               std::vector<std::vector<unsigned int>> const &dims,
                               double const *lower, double const *upper,
                               RNG *rng) const {
    // extract parameters
    const double *alphaMat = par[0];
    const double *indexMat = par[1];
    const double *indexMax = par[2];
    double rowIndicator = par[3][0];

    // extract dimensions
    int ncol = dims[1][0];
    //unsigned int nrow = dims[1][1];

    // extract the selected row and the maximum index
    int selectedRow      = static_cast<int>(rowIndicator - 1);
    int selectedIndexMax = static_cast<int>(indexMax[selectedRow]);

    // ---- deal with non weightfunction cases ---- //
    if (selectedIndexMax == 0) {
        for (int i = 0; i < ncol; ++i) {
            x[i] = 1.0;
        }
        return;
    }

    // --- deal with fixed weightfunctions --- //
    if (selectedIndexMax == -1) {
        for (int i = 0; i < ncol; ++i) {
            x[i] = alphaMat[static_cast<int>(selectedRow * ncol + indexMat[selectedRow * ncol + i] - 1)];
        }
        return;
    }

    // --- sample etas from the gamma distribution --- //
    std::vector<double> eta(selectedIndexMax);
    for (int i = 0; i < selectedIndexMax; ++i) {
        eta[i] = rgamma(alphaMat[selectedRow * ncol + i], 1.0, rng);
    }

    // --- normalized etas --- //
    std::vector<double> std_eta(selectedIndexMax);
    double eta_sum = std::accumulate(eta.begin(), eta.end(), 0.0); // Sum of eta
    for (int i = 0; i < selectedIndexMax; ++i) {
        std_eta[i] = eta[i] / eta_sum;
    }

    // --- transform to cummulative dirichlet distribution --- //
    std::vector<double> omega(selectedIndexMax);
    omega[0] = std_eta[0];
    for (int i = 1; i < selectedIndexMax; ++i) {
        omega[i] = omega[i - 1] + std_eta[i];
    }

    // --- map using the indexMat to the correct index --- //
    for (int i = 0; i < ncol; ++i) {
        x[i] = omega[static_cast<int>(indexMat[selectedRow * ncol + i] - 1)];
    }
}

// Support
void DWEIGHTSMIX::support(double *lower, double *upper, unsigned int length,
                          std::vector<double const *> const &par,
                          std::vector<std::vector<unsigned int>> const &dims) const {
    for (unsigned int i = 0; i < length; ++i) {
        lower[i] = 0.0;
        upper[i] = JAGS_POSINF;
    }
}

// Typical Value
void DWEIGHTSMIX::typicalValue(double *x, unsigned int length,
                               std::vector<double const *> const &par,
                               std::vector<std::vector<unsigned int>> const &dims,
                               double const *lower, double const *upper) const {
    // Not implemented
}

// Support Fixed
bool DWEIGHTSMIX::isSupportFixed(std::vector<bool> const &fixmask) const {
    return true;
}

}} // namespace jags::RoBMA
