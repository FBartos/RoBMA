#include "DWEIGHTSMIX.h"

#include <rng/RNG.h>
#include <util/dim.h>
#include <util/nainf.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <JRmath.h>
#include <numeric>

#include "../source/tools.h"

namespace jags {
namespace RoBMA {

DWEIGHTSMIX::DWEIGHTSMIX() : ArrayDist("weights_mix", 3) {}

// Dimension of the output
std::vector<unsigned int> DWEIGHTSMIX::dim(std::vector<std::vector<unsigned int>> const &dims) const {
    return std::vector<unsigned int>(1,dims[1][0]);
}

// Checking parameter dimensions
bool DWEIGHTSMIX::checkParameterDim(std::vector<std::vector<unsigned int>> const &dims) const {
    return true;
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


    // extract parameters
    const double *alphaMat = par[0];
    const double *indexMat = par[1];
    const int indexMax  = static_cast<int>(*par[2]);

    // extract dimensions
    int ncol = dims[1][0];

    // ---- deal with non weightfunction cases ---- //
    if (indexMax == 0) {
        return 0.0;
    }

    // --- deal with fixed weightfunctions --- //
    if (indexMax == -1) {
        return 0.0;
    }

    // --- extract omegas from x using the indexMat --- //
    std::vector<double> omega(indexMax);
    for (int i = 0; i < ncol; ++i) {
        omega[static_cast<int>(indexMat[i]) - 1] = x[i];
    }

    // --- extract normalized etas from the omegas --- //
    std::vector<double> std_eta(indexMax);
    std_eta[0] = omega[0];
    for (int i = 1; i < indexMax; ++i) {
        std_eta[i] = omega[i] - omega[i - 1];
    }

    // --- extract the corresponding alpha parameters --- //
    std::vector<double> alpha(indexMax);
    for (int i = 0; i < indexMax; ++i) {
        alpha[i] = alphaMat[i];
    }

    // --- use density of dirichlet distribution to compute the log density --- //
    double log_lik = ddirichlet(std_eta, alpha);

    return log_lik;
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
    const int indexMax  = static_cast<int>(*par[2]);

    // extract dimensions
    int ncol = dims[1][0];
    //unsigned int nrow = dims[1][1];

    // ---- deal with non weightfunction cases ---- //
    if (indexMax == 0) {
        for (int i = 0; i < ncol; ++i) {
            x[i] = 1.0;
        }
        return;
    }

    // --- deal with fixed weightfunctions --- //
    if (indexMax == -1) {
        for (int i = 0; i < ncol; ++i) {
            x[i] = alphaMat[ static_cast<int>(indexMat[i]) - 1 ];
        }
        return;
    }

    // --- sample etas from the gamma distribution --- //
    std::vector<double> eta(indexMax);
    for (int i = 0; i < indexMax; ++i) {
        eta[i] = rgamma(alphaMat[i], 1.0, rng);
    }

    // --- normalized etas --- //
    std::vector<double> std_eta(indexMax);
    double eta_sum = std::accumulate(eta.begin(), eta.end(), 0.0); // Sum of eta
    for (int i = 0; i < indexMax; ++i) {
        std_eta[i] = eta[i] / eta_sum;
    }

    // --- transform to cummulative dirichlet distribution --- //
    std::vector<double> omega(indexMax);
    omega[0] = std_eta[0];
    for (int i = 1; i < indexMax; ++i) {
        omega[i] = omega[i - 1] + std_eta[i];
    }

    // --- map using the indexMat to the correct index --- //
    for (int i = 0; i < ncol; ++i) {
        x[i] = omega[ static_cast<int>(indexMat[i]) - 1 ];
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
    // extract parameters
    const double *alphaMat = par[0];
    const double *indexMat = par[1];
    const int indexMax  = static_cast<int>(*par[2]);

    // extract dimensions
    int ncol = dims[1][0];

    // ---- deal with non weightfunction cases ---- //
    if (indexMax == 0) {
        for (int i = 0; i < ncol; ++i) {
            x[i] = 1.0;
        }
        return;
    }

    // --- deal with fixed weightfunctions --- //
    if (indexMax == -1) {
        for (int i = 0; i < ncol; ++i) {
            x[i] = alphaMat[ static_cast<int>(indexMat[i]) - 1 ];
        }
        return;
    }

    // --- sample etas from the gamma distribution --- //
    std::vector<double> eta(indexMax);
    for (int i = 0; i < indexMax; ++i) {
        eta[i] = alphaMat[i];
    }

    // --- normalized etas --- //
    std::vector<double> std_eta(indexMax);
    double eta_sum = std::accumulate(eta.begin(), eta.end(), 0.0); // Sum of eta
    for (int i = 0; i < indexMax; ++i) {
        std_eta[i] = eta[i] / eta_sum;
    }

    // --- transform to cummulative dirichlet distribution --- //
    std::vector<double> omega(indexMax);
    omega[0] = std_eta[0];
    for (int i = 1; i < indexMax; ++i) {
        omega[i] = omega[i - 1] + std_eta[i];
    }

    // --- map using the indexMat to the correct index --- //
    for (int i = 0; i < ncol; ++i) {
        x[i] = omega[ static_cast<int>(indexMat[i]) - 1 ];
    }
}

// Support Fixed
bool DWEIGHTSMIX::isSupportFixed(std::vector<bool> const &fixmask) const {
    return false;
}

}} // namespace jags::RoBMA
