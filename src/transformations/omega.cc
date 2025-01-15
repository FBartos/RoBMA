#include "omega.h"

#include <cmath>
#include <algorithm>
#include <numeric>


namespace jags {
  namespace RoBMA {

    eta2omega::eta2omega() : VectorFunction("eta2omega", 4) // Function name and number of arguments
    {
    }

    void eta2omega::evaluate(double *value,
                         std::vector<double const *> const &args,
                         std::vector<unsigned int> const &lengths) const
    {
        // extract parameters
        const double *eta_vector  = args[0];
        const double *omega_index = args[1];
        const double *eta_index   = args[2];
        const int eta_index_max   = static_cast<int>(*args[3]);

        // extract dimensions
        int K = lengths[1];

        // ---- deal with non weightfunction cases ---- //
        if (eta_index_max == 0) {
            for (int i = 0; i < K; ++i) {
                value[i] = 1.0;
            }
            return;
        }

        // --- deal with fixed weightfunctions --- //
        if (eta_index_max == -1) {
            for (int i = 0; i < K; ++i) {
              value[i] = omega_index[i];
            }
            return;
        }

        // ---- select the correct etas ---- //
        std::vector<double> eta(eta_index_max);
        for (int i = 0; i < eta_index_max; ++i) {
            eta[i] = eta_vector[ static_cast<int>(eta_index[i]) - 1 ];
        }

        // --- normalized etas --- //
        std::vector<double> std_eta(eta_index_max);
        double eta_sum = std::accumulate(eta.begin(), eta.end(), 0.0); // Sum of eta
        for (int i = 0; i < eta_index_max; ++i) {
            std_eta[i] = eta[i] / eta_sum;
        }

        // --- transform to cummulative dirichlet distribution --- //
        std::vector<double> omega(eta_index_max);
        omega[0] = std_eta[0];
        for (int i = 1; i < eta_index_max; ++i) {
            omega[i] = omega[i - 1] + std_eta[i];
        }

        // --- map using the indexMat to the correct index --- //
        for (int i = 0; i < K; ++i) {
            value[i] = omega[ static_cast<int>(omega_index[i]) - 1 ];
        }
    }

    unsigned int eta2omega::length(std::vector<unsigned int> const &lengths,
                                  std::vector<double const *> const &args) const
    {
        // Output length matches the length of the second std::vector
        return lengths[1];
    }

    bool eta2omega::checkParameterLength(std::vector<unsigned int> const &len) const
    {
        return true;
    }

    bool eta2omega::isDiscreteValued(std::vector<bool> const &mask) const
    {
        return false;
    }

    bool eta2omega::checkParameterDiscrete(std::vector<bool> const &mask) const
    {
        return true;
    }

    bool eta2omega::checkParameterFixed(std::vector<bool> const &mask) const
    {
        return true;
    }

    bool eta2omega::checkParameterValue(std::vector<double const *> const &args,
                                       std::vector<unsigned int> const &lens) const
    {
        return true;
    }
}}
