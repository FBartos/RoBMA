// Standalone numerical check for the analytic Gaussian block.
// Compile as C++17 or later and run the resulting executable.
#include "../../../src/samplers/GaussianJointKernel.h"

#include <random>
#include <stdexcept>

namespace {

namespace kernel = robma_gaussian_joint;
using Matrix = std::vector<std::vector<long double>>;
using Vector = std::vector<long double>;

constexpr long double log_two_pi = 1.837877066409345483560659472811235L;
using Coefficients = std::array<long double, 2>;

void check(long double actual, long double expected)
{
  if (!std::isfinite(actual) ||
      std::fabs(actual - expected) > 3e-11L * std::max(1.0L, std::fabs(expected))) {
    throw std::runtime_error("Gaussian joint kernel disagrees with the observation-space reference.");
  }
}

Vector forward_solve(Matrix const &lower, Vector rhs)
{
  for (std::size_t i = 0; i < rhs.size(); ++i) {
    for (std::size_t j = 0; j < i; ++j) rhs[i] -= lower[i][j] * rhs[j];
    rhs[i] /= lower[i][i];
  }
  return rhs;
}

long double dot(Vector const &a, Vector const &b)
{
  long double sum = 0;
  for (std::size_t i = 0; i < a.size(); ++i) sum += a[i] * b[i];
  return sum;
}

struct Reference {
  long double log_likelihood;
  Coefficients mean;
  std::array<Coefficients, 2> covariance;
};

Reference observation_space(kernel::Model const &model, unsigned int cell, double tau)
{
  const std::size_t n = model.rows.size();
  Matrix covariance(n, Vector(n));
  std::array<Vector, 2> cross_covariance{Vector(n), Vector(n)};
  Vector residual(n);
  for (std::size_t i = 0; i < n; ++i) {
    const auto &row = model.rows[i];
    const Coefficients design {{
      kernel::active(cell, 0) ? model.effect_sign : 0,
      kernel::active(cell, 1) ? model.effect_sign * row.moderator : 0
    }};
    residual[i] = row.y - design[0] * model.mean[0] - design[1] * model.mean[1];
    for (unsigned int k = 0; k < 2; ++k) {
      cross_covariance[k][i] = design[k] * model.sd[k] * model.sd[k];
    }
    covariance[i][i] = row.se * row.se + tau * tau;
    for (std::size_t j = 0; j <= i; ++j) {
      covariance[i][j] += cross_covariance[0][i] * design[0] +
        cross_covariance[1][i] * model.effect_sign *
          (kernel::active(cell, 1) ? model.rows[j].moderator : 0);
    }
  }
  Reference out;
  out.log_likelihood = -0.5L * n * log_two_pi;
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      for (std::size_t k = 0; k < j; ++k) covariance[i][j] -= covariance[i][k] * covariance[j][k];
      if (i == j) {
        covariance[i][j] = std::sqrt(covariance[i][j]);
        out.log_likelihood -= std::log(covariance[i][j]);
      } else {
        covariance[i][j] /= covariance[j][j];
      }
    }
  }
  const auto solved_residual = forward_solve(covariance, residual);
  out.log_likelihood -= 0.5L * dot(solved_residual, solved_residual);
  for (unsigned int k = 0; k < 2; ++k) {
    cross_covariance[k] = forward_solve(covariance, cross_covariance[k]);
    out.mean[k] = model.mean[k] + dot(cross_covariance[k], solved_residual);
  }
  for (unsigned int i = 0; i < 2; ++i) {
    for (unsigned int j = 0; j < 2; ++j) {
      out.covariance[i][j] = (i == j ? model.sd[i] * model.sd[i] : 0)
        - dot(cross_covariance[i], cross_covariance[j]);
    }
  }
  return out;
}

void check_conditionals(kernel::Model const &model)
{
  for (unsigned int cell = 0; cell < 4; ++cell) {
    for (double tau : {0., 1e-12, .2, 1.8, 9.}) {
      const auto exact = observation_space(model, cell, tau);
      const auto at = kernel::evaluate(model, cell, tau);
      check(at.log_likelihood, exact.log_likelihood);
      const auto center = kernel::restore_coefficients(model, at, 0, 0);
      const auto one = kernel::restore_coefficients(model, at, 1, 0);
      const auto two = kernel::restore_coefficients(model, at, 0, 1);
      for (unsigned int i = 0; i < 2; ++i) {
        check(center[i], exact.mean[i]);
        for (unsigned int j = 0; j < 2; ++j) {
          const auto value = (one[i] - center[i]) * (one[j] - center[j])
            + (two[i] - center[i]) * (two[j] - center[j]);
          check(value, exact.covariance[i][j]);
        }
      }
    }
  }
}

struct SeededUniform {
  std::mt19937 generator{7186};
  std::uniform_real_distribution<double> draws{0, 1};
  double uniform() { return draws(generator); }
};

void check_proposal(kernel::Model const &model)
{
  const auto quantile = [&](double probability) {
    double lower = 0;
    double upper = 20;
    for (unsigned int i = 0; i < 70; ++i) {
      const double mid = (lower + upper) / 2;
      if (std::erfc(mid / std::sqrt(2.0)) > probability) lower = mid;
      else upper = mid;
    }
    return (lower + upper) / 2 * model.sd[2];
  };
  kernel::Proposal proposal;
  if (!proposal.prepare(model, quantile)) throw std::runtime_error("Proposal preparation failed.");
  long double total = 0;
  for (unsigned int cell = 0; cell < kernel::Proposal::coefficient_cells; ++cell) {
    for (unsigned int bin = 0; bin < kernel::Proposal::scale_bins; ++bin) {
      total += std::exp(proposal.log_density(cell, bin)) / (bin == 0 ? 1 : kernel::Proposal::slab_bins);
    }
  }
  check(total, 1);

  SeededUniform rng;
  for (unsigned int rep = 0; rep < 400; ++rep) {
    const auto candidate = proposal.draw(model, quantile, rng);
    const auto current = proposal.state(model, rep % 4, 0, 0);
    const auto physical = [&](const kernel::Proposal::State &state) {
      const bool scale_on = state.bin > 0;
      if (scale_on && state.bin != proposal.bucket(state.tau)) {
        throw std::runtime_error("Scale draw left its proposal bin.");
      }
      const auto exact = observation_space(model, state.coefficients, state.tau);
      const long double z = state.tau / model.sd[2];
      const long double log_half_normal = scale_on ?
        std::log(2.0L) - std::log(model.sd[2]) - 0.5L * (log_two_pi + z * z) : 0;
      const long double target = exact.log_likelihood + log_half_normal +
        model.log_prior[0][(state.coefficients % 2) == 1] +
        model.log_prior[1][state.coefficients > 1] + model.log_prior[2][scale_on];
      const long double proposal_density =
        proposal.log_density(state.coefficients, state.bin) + log_half_normal;
      check(state.log_ratio, target - proposal_density);
      return std::array<long double, 2>{target, proposal_density};
    };
    const auto previous = physical(current);
    const auto proposed = physical(candidate);
    check(previous[0] + proposed[1] + std::min(0.0L, candidate.log_ratio - current.log_ratio),
          proposed[0] + previous[1] + std::min(0.0L, current.log_ratio - candidate.log_ratio));
  }
}

void check_location_cancellation()
{
  kernel::Model model;
  model.rows = {{0.7, 0.2, 0}};
  model.effect_sign = 1;
  model.mean = {{1e50, 0}};
  model.sd = {{1e50, 1, 1}};
  const auto at = kernel::evaluate(model, 1, 0);
  // At this scale the prior's shift of the conditional mean is 4e-52 and
  // its effect on the conditional SD is below the precision of the check.
  check(kernel::restore_coefficients(model, at, 0, 0)[0], 0.7L);
  check(kernel::restore_coefficients(model, at, 1, 0)[0], 0.9L);
  check(at.log_likelihood, -0.5L * (log_two_pi + 1) - std::log(1e50L));
}

} // namespace

int main()
{
  kernel::Model model;
  model.rows = {{-.87, .09, -1.3}, {-.5, .4, -.77}, {1.02, .2, .42}, {.21, .51, .74}};
  model.mean = {{.31, -.19}};
  model.sd = {{.8, 1.4, .61}};
  model.log_prior = {{
    {{std::log(.23L), std::log(.77L)}},
    {{std::log(.62L), std::log(.38L)}},
    {{std::log(.88L), std::log(.12L)}}
  }};
  for (double sign : {1.0, -1.0}) {
    model.effect_sign = sign;
    check_conditionals(model);
    check_proposal(model);
  }
  check_location_cancellation();
}
