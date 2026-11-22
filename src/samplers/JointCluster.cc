#include "JointCluster.h"

#include <BayesTools/sampler_declarations.h>

#include <JRmath.h>
#include <distribution/Distribution.h>
#include <graph/StochasticNode.h>
#include <module/ModuleError.h>
#include <rng/RNG.h>
#include <sampler/GraphView.h>
#include <sampler/MutableSampler.h>
#include <sampler/SamplerFactory.h>
#include <sampler/SingletonGraphView.h>
#include <sampler/Slicer.h>
#include <util/nainf.h>

#include "../distributions/selnorm-jags-bounds.h"
#include "../selnorm/selnorm.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <vector>

namespace jags {
namespace RoBMA {
namespace {

constexpr char member_name[] = "robma_cluster_member";
constexpr char setup_name[] = "robma_cluster_setup";
constexpr char sampler_name[] = "RoBMA::JointCluster";
constexpr double log_two_pi = 1.8378770664093455;

enum Role { gate = 1, intercept, moderator, study, total_sd, allocation };
enum Setup { group_gate, design, groups, response, sampling_se, sampling_variance,
             observed_bin, lower, upper, omega, mode, telescope, intercept_gate,
             intercept_label, moderator_label, physical_tau, physical_rho,
             heterogeneity_gate, heterogeneity_label, pet, peese, study_count,
             cadence, scale_enabled, setup_count };

using Params = std::vector<Node const *>;
using Vec = std::vector<double>;
using BayesTools::sampler_declarations::depends_on;

bool integer_between(double x, unsigned long first, unsigned long last)
{
  return x >= first && x <= last && x == std::floor(x);
}

bool valid_lengths(std::vector<unsigned long> const &size)
{
  if (size.size() != setup_count || size[groups] < 2 ||
      size[design] == 0 || size[design] % size[groups] != 0 ||
      size[omega] == 0) return false;
  for (unsigned int i : {response, sampling_se, sampling_variance, observed_bin}) {
    if (size[i] != size[groups]) return false;
  }
  if (size[lower] != size[omega] || size[upper] != size[omega]) return false;
  for (unsigned int i : {group_gate, mode, telescope, intercept_gate, intercept_label,
                         moderator_label, physical_tau, physical_rho, heterogeneity_gate,
                         heterogeneity_label, pet, peese,
                         study_count, cadence, scale_enabled}) {
    if (size[i] != 1) return false;
  }
  return true;
}

bool normal_parameters(StochasticNode const *node, Vec &mean, Vec &precision,
                       bool half_normal = false)
{
  unsigned long width = node->length();
  auto const &parents = node->parents();
  bool scalar = node->distribution()->name() == "dnorm" && width == 1;
  bool vector = node->distribution()->name() == "dmnorm" && !half_normal;
  if ((!scalar && !vector) || node->upperBound() != nullptr ||
      (half_normal ? !fixed_scalar(node->lowerBound()) || node->lowerBound()->value(0)[0] != 0
                   : node->lowerBound() != nullptr) ||
      parents.size() != (half_normal ? 3 : 2) || !parents[0]->isFixed() ||
      !parents[1]->isFixed() || parents[0]->length() != width ||
      parents[1]->length() != width * width) return false;

  for (unsigned long i = 0; i < width; ++i) {
    double location = parents[0]->value(0)[i];
    double diagonal = parents[1]->value(0)[i * (width + 1)];
    if (!std::isfinite(location) || (half_normal && location != 0) ||
        !(diagonal > 0 && std::isfinite(diagonal) &&
          std::isfinite(1 / std::sqrt(diagonal)))) return false;
    for (unsigned long j = 0; j < width; ++j) {
      if (i != j && parents[1]->value(0)[i + j * width] != 0) return false;
    }
    mean.push_back(location);
    precision.push_back(diagonal);
  }
  return true;
}

bool beta_parameters(StochasticNode const *node)
{
  auto const &p = node->parents();
  return node->length() == 1 && node->distribution()->name() == "dbeta" &&
    !node->lowerBound() && !node->upperBound() && p.size() == 2 &&
    fixed_scalar(p[0]) && fixed_scalar(p[1]) && std::isfinite(p[0]->value(0)[0]) &&
    std::isfinite(p[1]->value(0)[0]) && p[0]->value(0)[0] > 0 && p[1]->value(0)[0] > 0;
}

double scalar(Node const *node, unsigned int chain)
{
  return node->value(chain)[0];
}

// Lower triangular Cholesky and triangular solves for the small fixed-effect
// block. Study coordinates are eliminated separately as a diagonal block.
bool cholesky(Vec &matrix, std::size_t width)
{
  for (std::size_t i = 0; i < width; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      long double value = matrix[i * width + j];
      for (std::size_t k = 0; k < j; ++k) {
        value -= static_cast<long double>(matrix[i * width + k]) * matrix[j * width + k];
      }
      if (i == j) {
        if (!(value > 0 && std::isfinite(value))) return false;
        matrix[i * width + j] = std::sqrt(value);
      } else {
        matrix[i * width + j] = value / matrix[j * width + j];
      }
    }
  }
  return true;
}

void solve_upper(Vec const &lower_chol, Vec &rhs)
{
  std::size_t width = rhs.size();
  for (std::size_t i = width; i-- > 0;) {
    for (std::size_t j = i + 1; j < width; ++j) rhs[i] -= lower_chol[j * width + i] * rhs[j];
    rhs[i] /= lower_chol[i * width + i];
  }
}

struct FitData {
  Params setup;
  std::vector<unsigned long> studies;
  Vec y, se, variance, x, prior_mean, precision;
  std::vector<int> bins;
  Vec z_lower, z_upper;
  std::array<double, 2> log_gate{};
  unsigned long n_studies = 0;
  unsigned long width = 0;
  unsigned long spacing = 0;
  bool scale = false;
};

struct Context {
  FitData const &data;
  const double *weights;
  SelNormKernelData kernel{};
  double between;
  Vec sd, offset;
  bool selected;
  bool intercept_on;

  Context(FitData const &data_, unsigned int chain)
    : data(data_), weights(data.setup[omega]->value(chain)),
      sd(data.y.size()), offset(data.y.size())
  {
    double tau = scalar(data.setup[physical_tau], chain);
    double rho = scalar(data.setup[physical_rho], chain);
    double within = tau * std::sqrt(1 - rho);
    between = tau * std::sqrt(rho);
    kernel.n_bins = static_cast<int>(data.z_lower.size());
    kernel.z_lower = data.z_lower.data();
    kernel.z_upper = data.z_upper.data();
    kernel.effect_sign = 1;
    kernel.trusted_step_partition = true;
    kernel.telescope_probabilities = scalar(data.setup[telescope], chain) == 1;
    selected = scalar(data.setup[mode], chain) == SELKERNEL_STEP;
    intercept_on = scalar(data.setup[intercept_gate], chain) ==
      scalar(data.setup[intercept_label], chain);
    double pet_value = scalar(data.setup[pet], chain);
    double peese_value = scalar(data.setup[peese], chain);
    for (std::size_t i = 0; i < data.y.size(); ++i) {
      sd[i] = std::hypot(std::sqrt(data.variance[i]), within);
      offset[i] = data.se[i] * (pet_value + data.se[i] * peese_value);
    }
  }

  double coefficient(std::size_t row, std::size_t column) const
  {
    return column == 0 ? 1 : data.x[(column - 1) * data.y.size() + row];
  }

  double row_score(std::size_t i, double mean, bool use_selection,
                   double *gradient, double *weight) const
  {
    double row_sd = sd[i];
    double y = data.y[i];
    double residual = (y - mean) / row_sd;
    if (!use_selection) {
      if (gradient) {
        *gradient = residual / row_sd;
        *weight = (1 / row_sd) / row_sd;
      }
      return -0.5 * (log_two_pi + residual * residual) - std::log(row_sd);
    }
    if (gradient) {
      double m = 0, second = 0, var = 0;
      cpp_selnorm_kernel_moments(mean, row_sd, data.se[i], weights, 0, 0,
                                SELKERNEL_STEP, kernel, &m, &second, 1, false, &var);
      *gradient = ((y - m) / row_sd) / row_sd;
      *weight = ((var / row_sd) / row_sd) / row_sd / row_sd;
    }
    return cpp_selnorm_kernel_lpdf(y, mean, row_sd, mean, row_sd, data.se[i], 1,
                                 weights, data.bins[i], 0, 0, SELKERNEL_STEP,
                                 kernel, 1, false);
  }
};

struct Reference {
  Vec beta, gamma, chol, study_precision, study_solve;
  std::vector<unsigned long> active;
  double log_mass = 0, log_det_sqrt = 0;

  explicit Reference(Context const &ctx, bool included)
    : gamma(ctx.data.n_studies, 0)
  {
    if (ctx.intercept_on) active.push_back(0);
    if (included) for (unsigned long i = 1; i <= ctx.data.width; ++i) active.push_back(i);
    for (unsigned long i : active) beta.push_back(ctx.data.prior_mean[i]);
  }

  double evaluate(Context const &ctx, Vec const &fixed, Vec const &latent,
                  bool use_selection, Vec *fixed_step = nullptr, Vec *study_step = nullptr)
  {
    const std::size_t p = active.size();
    const std::size_t n = gamma.size();
    Vec fixed_score, study_score;
    if (fixed_step) {
      chol.assign(p * p, 0);
      study_precision.assign(n, 1);
      study_solve.assign(n * p, 0);
      fixed_score.resize(p);
      study_score.resize(n);
    }

    double target = 0;
    for (std::size_t a = 0; a < p; ++a) {
      unsigned long column = active[a];
      double precision = ctx.data.precision[column];
      double delta = fixed[a] - ctx.data.prior_mean[column];
      target += 0.5 * (std::log(precision) - log_two_pi - precision * delta * delta);
      if (fixed_step) {
        fixed_score[a] = -precision * delta;
        chol[a * p + a] = precision;
      }
    }
    for (std::size_t j = 0; j < n; ++j) {
      target -= 0.5 * (log_two_pi + latent[j] * latent[j]);
      if (fixed_step) study_score[j] = -latent[j];
    }

    for (std::size_t i = 0; i < ctx.data.y.size(); ++i) {
      std::size_t j = ctx.data.studies[i];
      double mean = ctx.offset[i] + ctx.between * latent[j];
      for (std::size_t a = 0; a < p; ++a) mean += ctx.coefficient(i, active[a]) * fixed[a];
      double gradient = 0, weight = 0;
      target += ctx.row_score(i, mean, use_selection, fixed_step ? &gradient : nullptr, &weight);
      if (!fixed_step) continue;
      if (!(std::isfinite(gradient) && std::isfinite(weight) && weight > 0)) return JAGS_NEGINF;
      study_precision[j] += ctx.between * (ctx.between * weight);
      study_score[j] += ctx.between * gradient;
      for (std::size_t a = 0; a < p; ++a) {
        double xa = ctx.coefficient(i, active[a]);
        fixed_score[a] += xa * gradient;
        study_solve[j * p + a] += xa * (ctx.between * weight);
        for (std::size_t b = 0; b <= a; ++b) {
          chol[a * p + b] += xa * (ctx.coefficient(i, active[b]) * weight);
        }
      }
    }
    if (!fixed_step || !std::isfinite(target)) return target;

    for (std::size_t j = 0; j < n; ++j) {
      const double h = study_precision[j];
      if (!(h > 0 && std::isfinite(h))) return JAGS_NEGINF;
      for (std::size_t a = 0; a < p; ++a) {
        double c = study_solve[j * p + a] / h;
        fixed_score[a] -= c * study_score[j];
        for (std::size_t b = 0; b <= a; ++b) {
          chol[a * p + b] -= c * study_solve[j * p + b];
        }
      }
      for (std::size_t a = 0; a < p; ++a) study_solve[j * p + a] /= h;
    }
    if (!cholesky(chol, p)) return JAGS_NEGINF;

    for (std::size_t a = 0; a < p; ++a) {
      for (std::size_t b = 0; b < a; ++b) fixed_score[a] -= chol[a * p + b] * fixed_score[b];
      fixed_score[a] /= chol[a * p + a];
    }
    solve_upper(chol, fixed_score);
    for (std::size_t j = 0; j < n; ++j) {
      study_score[j] /= study_precision[j];
      for (std::size_t a = 0; a < p; ++a) study_score[j] -= study_solve[j * p + a] * fixed_score[a];
    }
    *fixed_step = std::move(fixed_score);
    *study_step = std::move(study_score);
    return target;
  }

  bool prepare(Context const &ctx, bool use_selection, unsigned int label)
  {
    Vec bstep, gstep;
    Vec next_beta(beta.size()), next_gamma(gamma.size());
    for (unsigned int iteration = 0; iteration < 40; ++iteration) {
      double target = evaluate(ctx, beta, gamma, use_selection, &bstep, &gstep);
      if (!std::isfinite(target)) return false;
      auto small = [](double step) { return std::isfinite(step) && std::abs(step) < 1e-7; };
      if ((!use_selection && iteration == 1) ||
          (std::all_of(bstep.begin(), bstep.end(), small) &&
           std::all_of(gstep.begin(), gstep.end(), small))) {
        for (std::size_t a = 0; a < beta.size(); ++a) log_det_sqrt += std::log(chol[a * (beta.size() + 1)]);
        for (double h : study_precision) log_det_sqrt += 0.5 * std::log(h);
        log_mass = ctx.data.log_gate[label - 1] + target +
          0.5 * (beta.size() + gamma.size()) * log_two_pi - log_det_sqrt;
        return std::isfinite(log_mass);
      }

      bool improved = false;
      for (double ratio = 1; ratio >= 1e-8; ratio *= 0.5) {
        for (std::size_t a = 0; a < beta.size(); ++a) next_beta[a] = beta[a] + ratio * bstep[a];
        for (std::size_t j = 0; j < gamma.size(); ++j) next_gamma[j] = gamma[j] + ratio * gstep[j];
        if (!use_selection || evaluate(ctx, next_beta, next_gamma, true) >= target - 1e-10) {
          beta.swap(next_beta);
          gamma.swap(next_gamma);
          improved = true;
          break;
        }
      }
      if (!improved) return false;
    }
    return false;
  }

  double log_density(FitData const &data, Vec const &state) const
  {
    double result = log_det_sqrt - 0.5 * (data.precision.size() + gamma.size()) * log_two_pi;
    Vec delta(beta.size());
    for (std::size_t a = 0; a < beta.size(); ++a) delta[a] = state[1 + active[a]] - beta[a];
    for (std::size_t a = 0; a < beta.size(); ++a) {
      double z = 0;
      for (std::size_t b = a; b < beta.size(); ++b) z += chol[b * beta.size() + a] * delta[b];
      result -= 0.5 * z * z;
    }
    for (std::size_t a = 0; a < data.precision.size(); ++a) {
      if (std::find(active.begin(), active.end(), a) != active.end()) continue;
      double z = (state[a + 1] - data.prior_mean[a]) * std::sqrt(data.precision[a]);
      result += 0.5 * (std::log(data.precision[a]) - z * z);
    }
    for (std::size_t j = 0; j < gamma.size(); ++j) {
      double z = state[1 + data.precision.size() + j] - gamma[j];
      for (std::size_t a = 0; a < beta.size(); ++a) z += study_solve[j * beta.size() + a] * delta[a];
      z *= std::sqrt(study_precision[j]);
      result -= 0.5 * z * z;
    }
    return result;
  }

  void draw(FitData const &data, Vec &state, RNG *rng) const
  {
    Vec delta(beta.size());
    for (double &z : delta) z = rng->normal();
    solve_upper(chol, delta);
    for (std::size_t a = 0; a < data.precision.size(); ++a) {
      if (std::find(active.begin(), active.end(), a) == active.end()) {
        state[a + 1] = data.prior_mean[a] + rng->normal() / std::sqrt(data.precision[a]);
      }
    }
    for (std::size_t a = 0; a < beta.size(); ++a) state[1 + active[a]] = beta[a] + delta[a];
    for (std::size_t j = 0; j < gamma.size(); ++j) {
      double center = gamma[j];
      for (std::size_t a = 0; a < beta.size(); ++a) center -= study_solve[j * beta.size() + a] * delta[a];
      state[1 + data.precision.size() + j] = center + rng->normal() / std::sqrt(study_precision[j]);
    }
  }
};

class ScalarRefresh : public Slicer {
  SingletonGraphView const &_view;
  unsigned int _chain;
public:
  ScalarRefresh(SingletonGraphView const &view, unsigned int chain)
    : Slicer(1, 10), _view(view), _chain(chain) {}
  double value() const override { return scalar(_view.node(), _chain); }
  void setValue(double x) override { _view.setValue(&x, 1, _chain); }
  void getLimits(double *lower_, double *upper_) const override
  {
    _view.node()->support(lower_, upper_, 1, _chain);
  }
  double logDensity() const override { return _view.logFullConditional(_chain); }
  void update(RNG *rng) override
  {
    const double previous = value();
    try {
      if (!updateStep(rng) && state() != SLICER_OK) {
        throwNodeError(_view.node(), "Joint cluster scalar refresh has a nonfinite density");
      }
    } catch (...) {
      setValue(previous);
      throw;
    }
  }
};

class JointCluster : public MutableSampleMethod {
  GraphView const &_view;
  const unsigned int _chain;
  std::shared_ptr<const FitData> _data;
  std::vector<std::unique_ptr<SingletonGraphView>> _singletons;
  std::vector<std::unique_ptr<ScalarRefresh>> _slices;
  bool _adapting = true;
  unsigned long _ticks = 0;
  Vec _current, _trial;

  void finite_refresh(RNG *rng)
  {
    auto const &view = *_singletons[0];
    double previous = scalar(view.node(), _chain);
    double labels[2] = {1, 2};
    std::array<double, 2> scores;
    try {
      for (unsigned int k = 0; k < 2; ++k) {
        view.setValue(&labels[k], 1, _chain);
        scores[k] = view.logFullConditional(_chain);
      }
      double maximum = std::max(scores[0], scores[1]);
      if (!std::isfinite(maximum) || std::isnan(scores[0]) || std::isnan(scores[1])) {
        throwNodeError(view.node(), "Joint cluster gate refresh has a nonfinite density");
      }
      double a = std::exp(scores[0] - maximum);
      double b = std::exp(scores[1] - maximum);
      view.setValue(&labels[rng->uniform() * (a + b) < a ? 0 : 1], 1, _chain);
    } catch (...) {
      view.setValue(&previous, 1, _chain);
      throw;
    }
  }

  void moderator_refresh(RNG *rng)
  {
    // Elliptical slice along the declared fixed diagonal Gaussian slab.
    auto const &view = *_singletons[2];
    Vec original(_data->width), trial(_data->width), noise(_data->width);
    view.getValue(original, _chain);
    bool active = scalar(_singletons[0]->node(), _chain) ==
      scalar(_data->setup[moderator_label], _chain);
    for (std::size_t j = 0; j < noise.size(); ++j) {
      noise[j] = rng->normal() / std::sqrt(_data->precision[j + 1]);
      trial[j] = _data->prior_mean[j + 1] + noise[j];
    }
    if (!active) {
      // The factor slab is dormant. Check the likelihood as well so mistaken
      // declarations cannot turn this shortcut into an unconditional redraw.
      double old_score = view.logLikelihood(_chain);
      try {
        view.setValue(trial, _chain);
        double new_score = view.logLikelihood(_chain);
        if (!std::isnan(new_score) && std::log(rng->uniform()) < new_score - old_score) return;
      } catch (...) {
        view.setValue(original, _chain);
        throw;
      }
      view.setValue(original, _chain);
      return;
    }

    constexpr double full_turn = 6.2831853071795865;
    double cutoff = view.logLikelihood(_chain) - rng->exponential();
    double angle = full_turn * rng->uniform();
    double lo = angle - full_turn, hi = angle;
    try {
      for (;;) {
        for (std::size_t j = 0; j < trial.size(); ++j) {
          double center = _data->prior_mean[j + 1];
          trial[j] = center + (original[j] - center) * std::cos(angle) + noise[j] * std::sin(angle);
        }
        view.setValue(trial, _chain);
        if (view.logLikelihood(_chain) >= cutoff) return;
        if (angle < 0) lo = angle;
        else hi = angle;
        double next = lo + (hi - lo) * rng->uniform();
        if (next == angle || next == 0 || next <= lo || next >= hi) break;
        angle = next;
      }
    } catch (...) {
      view.setValue(original, _chain);
      throw;
    }
    view.setValue(original, _chain);
  }

  void metropolis(RNG *rng, double adjustment)
  {
    if (!std::isfinite(adjustment) ||
        !std::all_of(_trial.begin(), _trial.end(), [](double x) { return std::isfinite(x); })) return;
    double before = _view.logFullConditional(_chain);
    try {
      _view.setValue(_trial, _chain);
      double after = _view.logFullConditional(_chain);
      if (std::isfinite(after) && std::log(rng->uniform()) < after - before + adjustment) return;
    } catch (...) {
      _view.setValue(_current, _chain);
      throw;
    }
    _view.setValue(_current, _chain);
  }

  void group_update(RNG *rng)
  {
    Context ctx(*_data, _chain);
    unsigned int on = static_cast<unsigned int>(scalar(_data->setup[moderator_label], _chain));
    std::array<Reference, 2> refs{Reference(ctx, on == 1), Reference(ctx, on == 2)};
    for (unsigned int k = 0; k < 2; ++k) {
      if (!refs[k].prepare(ctx, ctx.selected, k + 1)) {
        // A Gaussian reference is still full support. The original selected
        // density supplies the MH correction. The fallback depends only on
        // the unchanged context and the proposed branch, never the old block.
        refs[k] = Reference(ctx, on == k + 1);
        if (!refs[k].prepare(ctx, false, k + 1)) return;
      }
    }
    _view.getValue(_current, _chain);
    std::size_t from = static_cast<std::size_t>(_current[0]) - 1;
    if (from >= refs.size()) return;
    double max_mass = std::max(refs[0].log_mass, refs[1].log_mass);
    // A small, fixed uniform component maintains numerical support for very
    // different branch masses and makes both labels available from every state.
    double a = std::exp(refs[0].log_mass - max_mass);
    double b = std::exp(refs[1].log_mass - max_mass);
    double p0 = 0.01 + 0.98 * a / (a + b);
    double p[2] = {p0, 1 - p0};
    std::size_t to = rng->uniform() < p[0] ? 0 : 1;
    _trial = _current;
    _trial[0] = static_cast<double>(to + 1);
    refs[to].draw(*_data, _trial, rng);
    double adjustment = std::log(p[from]) + refs[from].log_density(*_data, _current) -
      std::log(p[to]) - refs[to].log_density(*_data, _trial);
    metropolis(rng, adjustment);
  }

  void scale_update(RNG *rng)
  {
    _view.getValue(_current, _chain);
    std::size_t start = 1 + _data->precision.size();
    std::size_t tau_index = start + _data->n_studies;
    std::size_t rho_index = tau_index + 1;
    double tau = _current[tau_index], rho = _current[rho_index];
    if (!(tau > 0 && rho > 0 && rho < 1) ||
        scalar(_data->setup[heterogeneity_gate], _chain) !=
          scalar(_data->setup[heterogeneity_label], _chain)) return;

    long double ss = 0;
    for (std::size_t j = start; j < tau_index; ++j) {
      ss += static_cast<long double>(_current[j]) * _current[j];
    }
    if (!(ss > 0 && std::isfinite(ss))) return;
    double chi = rchisq(static_cast<double>(_data->n_studies), rng);
    double log_change = 0.5 * (std::log(ss) - std::log(chi));
    double new_between = std::exp(std::log(tau) + 0.5 * std::log(rho) + log_change);
    double within = tau * std::sqrt(1 - rho);
    double new_tau = std::hypot(within, new_between);
    double new_rho = new_between / new_tau;
    new_rho *= new_rho;
    if (!(std::isfinite(log_change) && std::isfinite(new_tau) && new_tau > 0 &&
          new_rho > 0 && new_rho < 1)) return;

    _trial = _current;
    _trial[tau_index] = new_tau;
    _trial[rho_index] = new_rho;
    long double new_ss = 0;
    double gamma_factor = std::exp(-log_change);
    for (std::size_t j = start; j < tau_index; ++j) {
      _trial[j] *= gamma_factor;
      new_ss += static_cast<long double>(_trial[j]) * _trial[j];
    }
    // Jacobian of (within, log(between), between * gamma) plus the reverse
    // chi-square proposal density in log(between) units. The GraphView
    // supplies the exact original density, including all stochastic priors.
    double adjustment = 2 * log_change - 3 * (std::log(new_tau) - std::log(tau)) +
      0.5 * (new_ss - ss);
    metropolis(rng, adjustment);
  }

public:
  JointCluster(GraphView const &view, Graph const &graph, unsigned int chain,
               std::shared_ptr<const FitData> data)
    : _view(view), _chain(chain), _data(std::move(data)), _current(view.length()), _trial(view.length())
  {
    for (StochasticNode *node : view.nodes()) {
      _singletons.push_back(std::make_unique<SingletonGraphView>(node, graph));
    }
    for (std::size_t i = 1; i < _singletons.size(); ++i) {
      if (i != 2) _slices.push_back(std::make_unique<ScalarRefresh>(*_singletons[i], chain));
    }
    _view.checkFinite(chain);
  }

  void update(RNG *rng) override
  {
    finite_refresh(rng);
    for (auto &slice : _slices) slice->update(rng);
    moderator_refresh(rng);
    if (!_adapting && ++_ticks == _data->spacing) {
      _ticks = 0;
      group_update(rng);
      if (_data->scale) scale_update(rng);
    }
  }

  bool isAdaptive() const override { return true; }
  bool checkAdaptation() const override
  {
    return std::all_of(_slices.begin(), _slices.end(), [](auto const &slice) {
      return slice->checkAdaptation();
    });
  }
  void adaptOff() override
  {
    for (auto &slice : _slices) slice->adaptOff();
    _adapting = false;
    _ticks = 0;
  }
};

bool collect_setup(Node const *tag, FitData &data)
{
  data.setup = tag->parents();
  std::vector<unsigned long> sizes;
  for (Node const *parent : data.setup) sizes.push_back(parent->length());
  if (!valid_lengths(sizes)) return false;
  for (unsigned int i : {design, groups, response, sampling_se, sampling_variance,
                        observed_bin, lower, upper, telescope, intercept_label,
                        moderator_label, heterogeneity_label, study_count, cadence, scale_enabled}) {
    if (!data.setup[i]->isFixed()) return false;
  }

  if (!integer_between(scalar(data.setup[study_count], 0), 2, sizes[groups]) ||
      !integer_between(scalar(data.setup[cadence], 0), 1, 100) ||
      !integer_between(scalar(data.setup[intercept_label], 0), 1, 2) ||
      !integer_between(scalar(data.setup[moderator_label], 0), 1, 2) ||
      !integer_between(scalar(data.setup[heterogeneity_label], 0), 1, 2) ||
      !integer_between(scalar(data.setup[scale_enabled], 0), 0, 1)) return false;

  data.n_studies = scalar(data.setup[study_count], 0);
  data.spacing = scalar(data.setup[cadence], 0);
  data.scale = scalar(data.setup[scale_enabled], 0) == 1;
  data.width = sizes[design] / sizes[groups];
  if (data.width > sizes[groups] || sizes[omega] > static_cast<unsigned long>(std::numeric_limits<int>::max())) {
    return false;
  }
  auto copy = [&](Vec &to, unsigned int i) {
    auto from = data.setup[i]->value(0);
    to.assign(from, from + sizes[i]);
  };
  copy(data.x, design);
  copy(data.y, response);
  copy(data.se, sampling_se);
  copy(data.variance, sampling_variance);
  SelNormJagsBounds lo(data.setup[lower]->value(0), sizes[lower]);
  SelNormJagsBounds hi(data.setup[upper]->value(0), sizes[upper]);
  data.z_lower.assign(lo.data(), lo.data() + sizes[lower]);
  data.z_upper.assign(hi.data(), hi.data() + sizes[upper]);
  for (double x : data.x) if (!std::isfinite(x)) return false;
  for (unsigned long i = 0; i < sizes[groups]; ++i) {
    double index = data.setup[groups]->value(0)[i];
    double bin = data.setup[observed_bin]->value(0)[i];
    if (!integer_between(index, 1, data.n_studies) || !integer_between(bin, 1, sizes[omega]) ||
        !std::isfinite(data.y[i]) || !(data.se[i] > 0 && std::isfinite(data.se[i])) ||
        !(data.variance[i] > 0 && std::isfinite(data.variance[i]))) return false;
    data.studies.push_back(static_cast<unsigned long>(index) - 1);
    data.bins.push_back(static_cast<int>(bin));
  }
  return true;
}

struct Member {
  StochasticNode *node;
  Role role;
  unsigned long slot;
};

bool collect_nodes(std::vector<Member> const &members, FitData &data,
                   std::vector<StochasticNode *> &nodes)
{
  std::array<StochasticNode *, 7> single{};
  std::vector<StochasticNode *> gamma(data.n_studies, nullptr);
  std::set<StochasticNode *> distinct;
  for (auto const &m : members) {
    if (!distinct.insert(m.node).second) return false;
    for (unsigned long i = 0; i < m.node->length(); ++i) {
      if (m.node->isObserved(i)) return false;
    }
    if (m.role == study) {
      if (m.slot > gamma.size() || gamma[m.slot - 1]) return false;
      gamma[m.slot - 1] = m.node;
    } else {
      if (m.slot != 1 || single[m.role]) return false;
      single[m.role] = m.node;
    }
  }
  if (!single[gate] || !single[intercept] || !single[moderator] ||
      single[gate] != data.setup[group_gate] || (single[total_sd] != nullptr) != data.scale ||
      (single[allocation] != nullptr) != data.scale ||
      std::any_of(gamma.begin(), gamma.end(), [](auto node) { return node == nullptr; })) return false;

  auto gate_node = single[gate];
  auto &parents = gate_node->parents();
  if (gate_node->distribution()->name() != "dcat" || gate_node->length() != 1 ||
      gate_node->lowerBound() || gate_node->upperBound() ||
      parents.size() != 1 || parents[0]->length() != 2 || !parents[0]->isFixed()) return false;
  double a = parents[0]->value(0)[0], b = parents[0]->value(0)[1];
  if (!(a > 0 && b > 0 && std::isfinite(a + b))) return false;
  data.log_gate = {std::log(a) - std::log(a + b), std::log(b) - std::log(a + b)};
  if (single[intercept]->length() != 1 || single[moderator]->length() != data.width ||
      !normal_parameters(single[intercept], data.prior_mean, data.precision) ||
      !normal_parameters(single[moderator], data.prior_mean, data.precision)) return false;

  nodes = {gate_node, single[intercept], single[moderator]};
  for (auto node : gamma) {
    Vec mean, precision;
    if (node->length() != 1 || !normal_parameters(node, mean, precision) ||
        mean[0] != 0 || precision[0] != 1) return false;
    nodes.push_back(node);
  }
  if (data.scale) {
    Vec mean, precision;
    if (!normal_parameters(single[total_sd], mean, precision, true) ||
        !beta_parameters(single[allocation])) return false;
    nodes.push_back(single[total_sd]);
    nodes.push_back(single[allocation]);
  }
  for (auto node : nodes) {
    for (unsigned int i : {intercept_gate, heterogeneity_gate, omega, mode, pet, peese}) {
      if (depends_on(data.setup[i], node)) return false;
    }
    if (node != single[total_sd] && node != single[allocation] &&
        (depends_on(data.setup[physical_tau], node) || depends_on(data.setup[physical_rho], node))) return false;
  }
  return true;
}

class JointClusterFactory : public SamplerFactory {
public:
  std::string name() const override { return sampler_name; }

  std::vector<Sampler *> makeSamplers(std::list<StochasticNode *> const &free,
                                      Graph const &graph) const override
  {
    std::vector<Member> members;
    Node const *setup = nullptr;
    for (StochasticNode *node : free) {
      for (Node const *tag : direct_terminal_tags(node, setup_name, setup_count, graph)) {
        if (setup) return {};
        setup = tag;
      }
      for (Node const *tag : direct_terminal_tags(node, member_name, 3, graph)) {
        auto const &p = tag->parents();
        if (!fixed_scalar(p[1]) || !fixed_scalar(p[2]) ||
            !integer_between(scalar(p[1], 0), gate, allocation) ||
            !integer_between(scalar(p[2], 0), 1, free.size())) return {};
        members.push_back({node, static_cast<Role>(static_cast<int>(scalar(p[1], 0))),
                           static_cast<unsigned long>(scalar(p[2], 0))});
      }
    }
    if (!setup) return {};
    auto data = std::make_shared<FitData>();
    std::vector<StochasticNode *> nodes;
    if (!collect_setup(setup, *data) || !collect_nodes(members, *data, nodes)) return {};
    for (StochasticNode *node : nodes) {
      for (StochasticNode *child : *node->stochasticChildren()) {
        auto const &name = child->distribution()->name();
        if (name == "dbt_gate" || name == "dbt_block") return {};
      }
    }

    auto view = std::make_unique<GraphView>(nodes, graph);
    std::size_t observations = 0;
    for (StochasticNode const *child : view->stochasticChildren()) {
      auto const &name = child->distribution()->name();
      if (!isObserved(child) || child->length() != 1 ||
          (name != "dnorm" && name != "dselnorm_step_switch")) return {};
      ++observations;
    }
    if (observations != data->y.size()) return {};

    std::vector<std::unique_ptr<MutableSampleMethod>> methods;
    std::vector<MutableSampleMethod *> raw_methods;
    for (unsigned int chain = 0; chain < nodes[0]->nchain(); ++chain) {
      methods.push_back(std::make_unique<JointCluster>(*view, graph, chain, data));
      raw_methods.push_back(methods.back().get());
    }
    std::string name(sampler_name);
    std::unique_ptr<Sampler> sampler(new MutableSampler(view.release(), raw_methods, name));
    for (auto &method : methods) method.release();
    std::vector<Sampler *> out{sampler.get()};
    sampler.release();
    return out;
  }
};
} // namespace

ClusterTag::ClusterTag(bool setup)
  : JointSamplerTag(setup ? setup_name : member_name, setup ? static_cast<unsigned long>(setup_count) : 3),
    _setup(setup) {}

bool ClusterTag::checkParameterLength(std::vector<unsigned long> const &sizes) const
{
  return _setup ? valid_lengths(sizes) : sizes.size() == 3 && sizes[0] >= 1 && sizes[1] == 1 && sizes[2] == 1;
}

SamplerFactory *joint_cluster_factory()
{
  static JointClusterFactory factory;
  return &factory;
}

} // namespace RoBMA
} // namespace jags
