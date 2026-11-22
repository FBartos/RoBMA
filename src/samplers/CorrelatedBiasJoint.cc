#include "CorrelatedBiasJoint.h"
#include "CoarseSurrogate.h"
#include "CorrectedSampler.h"
#include "IndependentBiasProposal.h"

#include "../distributions/selnorm-jags-bounds.h"

#include <BayesTools/sampler_declarations.h>

#include <graph/Graph.h>
#include <graph/StochasticNode.h>
#include <rng/RNG.h>
#include <sampler/GraphView.h>
#include <sampler/Sampler.h>
#include <sampler/SingletonGraphView.h>
#include <util/nainf.h>

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace jags::RoBMA {
namespace {

namespace declarations = BayesTools::sampler_declarations;
using joint_bias::Family;

constexpr char marker_name[] = "drobma_joint_correlated";
constexpr char sampler_name[] = "RoBMA::CorrelatedBiasJoint";

// The marker is declared by sampler-joint-correlated.R. Branch arrays use the
// original dcat order, omega is branch-major, and geometry matrices are column-major.
enum Parent {
  bias, raw_intercept, effect_gate, active_label, within_sd, retained_offset,
  moderator, moderator_design, family, slopes, omega, row_data,
  eigenvalues, eigenvectors, lower, upper, parent_count
};

bool fixed_scalar(Node const *node)
{
  return node->isFixed() && node->length() == 1;
}

unsigned long categorical_size(StochasticNode const *node)
{
  const auto &parents = node->parents();
  if (node->length() != 1 || node->distribution()->name() != "dcat" ||
      parents.size() != 1 || !isSupportFixed(node) || !parents.front()->isFixed()) return 0;
  const double *p = parents.front()->value(0);
  const auto n = parents.front()->length();
  if (!std::all_of(p, p + n, [](double w) { return std::isfinite(w) && w > 0; }) ||
      !std::isfinite(std::accumulate(p, p + n, 0.0))) return 0;
  return n;
}

struct Normal {
  double mean = 0;
  double sd = 1;

  double draw(RNG *rng) const { return mean + sd * rng->normal(); }

  double logDensity(double x) const { return joint_bias::log_normal(x, mean, sd); }
};

struct Reference {
  Normal base, shifted;
  bool mixture = false;

  double draw(RNG *rng) const
  {
    return (mixture && rng->uniform() < 0.85 ? shifted : base).draw(rng);
  }

  double logDensity(double x) const
  {
    const double a = base.logDensity(x);
    if (!mixture) return a;
    const double b = shifted.logDensity(x);
    const double top = std::max(a, b);
    return std::isfinite(top) ?
      top + std::log(0.15 * std::exp(a - top) + 0.85 * std::exp(b - top)) : top;
  }
};

struct Geometry {
  Normal prior;
  double prior_precision = 0;
  std::vector<Family> families;
  std::vector<double> cuts, yi, sei, variance;
  std::vector<double> lambda, basis, projected_one, projected_y, projected_se, projected_se2;
  std::vector<double> design;
};

bool normal_prior(StochasticNode const *node, Geometry &geometry)
{
  const auto &parents = node->parents();
  if (node->length() != 1 || node->distribution()->name() != "dnorm" ||
      !node->fullRank() || parents.size() != 2 || !isSupportFixed(node) ||
      !fixed_scalar(parents[0]) || !fixed_scalar(parents[1])) return false;
  double lo, hi;
  node->support(&lo, &hi, 1, 0);
  const double mean = parents[0]->value(0)[0];
  const double precision = parents[1]->value(0)[0];
  if (lo != JAGS_NEGINF || hi != JAGS_POSINF || !std::isfinite(mean) ||
      !std::isfinite(precision) || !(precision > 0)) return false;
  geometry.prior = {mean, 1 / std::sqrt(precision)};
  geometry.prior_precision = precision;
  return true;
}

bool guard_children(StochasticNode const *marker, const GraphView &view, unsigned long n)
{
  const auto &owned = view.nodes();
  unsigned long observed = 0;
  for (StochasticNode const *child : view.stochasticChildren()) {
    if (child == marker) continue;
    const auto &name = child->distribution()->name();
    const auto &parents = child->parents();
    if (name == declarations::gate_marker_name()) {
      // The bias label controls independent nuisance samplers' activity. An
      // actual declaration for either owned node is an ownership conflict.
      if (!declarations::observed_zero(child) || parents.size() != 2 ||
          !dynamic_cast<StochasticNode const *>(parents[0]) || parents[1]->length() != 1 ||
          view.isDependent(parents[0]) || declarations::depends_on(parents[1], owned[1])) return false;
      continue;
    }
    // A foreign block/readout/owner marker (even malformed) also fails this
    // allowlist. No sampling child or unobserved row should be shared.
    if (!isObserved(child) || parents.empty() || !view.isDependent(parents.front())) return false;
    unsigned int sign_index;
    if (name == "dselnorm_step_switch") sign_index = 8;
    else if (name == "dselnorm_mnorm_step") sign_index = 7;
    else if (name == "dselnorm_cluster_step" || name == "dselnorm_factor_step") sign_index = 8;
    else return false;
    if (parents.size() <= sign_index || !fixed_scalar(parents[sign_index]) ||
        parents[sign_index]->value(0)[0] != 1 || child->length() > n - observed) return false;
    observed += child->length();
  }
  return observed == n;
}

bool configure(StochasticNode const *marker, const GraphView &view, Geometry &geometry)
{
  const auto &par = marker->parents();
  const auto &roots = view.nodes();
  const auto *gate = dynamic_cast<StochasticNode const *>(par[effect_gate]);
  if (!gate || categorical_size(gate) != 2 || !fixed_scalar(par[active_label]) ||
      (par[active_label]->value(0)[0] != 1 && par[active_label]->value(0)[0] != 2) ||
      !normal_prior(roots[1], geometry)) return false;

  // Check actual ancestry through stochastic nodes as well as deterministic
  // descendants. Every live proposal input must be fixed during this move.
  for (auto role : {effect_gate, within_sd, retained_offset, moderator, slopes, omega}) {
    for (const auto *root : roots) {
      if (declarations::depends_on(par[role], root)) return false;
    }
  }
  if (par[within_sd]->length() != 1) return false;
  for (auto role : {moderator_design, family, row_data, eigenvalues, eigenvectors, lower, upper}) {
    if (!par[role]->isFixed()) return false;
  }

  const auto n_branches = categorical_size(roots[0]);
  const auto n_bins = par[lower]->length();
  const auto n = par[eigenvalues]->length();
  if (n_branches < 2 || n_branches > static_cast<unsigned long>(std::numeric_limits<int>::max()) ||
      par[family]->length() != n_branches || par[slopes]->length() != n_branches ||
      n_bins < 2 || par[upper]->length() != n_bins ||
      par[omega]->length() % n_branches || par[omega]->length() / n_branches != n_bins ||
      n < 2 || par[retained_offset]->length() != n ||
      par[row_data]->length() % n || par[row_data]->length() / n != 3 ||
      par[eigenvectors]->length() % n || par[eigenvectors]->length() / n != n ||
      par[moderator_design]->length() % n ||
      par[moderator_design]->length() / n != par[moderator]->length() ||
      !guard_children(marker, view, n)) return false;

  const double *codes = par[family]->value(0);
  for (unsigned long b = 0; b < n_branches; ++b) {
    if (!selnorm_jags_integer_in_range(codes[b], 0, 3)) return false;
    geometry.families.push_back(static_cast<Family>(static_cast<int>(codes[b])));
  }
  const auto selected = std::count(geometry.families.begin(), geometry.families.end(), Family::selection);
  if (selected == 0 || static_cast<unsigned long>(selected) == n_branches) return false;

  const double *z_lower = par[lower]->value(0), *z_upper = par[upper]->value(0);
  geometry.cuts.push_back(JAGS_POSINF);
  for (unsigned long j = 0; j < n_bins; ++j) {
    const double lo = selnorm_jags_bound_to_ieee(z_lower[j]);
    const double hi = selnorm_jags_bound_to_ieee(z_upper[j]);
    if (!(lo < hi) || hi != geometry.cuts.back()) return false;
    geometry.cuts.push_back(lo);
  }
  if (geometry.cuts.back() != JAGS_NEGINF) return false;

  const double *rows = par[row_data]->value(0);
  geometry.yi.assign(rows, rows + n);
  geometry.sei.assign(rows + n, rows + 2 * n);
  geometry.variance.assign(rows + 2 * n, rows + 3 * n);

  const auto copy_fixed = [&par](Parent role, std::vector<double> &out) {
    const double *x = par[role]->value(0);
    out.assign(x, x + par[role]->length());
  };
  copy_fixed(eigenvalues, geometry.lambda);
  copy_fixed(eigenvectors, geometry.basis);
  copy_fixed(moderator_design, geometry.design);
  for (auto values : {&geometry.yi, &geometry.basis, &geometry.design}) {
    if (!std::all_of(values->begin(), values->end(), [](double x) { return std::isfinite(x); })) return false;
  }
  for (auto values : {&geometry.sei, &geometry.variance, &geometry.lambda}) {
    if (!std::all_of(values->begin(), values->end(), [](double x) { return std::isfinite(x) && x > 0; })) return false;
  }

  for (unsigned long j = 0; j < n; ++j) {
    const double *u = geometry.basis.data() + j * n;
    geometry.projected_one.push_back(std::accumulate(u, u + n, 0.0));
    geometry.projected_y.push_back(std::inner_product(u, u + n, geometry.yi.begin(), 0.0));
    geometry.projected_se.push_back(std::inner_product(u, u + n, geometry.sei.begin(), 0.0));
    double se2 = 0;
    for (unsigned long i = 0; i < n; ++i) se2 += u[i] * geometry.sei[i] * geometry.sei[i];
    geometry.projected_se2.push_back(se2);
  }
  return true;
}

class CorrelatedBiasJoint : public Sampler {
  struct Chain {
    std::vector<double> offset, mass, derivative;
    std::vector<double> score, candidate, target, weight;
    std::vector<Reference> references;

    explicit Chain(const Geometry &geometry)
      : offset(geometry.yi.size()), mass(geometry.cuts.size() - 1), derivative(mass.size()),
        score(geometry.families.size()), candidate(score.size()), target(score.size()),
        weight(score.size()), references(score.size()) {}
  };

  const GraphView *_view;
  const std::vector<Node const *> _roles;
  const Geometry _geometry;
  std::vector<Chain> _chains;

  void propose(unsigned int chain, Chain &work) const
  {
    const auto &g = _geometry;
    const auto n = g.yi.size();
    const double sd = _roles[within_sd]->value(chain)[0];
    const double variance = sd * sd;
    if (!std::isfinite(variance) || sd < 0) return;
    const double *offset = _roles[retained_offset]->value(chain);
    const double *beta = _roles[moderator]->value(chain);
    std::copy(offset, offset + n, work.offset.begin());
    for (unsigned long j = 0; j < _roles[moderator]->length(); ++j) {
      for (std::size_t i = 0; i < n; ++i) work.offset[i] += beta[j] * g.design[i + n * j];
    }

    // In the rotated basis the exact Gaussian reference covariance is
    // diagonal: V + within^2 I. The study effects stay in the mean.
    double precision = g.prior_precision, location = precision * g.prior.mean;
    double pet = 0, peese = 0;
    for (std::size_t j = 0; j < n; ++j) {
      const double *u = g.basis.data() + j * n;
      const double residual = g.projected_y[j] -
        std::inner_product(u, u + n, work.offset.begin(), 0.0);
      const double info = g.projected_one[j] / (g.lambda[j] + variance);
      precision += info * g.projected_one[j];
      location += info * residual;
      pet += info * g.projected_se[j];
      peese += info * g.projected_se2[j];
    }
    const Normal base{location / precision, 1 / std::sqrt(precision)};
    if (!std::isfinite(base.mean) || !std::isfinite(base.sd) || !(base.sd > 0)) return;

    // The row-marginal normalizer derivative only changes the reference.
    // Invalid proposal geometry keeps the full-support original prior reference.
    std::fill(work.score.begin(), work.score.end(), 0);
    const double *weights = _roles[omega]->value(chain);
    const auto bins = work.mass.size();
    for (std::size_t i = 0; i < n; ++i) {
      const double row_sd = std::sqrt(g.variance[i] + variance);
      const double mean = base.mean + work.offset[i];
      if (!std::isfinite(mean) || !std::isfinite(row_sd)) return;
      auto hi = joint_bias::boundary(JAGS_POSINF);
      for (std::size_t j = 0; j < bins; ++j) {
        const auto lo = joint_bias::boundary((g.cuts[j + 1] * g.sei[i] - mean) / row_sd);
        work.mass[j] = joint_bias::interval_area(lo, hi);
        work.derivative[j] = (lo.density - hi.density) / row_sd;
        hi = lo;
      }
      for (std::size_t b = 0; b < g.families.size(); ++b) {
        if (g.families[b] != Family::selection) continue;
        const double *w = weights + b * bins;
        const double mass = std::inner_product(w, w + bins, work.mass.begin(), 0.0);
        const double derivative = std::inner_product(w, w + bins, work.derivative.begin(), 0.0);
        if (mass > 0 && std::isfinite(derivative / mass)) work.score[b] -= derivative / mass;
      }
    }

    const double *slope = _roles[slopes]->value(chain);
    for (std::size_t b = 0; b < g.families.size(); ++b) {
      if (g.families[b] == Family::selection) {
        // Smoothly bound the shift to three base SDs.
        const double standardized = work.score[b] * base.sd;
        const double shift = standardized / (1 + std::abs(standardized) / 3) * base.sd;
        const Normal shifted{base.mean + shift, 1.25 * base.sd};
        if (std::isfinite(shifted.mean) && std::isfinite(shifted.sd)) {
          work.references[b] = {base, shifted, true};
        }
      } else {
        Normal gaussian = base;
        if (g.families[b] == Family::pet) gaussian.mean -= pet / precision * slope[b];
        if (g.families[b] == Family::peese) gaussian.mean -= peese / precision * slope[b];
        if (std::isfinite(gaussian.mean)) work.references[b] = {gaussian, {}, false};
      }
    }
  }

  void set(double label, double raw, unsigned int chain) const
  {
    const std::array<double, 2> value{label, raw};
    _view->setValue(value.data(), value.size(), chain);
  }

  double logTarget(unsigned int chain) const
  {
    const double log_target = _view->logFullConditional(chain);
    if (std::isnan(log_target) || log_target == JAGS_POSINF) {
      throw std::runtime_error("Correlated bias sampler could not evaluate its full conditional.");
    }
    return log_target;
  }

public:
  CorrelatedBiasJoint(std::unique_ptr<GraphView> view,
                      std::vector<Node const *> roles, Geometry geometry)
    : Sampler(view.get()), _view(view.release()), _roles(std::move(roles)), _geometry(std::move(geometry))
  {
    const auto chains = _view->nodes()[0]->nchain();
    _chains.reserve(chains);
    for (unsigned int chain = 0; chain < chains; ++chain) _chains.emplace_back(_geometry);
  }

  void update(unsigned int chain, RNG *rng) override
  {
    auto &work = _chains[chain];
    const auto &roots = _view->nodes();
    const double label = roots[0]->value(chain)[0], raw = roots[1]->value(chain)[0];
    const auto n = _geometry.families.size();
    const double gate = _roles[effect_gate]->value(chain)[0];
    const bool active = gate == _roles[active_label]->value(0)[0];
    if (!selnorm_jags_integer_in_range(label, 1, static_cast<int>(n)) || !std::isfinite(raw) ||
        (gate != 1 && gate != 2)) {
      throw std::runtime_error("Correlated bias sampler encountered an invalid current state.");
    }
    const auto incumbent = static_cast<std::size_t>(label - 1);
    std::fill(work.references.begin(), work.references.end(), Reference{_geometry.prior, {}, false});
    if (active) propose(chain, work);

    try {
      // When active, the temporary joint law pairs pi(b,x_b) with the
      // independent q_j auxiliary draws for every j != b. Its finite update
      // has mass pi(b,x_b) / q_b(x_b). When dormant, draw finite weights at x.
      for (std::size_t b = 0; b < n; ++b) {
        const auto &q = work.references[b];
        work.candidate[b] = active && b != incumbent ? q.draw(rng) : raw;
        if (!std::isfinite(work.candidate[b])) {
          throw std::runtime_error("Correlated bias sampler generated a non-finite reference draw.");
        }
        set(static_cast<double>(b + 1), work.candidate[b], chain);
        work.target[b] = logTarget(chain);
        work.weight[b] = work.target[b];
        if (active && work.target[b] != JAGS_NEGINF) {
          work.weight[b] -= q.logDensity(work.candidate[b]);
        }
        if (std::isnan(work.weight[b]) || work.weight[b] == JAGS_POSINF) {
          throw std::runtime_error("Correlated bias sampler could not evaluate an ensemble weight.");
        }
      }
      if (!std::isfinite(work.weight[incumbent])) {
        throw std::runtime_error("Correlated bias sampler requires a finite current weight.");
      }
      const double max_weight = *std::max_element(work.weight.begin(), work.weight.end());
      for (double &value : work.weight) value = std::exp(value - max_weight);
      std::partial_sum(work.weight.begin(), work.weight.end(), work.weight.begin());
      const double u = rng->uniform() * work.weight.back();
      auto pick = std::upper_bound(work.weight.begin(), work.weight.end(), u);
      if (pick == work.weight.end()) {
        pick = std::lower_bound(work.weight.begin(), work.weight.end(), work.weight.back());
      }
      const auto chosen = static_cast<std::size_t>(pick - work.weight.begin());
      const double chosen_label = static_cast<double>(chosen + 1);
      const double chosen_raw = work.candidate[chosen];
      const auto &q = work.references[chosen];

      // Always refresh the chosen location. This independence MH step reuses
      // its exact ensemble score. On the spike, q is the raw prior and the
      // correction also protects against any unknown dependence on the slab.
      const double refreshed = q.draw(rng);
      if (std::isfinite(refreshed)) {
        set(chosen_label, refreshed, chain);
        const double log_target = logTarget(chain);
        if (log_target != JAGS_NEGINF) {
          const double ratio = (log_target - work.target[chosen]) +
            (q.logDensity(chosen_raw) - q.logDensity(refreshed));
          if (std::isnan(ratio)) {
            throw std::runtime_error("Correlated bias sampler could not evaluate its refresh acceptance.");
          }
          if (std::log(rng->uniform()) < ratio) return;
        }
      }
      set(chosen_label, chosen_raw, chain);
    } catch (...) {
      set(label, raw, chain);
      throw;
    }
  }

  void adaptOff() override {}
  bool checkAdaptation() const override { return true; }
  bool isAdaptive() const override { return false; }
  std::string name() const override { return sampler_name; }
};

} // namespace

DCorrelatedBiasJointMarker::DCorrelatedBiasJointMarker() : VectorDist(marker_name, parent_count) {}

double DCorrelatedBiasJointMarker::logDensity(
    double const *x, PDFType, std::vector<double const *> const &,
    std::vector<unsigned long> const &) const
{
  return *x == 0 ? 0 : JAGS_NEGINF;
}

void DCorrelatedBiasJointMarker::randomSample(double *x, std::vector<double const *> const &,
    std::vector<unsigned long> const &, RNG *) const
{
  *x = 0;
}

void DCorrelatedBiasJointMarker::support(double *lo, double *hi,
    std::vector<double const *> const &, std::vector<unsigned long> const &) const
{
  *lo = *hi = 0;
}

bool DCorrelatedBiasJointMarker::checkParameterLength(std::vector<unsigned long> const &len) const
{
  return len.size() == parent_count && std::all_of(len.begin(), len.begin() + within_sd + 1,
    [](unsigned long n) { return n == 1; });
}

bool DCorrelatedBiasJointMarker::checkParameterValue(std::vector<double const *> const &,
    std::vector<unsigned long> const &) const { return true; }
bool DCorrelatedBiasJointMarker::checkParameterDiscrete(std::vector<bool> const &) const { return true; }
bool DCorrelatedBiasJointMarker::isDiscreteValued(std::vector<bool> const &) const { return true; }
bool DCorrelatedBiasJointMarker::isSupportFixed(std::vector<bool> const &) const { return true; }
unsigned long DCorrelatedBiasJointMarker::length(std::vector<unsigned long> const &) const { return 1; }

std::vector<Sampler *> CorrelatedBiasJointFactory::makeSamplers(
    std::list<StochasticNode *> const &nodes, Graph const &graph) const
{
  const auto config = corrected_sampler_configuration();
  if (!config.enabled) return {};

  for (auto *node : nodes) {
    if (node->distribution()->name() != "dcat") continue;
    const auto markers = declarations::markers_of(node, marker_name);
    if (markers.size() != 1 || !graph.contains(markers[0]) ||
        !declarations::observed_zero(markers[0])) continue;
    const auto &roles = markers[0]->parents();
    if (roles.size() != parent_count) continue;
    auto raw = std::find(nodes.begin(), nodes.end(), roles[raw_intercept]);
    if (raw == nodes.end() || *raw == node || (*raw)->distribution()->name() != "dnorm" ||
        declarations::has_sampler_marker_child(*raw)) continue;

    // Replace only the known corrected scalar owner. Use its same eligibility
    // predicate on the actual likelihood children, without creating a second owner.
    const SingletonGraphView scalar(*raw, graph);
    const auto &children = scalar.stochasticChildren();
    if (std::none_of(children.begin(), children.end(), [&config](const auto *child) {
      return child_kind(child, config.settings) != ChildKind::other;
    })) continue;

    auto view = std::make_unique<GraphView>(std::vector<StochasticNode *>{node, *raw}, graph, true);
    Geometry geometry;
    if (!configure(markers[0], *view, geometry)) continue;
    auto sampler = std::make_unique<CorrelatedBiasJoint>(std::move(view), roles, std::move(geometry));
    std::vector<Sampler *> result{sampler.get()};
    sampler.release();
    return result;
  }
  return {};
}

std::string CorrelatedBiasJointFactory::name() const { return sampler_name; }

} // namespace jags::RoBMA
