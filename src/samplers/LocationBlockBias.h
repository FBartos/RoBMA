// The bias-label members of the location block (RoBMA::LocationBlock): the
// declaration, the verified structure and the numerical helpers of the bias
// indicator and its PET and PEESE slope sources (design: wp3-spec.md of the
// workspace). A source fragment of LocationBlock.cc, included there after the
// block's own structure types, inside its anonymous namespace; it is not a
// stand-alone header.
//
// R declares the bias members with one terminal tag on the indicator:
//   robma_location_bias(indicator, codes, PET source, PEESE source, intercept slot)
// codes gives each label's family (0 none, 1 step selection, 2 PET, 3 PEESE);
// an absent slope passes a fixed scalar. The factory verifies on the compiled
// graph (verify_bias, probe_bias): the indicator is a dcat with fixed
// probabilities; each slope source is a scalar continuous node with fixed
// prior parameters and support exactly [0, Inf), entering the graph only
// through its chart physical = source * equals(indicator, label) and that
// chart only through the children's means, linearly with a fixed per-row
// coefficient (probed); the children's kernel mode is NORMAL exactly on the
// labels coded none, PET and PEESE; every other stochastic descendant of the
// indicator is a BayesTools gate marker or the nested weight node
// (BayesTools_nested_a, a beta child whose parameters depend on the label
// only); the label enters the children otherwise only through parameters the
// block's Gaussian working model does not read (weights, kernel mode, rules).

constexpr char bias_tag_name[] = "robma_location_bias";

// A PET or PEESE slope source and its chart.
struct SlopeSpec {
  StochasticNode *source = nullptr;   // null: the label family is absent
  Node const *physical = nullptr;     // source * equals(indicator, label)
  int label = 0;                      // the indicator level that activates the slope
  unsigned pos = 0;                   // position in the block's member values
  ScalarDist const *dist = nullptr;   // the source's prior (with its bounds)
  Vec column;                         // per row: mean change per unit physical slope (probed)
  std::vector<unsigned> children;     // children with a nonzero column entry
};

struct BiasSpec {
  StochasticNode *indicator = nullptr;
  unsigned pos = 0;                   // position of the indicator in the member values
  int levels = 0;                     // labels 1..levels
  Vec log_prior;                      // normalized log prior mass of each label (index label - 1)
  std::vector<int> family;            // 0 none, 1 step selection, 2 PET, 3 PEESE
  std::vector<int> slope_of;          // per label: -1, or the active slope (0 PET, 1 PEESE)
  SlopeSpec slope[2];                 // PET, PEESE
  unsigned intercept = 0;             // the intercept group (scalar)
  // the deterministic descendants of the indicator outside the slope charts,
  // in the view's (topological) order: recomputed after a raw write of the
  // indicator so that the children's weights, kernel modes and rules and the
  // nested node's parameters follow the written label
  std::vector<DeterministicNode *> label_nodes;
  StochasticNode *nested = nullptr;   // the nested weight child, if any
  std::set<StochasticNode const *> skip;  // gate markers and the nested node: not likelihood children

  bool gaussian(int label) const { return family[static_cast<std::size_t>(label - 1)] != 1; }
};

// The raw declaration of collect().
struct BiasDeclared {
  Node const *tag = nullptr;
  StochasticNode *indicator = nullptr;
  std::vector<int> codes;
  StochasticNode *sources[2] = {nullptr, nullptr};
  unsigned intercept_slot = 0;
};

// ---- numerical helpers

// log Phi(x), tail-safe.
inline double log_phi_upper(double x) { return pnorm(x, 0.0, 1.0, 1, 1); }

// A draw from N(mean, sd^2) truncated to [0, Inf): rejection from the normal
// when the mass above zero is at least 0.3, otherwise Robert's (1995)
// exponential rejection sampler for the tail beyond a = -mean / sd > 0.
inline double truncated_normal_positive(double mean, double sd, RNG *rng)
{
  const double a = -mean / sd;
  if (a < 0.5) {
    for (int attempt = 0; attempt < 1000; ++attempt) {
      const double x = mean + sd * rng->normal();
      if (x >= 0) return x;
    }
  }
  const double alpha = 0.5 * (a + std::sqrt(a * a + 4.0));
  for (int attempt = 0; attempt < 100000; ++attempt) {
    const double z = a + rng->exponential() / alpha;
    const double rho = std::exp(-0.5 * (z - alpha) * (z - alpha));
    if (rng->uniform() <= rho) return mean + sd * z;
  }
  return std::numeric_limits<double>::quiet_NaN();
}

// log density of N(mean, sd^2) truncated to [0, Inf) at x (-Inf below zero).
inline double log_truncated_normal_positive(double x, double mean, double sd)
{
  if (!(x >= 0)) return -std::numeric_limits<double>::infinity();
  const double z = (x - mean) / sd;
  return -0.5 * z * z - 0.91893853320467274178 - std::log(sd) - log_phi_upper(mean / sd);
}

// ---- the slope prior

struct SlopePrior {
  ScalarDist const *dist = nullptr;
  std::vector<double const *> par;
  double const *lower = nullptr, *upper = nullptr;
};

inline SlopePrior slope_prior(SlopeSpec const &spec, unsigned chain)
{
  SlopePrior out;
  out.dist = spec.dist;
  for (Node const *parent : spec.source->parents()) out.par.push_back(parent->value(chain));
  if (spec.source->lowerBound() != nullptr) out.lower = spec.source->lowerBound()->value(chain);
  if (spec.source->upperBound() != nullptr) out.upper = spec.source->upperBound()->value(chain);
  return out;
}

// The normalized prior log density (PDF_FULL: truncation constants included).
inline double slope_log_prior(SlopePrior const &prior, double x)
{
  return prior.dist->logDensity(x, PDF_FULL, prior.par, prior.lower, prior.upper);
}

// ---- structure (factory)

// equals(indicator, k) with a fixed integer k: k, or 0.
inline int indicator_test(Node const *node, StochasticNode const *indicator)
{
  auto const *logical = dynamic_cast<LogicalNode const *>(node);
  if (logical == nullptr || node->parents().size() != 2 || node->parents()[0] != indicator) return 0;
  const std::string text = logical->deparse({"a", "b"});
  if (text != "equals(a,b)" && text != "(equals(a,b))") return 0;
  Node const *k = node->parents()[1];
  if (!fixed_scalar(k)) return 0;
  const double v = k->value(0)[0];
  if (!(v >= 1 && v == std::floor(v) && v <= 1e6)) return 0;
  return static_cast<int>(v);
}

// The chart physical = source * equals(indicator, k) in either order: k, or 0.
inline int slope_chart(Node const *physical, StochasticNode const *source, StochasticNode const *indicator)
{
  auto const *logical = dynamic_cast<LogicalNode const *>(physical);
  if (logical == nullptr || physical->length() != 1 || physical->parents().size() != 2 ||
      logical->deparse({"a", "b"}) != "(a*b)") return 0;
  for (unsigned order = 0; order < 2; ++order) {
    if (physical->parents()[order] != source) continue;
    const int k = indicator_test(physical->parents()[1 - order], indicator);
    if (k > 0) return k;
  }
  return 0;
}

// A slope source: scalar, continuous, a scalar distribution with fixed
// parameters and support exactly [0, Inf), unobserved.
inline bool slope_source_law(StochasticNode const *node)
{
  if (node->length() != 1 || node->isDiscreteValued() || !isSupportFixed(node) || !unobserved(node) ||
      dynamic_cast<ScalarDist const *>(node->distribution()) == nullptr) return false;
  for (Node const *parent : node->parents()) if (!parent->isFixed()) return false;
  if (node->lowerBound() != nullptr && !node->lowerBound()->isFixed()) return false;
  if (node->upperBound() != nullptr && !node->upperBound()->isFixed()) return false;
  double lower, upper;
  node->support(&lower, &upper, 1, 0);
  return lower == 0.0 && upper == std::numeric_limits<double>::infinity();
}

// The parameter slots of a likelihood child that the bias label may enter:
// none of those the block's Gaussian working model and row models read (mean,
// variances and loadings, standard errors, power, bounds, observed bin, sign,
// telescoping); the weights, kernel modes and rules are read from the graph at
// each label.
inline bool label_slot_allowed(ChildType type, std::size_t slot)
{
  switch (type) {
  case ChildType::normal:
    return false;
  case ChildType::step_switch:
    return slot == 4 || slot == 9;
  case ChildType::cluster:
  case ChildType::factor:
    return slot != 0 && !(slot >= 1 && slot <= 3) && !(slot >= 5 && slot <= 8);
  case ChildType::mv:
    return slot != 0 && slot != 1 && slot != 2 && !(slot >= 4 && slot <= 7);
  }
  return false;
}

// The nested weight child of the indicator: a scalar unobserved dbeta whose
// parameters are scalars depending on the indicator only (the shared weight
// coordinate of the nested layout, BayesTools_nested_a).
inline bool nested_weight_node(StochasticNode const *node, StochasticNode const *indicator)
{
  if (node->distribution()->name() != "dbeta" || node->length() != 1 || !unobserved(node) ||
      node->lowerBound() != nullptr || node->upperBound() != nullptr || node->parents().size() != 2) return false;
  bool dependent = false;
  for (Node const *parent : node->parents()) {
    if (parent->length() != 1) return false;
    // the stochastic inputs of the parameter: the indicator only
    std::vector<Node const *> pending = {parent};
    std::set<Node const *> seen;
    while (!pending.empty()) {
      Node const *current = pending.back();
      pending.pop_back();
      if (!seen.insert(current).second || current->isFixed()) continue;
      if (current->isStochastic()) {
        if (current != indicator) return false;
        dependent = true;
        continue;
      }
      pending.insert(pending.end(), current->parents().begin(), current->parents().end());
    }
  }
  return dependent;
}

// A BayesTools gate marker: an observed dbt_gate with value zero (constant
// density in its parents).
inline bool gate_marker(StochasticNode const *node)
{
  if (node->distribution()->name() != "dbt_gate" || node->length() != 1) return false;
  for (unsigned long i = 0; i < node->length(); ++i) if (!node->isObserved(i)) return false;
  return node->value(0)[0] == 0.0;
}

// The bias tag's declaration (collect()): the indicator free and unobserved,
// the codes a fixed vector of family codes (0 to 3) with at most one PET and
// one PEESE label, each declared source free and unobserved or a fixed scalar
// (absent), and the intercept slot among the groups. The reason of a refusal
// in `reason` (the location members stay owned).
inline bool read_bias_declaration(Node const *tag, std::set<StochasticNode *> const &available,
                                  std::size_t n_groups, BiasDeclared &out, std::string &reason)
{
  auto const &par = tag->parents();
  out.tag = tag;
  out.indicator = dynamic_cast<StochasticNode *>(const_cast<Node *>(par[0]));
  if (out.indicator == nullptr || !available.count(out.indicator) || !unobserved(out.indicator)) {
    reason = "indicator not free or observed";
    return false;
  }
  Node const *codes = par[1];
  if (!codes->isFixed() || codes->length() < 2) {
    reason = "codes";
    return false;
  }
  int pet = 0, peese = 0;
  out.codes.clear();
  for (unsigned long i = 0; i < codes->length(); ++i) {
    const double v = codes->value(0)[i];
    if (!(v == 0 || v == 1 || v == 2 || v == 3)) {
      reason = "codes";
      return false;
    }
    out.codes.push_back(static_cast<int>(v));
    if (v == 2) ++pet;
    if (v == 3) ++peese;
  }
  if (pet > 1 || peese > 1) {
    reason = "codes";
    return false;
  }
  for (int b = 0; b < 2; ++b) {
    Node const *source = par[2 + b];
    const bool present = (b == 0 ? pet : peese) == 1;
    auto *node = dynamic_cast<StochasticNode *>(const_cast<Node *>(source));
    if (present) {
      if (node == nullptr || !available.count(node) || !unobserved(node)) {
        reason = "slope source not free or observed";
        return false;
      }
      out.sources[b] = node;
    } else {
      if (!fixed_scalar(source)) {
        reason = "absent slope is not a constant";
        return false;
      }
      out.sources[b] = nullptr;
    }
  }
  Node const *slot = par[4];
  if (!fixed_scalar(slot)) {
    reason = "intercept slot";
    return false;
  }
  const double v = slot->value(0)[0];
  if (!(v >= 1 && v == std::floor(v) && v <= static_cast<double>(n_groups))) {
    reason = "intercept slot";
    return false;
  }
  out.intercept_slot = static_cast<unsigned>(v);
  return true;
}

// Whether a likelihood child is a selection child in step mode at the graph's
// current values of the chain (LocationBlock::is_selected).
inline bool kernel_selected(ChildSpec const &c, unsigned chain)
{
  auto const &par = c.node->parents();
  switch (c.type) {
  case ChildType::normal:
    return false;
  case ChildType::step_switch:
    return par[9]->value(chain)[0] != SELKERNEL_NORMAL && par[7]->value(chain)[0] != 0;
  case ChildType::cluster:
  case ChildType::factor:
    return par[10]->value(chain)[0] != SELKERNEL_NORMAL;
  case ChildType::mv:
    return par[9]->value(chain)[0] != SELKERNEL_NORMAL;
  }
  return true;
}

// A finite readout tag of BayesTools on the indicator (sampler.probabilities):
// its recorded probabilities are exact conditionals, which a joint label
// update does not give.
inline bool finite_readout_tag(StochasticNode *node)
{
  for (DeterministicNode const *child : *node->deterministicChildren()) {
    auto const *logical = dynamic_cast<LogicalNode const *>(child);
    if (logical == nullptr || child->parents().empty() || child->parents()[0] != node) continue;
    std::vector<std::string> names(child->parents().size(), "w");
    names[0] = "x";
    const std::string text = logical->deparse(names);
    if (text.find("bt_finite_tag(") != std::string::npos ||
        text.find("bt_finite_witness_tag(") != std::string::npos) return true;
  }
  return false;
}

// The bias members from their declaration: the indicator's law and label
// families, each slope source's law and chart, the intercept group.
inline bool init_bias(BiasDeclared const &declared, Shared const &S, BiasSpec &B)
{
  B.indicator = declared.indicator;
  const CategoricalWeights weights = categorical_weights(B.indicator);
  if (!weights || !isSupportFixed(B.indicator) || weights.size != declared.codes.size()) {
    return fail("bias: indicator law");
  }
  if (finite_readout_tag(B.indicator)) return fail("bias: exact probability readout of the indicator");
  B.levels = static_cast<int>(weights.size);
  B.family = declared.codes;
  B.log_prior.assign(weights.size, 0.0);
  const double log_total = std::log(weights.total);
  for (unsigned long i = 0; i < weights.size; ++i) B.log_prior[i] = std::log(weights.values[i]) - log_total;
  bool selection = false;
  for (int code : B.family) selection = selection || code == 1;
  if (!selection) return fail("bias: no selection label");
  B.slope_of.assign(weights.size, -1);
  for (int b = 0; b < 2; ++b) {
    SlopeSpec &spec = B.slope[b];
    spec.source = declared.sources[b];
    if (spec.source == nullptr) continue;
    const int code = b == 0 ? 2 : 3;
    for (std::size_t l = 0; l < B.family.size(); ++l) if (B.family[l] == code) spec.label = static_cast<int>(l) + 1;
    if (spec.label == 0 || !slope_source_law(spec.source)) return fail("bias: slope source law or support");
    spec.dist = dynamic_cast<ScalarDist const *>(spec.source->distribution());
    for (DeterministicNode *child : *spec.source->deterministicChildren()) {
      if (slope_chart(child, spec.source, B.indicator) == 0) continue;
      if (spec.physical != nullptr) return fail("bias: two slope charts");
      spec.physical = child;
    }
    if (spec.physical == nullptr || slope_chart(spec.physical, spec.source, B.indicator) != spec.label) {
      return fail("bias: slope chart");
    }
    B.slope_of[static_cast<std::size_t>(spec.label - 1)] = b;
  }
  for (std::size_t l = 0; l < B.family.size(); ++l) {
    if ((B.family[l] == 2 || B.family[l] == 3) && B.slope_of[l] < 0) return fail("bias: slope label without a source");
  }
  B.intercept = declared.intercept_slot - 1;
  if (B.intercept >= S.groups.size() || S.groups[B.intercept].d != 1) return fail("bias: intercept group");
  return true;
}

// The stochastic children of an owned bias label that are not likelihood
// children: gate markers (constant densities) and at most one nested weight
// node. Every other stochastic child must pass the block's child checks.
inline bool classify_bias_children(GraphView const &view, BiasSpec &B)
{
  B.skip.clear();
  B.nested = nullptr;
  for (StochasticNode *child : view.stochasticChildren()) {
    if (gate_marker(child)) {
      B.skip.insert(child);
    } else if (nested_weight_node(child, B.indicator)) {
      if (B.nested != nullptr) return fail("bias: two nested weight nodes");
      B.nested = child;
      B.skip.insert(child);
    }
  }
  return true;
}

// The label families and the slope columns, probed on chain 0 (restored):
// at each label every child's kernel mode (NORMAL exactly on the labels
// coded none, PET and PEESE; some child in step mode on a selection label)
// and the nested node's density (finite); per slope its column (the change of
// the means per unit source at its label, affine, no effect at another
// label, additive with the location coordinates and latents).
inline bool probe_bias(GraphView const &view, Shared &S, BiasSpec &B)
{
  Vec saved(view.length());
  view.getValue(saved, 0);
  bool ok = true;
  Vec v;
  auto means = [&](Vec &out) {
    view.setValue(v, 0);
    out.assign(S.rows, 0.0);
    for (ChildSpec const &c : S.children) {
      double const *m = c.node->parents()[0]->value(0);
      for (unsigned i = 0; i < c.n; ++i) out[c.row0 + i] = m[i];
    }
  };
  auto zero_location = [&]() {
    for (GroupSpec const &g : S.groups) {
      if (g.gate != nullptr) v[g.gate_pos] = g.label;
      for (unsigned i = 0; i < g.d; ++i) v[g.slab_pos + i] = 0.0;
    }
    for (unsigned pos : S.latent_pos) v[pos] = 0.0;
  };
  try {
    for (int l = 1; ok && l <= B.levels; ++l) {
      v = saved;
      v[B.pos] = l;
      view.setValue(v, 0);
      bool any = false;
      for (ChildSpec const &c : S.children) any = any || kernel_selected(c, 0);
      if (any != !B.gaussian(l)) ok = fail("bias: kernel mode does not follow the label families");
      if (ok && B.nested != nullptr && !std::isfinite(B.nested->logDensity(0, PDF_LIKELIHOOD))) {
        ok = fail("bias: nested node density");
      }
    }
    for (int b = 0; ok && b < 2; ++b) {
      SlopeSpec &spec = B.slope[b];
      if (spec.source == nullptr) continue;
      Vec mu0, mu1, mu2, latent_one;
      v = saved;
      v[B.pos] = spec.label;
      v[spec.pos] = 0.0;
      means(mu0);
      v[spec.pos] = 1.0;
      means(mu1);
      v[spec.pos] = 2.0;
      means(mu2);
      spec.column.assign(S.rows, 0.0);
      for (unsigned r = 0; ok && r < S.rows; ++r) {
        spec.column[r] = mu1[r] - mu0[r];
        const double scale = 1.0 + std::fabs(mu0[r]);
        if (!std::isfinite(spec.column[r]) ||
            !(std::fabs(mu2[r] - (mu0[r] + 2.0 * spec.column[r])) <= 1e-9 * (scale + std::fabs(mu2[r])))) {
          ok = fail("bias: slope column not affine");
        }
      }
      // no effect at another label
      int other = 0;
      for (int l = 1; l <= B.levels && other == 0; ++l) if (l != spec.label) other = l;
      v = saved;
      v[B.pos] = other;
      v[spec.pos] = 0.0;
      means(mu0);
      v[spec.pos] = 1.0;
      means(mu1);
      for (unsigned r = 0; ok && r < S.rows; ++r) if (mu1[r] != mu0[r]) ok = fail("bias: slope acts at another label");
      // additive with the location coordinates and latents: every
      // coordinate, the latents and the source at one against the latents at
      // one alone plus the probed columns
      if (ok) {
        v = saved;
        v[B.pos] = spec.label;
        zero_location();
        v[spec.pos] = 0.0;
        for (unsigned pos : S.latent_pos) v[pos] = 1.0;
        means(latent_one);
        for (GroupSpec const &g : S.groups) for (unsigned i = 0; i < g.d; ++i) v[g.slab_pos + i] = 1.0;
        v[spec.pos] = 1.0;
        means(mu1);
        for (unsigned r = 0; ok && r < S.rows; ++r) {
          long double x = latent_one[r] + static_cast<long double>(spec.column[r]);
          for (unsigned k = 0; k < S.p; ++k) x += S.X[r * S.p + k];
          const double scale = 1.0 + std::fabs(latent_one[r]);
          if (!(std::fabs(mu1[r] - static_cast<double>(x)) <= 1e-9 * (scale + std::fabs(mu1[r])))) {
            ok = fail("bias: slope column not additive");
          }
        }
      }
      spec.children.clear();
      for (unsigned c = 0; c < S.children.size(); ++c) {
        ChildSpec const &child = S.children[c];
        for (unsigned i = 0; i < child.n; ++i) {
          if (spec.column[child.row0 + i] != 0.0) {
            spec.children.push_back(c);
            break;
          }
        }
      }
    }
  } catch (...) {
    view.setValue(saved, 0);
    throw;
  }
  view.setValue(saved, 0);
  return ok;
}

// ---- the evidence record of the bias label (bias_evidence): per chain,
// over the updates after adaptation, the sweeps ending in each label, the
// proposals and acceptances of L0, L1, S and L2 (L2 also its stage-1 passes),
// and by source family (0: leaving the selection family, 1: leaving the
// Gaussian family) for L2 the count, sum and maximum of the stage-1 log ratio
// and the sum of its pass probabilities min(1, exp(stage 1)) over every
// attempt, the count, sum and maximum of the stage-1 log ratio of the swaps
// retained after the sweep's exact correction, and for L0 (from the label
// scores it computes, no further evaluation) the updates made, the sum of
// their probabilities of moving to the other family, the moves it made there
// and the maximum of log P(other family) - log P(own family) of its
// surrogate full conditional. Never destroyed, like the counters.
struct BiasEvidence {
  std::vector<double> label_sweeps;
  double proposed[4] = {0, 0, 0, 0}, stage1[4] = {0, 0, 0, 0}, accepted[4] = {0, 0, 0, 0};
  double attempts[2] = {0, 0}, attempt_sum[2] = {0, 0}, pass_sum[2] = {0, 0};
  double attempt_max[2] = {-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
  double kept[2] = {0, 0}, kept_sum[2] = {0, 0};
  double kept_max[2] = {-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};
  double label_updates[2] = {0, 0}, label_cross_sum[2] = {0, 0}, label_crossings[2] = {0, 0};
  double label_log_ratio_max[2] = {-std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()};

  void add(BiasEvidence const &other)
  {
    if (label_sweeps.size() < other.label_sweeps.size()) label_sweeps.resize(other.label_sweeps.size(), 0.0);
    for (std::size_t i = 0; i < other.label_sweeps.size(); ++i) label_sweeps[i] += other.label_sweeps[i];
    for (int i = 0; i < 4; ++i) {
      proposed[i] += other.proposed[i];
      stage1[i] += other.stage1[i];
      accepted[i] += other.accepted[i];
    }
    for (int d = 0; d < 2; ++d) {
      attempts[d] += other.attempts[d];
      attempt_sum[d] += other.attempt_sum[d];
      pass_sum[d] += other.pass_sum[d];
      attempt_max[d] = std::max(attempt_max[d], other.attempt_max[d]);
      kept[d] += other.kept[d];
      kept_sum[d] += other.kept_sum[d];
      kept_max[d] = std::max(kept_max[d], other.kept_max[d]);
      label_updates[d] += other.label_updates[d];
      label_cross_sum[d] += other.label_cross_sum[d];
      label_crossings[d] += other.label_crossings[d];
      label_log_ratio_max[d] = std::max(label_log_ratio_max[d], other.label_log_ratio_max[d]);
    }
  }
};

struct BiasEvidenceStore {
  std::mutex mutex;
  std::map<unsigned, BiasEvidence> records;
};

inline BiasEvidenceStore &bias_evidence_store()
{
  static BiasEvidenceStore *store = new BiasEvidenceStore();
  return *store;
}

// The component index of the evidence record: L0, L1, S, L2.
inline int bias_record_index(int component)
{
  switch (component) {
  case 1: return 0;   // bias_label
  case 2: return 1;   // bias_gaussian
  case 3: return 2;   // bias_slope
  case 4: return 3;   // bias_swap
  default: return -1;
  }
}
