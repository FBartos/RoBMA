// The bias-label components of the location block sweep (wp3-spec.md): the
// per-update state of the owned bias indicator and slope sources and the
// components D (dormant slope draws) and L0 (the screened label update). A
// source fragment of LocationBlock.cc, included inside the class body of
// LocationBlock; it is not a stand-alone header.
//
// The label is evaluated by writing it to the graph without recomputing the
// children's means: the indicator's value is set directly and only its
// deterministic descendants outside the slope charts are recomputed
// (BiasSpec::label_nodes: weights, kernel modes, rules, the nested node's
// parameters). The block evaluates every child at its own working mean, so the
// graph's means, which still hold the update's entry state, are never read
// during the sweep; correct() writes the indicator back to its entry value
// before it reads the entry's exact densities and then writes the final state
// in full.

  // ---- per-update state of the bias members
  int _label = 0, _label_entry = 0, _graph_label = 0;
  double _theta[2] = {0.0, 0.0};      // the slope sources (PET, PEESE)
  bool _selection_entry = false;      // the entry label is a selection label
  std::vector<char> _entry_valued;    // _value_entry holds the entry label's term
  Vec _label_terms;                   // levels x children: the terms of the last label scores
  Vec _label_score;                   // per label: its score (log prior + nested + children)
  std::vector<char> _step_children;   // per child: a step-switch row
  Vec _shift_mean;                    // scratch: a child's mean at another label

  bool bias() const { return s().bias != nullptr; }
  // An endpoint of the sweep is in a selection label: the exact correction
  // can differ from one only then (a sweep whose entry and final labels are
  // Gaussian has equal inner and exact targets at both ends).
  bool selection_endpoint() const { return _selection_entry || _selection; }
  BiasSpec const &bspec() const { return *s().bias; }

  // The physical slope b at label l with the current sources.
  double physical_at(int l, int b) const
  {
    SlopeSpec const &spec = bspec().slope[b];
    return spec.source != nullptr && spec.label == l ? _theta[b] : 0.0;
  }

  // The step-switch row's validity at the graph's current label: the
  // distribution's parameter check, with its label-free parts taken from the
  // update's entry (StepChild::free_ok, step_ok) and the kernel mode and
  // weights read now.
  void refresh_step(unsigned c)
  {
    StepChild &st = _step[c];
    auto const &par = s().children[c].node->parents();
    st.mode = static_cast<int>(par[9]->value(_chain)[0]);
    st.omega = par[4]->value(_chain);
    bool ok = st.free_ok && (st.mode == SELKERNEL_NORMAL || st.mode == SELKERNEL_STEP);
    if (ok && st.mode == SELKERNEL_STEP) {
      ok = st.step_ok;
      for (int b = 0; ok && b < st.bins; ++b) ok = std::isfinite(st.omega[b]) && st.omega[b] >= 0;
    }
    st.valid = ok;
  }

  // Writes label l to the indicator and its label nodes (not the means).
  void set_graph_label(int l)
  {
    if (l == _graph_label) return;
    BiasSpec const &B = bspec();
    const double value = l;
    B.indicator->setValue(&value, 1, _chain);
    for (DeterministicNode *node : B.label_nodes) node->randomSample(nullptr, _chain);
    _graph_label = l;
    for (unsigned c = 0; c < s().children.size(); ++c) if (_step_children[c]) refresh_step(c);
  }

  // After prepare(): the bias members' entry values.
  void bias_prepare()
  {
    if (!bias()) return;
    BiasSpec const &B = bspec();
    _label = _label_entry = _graph_label = static_cast<int>(_entry[B.pos]);
    for (int b = 0; b < 2; ++b) _theta[b] = B.slope[b].source != nullptr ? _entry[B.slope[b].pos] : 0.0;
  }

  // The nested node's term at the graph's label (PDF_LIKELIHOOD, as the
  // indicator's full conditional counts it).
  double nested_term() const
  {
    StochasticNode const *nested = bspec().nested;
    return nested == nullptr ? 0.0 : nested->logDensity(_chain, PDF_LIKELIHOOD);
  }

  // The inner target's term of child c at the graph's label with the slopes
  // of label l (its mean shifted from the working mean by the change of the
  // physical slopes): the surrogate or exact selected-normal term when the
  // child is a selection child at l, its exact (normal) density otherwise.
  // *selected receives the child's status at l.
  double term_at_label(unsigned c, int l, bool *selected)
  {
    ChildSpec const &child = s().children[c];
    BiasSpec const &B = bspec();
    _shift_mean.assign(_mean.begin() + child.row0, _mean.begin() + child.row0 + child.n);
    for (int b = 0; b < 2; ++b) {
      SlopeSpec const &spec = B.slope[b];
      if (spec.source == nullptr) continue;
      const double change = physical_at(l, b) - physical_at(_label, b);
      if (change == 0.0) continue;
      for (unsigned i = 0; i < child.n; ++i) _shift_mean[i] += change * spec.column[child.row0 + i];
    }
    *selected = is_selected(child);
    if (*selected) return selected_value(c, _shift_mean.data());
    const double v = exact_with_mean(c, _shift_mean.data());
    if (std::isnan(v) || v == std::numeric_limits<double>::infinity()) {
      throwNodeError(child.node, "Location block sampler: the log density of this likelihood node is not defined at a label");
    }
    return v;
  }

  // The scores of every label with prior mass at the current working state:
  // log prior + the nested node's term + every child's term (L0, sweep
  // condition 5: the block's buffers, the graph only for the label's own
  // parameters). The terms are kept in _label_terms; the graph is left at the
  // working label.
  void score_labels()
  {
    BiasSpec const &B = bspec();
    const std::size_t n = s().children.size();
    _label_terms.assign(static_cast<std::size_t>(B.levels) * n, 0.0);
    _label_score.assign(static_cast<std::size_t>(B.levels), negative_infinity);
    for (int l = 1; l <= B.levels; ++l) {
      const double prior = B.log_prior[static_cast<std::size_t>(l - 1)];
      if (!(prior > negative_infinity)) continue;
      set_graph_label(l);
      long double score = prior + nested_term();
      double *terms = &_label_terms[static_cast<std::size_t>(l - 1) * n];
      for (unsigned c = 0; c < n; ++c) {
        bool selected;
        terms[c] = term_at_label(c, l, &selected);
        score += terms[c];
      }
      _label_score[static_cast<std::size_t>(l - 1)] = static_cast<double>(score);
    }
    set_graph_label(_label);
  }

  // Moves the working state to label l (the sources unchanged): the physical
  // slopes, means and residuals, the selection status of every child, the
  // inner-target terms (terms: the children's terms at l, or null to
  // evaluate them lazily) and the Gaussian working model's linear term.
  void apply_label(int l, double const *terms)
  {
    BiasSpec const &B = bspec();
    Shift sh;
    for (int b = 0; b < 2; ++b) {
      if (B.slope[b].source == nullptr) continue;
      const double change = physical_at(l, b) - physical_at(_label, b);
      if (change != 0.0) {
        sh.slopes.push_back(b);
        sh.dslope.push_back(change);
      }
    }
    set_graph_label(l);
    _label = l;
    if (!sh.slopes.empty()) {
      prepare_shift(sh, {});
      move_means(sh, 1.0);
    }
    _selection = false;
    for (unsigned c = 0; c < s().children.size(); ++c) {
      _selected[c] = is_selected(s().children[c]);
      if (_selected[c]) _selection = true;
    }
    for (unsigned c = 0; c < s().children.size(); ++c) {
      _fit_ok[c] = 0;
      if (terms != nullptr && !gaussian_child(c)) {
        _value[c] = terms[c];
        _value_ok[c] = 1;
      } else {
        _value_ok[c] = 0;
      }
    }
    _line.valid = false;
    if (!sh.slopes.empty()) rebuild_h();
  }

  // The means, residuals and weighted residuals moved by t times a prepared
  // shift (no member coordinates; the children become dirty).
  void move_means(Shift const &sh, double t)
  {
    Shared const &S = s();
    for (std::size_t q = 0; q < sh.children.size(); ++q) {
      const unsigned c = sh.children[q];
      ChildSpec const &child = S.children[c];
      double const *delta = &sh.delta[sh.offset[q]];
      double const *wdelta = &sh.wdelta[sh.offset[q]];
      for (unsigned i = 0; i < child.n; ++i) {
        const unsigned row = child.row0 + i;
        _mean[row] += t * delta[i];
        _r[row] -= t * delta[i];
        _wr[row] -= t * wdelta[i];
      }
      _fit_ok[c] = 0;
      _dirty[c] = 1;
    }
  }

  // U' S^-1 V for row vectors (S = Sigma + sum_s sd_s^2 h_s h_s' over the
  // latents), per latent in the within/between form of prepare().
  void cross_S(std::vector<double const *> const &U, std::vector<double const *> const &V, Vec &out)
  {
    Shared const &S = s();
    const std::size_t nu = U.size(), nv = V.size();
    std::vector<long double> M(nu * nv, 0.0L);
    Vec wv, wd;
    auto products = [&](ChildSpec const &child) {
      const unsigned n = child.n;
      double const *W = &_W[child.woff];
      wv.assign(static_cast<std::size_t>(n) * nv, 0.0);
      wd.assign(n, 0.0);
      for (unsigned i = 0; i < n; ++i) {
        for (unsigned j = 0; j < n; ++j) {
          const double w = W[i * n + j];
          if (w == 0.0) continue;
          const unsigned row = child.row0 + j;
          for (std::size_t k = 0; k < nv; ++k) wv[i * nv + k] += w * V[k][row];
          wd[i] += w * _d[row];
        }
      }
    };
    for (ChildSpec const &child : S.children) {
      if (child.latent >= 0) continue;
      products(child);
      for (unsigned i = 0; i < child.n; ++i) {
        const unsigned row = child.row0 + i;
        for (std::size_t a = 0; a < nu; ++a) {
          const double u = U[a][row];
          if (u == 0.0) continue;
          for (std::size_t k = 0; k < nv; ++k) M[a * nv + k] += static_cast<long double>(u) * wv[i * nv + k];
        }
      }
    }
    for (unsigned lat = 0; lat < S.latents.size(); ++lat) {
      std::vector<unsigned> const &kids = S.latent_children[lat];
      const double a = _a[lat];
      const bool between = a > 0;
      std::vector<long double> ubar(nu, 0.0L), vbar(nv, 0.0L);
      if (between) {
        for (unsigned c : kids) {
          ChildSpec const &child = S.children[c];
          products(child);
          for (unsigned i = 0; i < child.n; ++i) {
            const unsigned row = child.row0 + i;
            for (std::size_t q = 0; q < nu; ++q) ubar[q] += static_cast<long double>(U[q][row]) * wd[i];
            for (std::size_t k = 0; k < nv; ++k) vbar[k] += static_cast<long double>(V[k][row]) * wd[i];
          }
        }
        for (auto &x : ubar) x /= a;
        for (auto &x : vbar) x /= a;
      }
      for (unsigned c : kids) {
        ChildSpec const &child = S.children[c];
        products(child);
        for (unsigned i = 0; i < child.n; ++i) {
          const unsigned row = child.row0 + i;
          const double di = _d[row];
          for (std::size_t q = 0; q < nu; ++q) {
            const long double uc = U[q][row] - (between ? di * ubar[q] : 0.0L);
            if (uc == 0.0L) continue;
            for (std::size_t k = 0; k < nv; ++k) {
              M[q * nv + k] += uc * (wv[i * nv + k] - (between ? wd[i] * vbar[k] : 0.0L));
            }
          }
        }
      }
      if (between) {
        const long double scale = a / _lambda[lat];
        for (std::size_t q = 0; q < nu; ++q) for (std::size_t k = 0; k < nv; ++k) M[q * nv + k] += scale * ubar[q] * vbar[k];
      }
    }
    out.assign(nu * nv, 0.0);
    for (std::size_t i = 0; i < nu * nv; ++i) out[i] = static_cast<double>(M[i]);
  }

  // The offset of the working model, r0 = y - offset = r + X b + d z (rows).
  void working_r0(Vec &r0) const
  {
    Shared const &S = s();
    const unsigned p = S.p;
    r0.assign(S.rows, 0.0);
    for (ChildSpec const &child : S.children) {
      for (unsigned i = 0; i < child.n; ++i) {
        const unsigned row = child.row0 + i;
        long double v = _r[row];
        for (unsigned k = 0; k < p; ++k) v += static_cast<long double>(S.X[row * p + k]) * _b[k];
        if (child.latent >= 0) v += static_cast<long double>(_d[row]) * _z[child.latent];
        r0[row] = static_cast<double>(v);
      }
    }
  }

  // h = B m0 + X' S^-1 r0 from the current offset (after a change of the
  // physical slopes, which are part of the offset).
  void rebuild_h()
  {
    Shared const &S = s();
    const unsigned p = S.p;
    Vec r0, cross;
    working_r0(r0);
    std::vector<Vec> columns(p, Vec(S.rows, 0.0));
    std::vector<double const *> U(p);
    for (unsigned k = 0; k < p; ++k) {
      for (unsigned row = 0; row < S.rows; ++row) columns[k][row] = S.X[row * p + k];
      U[k] = columns[k].data();
    }
    cross_S(U, {r0.data()}, cross);
    for (unsigned k = 0; k < p; ++k) _h[k] = cross[k];
    for (GroupSpec const &g : S.groups) {
      for (unsigned i = 0; i < g.d; ++i) {
        long double bm = 0;
        for (unsigned j = 0; j < g.d; ++j) bm += static_cast<long double>(g.B[i * g.d + j]) * g.m0[j];
        _h[g.coord0 + i] += static_cast<double>(bm);
      }
    }
  }

  // ---- D: the dormant slopes, exact draws from their priors (a Gibbs step:
  // a dormant slope enters nothing but its prior)
  void dormant_draws(RNG *rng)
  {
    BiasSpec const &B = bspec();
    for (int b = 0; b < 2; ++b) {
      SlopeSpec const &spec = B.slope[b];
      if (spec.source == nullptr || spec.label == _label) continue;
      count(counters().bias_proposed[bias_dormant]);
      const SlopePrior prior = slope_prior(spec, _chain);
      const double draw = prior.dist->randomSample(prior.par, prior.lower, prior.upper, rng);
      const double log_new = slope_log_prior(prior, draw);
      const double log_old = slope_log_prior(prior, _theta[b]);
      if (!std::isfinite(draw) || !std::isfinite(log_new) || !std::isfinite(log_old)) continue;
      count(counters().bias_accepted[bias_dormant]);
      _repr_delta += static_cast<long double>(log_new) - log_old;
      _theta[b] = draw;
    }
  }

  // ---- L0: the screened label update. Liu's Metropolized Gibbs step on the
  // label's surrogate full conditional at the current working state (every
  // label with prior mass: the surrogate terms of the coarse-capable children,
  // the exact terms of the others, the nested node's term), reversible for
  // the inner target; the sweep's exact correction covers it.
  void label_update(RNG *rng)
  {
    BiasSpec const &B = bspec();
    score_labels();
    const std::size_t n = s().children.size();
    const int k = _label;
    // the entry label's terms are the children's entry values
    {
      double const *terms = &_label_terms[static_cast<std::size_t>(k - 1) * n];
      for (unsigned c = 0; c < n; ++c) {
        if (gaussian_child(c) || _value_ok[c]) continue;
        _value[c] = terms[c];
        _value_ok[c] = 1;
        if (_label == _label_entry && !_dirty[c]) {
          _value_entry[c] = terms[c];
          _entry_valued[c] = 1;
        }
      }
    }
    double maximum = negative_infinity;
    for (double v : _label_score) maximum = std::max(maximum, v);
    if (!std::isfinite(maximum) || !std::isfinite(_label_score[static_cast<std::size_t>(k - 1)])) return;
    std::vector<long double> mass(static_cast<std::size_t>(B.levels), 0.0L);
    long double total = 0;
    for (int l = 1; l <= B.levels; ++l) {
      const double v = _label_score[static_cast<std::size_t>(l - 1)];
      mass[static_cast<std::size_t>(l - 1)] = v > negative_infinity ? std::exp(static_cast<long double>(v - maximum)) : 0.0L;
      total += mass[static_cast<std::size_t>(l - 1)];
    }
    const long double others = total - mass[static_cast<std::size_t>(k - 1)];
    if (!(others > 0)) return;
    count(counters().bias_proposed[bias_label]);
    long double draw = static_cast<long double>(rng->uniform()) * others;
    int proposed = 0;
    for (int l = 1; l <= B.levels; ++l) {
      if (l == k || mass[static_cast<std::size_t>(l - 1)] == 0.0L) continue;
      proposed = l;
      if (draw < mass[static_cast<std::size_t>(l - 1)]) break;
      draw -= mass[static_cast<std::size_t>(l - 1)];
    }
    if (proposed == 0) return;
    const long double others_at = total - mass[static_cast<std::size_t>(proposed - 1)];
    const long double ratio = others / others_at;
    if (!(ratio >= 1.0L || static_cast<long double>(rng->uniform()) < ratio)) return;
    count(counters().bias_accepted[bias_label]);
    _repr_delta += static_cast<long double>(_label_score[static_cast<std::size_t>(proposed - 1)]) -
      _label_score[static_cast<std::size_t>(k - 1)];
    apply_label(proposed, &_label_terms[static_cast<std::size_t>(proposed - 1) * n]);
  }

  // ---- the bias slots of the sweep

  // Decided once per sweep before any component runs, independently of the
  // state (sweep conditions; the schedule of L2).
  bool _swap_this_sweep = false;

  void bias_sweep_start(RNG *rng)
  {
    _swap_this_sweep = false;
    (void) rng;
  }

  Category bias_category(unsigned k) const
  {
    switch (k) {
    case bias_dormant: return cat_dormant;
    case bias_label: return cat_label;
    case bias_gaussian: return cat_gaussian;
    case bias_slope: return cat_slope;
    default: return cat_swap;
    }
  }

  void bias_component(unsigned k, RNG *rng)
  {
    switch (k) {
    case bias_dormant:
      dormant_draws(rng);
      break;
    case bias_label:
      label_update(rng);
      break;
    default:
      break;
    }
  }
