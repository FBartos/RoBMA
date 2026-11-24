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

  // probe record types of the bias components (LocationBlockReport.h)
  static constexpr int probe_family = 7, probe_slope = 8, probe_swap = 9;

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

  // ---- L1 and S: the Gaussian-family moves (wp3-spec.md, L1 and S)
  //
  // In a label none, PET or PEESE every child is Gaussian. With the other
  // location coefficients fixed, r0 = y - X_rest b_rest (no intercept, slope or
  // latent), S = Sigma + sum_s d_s d_s' (the latents integrated) and the
  // intercept's slab N(m0, 1/B0), the likelihood of cell c = (label, gate) in
  // the active slope theta is k_c + r_c theta - q_c theta^2 / 2 (the intercept
  // integrated when its gate is on).

  struct Collapse {
    bool valid = false;
    int intercept = 0;                // the intercept group
    bool gated = false;
    double B0 = 0, m0 = 0, log_on = 0, log_off = 0;
    double xx = 0, xr = 0, rr = 0;    // x0'S^-1 x0, x0'S^-1 r0, r0'S^-1 r0
    double xs[2] = {0, 0}, ss[2] = {0, 0}, sr[2] = {0, 0};
    Vec r0;                           // rows
    // per latent: d'W r0, d'W x0, d'W s_b (the latents' conditional means)
    Vec dr, dx, ds[2];
  };

  struct Cell {
    int label = 0, slope = -1;        // the slope active in the label (-1 none)
    bool gate = true;
    double q = 0, r = 0, k = 0;       // the likelihood in theta
    double P0 = 0, a = 0, h = 0;      // the intercept's law given theta (gate on)
    bool prior_cell = false;
    double anchor = 0;                // the truncated-normal mean (regular cells)
    double log_weight = 0;            // log proposal weight (unnormalized)
  };

  // The collapsed Gaussian quantities at the current state (state the moves
  // do not change: the other coefficients, the precisions, the latent column).
  bool collapse(Collapse &C)
  {
    Shared const &S = s();
    BiasSpec const &B = bspec();
    C.valid = false;
    C.intercept = static_cast<int>(B.intercept);
    GroupSpec const &g0 = S.groups[B.intercept];
    C.gated = g0.gate != nullptr;
    C.B0 = g0.B[0];
    C.m0 = g0.m0[0];
    C.log_on = C.gated ? g0.log_on : 0.0;
    C.log_off = C.gated ? g0.log_off : negative_infinity;
    const unsigned p = S.p, k0 = g0.coord0;
    // r0 = r + x0 t + sum_b s_b physical_b + d z
    C.r0.assign(S.rows, 0.0);
    Vec x0(S.rows, 0.0);
    for (ChildSpec const &child : S.children) {
      for (unsigned i = 0; i < child.n; ++i) {
        const unsigned row = child.row0 + i;
        x0[row] = S.X[row * p + k0];
        long double v = _r[row];
        v += static_cast<long double>(x0[row]) * _b[k0];
        for (int b = 0; b < 2; ++b) {
          if (B.slope[b].source != nullptr) v += static_cast<long double>(B.slope[b].column[row]) * physical_at(_label, b);
        }
        if (child.latent >= 0) v += static_cast<long double>(_d[row]) * _z[child.latent];
        C.r0[row] = static_cast<double>(v);
      }
    }
    Vec zero(S.rows, 0.0);
    std::vector<double const *> U = {x0.data(), C.r0.data(),
      B.slope[0].source != nullptr ? B.slope[0].column.data() : zero.data(),
      B.slope[1].source != nullptr ? B.slope[1].column.data() : zero.data()};
    Vec G;
    cross_S(U, U, G);
    auto at = [&](int i, int j) { return 0.5 * (G[i * 4 + j] + G[j * 4 + i]); };
    C.xx = at(0, 0);
    C.xr = at(0, 1);
    C.rr = at(1, 1);
    for (int b = 0; b < 2; ++b) {
      C.xs[b] = at(0, 2 + b);
      C.ss[b] = at(2 + b, 2 + b);
      C.sr[b] = at(1, 2 + b);
    }
    // the latents' conditional means: d'W v over each latent's children
    const unsigned L = static_cast<unsigned>(S.latents.size());
    C.dr.assign(L, 0.0);
    C.dx.assign(L, 0.0);
    C.ds[0].assign(L, 0.0);
    C.ds[1].assign(L, 0.0);
    for (unsigned lat = 0; lat < L; ++lat) {
      long double vr = 0, vx = 0, v0 = 0, v1 = 0;
      for (unsigned c : S.latent_children[lat]) {
        ChildSpec const &child = S.children[c];
        double const *W = &_W[child.woff];
        for (unsigned i = 0; i < child.n; ++i) {
          long double wd = 0;
          for (unsigned j = 0; j < child.n; ++j) wd += static_cast<long double>(W[i * child.n + j]) * _d[child.row0 + j];
          const unsigned row = child.row0 + i;
          vr += wd * C.r0[row];
          vx += wd * x0[row];
          if (B.slope[0].source != nullptr) v0 += wd * B.slope[0].column[row];
          if (B.slope[1].source != nullptr) v1 += wd * B.slope[1].column[row];
        }
      }
      C.dr[lat] = static_cast<double>(vr);
      C.dx[lat] = static_cast<double>(vx);
      C.ds[0][lat] = static_cast<double>(v0);
      C.ds[1][lat] = static_cast<double>(v1);
    }
    C.valid = std::isfinite(C.xx) && std::isfinite(C.rr) && C.B0 > 0;
    return C.valid;
  }

  // The cell's likelihood in theta (spec L1) and its proposal weight: label
  // prior and nested term (lnest), gate prior, k_c and I_c (the slope's
  // normalized prior at the anchor times the truncated Gaussian integral; 1
  // for none and for a prior cell, whose slope is proposed from its prior).
  Cell make_cell(Collapse const &C, int label, bool gate, double lnest)
  {
    BiasSpec const &B = bspec();
    Cell cell;
    cell.label = label;
    cell.gate = gate;
    cell.slope = B.slope_of[static_cast<std::size_t>(label - 1)];
    const int b = cell.slope;
    const double ss = b >= 0 ? C.ss[b] : 0.0, sr = b >= 0 ? C.sr[b] : 0.0, xs = b >= 0 ? C.xs[b] : 0.0;
    if (gate) {
      cell.P0 = C.B0 + C.xx;
      cell.a = xs;
      cell.h = C.B0 * C.m0 + C.xr;
      cell.q = ss - cell.a * cell.a / cell.P0;
      cell.r = sr - cell.a * cell.h / cell.P0;
      cell.k = -0.5 * C.rr + cell.h * cell.h / (2.0 * cell.P0) - 0.5 * C.B0 * C.m0 * C.m0 + 0.5 * std::log(C.B0 / cell.P0);
    } else {
      cell.q = ss;
      cell.r = sr;
      cell.k = -0.5 * C.rr;
    }
    cell.log_weight = B.log_prior[static_cast<std::size_t>(label - 1)] + lnest + (gate ? C.log_on : C.log_off) + cell.k;
    if (b < 0) return cell;
    const SlopePrior prior = slope_prior(B.slope[b], _chain);
    double log_anchor = negative_infinity;
    if (cell.q > 0 && cell.q >= 1e-12 * ss) {
      const double mean = cell.r / cell.q, sd = 1.0 / std::sqrt(cell.q);
      // the mean of N(mean, sd^2) truncated to [0, Inf): mean + sd lambda(-mean / sd)
      const double alpha = -mean / sd;
      const double mills = std::exp(-0.5 * alpha * alpha - 0.91893853320467274178 - log_phi_upper(-alpha));
      cell.anchor = mean + sd * mills;
      if (cell.anchor > 0 && std::isfinite(cell.anchor)) log_anchor = slope_log_prior(prior, cell.anchor);
    }
    if (!std::isfinite(log_anchor)) {
      cell.prior_cell = true;
      return cell;
    }
    const double sd = 1.0 / std::sqrt(cell.q);
    cell.log_weight += log_anchor + 0.5 * std::log(2.0 * 3.14159265358979323846 / cell.q) +
      cell.r * cell.r / (2.0 * cell.q) + log_phi_upper(cell.r * sd);
    return cell;
  }

  // The cells of the Gaussian family with prior mass.
  void family_cells(Collapse const &C, std::vector<Cell> &cells)
  {
    BiasSpec const &B = bspec();
    cells.clear();
    for (int l = 1; l <= B.levels; ++l) {
      if (!B.gaussian(l) || !(B.log_prior[static_cast<std::size_t>(l - 1)] > negative_infinity)) continue;
      set_graph_label(l);
      const double lnest = nested_term();
      for (int gate = C.gated ? 0 : 1; gate <= 1; ++gate) cells.push_back(make_cell(C, l, gate == 1, lnest));
    }
    set_graph_label(_label);
  }

  // log density of the active slope under its cell's proposal: the truncated
  // normal of a regular cell, the prior of a prior cell.
  double cell_slope_density(Cell const &cell, double theta)
  {
    if (cell.slope < 0) return 0.0;
    if (cell.prior_cell) return slope_log_prior(slope_prior(bspec().slope[cell.slope], _chain), theta);
    return log_truncated_normal_positive(theta, cell.r / cell.q, 1.0 / std::sqrt(cell.q));
  }

  // The latents' exact Gaussian conditional means given the intercept's
  // physical value t and the physical slope theta of slope b (-1 none).
  void latent_means(Collapse const &C, double t, int b, double theta, Vec &m) const
  {
    const std::size_t L = C.dr.size();
    m.assign(L, 0.0);
    for (std::size_t s = 0; s < L; ++s) {
      long double v = C.dr[s] - static_cast<long double>(C.dx[s]) * t;
      if (b >= 0) v -= static_cast<long double>(C.ds[b][s]) * theta;
      m[s] = static_cast<double>(v / _lambda[s]);
    }
  }

  double latents_log_density(Vec const &z, Vec const &m) const
  {
    long double out = 0;
    for (std::size_t s = 0; s < z.size(); ++s) {
      const double e = (z[s] - m[s]) * std::sqrt(_lambda[s]);
      out += -0.5L * e * e - 0.91893853320467274178L + 0.5L * std::log(static_cast<long double>(_lambda[s]));
    }
    return static_cast<double>(out);
  }

  // A state of the L1 block: label, intercept gate and slab, both slopes,
  // the latents.
  struct FamilyState {
    int label = 0;
    bool gate = true;
    double beta0 = 0;
    double theta[2] = {0, 0};
    Vec z;
  };

  FamilyState family_state() const
  {
    FamilyState x;
    BiasSpec const &B = bspec();
    GroupSpec const &g0 = s().groups[B.intercept];
    x.label = _label;
    x.gate = active(B.intercept);
    x.beta0 = _beta[g0.coord0];
    x.theta[0] = _theta[0];
    x.theta[1] = _theta[1];
    x.z = _z;
    return x;
  }

  // log q(x) of the L1 independence proposal (log_total: the log of the sum
  // of the cell weights).
  double family_log_q(Collapse const &C, std::vector<Cell> const &cells, double log_total, FamilyState const &x)
  {
    BiasSpec const &B = bspec();
    GroupSpec const &g0 = s().groups[B.intercept];
    Cell const *cell = nullptr;
    for (Cell const &c : cells) if (c.label == x.label && c.gate == x.gate) cell = &c;
    if (cell == nullptr) return negative_infinity;
    long double out = cell->log_weight - log_total;
    const int b = cell->slope;
    const double theta = b >= 0 ? x.theta[b] : 0.0;
    out += cell_slope_density(*cell, theta);
    double t = 0.0;
    if (x.gate) {
      const double mean = (cell->h - cell->a * theta) / cell->P0;
      const double e = (x.beta0 - mean) * std::sqrt(cell->P0);
      out += -0.5 * e * e - 0.91893853320467274178 + 0.5 * std::log(cell->P0);
      t = x.beta0;
    } else {
      out += log_slab(g0, &x.beta0);
    }
    Vec m;
    latent_means(C, t, b, theta, m);
    out += latents_log_density(x.z, m);
    for (int d = 0; d < 2; ++d) {
      if (B.slope[d].source == nullptr || d == b) continue;
      out += slope_log_prior(slope_prior(B.slope[d], _chain), x.theta[d]);
    }
    return static_cast<double>(out);
  }

  // The change of the inner target from the current state to x within the
  // L1 block (every child Gaussian at both ends): label prior and nested
  // term, gate prior, intercept slab, both slope priors, latent priors and the
  // children; sh receives the prepared shift of the means.
  long double family_delta(FamilyState const &x, Shift &sh, bool *ok)
  {
    Shared const &S = s();
    BiasSpec const &B = bspec();
    GroupSpec const &g0 = S.groups[B.intercept];
    *ok = true;
    sh.clear();
    const double t_old = _b[g0.coord0], t_new = x.gate ? x.beta0 : 0.0;
    sh.coords.push_back(g0.coord0);
    sh.dcoef.push_back(t_new - t_old);
    for (unsigned lat = 0; lat < S.latents.size(); ++lat) {
      sh.lats.push_back(lat);
      sh.dlat.push_back(x.z[lat] - _z[lat]);
    }
    for (int b = 0; b < 2; ++b) {
      if (B.slope[b].source == nullptr) continue;
      const double now = physical_at(_label, b);
      const double next = B.slope[b].label == x.label ? x.theta[b] : 0.0;
      sh.slopes.push_back(b);
      sh.dslope.push_back(next - now);
    }
    prepare_shift(sh, {});
    long double out = children_delta(sh, 1.0) + prior_delta(sh, 1.0);
    out += B.log_prior[static_cast<std::size_t>(x.label - 1)] - B.log_prior[static_cast<std::size_t>(_label - 1)];
    if (x.label != _label && B.nested != nullptr) {
      set_graph_label(x.label);
      out += nested_term();
      set_graph_label(_label);
      out -= nested_term();
    }
    if (g0.gate != nullptr && x.gate != active(B.intercept)) out += x.gate ? g0.log_on - g0.log_off : g0.log_off - g0.log_on;
    out += log_slab(g0, &x.beta0) - log_slab(g0, &_beta[g0.coord0]);
    for (int b = 0; b < 2; ++b) {
      if (B.slope[b].source == nullptr || x.theta[b] == _theta[b]) continue;
      const SlopePrior prior = slope_prior(B.slope[b], _chain);
      out += static_cast<long double>(slope_log_prior(prior, x.theta[b])) - slope_log_prior(prior, _theta[b]);
    }
    if (!std::isfinite(static_cast<double>(out))) *ok = false;
    return out;
  }

  // Moves the working state to x (the prepared shift of family_delta).
  void family_apply(FamilyState const &x, Shift &sh, long double delta)
  {
    Shared const &S = s();
    BiasSpec const &B = bspec();
    GroupSpec const &g0 = S.groups[B.intercept];
    apply_shift(sh, 1.0);
    _g[B.intercept] = x.gate ? 1 : 0;
    _beta[g0.coord0] = x.beta0;
    _b[g0.coord0] = x.gate ? x.beta0 : 0.0;
    _theta[0] = x.theta[0];
    _theta[1] = x.theta[1];
    if (x.label != _label) {
      set_graph_label(x.label);
      _label = x.label;
      _selection = false;
      for (unsigned c = 0; c < S.children.size(); ++c) {
        _selected[c] = is_selected(S.children[c]);
        if (_selected[c]) _selection = true;
        _value_ok[c] = 0;
      }
    }
    _line.valid = false;
    rebuild_h();
    _repr_delta += delta;
  }

  // The analytic L1 log ratio F_c'(theta') / F_c(theta) (spec L1, for the
  // frozen-state test): p_b(theta) / p_b(anchor) for a regular cell,
  // exp(r theta - q theta^2 / 2) for a prior cell, 1 for none.
  double family_log_F(Cell const &cell, double theta)
  {
    if (cell.slope < 0) return 0.0;
    if (cell.prior_cell) return cell.r * theta - 0.5 * cell.q * theta * theta;
    const SlopePrior prior = slope_prior(bspec().slope[cell.slope], _chain);
    return slope_log_prior(prior, theta) - slope_log_prior(prior, cell.anchor);
  }

  // ---- L1: the Gaussian-family independence move (label in none, PET,
  // PEESE x intercept gate, the slopes, intercept and latents redrawn), a
  // Metropolis-Hastings step for the inner target computed from its values
  // and the complete proposal densities at both states.
  void family_move(RNG *rng, bool dry = false)
  {
    BiasSpec const &B = bspec();
    if (_selection || !B.gaussian(_label)) return;
    Collapse C;
    if (!collapse(C)) return;
    std::vector<Cell> cells;
    family_cells(C, cells);
    if (cells.empty()) return;
    double maximum = negative_infinity;
    for (Cell const &c : cells) maximum = std::max(maximum, c.log_weight);
    if (!std::isfinite(maximum)) return;
    long double total = 0;
    for (Cell const &c : cells) total += std::exp(static_cast<long double>(c.log_weight - maximum));
    const double log_total = maximum + std::log(static_cast<double>(total));
    if (!dry) count(counters().bias_proposed[bias_gaussian]);
    // the proposal
    long double draw = static_cast<long double>(rng->uniform()) * total;
    std::size_t pick = cells.size() - 1;
    for (std::size_t i = 0; i < cells.size(); ++i) {
      const long double w = std::exp(static_cast<long double>(cells[i].log_weight - maximum));
      if (draw < w) {
        pick = i;
        break;
      }
      draw -= w;
    }
    Cell const &cell = cells[pick];
    GroupSpec const &g0 = s().groups[B.intercept];
    FamilyState next;
    next.label = cell.label;
    next.gate = cell.gate;
    for (int d = 0; d < 2; ++d) {
      if (B.slope[d].source == nullptr) continue;
      const SlopePrior prior = slope_prior(B.slope[d], _chain);
      // the active slope of a regular cell from its truncated normal, every
      // other slope (dormant, or the active one of a prior cell) from its prior
      next.theta[d] = (d == cell.slope && !cell.prior_cell) ?
        truncated_normal_positive(cell.r / cell.q, 1.0 / std::sqrt(cell.q), rng) :
        prior.dist->randomSample(prior.par, prior.lower, prior.upper, rng);
    }
    const double theta = cell.slope >= 0 ? next.theta[cell.slope] : 0.0;
    if (cell.gate) {
      next.beta0 = (cell.h - cell.a * theta) / cell.P0 + rng->normal() / std::sqrt(cell.P0);
    } else {
      draw_prior(g0, rng, &next.beta0);
    }
    Vec m;
    latent_means(C, cell.gate ? next.beta0 : 0.0, cell.slope, theta, m);
    next.z.assign(m.size(), 0.0);
    for (std::size_t s = 0; s < m.size(); ++s) next.z[s] = m[s] + rng->normal() / std::sqrt(_lambda[s]);
    if (!std::isfinite(next.beta0) || !std::isfinite(next.theta[0]) || !std::isfinite(next.theta[1]) ||
        !all_finite(next.z.data(), next.z.size())) return;
    // the acceptance from the inner target and the proposal densities
    const FamilyState current = family_state();
    Shift sh;
    bool ok;
    const long double delta = family_delta(next, sh, &ok);
    const double log_q_current = family_log_q(C, cells, log_total, current);
    const double log_q_next = family_log_q(C, cells, log_total, next);
    const double log_ratio = static_cast<double>(delta + log_q_current - log_q_next);
    if (dry) {
      Cell const *from = nullptr;
      for (Cell const &c : cells) if (c.label == current.label && c.gate == current.gate) from = &c;
      const double analytic = from == nullptr ? negative_infinity :
        family_log_F(cell, theta) - family_log_F(*from, from->slope >= 0 ? current.theta[from->slope] : 0.0);
      probe_record(probe_family, analytic, log_ratio, (from != nullptr && from->prior_cell ? 1.0 : 0.0) +
                   (cell.prior_cell ? 2.0 : 0.0), static_cast<double>(cell.label));
      return;
    }
    if (!ok || !std::isfinite(log_ratio) || !metropolis(rng, log_ratio)) return;
    count(counters().bias_accepted[bias_gaussian]);
    family_apply(next, sh, delta);
  }

  // ---- S: one slice update of the active slope on its collapsed conditional
  // log p_b(theta) + r theta - q theta^2 / 2 (the current cell's q and r, the
  // latents and, gate on, the intercept integrated), then the exact
  // conditional draws of the intercept (gate on) and the latents; accepted
  // with the ratio of the inner target to the collapsed law times the
  // conditional densities (one when the Gaussian algebra is exact).
  void slope_move(RNG *rng, bool dry = false)
  {
    BiasSpec const &B = bspec();
    if (_selection || !B.gaussian(_label)) return;
    const int b = B.slope_of[static_cast<std::size_t>(_label - 1)];
    if (b < 0) return;
    Collapse C;
    if (!collapse(C)) return;
    const bool gate = active(B.intercept);
    const Cell cell = make_cell(C, _label, gate, 0.0);
    const SlopePrior prior = slope_prior(B.slope[b], _chain);
    // the collapsed law (q floored at zero against rounding: it only shapes the
    // proposal, the acceptance below uses the inner target)
    const double q = std::max(cell.q, 0.0);
    auto f = [&](double theta) {
      if (!(theta >= 0) || !std::isfinite(theta)) return negative_infinity;
      const double lp = slope_log_prior(prior, theta);
      return std::isfinite(lp) ? lp + cell.r * theta - 0.5 * q * theta * theta : negative_infinity;
    };
    const double theta_old = _theta[b];
    const double f_old = f(theta_old);
    if (!std::isfinite(f_old)) return;
    if (!dry) count(counters().bias_proposed[bias_slope]);
    // the width: 2.5 SDs of the collapsed Gaussian, the prior's scale when the
    // slope is weakly identified (a function of the cell, not of theta)
    double width = q > 0 ? 2.5 / std::sqrt(q) : 2.5;
    if (!std::isfinite(width) || !(width > 0) || width > 1e6) width = 2.5;
    std::vector<double> u(1, theta_old), direction(1, 1.0);
    const slice::SliceResult result = slice::slice_along(
      u, direction, width, f_old, [&](std::vector<double> const &point) { return f(point[0]); },
      [rng]() { return rng->uniform(); }, [rng]() { return rng->exponential(); });
    if (!result.accepted) return;
    const double theta_new = u[0];
    FamilyState next = family_state();
    next.theta[b] = theta_new;
    // the conditional draws
    auto beta_mean = [&](double theta) { return (cell.h - cell.a * theta) / cell.P0; };
    if (gate) next.beta0 = beta_mean(theta_new) + rng->normal() / std::sqrt(cell.P0);
    Vec m_new, m_old;
    latent_means(C, gate ? next.beta0 : 0.0, b, theta_new, m_new);
    for (std::size_t s = 0; s < m_new.size(); ++s) next.z[s] = m_new[s] + rng->normal() / std::sqrt(_lambda[s]);
    if (!std::isfinite(next.beta0) || !all_finite(next.z.data(), next.z.size())) return;
    FamilyState const current = family_state();
    latent_means(C, gate ? current.beta0 : 0.0, b, theta_old, m_old);
    auto conditional = [&](FamilyState const &x, double theta, Vec const &m) {
      long double out = latents_log_density(x.z, m);
      if (gate) {
        const double e = (x.beta0 - beta_mean(theta)) * std::sqrt(cell.P0);
        out += -0.5 * e * e - 0.91893853320467274178 + 0.5 * std::log(cell.P0);
      }
      return out;
    };
    Shift sh;
    bool ok;
    const long double delta = family_delta(next, sh, &ok);
    const double log_ratio = static_cast<double>(delta + (f_old + conditional(current, theta_old, m_old)) -
                                                 (result.log_density + conditional(next, theta_new, m_new)));
    if (dry) {
      probe_record(probe_slope, 0.0, log_ratio, 0.0, static_cast<double>(_label));
      return;
    }
    record_max(counters().max_stage2_gaussian, std::fabs(log_ratio));
    if (!ok || !std::isfinite(log_ratio) || !metropolis(rng, log_ratio)) return;
    count(counters().bias_accepted[bias_slope]);
    family_apply(next, sh, delta);
  }

  // ---- the bias slots of the sweep (L1, S and L2 above this point)

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
    case bias_gaussian:
      family_move(rng);
      break;
    case bias_slope:
      slope_move(rng);
      break;
    default:
      break;
    }
  }
