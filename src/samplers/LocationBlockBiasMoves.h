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
  static constexpr int probe_family = 7, probe_slope = 8, probe_swap = 9, probe_label = 10;

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
      bias_tally(bias_dormant, 0);
      const SlopePrior prior = slope_prior(spec, _chain);
      const double draw = prior.dist->randomSample(prior.par, prior.lower, prior.upper, rng);
      const double log_new = slope_log_prior(prior, draw);
      const double log_old = slope_log_prior(prior, _theta[b]);
      if (!std::isfinite(draw) || !std::isfinite(log_new) || !std::isfinite(log_old)) continue;
      bias_tally(bias_dormant, 2);
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
    bias_tally(bias_label, 0);
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
    bias_tally(bias_label, 2);
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
    if (!dry) bias_tally(bias_gaussian, 0);
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
    bias_tally(bias_gaussian, 2);
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
    if (!dry) bias_tally(bias_slope, 0);
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
    bias_tally(bias_slope, 2);
    family_apply(next, sh, delta);
  }

  // ---- L2: the cross-family swap (wp3-spec.md, L2)
  //
  // Coordinates: the intercept's physical value t, the active slope theta_b
  // of a PET or PEESE label, and E = (the intercept line's partners, every
  // latent). The shears K_t (the intercept's line, LocationBlock::build_line)
  // and K_b (the conditional-mean shear of E per unit physical slope b, the
  // same partner set) are functions of the precisions, the other gates and
  // the non-members only; E_0 = E - K_t t - K_b theta_b is what the swap holds
  // fixed (a unit-Jacobian bijection). Both directions evaluate the label
  // proposal, the evidences and the laws at the base point (t = 0, no slope,
  // E = E_0) from unchanged state.

  struct SwapGeometry {
    Line L;                           // the intercept's line (partners, latent shear)
    Vec Kt, Kb[2];                    // shears over E = (partners, latents)
    bool slope[2] = {false, false};
    std::size_t nq = 0, nz = 0;
  };

  // The base-point quantities of one attempt (functions of E_0 and unchanged
  // state only).
  struct SwapBase {
    Vec E0;                           // partners then latents
    double G0 = 0;                    // sum of the children's normal densities at the base point
    // the Gaussian plane in (t, theta_b): G(x) = G0 + eta'x - x'Pi x / 2
    double eta_t = 0, Pi_tt = 0, eta_b[2] = {0, 0}, Pi_bb[2] = {0, 0}, Pi_tb[2] = {0, 0};
    Vec log_e_off, log_e_on;          // per label: log evidences (-Inf without mass or gate)
    Vec screen;                       // per selection label: the row-marginal screen
    // the Laplace law of t at the evaluated selection labels
    Vec fit_mean, fit_precision;
    std::vector<char> fit_ok;
  };

  // The swap geometry at the current state; false when a law is singular.
  bool swap_geometry(SwapGeometry &G)
  {
    Shared const &S = s();
    BiasSpec const &B = bspec();
    const unsigned p = S.p;
    if (!build_line({B.intercept}, G.L)) return false;
    Line const &L = G.L;
    G.nq = L.partners.size();
    G.nz = S.latents.size();
    const std::size_t ne = G.nq + G.nz;
    G.Kt.assign(ne, 0.0);
    for (std::size_t q = 0; q < G.nq; ++q) G.Kt[q] = L.KQ[q];
    for (std::size_t z = 0; z < L.lats.size(); ++z) G.Kt[G.nq + L.lats[z]] = L.KZ[z];
    // the slope shears: K_Q = -P_QQ^-1 X_Q' S^-1 s_b, K_Z = -(d'W s_b + u_Q K_Q) / lambda
    Vec PQ;
    if (G.nq > 0) {
      PQ.assign(G.nq * G.nq, 0.0);
      for (std::size_t i = 0; i < G.nq; ++i) for (std::size_t j = 0; j < G.nq; ++j) PQ[i * G.nq + j] = _P[L.partners[i] * p + L.partners[j]];
      if (!cholesky(PQ, G.nq)) return false;
    }
    std::vector<Vec> columns(G.nq, Vec(S.rows, 0.0));
    std::vector<double const *> U(G.nq);
    for (std::size_t q = 0; q < G.nq; ++q) {
      for (unsigned row = 0; row < S.rows; ++row) columns[q][row] = S.X[row * p + L.partners[q]];
      U[q] = columns[q].data();
    }
    for (int b = 0; b < 2; ++b) {
      G.slope[b] = B.slope[b].source != nullptr;
      G.Kb[b].assign(ne, 0.0);
      if (!G.slope[b]) continue;
      Vec cross;
      if (G.nq > 0) {
        cross_S(U, {B.slope[b].column.data()}, cross);
        for (std::size_t q = 0; q < G.nq; ++q) cross[q] = -cross[q];
        solve_lower(PQ, G.nq, cross.data());
        solve_lower_transpose(PQ, G.nq, cross.data());
        for (std::size_t q = 0; q < G.nq; ++q) G.Kb[b][q] = cross[q];
      }
      for (unsigned lat = 0; lat < G.nz; ++lat) {
        long double v = 0;
        for (unsigned c : S.latent_children[lat]) {
          ChildSpec const &child = S.children[c];
          double const *W = &_W[child.woff];
          for (unsigned i = 0; i < child.n; ++i) {
            long double wd = 0;
            for (unsigned j = 0; j < child.n; ++j) wd += static_cast<long double>(W[i * child.n + j]) * _d[child.row0 + j];
            v += wd * B.slope[b].column[child.row0 + i];
          }
        }
        for (std::size_t q = 0; q < G.nq; ++q) v += static_cast<long double>(_u[lat * p + L.partners[q]]) * G.Kb[b][q];
        G.Kb[b][G.nq + lat] = static_cast<double>(-v / _lambda[lat]);
      }
    }
    return all_finite(G.Kt.data(), ne) && all_finite(G.Kb[0].data(), ne) && all_finite(G.Kb[1].data(), ne);
  }

  // The current E = (partners, latents).
  Vec swap_E(SwapGeometry const &G) const
  {
    Vec E(G.nq + G.nz, 0.0);
    for (std::size_t q = 0; q < G.nq; ++q) E[q] = _b[G.L.partners[q]];
    for (std::size_t z = 0; z < G.nz; ++z) E[G.nq + z] = _z[z];
    return E;
  }

  // A shift of the intercept by dt, the physical slope b by dth (b < 0: none)
  // and E by dE, prepared (the partners' slab priors and the latents' priors
  // in p1, p2).
  void swap_shift(SwapGeometry const &G, double dt, int b, double dth, Vec const &dE, Shift &sh)
  {
    BiasSpec const &B = bspec();
    sh.clear();
    sh.coords.push_back(s().groups[B.intercept].coord0);
    sh.dcoef.push_back(dt);
    for (std::size_t q = 0; q < G.nq; ++q) {
      sh.coords.push_back(G.L.partners[q]);
      sh.dcoef.push_back(dE[q]);
    }
    for (std::size_t z = 0; z < G.nz; ++z) {
      sh.lats.push_back(static_cast<unsigned>(z));
      sh.dlat.push_back(dE[G.nq + z]);
    }
    if (b >= 0) {
      sh.slopes.push_back(b);
      sh.dslope.push_back(dth);
    }
    prepare_shift(sh, G.L.partners);
  }

  // Moves the working coordinates (not the label) by a prepared swap shift.
  void swap_move_state(Shift const &sh)
  {
    Shared const &S = s();
    move_means(sh, 1.0);
    for (std::size_t i = 0; i < sh.coords.size(); ++i) {
      const unsigned k = sh.coords[i];
      _b[k] += sh.dcoef[i];
      if (active(S.coord_group[k])) _beta[k] = _b[k];
    }
    for (std::size_t i = 0; i < sh.lats.size(); ++i) _z[sh.lats[i]] += sh.dlat[i];
  }

  // Sets the working label (its graph nodes, selection status; values lazily).
  void swap_set_label(int l)
  {
    Shared const &S = s();
    set_graph_label(l);
    _label = l;
    _selection = false;
    for (unsigned c = 0; c < S.children.size(); ++c) {
      _selected[c] = is_selected(S.children[c]);
      if (_selected[c]) _selection = true;
      _value_ok[c] = 0;
      _fit_ok[c] = 0;
    }
    _line.valid = false;
  }

  // The working state, saved and restored around the base-point evaluations.
  struct Snapshot {
    Vec b, beta, z, mean, r, wr, h, value, fit_value;
    std::vector<char> g, value_ok, fit_ok, selected, dirty;
    bool selection = false;
    int label = 0;
  };

  void take(Snapshot &x) const
  {
    x.b = _b; x.beta = _beta; x.z = _z; x.mean = _mean; x.r = _r; x.wr = _wr; x.h = _h;
    x.value = _value; x.fit_value = _fit_value; x.g = _g; x.value_ok = _value_ok; x.fit_ok = _fit_ok;
    x.selected = _selected; x.dirty = _dirty; x.selection = _selection; x.label = _label;
  }

  void put(Snapshot const &x)
  {
    _b = x.b; _beta = x.beta; _z = x.z; _mean = x.mean; _r = x.r; _wr = x.wr; _h = x.h;
    _value = x.value; _fit_value = x.fit_value; _g = x.g; _value_ok = x.value_ok; _fit_ok = x.fit_ok;
    _selected = x.selected; _dirty = x.dirty; _selection = x.selection; _label = x.label;
    set_graph_label(_label);
    _line.valid = false;
  }

  // The inner target's term of child c at the graph's label and the
  // working mean (selection children their surrogate or exact term, the
  // others their exact normal density).
  double term_now(unsigned c)
  {
    if (!gaussian_child(c)) {
      ensure_value(c);
      return _value[c];
    }
    ChildSpec const &child = s().children[c];
    const double v = exact_with_mean(c, &_mean[child.row0]);
    if (std::isnan(v) || v == std::numeric_limits<double>::infinity()) {
      throwNodeError(child.node, "Location block sampler: the log density of this likelihood node is not defined at a label");
    }
    return v;
  }

  // sum over the rows of a prepared shift's children of delta_a' W delta_b
  // for the shift sa's deltas and sb's weighted deltas.
  double shift_cross(Shift const &sa, Shift const &sb)
  {
    std::vector<int> &position = _tilt_position;
    for (std::size_t q = 0; q < sb.children.size(); ++q) position[sb.children[q]] = static_cast<int>(q);
    long double out = 0;
    for (std::size_t q = 0; q < sa.children.size(); ++q) {
      const int other = position[sa.children[q]];
      if (other < 0) continue;
      const unsigned n = s().children[sa.children[q]].n;
      for (unsigned i = 0; i < n; ++i) {
        out += static_cast<long double>(sa.delta[sa.offset[q] + i]) * sb.wdelta[sb.offset[static_cast<std::size_t>(other)] + i];
      }
    }
    for (unsigned c : sb.children) position[c] = -1;
    return static_cast<double>(out);
  }

  // K_a' B_E K_b over E (the partners' slab precisions, unit latents).
  double prior_cross(SwapGeometry const &G, Vec const &Ka, Vec const &Kb) const
  {
    Shared const &S = s();
    long double out = 0;
    for (std::size_t i = 0; i < G.nq; ++i) {
      const unsigned ki = G.L.partners[i];
      GroupSpec const &gi = S.groups[S.coord_group[ki]];
      for (std::size_t j = 0; j < G.nq; ++j) {
        const unsigned kj = G.L.partners[j];
        if (S.coord_group[kj] != S.coord_group[ki]) continue;
        out += static_cast<long double>(Ka[i]) * gi.B[(ki - gi.coord0) * gi.d + (kj - gi.coord0)] * Kb[j];
      }
    }
    for (std::size_t z = 0; z < G.nz; ++z) out += static_cast<long double>(Ka[G.nq + z]) * Kb[G.nq + z];
    return static_cast<double>(out);
  }

  // The law of a Gaussian label's cell along the plane: the slope's
  // quadratic (k + r theta - q theta^2 / 2 after integrating t when the gate
  // is on) and t's conditional law given theta (gate on).
  struct PlaneCell {
    int slope = -1;
    bool gate = true;
    double q = 0, r = 0, k = 0;       // in theta
    double Pt = 0, mt = 0, ct = 0;    // gate on: t | theta ~ N(mt - ct theta, 1 / Pt)
    bool prior_cell = false;
    double anchor = 0, log_e = negative_infinity;  // log evidence (label prior etc. excluded)
  };

  PlaneCell plane_cell(SwapBase const &base, int slope, bool gate)
  {
    BiasSpec const &B = bspec();
    GroupSpec const &g0 = s().groups[B.intercept];
    const double B0 = g0.B[0], m0 = g0.m0[0];
    PlaneCell cell;
    cell.slope = slope;
    cell.gate = gate;
    const double eta_b = slope >= 0 ? base.eta_b[slope] : 0.0;
    const double Pi_bb = slope >= 0 ? base.Pi_bb[slope] : 0.0;
    const double Pi_tb = slope >= 0 ? base.Pi_tb[slope] : 0.0;
    if (gate) {
      cell.Pt = base.Pi_tt + B0;
      const double lin = base.eta_t + B0 * m0;
      cell.mt = lin / cell.Pt;
      cell.ct = Pi_tb / cell.Pt;
      cell.q = Pi_bb - Pi_tb * Pi_tb / cell.Pt;
      cell.r = eta_b - Pi_tb * lin / cell.Pt;
      cell.k = 0.5 * std::log(B0 / cell.Pt) - 0.5 * B0 * m0 * m0 + lin * lin / (2.0 * cell.Pt);
    } else {
      cell.q = Pi_bb;
      cell.r = eta_b;
      cell.k = 0.0;
    }
    cell.log_e = base.G0 + cell.k;
    if (slope < 0) return cell;
    const SlopePrior prior = slope_prior(B.slope[slope], _chain);
    double log_anchor = negative_infinity;
    if (cell.q > 0 && cell.q >= 1e-12 * Pi_bb) {
      const double mean = cell.r / cell.q, sd = 1.0 / std::sqrt(cell.q);
      const double alpha = -mean / sd;
      const double mills = std::exp(-0.5 * alpha * alpha - 0.91893853320467274178 - log_phi_upper(-alpha));
      cell.anchor = mean + sd * mills;
      if (cell.anchor > 0 && std::isfinite(cell.anchor)) log_anchor = slope_log_prior(prior, cell.anchor);
    }
    if (!std::isfinite(log_anchor)) {
      cell.prior_cell = true;
      return cell;
    }
    cell.log_e += log_anchor + 0.5 * std::log(2.0 * 3.14159265358979323846 / cell.q) + cell.r * cell.r / (2.0 * cell.q) +
      log_phi_upper(cell.r / std::sqrt(cell.q));
    return cell;
  }

  // log density of a Gaussian cell's coordinates: the slope (truncated
  // normal or prior) and t given the slope (gate on) or the dormant slab.
  double plane_log_density(PlaneCell const &cell, double theta, double beta0)
  {
    BiasSpec const &B = bspec();
    long double out = 0;
    if (cell.slope >= 0) {
      out += cell.prior_cell ? slope_log_prior(slope_prior(B.slope[cell.slope], _chain), theta) :
        log_truncated_normal_positive(theta, cell.r / cell.q, 1.0 / std::sqrt(cell.q));
    }
    if (cell.gate) {
      const double e = (beta0 - (cell.mt - cell.ct * theta)) * std::sqrt(cell.Pt);
      out += -0.5 * e * e - 0.91893853320467274178 + 0.5 * std::log(cell.Pt);
    } else {
      out += log_slab(s().groups[B.intercept], &beta0);
    }
    return static_cast<double>(out);
  }

  // log density of a selection cell's coordinates: t from the Laplace law
  // (gate on) or the dormant slab.
  double fit_log_density(SwapBase const &base, int l, bool gate, double beta0)
  {
    BiasSpec const &B = bspec();
    if (!gate) return log_slab(s().groups[B.intercept], &beta0);
    const std::size_t i = static_cast<std::size_t>(l - 1);
    const double e = (beta0 - base.fit_mean[i]) * std::sqrt(base.fit_precision[i]);
    return -0.5 * e * e - 0.91893853320467274178 + 0.5 * std::log(base.fit_precision[i]);
  }

  // The row-marginal screen of selection label l at the working means (the
  // graph at l): sum over the selection rows of power (log omega_obs - log A).
  double label_screen()
  {
    Shared const &S = s();
    long double total = 0;
    RowModel row;
    for (unsigned c = 0; c < S.children.size(); ++c) {
      if (!_selected[c]) continue;
      ChildSpec const &child = S.children[c];
      for (unsigned i = 0; i < child.n; ++i) {
        const double mu = _mean[child.row0 + i];
        double log_a, d1, d2, power;
        int bin;
        double omega_obs;
        if (child.type == ChildType::step_switch && _step[c].valid) {
          StepChild const &st = _step[c];
          bin = st.bin;
          if (bin == 0) continue;
          power = st.power;
          omega_obs = st.omega[bin - 1];
          if (!step_row_log_normalizer(step_law(c), mu, &log_a, &d1, &d2)) return negative_infinity;
        } else {
          if (!row_model(child, i, row) || row.bin == 0) continue;
          bin = row.bin;
          power = row.power;
          omega_obs = row.omega[bin - 1];
          if (!row_log_normalizer(mu, row.sd, row.se, row.sign, row.omega, row.lower.data(), row.upper.data(), row.bins,
                                  &log_a, &d1, &d2)) return negative_infinity;
        }
        if (!(omega_obs > 0)) return negative_infinity;
        total += static_cast<long double>(power) * (std::log(omega_obs) - log_a);
      }
    }
    return static_cast<double>(total);
  }

  // The base-point evaluations: moves the working state to the base point,
  // evaluates the Gaussian plane, every Gaussian label's evidences, every
  // selection label's screen and the absolute evidences and Laplace laws of
  // the selection labels in `fits`, then restores the working state.
  bool swap_base(SwapGeometry const &G, std::vector<int> const &fits, SwapBase &base)
  {
    Shared const &S = s();
    BiasSpec const &B = bspec();
    GroupSpec const &g0 = S.groups[B.intercept];
    const int levels = B.levels;
    const int b_now = B.slope_of[static_cast<std::size_t>(_label - 1)];
    const double t = _b[g0.coord0];
    const double th = b_now >= 0 ? _theta[b_now] : 0.0;
    const Vec E = swap_E(G);
    base.E0.assign(E.size(), 0.0);
    for (std::size_t e = 0; e < E.size(); ++e) {
      base.E0[e] = E[e] - G.Kt[e] * t - (b_now >= 0 ? G.Kb[b_now][e] * th : 0.0);
    }
    Vec dE(E.size());
    for (std::size_t e = 0; e < E.size(); ++e) dE[e] = base.E0[e] - E[e];
    Snapshot saved;
    take(saved);
    bool ok = true;
    try {
      Shift to_base;
      swap_shift(G, -t, b_now, -th, dE, to_base);
      swap_move_state(to_base);
      for (unsigned c = 0; c < S.children.size(); ++c) _value_ok[c] = 0;
      rebuild_h();
      // the children's normal densities at the base point (a Gaussian label)
      int gaussian_label = 0;
      for (int l = 1; l <= levels && gaussian_label == 0; ++l) {
        if (B.gaussian(l) && B.log_prior[static_cast<std::size_t>(l - 1)] > negative_infinity) gaussian_label = l;
      }
      swap_set_label(gaussian_label);
      long double g0_sum = 0;
      for (unsigned c = 0; c < S.children.size(); ++c) g0_sum += term_now(c);
      base.G0 = static_cast<double>(g0_sum);
      // the plane: t and each slope
      Shift st;
      swap_shift(G, 1.0, -1, 0.0, G.Kt, st);
      long double qt1 = 0, qt2 = 0;
      for (std::size_t q = 0; q < st.children.size(); ++q) { qt1 += st.q1[q]; qt2 += st.q2[q]; }
      base.eta_t = static_cast<double>(qt1 + st.p1);
      base.Pi_tt = static_cast<double>(qt2 + st.p2);
      for (int b = 0; b < 2; ++b) {
        if (!G.slope[b]) continue;
        Shift sb;
        swap_shift(G, 0.0, b, 1.0, G.Kb[b], sb);
        long double q1 = 0, q2 = 0;
        for (std::size_t q = 0; q < sb.children.size(); ++q) { q1 += sb.q1[q]; q2 += sb.q2[q]; }
        base.eta_b[b] = static_cast<double>(q1 + sb.p1);
        base.Pi_bb[b] = static_cast<double>(q2 + sb.p2);
        base.Pi_tb[b] = shift_cross(st, sb) + prior_cross(G, G.Kt, G.Kb[b]);
      }
      // test hook (probe): the fitted quadratic of the plane against the
      // children's own normal densities and the priors of E at a few points
      if (S.probe) {
        double worst = 0.0, scale = 1.0;
        const double points[3][2] = {{0.1, 0.2}, {-0.2, 0.4}, {0.3, 0.05}};
        for (int b = -1; b < 2; ++b) {
          if (b >= 0 && !G.slope[b]) continue;
          for (auto const &point : points) {
            const double t = point[0], th = b >= 0 ? point[1] : 0.0;
            Vec dE(G.Kt.size());
            for (std::size_t e = 0; e < dE.size(); ++e) dE[e] = G.Kt[e] * t + (b >= 0 ? G.Kb[b][e] * th : 0.0);
            Shift move;
            swap_shift(G, t, b, th, dE, move);
            const long double prior = prior_delta(move, 1.0);
            Snapshot at_base;
            take(at_base);
            swap_move_state(move);
            for (unsigned c = 0; c < S.children.size(); ++c) _value_ok[c] = 0;
            long double direct = prior - base.G0;
            for (unsigned c = 0; c < S.children.size(); ++c) direct += term_now(c);
            put(at_base);
            const double eta_b = b >= 0 ? base.eta_b[b] : 0.0, Pi_bb = b >= 0 ? base.Pi_bb[b] : 0.0;
            const double Pi_tb = b >= 0 ? base.Pi_tb[b] : 0.0;
            const double quadratic = base.eta_t * t + eta_b * th -
              0.5 * (base.Pi_tt * t * t + 2.0 * Pi_tb * t * th + Pi_bb * th * th);
            worst = std::max(worst, std::fabs(static_cast<double>(direct) - quadratic));
            scale = std::max(scale, std::fabs(quadratic));
          }
        }
        probe_raw(probe_swap, 2.0, worst, scale, 0.0, 0.0);
      }
      // evidences of every label (log prior, nested term and gate prior included)
      base.log_e_off.assign(static_cast<std::size_t>(levels), negative_infinity);
      base.log_e_on.assign(static_cast<std::size_t>(levels), negative_infinity);
      base.screen.assign(static_cast<std::size_t>(levels), negative_infinity);
      base.fit_mean.assign(static_cast<std::size_t>(levels), 0.0);
      base.fit_precision.assign(static_cast<std::size_t>(levels), 0.0);
      base.fit_ok.assign(static_cast<std::size_t>(levels), 0);
      const bool gated = g0.gate != nullptr;
      for (int l = 1; l <= levels; ++l) {
        const std::size_t i = static_cast<std::size_t>(l - 1);
        if (!(B.log_prior[i] > negative_infinity)) continue;
        set_graph_label(l);
        const double head = B.log_prior[i] + (S.mutant == 9 ? 0.0 : nested_term());
        if (B.gaussian(l)) {
          const int b = B.slope_of[i];
          if (gated) base.log_e_off[i] = head + g0.log_off + plane_cell(base, b, false).log_e;
          base.log_e_on[i] = head + (gated ? g0.log_on : 0.0) + plane_cell(base, b, true).log_e;
        } else {
          swap_set_label(l);
          base.screen[i] = label_screen();
        }
      }
      for (int l : fits) {
        const std::size_t i = static_cast<std::size_t>(l - 1);
        swap_set_label(l);
        long double absolute = 0;
        for (unsigned c = 0; c < S.children.size(); ++c) absolute += term_now(c);
        const double head = B.log_prior[i] + (S.mutant == 9 ? 0.0 : nested_term());
        // the absolute evidence at the base point with the gate off, and the
        // Laplace odds of the gate along the intercept's line (prior odds
        // included; LocationBlock::group_fit)
        const double off = S.mutant == 10 ? 0.0 : head + static_cast<double>(absolute);
        if (!group_line(B.intercept)) { ok = false; break; }
        group_fit();
        if (!_line.fit.ok || !std::isfinite(_line.log_odds)) { ok = false; break; }
        if (gated) {
          base.log_e_off[i] = off + g0.log_off;
          base.log_e_on[i] = off + g0.log_off + _line.log_odds;
        } else {
          base.log_e_on[i] = off + _line.log_odds;
        }
        base.fit_mean[i] = _line.fit.xhat[0];
        base.fit_precision[i] = _line.fit.Pt[0];
        base.fit_ok[i] = 1;
      }
    } catch (...) {
      put(saved);
      throw;
    }
    put(saved);
    return ok;
  }

  static double log_sum(double a, double b)
  {
    if (!(a > negative_infinity)) return b;
    if (!(b > negative_infinity)) return a;
    const double m = std::max(a, b);
    return m + std::log(std::exp(a - m) + std::exp(b - m));
  }

  // The label proposal of each family: 5% uniform and 95% in proportion to
  // the Gaussian labels' evidences (q_G) or the selection labels' screens
  // (q_S); log probabilities, -Inf outside the family.
  void swap_label_proposals(SwapBase const &base, Vec &log_qG, Vec &log_qS)
  {
    BiasSpec const &B = bspec();
    const int levels = B.levels;
    log_qG.assign(static_cast<std::size_t>(levels), negative_infinity);
    log_qS.assign(static_cast<std::size_t>(levels), negative_infinity);
    for (int family = 0; family < 2; ++family) {
      Vec &out = family == 0 ? log_qG : log_qS;
      std::vector<std::size_t> members;
      Vec score;
      for (int l = 1; l <= levels; ++l) {
        const std::size_t i = static_cast<std::size_t>(l - 1);
        if (!(B.log_prior[i] > negative_infinity) || B.gaussian(l) != (family == 0)) continue;
        members.push_back(i);
        score.push_back(family == 0 ? log_sum(base.log_e_off[i], base.log_e_on[i]) : base.screen[i]);
      }
      if (members.empty()) continue;
      double maximum = negative_infinity;
      for (double v : score) maximum = std::max(maximum, v);
      long double total = 0;
      for (double v : score) total += std::isfinite(maximum) && v > negative_infinity ? std::exp(static_cast<long double>(v - maximum)) : 0.0L;
      const double uniform = 0.05 / static_cast<double>(members.size());
      for (std::size_t m = 0; m < members.size(); ++m) {
        const double share = total > 0 && score[m] > negative_infinity ?
          static_cast<double>(std::exp(static_cast<long double>(score[m] - maximum)) / total) : 0.0;
        out[members[m]] = std::log(uniform + (total > 0 ? 0.95 * share : 0.95 / static_cast<double>(members.size())));
      }
    }
  }

  // A state of the swap's coordinates: the label, the intercept's gate and
  // slab value (its physical value when the gate is on), both slopes, E.
  struct SwapState {
    int label = 0;
    bool gate = true;
    double beta0 = 0;
    double theta[2] = {0, 0};
    Vec E;
  };

  SwapState swap_current(SwapGeometry const &G) const
  {
    BiasSpec const &B = bspec();
    SwapState x;
    x.label = _label;
    x.gate = active(B.intercept);
    x.beta0 = _beta[s().groups[B.intercept].coord0];
    x.theta[0] = _theta[0];
    x.theta[1] = _theta[1];
    x.E = swap_E(G);
    return x;
  }

  // The stage-1 log ratio of a move from label `from` to label `to`: the
  // ratio of the cells' summed evidences times the label Hastings factor.
  double swap_stage1(SwapBase const &base, Vec const &log_qG, Vec const &log_qS, int from, int to) const
  {
    BiasSpec const &B = bspec();
    const std::size_t it = static_cast<std::size_t>(to - 1), ifr = static_cast<std::size_t>(from - 1);
    const bool from_selection = !B.gaussian(from);
    const double hastings = s().mutant == 6 ? 0.0 :
      (from_selection ? log_qS : log_qG)[ifr] - (from_selection ? log_qG : log_qS)[it];
    return log_sum(base.log_e_off[it], base.log_e_on[it]) - log_sum(base.log_e_off[ifr], base.log_e_on[ifr]) + hastings;
  }

  // log q(to | from): the label proposal of from's family, the gate share of
  // the target label, the target cell's coordinates and the prior of the
  // slope that becomes dormant (mutants 6, 7 and 8 omit one of them).
  double swap_log_q(SwapBase const &base, Vec const &log_qG, Vec const &log_qS, SwapState const &from,
                    SwapState const &to)
  {
    Shared const &S = s();
    BiasSpec const &B = bspec();
    GroupSpec const &g0 = S.groups[B.intercept];
    const std::size_t it = static_cast<std::size_t>(to.label - 1), ifr = static_cast<std::size_t>(from.label - 1);
    long double out = 0;
    if (S.mutant != 6) out += (B.gaussian(from.label) ? log_qS : log_qG)[it];
    if (g0.gate != nullptr && S.mutant != 7) {
      out += (to.gate ? base.log_e_on[it] : base.log_e_off[it]) - log_sum(base.log_e_off[it], base.log_e_on[it]);
    }
    if (B.gaussian(to.label)) {
      const int b = B.slope_of[it];
      out += plane_log_density(plane_cell(base, b, to.gate), b >= 0 ? to.theta[b] : 0.0, to.beta0);
    } else {
      out += fit_log_density(base, to.label, to.gate, to.beta0);
    }
    const int slope_out = B.slope_of[ifr];
    if (slope_out >= 0 && S.mutant != 8) out += slope_log_prior(slope_prior(B.slope[slope_out], _chain), to.theta[slope_out]);
    return static_cast<double>(out);
  }

  // Draws the proposal from the current state: the target label's gate and
  // cell coordinates (a Gaussian label's slope and then t, a selection label's
  // t from its Laplace law, the dormant slab from its prior), the leaving
  // slope from its prior, and E = E_0 + K_t t + K_b theta.
  bool swap_draw(SwapGeometry const &G, SwapBase const &base, SwapState const &from, int target, RNG *rng,
                 SwapState &to)
  {
    BiasSpec const &B = bspec();
    GroupSpec const &g0 = s().groups[B.intercept];
    const std::size_t it = static_cast<std::size_t>(target - 1);
    to = from;
    to.label = target;
    to.gate = true;
    if (g0.gate != nullptr) {
      to.gate = rng->uniform() < std::exp(base.log_e_on[it] - log_sum(base.log_e_off[it], base.log_e_on[it]));
    }
    const int slope_in = B.slope_of[it];
    const int slope_out = B.slope_of[static_cast<std::size_t>(from.label - 1)];
    if (B.gaussian(target)) {
      const PlaneCell cell = plane_cell(base, slope_in, to.gate);
      double theta = 0.0;
      if (slope_in >= 0) {
        const SlopePrior prior = slope_prior(B.slope[slope_in], _chain);
        theta = cell.prior_cell ? prior.dist->randomSample(prior.par, prior.lower, prior.upper, rng) :
          truncated_normal_positive(cell.r / cell.q, 1.0 / std::sqrt(cell.q), rng);
        to.theta[slope_in] = theta;
      }
      if (to.gate) to.beta0 = cell.mt - cell.ct * theta + rng->normal() / std::sqrt(cell.Pt);
      else draw_prior(g0, rng, &to.beta0);
    } else {
      if (to.gate) to.beta0 = base.fit_mean[it] + rng->normal() / std::sqrt(base.fit_precision[it]);
      else draw_prior(g0, rng, &to.beta0);
    }
    if (slope_out >= 0) {
      const SlopePrior prior = slope_prior(B.slope[slope_out], _chain);
      to.theta[slope_out] = prior.dist->randomSample(prior.par, prior.lower, prior.upper, rng);
    }
    const double t = to.gate ? to.beta0 : 0.0;
    for (std::size_t e = 0; e < to.E.size(); ++e) {
      to.E[e] = base.E0[e] + G.Kt[e] * t + (slope_in >= 0 ? G.Kb[slope_in][e] * to.theta[slope_in] : 0.0);
    }
    return std::isfinite(to.beta0) && std::isfinite(to.theta[0]) && std::isfinite(to.theta[1]) &&
      all_finite(to.E.data(), to.E.size());
  }

  // Moves the working state from `from` (the current state) to `to` and
  // returns the change of the inner target: every child's absolute term at
  // both labels and means, the label prior, the nested node's term (mutant 9
  // omits it), the gate prior, the intercept's slab, both slope priors and the
  // priors of E.
  long double swap_move(SwapGeometry const &G, SwapState const &from, SwapState const &to)
  {
    Shared const &S = s();
    BiasSpec const &B = bspec();
    GroupSpec const &g0 = S.groups[B.intercept];
    long double old_terms = 0, new_terms = 0;
    for (unsigned c = 0; c < S.children.size(); ++c) old_terms += term_now(c);
    const double nested_old = S.mutant == 9 ? 0.0 : nested_term();
    const double t_old = from.gate ? from.beta0 : 0.0, t_new = to.gate ? to.beta0 : 0.0;
    Vec dE(to.E.size());
    for (std::size_t e = 0; e < dE.size(); ++e) dE[e] = to.E[e] - from.E[e];
    Shift sh;
    swap_shift(G, t_new - t_old, -1, 0.0, dE, sh);
    sh.slopes.clear();
    sh.dslope.clear();
    for (int b = 0; b < 2; ++b) {
      if (B.slope[b].source == nullptr) continue;
      const double now = physical_at(from.label, b);
      const double next = B.slope[b].label == to.label ? to.theta[b] : 0.0;
      if (next != now) {
        sh.slopes.push_back(b);
        sh.dslope.push_back(next - now);
      }
    }
    prepare_shift(sh, G.L.partners);
    long double delta = prior_delta(sh, 1.0);
    swap_move_state(sh);
    swap_set_label(to.label);
    _theta[0] = to.theta[0];
    _theta[1] = to.theta[1];
    _g[B.intercept] = to.gate ? 1 : 0;
    _beta[g0.coord0] = to.beta0;
    _b[g0.coord0] = t_new;
    for (unsigned c = 0; c < S.children.size(); ++c) new_terms += term_now(c);
    const double nested_new = S.mutant == 9 ? 0.0 : nested_term();
    delta += new_terms - old_terms + nested_new - nested_old;
    delta += B.log_prior[static_cast<std::size_t>(to.label - 1)] - B.log_prior[static_cast<std::size_t>(from.label - 1)];
    if (g0.gate != nullptr && to.gate != from.gate) delta += to.gate ? g0.log_on - g0.log_off : g0.log_off - g0.log_on;
    delta += log_slab(g0, &to.beta0) - log_slab(g0, &from.beta0);
    for (int b = 0; b < 2; ++b) {
      if (B.slope[b].source == nullptr || to.theta[b] == from.theta[b]) continue;
      const SlopePrior prior = slope_prior(B.slope[b], _chain);
      delta += static_cast<long double>(slope_log_prior(prior, to.theta[b])) - slope_log_prior(prior, from.theta[b]);
    }
    return delta;
  }

  // Restores the working state saved before swap_move.
  void swap_restore(Snapshot const &saved, SwapState const &from)
  {
    put(saved);
    _theta[0] = from.theta[0];
    _theta[1] = from.theta[1];
  }

  // One attempt of L2: propose a label of the other family, a gate, the
  // intercept and the entering slope from the target cell's law at the base
  // point (the leaving slope from its prior), E = E_0 + K_t t + K_b theta;
  // stage 1 on the evidences and the label Hastings factor, stage 2 on the
  // inner target with the complete proposal densities of both directions. A
  // dry run (probe) records the frozen-state identities of the attempt and
  // changes nothing.
  void label_swap(RNG *rng, bool dry = false)
  {
    struct SwapScope {
      bool &flag;
      explicit SwapScope(bool &f) : flag(f) { flag = true; }
      ~SwapScope() { flag = false; }
    } scope(_in_swap);
    BiasSpec const &B = bspec();
    int n_selection = 0, n_gaussian = 0;
    for (int l = 1; l <= B.levels; ++l) {
      if (!(B.log_prior[static_cast<std::size_t>(l - 1)] > negative_infinity)) continue;
      (B.gaussian(l) ? n_gaussian : n_selection) += 1;
    }
    if (n_selection == 0 || n_gaussian == 0) return;
    const int from = _label;
    const bool from_selection = !B.gaussian(from);
    SwapGeometry G;
    if (!swap_geometry(G)) return;
    if (!dry) bias_tally(bias_swap, 0);
    // the base point: the Gaussian labels' evidences, the selection labels'
    // screens and, leaving a selection label, its own Laplace evidence
    SwapBase base;
    std::vector<int> fits;
    if (from_selection) fits.push_back(from);
    if (!swap_base(G, fits, base)) return;
    Vec log_qG, log_qS;
    swap_label_proposals(base, log_qG, log_qS);
    Vec const &forward_q = from_selection ? log_qG : log_qS;
    int target = 0;
    {
      double u = rng->uniform();
      for (int l = 1; l <= B.levels; ++l) {
        const double q = forward_q[static_cast<std::size_t>(l - 1)];
        if (!(q > negative_infinity)) continue;
        target = l;
        const double p = std::exp(q);
        if (u < p) break;
        u -= p;
      }
    }
    if (target == 0) return;
    // entering a selection label: its Laplace evidence at the base point
    if (!from_selection && !swap_base(G, {target}, base)) return;
    const double stage1 = swap_stage1(base, log_qG, log_qS, from, target);
    if (!std::isfinite(stage1)) return;
    if (!dry) {
      swap_record_attempt(from_selection, stage1);
      if (!metropolis(rng, stage1)) {
        swap_trace(from, target, false, stage1, 0.0, false);
        return;
      }
      bias_tally(bias_swap, 1);
    }
    const SwapState current = swap_current(G);
    SwapState proposal;
    if (!swap_draw(G, base, current, target, rng, proposal)) return;
    Snapshot saved;
    take(saved);
    bool accepted = false;
    double full = negative_infinity;
    long double delta = 0;
    try {
      delta = swap_move(G, current, proposal);
      const double forward = swap_log_q(base, log_qG, log_qS, current, proposal);
      const double reverse = swap_log_q(base, log_qG, log_qS, proposal, current);
      full = static_cast<double>(delta + reverse - forward);
      if (dry) {
        swap_probe(G, base, log_qG, log_qS, current, proposal, stage1, full, delta, forward);
      } else {
        accepted = std::isfinite(full) && metropolis(rng, full - stage1);
      }
      if (accepted) {
        rebuild_h();
        _repr_delta += delta;
        bias_tally(bias_swap, 2);
        swap_record_accepted(from_selection, stage1);
      }
    } catch (...) {
      swap_restore(saved, current);
      throw;
    }
    if (!accepted) swap_restore(saved, current);
    if (!dry) swap_trace(from, target, true, stage1, full - stage1, accepted);
  }

  // Test hook (probe) of L2 at a proposal (the working state at it): the
  // geometry, E_0, every label's evidences and the label proposals rebuilt
  // from the proposal agree with those from the current state; the reverse
  // move's stage-1 ratio and full ratio are the inverses; and the accepted
  // probability flows of the two directions, from the target values and the
  // proposal densities computed at each end, are equal.
  void swap_probe(SwapGeometry const &G, SwapBase const &base, Vec const &log_qG, Vec const &log_qS,
                  SwapState const &current, SwapState const &proposal, double stage1, double full,
                  long double delta, double forward)
  {
    BiasSpec const &B = bspec();
    SwapGeometry G2;
    if (!swap_geometry(G2)) return;
    double k_diff = 0.0;
    for (std::size_t e = 0; e < G.Kt.size(); ++e) {
      k_diff = std::max(k_diff, std::fabs(G2.Kt[e] - G.Kt[e]));
      for (int b = 0; b < 2; ++b) k_diff = std::max(k_diff, std::fabs(G2.Kb[b][e] - G.Kb[b][e]));
    }
    std::vector<int> fits;
    if (!B.gaussian(current.label)) fits.push_back(current.label);
    if (!B.gaussian(proposal.label)) fits.push_back(proposal.label);
    SwapBase base2;
    if (!swap_base(G2, fits, base2)) return;
    Vec qG2, qS2;
    swap_label_proposals(base2, qG2, qS2);
    double e_diff = 0.0, q_diff = 0.0, e0_diff = 0.0;
    for (std::size_t e = 0; e < base.E0.size(); ++e) e0_diff = std::max(e0_diff, std::fabs(base2.E0[e] - base.E0[e]));
    for (int l = 1; l <= B.levels; ++l) {
      const std::size_t i = static_cast<std::size_t>(l - 1);
      auto diff = [](double a, double b) { return (a == b) ? 0.0 : std::fabs(a - b); };
      e_diff = std::max(e_diff, std::max(diff(base.log_e_off[i], base2.log_e_off[i]), diff(base.log_e_on[i], base2.log_e_on[i])));
      q_diff = std::max(q_diff, std::max(diff(log_qG[i], qG2[i]), diff(log_qS[i], qS2[i])));
    }
    const double stage1_reverse = swap_stage1(base2, qG2, qS2, proposal.label, current.label);
    const double forward2 = swap_log_q(base2, qG2, qS2, proposal, current);   // the reverse move's proposal
    const double reverse2 = swap_log_q(base2, qG2, qS2, current, proposal);
    const double full_reverse = static_cast<double>(-delta + reverse2 - forward2);
    // accepted flows: pi(x) q(x'|x) a(x, x') against pi(x') q(x|x') a(x', x)
    auto log_accept = [](double s1, double total) { return std::min(0.0, s1) + std::min(0.0, total - s1); };
    const double flow_forward = forward + log_accept(stage1, full);
    const double flow_reverse = static_cast<double>(delta) + forward2 + log_accept(stage1_reverse, full_reverse);
    probe_raw(probe_swap, 0.0, stage1, stage1 + stage1_reverse, full, full + full_reverse);
    probe_raw(probe_swap, 1.0, e0_diff, std::max(e_diff, q_diff), k_diff, flow_forward - flow_reverse);
  }

  // A raw probe record (type, five values).
  void probe_raw(int type, double a, double b, double c, double d, double e)
  {
    ProbeStore &store = probe_store();
    const std::lock_guard<std::mutex> lock(store.mutex);
    store.records.insert(store.records.end(), {static_cast<double>(type), a, b, c, d, e});
  }

  // ---- the evidence record of the bias label (LocationBlockBias.h,
  // BiasEvidence), over the updates after adaptation, flushed after each update
  BiasEvidence _bias_evidence;
  bool _bias_dirty = false;
  bool _swap_kept = false;            // L2 accepted a swap in this sweep
  bool _swap_kept_from_selection = false;
  double _swap_kept_stage1 = 0.0;

  // Counts a bias component's proposal (what 0), stage-1 pass (1) or
  // acceptance (2) in the per-fit counters and the evidence record.
  void bias_tally(int component, int what)
  {
    Counters &c = counters();
    count(what == 0 ? c.bias_proposed[component] : what == 1 ? c.bias_stage1[component] : c.bias_accepted[component]);
    const int i = bias_record_index(component);
    if (_adapting || i < 0) return;
    (what == 0 ? _bias_evidence.proposed[i] : what == 1 ? _bias_evidence.stage1[i] : _bias_evidence.accepted[i]) += 1;
    _bias_dirty = true;
  }

  // An L2 attempt's stage-1 log ratio by direction (0 leaving the selection
  // family, 1 leaving the Gaussian family).
  void swap_record_attempt(bool from_selection, double stage1)
  {
    if (_adapting) return;
    const int d = from_selection ? 0 : 1;
    _bias_evidence.attempts[d] += 1;
    _bias_evidence.attempt_sum[d] += stage1;
    _bias_evidence.attempt_max[d] = std::max(_bias_evidence.attempt_max[d], stage1);
    _bias_dirty = true;
  }

  void swap_record_accepted(bool from_selection, double stage1)
  {
    _swap_kept = true;
    _swap_kept_from_selection = from_selection;
    _swap_kept_stage1 = stage1;
  }

  // After the exact correction: the sweep's final label and, when the
  // correction kept the sweep, its accepted swap.
  void bias_after_correction(bool kept)
  {
    if (!bias() || _adapting) return;
    const int label = kept ? _label : _label_entry;
    if (_bias_evidence.label_sweeps.size() < static_cast<std::size_t>(bspec().levels)) {
      _bias_evidence.label_sweeps.resize(static_cast<std::size_t>(bspec().levels), 0.0);
    }
    _bias_evidence.label_sweeps[static_cast<std::size_t>(label - 1)] += 1;
    if (kept && _swap_kept) {
      const int d = _swap_kept_from_selection ? 0 : 1;
      _bias_evidence.kept[d] += 1;
      _bias_evidence.kept_sum[d] += _swap_kept_stage1;
      _bias_evidence.kept_max[d] = std::max(_bias_evidence.kept_max[d], _swap_kept_stage1);
    }
    _bias_dirty = true;
  }

  void flush_bias_evidence()
  {
    if (!_bias_dirty) return;
    BiasEvidenceStore &store = bias_evidence_store();
    const std::lock_guard<std::mutex> lock(store.mutex);
    store.records[_chain].add(_bias_evidence);
    _bias_evidence = BiasEvidence();
    _bias_dirty = false;
  }

  // ---- the L2 records (evidence record and per-attempt trace)
  // The per-attempt trace of L2 (test switch): chain, source and target
  // labels, stage 1 passed, stage-1 log ratio, stage-2 log ratio, accepted,
  // the first latent column entry (the random-effect SD) and the target
  // label's first six weights at the attempt.
  void swap_trace(int from, int target, bool passed, double stage1, double stage2, bool accepted)
  {
    if (!s().trace) return;
    Shared const &S = s();
    double sd = 0.0;
    for (unsigned c = 0; c < S.children.size() && sd == 0.0; ++c) {
      if (S.children[c].latent >= 0) sd = _d[S.children[c].row0];
    }
    double omega[6] = {0, 0, 0, 0, 0, 0};
    set_graph_label(target);
    for (unsigned c = 0; c < S.children.size(); ++c) {
      ChildSpec const &child = S.children[c];
      if (child.type == ChildType::normal) continue;
      const unsigned slot = child.type == ChildType::mv ? 3 : 4;
      Node const *weights = child.node->parents()[slot];
      for (unsigned long i = 0; i < weights->length() && i < 6; ++i) omega[i] = weights->value(_chain)[i];
      break;
    }
    set_graph_label(_label);
    ProbeStore &store = trace_store();
    const std::lock_guard<std::mutex> lock(store.mutex);
    store.records.insert(store.records.end(), {static_cast<double>(_chain), static_cast<double>(from),
      static_cast<double>(target), passed ? 1.0 : 0.0, stage1, stage2, accepted ? 1.0 : 0.0, sd,
      omega[0], omega[1], omega[2], omega[3], omega[4], omega[5]});
  }

  // Test hook (probe) of sweep condition 5 at L0's slot: every label's exact
  // score from the block's buffers (the working means shifted for the label's
  // slopes, the label's parameters written to the graph) against the graph's
  // own densities with the whole working state written to the graph (the
  // children's PDF_LIKELIHOOD, the nested node, the label prior); the graph is
  // restored to the update's entry and the working label afterwards.
  void label_probe()
  {
    Shared const &S = s();
    BiasSpec const &B = bspec();
    Vec working = _entry;
    for (unsigned j = 0; j < S.groups.size(); ++j) {
      GroupSpec const &g = S.groups[j];
      if (g.gate != nullptr) working[g.gate_pos] = _g[j] ? g.label : 3.0 - g.label;
      for (unsigned i = 0; i < g.d; ++i) working[g.slab_pos + i] = _beta[g.coord0 + i];
    }
    for (unsigned lat = 0; lat < S.latents.size(); ++lat) working[S.latent_pos[lat]] = _z[lat];
    for (int b = 0; b < 2; ++b) if (B.slope[b].source != nullptr) working[B.slope[b].pos] = _theta[b];
    double worst = 0.0, scale = 1.0;
    for (int l = 1; l <= B.levels; ++l) {
      const double prior = B.log_prior[static_cast<std::size_t>(l - 1)];
      if (!(prior > negative_infinity)) continue;
      set_graph_label(l);
      long double block = prior + nested_term();
      for (unsigned c = 0; c < S.children.size(); ++c) {
        ChildSpec const &child = S.children[c];
        _shift_mean.assign(_mean.begin() + child.row0, _mean.begin() + child.row0 + child.n);
        for (int b = 0; b < 2; ++b) {
          if (B.slope[b].source == nullptr) continue;
          const double change = physical_at(l, b) - physical_at(_label, b);
          for (unsigned i = 0; i < child.n; ++i) _shift_mean[i] += change * B.slope[b].column[child.row0 + i];
        }
        block += exact_with_mean(c, _shift_mean.data());
      }
      working[B.pos] = l;
      _view->setValue(working, _chain);
      _graph_label = l;
      long double graph = prior + nested_term();
      for (ChildSpec const &child : S.children) graph += child.node->logDensity(_chain, PDF_LIKELIHOOD);
      worst = std::max(worst, std::fabs(static_cast<double>(graph - block)));
      scale = std::max(scale, std::fabs(static_cast<double>(graph)));
    }
    _view->setValue(_entry, _chain);
    _graph_label = _label_entry;
    for (unsigned c = 0; c < S.children.size(); ++c) if (_step_children[c]) refresh_step(c);
    set_graph_label(_label);
    probe_raw(probe_label, 0.0, worst, scale, 0.0, 0.0);
  }

  // ---- the bias slots of the sweep (L1, S and L2 above this point)

  // Decided once per sweep before any component runs, independently of the
  // state (sweep conditions; the schedule of L2).
  bool _swap_this_sweep = false;
  // L2 is running: its fits and evaluations count as label_swap (LocationBlock::enter)
  bool _in_swap = false;

  void bias_sweep_start(RNG *rng)
  {
    _swap_kept = false;
    // the schedule of L2: attempted with a fixed probability per sweep
    _swap_this_sweep = s().swap_probability >= 1.0 || rng->uniform() < s().swap_probability;
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

  // One bias slot. With the probe (test hook) each move first records its
  // frozen-state identities at the state where it runs: in the reversed sweep
  // after the location components have moved the state (sweep condition 5).
  void bias_component(unsigned k, RNG *rng)
  {
    const bool probe = s().probe;
    switch (k) {
    case bias_dormant:
      dormant_draws(rng);
      break;
    case bias_label:
      if (probe) label_probe();
      label_update(rng);
      break;
    case bias_gaussian:
      if (probe) family_move(rng, true);
      family_move(rng);
      break;
    case bias_slope:
      if (probe) slope_move(rng, true);
      slope_move(rng);
      break;
    case bias_swap:
      if (probe) label_swap(rng, true);
      if (_swap_this_sweep) label_swap(rng);
      break;
    default:
      break;
    }
  }
