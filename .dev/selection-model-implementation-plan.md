# Selection-model implementation specification

Status: implemented and verified by independent source-cell references, actual JAGS checks and the standard/package checks. The 17 preexisting reference comparisons remain for maintainer review; scenario regeneration and extended certification were not rerun. Updated 2026-09-10. This is the authoritative implementation specification.

## 1. Scope and accepted interface

Use one model specification, separately from numerical controls. The supported selection constructors accept `selection = selection_model(...)` for automatically generated weightfunction priors. Explicit priors own their model:

~~~r
prior_bias <- prior_weightfunction(
  steps = .05,
  weights = wf_cumulative(),
  model = selection_model(
    estimate_random_effects = "integrate",
    other_random_effects    = "condition",
    known_sampling_variance = "condition",
    weight_rule             = "product",
    group                   = paper_id
  )
)
~~~

| Field | Accepted values | Default |
|---|---|---|
| estimate_random_effects | "condition", "integrate" | "integrate" |
| other_random_effects | "condition", "integrate" | "condition" |
| known_sampling_variance | "condition", "integrate" | "condition" |
| weight_rule | "product", "best" | "product" |
| group | Deferred data-column reference or NULL | NULL; resolve by section 5 |

Omitting the model or fields gives **integrate/condition/condition/product**, in the order above. Estimate random effects are integrated by default and can now be explicitly conditioned upon. All-integrated product selection is **integrate/integrate/integrate/product**. An absent source makes its two settings equivalent without creating a latent variable.

Keep these exact fields and the single constructor input `selection = selection_model()` in `bselmodel()`, `bselmodel.mv()`, `RoBMA()` and `RoBMA.mv()`. It configures generated default weightfunction priors only; every explicit prior retains its own complete `model`. Do not retain superseded mode arguments, aliases, object adapters, a second `condition_on` interface or a constructor override of explicit priors. `prior_weights` remains model prior odds; weight-height priors, mixture probabilities and outcome weights retain their separate meanings.

BayesTools owns the generic specification, prior validation/composition and declared estimate/other role resolution through `random_effects_level_roles()`. RoBMA owns data/group binding, covariance compilation, likelihood routing and summaries. Reuse compiled metadata and existing numerical kernels.

Current scope is three binary source axes and product/best weights. Estimated fractions, fraction-sharing controls, focal-result selection, combined reporting stages and finite attempted pools remain outside scope. Section 8 records a possible future extension only.

## 2. Statistical target and source roles

Conditional on population parameters, a representative hierarchy is

\[
Y_{ij}=\Gamma_j+u_{ij}+(B_jc_j)_i+\epsilon_{ij},\qquad
c_j\sim N(0,I),\quad \epsilon_j\sim N(0,D_j).
\]

`estimate_random_effects` controls the estimate-level true-effect source \(u\); `other_random_effects` controls the remaining true-effect sources, represented by \(\Gamma\). `known_sampling_variance` controls declared sampling context \(Bc\), not the numerical covariance matrix or its population parameters. Residual sampling variation \(D\) is **always integrated**. No axis alters variance or redistributes it between \(BB^\top\) and \(D\). The known baseline sampling covariance remains \(V=BB^\top+D\); selected outcomes generally do not have covariance \(V\).

### Omitted multivariate random structure

All six .mv constructors treat omitted or NULL `random` as no random-effects structure, like omitted `mods` leaves no moderator terms. The model is fixed-effect with respect to heterogeneity: no implicit `tau` source or heterogeneity mixture is introduced. Known sampling covariance `V` and any declared selection sampling context remain intact. Heterogeneity and `scale` formulas require an explicit `random` declaration. Do not supply `prior_heterogeneity` or `prior_heterogeneity_null` without that declaration, including explicit zero/NULL priors; omission represents the absent component.

Consequently, selection's estimate and other random-effect axes are inapplicable when .mv `random` is absent. They cannot manufacture an implicit source. This change is confined to the multivariate interfaces. Specialized univariate models keep their implicit estimate heterogeneity and, with `cluster`, their existing `tau`/`rho`/`I2` parameterization.

### Declared role resolution

For selected `brma.mv()` models and their wrappers, an estimate-level factor has a declared grouping factor whose levels are one-to-one with retained estimate rows. Allow zero or one such random-effect factor; more than one is a hard error naming the conflicting factors. With none, the estimate role is absent. All other random-effect terms use `other_random_effects`.

Resolve roles from BayesTools' declared grouping metadata after common input row selection. Keep each declared family intact: do not split a factor, slope, coefficient basis or covariance family according to realized loadings or numerical covariance. Non-diagonal known group `R` does not change the role of a one-to-one estimate factor. Preserve its full covariance and dependencies. Only diagonal cases use the scalar estimate-variance fast path; general covariance uses the existing compiled path. Sampling a latent node, numerical rank and a loading affecting one row do not determine statistical role.

In unclustered univariate models, ordinary heterogeneity is the estimate source. In the specialized `brma(..., cluster = ...)` interface, within-cluster heterogeneity is the estimate source and the shared cluster effect is the other source. Preserve specialized total `tau`/`tau2`, allocation `rho` and within-/between-cluster `I2` summaries.

Plain `vi`/`sei` and undeclared diagonal `V` have no sampling context. Explicit or certified sampling factors, including singleton columns, remain sampling context and never become estimate random effects. Random-effect roles do not define publication groups.

| estimate_random_effects | other_random_effects | known_sampling_variance | Retained \(H\) | Integrated before normalization |
|---|---|---|---|---|
| integrate | condition | condition | \(\Gamma,c\) | \(u,\epsilon\) |
| condition | condition | condition | \(u,\Gamma,c\) | \(\epsilon\) |
| integrate | condition | integrate | \(\Gamma\) | \(u,c,\epsilon\) |
| condition | condition | integrate | \(u,\Gamma\) | \(c,\epsilon\) |
| integrate | integrate | condition | \(c\) | \(u,\Gamma,\epsilon\) |
| condition | integrate | condition | \(u,c\) | \(\Gamma,\epsilon\) |
| integrate | integrate | integrate | None | \(u,\Gamma,c,\epsilon\) |
| condition | integrate | integrate | \(u\) | \(\Gamma,c,\epsilon\) |

“Condition” means normalize at each retained latent value, then average over its original mixing law. The latent source remains unknown and estimated, and its posterior can update with data. Selection can reweight the population mixing law of integrated sources. Constant weights or absent sources give coincident cases. Correlated retained/integrated sources require their actual conditional Gaussian law.

These are different sampling experiments, not numerical approximations to one model. One-to-one estimate grouping does not imply independence.

### One likelihood for all eight cells

Let \(H\) contain retained sources and \(U\) those integrated before normalization, always including residual sampling variation. Define

\[
g_H(\mathbf y\mid h)=\int p_0(u\mid h)f_0(\mathbf y\mid h,u)\,du,
\qquad
A_H(h)=\int g_H(\mathbf z\mid h)W(\mathbf z)\,d\mathbf z.
\]

The likelihood is

\[
L_H(\mathbf y)=\int p_0(h)
\frac{g_H(\mathbf y\mid h)W(\mathbf y)}{A_H(h)}\,dh.
\]

With no retained source, normalize the marginalized Gaussian vector once. Otherwise keep \(A_H(h)\) inside the integral over \(h\). At fixed \(h\), the selected integrated-source law is proportional to \(p_0(u\mid h)E_0[W(Y)\mid h,u]\); retained sources keep their original mixing law \(p_0(h)\).

Numerical augmentation may sample \(U\), but must use \(p_0(h)p_0(u\mid h)f_0(y\mid h,u)W(y)/A_H(h)\), not a denominator conditional on \(u\). Backend sampling versus analytic marginalization must not move this boundary.

Product normalizers factor by row only when rows are conditionally independent given \(H\). Integrating a shared source can require a joint normalizer even when another source is retained. Fixed role choices add no model parameters.

### Preserve existing scenario targets

Existing `cond` labels mean **integrate/condition/condition/product**; `marg` means **integrate/integrate/integrate/product**. Neither varies the estimate axis. Preserve budgets, seeds, source declarations and deferred work when migrating syntax. Do not relabel results as equivalent when source-role semantics or the sampling law differ, or treat old caches/snapshots as refreshed evidence.

For the specialized clustered `cond` model,

\[
\Gamma_j\sim N(\mu,\tau_b^2),\qquad
Y_{ij}\mid\Gamma_j\sim SN(\Gamma_j,s_{ij}^2+\tau_w^2,\omega).
\]

Its cluster integral is outside the product of row normalizers; within-study variance is inside each row kernel. Its `marg` counterpart also integrates the cluster source before normalization. With no other or sampling context, the default stays \(SN(\mu,s_i^2+\tau^2,\omega)\). Explicit estimate conditioning is a newly available, different target.

## 3. Vector weights and the observation design

Conditioning specifies the normalization boundary; weight_rule specifies \(W\). Support both independently, including when all applicable sources are conditioned upon:

- **Product:** \(W(\mathbf y)=\prod_i w_i(y_i)\). One relative weight factor per estimate. With positive weights it does not require all results to be significant. Independent reporting gates, normalized at the chosen conditioning level for the specified retained vector, give this rule. With independent conditional kernels it also admits repeated candidate draws to a prescribed reported count.
- **Best:** \(W(\mathbf y)=w(\min_i p_i(y_i))\). Use the bin of the smallest p-value, not the largest weight. The one-cutoff, fully integrated version is the paper's relaxed, report-all rule: a favorable result can qualify the vector, whose measured outcomes are then all reported. Other cells and the multiple-cutoff version are extensions of that specification.

There are eight conditioning cells for each rule, with coincident special cases. Do not automatically expand mixture model odds over these choices. Product remains the default; best does not reproduce the conditional per-row reporting model.

For two estimates with favorable weight 1 and unfavorable weight .2:

| Pattern | Product relative weight | Best relative weight |
|---|---|---|
| Both unfavorable | .04 | .2 |
| One favorable | .2 | 1 |
| Both favorable | 1 | 1 |

These are unnormalized weights, not pattern probabilities. Even with independent residual rows after conditioning, best generally couples the selected residual outcomes. It can produce weak marginal z distortion when nearly every large study contains a favorable result. The paper's strict rule differs from product: its one penalty for any unfavorable result is not a penalty raised to the number of unfavorable results.

Use original sampling SEs, \(\sqrt{\operatorname{diag}(V)}\), for selection thresholds; never substitute conditional SDs or heterogeneity-inflated SDs. Preserve existing weightfunction side, reference-bin, and exact-zero contracts. Support unequal SEs, moderators, positive/negative one-sided and two-sided geometry, and valid nonmonotone or above-one relative weights.

For cutpoints \(0=a_0<\cdots<a_K=1\), define \(G_H(a\mid h)=P_0(\text{all }p_i>a\mid h)\) under \(g_H\). Then

\[
A_H^{\mathrm{best}}(h)=\sum_{k=1}^K\omega_k
\{G_H(a_{k-1}\mid h)-G_H(a_k\mid h)\}.
\]

Evaluate the corresponding conditional Gaussian rectangle/tail probabilities. Independent covariance subblocks within one publication group factor probabilities inside \(G_H\); they do not create separate publication events. Use continuous-law endpoint identities only where valid; handle supported singular boundaries explicitly. A zero-probability selection event cannot be normalized.

Condition on the specified reported/measured vector and its design, including count and SEs. Do not infer an attempted count, an absolute publication probability, or a finite missing pool from relative weights. Conditional reporting can use an effectively unbounded candidate distribution while keeping context fixed; this neither implies fresh datasets nor establishes candidate independence. A finite-pool/count likelihood is a separate model outside this change.

## 4. Known sampling context: resolve once, then evaluate

For sampling-conditioned targets, the split \(V=S+D\), \(S=BB^\top\), is part of the statistical model. \(V\) alone does not uniquely identify which sampling variation persists during hypothetical reporting/search. Equivalent factor rotations preserving the declared retained law and residual law are computational changes; moving variance between \(S\) and \(D\) changes the conditional likelihood.

Resolve the structure after canonical input row selection and before backend choice, using existing known-V parsing and factor helpers:

| Input/provenance | Required resolution |
|---|---|
| Plain vi/sei or undeclared diagonal V | No sampling context; all sampling variance in \(D\) |
| Explicit known_v_factor(diagonal = d, loading = B) | Honor its declared factor and residual, even if reconstructed V happens to be diagonal |
| vcalc2() construction with factors certified by the current implementation | Honor its existing declared factor and residual |
| Correlated matrix without a valid declared/certified structure, including an uncertified vcalc2() construction | Preserve the current versioned generic decomposition; warn as specified below |
| Fit created by this implementation | Restore its stored target and decomposition; do not resolve using new defaults |

For a generic nonsingular correlated block, the existing default is

\[
D=\alpha\,\operatorname{diag}(\operatorname{diag}(V)),\qquad
\alpha=\min\{0.10,\;0.99\lambda_{\min}(\operatorname{cor}(V))\}.
\]

Here \(D\) is a diagonal matrix. Preserve existing block planning, diagonal/singleton special cases, supported exact rank-one factors with zero \(D\), and explicit failure for unsupported near-singular cases. This generic split is a versioned modelling convention chosen to preserve the existing conditional target; it is not a uniquely inferred physical overlap structure.

Persist the residual, retained covariance/factor structure, source identities, row alignment, provenance, and decomposition policy/version. A structured representation suffices; do not materialize a dense retained covariance solely for storage. Record generic requested/effective fractions as reconstruction metadata without introducing a new estimated parameter.

Numerical backends may exploit blocks, sparsity, or factor rotations while preserving the declared statistical structure. They may not change \(D\), drop covariance, or move a source across the normalizer to gain speed. Fitting and every post-fit consumer must use the same resolved structure. Within-fit deletion and partial-vector operations keep the original selection experiment; do not rerun the generic split on the remaining rows. A genuinely new fit/design resolves its own supplied inputs under the recorded policy.

Validate row identity, metadata version/hash, factor reconstruction, positive semidefiniteness, residual nonnegativity, and supported conditional-kernel boundaries. Invalid or stale declarations are errors, not permission to fall back silently. Zero residual sampling variance can be valid when integrated estimate-level variance remains. Do not invent residual noise or repair covariance to avoid a boundary.

When known_sampling_variance="integrate", this axis depends only on the full \(V\) law; no statistical contextual split is required. Gaussian brma.mv is likewise invariant to valid numerical decompositions. Keep known_v_residual_fraction out of the new public modelling interface and remove it from unreleased Gaussian examples/API.

### Required fallback warning

Emit a standard R warning, once per fit, when **at least one active weightfunction branch uses known_sampling_variance="condition" and correlated sampling input requires the generic undeclared decomposition**. Continue fitting with that documented decomposition; do not switch cells or reject an otherwise supported matrix.

Use this message for a plain matrix:

~~~text
'known_sampling_variance = "condition"' is using the default decomposition of a correlated 'V' because no declared sampling factors are available. The conditional likelihood depends on this decomposition. Use a supported 'vcalc2()' construction that supplies certified factors, or 'known_v_factor()' to specify them explicitly.
~~~

If the input already came from an unsupported/uncertified vcalc2() construction, replace the first sentence with:

~~~text
The supplied 'vcalc2()' construction does not provide certified sampling factors; 'known_sampling_variance = "condition"' is using the default matrix decomposition.
~~~

Retain the remaining explanation and remedies. Use call. = FALSE. Deduplicate across blocks and mixture branches, and do not reissue it during MCMC or ordinary post-fit operations. Store/report the fallback provenance independently of warning display.

Do not emit this warning for an undeclared diagonal input with no sampling context, a valid explicit factor, an already certified construction, Gaussian/PET/PEESE-only fits, or fits whose active selection branches all integrate the sampling context. Here an active branch is one with nonzero prior model probability whose outcome weights can affect the likelihood; a structurally constant weightfunction alone does not require this warning. A no-selection mixture component does not suppress the warning needed by an active conditional selection component. Invalid/stale metadata remains an error.

Help must explain that merely wrapping an arbitrary matrix construction in vcalc2() is insufficient, and explicitly choosing different structural factors can change the conditional target. Recommending vcalc2() is not a promise that every call supplies factors or reproduces the old generic split.

### Improve metadata without changing the default model

Current factor certification covers only supported constructions. In particular, the economics example's vcalc2(vi = sez^2, cluster = newid, obs = obs, rho = .5, ...) lacks type and currently uses the generic split. Improve basic cluster/row metadata for scalar-rho/no-type calls independently of factor certification. Preserve original clusters and subgroup identities, unequal SEs, and duplicate observation identities. Reuse evaluated inputs instead of evaluating user expressions twice with different side effects.

Do not silently give a previously uncertified call a new structural split in this refactor. Retain its existing conditional target and issue the applicable warning. Adoption of a new factor representation can be explicit through known_v_factor(); broader automatic certification with a changed default target requires a separate model/default migration.

## 5. Publication groups, binding, and dependencies

Publication groups, random-effect groups, and independent covariance blocks are distinct. Resolve the publication partition in this order:

1. Explicit selection_model(group = ...) column reference.
2. Validated vcalc2() publication-cluster metadata.
3. Explicit cluster from a supported specialized/univariate constructor.
4. Singleton units only for a genuinely independent, unambiguous univariate input.

Otherwise require an explicit group. Do not infer publications from covariance connected components, and do not replace original cluster by cluster × subgroup. Explicit grouping can combine cohorts into their publication unit. Group metadata and covariance-source metadata have separate validation and provenance.

### Deferred column binding

Capture group = paper_id unevaluated when the prior/model is constructed. It must work before paper_id or the fitting data exists. Store a serializable unresolved column reference, not a vector, data frame, or arbitrary captured environment. Bind it to data when bselmodel.mv, RoBMA.mv, or another supported fitting interface is called.

Initially support a bare column name, including a backticked non-syntactic name, and a single character column name for programmatic construction. Normalize both to the same reference; omission/NULL requests automatic resolution. Bare symbols denote data columns, not same-named global vectors. Wrappers must forward the captured reference or construct the call explicitly, without capturing their local formal-argument name. Reuse existing NSE/input-binding helpers; do not add a general expression interpreter. Derived groups should be columns in data.

Apply the same subset/NA row map as outcomes and V. Validate missing columns, length, and missing identifiers among retained rows; do not silently fall back to another environment or singleton groups. Preserve requested reference separately from resolved values/provenance through mixtures, serialization, fitted reconstruction, and newdata. Reusing a prior with another data frame binds afresh. Preserve/revalidate metadata under row filtering and permutations.

### Normalization and initial support boundary

Plan dependencies **after conditioning on retained sources**. Shared retained effects may connect publication groups while their selected kernels factor conditional on those effects; integrating them afterwards can still couple the marginal likelihood. Integrated sources connecting publication groups generally require joint normalization. Shared population parameters alone do not break conditional independence.

For jointly handled selection events, the vector weight is the product of the declared publication-group weights; do not reinterpret independent Gaussian subblocks as new publication decisions. Inspect both known sampling covariance and random-effect sources when building the dependency plan.

Deliver all eight cells for product and best on supported Gaussian sources confined to publication groups. Preserve the existing broader supported product paths. Additional crossed-source combinations require a valid dependency/normalization implementation; if unavailable, reject explicitly, naming the source and affected groups. Do not silently select the conditional cell, change weights, zero cross-group covariance, or merge publication units.

For the initial ensemble implementation, require a common resolved conditioning cell and publication partition across active weightfunction branches, allowing the existing different bin/weight priors and supported weight rules. Preserve each child's requested specification even when a mixed-cell fit is unsupported; never copy the first child's settings over other children. Extending mixed-cell fitting is separate from sharing fixed settings.

## 6. Compatibility, inference, and user-visible output

### Target preservation

The unreleased `exact`/`approximate` selection-target API is replaced outright by the resolved `selection_model()` specification. New fits default to **integrate/condition/condition/product**; all-integrated product selection is **integrate/integrate/integrate/product**. Preserve prior targets only when declared roles and sampling decomposition match. Newly fitted objects record explicit modes, source roles, selection groups, applicability, and decomposition policy/version. Their log likelihood, evidence, prediction, and other post-fit operations preserve that recorded target and decomposition.

Historical prior and fitted objects from superseded iterations of this unreleased API are not migrated. Recreate their priors with `selection_model()` and refit when the current resolved specification or execution metadata are unavailable; do not add aliases, schema adapters, inferred defaults, or compatibility likelihoods for those objects. This full replacement does not authorize breaking released interfaces: preserve released contracts under repository policy, with a separate maintainer decision for any released incompatibility.

Gaussian no-selection, PET, and PEESE branches retain their existing targets. Constant outcome weights must recover the same Gaussian observed-data law in every cell. Binary settings add no parameters to no-selection or selected branches.

### One model specification across inference

Route JAGS/native fitting, joint log likelihood, marginal likelihood/evidence, density/CDF/RNG methods, latent reconstruction, prediction, residuals, z plots, and diagnostics through the same resolved model. Numerical budgets remain in evaluator controls such as selection_control. Use analytic rows only when justified, then existing quadrature or diagnosed multivariate integration. Integration failure must not trigger a different cell or weight rule.

For selected new-group predictions, draw retained context from its original mixing law and then draw the selected integrated sources/outcomes conditional on that context. For existing-group predictions use the appropriate inferred context. Label selected versus unselected-population and new-group versus existing-group predictions separately.

Given complete outcomes and fixed retained context, weights cancel from the posterior of integrated Gaussian sources (for an observation with positive selected density). They do not cancel the retained-context-dependent normalizer from inference about retained sources. Implement latent summaries accordingly.

Group predictive scores must use the full specified selection event. Partial-vector prediction and row deletion integrate missing outcomes under that event; they do not define a new event from the observed subset. Do not present sums of conditional row scores as interchangeable with joint group scores. Simulation recovery and held-out selected-data prediction serve different purposes; predictive fit to the published sample does not identify the unseen population.

Reframe discrepancies between the former exact and approximate targets as **conditioning/weight-rule sensitivity**, separately from integration accuracy. They are not numerical error bounds or an automatic rule for choosing the correct reporting mechanism. Avoid creating a new diagnostic API name unless needed by existing conventions.

Known covariance still protects against treating overlapping data as independent. Before selection, the equicorrelated example has

\[
\operatorname{Var}(\bar Y\mid\Gamma)=rs^2+\{(1-r)s^2+\tau_w^2\}/n.
\]

Do not promise this Gaussian variance formula or an unchanged scalar study weight after selection. For the retained-common-error construction \(L(\Gamma)=\int\phi(t;\Gamma,\sigma_C^2)K(y\mid t)\,dt\), with a normalized kernel \(K\) and other parameters fixed,

\[
-\ell''(\Gamma)=\frac{1}{\sigma_C^2}
-\frac{\operatorname{Var}(T\mid y,\Gamma)}{\sigma_C^4}
\leq\frac{1}{\sigma_C^2}.
\]

This is a scoped likelihood-curvature bound, not a universal posterior-variance floor or a guarantee for sampling-integrated targets. Keep the existing numerical information check tied to those assumptions.

### Printing and documentation

Print the modes and their meaning, vector rule, group reference/provenance, and sampling-structure provenance. Before binding, automatic sources/groups are unresolved; afterwards print the resolved result. Report absent contexts as absent with coincident choices. Preserve requested settings separately from resolved applicability.

Use human labels “Estimate random effects”, “Other random effects” and “Known sampling context (from V)”, with source details appropriate to the constructor, and “Product of estimate weights” / “Weight of the best p-value”. Explain that conditioned sources remain unknown and estimated. State which population mixing laws can change, without claiming their posteriors cannot update. Do not label product “all estimates significant”, best “maximum weight”, or statistical cells “exact”/“approximate”.

## 7. Prior composition now

Preserve the model specification through raw lists, prior_mixture, bias wrappers, branch maps, serialization, reconstruction, and fit-time binding. Repeated fixed choices add no prior density, conditioning parameter, monitor, sampler coordinate, or sharing argument. Reusing the same R object is a convenience and does not establish ownership of estimated parameters.

Keep existing weight-height priors, mixture indicators, components/is_null, and model odds unchanged. Do not implement a new bias prior_spike_and_slab workflow incidentally; retain existing supported selection/no-selection mixture behavior. Unsupported mixtures must fail without overwriting child metadata.

## 8. Future extension boundary only

Keep estimate random effects, other random effects and sampling context separate and preserve declared covariance-source identity independently of numerical factors. A future retained fraction for each role could split its source between retained and integrated covariance. Its endpoints would be condition = 1 and integrate = 0. Current public arguments accept only the two strings; compile constants directly, without point-prior parameters or unused ownership machinery.

If estimated fractions are added later, the intended default is one owner per compatible role, source/group meaning, fraction scale, and prior within a containing bias mixture, with explicit opt-out. Distinct roles are not shared merely because their prior distributions match. Repeated references contribute their prior once; nested scope, branch odds, and serialization must remain stable. All inference consumers must use the same owner references. No selection makes a fraction inactive, not zero; zero means full integration.

Estimability requires a separate assessment. For aligned common random and sampling shifts, the observed likelihood can identify only a combined retained variance such as \(h_\Gamma\tau_b^2+h_c v_c\), not the two fractions separately. Proper priors do not remove exact observational equivalence. No Uniform default, fraction API, or general sharing registry is required now.

## 9. Implementation sequence and completion checks

1. Read current repository AGENTS.md and only the relevant detailed guides. Audit existing target/source metadata and released compatibility. Add the generic prior/model object, deferred group capture, validation, printing, and propagation through composition.
2. Resolve/persist source roles, publication groups, and conditional sampling structure before backend selection. Preserve the specified generic conditional decomposition, add the fallback warning, and improve vcalc2() grouping metadata independently of factor certification.
3. Wire integrate/condition/condition/product defaults and explicit integrate/integrate/integrate/product, checking preserved targets at the same declared roles and decomposition. Implement estimate conditioning and the other supported cells and best-weight normalizers against independent small-dimensional references. Publish a precise capability matrix; do not advertise a fitting target until its required inference operations work.
4. Route all inference and post-fit consumers through the resolved specification. Preserve the recorded targets of newly fitted objects, selection events for partial vectors, and existing nonselection branches. Reclassify target differences separately from numerical diagnostics.
5. Update help, examples, and release metadata according to repository policy; remove superseded unreleased arguments and terminology. Run focused checks, then required standard checks. Put expensive recovery/integration certification in the certification profile.

Completion requires these focused checks; use independent references rather than tests that repeat the implementation:

| Area | Required evidence |
|---|---|
| Statistical targets | Normalization for every supported cell/rule; constant weights recover Gaussian inference; independent low-dimensional references for conditional/integrated numerators and denominators |
| Defaults and full replacement | Omitted model/modes equal explicit integrate/condition/condition/product; constructor selection applies only to generated priors; explicit priors retain all three modes/group; all-integrated settings reproduce former marginal references only when declared roles and decomposition match; newly fitted objects retain their resolved targets; superseded unreleased objects are not migrated |
| Source roles | Default ordinary marginal SN unchanged; independent estimate-conditioning references; clustered vi/sei retains tau/rho/I2 and its default model; zero/one one-to-one factor supported and multiple factors fail; non-diagonal known R retains full estimate-source covariance; no mixed-family splitting; absent roles collapse dimensions |
| Decomposition | Integrated targets invariant to equivalent representations of V; conditional targets invariant to refactorization of fixed retained/residual laws; plain diagonal differs from an intentional explicit contextual declaration when appropriate |
| Source resolution | Explicit/certified precedence; generic .10/capped block convention and supported singular cases; stale metadata errors; scalar-rho/no-type grouping improvements leave the old conditional split unchanged |
| Warning | Complete messages for plain matrix and uncertified vcalc2(); once per fit across blocks/branches; fitting continues with stored generic provenance; no warning in the excluded cases from section 4 |
| Weight rules | Product conditional row factorization where justified; best-induced dependence and the toy weights above; multiple bins, unequal SEs, moderators, both sides, exact zeros, valid nonmonotone/above-one weights, and subblocks within one publication |
| Group binding | Prior construction before data exists; bare/string/backticked names; wrapper forwarding; no global-vector capture; missing columns/IDs; subset/NA alignment; duplicates; permutation; explicit precedence; serialization, prior reuse, and newdata |
| Dependencies and mixtures | Cross-publication sampling/random sources checked; supported product paths preserved; unsupported combinations identified; child specifications/odds preserved; repeated fixed settings add no parameters or priors |
| Inference consistency | Fit/log density/evidence/predictive/latent consumers target the same normalized law; selected versus unselected predictions distinguished; partial vectors retain the full selection event |
| Numerical validation | Diagnosed integration convergence against independent references; no target-changing fallback, silent covariance repair, or invented residual noise; selected-information checks keep their stated scope |

### Source map and reference material

The following paths were inspected during planning; verify the current tree before editing. This map is a starting point, not a requirement to modify every file.

| Repository | Starting points |
|---|---|
| C:/R-Packages/BayesTools | R/priors-weightfunction.R, R/priors.R, R/selection-kernels.R, existing prior/NSE validators and random-effect source metadata |
| C:/R-Packages/RoBMA | R/bselmodel.R, R/RoBMA.R, R/bselmodel.mv.R, R/RoBMA.mv.R, R/fit.R, prior/input constructors |
| C:/R-Packages/RoBMA | R/input-data-mv.R (.known_v_decompose_block), R/known-v-representation.R, R/vcalc2.R, R/selection-likelihood.R, R/selection-mapping.R, src/selnorm/ |
| C:/R-Packages/RoBMA | R/marglik.R, R/iwmde-likelihood.R, density/prediction/RNG/zplot/residual methods, R/selection-sensitivity-diagnostics.R |

Paper context: van Aert, Riley, and Jackson, *Correcting for publication bias in multivariate and multilevel meta-analysis: A multivariate step function selection model approach*, Figure 1 and printed pp. 8–9. Attribute its relaxed rule accurately and label additional conditioning/multicutoff constructions as extensions.

This specification covers RoBMA and BayesTools implementation, including their tests, maintainer scenarios, documentation, and scenario caches. It does not request refitting the separate economics analysis, redeploying the demonstration site, or performing Git operations.
