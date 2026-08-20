#' Random-effect prior constructors
#'
#' @description
#' Constructors for specifying priors and policies for the random effects in
#' [brma.mv()]. These functions are re-exported from BayesTools so a partial or
#' complete random-effect prior can be written with the RoBMA namespace alone.
#'
#' `prior_random()` collects the specification supplied through
#' `prior_heterogeneity`. `random_block()` overrides priors or policies for one
#' named random-effect block.
#' `random_covariance()` supplies correlation priors, with `prior_lkj()` used
#' for an unstructured correlation matrix. `random_variance_allocation()`
#' distributes an aggregate random-effect standard deviation across blocks or
#' across the SD components of one heterogeneous block. Its required `name` is
#' a stable internal identifier, while `display_name` and `component_names`
#' control public semantic labels. `allocation_ref()` connects a child
#' allocation to a component of an earlier named allocation.
#' `random_monitor()` controls saved random-effect quantities, and
#' `random_new_levels()` controls prediction for previously unseen grouping
#' levels. `random_sd_source()` is the advanced adapter for an existing scalar
#' or row-shaped JAGS node that supplies a random-effect SD without owning a new
#' prior for that node.
#'
#' @details
#' ## Formula structure is separate from prior structure
#'
#' The `random` formula owns the grouping variables, coefficient or index
#' design, and covariance family. These prior constructors configure the
#' already parsed blocks; they do not reinterpret the formula. The `id()`,
#' `diag()`, and `us()` / `un()` families use a coefficient formula: `1`, `0`,
#' and `-1` control the random intercept, and continuous slopes, factor slopes,
#' and interactions are supported. Plain `|` syntax defaults to US and `||` to
#' DIAG.
#'
#' Factor coding in a random-coefficient block is owned by resolved contrast
#' metadata, not by the presence of an intercept. Existing fixed-factor
#' contrasts are reused when available; `random_block(contrasts = ...)` supplies
#' an explicit block-specific override. Thus `us(0 + group | study)` combined
#' with `contrasts = c(group = "independent")` gives one correlated coefficient
#' per factor level and no random intercept. Omitting the intercept alone does
#' not force a level-indicator basis.
#'
#' The `cs()` / `hcs()`, `ar1()` / `ar()` / `har()`, and `car()` structures
#' instead own an index basis. Their blocks reject intercept controls and
#' `random_block(contrasts = ...)`. HCS estimates one common pairwise
#' correlation; US estimates an unrestricted correlation matrix among the
#' selected coefficient columns. See [random_effect_formula_tags()] for the
#' full formula and index-type contract.
#'
#' ## RoBMA completion of partial specifications
#'
#' A named `random_block()` must match the block name produced by the parsed
#' formula. Use `only_priors = TRUE` when constructing a [brma.mv()] object and
#' then `print_prior(object, parameter = "random")` to inspect the resolved
#' names and completed specification. An ordinary positive prior supplied to
#' `prior_heterogeneity` changes RoBMA's base heterogeneity prior while retaining
#' its automatic allocation and covariance rules. A `prior_random()` object with
#' no SD or allocation specification is completed in the same way: RoBMA adds
#' its UISD-scaled base SD and any required Dirichlet allocation while retaining
#' the supplied block contrasts, covariance priors, and policies. If any global
#' or block SD, SD source, term-specific SD override, covariance-owned SD, or
#' allocation is supplied, the user owns the complete scale architecture and
#' RoBMA does not merge in scale defaults.
#'
#' Consequently, covariance-only, contrast-only, monitoring-only,
#' prediction-only, and parameterization-only objects are partial and
#' self-finishing in
#' RoBMA. For example,
#' `prior_random(study = random_block(contrasts = c(group = "independent")))`
#' needs no explicit SD, allocation, or LKJ prior for a US block. In contrast,
#' adding `sd = ...` makes the object complete; every other required scale or
#' allocation must then be supplied explicitly.
#'
#' ## Standard-deviation and correlation priors
#'
#' Standard-deviation priors must have nonnegative support. Scalar-correlation
#' structures use `cor` in `random_covariance()`; `cor_scale = "fisher_z"`
#' places the prior on Fisher's z, whereas `cor_scale = "cor"` places it on the
#' raw correlation. When omitted, CS/HCS receive a uniform prior over
#' `(-1 / (K - 1), 1)`, AR1/HAR over `(-1, 1)`, and CAR over `(0, 1)`.
#' Unstructured `us()` blocks default to `prior_lkj(eta = 1)`. Values of `eta`
#' above one shrink correlations toward zero; `eta = 1` is uniform over
#' correlation matrices.
#'
#' The covariance family remains formula-owned. `random_covariance(structure =
#' ...)` is mainly useful in direct BayesTools workflows; when used with
#' `brma.mv()`, an explicit structure must match the parsed formula family.
#' Use `prior_lkj()` as `cor` for an unstructured matrix prior and an ordinary
#' scalar prior as `cor` for a scalar-correlation structure.
#'
#' ## Variance allocation
#'
#' A block allocation with `scale = "total_variance"` exposes `sd_total` and
#' `var_total` and splits that total variance across multiple random-effect
#' blocks:
#' \deqn{\sigma_j = \sigma_{\mathrm{total}}\sqrt{w_j}, \qquad
#'       \sum_j \sigma_j^2 = \sigma_{\mathrm{total}}^2.}
#' An SD-component allocation with
#' `scale = "mean_variance"` instead exposes `sd_common` and `var_common` and
#' gives each component the common-SD scale at equal weights:
#' \deqn{\sigma_j = \sigma_{\mathrm{common}}\sqrt{K w_j}, \qquad
#'       K^{-1}\sum_j \sigma_j^2 = \sigma_{\mathrm{common}}^2.}
#' This is useful for heterogeneous level-specific SDs. Root allocations own an
#' SD prior or external SD source. Child allocations inherit one named parent
#' component through `allocation_ref()` and therefore must not supply another
#' SD source. `inclusion` optionally places independent Bernoulli gates on
#' named block-allocation components.
#'
#' ## Public parameter names
#'
#' Random-effect quantities are exposed through the semantic catalog view of
#' BayesTools' fitted parameter map. Concrete posterior coordinates are a linked
#' view of the same map and are not additional public names. BayesTools
#' canonical names use `(formula) owner: quantity(arguments)`. RoBMA prints and
#' accepts simplified aliases: a sole random intercept is `sd` for a bare block
#' and `study: sd` for a named block; a non-intercept coefficient remains
#' explicit, for example `study: sd(x)`. An owner-free `sd` is accepted only
#' when it resolves uniquely. Public correlations use `cor`; compact scalar
#' `rho` and LKJ construction coordinates remain internal backend dependencies.
#'
#' A bare formula or unnamed one-entry list omits a redundant component owner.
#' An explicitly named one-entry list retains its owner. For multiple
#' components, public component names come from the formula list; missing names
#' are generated as `component 1`, `component 2`, and so on.
#' Total-variance allocations expose `sd_total`, `var_total`, and
#' `var_prop(...)`. Mean-variance allocations expose `sd_common`, `var_common`,
#' `var_ratio(...)`, and `sd_ratio(...)`.
#'
#' ## Monitoring, prediction, and parameterization
#'
#' Keep the default `random_monitor(latent = TRUE)` when bridge sampling or
#' conditional prediction is required. For new grouping levels,
#' `random_new_levels("error")` is the safe default; `"zero"` removes their
#' random contribution and `"sample"` draws from the fitted random-effect
#' distribution.
#'
#' RoBMA creates external SD sources automatically for its location-scale
#' random-effect models. Use `random_sd_source()` directly only when a custom
#' complete specification intentionally references a JAGS node created elsewhere
#' in the model being fitted.
#'
#' `parameterization` changes only the sampled backend representation.
#' `"noncentered"` samples standardized latent effects, `"centered"` samples
#' realized group coefficients, and `"auto"` asks BayesTools to choose from the
#' resolved design. It does not change any SD, allocation, or correlation
#' prior. Unsupported centered configurations are rejected; `"auto"` can fall
#' back to noncentering.
#'
#' @param ... named `random_block()` overrides and, optionally,
#'   `random_variance_allocation()` objects. Override names must match resolved
#'   block names.
#' @param sd an SD prior with nonnegative support. In `prior_random()`, it is the
#'   top-level default inherited by blocks. In `random_block()` and
#'   `random_covariance()`, it is a more specific block or covariance SD prior.
#'   In a root `random_variance_allocation()`, it is the aggregate-SD prior.
#' @param covariance a `random_covariance()` object. A top-level value is
#'   inherited by blocks without a block-specific covariance override.
#' @param cor an optional correlation prior. Use `prior_lkj()` for an
#'   unstructured correlation matrix and an ordinary scalar prior for CS, HCS,
#'   AR1, HAR, or CAR. It is a convenience alternative to
#'   `covariance = random_covariance(cor = ...)`.
#' @param monitor a `random_monitor()` object. In `prior_random()`, it supplies
#'   the top-level default; in `random_block()`, `NULL` inherits that default.
#' @param new_levels a `random_new_levels()` object controlling prediction for
#'   grouping levels absent from the fitted data. In `random_block()`, `NULL`
#'   inherits the top-level policy.
#' @param allocation a `random_variance_allocation()` object or list of such
#'   objects. Allocations can equivalently be supplied through `...`, but the
#'   two routes cannot be combined.
#' @param parameterization requested sampled representation:
#'   `"noncentered"`, `"centered"`, or `"auto"`. The top-level value is
#'   inherited by blocks whose `random_block()` value is `NULL`.
#' @param terms in `random_block()`, an optional named list of term-specific SD
#'   overrides. Each entry is an SD prior or `random_block(sd = ...)`, and its
#'   name must match a resolved term such as `"intercept"` or a slope label. In
#'   `random_variance_allocation()`, a character vector of targeted block names;
#'   `NULL` lets a single root block allocation target all remaining blocks at
#'   resolution time.
#' @param contrasts for `random_block()`, an optional named character vector
#'   assigning contrast families to factor predictors in that block. Accepted
#'   values are `"treatment"`, `"independent"`, `"orthonormal"`, `"meandif"`,
#'   `"cumulative"`, and `"cumulative_levels"`, with or without the `contr.`
#'   prefix. Structure-owned index families do not accept this argument.
#' @param sd_source an external SD source created with `random_sd_source()`.
#'   In `random_block()`, it conflicts with block-local `sd`, covariance-owned
#'   `sd`, and term-specific SD overrides. In a root allocation, exactly one of
#'   `sd` and `sd_source` is required; child allocations inherit from `parent`.
#' @param name required stable allocation identifier used internally and by
#'   `allocation_ref()`.
#' @param display_name public semantic owner used as the parameter-name prefix.
#'   It defaults to `name`; use `""` when a single model component needs no
#'   disambiguating prefix.
#' @param component_names optional public names for allocation components, in
#'   the same order as `terms`. Internal component identifiers remain the
#'   sanitized names of `terms`.
#' @param parent an optional `allocation_ref()` selecting one component of an
#'   earlier named allocation. A child allocation must not specify `sd` or
#'   `sd_source`.
#' @param target allocation target: `"block"` divides a variance budget among
#'   blocks or symbolic components; `"sd_component"` divides it among the
#'   resolved SD components of exactly one heterogeneous block.
#' @param scale variance normalization. `"total_variance"` uses
#'   `sd_child = sd_parent * sqrt(w)`. `"mean_variance"` uses
#'   `sd_child = sd_parent * sqrt(K * w)` and is available only for
#'   `target = "sd_component"`.
#' @param weights an optional Dirichlet prior over allocation weights. For a
#'   block allocation with explicit `terms`, omission creates a symmetric
#'   `Dirichlet(1, ..., 1)` prior.
#' @param inclusion an optional named list of scalar probability priors. For a
#'   block allocation, each named component receives an independent Bernoulli
#'   inclusion gate.
#' @param allocation_name name of an earlier allocation referred to by
#'   `allocation_ref()`.
#' @param component component label within the referenced allocation. Named
#'   allocation `terms` define labels; unnamed terms use sanitized block names.
#' @param structure optional covariance-family label for
#'   `random_covariance()`: `"ID"`, `"DIAG"`, `"CS"`, `"HCS"`, `"AR"` /
#'   `"AR1"`, `"HAR"`, `"CAR"`, or `"UN"` / `"US"`. Matching is
#'   case-insensitive. In RoBMA it must agree with the formula-owned family.
#' @param eta positive LKJ concentration. Values above one increasingly favor
#'   correlations near zero; `eta = 1` is uniform over correlation matrices.
#' @param cor_scale scale for an explicitly supplied scalar `cor` prior.
#'   `"fisher_z"` places it on `atanh(cor)` and is the default; `"cor"` uses
#'   raw correlation; `"logit"` uses a logit transform of the admissible raw
#'   interval. An omitted `cor` always receives the raw-uniform default after
#'   the structure and dimension are resolved.
#' @param include_correlation whether `prior_lkj()` requests monitoring of the
#'   semantic correlation matrix.
#' @param include_primitives whether `prior_lkj()` requests monitoring of its
#'   primitive beta coordinates. These coordinates remain internal and are not
#'   added to the map's public quantities.
#' @param latent whether to monitor standardized latent random effects. These
#'   are required for bridge sampling and random-coefficient reconstruction.
#' @param coefficients whether to monitor realized group-level coefficients.
#'   This can be convenient for inspection but memory intensive.
#' @param correlation whether semantic correlation matrices are included in
#'   monitored output. Scalar structures may retain an internal backend
#'   coordinate during sampling and derive the requested correlation matrix
#'   afterward.
#' @param lkj_primitives whether to retain LKJ primitive beta coordinates for
#'   backend inspection. These coordinates remain internal and are not public
#'   plotting, density, or hypothesis targets.
#' @param method new-level prediction method. `"error"` rejects unseen levels,
#'   `"zero"` assigns zero random contribution, and `"sample"` draws from the
#'   fitted random-effect distribution.
#' @param source existing JAGS node name, or a BayesTools `parameter_source`
#'   object, used by `random_sd_source()`.
#' @param shape source shape for `random_sd_source()`: `"scalar"` for one node
#'   or `"row"` for a vector aligned with the formula rows.
#'
#' @return
#' `prior_random()` returns a `prior_random` object. The remaining functions
#' return the corresponding list-like helper objects consumed by
#' `prior_random()`.
#'
#' @examples
#' sd_prior <- prior(
#'   "normal",
#'   parameters = list(mean = 0, sd = 0.5),
#'   truncation = list(lower = 0, upper = Inf)
#' )
#'
#' # One random coefficient per group level, with unrestricted correlation.
#' random_formula <- ~ us(0 + group | study)
#' unstructured_prior <- prior_random(
#'   study = random_block(
#'     contrasts = c(group = "independent")
#'   )
#' )
#'
#' # The same coefficient basis can instead have a diagonal covariance.
#' diagonal_formula <- ~ diag(0 + group | study)
#' diagonal_prior <- prior_random(
#'   study = random_block(
#'     contrasts = c(group = "independent")
#'   )
#' )
#'
#' # Correlation-only customization remains partial: RoBMA supplies the SDs.
#' cor_prior <- prior("normal", list(mean = 0, sd = 1))
#' hcs_prior <- prior_random(
#'   study = random_block(
#'     covariance = random_covariance(
#'       cor = cor_prior,
#'       cor_scale = "fisher_z"
#'     )
#'   )
#' )
#'
#' # Allocate one common SD scale across a block's level-specific SD components.
#' component_allocation <- random_variance_allocation(
#'   name = "heterogeneity",
#'   display_name = "",
#'   terms = "study",
#'   sd = sd_prior,
#'   weights = prior("dirichlet", list(alpha = c(1, 1))),
#'   target = "sd_component",
#'   scale = "mean_variance"
#' )
#' allocated_prior <- prior_random(
#'   study = random_block(
#'     covariance = random_covariance(cor = prior_lkj(eta = 1)),
#'     contrasts = c(group = "independent")
#'   ),
#'   allocation = component_allocation
#' )
#'
#' # A complete AR prior owns both its SD and its scalar correlation prior.
#' ar_cor_prior <- prior("normal", list(mean = 0, sd = 0.5))
#' autoregressive_prior <- prior_random(
#'   study = random_block(
#'     sd = sd_prior,
#'     covariance = random_covariance(
#'       cor = ar_cor_prior,
#'       cor_scale = "fisher_z"
#'     )
#'   )
#' )
#'
#' # Monitoring and new-level prediction policies are explicit when customized.
#' prediction_prior <- prior_random(
#'   sd = sd_prior,
#'   monitor = random_monitor(latent = TRUE, coefficients = FALSE),
#'   new_levels = random_new_levels(method = "sample")
#' )
#'
#' # Advanced: refer to an existing row-shaped model node as the SD source.
#' external_sd <- random_sd_source("tau", shape = "row")
#'
#' @seealso
#' [random_effect_prior_specification()], [random_effect_formula_tags()],
#' [prior()], [brma.mv()]
#'
#' @name prior_random
NULL


#' @rdname prior_random
#' @importFrom BayesTools prior_random
#' @export
BayesTools::prior_random

#' @rdname prior_random
#' @importFrom BayesTools random_block
#' @export
BayesTools::random_block

#' @rdname prior_random
#' @importFrom BayesTools random_variance_allocation
#' @export
BayesTools::random_variance_allocation

#' @rdname prior_random
#' @importFrom BayesTools allocation_ref
#' @export
BayesTools::allocation_ref

#' @rdname prior_random
#' @importFrom BayesTools random_covariance
#' @export
BayesTools::random_covariance

#' @rdname prior_random
#' @importFrom BayesTools prior_lkj
#' @export
BayesTools::prior_lkj

#' @rdname prior_random
#' @importFrom BayesTools random_monitor
#' @export
BayesTools::random_monitor

#' @rdname prior_random
#' @importFrom BayesTools random_new_levels
#' @export
BayesTools::random_new_levels

#' @rdname prior_random
#' @importFrom BayesTools random_sd_source
#' @export
BayesTools::random_sd_source
