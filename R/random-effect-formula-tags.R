#' @title Random-Effect Formula Structure Tags
#' @name random_effect_formula_tags
#'
#' @description
#' `brma.mv()` random-effect formulas use BayesTools random-effect formula
#' structure tags to describe the covariance structure of random coefficients.
#' These tags are formula specials parsed by
#' `BayesTools::random_effects_formula()`, not ordinary user-facing functions.
#'
#' @details
#' Use tags as bare names inside the `random` formula supplied to `brma.mv()`,
#' for example `random = ~ cs(comp | study)`.
#' They are parsed from the formula call and are not called as normal R
#' functions, so users should not write `RoBMA::cs()` or `BayesTools::cs()`.
#' RoBMA forwards the formula to `BayesTools::random_effects_formula()`, stores
#' the parsed BayesTools random-effect metadata, and uses that metadata for
#' JAGS compilation, prior matching, prediction, and summary labels.
#'
#' BayesTools distinguishes random-coefficient formulas from structured index
#' specifications. For `id()`, `diag()`, and `us()` / `un()`, the left side of
#' `|` is a coefficient formula. It can contain an intercept, continuous or
#' factor slopes, and interactions. As in an ordinary random-slope formula,
#' `1`, `0`, and `-1` include or remove the random intercept. Plain
#' `(expr | group)` defaults to `us(expr | group)`, while `expr || group` is the
#' diagonal shorthand.
#'
#' Factor coding inside a random-coefficient block is not selected by removing
#' the intercept. BayesTools reuses an already established fixed-factor
#' contrast by default and lets a
#' \code{\link{random_block}} set an explicit block-specific basis through
#' `contrasts`. For example, `us(0 + group | study)` together with
#' `random_block(contrasts = c(group = "independent"))` creates one correlated
#' random coefficient per `group` level and no random intercept. Without the
#' independent contrast override, `0 + group` preserves the resolved contrast
#' basis instead of forcing raw level indicators.
#'
#' In contrast, `cs()`, `hcs()`, `ar1()` / `ar()`, `har()`, and `car()` own a
#' structure-specific index basis. Their left side names index columns rather
#' than coefficients. These tags reject explicit `1`, `0`, and `-1`, and their
#' blocks do not accept `random_block(contrasts = ...)`.
#'
#' Supported structure tags:
#' \itemize{
#'   \item `id()`: independent random-coefficient columns with one shared
#'     standard deviation.
#'   \item `diag()`: independent random-coefficient columns with one standard
#'     deviation per column. `||` is shorthand for this structure.
#'   \item `us()` / `un()`: one standard deviation per random-coefficient column
#'     and an unrestricted correlation matrix. Plain `|` syntax defaults to
#'     this structure.
#'   \item `cs()`: compound symmetry over one or more discrete index columns;
#'     index levels share one standard deviation and one common pairwise
#'     correlation.
#'   \item `hcs()`: heterogeneous compound symmetry over one or more discrete
#'     index columns; levels have separate standard deviations and one common
#'     pairwise correlation.
#'   \item `ar1()` / `ar()`: homogeneous autoregressive order-1 covariance over
#'     exactly one discrete index; levels share one standard deviation and
#'     correlations decay as \eqn{\mathrm{cor}^{|i-j|}} over the resolved order.
#'   \item `har()`: heterogeneous autoregressive order-1 covariance;
#'     index levels have separate standard deviations and AR(1) correlations.
#'   \item `car()`: continuous-time autoregressive covariance; use exactly one
#'     finite numeric/integer coordinate or an ordered factor with numeric level
#'     labels, as in `car(time | study)`. Correlations are
#'     \eqn{\mathrm{cor}^d} based on actual absolute distances, and
#'     \eqn{\mathrm{cor}} is constrained to \eqn{[0, 1)}.
#'   \item `random()` / `re()`: explicit random-effect block wrapper. By
#'     default this uses an unstructured covariance, and it can be useful when
#'     specifying an explicit `name` for prior matching and summaries.
#' }
#'
#' The discrete structured tags `cs()`, `hcs()`, `ar1()`, and `har()` accept
#' factor, character, numeric/integer, or logical index columns. Existing factor
#' level order is retained; otherwise sorted unique values determine the levels
#' and, for AR(1), their order. `cs()` and `hcs()` may combine multiple columns
#' with `+`; BayesTools forms the observed interaction in lexicographic order.
#' `ar1()` and `har()` require exactly one index column.
#'
#' `hcs(index | group)` and `us(0 + index | group)` are different models even
#' when both produce one standard deviation per displayed level. HCS constrains
#' every pair to the same correlation. US estimates an unrestricted correlation
#' matrix among the selected coefficient columns and its factor basis must be
#' specified through contrast metadata when it matters.
#' Alias tags are normalized before model compilation: `un()` to `us()`,
#' `ar()` to `ar1()`, and unqualified `random()` / `re()` to `us()` unless a
#' covariance override changes the structure.
#'
#' `brma.mv()` also accepts metafor-style known group covariance through its
#' `R` and `Rscale` arguments. This covariance is attached to random-effect
#' grouping levels and is distinct from the known sampling covariance `V`.
#' RoBMA translates `R` internally to BayesTools random-effect group covariance
#' metadata; users should continue to use `R`/`Rscale` in `brma.mv()` rather
#' than BayesTools-native `group_covariance` arguments. Current support is
#' limited to random-intercept blocks such as `random = ~ 1 | study`. Supported
#' one-to-one known-`R` blocks in known-`V` normal models can be compiled as
#' marginalized diagonal row-variance terms; other supported known-`R` blocks
#' remain sampled.
#'
#' The scalar-correlation structures `cs()`, `hcs()`, `ar1()` / `ar()`,
#' and `har()` use uniform defaults over their complete valid raw-correlation
#' intervals: `(-1 / (K - 1), 1)` for CS/HCS and `(-1, 1)` for AR1/HAR, where
#' `K` is the number of resolved index levels. CAR uses `Uniform(0, 1)` because
#' its correlation-decay parameter is one-sided. These defaults therefore do
#' not silently exclude negative correlations when the structure permits them.
#'
#' A list of formulas defines top-level random-effect components:
#'
#' ```r
#' random = list(
#'   study  = ~ cs(comp | study),
#'   design = ~ cs(comp | design)
#' )
#' ```
#'
#' A bare formula or unnamed one-entry list has no redundant user-facing
#' component prefix. An explicitly named one-entry list retains its name. In a
#' list with two or more formulas, missing names are generated as
#' `component 1`, `component 2`, and so on. These component names are used by
#' RoBMA summaries and by component-specific scale models.
#' Explicit `name =` arguments inside `random()` or `re()` define the final
#' random-effect block name used by `prior_random()` block
#' overrides.
#'
#' ## Public parameter names
#'
#' BayesTools stores concrete coordinates, selectable semantic quantities, and
#' aliases in one versioned parameter map. Random-effect canonical names use
#' `(formula) owner: quantity(arguments)`, with the formula prefix optional when
#' the remaining name is unambiguous. Parentheses contain coefficient or
#' parameter names and square brackets contain factor or index levels. Examples
#' include `study: sd(intercept)` and
#' `study: cor(group[sensitivity],group[specificity])`.
#'
#' Public correlations are named `cor`; compact scalar `rho` and LKJ
#' construction coordinates remain internal. Automatic total-variance
#' allocations expose `sd_total`, `var_total`, and `var_prop(...)`.
#' Mean-variance allocations within heterogeneous blocks expose `sd_common`,
#' `var_common`, `var_ratio(...)`, and `sd_ratio(...)`. These semantic names are
#' the supported selectors for summaries, plotting, density estimation, and
#' hypotheses.
#'
#' @seealso [brma.mv()], \code{\link{random_effect_prior_specification}},
#'   \code{\link{prior_random}}, \code{\link{random_covariance}},
#'   \code{BayesTools::random_effects_formula()}
#'
#' @examples
#' \dontrun{
#' fit <- brma.mv(
#'   yi      = yi,
#'   V       = V,
#'   random  = ~ cs(outcome | study),
#'   data    = dat,
#'   measure = "GEN",
#'   prior_unit_information_sd = 1
#' )
#'
#' fit_network <- brma.mv(
#'   yi = yi,
#'   V  = V,
#'   random = list(
#'     study  = ~ cs(comp | study),
#'     design = ~ cs(comp | design)
#'   ),
#'   data    = dat,
#'   measure = "OR"
#' )
#'
#' group_prior <- prior_factor(
#'   distribution = "normal",
#'   parameters   = list(mean = 0, sd = 1),
#'   contrast     = "independent"
#' )
#' random_prior <- prior_random(
#'   study = random_block(
#'     contrasts = c(group = "independent")
#'   )
#' )
#' fit_group <- brma.mv(
#'   yi = yi, V = V, ni = ni,
#'   mods = ~ 0 + group,
#'   random = ~ us(0 + group | study),
#'   prior_mods = group_prior,
#'   prior_heterogeneity = random_prior,
#'   data = dat, measure = "GEN",
#'   prior_unit_information_sd = 1
#' )
#' }
NULL
