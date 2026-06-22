#' @title Random-Effect Formula Structure Tags
#' @name random_effect_formula_tags
#'
#' @description
#' `brma.mv()` random-effect formulas use BayesTools random-effect formula
#' structure tags to describe the covariance structure of random coefficients.
#' These tags are formula specials parsed by
#' [BayesTools::random_effects_formula()], not ordinary user-facing functions.
#'
#' @details
#' Use tags as bare names inside the `random` formula supplied to `brma.mv()`,
#' for example `random = ~ cs(comp | study)`.
#' They are parsed from the formula call and are not called as normal R
#' functions, so users should not write `RoBMA::cs()` or `BayesTools::cs()`.
#' RoBMA forwards the formula to [BayesTools::random_effects_formula()], stores
#' the parsed BayesTools random-effect metadata, and uses that metadata for
#' JAGS compilation, prior matching, prediction, and summary labels.
#'
#' Basic syntax is `random = ~ tag(term | group)`.
#' The left side of `|` defines the random coefficient design and the right
#' side defines the grouping factor.
#' Plain random-intercept syntax, such as `random = ~ 1 | study`, is accepted
#' for intercept-only blocks.
#' Random slopes require an explicit structure tag, such as
#' `random = ~ diag(1 + x | study)` or `random = ~ us(1 + x | study)`.
#' The lme4-style `||` shorthand is accepted as a diagonal covariance shortcut.
#'
#' Supported structure tags:
#' \itemize{
#'   \item `id()`: identity covariance; random coefficients are independent
#'     with a common standard deviation.
#'   \item `diag()`: diagonal covariance; random coefficients are independent
#'     with coefficient-specific standard deviations.
#'   \item `us()` / `un()`: unstructured covariance; standard deviations and
#'     correlations are estimated freely.
#'   \item `cs()`: compound-symmetric covariance; coefficients share one
#'     standard deviation and one common correlation.
#'   \item `hcs()`: heterogeneous compound symmetry; coefficients have
#'     coefficient-specific standard deviations and one common correlation.
#'   \item `ar1()` / `ar()`: homogeneous autoregressive order-1 covariance;
#'     coefficients share one standard deviation and correlations decay as
#'     \eqn{\rho^{|i-j|}} over the ordered coefficient levels.
#'   \item `har()`: heterogeneous autoregressive order-1 covariance;
#'     coefficients have coefficient-specific standard deviations and AR(1)
#'     correlations.
#'   \item `car()`: continuous-time autoregressive covariance; use exactly one
#'     untransformed time coordinate, as in `car(time | study)`, with
#'     correlations \eqn{\rho^d} based on absolute time distances. The
#'     correlation parameter is constrained to \eqn{[0, 1)}.
#'   \item `random()` / `re()`: explicit random-effect block wrapper. By
#'     default this uses an unstructured covariance, and it can be useful when
#'     specifying an explicit `name` for prior matching and summaries.
#' }
#' Alias tags are normalized before model compilation: `un()` to `us()`,
#' `ar()` to `ar1()`, and unqualified `random()` / `re()` to `us()` unless a
#' covariance override changes the structure.
#'
#' The scalar-correlation structures `cs()`, `hcs()`, `ar1()` / `ar()`,
#' `har()`, and `car()` use a default `Beta(1, 1)` prior distribution directly
#' on the raw correlation parameter \eqn{\rho}. This default prior distribution
#' is bounded to non-negative correlations. To allow negative correlations,
#' supply an explicit [BayesTools::prior_random()] object with a
#' [BayesTools::random_covariance()] prior on either the raw or Fisher-z scale.
#'
#' A named list of formulas defines top-level random-effect components:
#'
#' ```r
#' random = list(
#'   study  = ~ cs(comp | study),
#'   design = ~ cs(comp | design)
#' )
#' ```
#'
#' List names define component names used by RoBMA summaries and by
#' component-specific scale models.
#' Explicit `name =` arguments inside `random()` or `re()` define the final
#' random-effect block name used by `BayesTools::prior_random()` block
#' overrides.
#'
#' @seealso [brma.mv()], [BayesTools::random_effects_formula()],
#'   [BayesTools::random_covariance()], [BayesTools::prior_random()]
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
#' }
NULL
