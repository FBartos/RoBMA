# Vignette Writing Instructions

Reference for writing vignettes in `vignettes/`. Canonical references for general functionality:

- **Prior distributions / default prior distributions**: `vignettes/v01-prior-distributions.Rmd`
- **Getting-started / simple workflow**: `vignettes/v02-bayesian-meta-analysis.Rmd`
- **Comprehensive feature list**: `vignettes/v03-feature-coverage.Rmd`

## Filename and Ordering

Vignettes are ordered alphabetically by filename in `browseVignettes()`, on the CRAN page, and in pkgdown. We use a tiered numeric prefix with gaps so new vignettes can be inserted without renumbering:

- `0X` — Foundations (`v00-introduction.Rmd`, `v01-prior-distributions.Rmd`, `v02-bayesian-meta-analysis.Rmd`)
- `1X` — Correspondence with `metafor` (`v10-metafor-parity-multilevel.Rmd`, ...)
- `2X` — Bayesian Model Averaging (`v20-bayesian-model-averaging.Rmd`, `21-publication-bias-adjustment.Rmd`)
- `3X` — Paper companions (`v30-tutorial.Rmd`, `31-meta-regression.Rmd`, ...)

Filenames are lowercase with hyphens. The YAML `title:` is the clean human-readable name, without the prefix.

When a vignette caches model fits, the cache directory under `models/` MUST match the filename without `.Rmd`. For example, `v02-bayesian-meta-analysis.Rmd` caches into `models/v02-bayesian-meta-analysis/`. The cache `name` argument inside the vignette must use the same string.

## YAML Header

```yaml
---
title: "Vignette Title"
author: "Author Name"
date: "27th of April 2026"
output:
  rmarkdown::html_vignette:
    self_contained: yes
bibliography: ../inst/REFERENCES.bib
csl: ../inst/apa.csl
link-citations: true
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown_notangle}
---
```

Use author `František Bartoš` unless instructed otherwise.

Use a spelled date, not `\`r Sys.Date()\``. `link-citations: true` and the `_notangle` engine line are required.

Immediately after the YAML header, include the shared print-output rule:

````markdown
```{r child = "_vignette-nowrap.md", echo = FALSE, eval = TRUE}
```
````

This child sets `options(width = 10000)` and CSS overrides for `pre` / `sourceCode` blocks, so printed output and source blocks use horizontal scrolling instead of browser wrapping.
Do not duplicate this CSS inside individual vignettes; edit `vignettes/_vignette-nowrap.md` if the rule changes.

## Setup Chunk

For prior-only or display-only vignettes (no MCMC fits), the setup is minimal:

```r
{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse   = TRUE,
  comment    = "#>",
  message    = FALSE,
  warning    = FALSE,
  fig.width  = 6,
  fig.height = 4
)

library(RoBMA)

has_metafor <- requireNamespace("metafor", quietly = TRUE)
has_bcg     <- has_metafor && requireNamespace("metadat", quietly = TRUE)
```

Vertically align `=` in both `opts_chunk$set()` and the `requireNamespace` flag block. Gate any chunk that depends on a Suggests-only package via `eval = has_<flag>`.

### Figures

For vignettes that render plots (most non-display-only vignettes), use `fig.retina = 3` and the cairo device on Windows:

```r
{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse    = TRUE,
  comment     = "#>",
  message     = FALSE,
  warning     = FALSE,
  fig.width   = 7,
  fig.height  = 3.5,
  dev         = "png",
  fig.retina  = 3
)
if (.Platform$OS.type == "windows") {
  knitr::opts_chunk$set(dev.args = list(type = "cairo"))
}
```

`fig.retina = 3` saves PNGs at 3× pixel density and inserts the `<img>` tag at the original CSS width, giving sharp figures on retina/high-DPI screens at the same physical display size as the default. **Do not bump `dpi`** to improve quality: `dpi` scales both the file resolution *and* the displayed CSS size, so the figures grow on the page. Use `fig.retina` for resolution; use per-chunk `out.width` only when explicitly needed.

The global `fig.height = 3.5` is for single-panel chunks. Override per chunk for double-column comparison plots (see *Side-by-Side `metafor` Parity* below).

## Caching (only when vignettes fit MCMC models)

When a vignette runs actual MCMC fits, use `vignettes/vignette-cache.R`. The cache `name` MUST match the vignette filename without `.Rmd`. The pattern is three chunks: setup with cache → load cached fits → re-fit code with `eval = FALSE`.

`vignettes/v02-bayesian-meta-analysis.Rmd` and `vignettes/v30-tutorial.Rmd` are working examples. If the vignette only assembles priors (`only_priors = TRUE`) or otherwise runs cheaply, skip caching entirely — see `v01-prior-distributions.Rmd`.

## Section Structure

- **Lede.** One to three short paragraphs before the first heading. Two patterns depending on vignette type:
  - *Methodological framing* (priors, distributions, theory). Open with what the topic is. If the topic has a substantive asymmetry or tension (e.g., priors behave differently in estimation vs. testing/BMA), state it up front so the rest of the vignette has a frame. Cite the umbrella reference at the end of the lede.
  - *Getting-started framing* (`metafor` parity, paper companions, tutorials). Open with "this vignette is a starting point for ...". Name the example dataset, the side-by-side device, and the simple-to-complex arc in one sentence. Cite the package and the comparison reference (e.g., `[@RoBMA]`, `[@metafor]`).
  - In both patterns, avoid hype-ish framings: do not list features as a "toolkit", do not write "what changes is X" sentences, do not use marketing language.
- **Headings.** `##` for top-level sections, `###` for subsections. Inside a subsection, finer divisions use bolded sentence-leading paragraphs like `**Categorical moderators.**` — not `####`.
- **Closer.** Substantive vignettes end with a short *General Considerations and Reporting* (or analogous) section that re-iterates the main framing and gives reporting guidance, in bolded principle paragraphs.
- **References.** `## References` is the final heading; the bibliography is auto-rendered.

## Citations

Pandoc syntax, not `\insertCite{}`:

- `[@bartos2021no]` → (Bartoš et al., 2021)
- `@bartos2021no` → Bartoš et al. (2021)
- `[@key1; @key2]` for multiple keys

Always check `inst/REFERENCES.bib` before citing. Don't invent keys and don't add keys to `inst/REFERENCES.bib`; ask the user to add the keys manually!

## Cross-vignette links

Reference other vignettes as clickable links to their rendered HTML, not as plain italics:

```markdown
see the [*Prior Distributions*](v01-prior-distributions.html) vignette
```

This works in both pkgdown sites and `browseVignettes()` HTML. The path is relative to the rendered HTML of the current vignette (same directory). Keep the title italicized inside the link text.

## Code Chunk Style

The single most important convention is the **argument layout in fitting calls**. The user has spent considerable time tuning this; respect it.

### Argument layout: highlight the demonstrated arg, suppress the boilerplate

`data` and `only_priors` are *boilerplate*: they appear in nearly every call and carry no per-example information. They are crammed onto a single line and always go at the bottom of the call. The arguments that this particular example is *demonstrating* get their own line above the boilerplate so the reader's eye lands on them.

The order of arguments inside the call is always:

1. **Input** — outcome variables and effect-size spec (`yi`, `vi`/`sei`, `ni`, `measure`).
2. **Formulas** — `mods`, scale formulas, etc.
3. **Other important arguments** — the per-example highlights: `rescale_priors`, `prior_unit_information_sd`, `set_contrast_factor_predictors`, `prior_effect`, `prior_informed_field`, `model_type`, custom `prior_*` lists.
4. **Boilerplate** — `data` and fitting-control args (`only_priors`, `seed`, `parallel`, ...). Always last.

**Basic call (input + boilerplate only):**

```r
bcg_priors <- brma(
  yi = yi, vi = vi, measure = "RR",
  data = bcg, only_priors = TRUE
)
```

**One demonstrated argument (a tweak / decoration):**

```r
bcg_wider_priors <- brma(
  yi = yi, vi = vi, measure = "RR",
  rescale_priors = 2,
  data = bcg, only_priors = TRUE
)
```

```r
Havrankova_manual_priors <- brma(
  yi = y, sei = se, measure = "GEN",
  prior_unit_information_sd = 10,
  data = Havrankova2025, only_priors = TRUE
)
```

**Formula plus a structural argument:**

```r
Havrankova_treatment <- brma(
  yi = y, sei = se, ni = N, measure = "GEN",
  mods = ~ facing_customer + N,
  set_contrast_factor_predictors = "treatment",
  data = Havrankova2025, only_priors = TRUE
)
```

**Multiple structural overrides:**

```r
bcg_informed_priors <- brma(
  yi = yi, vi = vi, measure = "RR",
  prior_informed_field    = "medicine",
  prior_informed_subfield = "Cochrane",
  data = bcg, only_priors = TRUE
)
```

### Vertical `=` alignment within multi-line blocks

When a run of single-line named arguments sits at the same indentation level, align the `=`:

```r
prior_informed_field    = "medicine",
prior_informed_subfield = "Cochrane",
```

```r
prior_effect        = prior_informed("Oosterwijk"),
prior_heterogeneity = prior_informed("van Erp", parameter = "heterogeneity"),
```

Do not align across blocks separated by multi-line `prior(...)` constructs. Inside a `prior(...)` call, do not align — the args read naturally without it.

### Named arguments in `prior()`, `prior_factor()`, etc.

Always pass `distribution = "..."` explicitly. Don't rely on positional matching:

```r
prior(
  distribution = "normal",
  parameters = list(mean = 0, sd = 0.5)
)

prior_factor(
  distribution = "mnormal",
  parameters = list(mean = 0, sd = 5),
  contrast   = "meandif"
)
```

### Effect-size input

Prefer `vi = vi` over `sei = sqrt(vi)` when `escalc()` has already produced `vi`. Use `sei = se` only when the dataset stores SE directly (e.g., Havrankova2025).

### `print_prior()` and `plot_prior()`

Default to the unargumented form for general inspection:

```r
print_prior(fit)
```

This prints all priors and is the right call after a fit when the whole model is interesting. Use the focused forms only when comparing one specific component or moderator across fits:

```r
print_prior(fit, parameter      = "mu")          # focus on a single component
print_prior(fit, parameter      = "bias")        # focus on bias-adjustment priors
print_prior(fit, parameter_mods = "x")           # focus on one moderator
```

For continuous-moderator priors on both scales, use two `plot_prior()` calls:

```r
plot_prior(fit, parameter_mods = "x")
plot_prior(fit, parameter_mods = "x",
  standardized_coefficients = FALSE)
```

`print_prior()` does not support `standardized_coefficients`; use `plot_prior()` for that toggle.

### `escalc()`

Group the count/rate input variables on one line, then put `measure` and `data` on the next:

```r
bcg <- metafor::escalc(
  ai = tpos, bi = tneg, ci = cpos, di = cneg,
  measure = "RR", data = dat.bcg
)

Kroupova2021 <- metafor::escalc(
  ri = r, ni = sample_size, measure = "ZCOR",
  data = Kroupova2021
)
```

## Side-by-Side `metafor` Parity

A `metafor`-parity (or getting-started) vignette pairs each `metafor::rma()` call with the matching `brma()` call, building from a simple random-effects model up to meta-regression with a categorical moderator. Do not dump all `brma()` fits in one chunk; introduce one model per section, each with its own summary, plot, and diagnostic before moving on.

**Do not duplicate generic inference helpers across parity vignettes.** A specialized parity vignette (multilevel, GLMM, publication-bias) should only walk through the features that are *specific* to the model class it covers — for the multilevel vignette, that is `cluster`, the `tau`/`rho` parameterization, and `summary_heterogeneity()`'s within/between decomposition. Generic helpers that work the same way for any `brma()` fit (`mods`, `regplot()`, `predict()`, `bf()` / `loo_compare()`, `rstudent()`, `qqnorm()`, `plot(loo(fit))`) belong in `v02-bayesian-meta-analysis.Rmd`. Close the specialized vignette with a short *Other Inference Helpers* section that points back to `v02-bayesian-meta-analysis.html` and notes that the only thing that changes is the new structural argument.

### Per-section pattern

1. **Lead with the metafor call**, evaluated inline (cheap, no caching):

```{r fit1-metafor}
fit1_metafor <- metafor::rma(yi, vi, data = dat, method = "REML")
fit1_metafor
```

2. **Then the brma counterpart**, displayed with `eval = FALSE`. The cached fit is loaded by the hidden `load-models` chunk, so all subsequent computation uses it:

```{r fit1-brma, eval = FALSE}
fit1_brma <- brma(
  yi = yi, vi = vi, measure = "RR",
  data = dat, seed = 1
)
```

3. **Then the summaries, plots, and diagnostics**, regular eval, using the cached fit:

```{r fit1-brma-summary}
summary(fit1_brma, include_MCMC_diagnostics = FALSE)
```

`add_loo()` and `add_marglik()` are also displayed with `eval = FALSE` at the point they're discussed (the cached fits already have them attached, so subsequent `loo()` / `bf()` calls work).

### Function pairings

Make the metafor ↔ RoBMA correspondence explicit when introducing each one:

| metafor | RoBMA | Notes |
|:--|:--|:--|
| `confint(fit)` | `summary_heterogeneity(fit)` | profile-likelihood vs posterior summaries for `tau`, `tau^2`, `I^2`, `H^2` |
| `funnel(fit)` | `funnel(fit, sampling_heterogeneity = FALSE)` | RoBMA shows the full sampling distribution under the model by default; turn it off to match metafor |
| `regplot(fit, mod = "...")` | `regplot(fit, mod = "...")` | RoBMA plots a categorical moderator as a single factor, not by dummy column |
| `predict(fit, newmods = ...)` | `predict(fit, newdata = data.frame(...), type = "terms")` | brma takes a data frame on the original moderator scale |
| `rstudent(fit)` | `rstudent(fit)` | metafor: studentized (deletion); RoBMA: LOO-PIT (requires `add_loo()`) |
| `qqnorm(fit)` | `qqnorm(fit)` | match `xlim` / `ylim` between panels for direct visual comparison |
| `influence(fit)` | `influence(fit)` | DFFITS, Cook's distance, COVRATIO, leverage, DFBETAS where available |
| `AIC(fit1, fit2)` | `loo_compare(fit1, fit2)` | requires `add_loo()` on every fit being compared |
| (no analogue) | `bf(fit2, fit1)` | Bayes factor; requires `add_marglik()` on both fits |
| (no analogue) | `marginal_means(fit)` | estimated marginal means averaged over other moderators |
| (no analogue) | `pooled_effect(fit, transform = "EXP")` | natural-scale (RR / OR / IRR / HR) effect summary |
| (no analogue) | `plot(fit, parameter = "mu", prior = TRUE)` | posterior with prior overlay |

### Side-by-side plots: `par()` recipe

Double-column comparison plots use tuned `par()` settings so axis tick labels, axis titles, and panel titles fit the smaller panels:

```{r regplot-continuous, fig.width = 7, fig.height = 3.2}
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1), cex.axis = 0.8, cex.lab = 0.8, cex.main = 0.75)
metafor::regplot(fit2_metafor, main = "metafor", xlim = c(10, 60), ylim = c(-2, 0.5))
regplot(fit2_brma,             main = "RoBMA",   xlim = c(10, 60), ylim = c(-2, 0.5))
```

Conventions:

- Set `fig.width = 7, fig.height = 3.2` per chunk; the global default `fig.height = 3.5` is for single-panel chunks.
- `cex.axis = 0.8, cex.lab = 0.8, cex.main = 0.75` keeps tick labels, axis titles, and panel titles legible without dominating the plotting area.
- Match `xlim` and `ylim` across both panels so visual differences come from the models, not the axis ranges.
- Panel titles are `main = "metafor"` and `main = "RoBMA"`. Expand only when the panel needs disambiguation, e.g., `"metafor: random dummy"` / `"RoBMA: alloc factor"`.

For single-panel posterior plots (e.g., `plot(fit, parameter = "mu", prior = TRUE)`), tighten margins instead of using the double-column recipe:

```r
par(mar = c(4, 4, 1, 1))
plot(fit, parameter = "mu", prior = TRUE, xlim = c(-3, 3))
```

## Prose Style

- Concise and direct. Simple sentences. No flowery or marketing language.
- One sentence per markdown line (soft wrap). Paragraphs are separated by blank lines.
- **Avoid em dashes.** Use commas, colons, parentheses, or sentence breaks.
- Avoid prescriptive "you want to..." or "What X buys you" framings; prefer neutral descriptions.
- Tables: neutral column headers (e.g., *Approach* / *Description* / *Main arguments*). Don't include "Use when" columns.
- Cross-reference subsections by their human-readable title in italics (`see *Per-moderator overrides* below`), not by filename or numeric prefix.
- Bolded sentence-leading paragraphs are the right format for re-iterating principles or warning notes (`**Wide priors are not "uninformative" in testing or model averaging.**`).
- For tone, the lede should establish the asymmetry that motivates the rest of the vignette (e.g., the role of priors differs between estimation and hypothesis testing). Re-iterate the same point in the closing section.

### Voice and word choice

- Use both **Subject + active verb stuctures** ("This vignette illustrates", ...) and **"we" as the narrator.** ("We start by fitting...", "we can model-average over the two models..."). Do not relly too much on passive tense. 
- The vignette is going through the analysis with the reader, not delivering a brochure.
- **Nouns   Not "the vignette walks through" or bare imperatives. 
- **Name the specific referent of nouns.** Prefer "mixture of *the two* models" over "mixture of models"; "the *prior distribution*" over "the prior"; "without fitting *the models*" over "without fitting". When the noun has a definite referent in context, name it; tightening this away reads as math jargon.
- **Open the lede with the practical situation, not the abstract claim.** "In meta-analysis, analysts often have to balance different assumptions about the data-generating process and decide what kind of model to fit, which covariates to include, and what publication-bias adjustment method to use" beats "Choosing one meta-analytic model out of several plausible specifications discards the evidence the discarded models carry." The conceptual claim follows the practical setup; do not lead with abstract framing.
- **Concrete, visible action over abstract metaphor.** "`regplot()` shrinks the difference towards zero" beats "shrinks the line toward a flat fit when the inclusion BF is modest". Describe what the figure or output literally shows, not what it metaphorically represents.
- **Plain words over math jargon when describing visualizations.** "the corresponding visualization for the prior distribution" beats "the same decomposition for the prior"; "drops the null component on the effect (spike at zero)" beats "drops the spike". When introducing a technical term, follow with a parenthetical that grounds it.
- **No clever paired phrasings.** "remove the component" beats "force the component in or out". List alternatives with "or"; do not stack parallels for symmetry.
- **When introducing a paired naming convention, name the example.** "Every component has a paired `_null` argument (i.e. `prior_effect` and `prior_effect_null`)" beats "Every component has a paired `_null` argument". The example does the lifting that the rule alone cannot.
- **Cross-references go in narrative voice.** "In a later vignette [*Robust Bayesian Meta-Analysis*](v21-robust-bayesian-meta-analysis.html), we illustrate..." or "Note that we do not demonstrate all the diagnostics; see [*Bayesian Meta-Analysis*](v02-bayesian-meta-analysis.html)". Bare terminal pointers ("for X, see Y") are reserved for the closer's see-also list.
- always use "prior distribution" / "default prior distribution" in full
- always use "`RoBMA` R package" / "`metafor` package" in full (note the "R" only for `RoBMA`)

## Common Pitfalls

- Do not use absolute paths.
- Do not use `library()` in package R/ source; vignettes can use it.
- Do not use `\insertCite{}`; use Pandoc `[@key]` syntax.
- Do not align `=` across blocks separated by multi-line constructs.
- Do not use em dashes in prose.
- Do not invent citation keys; check `inst/REFERENCES.bib` first.
- Do not use `\`r Sys.Date()\`` in the YAML date; use a spelled date.
- Do not let printed chunk output wrap in the browser; include `_vignette-nowrap.md` immediately after the YAML header.
- Do not let the cache `name` and the vignette filename drift apart.
- Do not bump `dpi` to improve figure quality. Use `fig.retina = 3` instead. `dpi` scales both file resolution and displayed CSS size, so the figures grow on the page.
- Do not list features as a "toolkit" or use "what changes is X" framings in the lede.
- Do not reference other vignettes as plain italics; use `[*Title*](filename.html)` so the link is clickable.
- Do not put all `brma()` fits in one chunk in `metafor`-parity vignettes. Build models step by step, one section per added complexity (intercept-only → continuous moderator → categorical moderator).
