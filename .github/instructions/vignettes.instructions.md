---
description: 'Vignette writing: Guidance for writing vignette documentation.'
applyTo: "**/vignettes/*.Rmd"
---

# Vignette Writing Instructions for RoBMA

Reference for writing vignettes in `vignettes/`. The canonical reference vignette is `vignettes/01-prior-distributions.Rmd`; mirror its conventions whenever uncertain.

## Filename and Ordering

Vignettes are ordered alphabetically by filename in `browseVignettes()`, on the CRAN page, and in pkgdown. We use a tiered numeric prefix with gaps so new vignettes can be inserted without renumbering:

- `0X` — Foundations (`00-introduction.Rmd`, `01-prior-distributions.Rmd`, `02-bayesian-meta-analysis.Rmd`)
- `1X` — Correspondence with `metafor` (`10-metafor-parity-multilevel.Rmd`, ...)
- `2X` — Bayesian Model Averaging (`20-bayesian-model-averaging.Rmd`, `21-publication-bias-adjustment.Rmd`)
- `3X` — Paper companions (`30-tutorial.Rmd`, `31-meta-regression.Rmd`, ...)

Filenames are lowercase with hyphens. The YAML `title:` is the clean human-readable name, without the prefix.

When a vignette caches model fits, the cache directory under `models/` MUST match the filename without `.Rmd`. For example, `02-bayesian-meta-analysis.Rmd` caches into `models/02-bayesian-meta-analysis/`. The cache `name` argument inside the vignette must use the same string.

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

Use a spelled date, not `\`r Sys.Date()\``. `link-citations: true` and the `_notangle` engine line are required.

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

## Caching (only when vignettes fit MCMC models)

When a vignette runs actual MCMC fits, use `vignettes/vignette-cache.R`. The cache `name` MUST match the vignette filename without `.Rmd`. The pattern is three chunks: setup with cache → load cached fits → re-fit code with `eval = FALSE`.

`vignettes/02-bayesian-meta-analysis.Rmd` and `vignettes/30-tutorial.Rmd` are working examples. If the vignette only assembles priors (`only_priors = TRUE`) or otherwise runs cheaply, skip caching entirely — see `01-prior-distributions.Rmd`.

## Section Structure

- **Lede.** One to three short paragraphs before the first heading. Open with what the topic is. If the topic has a substantive asymmetry or tension (e.g., priors behave differently in estimation vs. testing/BMA), state it up front so the rest of the vignette has a frame. Cite the umbrella reference at the end of the lede.
- **Headings.** `##` for top-level sections, `###` for subsections. Inside a subsection, finer divisions use bolded sentence-leading paragraphs like `**Categorical moderators.**` — not `####`.
- **Closer.** Substantive vignettes end with a short *General Considerations and Reporting* (or analogous) section that re-iterates the main framing and gives reporting guidance, in bolded principle paragraphs.
- **References.** `## References` is the final heading; the bibliography is auto-rendered.

## Citations

Pandoc syntax, not `\insertCite{}`:

- `[@bartos2021no]` → (Bartoš et al., 2021)
- `@bartos2021no` → Bartoš et al. (2021)
- `[@key1; @key2]` for multiple keys

Always check `inst/REFERENCES.bib` before citing. Don't invent keys; if the desired entry is missing, add it to the bib file or drop the citation.

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

## Prose Style

- Concise and direct. Simple sentences. No flowery or marketing language.
- One sentence per markdown line (soft wrap). Paragraphs are separated by blank lines.
- **Avoid em dashes.** Use commas, colons, parentheses, or sentence breaks.
- Avoid prescriptive "you want to..." or "What X buys you" framings; prefer neutral descriptions.
- Tables: neutral column headers (e.g., *Approach* / *Description* / *Main arguments*). Don't include "Use when" columns.
- Cross-reference subsections by their human-readable title in italics (`see *Per-moderator overrides* below`), not by filename or numeric prefix.
- Bolded sentence-leading paragraphs are the right format for re-iterating principles or warning notes (`**Wide priors are not "uninformative" in testing or model averaging.**`).
- For tone, the lede should establish the asymmetry that motivates the rest of the vignette (e.g., the role of priors differs between estimation and hypothesis testing). Re-iterate the same point in the closing section.

## Common Pitfalls

- Do not use absolute paths.
- Do not use `library()` in package R/ source; vignettes can use it.
- Do not use `\insertCite{}`; use Pandoc `[@key]` syntax.
- Do not align `=` across blocks separated by multi-line constructs.
- Do not use em dashes in prose.
- Do not invent citation keys; check `inst/REFERENCES.bib` first.
- Do not use `\`r Sys.Date()\`` in the YAML date; use a spelled date.
- Do not let the cache `name` and the vignette filename drift apart.
