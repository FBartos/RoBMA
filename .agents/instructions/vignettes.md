# Vignettes

Use this guide for files under `vignettes/` and their cached fits under
`models/`. Follow the nearest current vignette when a local pattern differs
from this guide.

## Files and Headers

Vignette filenames use ordered lowercase prefixes such as `v02-`, `v14-`,
and `v22-`. Inspect the current directory instead of copying a fixed inventory.
A cache directory under `models/` must match the vignette filename without
`.Rmd`.

Keep sources UTF-8. Use a literal human-readable date, the actual author name,
and the standard `html_vignette` header used by neighboring files. Do not copy
an old date from a template.

Immediately after YAML, include:

````markdown
```{r child = "_vignette-nowrap.md", echo = FALSE, eval = TRUE}
```
````

The shared child owns print-width and no-wrap CSS.

## Setup and Dependencies

Use `knitr::opts_chunk$set()` consistently with neighboring vignettes. For
plots, use `fig.retina = 3`; do not increase `dpi` to sharpen figures. Use
the existing Windows cairo setup where applicable.

Vignettes may attach packages because they are user-facing examples. Gate
Suggests-only packages with `requireNamespace(..., quietly = TRUE)` flags.

## Cached Fits

Use `vignettes/vignette-cache.R` only for expensive fitted models. The cache
name must match the vignette filename. Follow the existing pattern:

1. declare the cache;
2. load cached objects in a hidden chunk;
3. show fitting code with `eval = FALSE`.

Do not add migrations for stale cached objects. Regenerate them with the current
RoBMA and BayesTools implementations. Cheap prior-only examples should not be
cached.

## Code Examples

Make the scientific point visually prominent:

1. outcome inputs;
2. formulas;
3. arguments demonstrated by the example;
4. boilerplate such as `data`, `only_priors`, `seed`, and fitting controls.

Keep compact related inputs on one line where that is the established vignette
style. Align runs of named arguments, but do not force alignment across nested
multi-line calls. Name `distribution` explicitly in prior constructors.

Use `vi = vi` when data already contain sampling variances. Use `sei` when
the dataset contains standard errors.

In metafor-parity vignettes, introduce one model at a time: metafor fit, RoBMA
counterpart, then summaries, plots, and diagnostics. Keep matching plot limits
for side-by-side comparisons. Avoid repeating generic helpers in every
specialized vignette; link to the main Bayesian meta-analysis vignette.

## Citations and Links

- Use Pandoc citations such as `[@key]`, not `\insertCite{}`.
- Verify every key in `inst/REFERENCES.bib`. Do not invent or silently add
  bibliography entries.
- Link other vignettes as
  `[*Bayesian Meta-Analysis*](v02-bayesian-meta-analysis.html)`.

## Prose

Write concise, active, non-marketing prose. Open with the practical problem,
then introduce the method. Use `##` and `###` headings; avoid unnecessary
fourth-level sections. End substantive guides with reporting considerations and
references.

Use specific referents and plain descriptions of visible output. Avoid hype,
clever parallel phrasing, and em dashes. Prefer complete names such as "prior
distribution," "the `RoBMA` R package," and "the `metafor` package."

When creating a new vignette, examine other existing vignette for the expected
sections, style, etc.
