---
applyTo: "**/vignettes/*.Rmd"
---

# Vignette Writing Instructions for RoBMA

This document provides guidance for writing and maintaining vignettes in the RoBMA package.

## Overview

RoBMA vignettes are R Markdown documents that demonstrate package functionality with real-world examples. They are pre-computed and cached to avoid CRAN check timeouts, as Bayesian model fitting is computationally intensive.

## Vignette Structure

### Current Vignettes
1. **Tutorial.Rmd** - Introduction to RoBMA-PSMA (publication bias adjustment)
2. **ReproducingBMA.Rmd** - Classic Bayesian model-averaged meta-analysis (no publication bias)
3. **MetaRegression.Rmd** - `RoBMA.reg()` with moderators
4. **HierarchicalRoBMA.Rmd** - Multilevel RoBMA
5. **HierarchicalRoBMARegression.Rmd** - Multilevel RoBMA with moderators
6. **HierarchicalBMA.Rmd** - Simpler multilevel models via `study_ids`
7. **MedicineBMA.Rmd** - Informed priors for medical meta-analysis (continuous outcomes)
8. **MedicineBiBMA.Rmd** - Informed priors for binary outcomes (log OR, RR, RD, HR)
9. **CustomEnsembles.Rmd** - Advanced ensemble customization
10. **FastRoBMA.Rmd** - Spike-and-slab algorithm (`algorithm = "ss"`)
11. **ZCurveDiagnostics.Rmd** - Meta-analytic z-curve publication bias diagnostics

## Standard YAML Header

```yaml
---
title: "Your Vignette Title"
author: "Author Name(s)"
date: "`r Sys.Date()`"  # or fixed year for published papers
output:
  rmarkdown::html_vignette:
    self_contained: yes
bibliography: ../inst/REFERENCES.bib
csl: ../inst/apa.csl
vignette: >
  %\VignetteIndexEntry{Your Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown_notangle}
---
```

**Important**: Use `../inst/REFERENCES.bib` (relative path) for bibliography, not absolute paths.

## Code Chunk Strategy (Pre-computation Pattern)

All vignettes follow a **three-chunk pattern** to handle computationally expensive model fitting:

### Chunk 1: Setup & Check Detection
```r
```{r setup, include = FALSE}
is_check <- ("CheckExEnv" %in% search()) ||
              any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv())) ||
              !file.exists("../models/YourVignette/your_model.RDS")
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  eval      = !is_check,
  dev       = "png")
if(.Platform$OS.type == "windows"){
  knitr::opts_chunk$set(dev.args = list(type = "cairo"))
}
```
```

**Purpose**: Detect CRAN checks or missing cached models and disable evaluation to avoid timeouts.

### Chunk 2: Load Pre-computed Models
```r
```{r include = FALSE}
library(RoBMA)
# Pre-load fitted models to avoid re-fitting during vignette build
fit_model <- readRDS(file = "../models/YourVignette/your_model.RDS")
```
```

**Purpose**: Load cached model results silently (not shown to user).

### Chunk 3: Model Fitting Code (Not Evaluated)
```r
```{r include = FALSE, eval = FALSE}
# R package version updating
library(RoBMA)

# Actual model fitting code that was used to create cached models
fit_model <- RoBMA(d = data$d, se = data$se, seed = 1, parallel = TRUE)

# Save for future vignette builds
saveRDS(fit_model, file = "../models/YourVignette/your_model.RDS")
```
```

**Purpose**: Document the exact code used to generate cached models. This is **never evaluated** during package checks but serves as a record for updating models when package versions change.

### Why This Pattern?

- **CRAN compliance**: Vignettes must build in < 10 minutes; MCMC fitting takes much longer
- **Reproducibility**: Exact fitting code is preserved but not executed
- **Version tracking**: When RoBMA updates, re-run chunk 3 to regenerate all cached models
- **User clarity**: Users see the actual fitting code in chunk 3 (via `include = FALSE` it doesn't clutter output)

## Model Caching Location

All pre-computed models are stored in `models/` directory:
```
models/
  Tutorial/
    fit_RoBMA_Lui2015.RDS
  ReproducingBMA/
    PowerPoseTest.RDS
  MetaRegression/
    fit_RoBMA.RDS
  ...
```

- **Naming convention**: Use descriptive names (dataset + model type)
- **Compression**: Use `compress = "xz"` for large models: `saveRDS(fit, file = "path.RDS", compress = "xz")`
- **Git tracking**: Models are committed to the repository (not gitignored)

## Code Presentation for Users

Code that **users should see and run** goes in regular chunks:

```r
```{r}
library(RoBMA)
data("Lui2015", package = "RoBMA")
head(Lui2015)
```
```

**Never show** the model loading code (`readRDS()`) to users. They should see the fitting code from chunk 3.

## Displaying Pre-computed Results

After loading cached models with `readRDS()`, display them normally:

```r
```{r}
# This uses the pre-loaded fit_model from chunk 2
summary(fit_model)
plot(fit_model, parameter = "mu")
```
```

Users see the output without knowing it came from cache.

## Citations

Use `\insertCite{key}{RoBMA}` for inline citations:
- `\insertCite{bartos2021no}{RoBMA}` → (Bartoš et al., 2021)
- `\insertCite{bartos2021no;textual}{RoBMA}` → Bartoš et al. (2021)

Add new references to `inst/REFERENCES.bib`. The bibliography is automatically rendered at the end.

## Code Style in Vignettes

- **Function calls**: Use full argument names for clarity (no abbreviations)
- **Seeds**: Always set `seed = 1` (or another fixed value) for reproducibility
- **Parallel processing**: Use `parallel = TRUE` when fitting to speed up model generation
- **Save argument**: Consider `save = "min"` to reduce model size if posterior samples aren't needed

### Example
```r
fit <- RoBMA(
  d       = data$effectSize,
  se      = data$SE,
  seed    = 1,
  parallel = TRUE,
  save    = "min"  # Reduces file size
)
```

## Figures

- **Captions**: Use `fig.cap` for meaningful captions
  ```r
  ```{r, fig.cap="Forest Plot of Effect Sizes"}
  forest(fit_model)
  ```
  ```
- **Size**: Let knitr use defaults; override only if necessary
- **Device**: The setup chunk handles Windows Cairo device automatically

## Updating Vignettes for New Package Versions

When RoBMA is updated and model structures change:

1. **Identify affected vignettes** (check NEWS.md for breaking changes)
2. **Re-run chunk 3** in each affected vignette:
   ```r
   # Set eval = TRUE temporarily in chunk 3 header
   ```{r include = FALSE, eval = TRUE}
   ```
3. **Verify outputs** match expectations
4. **Commit updated .RDS files** to `models/`
5. **Reset chunk 3** back to `eval = FALSE`
6. **Rebuild vignettes**: `devtools::build_vignettes()`

## Testing Vignettes Locally

```r
# Build all vignettes
devtools::build_vignettes()

# Preview specific vignette
rmarkdown::render("vignettes/Tutorial.Rmd")

# Check if vignettes build during R CMD check
devtools::check()
```

## Common Pitfalls

❌ **Don't** use `library()` or `require()` in package functions (only in vignettes is OK)
❌ **Don't** use absolute paths (`C:/Users/...`)
❌ **Don't** commit temporary files (`.html` vignette outputs go to `doc/`)
❌ **Don't** use `eval = TRUE` in chunk 3 (model fitting) unless intentionally regenerating
✅ **Do** use relative paths (`../models/`, `../inst/`)
✅ **Do** compress models (`compress = "xz"`)
✅ **Do** test that vignettes build with `is_check = TRUE` condition (simulates CRAN)

## Example Vignette Skeleton

```rmd
---
title: "My New RoBMA Vignette"
author: "Your Name"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    self_contained: yes
bibliography: ../inst/REFERENCES.bib
csl: ../inst/apa.csl
vignette: >
  %\VignetteIndexEntry{My New RoBMA Vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
is_check <- ("CheckExEnv" %in% search()) ||
              any(c("_R_CHECK_TIMINGS_", "_R_CHECK_LICENSE_") %in% names(Sys.getenv())) ||
              !file.exists("../models/MyVignette/my_model.RDS")
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  eval      = !is_check,
  dev       = "png")
if(.Platform$OS.type == "windows"){
  knitr::opts_chunk$set(dev.args = list(type = "cairo"))
}
```

```{r include = FALSE}
library(RoBMA)
my_model <- readRDS(file = "../models/MyVignette/my_model.RDS")
```

```{r include = FALSE, eval = FALSE}
library(RoBMA)
data("MyData", package = "RoBMA")

my_model <- RoBMA(d = MyData$d, se = MyData$se, seed = 1, parallel = TRUE)
saveRDS(my_model, file = "../models/MyVignette/my_model.RDS")
```

## Introduction

This vignette demonstrates...

```{r}
library(RoBMA)
data("MyData", package = "RoBMA")
head(MyData)
```

## Analysis

```{r}
summary(my_model)
```

## References
```

## Prose Editing Guidelines

When editing vignette prose, follow the Eric-Jan Wagenmakers style: concise, direct, and logically structured. Clarify meaning, tighten flow, and preserve all scientific content.

### Writing Style & Formatting
- **Concise and Direct**: Use simple sentences to describe outputs. Avoid flowery language or filler phrases.
- **No Excessive Bold**: Use bold text sparingly. Do not bold every list item or emphasis point. Use it only for headers or defining key terms.
- **Flowing Text**: Prefer paragraphs over bulleted lists when describing plots or outputs. Integrate the description into a narrative flow.
- **Interpretation Focused**: Focus on what the output *means* (interpretation) rather than just listing what is displayed.
- **Concrete Examples**: Use specific values from the example to illustrate points (e.g., "In our example, we find...").
- **Technical but Accessible**: Use correct terminology (e.g., "heterogeneity allocation parameter") but explain it simply.

### Non-Negotiables
- **Do not** add/remove references, change results, mathematical notation, or variable names
- **Preserve UI specifics exactly**: argument names like `priors_effect`, function names like `RoBMA.reg()`, parameter names like `mu`, `tau`, `omega`
- **Keep defined abbreviations**; spell out on first use (e.g., "Markov Chain Monte Carlo (MCMC)")
- **Prefer full terms**: "prior distributions" over "priors"; spell out "null hypothesis" and "alternative hypothesis"
- **Do not omit technical details**: exact argument labels, full file paths, figure references

### Voice & Rhythm
- **Prefer passive tense** for objectivity, but use collaborative first-person plural ("we set...", "we estimate...") when it improves flow
- **Avoid "we... we..." runs**: vary sentence structure to maintain rhythm
- **Keep tone precise and readable**: cut redundancy, avoid filler phrases, use commas for disambiguation only

### Editing Passes (Apply in Order)

1. **Meaning**: Remove clutter; define key terms briefly when first introduced; add a concrete example if needed for clarity
2. **Structure**: Smooth transitions between paragraphs; align parallel or contrasting ideas; keep section logic tight
3. **Emphasis & Rhythm**: Place key words in strong positions (sentence start/end); use light anaphora/epistrophe only if it clarifies
4. **Style**: One tasteful rhetorical device per sentence at most (e.g., parallelism, anticipating objections); maintain EJW tone
5. **Polish**: Fix punctuation for disambiguation; correct typos quietly

### Clarity Techniques (Use Sparingly)
- **Parallelism**: Align list items or related sentences for easier comparison
- **Procatalepsis**: Anticipate and answer likely reader objections in one sentence when helpful
- **Selective repetition**: Repeat key terms for emphasis, but avoid redundancy

### Examples

❌ **Verbose**: "In this section, we are going to discuss how to fit models using the RoBMA package"
✅ **Concise**: "We fit models using the `RoBMA()` function"

❌ **Vague**: "We can use different priors for the analysis"
✅ **Specific**: "We specify prior distributions via the `priors_effect` and `priors_heterogeneity` arguments"

❌ **Redundant**: "The results show that the effect is significant and statistically significant"
✅ **Tight**: "The effect is statistically significant"

❌ **Cluttered**: "We can see from the output that..."
✅ **Direct**: "The output shows..."

❌ **Excessive Bold/Lists**:
> This plot displays:
> - **x-axis**: One-sided *p*-value cutoffs
> - **y-axis**: Relative probability of publication

✅ **Flowing Description**:
> The plot displays one-sided *p*-value cutoffs (x-axis) against relative publication probability (y-axis).

## Additional Resources

- [R Markdown Guide](https://rmarkdown.rstudio.com/articles_intro.html)
- [Vignette Best Practices](https://r-pkgs.org/vignettes.html)
- RoBMA paper: \insertCite{bartos2022adjusting}{RoBMA}
