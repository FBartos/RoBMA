# Vignette Writing Instructions

Reference this file when working with vignettes in `vignettes/`.

## Overview

RoBMA vignettes are R Markdown documents that demonstrate package functionality with real-world examples. They are pre-computed and cached to avoid CRAN check timeouts, as Bayesian model fitting is computationally intensive.

## Current Vignettes

1. **Tutorial.Rmd** - Introduction to RoBMA-PSMA (publication bias adjustment)
2. **ReproducingBMA.Rmd** - Classic Bayesian model-averaged meta-analysis
3. **MetaRegression.Rmd** - `RoBMA.reg()` with moderators
4. **HierarchicalRoBMA.Rmd** - Multilevel RoBMA
5. **HierarchicalRoBMARegression.Rmd** - Multilevel RoBMA with moderators
6. **HierarchicalBMA.Rmd** - Simpler multilevel models via `cluster`
7. **MedicineBMA.Rmd** - Informed priors for medical meta-analysis
8. **MedicineBiBMA.Rmd** - Informed priors for binary outcomes
9. **CustomEnsembles.Rmd** - Advanced ensemble customization
10. **FastRoBMA.Rmd** - Spike-and-slab algorithm (`algorithm = "ss"`)
11. **ZCurveDiagnostics.Rmd** - Meta-analytic z-curve publication bias diagnostics

## Standard YAML Header

```yaml
---
title: "Your Vignette Title"
author: "Author Name(s)"
date: "`r Sys.Date()`"
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

**Important**: Use `../inst/REFERENCES.bib` (relative path) for bibliography.

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

**Purpose**: Detect CRAN checks or missing cached models and disable evaluation.

### Chunk 2: Load Pre-computed Models

```r
```{r include = FALSE}
library(RoBMA)
my_model <- readRDS(file = "../models/YourVignette/my_model.RDS")
```
```

**Purpose**: Load cached model results silently (not shown to user).

### Chunk 3: Model Fitting Code (Not Evaluated)

```r
```{r include = FALSE, eval = FALSE}
library(RoBMA)
data("MyData", package = "RoBMA")

my_model <- RoBMA(d = MyData$d, se = MyData$se, seed = 1, parallel = TRUE)
saveRDS(my_model, file = "../models/YourVignette/my_model.RDS")
```
```

**Purpose**: Document the exact code used to generate cached models. Never evaluated during package checks.

### Why This Pattern?

- **CRAN compliance**: Vignettes must build in < 10 minutes; MCMC fitting takes much longer
- **Reproducibility**: Exact fitting code is preserved but not executed
- **Version tracking**: When RoBMA updates, re-run chunk 3 to regenerate cached models

## Model Caching Location

All pre-computed models are stored in `models/` directory:

```
models/
  Tutorial/
    fit_RoBMA_Lui2015.RDS
  ReproducingBMA/
    PowerPoseTest.RDS
  ...
```

- **Naming convention**: Use descriptive names (dataset + model type)
- **Compression**: Use `compress = "xz"` for large models
- **Git tracking**: Models are committed to the repository

## Citations

Use `\insertCite{key}{RoBMA}` for inline citations:

- `\insertCite{bartos2021no}{RoBMA}` -> (Bartos et al., 2021)
- `\insertCite{bartos2021no;textual}{RoBMA}` -> Bartos et al. (2021)

Add new references to `inst/REFERENCES.bib`.

## Code Style in Vignettes

- **Function calls**: Use full argument names for clarity
- **Seeds**: Always set `seed = 1` (or another fixed value) for reproducibility
- **Parallel processing**: Use `parallel = TRUE` when fitting
- **Save argument**: Consider `save = "min"` to reduce model size

## Prose Editing Guidelines

Follow concise, direct, and logically structured style:

- **Concise and Direct**: Use simple sentences. Avoid flowery language.
- **Flowing Text**: Prefer paragraphs over bulleted lists when describing outputs.
- **Interpretation Focused**: Focus on what the output *means* rather than just what is displayed.
- **Technical but Accessible**: Use correct terminology but explain simply.

### Non-Negotiables

- **Do not** add/remove references, change results, mathematical notation, or variable names
- **Preserve UI specifics exactly**: argument names, function names, parameter names
- **Keep defined abbreviations**; spell out on first use

### Examples

- Verbose: "In this section, we are going to discuss how to fit models"
- Concise: "We fit models using the `RoBMA()` function"
- Vague: "We can use different priors for the analysis"
- Specific: "We specify prior distributions via the `priors_effect` argument"

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

- Do not use `library()` in package functions (only in vignettes is OK)
- Do not use absolute paths
- Do not commit temporary `.html` outputs
- Do not use `eval = TRUE` in chunk 3 unless regenerating models
- Do use relative paths (`../models/`, `../inst/`)
- Do compress models (`compress = "xz"`)
