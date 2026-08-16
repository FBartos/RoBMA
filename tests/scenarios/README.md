# Maintainer scenarios

Scenario files are deliberately excluded from routine package tests. Name each
file `test-<scenario>.R` and run it with:

```text
Rscript tools/test-scenario.R <scenario>
```

Use `--regenerate` to refit every wrapped model and replace every wrapped text
and plot baseline exercised by the script. Review all resulting Git diffs.
Use `scenario_agreement_plot()` for reference-versus-RoBMA difference plots so
their finite-value filtering, agreement band, and axes stay consistent.

Minimal scenario:

```r
REGENERATE_SCENARIO_FILES <- FALSE
SHOW_SCENARIO_OUTPUT      <- FALSE

scenario_start("bcg")

testthat::test_that("BCG analysis remains stable", {

  fit <- scenario_fit("meta-analysis", {
    brma(y = yi, se = sqrt(vi), data = dat)
  })

  scenario_text("summary", {
    summary(fit)
  })

  scenario_plot("forest", {
    forest(fit)
  })
})
```

Fit caches under `cache/` are local and ignored. Text files under `results/`
and SVG files under `_snaps/` are regression baselines and must be committed.
Missing baselines are added on a local first run and rejected on CI.
Scenario and artifact names may use ASCII letters in either case, numbers,
underscores, hyphens, and internal periods.

Each fit cache records its normalized fitting call and refits automatically when
that call changes. `scenario_text()` automatically prints a visible returned
value, so summary and table calls do not need an explicit `print()`. Set
`SHOW_SCENARIO_OUTPUT <- TRUE` while developing so captured text is echoed and
`scenario_plot()` draws on the active graphics device; leave it `FALSE` for
quiet `test_file()` / `tools/test-scenario.R` runs. Outside a scenario test
runner, `scenario_plot()` still draws when output is enabled (or when the
session is interactive) and stops there; runner-backed calls also perform the
vdiffr comparison. Base graphics parameters are restored after each drawn
evaluation, so settings such as `par(mar = ...)` do not leak into subsequent
plots. Base graphics overlay mode is reset before and after each drawn block,
preventing one scenario figure from being drawn over the previous figure.
