# Maintainer scenarios

Scenario files are deliberately excluded from routine package tests. Name each
file `test-<scenario>.R` and run it with:

```text
Rscript tools/test-scenario.R <scenario>
```

Use `--regenerate` to refit every wrapped model and replace every wrapped text
and plot baseline exercised by the script. Review all resulting Git diffs.

Minimal scenario:

```r
REGENERATE_SCENARIO_FILES <- FALSE

scenario_start("bcg")

testthat::test_that("BCG analysis remains stable", {

  fit <- scenario_fit("meta-analysis", {
    brma(y = yi, se = sqrt(vi), data = dat)
  })

  scenario_text("summary", {
    print(summary(fit))
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
that call changes. In an interactive session, `scenario_text()` also echoes its
captured output. `scenario_plot()` draws on the active graphics device; outside
a scenario test runner it stops there, while runner-backed calls also perform
the vdiffr comparison. Base graphics parameters are restored after each
interactive evaluation, so settings such as `par(mar = ...)` do not leak into
subsequent plots. Base graphics overlay mode is reset before and after each
interactive block, preventing one scenario figure from being drawn over the
previous figure.
