# Remove a fitted model's selection cache

Removes retained selection-normalizer cache entries from a fitted object
without changing its draws, inference, or ability to obtain more draws.

## Usage

``` r
remove_selection_cache(object)
```

## Arguments

- object:

  a fitted `brma` object, including
  [`bselmodel.mv()`](https://fbartos.github.io/RoBMA/dev/reference/bselmodel.mv.md)
  and
  [`RoBMA()`](https://fbartos.github.io/RoBMA/dev/reference/RoBMA.md)
  fits.

## Value

The fitted object with its retained selection cache removed.

## Details

Set `RoBMA.options(selection.cache_retain = TRUE)` before fitting or
extending to retain exact and coarse cache entries. The default is
`FALSE`. This option is checked for each RoBMA fit or
[`update()`](https://rdrr.io/r/stats/update.html) call. Existing
retained entries are used during extension even when retention is turned
off; in that case the returned fit does not retain a new cache.

Retained entries are ordinary serialized data and survive
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html). They are restored
only with matching package code, native builds, R and JAGS. Incompatible
or damaged cache data are discarded with a warning; extension proceeds
with ordinary likelihood evaluation. Cache contents never replace the
model, numerical integration controls, or JAGS random-number state.
Reconstructed JAGS models can still require adaptation before collecting
new samples; retaining the normalizer cache does not retain sampler
tuning.

Each worker saves its own entries. Extension with the same worker count
restores the corresponding worker cache. With a different worker count,
previous worker caches are assigned round-robin without duplication, and
entries that do not fit the destination's share of the memory budget are
dropped. Increasing the worker count can leave some workers with an
empty cache initially. Cache reuse can improve speed but is not
guaranteed.

Successful capture releases the native cache. A retained object's
compact cache payload is bounded by the fitted total native budget,
apart from small R metadata overhead. The byte budget is not a limit on
total R memory: capture, serialization and extension can temporarily
hold both a snapshot and native storage, and other copies or saved fits
can retain their own snapshots. Save or load operations can therefore
require substantial memory.

Assign the result back to the object to release its reference to the
retained data. Memory remains in use if another object references the
same snapshot. This helper does not clear the current process's
independent native cache; use `selection_cache_info(clear = TRUE)` for
that. Future retention is controlled by `selection.cache_retain`.

## Examples

``` r
if (FALSE) { # \dontrun{
RoBMA.options(selection.cache_retain = TRUE)
fit <- bselmodel.mv(yi, V = V, random = ~ 1 | study / esid, data = dat)
saveRDS(fit, "fit.rds")
fit <- update(readRDS("fit.rds"), sample_extend = 1000)
fit <- remove_selection_cache(fit)
RoBMA.options(selection.cache_retain = FALSE)
} # }
```
