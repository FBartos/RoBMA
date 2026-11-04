# Inspect or clear the native selection cache

Reports the bounded native selection-normalizer cache in the current R
process. Exact caching reuses identical numerical inputs and does not
change the likelihood or its integration diagnostics. Post-fit requests
using bounded CDF interpolation bypass these caches, including when they
fall back to the ordinary primitive. Fitting retains its
ordinary-primitive cache.

## Usage

``` r
selection_cache_info(clear = FALSE, component = c("all", "exact", "coarse"))
```

## Arguments

- clear:

  clear the requested entries and reset their counters before reporting
  them. Clearing `"all"` also releases the shared storage.

- component:

  report the `"exact"` cache, the `"coarse"` proposal cache, or `"all"`
  caches.

## Value

A data frame with the cache component, byte counts, entry count, and
lookup, eviction and allocation counters. The `"total"` row reports the
shared capacity and directly allocated blocks, including retained
controls and index buckets; component byte counts are `NA` because
storage is shared. Allocator metadata, the surrounding R process and
JAGS model are excluded; this is not a process-RSS limit. Eviction
releases entry blocks, while controls and index buckets can remain
allocated. Oversize lookup keys increment misses without incrementing
allocation failures. Failed insertion, capacity or allocation attempts
can increment allocation failures; these counts do not necessarily
indicate exhaustion of system memory.

## Details

The default `selection.cache_max_bytes = "auto"` chooses the smaller of
4 GiB and one quarter of currently available system RAM. Available RAM
is reported by
[`ps::ps_system_memory()`](https://ps.r-lib.org/reference/ps_system_memory.html),
excluding swap. A numeric setting supplies an explicit total byte
budget; `0` disables storage. If the available-memory query fails, a
warning is issued and caching stays disabled until a budget can be
resolved or supplied explicitly. Loading the package does not query
available memory. The first selection computation or cache inspection
initializes the current process's cache, including computations on a
saved fit. Explicit option changes and fitting setup also initialize it.
Exact and coarse entries share one byte budget for directly allocated
key, control, list-node and index blocks, with global
least-recently-used eviction. Storage is allocated lazily. The ordinary
sampler does not use coarse entries.

The budget is resolved once before fitting and divided equally among
chains. Shares are pooled for chains assigned to the same process: ten
parallel workers receive one tenth each, while ten sequential chains
share the total. Integer-byte rounding can leave a few bytes unused. The
coordinator releases its cache while workers run and restores its local
budget when they stop. Separate simultaneous fits are separate budgets;
this is not a machine-wide memory manager. Stored fit extensions retain
the resolved budget and sampler settings, redistributing the budget for
their current worker topology.

This function reports only the current process, not separate workers or
saved fitted-object caches. By default entries live in the native
process. Set `RoBMA.options(selection.cache_retain = TRUE)` to retain
them with the fit; see
[`remove_selection_cache`](https://fbartos.github.io/RoBMA/reference/remove_selection_cache.md).
Changing the budget or clearing all entries releases the shared storage.
Clearing one component releases its keys, nodes and controls; shared
index buckets can remain allocated. Switching sampler options does not
clear existing entries or replace already compiled samplers.

The exact cache accelerates successful deterministic covariance-envelope
normalizers for multivariate step selection. Unsupported numerical
routes continue through their ordinary evaluator. The coarse cache can
also store its prescribed unit-normalizer fallback. If cache storage
cannot be allocated, the same exact or surrogate evaluator runs without
storage; memory capacity never determines the returned normalizer or its
numerical diagnostics.
