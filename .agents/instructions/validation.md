# Validation: RoBMA specifics

The shared [validation policy](../../../.agents/instructions/validation.md) owns
the three layers, the evidence rules and the runtime expectations. This guide
adds only what is specific to this package's native code, and the runners are
documented in [testing](testing.md) and [scenarios](scenarios.md).

## R native error boundaries

R-facing C++ entry points use `ROBMA_NATIVE_BEGIN`/`ROBMA_NATIVE_END` and include
`r-native-api.h` after external headers. Individual allocating/error-capable R
calls need the protected adapter; wrapping an entire C++ computation in
`R_UnwindProtect` would still skip its destructors. Evaluate native expressions
before entering the adapter, and preserve R continuation conditions separately
from ordinary C++ exceptions. Keep destructors nonthrowing.

Workers consume prepared numeric pointers and never call error-capable R APIs.
The row executor materializes ALTREP views on the main thread after validation
and transports worker exceptions to the outer boundary. Shared JAGS kernels do
not include the R API remapping. Changes to this boundary require the focused
native-unwind tests as well as the relevant numerical checks.

## Performance investigation

### The Poisson GLMM AGHQ layout collision

`tools/bench-glmm-aghq.R` is the check to run after **any** native change (a
source edit under `src/`, a `Makevars` edit, a compiler or Rtools change, an
object added to or removed from `OBJECTS`):

```
Rscript tools/bench-glmm-aghq.R                       # all three fits, 5 reps
Rscript tools/bench-glmm-aghq.R --fits=fit_BMA --reps=3
Rscript tools/bench-glmm-aghq.R --root=<other tree>   # a comparison tree
```

**What it detects.** Historical builds of the Poisson AGHQ kernel `run_poisson`
(`src/glmm-aghq.cc`, reached through `.Call("RoBMA_glmm_pois_aghq", ...)`) differed
by about 2.2x in runtime despite a byte-identical compiled object and unchanged
computation and entry address. The statically linked mingw libm functions it
calls (`exp`, `log`, `lgamma`) moved relative to it; linked call and data-reference
operands therefore differed. These experiments establish a code-placement
effect. A specific CPU predictor/cache mechanism remains unproven.

The measured builds returned bit-identical numbers, so numerical tests,
snapshots and scenario fingerprints did not detect the slowdown. Every native
relink can change placement. Tested function alignment, page offsets and moving
`glmm-aghq.o` to the end of the link order did not reliably prevent the slowdown.

**Reading it.** The tool prints a median per cached fit. On the development
machine (AMD Ryzen 9 9950X, R 4.6.0 UCRT, Rtools45 g++ 14.2.0, one thread)
the historical `fit_BMA` runs clustered near 12.7 s and 28 s in fast and slow
layouts. Absolute seconds are machine-specific:
measure the tree you compare against in the same session, on the same machine,
and compare the medians.

**What to do with a slow reading.** Report the measured build and comparison.
Past experiments restored speed by changing code placement, but no padding or
link-order workaround is shipped. Keep the numerical implementation and settings unchanged;
any workaround needs evidence for that linked image and a maintainer decision.

**Note on `src/Makevars.win`.** R 4.6.0 is UCRT and reads `src/Makevars.ucrt`,
so `src/Makevars.win` does not affect the build on this machine (its
`-D_GLIBCXX_USE_CXX11_ABI=0` is therefore not applied here). Keep `OBJECTS`
consistent across all three files anyway; other R builds read the other two.
