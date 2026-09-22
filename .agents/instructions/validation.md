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

**What it detects.** The Poisson AGHQ kernel `run_poisson`
(`src/glmm-aghq.cc`, reached through `.Call("RoBMA_glmm_pois_aghq", ...)`) runs
2.2x slower in some linked images than in others, with its compiled object
byte-identical and its computation and entry address unchanged. Linked call
and data-reference operands differ because their targets are relocated.
What moves is the placement of
the statically linked mingw libm bodies it calls (`exp`, `log`, `lgamma`)
relative to it: an exact placement collides (branch-target aliasing is the
consistent hypothesis; unproven, and naming the resource needs hardware
counters that are not installed here). Any native change relinks the image and
can land on such a placement; the measured chance is a few percent.

The collision is invisible to every other check. All repetitions return
bit-identical numbers, so no test, no snapshot and no scenario fingerprint
moves; only `add_loo()` on a Poisson GLMM gets slower, and the GLMM scenario
timings are the only place it would eventually surface. Neither 64-byte
function alignment, nor the page offset of `exp`, nor a link order that pins
`glmm-aghq.o` last prevents it (pinning it last was itself slow on both trees
measured).

**Reading it.** The tool prints a median per cached fit. On the development
machine (AMD Ryzen 9 9950X, R 4.6.0 UCRT, Rtools45 g++ 14.2.0, one thread)
`fit_BMA` takes about 12.7 s in a fast layout and about 28 s in the collided
one; there is no intermediate value. Absolute seconds are machine-specific:
measure the tree you compare against in the same session, on the same machine,
and compare the medians.

**What to do with a slow reading.** Report it - it is a build-layout property,
not a regression in the change under test, and the maintainer decides whether
to carry a workaround. The available results-identical mitigation is a
never-called pad object that shifts the image tail: a `used`, `noinline`
function whose body is `.skip N, 0x90` in its own translation unit, listed in
`OBJECTS` of `src/Makevars.in`, `src/Makevars.win` and `src/Makevars.ucrt`.
Any `N >= 0x100` restored the fast layout in all 15 padded builds measured, but
a pad is a property of one image: a later native change relinks and re-rolls
the placement, so re-run the benchmark rather than trusting the pad.

**Note on `src/Makevars.win`.** R 4.6.0 is UCRT and reads `src/Makevars.ucrt`,
so `src/Makevars.win` does not affect the build on this machine (its
`-D_GLIBCXX_USE_CXX11_ABI=0` is therefore not applied here). Keep `OBJECTS`
consistent across all three files anyway; other R builds read the other two.
