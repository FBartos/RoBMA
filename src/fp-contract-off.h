#ifndef ROBMA_FP_CONTRACT_OFF_H_
#define ROBMA_FP_CONTRACT_OFF_H_

// Included first by a translation unit whose arithmetic decides seeded draws.
//
// A compiler may contract a * b + c into one fused multiply-add wherever the
// target has the instruction: always on aarch64, and on x86-64 from
// -march=x86-64-v3 (the GCC default on Ubuntu 26.04). The fused result keeps
// the rounding error of the product, so such a build accepts and rejects
// proposals differently in the last bits and its chains leave those of the
// other builds: the same seed then gives different draws on different
// platforms (equally valid, but not reproducible). With contraction off, every
// build evaluates the same double operations. The pragmas act from here to the
// end of the translation unit, on the headers included after this one too.
#if defined(__clang__)
#pragma STDC FP_CONTRACT OFF
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

#endif
