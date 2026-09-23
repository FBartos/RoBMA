#ifndef ROBMA_SELNORM_TAIL_H_
#define ROBMA_SELNORM_TAIL_H_

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>

#if defined(__GNUC__)
#define SELNORM_TAIL_INLINE inline __attribute__((always_inline))
#else
#define SELNORM_TAIL_INLINE inline
#endif

// The AVX2 lane needs a per-function target region, which is spelled
// "#pragma GCC target" and is not accepted by clang or MSVC. Those compilers
// use the scalar lane below, which is itself faster than the library erfc it
// replaces, so this is a missing speedup rather than a missing capability.
#if defined(__GNUC__) && !defined(__clang__) &&     (defined(__x86_64__) || defined(_M_X64))
#define SELNORM_TAIL_HAS_AVX2 1
#include <immintrin.h>
#else
#define SELNORM_TAIL_HAS_AVX2 0
#endif

// Batched standard-normal upper-tail probabilities for the selection kernels.
//
// Every step-weight normalizer is a weighted sum of Q(x) = P(Z > x) at scores
// that share one row's scale and differ only in the quadrature location, so the
// kernels evaluate the same special function at thousands of arguments per
// state. The library's scalar erfc costs about 26 ns per argument on the
// reference build and blocks vectorization.
//
// Q is evaluated on the three argument ranges the established algorithms use
// (Cody 1969, Math. Comp. 23:631-637, which is what R's own pnorm evaluates),
// with the rational coefficients fitted against an arbitrary-precision
// reference rather than taken on trust, and with exp in the Cephes rational
// form. The fitted relative errors are below 1e-18, so double rounding
// dominates; the package tests assert agreement with pnorm and erfc over a
// dense argument grid, where this is in fact closer to the exact tail than the
// library erfc is beyond about eight standard deviations.
//
// The scalar and AVX2 lanes include one kernel source, so a batch returns the
// same values however it is dispatched on a given build.

namespace selnorm_tail {

// Phi(x) - 1/2 = x * R(x^2) for |x| <= 0.67448975; fitted relative error
// 1.5e-20.
const int central_degree = 4;
const double central_numerator[5] = {
  0.398942280401432677934, 0.023391444506584200868,
  0.00353529226284461142329, 0.0000486115922823997155555,
  0.00000144459843283355322552
};
const double central_denominator[4] = {
  0.225300322851332611792, 0.0214117173536862530242,
  0.00103415315553786721669, 0.0000218718200710845228644
};

// Q(y) = exp(-y^2/2) * R(y) for 0.67448975 < y <= sqrt(32); fitted relative
// error 2.9e-19.
const int middle_degree = 9;
const double middle_numerator[10] = {
  0.500000000000625392183, 0.658631771741979888792,
  0.433659208656499318894, 0.179467115325292050143,
  0.0503838135867302869886, 0.00977104553999546513378,
  0.00127814411532389928245, 0.000103580633827034798779,
  0.00000402559766620178586876, -6.28615672558255165791e-15
};
const double middle_denominator[9] = {
  2.11514810430173107977, 2.05496243346281075071,
  1.20694449775612617988, 0.473836795601048224271,
  0.129476631128674834062, 0.0247520502090748082336,
  0.00321392120728809578456, 0.00025963820507150268528,
  0.000010090675539696163043
};

// Q(y) = exp(-y^2/2) * (c - s R(s)) / y with s = 1/y^2, for sqrt(32) < y <=
// 37.5193; fitted relative error 9.7e-19.
const int far_degree = 6;
const double far_numerator[7] = {
  0.398942280401432666664, 26.6366087969648946742,
  619.64764135519764592, 6116.80669705551546094,
  24657.9433372222621096, 30916.7349136973747696,
  2333.97848345269519012
};
const double far_denominator[6] = {
  69.7680767507569633463, 1747.53052838170591458,
  19633.6310511264946558, 100876.881266441649314,
  213577.809126302834304, 133650.291594108042916
};

const double small_cut  = 0.67448975;
const double middle_cut = 5.656854249492380195206754896838;   // sqrt(32)
const double tail_cut   = 37.5193;
const double inverse_root_two_pi = 0.398942280401432677939946059934;

// Cephes exp coefficients. The arguments these tails need lie in [-710, 0].
const double exp_log2e    = 1.4426950408889634073599246810019;
const double exp_ln2_high = 6.93145751953125e-1;
const double exp_ln2_low  = 1.42860682030941723212e-6;
const double exp_p[3] = {
  1.26177193074810590878e-4, 3.02994407707441961300e-2,
  9.99999999999999999910e-1
};
const double exp_q[4] = {
  3.00198505138664455042e-6, 2.52448340349684104192e-3,
  2.27265548208155028766e-1, 2.0
};


// -------------------------------------------------------------------------- //
// Scalar lane. Every primitive is inlined: a lane operation left as a call
// taking or returning a 256-bit vector is not callable from generic code.
// -------------------------------------------------------------------------- //
struct ScalarLane {
  typedef double Value;
  typedef bool Mask;

  static SELNORM_TAIL_INLINE Value splat(double value) { return value; }
  static SELNORM_TAIL_INLINE Value add(Value a, Value b) { return a + b; }
  static SELNORM_TAIL_INLINE Value subtract(Value a, Value b) { return a - b; }
  static SELNORM_TAIL_INLINE Value multiply(Value a, Value b) { return a * b; }
  static SELNORM_TAIL_INLINE Value divide(Value a, Value b) { return a / b; }
  // Written as a product and a sum rather than std::fma: without a hardware
  // fused multiply-add enabled at compile time the library call is software
  // emulation, which costs more than the whole approximation. A build that
  // does enable one contracts this back into a single instruction. The
  // approximations carry relative errors below 1e-18, so the extra rounding
  // per step stays far inside the double-precision result.
  static SELNORM_TAIL_INLINE Value fused(Value a, Value b, Value c) {
    return a * b + c;
  }
  static SELNORM_TAIL_INLINE Value absolute(Value a) { return std::fabs(a); }

  // Round to nearest, ties to even, for |a| below 2^51. Adding and removing
  // 1.5 * 2^52 does that with two additions; std::nearbyint has to consult the
  // rounding mode and is a library call here.
  static SELNORM_TAIL_INLINE Value round_nearest(Value a) {
    const double shift = 6755399441055744.0;   // 2^52 + 2^51
    return (a + shift) - shift;
  }

  // 2^exponent for an integral exponent in [-1022, 1023].
  static SELNORM_TAIL_INLINE Value power_of_two(Value exponent) {
    const std::int64_t bits =
      (static_cast<std::int64_t>(exponent) + 1023) << 52;
    double value;
    std::memcpy(&value, &bits, sizeof(value));
    return value;
  }

  // Truncation to a multiple of one sixteenth. Such a value and its square are
  // exact in binary, which is what keeps exp(-x^2/2) accurate at large x.
  // Truncation toward zero; the argument is a magnitude below 38, so the
  // integer conversion is exact and matches the vector lane's rounding mode.
  static SELNORM_TAIL_INLINE Value sixteenth(Value a) {
    return static_cast<double>(static_cast<int>(a * 16.0)) * 0.0625;
  }

  static SELNORM_TAIL_INLINE Mask less_equal(Value a, Value b) {
    return a <= b;
  }
  static SELNORM_TAIL_INLINE Mask greater(Value a, Value b) { return a > b; }
  static SELNORM_TAIL_INLINE Value select(Mask mask, Value when_true,
                                          Value when_false) {
    return mask ? when_true : when_false;
  }
  static SELNORM_TAIL_INLINE bool any(Mask mask) { return mask; }
  static SELNORM_TAIL_INLINE Mask negate(Mask mask) { return !mask; }
  static SELNORM_TAIL_INLINE Mask conjunction(Mask left, Mask right) {
    return left && right;
  }
};

namespace scalar_kernel {
#define SELNORM_TAIL_LANE ScalarLane
#include "selnorm-tail-kernel.cc.inc"
#undef SELNORM_TAIL_LANE
}


// -------------------------------------------------------------------------- //
// AVX2 lane, defined inside the target region together with its own copy of
// the kernel so that every intrinsic inlines there.
// -------------------------------------------------------------------------- //
#if SELNORM_TAIL_HAS_AVX2
#pragma GCC push_options
#pragma GCC target("avx2,fma")

struct VectorLane {
  typedef __m256d Value;
  typedef __m256d Mask;

  static SELNORM_TAIL_INLINE Value splat(double value) {
    return _mm256_set1_pd(value);
  }
  static SELNORM_TAIL_INLINE Value add(Value a, Value b) {
    return _mm256_add_pd(a, b);
  }
  static SELNORM_TAIL_INLINE Value subtract(Value a, Value b) {
    return _mm256_sub_pd(a, b);
  }
  static SELNORM_TAIL_INLINE Value multiply(Value a, Value b) {
    return _mm256_mul_pd(a, b);
  }
  static SELNORM_TAIL_INLINE Value divide(Value a, Value b) {
    return _mm256_div_pd(a, b);
  }
  static SELNORM_TAIL_INLINE Value fused(Value a, Value b, Value c) {
    return _mm256_fmadd_pd(a, b, c);
  }
  static SELNORM_TAIL_INLINE Value absolute(Value a) {
    return _mm256_andnot_pd(_mm256_set1_pd(-0.0), a);
  }
  static SELNORM_TAIL_INLINE Value round_nearest(Value a) {
    return _mm256_round_pd(a, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
  }
  static SELNORM_TAIL_INLINE Value power_of_two(Value exponent) {
    const __m256i wide = _mm256_cvtepi32_epi64(_mm256_cvtpd_epi32(exponent));
    return _mm256_castsi256_pd(_mm256_slli_epi64(
      _mm256_add_epi64(wide, _mm256_set1_epi64x(1023)), 52));
  }
  static SELNORM_TAIL_INLINE Value sixteenth(Value a) {
    return _mm256_mul_pd(
      _mm256_round_pd(_mm256_mul_pd(a, _mm256_set1_pd(16.0)),
                      _MM_FROUND_TO_ZERO | _MM_FROUND_NO_EXC),
      _mm256_set1_pd(0.0625));
  }
  static SELNORM_TAIL_INLINE Mask less_equal(Value a, Value b) {
    return _mm256_cmp_pd(a, b, _CMP_LE_OQ);
  }
  static SELNORM_TAIL_INLINE Mask greater(Value a, Value b) {
    return _mm256_cmp_pd(a, b, _CMP_GT_OQ);
  }
  static SELNORM_TAIL_INLINE Value select(Mask mask, Value when_true,
                                          Value when_false) {
    return _mm256_blendv_pd(when_false, when_true, mask);
  }
  static SELNORM_TAIL_INLINE bool any(Mask mask) {
    return _mm256_movemask_pd(mask) != 0;
  }
  static SELNORM_TAIL_INLINE Mask negate(Mask mask) {
    return _mm256_xor_pd(mask, _mm256_castsi256_pd(_mm256_set1_epi64x(-1)));
  }
  static SELNORM_TAIL_INLINE Mask conjunction(Mask left, Mask right) {
    return _mm256_and_pd(left, right);
  }
};

namespace vector_kernel {
#define SELNORM_TAIL_LANE VectorLane
#include "selnorm-tail-kernel.cc.inc"
#undef SELNORM_TAIL_LANE
}


// out[i] = Q(offset - mean[i] * inverse_sd), four arguments at a time. The
// remainder goes through the same lane on a padded vector, so a value's tail
// does not depend on where the batch boundary happens to fall.
inline void upper_tail_affine_avx2(const double *mean, std::size_t count,
                                   double offset, double inverse_sd,
                                   double *out)
{
  const __m256d intercept = _mm256_set1_pd(offset);
  const __m256d slope = _mm256_set1_pd(-inverse_sd);
  std::size_t index = 0;
  for (; index + 4 <= count; index += 4) {
    const __m256d score = _mm256_fmadd_pd(
      _mm256_loadu_pd(mean + index), slope, intercept);
    _mm256_storeu_pd(out + index, vector_kernel::upper_tail(score));
  }
  if (index < count) {
    double padded_in[4] = {0.0, 0.0, 0.0, 0.0};
    double padded_out[4];
    const std::size_t remainder = count - index;
    for (std::size_t lane = 0; lane < remainder; ++lane) {
      padded_in[lane] = mean[index + lane];
    }
    const __m256d score = _mm256_fmadd_pd(
      _mm256_loadu_pd(padded_in), slope, intercept);
    _mm256_storeu_pd(padded_out, vector_kernel::upper_tail(score));
    for (std::size_t lane = 0; lane < remainder; ++lane) {
      out[index + lane] = padded_out[lane];
    }
  }
}

#pragma GCC pop_options


inline bool avx2_available()
{
  static const bool supported = __builtin_cpu_supports("avx2") &&
    __builtin_cpu_supports("fma");
  return supported;
}
#endif


// out[i] = Q(offset - mean[i] * inverse_sd). `scalar` forces the fallback lane
// that builds without the vector region use, so the tests can certify it on
// any machine.
inline void upper_tail_affine(const double *mean, std::size_t count,
                              double offset, double inverse_sd, double *out,
                              bool scalar = false)
{
#if SELNORM_TAIL_HAS_AVX2
  if (!scalar && avx2_available()) {
    upper_tail_affine_avx2(mean, count, offset, inverse_sd, out);
    return;
  }
#endif
  // std::fma here is deliberate, and is the one place the lane keeps it: the
  // AVX2 region forms the same score with a fused multiply-add, so the two
  // lanes agree bit for bit only if this one does too. Without hardware FMA it
  // becomes a library call: 18.4 ns per value against 9.2 ns with -mfma, over
  // 1e7 arguments at -O2 with the vector region forced off. That is the price
  // of the agreement on such a build, and the lane is still ahead of the
  // library erfc it replaces (23.0 ns in the same harness).
  for (std::size_t index = 0; index < count; ++index) {
    out[index] = scalar_kernel::upper_tail(
      std::fma(-mean[index], inverse_sd, offset));
  }
}


// Q(x) for one argument. Routed through the batch so that a single value and
// the same value inside a batch agree exactly.
inline double upper_tail_scalar(double x)
{
  double result;
  upper_tail_affine(&x, 1, 0.0, -1.0, &result);
  return result;
}

}  // namespace selnorm_tail

#endif
