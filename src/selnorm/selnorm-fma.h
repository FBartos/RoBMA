#ifndef ROBMA_SELNORM_FMA_H_
#define ROBMA_SELNORM_FMA_H_

#include <cmath>

// MinGW's software fma is expensive inside the selected-normal quadrature.
// Use hardware only when GCC reports CPU and OS AVX/FMA support; other
// platforms retain their standard-library implementation.
#if defined(__MINGW32__) && defined(__GNUC__) && defined(__x86_64__)
namespace selnorm_fma_detail {
using Operation = double (*)(double, double, double);

__attribute__((target("fma"), noinline))
inline double hardware(double a, double b, double c)
{
  return std::fma(a, b, c);
}

inline double standard(double a, double b, double c)
{
  return std::fma(a, b, c);
}

inline bool supported()
{
  return __builtin_cpu_supports("avx") && __builtin_cpu_supports("fma");
}
}

inline double selnorm_fma(double a, double b, double c)
{
  static const selnorm_fma_detail::Operation operation =
    selnorm_fma_detail::supported() ? selnorm_fma_detail::hardware :
                                    selnorm_fma_detail::standard;
  return operation(a, b, c);
}
#else
inline double selnorm_fma(double a, double b, double c)
{
  return std::fma(a, b, c);
}
#endif

#endif
