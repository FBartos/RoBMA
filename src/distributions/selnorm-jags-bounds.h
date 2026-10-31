#ifndef SELNORM_JAGS_BOUNDS_H_
#define SELNORM_JAGS_BOUNDS_H_

#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

// Validate controls before any floating-to-integer conversion.
inline bool selnorm_jags_integer_in_range(double value, int minimum, int maximum)
{
  return std::isfinite(value) && value >= minimum && value <= maximum &&
    value == std::floor(value);
}

// All operands are positive counts. Bound multiplication by the actual array
// length, which the distribution's length gate keeps within native int indexing.
inline bool selnorm_jags_qmc_length(unsigned int length, unsigned int dimensions,
    unsigned int points, unsigned int scrambles)
{
  if (dimensions == 0 || points == 0 || scrambles == 0 ||
      points > length / dimensions) return false;
  const unsigned int partial = dimensions * points;
  return length % partial == 0 && length / partial == scrambles;
}

// JAGS serializes infinite selection bounds as exact finite sentinels.
// Decode them only at distribution ingress; observations remain finite.
inline double selnorm_jags_bound_to_ieee(double value)
{
  if (value == -1e300) {
    return -std::numeric_limits<double>::infinity();
  }
  if (value == 1e300) {
    return std::numeric_limits<double>::infinity();
  }
  return value;
}

class SelNormJagsBounds
{
  enum { stack_capacity = 64 };

public:
  SelNormJagsBounds(const double *values, unsigned int length)
    : heap_(length > stack_capacity ? length : 0),
      decoded_(length > stack_capacity ? heap_.data() : stack_)
  {
    for (unsigned int i = 0; i < length; ++i) {
      decoded_[i] = selnorm_jags_bound_to_ieee(values[i]);
    }
  }

  const double *data() const
  {
    return decoded_;
  }

  double operator[](unsigned int index) const
  {
    return decoded_[index];
  }

private:
  double stack_[stack_capacity];
  std::vector<double> heap_;
  double *decoded_;

  SelNormJagsBounds(const SelNormJagsBounds &);
  SelNormJagsBounds &operator=(const SelNormJagsBounds &);
};

#endif
