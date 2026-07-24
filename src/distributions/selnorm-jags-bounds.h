#ifndef SELNORM_JAGS_BOUNDS_H_
#define SELNORM_JAGS_BOUNDS_H_

#include <limits>
#include <vector>

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
