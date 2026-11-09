#ifndef ROBMA_SELNORM_NODE_PRODUCT_H_
#define ROBMA_SELNORM_NODE_PRODUCT_H_

#include <cfloat>
#include <cmath>
#include <cstdint>

// Products of the row normalizers at a quadrature node reach omega^K, far below
// the range of double once one latent axis owns many rows. The nested rule of a
// forest has no log-scale fallback at the orders where the tensor rule is over
// the node budget, so its products need that whole range.
//
// Where long double has more range than double (x87 on x86-64, IEEE quad), a
// SelNormNodeProduct is a long double and these helpers are its plain checks, so
// the results are bit for bit those of plain long double arithmetic. Where it
// has the range of double (Apple silicon), the type carries a binary exponent
// per node, and only the mantissa rounds: a product has the relative accuracy
// of double arithmetic at every scale.
//
// The other direct rules keep long double products, because a log-scale fallback
// follows them: where long double has only the range of double, a product that
// reaches the subnormal range has lost bits and goes to that fallback.

#if LDBL_MIN_EXP < DBL_MIN_EXP

#define SELNORM_LONG_DOUBLE_HAS_RANGE 1

typedef long double SelNormNodeProduct;

inline bool selnorm_direct_product_usable(long double product)
{
  return product > 0.0L && std::isfinite(product);
}

// A term of a direct sum: on this range an exhausted term has no mass to lose.
inline bool selnorm_direct_term_usable(long double)
{
  return true;
}

inline bool selnorm_node_product_usable(long double product)
{
  return selnorm_direct_product_usable(product);
}

inline double selnorm_node_product_log(long double product)
{
  return static_cast<double>(std::log(product));
}

#else

#define SELNORM_LONG_DOUBLE_HAS_RANGE 0

inline bool selnorm_direct_product_usable(long double product)
{
  return product >= LDBL_MIN && std::isfinite(product);
}

// A subnormal term has lost bits and may grow again under a weight above one.
inline bool selnorm_direct_term_usable(long double term)
{
  return term >= LDBL_MIN;
}

// mantissa * 2^exponent, signed (the moments sum signed nodes). The mantissa
// stays within [2^-256, 2^256], so products of two mantissas and aligned sums
// neither overflow nor reach the subnormal range; zero and non-finite values
// are carried in the mantissa.
class SelNormScaledProduct {
 public:
  SelNormScaledProduct(long double value = 0.0L)
    : mantissa_(static_cast<double>(value)), exponent_(0)
  {
    rescale();
  }

  SelNormScaledProduct &operator*=(const SelNormScaledProduct &other)
  {
    mantissa_ *= other.mantissa_;
    exponent_ += other.exponent_;
    rescale();
    return *this;
  }

  SelNormScaledProduct &operator+=(const SelNormScaledProduct &other)
  {
    if (other.mantissa_ == 0.0) return *this;
    if (mantissa_ == 0.0) return *this = other;
    // Align on the larger exponent. What ldexp then loses is below 2^-800 of
    // the larger term, far under the rounding of the sum.
    const std::int64_t gap = exponent_ - other.exponent_;
    if (gap >= 0) {
      mantissa_ += std::ldexp(other.mantissa_, -bounded(gap));
    } else {
      mantissa_ = std::ldexp(mantissa_, -bounded(-gap)) + other.mantissa_;
      exponent_ = other.exponent_;
    }
    rescale();
    return *this;
  }

  friend SelNormScaledProduct operator*(SelNormScaledProduct left,
                                        const SelNormScaledProduct &right)
  {
    left *= right;
    return left;
  }

  // Division by a plain number keeps the scale.
  friend SelNormScaledProduct operator/(SelNormScaledProduct left, long double right)
  {
    left.mantissa_ /= static_cast<double>(right);
    left.rescale();
    return left;
  }

  // A ratio of two scaled values is a plain number (probabilities and moments
  // of the rule); a ratio below the range of double is zero mass.
  friend long double operator/(const SelNormScaledProduct &left,
                               const SelNormScaledProduct &right)
  {
    const std::int64_t gap = left.exponent_ - right.exponent_;
    return std::ldexp(left.mantissa_ / right.mantissa_,
                      gap >= 0 ? bounded(gap) : -bounded(-gap));
  }

  bool usable() const { return mantissa_ > 0.0 && std::isfinite(mantissa_); }

  double log() const
  {
    return std::log(mantissa_) + static_cast<double>(exponent_) * 0.693147180559945309417232121458;
  }

 private:
  static int bounded(std::int64_t gap)
  {
    return gap > 4096 ? 4096 : static_cast<int>(gap);
  }

  void rescale()
  {
    const double magnitude = std::fabs(mantissa_);
    if ((magnitude >= 0x1p-256 && magnitude <= 0x1p256) ||
        mantissa_ == 0.0 || !std::isfinite(mantissa_)) {
      return;
    }
    int shift = 0;
    mantissa_ = std::frexp(mantissa_, &shift);
    exponent_ += shift;
  }

  double mantissa_;
  std::int64_t exponent_;
};

typedef SelNormScaledProduct SelNormNodeProduct;

inline bool selnorm_node_product_usable(const SelNormScaledProduct &product)
{
  return product.usable();
}

inline double selnorm_node_product_log(const SelNormScaledProduct &product)
{
  return product.log();
}

#endif

#endif
