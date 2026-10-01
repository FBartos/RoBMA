// Standard-normal tail kernel, instantiated once per lane.
//
// Included from selnorm-tail.h with SELNORM_TAIL_LANE naming the lane, inside a
// namespace of its own. The AVX2 copy is included inside that header's target
// region so the intrinsics inline there; both copies come from this one source,
// so the two lanes evaluate the same expressions in the same order.

typedef SELNORM_TAIL_LANE Lane;
typedef Lane::Value Value;
typedef Lane::Mask Mask;


// exp(value) for value in [-710, 1]; the Cephes rational form.
SELNORM_TAIL_INLINE Value exponential(Value value)
{
  const Value scaled = Lane::round_nearest(
    Lane::multiply(value, Lane::splat(exp_log2e)));
  Value reduced = Lane::fused(scaled, Lane::splat(-exp_ln2_high), value);
  reduced = Lane::fused(scaled, Lane::splat(-exp_ln2_low), reduced);
  const Value square = Lane::multiply(reduced, reduced);
  Value numerator = Lane::fused(Lane::splat(exp_p[0]), square,
                                Lane::splat(exp_p[1]));
  numerator = Lane::fused(numerator, square, Lane::splat(exp_p[2]));
  numerator = Lane::multiply(numerator, reduced);
  Value denominator = Lane::fused(Lane::splat(exp_q[0]), square,
                                  Lane::splat(exp_q[1]));
  denominator = Lane::fused(denominator, square, Lane::splat(exp_q[2]));
  denominator = Lane::fused(denominator, square, Lane::splat(exp_q[3]));
  const Value ratio = Lane::divide(numerator,
                                   Lane::subtract(denominator, numerator));
  const Value mantissa = Lane::fused(Lane::splat(2.0), ratio, Lane::splat(1.0));
  return Lane::multiply(mantissa, Lane::power_of_two(scaled));
}


// exp(-x^2/2) on Cody's split, so the squared argument's rounding error does
// not grow with x^2. Requires 0 <= x <= 37.6.
SELNORM_TAIL_INLINE Value gaussian_ordinate(Value x)
{
  const Value coarse = Lane::sixteenth(x);
  const Value difference = Lane::multiply(
    Lane::subtract(x, coarse), Lane::add(x, coarse));
  const Value coarse_exponent = Lane::multiply(
    Lane::multiply(coarse, coarse), Lane::splat(-0.5));
  const Value fine_exponent = Lane::multiply(difference, Lane::splat(-0.5));
  return Lane::multiply(exponential(coarse_exponent),
                        exponential(fine_exponent));
}


// R(x) = numerator(x) / (1 + x * denominator(x)); two independent Horner
// chains so the rational costs about one multiply-add of latency per degree.
SELNORM_TAIL_INLINE Value rational(
  Value x, const double *numerator_coefficients,
  int numerator_degree, const double *denominator_coefficients,
  int denominator_degree)
{
  Value numerator = Lane::splat(numerator_coefficients[numerator_degree]);
  for (int index = numerator_degree - 1; index >= 0; --index) {
    numerator = Lane::fused(numerator, x,
                            Lane::splat(numerator_coefficients[index]));
  }
  Value denominator =
    Lane::splat(denominator_coefficients[denominator_degree - 1]);
  for (int index = denominator_degree - 2; index >= 0; --index) {
    denominator = Lane::fused(denominator, x,
                              Lane::splat(denominator_coefficients[index]));
  }
  denominator = Lane::fused(denominator, x, Lane::splat(1.0));
  return Lane::divide(numerator, denominator);
}


// Q(x) = P(Z > x) for finite x.
SELNORM_TAIL_INLINE Value upper_tail(Value x)
{
  const Value magnitude = Lane::absolute(x);
  const Mask central = Lane::less_equal(magnitude, Lane::splat(small_cut));
  const Mask far = Lane::greater(magnitude, Lane::splat(middle_cut));

  Value result = Lane::splat(0.0);

  if (Lane::any(central)) {
    const Value square = Lane::multiply(x, x);
    result = Lane::select(central,
      Lane::subtract(Lane::splat(0.5), Lane::multiply(x,
        rational(square, central_numerator, central_degree,
                       central_denominator, central_degree))),
      result);
  }

  const Mask outer = Lane::negate(central);
  if (Lane::any(outer)) {
    // Both outer branches evaluate the small tail at |x| and then reflect.
    Value small = Lane::splat(0.0);
    const Mask middle_region = Lane::conjunction(outer, Lane::negate(far));
    if (Lane::any(middle_region)) {
      // Lanes outside this region are still evaluated; the rational stays
      // finite there and the blend below discards the value.
      const Value guarded = Lane::select(middle_region, magnitude,
                                         Lane::splat(middle_cut));
      small = Lane::select(middle_region,
        Lane::multiply(gaussian_ordinate(guarded),
          rational(guarded, middle_numerator, middle_degree,
                         middle_denominator, middle_degree)),
        small);
    }
    if (Lane::any(far)) {
      // Beyond the tail cut the tail underflows the double range; the
      // established algorithms return an exact zero there. The guarded
      // magnitude keeps the unused lanes inside the fitted range.
      const Mask usable = Lane::conjunction(
        far, Lane::less_equal(magnitude, Lane::splat(tail_cut)));
      const Value guarded = Lane::select(usable, magnitude,
                                         Lane::splat(middle_cut + 1.0));
      const Value inverse_square = Lane::divide(Lane::splat(1.0),
        Lane::multiply(guarded, guarded));
      const Value series = Lane::divide(
        Lane::subtract(Lane::splat(inverse_root_two_pi),
          Lane::multiply(inverse_square,
            rational(inverse_square, far_numerator, far_degree,
                           far_denominator, far_degree))),
        guarded);
      small = Lane::select(usable,
        Lane::multiply(gaussian_ordinate(guarded), series), small);
      small = Lane::select(Lane::conjunction(far, Lane::negate(usable)),
                           Lane::splat(0.0), small);
    }
    const Mask reflect = Lane::conjunction(outer,
      Lane::less_equal(x, Lane::splat(0.0)));
    result = Lane::select(outer,
      Lane::select(reflect, Lane::subtract(Lane::splat(1.0), small), small),
      result);
  }

  return result;
}


