bool selnorm_bound_is_negative_infinity(double x)
{
  return x <= -1e290 || (std::isinf(x) && x < 0);
}

bool selnorm_bound_is_positive_infinity(double x)
{
  return x >= 1e290 || (std::isinf(x) && x > 0);
}

bool selnorm_same_finite_bound(double x, double y)
{
  return std::isfinite(x) && std::isfinite(y) && x == y;
}

bool selnorm_is_descending_step_partition(const double *z_lower,
                                          const double *z_upper,
                                          int n_bins)
{
  if (n_bins < 1 ||
      !selnorm_bound_is_positive_infinity(z_upper[0]) ||
      !selnorm_bound_is_negative_infinity(z_lower[n_bins - 1])) {
    return false;
  }

  for (int b = 0; b < n_bins; ++b) {
    if (!std::isfinite(z_lower[b]) &&
        !selnorm_bound_is_negative_infinity(z_lower[b]) &&
        !selnorm_bound_is_positive_infinity(z_lower[b])) {
      return false;
    }
    if (!std::isfinite(z_upper[b]) &&
        !selnorm_bound_is_negative_infinity(z_upper[b]) &&
        !selnorm_bound_is_positive_infinity(z_upper[b])) {
      return false;
    }
    if (!(z_lower[b] < z_upper[b])) {
      return false;
    }
  }

  for (int b = 1; b < n_bins; ++b) {
    if (!selnorm_same_finite_bound(z_upper[b], z_lower[b - 1])) {
      return false;
    }
  }

  return true;
}
