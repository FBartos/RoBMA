double phack_beta(double alpha, int q, const SelNormKernelData &data,
                  double mean_z, double sd_z)
{
  if (!(alpha > 0) || q <= 0) {
    return 0;
  }

  const double source_lower = data.phack_z_source[0];
  const double source_upper = data.phack_z_source[1];
  const double dest_lower   = data.phack_z_dest[0];
  const double dest_upper   = data.phack_z_dest[1];

  const double source_moment = source_power_moment(
    source_lower, source_upper, mean_z, sd_z, source_lower, source_upper, q
  );
  const double dest_moment = dest_power_moment(
    dest_lower, dest_upper, mean_z, sd_z, dest_lower, dest_upper, q
  );

  if (!(dest_moment > 0)) {
    return 0;
  }

  return alpha * source_moment / dest_moment;
}

double phack_kernel_weight(double z, double alpha, double beta, int q,
                           const SelNormKernelData &data)
{
  if (q <= 0 || !(alpha > 0)) {
    return 1;
  }

  const double source_lower = data.phack_z_source[0];
  const double source_upper = data.phack_z_source[1];
  const double dest_lower   = data.phack_z_dest[0];
  const double dest_upper   = data.phack_z_dest[1];

  if (z >= source_lower && z <= source_upper) {
    const double x = (z - source_lower) / (source_upper - source_lower);
    const double depletion = q == 1 ? x : x * x;
    return 1 - alpha * depletion;
  }

  if (z > dest_lower && z <= dest_upper) {
    const double x = (dest_upper - z) / (dest_upper - dest_lower);
    const double boost = q == 1 ? x : x * x;
    return 1 + beta * boost;
  }

  return 1;
}

double combined_normalizer(const double *omega, double alpha, double beta, int q,
                           const SelNormKernelData &data, double mean_z,
                           double sd_z, int omega_stride)
{
  double out = 0;

  for (int s = 0; s < data.n_segments; ++s) {
    const double lower = data.segment_bounds[s];
    const double upper = data.segment_bounds[s + 1];
    const int step_bin = segment_step_bin_at(data, s);

    if (step_bin < 1 || step_bin > data.n_bins) {
      continue;
    }

    const double w = omega_at(omega, step_bin, omega_stride);
    if (w <= 0) {
      continue;
    }

    double contribution = interval_prob(lower, upper, mean_z, sd_z);
    const int region = segment_phack_region_at(data, s);

    if (region == 1 && alpha > 0 && q > 0) {
      contribution -= alpha * source_power_moment(
        lower, upper, mean_z, sd_z, data.phack_z_source[0],
        data.phack_z_source[1], q
      );
    } else if (region == 2 && beta > 0 && q > 0) {
      contribution += beta * dest_power_moment(
        lower, upper, mean_z, sd_z, data.phack_z_dest[0],
        data.phack_z_dest[1], q
      );
    }

    if (contribution < 0 && contribution > -1e-12) {
      contribution = 0;
    }
    out += w * contribution;
  }

  return out;
}
