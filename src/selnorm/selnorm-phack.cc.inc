double phack_power_weight(double z, double lower, double upper, int q,
                          bool source)
{
  if (!std::isfinite(z) || !std::isfinite(lower) || !std::isfinite(upper) ||
      !(lower < upper) || (q != 1 && q != 2) ||
      z < lower || z > upper) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  if (z == lower) {
    return source ? 0.0 : 1.0;
  }
  if (z == upper) {
    return source ? 1.0 : 0.0;
  }

  const double width = upper - lower;
  double x;
  if (std::isfinite(width)) {
    x = source ? (z - lower) / width : (upper - z) / width;
  } else {
    const double scale        = std::max(std::abs(lower), std::abs(upper));
    const double scaled_lower = lower / scale;
    const double scaled_upper = upper / scale;
    const double scaled_z     = z / scale;
    const double scaled_width = scaled_upper - scaled_lower;
    x = source ?
      (scaled_z - scaled_lower) / scaled_width :
      (scaled_upper - scaled_z) / scaled_width;
  }

  if (!std::isfinite(x) || !(x > 0) || !(x < 1)) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  return q == 1 ? x : x * x;
}

double phack_log_power_moment_tail(double lower, double upper,
                                   double mean, double sd,
                                   double region_lower,
                                   double region_upper, int q,
                                   bool source, bool lower_tail)
{
  static const double nodes[8] = {
    0.09501250983763744, 0.28160355077925891,
    0.45801677765722739, 0.61787624440264375,
    0.75540440835500303, 0.86563120238783174,
    0.94457502307323258, 0.98940093499164993
  };
  static const double weights[8] = {
    0.18945061045506850, 0.18260341504492359,
    0.16915651939500254, 0.14959598881657673,
    0.12462897125553388, 0.09515851168249278,
    0.06225352393864789, 0.02715245941175409
  };

  double log_p_lower;
  double log_p_upper;
  if (lower_tail) {
    log_p_lower = pnorm(lower, mean, sd, true, true);
    log_p_upper = pnorm(upper, mean, sd, true, true);
  } else {
    log_p_lower = pnorm(upper, mean, sd, false, true);
    log_p_upper = pnorm(lower, mean, sd, false, true);
  }

  if (!std::isfinite(log_p_lower) || !std::isfinite(log_p_upper) ||
      !(log_p_upper > log_p_lower)) {
    return -INFINITY;
  }

  const double mid  = 0.5 * (log_p_lower + log_p_upper);
  const double half = 0.5 * (log_p_upper - log_p_lower);
  if (!(half > 0)) {
    return -INFINITY;
  }

  double out = -INFINITY;
  for (int i = 0; i < 8; ++i) {
    for (int sign = -1; sign <= 1; sign += 2) {
      const double v = mid + sign * half * nodes[i];
      const double z = qnorm(v, mean, sd, lower_tail, true);
      const double power = phack_power_weight(
        z, region_lower, region_upper, q, source
      );
      if (!(power > 0)) {
        continue;
      }
      out = add_logspace(
        out,
        std::log(weights[i]) + v + std::log(power)
      );
    }
  }

  return out == -INFINITY ? out : std::log(half) + out;
}

double phack_log_power_moment(double lower, double upper, double mean,
                              double sd, double region_lower,
                              double region_upper, int q, bool source)
{
  lower = std::max(lower, region_lower);
  upper = std::min(upper, region_upper);
  if (!(sd > 0) || lower >= upper || q <= 0) {
    return -INFINITY;
  }

  if (upper <= mean) {
    return phack_log_power_moment_tail(
      lower, upper, mean, sd, region_lower, region_upper, q, source, true
    );
  }
  if (lower >= mean) {
    return phack_log_power_moment_tail(
      lower, upper, mean, sd, region_lower, region_upper, q, source, false
    );
  }

  const double left = phack_log_power_moment_tail(
    lower, mean, mean, sd, region_lower, region_upper, q, source, true
  );
  const double right = phack_log_power_moment_tail(
    mean, upper, mean, sd, region_lower, region_upper, q, source, false
  );
  return add_logspace(left, right);
}

double phack_log_beta(double alpha, int q, const SelNormKernelData &data,
                      double mean_z, double sd_z)
{
  if (!(alpha > 0) || q <= 0) {
    return -INFINITY;
  }

  const double log_source = phack_log_power_moment(
    data.phack_z_source[0], data.phack_z_source[1],
    mean_z, sd_z, data.phack_z_source[0], data.phack_z_source[1],
    q, true
  );
  const double log_dest = phack_log_power_moment(
    data.phack_z_dest[0], data.phack_z_dest[1],
    mean_z, sd_z, data.phack_z_dest[0], data.phack_z_dest[1],
    q, false
  );

  if (!std::isfinite(log_source) || !std::isfinite(log_dest)) {
    return -INFINITY;
  }

  return std::log(alpha) + log_source - log_dest;
}

double phack_log_kernel_weight(double z, double alpha, double log_beta, int q,
                               const SelNormKernelData &data)
{
  if (q <= 0 || !(alpha > 0)) {
    return 0;
  }

  const double source_lower = data.phack_z_source[0];
  const double source_upper = data.phack_z_source[1];
  const double dest_lower   = data.phack_z_dest[0];
  const double dest_upper   = data.phack_z_dest[1];

  if (z >= source_lower && z <= source_upper) {
    const double depletion = phack_power_weight(
      z, source_lower, source_upper, q, true
    );
    if (!(depletion > 0)) {
      return 0;
    }
    return std::log1p(-alpha * depletion);
  }

  if (z > dest_lower && z <= dest_upper) {
    const double boost = phack_power_weight(
      z, dest_lower, dest_upper, q, false
    );
    if (!(boost > 0) || !std::isfinite(log_beta)) {
      return 0;
    }
    return log1pexp(log_beta + std::log(boost));
  }

  return 0;
}

double combined_log_normalizer(const double *omega, double alpha,
                               double log_beta, int q,
                               const SelNormKernelData &data,
                               double mean_z, double sd_z,
                               int omega_stride)
{
  double out = -INFINITY;

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

    double log_contribution = interval_log_prob(lower, upper, mean_z, sd_z);
    const int region = segment_phack_region_at(data, s);

    if (region == 1 && alpha > 0 && q > 0) {
      const double log_moment = phack_log_power_moment(
        lower, upper, mean_z, sd_z, data.phack_z_source[0],
        data.phack_z_source[1], q, true
      );
      if (std::isfinite(log_moment)) {
        log_contribution = sub_logspace(
          log_contribution, std::log(alpha) + log_moment
        );
      }
    } else if (region == 2 && std::isfinite(log_beta) && q > 0) {
      const double log_moment = phack_log_power_moment(
        lower, upper, mean_z, sd_z, data.phack_z_dest[0],
        data.phack_z_dest[1], q, false
      );
      if (std::isfinite(log_moment)) {
        log_contribution = add_logspace(
          log_contribution, log_beta + log_moment
        );
      }
    }

    if (std::isfinite(log_contribution)) {
      out = add_logspace(out, std::log(w) + log_contribution);
    }
  }

  return out;
}
