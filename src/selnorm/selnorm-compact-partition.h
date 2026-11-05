#ifndef ROBMA_SELNORM_COMPACT_PARTITION_H
#define ROBMA_SELNORM_COMPACT_PARTITION_H

#include <vector>

// A posterior branch can use only a subset of the shared prior-union cuts.
// This call-owned view removes exact adjacent repeats after validation, without
// changing the context, branch routing, or the weight at an observed boundary.
struct CompactStepPartition {
  std::vector<double> omega, lower, upper;
  std::vector<int> obs_bin, bin_map;

  bool prepare(const double *weights, const double *z_lower,
      const double *z_upper, int bins, const int *observed = nullptr,
      int observations = 0)
  {
    if (bins <= 2) return false;
    bool repeated = false;
    for (int bin = 1; bin < bins; ++bin)
      repeated = repeated || weights[bin] == weights[bin - 1];
    if (!repeated) return false;
    omega.clear(); lower.clear(); upper.clear();
    if (observed != nullptr) bin_map.resize(bins);
    for (int first = 0; first < bins;) {
      int last = first;
      while (last + 1 < bins && weights[last + 1] == weights[first]) ++last;
      omega.push_back(weights[first]);
      lower.push_back(z_lower[last]);
      upper.push_back(z_upper[first]);
      if (observed != nullptr)
        for (int bin = first; bin <= last; ++bin)
          bin_map[bin] = static_cast<int>(omega.size());
      first = last + 1;
    }
    if (observed != nullptr) {
      obs_bin.resize(observations);
      for (int row = 0; row < observations; ++row)
        obs_bin[row] = bin_map[observed[row] - 1];
    }
    return true;
  }
};

#endif
