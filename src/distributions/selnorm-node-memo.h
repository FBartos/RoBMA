#ifndef SELNORM_NODE_MEMO_H_
#define SELNORM_NODE_MEMO_H_

#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <vector>

// Last-evaluation memo of a selection step node (dselnorm_factor_step and
// dselnorm_cluster_step use it).
//
// JAGS's slice sampler evaluates each update first at the current value, whose
// inputs equal the node's previous evaluation, so a node is often asked for a
// log density it has just returned (about 17% of the calls on the reduced
// assink2016 CCI fit). The memo keeps, per node, owned copies of
// every input that determines a deterministic result and the density returned
// for them. A call is answered from it only when every stored value equals the
// call's value. JAGS mutates parent values in place, so a pointer never stands
// for its contents: the node's `x` pointer only selects the entry.
//
// The QMC design array is not an input. Only results the kernel flags as
// deterministic are stored (the ladder, the refinement tier and the sparse
// grid); those never read the design, and the design is a fixed data array of
// the node. A randomized result is never stored, even when its relative Monte
// Carlo standard error happens to be exactly 0.
//
// JAGS runs a process single-threaded, so the memo has no lock. It lives behind
// a function-local pointer that is intentionally never deleted, which avoids
// a destructor running while the DLL unloads.
//
// A kernel gains only when the density it skips costs more than the comparison
// of every input, which includes the quadrature tables. On that fit the factor
// kernel (about 140 us a call) saves 17% of its time for 0.2-0.3 us a call. The
// cluster kernel (about 8 us a call) saves 17% as well and pays about 0.85 us
// for the comparison of its 32 KB of tables on a hit and 0.8 us for the store
// after a miss, a net gain of about 7% of its time: a cluster kernel that costs
// less than about 5 us a call would lose by it.

// ROBMA_SELNORM_NODE_MEMO=0 disables the memo. The variable is read once per
// process, so on/off comparisons run in separate processes.
inline bool selnorm_node_memo_enabled()
{
  static const bool enabled = []() {
    const char *value = std::getenv("ROBMA_SELNORM_NODE_MEMO");
    return value == nullptr || std::strcmp(value, "0") != 0;
  }();
  return enabled;
}

struct SelNormMemoInput {
  const double *values;
  std::size_t count;
};

class SelNormNodeMemo {
public:
  // The map is cleared before it would hold more keys than this.
  enum { maximum_keys = 4096 };

  // True, with the remembered density, when the node's last stored call had
  // exactly these inputs. Lengths are compared first, then the inputs in the
  // order given, so the caller lists the ones that change most often first.
  bool find(const void *node, const SelNormMemoInput *inputs,
            std::size_t count, double *log_density) const
  {
    const auto found = entries_.find(node);
    if (found == entries_.end()) return false;
    const Entry &entry = found->second;
    if (entry.counts.size() != count) return false;
    for (std::size_t i = 0; i < count; ++i) {
      if (entry.counts[i] != inputs[i].count) return false;
    }
    const double *stored = entry.values.data();
    for (std::size_t i = 0; i < count; ++i) {
      if (inputs[i].count != 0 &&
          std::memcmp(stored, inputs[i].values,
                      inputs[i].count * sizeof(double)) != 0) return false;
      stored += inputs[i].count;
    }
    *log_density = entry.log_density;
    return true;
  }

  // Remembers the density for these inputs, copying only the inputs that differ
  // from the node's previous entry.
  void store(const void *node, const SelNormMemoInput *inputs,
             std::size_t count, double log_density)
  {
    if (entries_.size() >= maximum_keys && entries_.find(node) == entries_.end()) {
      entries_.clear();
    }
    Entry &entry = entries_[node];
    bool same_layout = entry.counts.size() == count;
    for (std::size_t i = 0; same_layout && i < count; ++i) {
      same_layout = entry.counts[i] == inputs[i].count;
    }
    try {
      if (same_layout) {
        double *stored = entry.values.data();
        for (std::size_t i = 0; i < count; ++i) {
          if (inputs[i].count != 0 &&
              std::memcmp(stored, inputs[i].values,
                          inputs[i].count * sizeof(double)) != 0) {
            std::memcpy(stored, inputs[i].values, inputs[i].count * sizeof(double));
          }
          stored += inputs[i].count;
        }
      } else {
        std::vector<std::size_t> counts(count);
        std::size_t total = 0;
        for (std::size_t i = 0; i < count; ++i) {
          counts[i] = inputs[i].count;
          total += inputs[i].count;
        }
        std::vector<double> values(total);
        double *stored = values.data();
        for (std::size_t i = 0; i < count; ++i) {
          if (inputs[i].count != 0) {
            std::memcpy(stored, inputs[i].values, inputs[i].count * sizeof(double));
          }
          stored += inputs[i].count;
        }
        entry.counts.swap(counts);
        entry.values.swap(values);
      }
      entry.log_density = log_density;
    } catch (...) {
      // An allocation failure must not leave an entry whose values and density
      // disagree.
      entries_.erase(node);
      throw;
    }
  }

private:
  struct Entry {
    std::vector<std::size_t> counts;
    std::vector<double> values;
    double log_density = 0.0;
  };
  std::unordered_map<const void *, Entry> entries_;
};

// Inputs of one call of the factor and cluster step distributions, which share
// their twenty-parameter layout: mean, residual SD, loading, selection SE,
// weights, step bounds, observed bins, effect sign, telescoping flag, kernel,
// quadrature nodes, log weights and orders, QMC design, point counts,
// scrambles, tolerance and weight rule. Everything but the QMC design is an
// input; the other scalar controls are packed with the density type. JAGS gives
// vector distributions no truncation bounds and the kernels do not read them;
// the callers bypass the memo for a call that has any, since their length is
// unknown here.
enum { selnorm_step_memo_inputs = 13, selnorm_step_memo_scalars = 9 };

inline void selnorm_step_memo_collect(
    const double *x, unsigned int length, int type,
    const std::vector<const double *> &par,
    const std::vector<unsigned int> &len,
    SelNormMemoInput *inputs, double *scalars)
{
  scalars[0] = static_cast<double>(type);
  scalars[1] = *par[8];
  scalars[2] = *par[9];
  scalars[3] = *par[10];
  scalars[4] = *par[15];
  scalars[5] = *par[16];
  scalars[6] = *par[17];
  scalars[7] = *par[18];
  scalars[8] = *par[19];
  // Per-row arrays and the weights change with the sampled parameters, so they
  // come first; the fixed tables, the largest, come last.
  inputs[0]  = SelNormMemoInput{par[0], len[0]};
  inputs[1]  = SelNormMemoInput{par[1], len[1]};
  inputs[2]  = SelNormMemoInput{par[2], len[2]};
  inputs[3]  = SelNormMemoInput{par[4], len[4]};
  inputs[4]  = SelNormMemoInput{par[3], len[3]};
  inputs[5]  = SelNormMemoInput{par[7], len[7]};
  inputs[6]  = SelNormMemoInput{x, length};
  inputs[7]  = SelNormMemoInput{par[5], len[5]};
  inputs[8]  = SelNormMemoInput{par[6], len[6]};
  inputs[9]  = SelNormMemoInput{scalars, selnorm_step_memo_scalars};
  inputs[10] = SelNormMemoInput{par[13], len[13]};
  inputs[11] = SelNormMemoInput{par[11], len[11]};
  inputs[12] = SelNormMemoInput{par[12], len[12]};
}

#endif
