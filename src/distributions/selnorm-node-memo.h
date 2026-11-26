#ifndef SELNORM_NODE_MEMO_H_
#define SELNORM_NODE_MEMO_H_

#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <unordered_map>
#include <vector>

// Evaluation memo of a selection step node (dselnorm_factor_step and
// dselnorm_cluster_step use it): the last two evaluations of each node.
//
// JAGS's slice sampler evaluates each update first at the current value, whose
// inputs equal the node's previous evaluation, so a node is often asked for a
// log density it has just returned (about 17% of the calls on the reduced
// assink2016 CCI fit). A sampler that ends an update by restoring an earlier
// state (a rejected Metropolis correction, an enumeration of an indicator's
// levels) leaves that state one evaluation back, so the memo keeps two entries
// per node and evicts the one used least recently. Each entry holds owned
// copies of every input that determines a deterministic result and the density
// returned for them. A call is answered from an entry only when every stored
// value equals the call's value. JAGS mutates parent values in place, so a
// pointer never stands for its contents: the node's `x` pointer only selects
// the node.
//
// The QMC design array is not an input. Only results the kernel flags as
// deterministic are stored (the ladders of the product and the best weight
// rule, the refinement tier and the sparse grid); those never read the design,
// and the design is a fixed data array of the node. A randomized result is never
// stored, even when its relative Monte Carlo standard error happens to be
// exactly 0. A remembered density is the one the kernel would return, so the
// draws are the same with any number of entries and without the memo.
//
// JAGS can evaluate different chains concurrently. The lock protects the shared
// map while copying or comparing an entry; density evaluation runs outside it.
//
// A kernel gains only when the density it skips costs more than the comparison
// of every input, which includes the quadrature tables. On that fit the factor
// kernel (about 140 us a call) saves 17% of its time for 0.2-0.3 us a call. The
// cluster kernel (about 8 us a call) saves 17% as well and pays about 0.85 us
// for the comparison of its 32 KB of tables on a hit and 0.8 us for the store
// after a miss, a net gain of about 7% of its time: a cluster kernel that costs
// less than about 5 us a call would lose by it. A miss compares the second
// entry too, which stops at the first differing input (the per-row arrays and
// weights come first).

// ROBMA_SELNORM_NODE_MEMO=0 disables the memo and =1 keeps one entry per node
// (the earlier memo); anything else, or no value, keeps two. The variable is
// read once per process, so on/off comparisons run in separate processes. The
// entry count can also be changed at run time by the test hook
// (selnorm_node_memo_set_slots), which clears nothing: a node's entries beyond
// the new count are ignored and then reused.
inline std::atomic<int> &selnorm_node_memo_slot_setting()
{
  static std::atomic<int> slots([]() {
    const char *value = std::getenv("ROBMA_SELNORM_NODE_MEMO");
    if (value != nullptr && std::strcmp(value, "0") == 0) return 0;
    if (value != nullptr && std::strcmp(value, "1") == 0) return 1;
    return 2;
  }());
  return slots;
}

inline int selnorm_node_memo_slots()
{
  return selnorm_node_memo_slot_setting().load(std::memory_order_relaxed);
}

inline bool selnorm_node_memo_enabled()
{
  return selnorm_node_memo_slots() > 0;
}

// Test hook: 0, 1 or 2 entries per node from now on.
inline void selnorm_node_memo_set_slots(int slots)
{
  selnorm_node_memo_slot_setting().store(slots < 0 ? 0 : (slots > 2 ? 2 : slots),
                                         std::memory_order_relaxed);
}

// Process counters of the memos (test hook): lookups and hits of the cluster
// (index 0) and factor (index 1) step nodes.
struct SelNormMemoCounters {
  std::atomic<std::uint64_t> lookups[2];
  std::atomic<std::uint64_t> hits[2];
};

inline SelNormMemoCounters &selnorm_node_memo_counters()
{
  static SelNormMemoCounters *counters = []() {
    SelNormMemoCounters *out = new SelNormMemoCounters();
    for (int i = 0; i < 2; ++i) {
      out->lookups[i].store(0, std::memory_order_relaxed);
      out->hits[i].store(0, std::memory_order_relaxed);
    }
    return out;
  }();
  return *counters;
}

struct SelNormMemoInput {
  const double *values;
  std::size_t count;
};

class SelNormNodeMemo {
public:
  // The map is cleared before it would hold more keys than this.
  enum { maximum_keys = 4096, maximum_slots = 2 };

  // `kind` selects the process counters (0 cluster, 1 factor).
  explicit SelNormNodeMemo(int kind) : kind_(kind) {}

  // True, with the remembered density, when one of the node's stored calls had
  // exactly these inputs (the most recently used entry is compared first).
  // Lengths are compared first, then the inputs in the order given, so the
  // caller lists the ones that change most often first.
  bool find(const void *node, const SelNormMemoInput *inputs,
            std::size_t count, double *log_density) const
  {
    const int slots = selnorm_node_memo_slots();
    SelNormMemoCounters &counters = selnorm_node_memo_counters();
    counters.lookups[kind_].fetch_add(1, std::memory_order_relaxed);
    const std::lock_guard<std::mutex> lock(mutex_);
    const auto found = entries_.find(node);
    if (found == entries_.end()) return false;
    Node const &stored = found->second;
    const int first = slots > 1 && stored.slots[1].used > stored.slots[0].used ? 1 : 0;
    for (int k = 0; k < slots; ++k) {
      const int slot = k == 0 ? first : 1 - first;
      Slot const &entry = stored.slots[slot];
      if (entry.used == 0 || !matches(entry, inputs, count)) continue;
      *log_density = entry.log_density;
      entry.used = ++clock_;
      counters.hits[kind_].fetch_add(1, std::memory_order_relaxed);
      return true;
    }
    return false;
  }

  // Remembers the density for these inputs in the node's unused or least
  // recently used entry (the only entry with one slot), copying only the inputs
  // that differ from what that entry holds.
  void store(const void *node, const SelNormMemoInput *inputs,
             std::size_t count, double log_density)
  {
    const int slots = selnorm_node_memo_slots();
    if (slots < 1) return;
    const std::lock_guard<std::mutex> lock(mutex_);
    if (entries_.size() >= maximum_keys && entries_.find(node) == entries_.end()) {
      entries_.clear();
    }
    Node &stored = entries_[node];
    int target = 0;
    if (slots > 1) {
      if (stored.slots[0].used == 0) target = 0;
      else if (stored.slots[1].used == 0) target = 1;
      else target = stored.slots[1].used < stored.slots[0].used ? 1 : 0;
    }
    Slot &entry = stored.slots[target];
    bool same_layout = entry.counts.size() == count;
    for (std::size_t i = 0; same_layout && i < count; ++i) {
      same_layout = entry.counts[i] == inputs[i].count;
    }
    try {
      if (same_layout) {
        double *values = entry.values.data();
        for (std::size_t i = 0; i < count; ++i) {
          if (inputs[i].count != 0 &&
              std::memcmp(values, inputs[i].values,
                          inputs[i].count * sizeof(double)) != 0) {
            std::memcpy(values, inputs[i].values, inputs[i].count * sizeof(double));
          }
          values += inputs[i].count;
        }
      } else {
        std::vector<std::size_t> counts(count);
        std::size_t total = 0;
        for (std::size_t i = 0; i < count; ++i) {
          counts[i] = inputs[i].count;
          total += inputs[i].count;
        }
        std::vector<double> values(total);
        double *copy = values.data();
        for (std::size_t i = 0; i < count; ++i) {
          if (inputs[i].count != 0) {
            std::memcpy(copy, inputs[i].values, inputs[i].count * sizeof(double));
          }
          copy += inputs[i].count;
        }
        entry.counts.swap(counts);
        entry.values.swap(values);
      }
      entry.log_density = log_density;
      entry.used = ++clock_;
    } catch (...) {
      // An allocation failure must not leave an entry whose values and density
      // disagree.
      entry.counts.clear();
      entry.values.clear();
      entry.used = 0;
      throw;
    }
  }

private:
  struct Slot {
    std::vector<std::size_t> counts;
    std::vector<double> values;
    double log_density = 0.0;
    // the memo's clock at the last store or hit; 0 for an unused entry
    mutable std::uint64_t used = 0;
  };
  struct Node {
    Slot slots[maximum_slots];
  };

  static bool matches(Slot const &entry, const SelNormMemoInput *inputs, std::size_t count)
  {
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
    return true;
  }

  const int kind_;
  mutable std::mutex mutex_;
  mutable std::uint64_t clock_ = 0;
  std::unordered_map<const void *, Node> entries_;
};

// Inputs of one call of the factor and cluster step distributions, which share
// their twenty-parameter layout: mean, residual SD, loading, selection SE,
// weights, step bounds, observed bins, effect sign, telescoping flag, kernel,
// quadrature nodes, log weights and orders, QMC design, point counts,
// scrambles, tolerance and weight rule. Everything but the QMC design is an
// input; the other scalar controls are packed with the density type.
enum { selnorm_step_memo_inputs = 13, selnorm_step_memo_scalars = 9 };

inline void selnorm_step_memo_collect(
    const double *x, unsigned long length, int type,
    const std::vector<const double *> &par,
    const std::vector<unsigned long> &len,
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
