// Probe of the key and lifetime rules of the selection-node memo
// (src/distributions/selnorm-node-memo.h). The test compiles and runs it: it
// mutates inputs of a recorded call in place, the way JAGS mutates parent
// values, which no R-level model can do. Each check prints "name 0|1".
#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

#include "selnorm-node-memo.h"

namespace {

// One call of the twenty-parameter step layout with a rank-two block of four
// rows. The arrays stay owned here, so a test mutates them in place.
struct Call {
  std::vector<double> data[20];
  std::vector<double> x;
  std::vector<const double *> par;
  std::vector<unsigned int> len;
  int type = 0;
  SelNormMemoInput inputs[selnorm_step_memo_inputs];
  double scalars[selnorm_step_memo_scalars];

  explicit Call(unsigned int rows = 4U)
  {
    const unsigned int rank = 2U;
    data[0].assign(rows, 0.1);
    data[1].assign(rows, 0.5);
    data[2].assign(rows * rank, 0.3);
    data[3].assign(rows, 0.2);
    data[4] = {1.0, 0.4, 0.2};
    data[5] = {0.0, -1.0, -1e300};
    data[6] = {1e300, 0.0, -1.0};
    data[7].assign(rows, 1.0);
    data[8] = {1.0};
    data[9] = {1.0};
    data[10] = {1.0};
    for (int i = 0; i < 12; ++i) data[11].push_back(-1.5 + 0.25 * i);
    for (int i = 0; i < 12; ++i) data[12].push_back(-0.5 - 0.1 * i);
    data[13] = {3.0, 4.0, 5.0};
    data[14].assign(32, 0.5);
    data[15] = {8.0};
    data[16] = {16.0};
    data[17] = {2.0};
    data[18] = {0.005};
    data[19] = {0.0};
    x.assign(rows, 0.25);
    for (unsigned int i = 0; i < rows; ++i) { x[i] += 0.01 * i; data[0][i] += 0.02 * i; }
    for (int i = 0; i < 20; ++i) {
      par.push_back(data[i].data());
      len.push_back(static_cast<unsigned int>(data[i].size()));
    }
  }

  const SelNormMemoInput *collect()
  {
    selnorm_step_memo_collect(x.data(), static_cast<unsigned int>(x.size()), type,
                              par, len, inputs, scalars);
    return inputs;
  }
};

// Every array of the call that is an input, with a label.
std::vector<std::vector<double> *> input_arrays(Call &call)
{
  std::vector<std::vector<double> *> arrays;
  arrays.push_back(&call.x);
  for (int i = 0; i < 20; ++i) if (i != 14) arrays.push_back(&call.data[i]);
  return arrays;
}

bool hit(SelNormNodeMemo &memo, const void *key, Call &call, double *value = nullptr)
{
  double found = 0.0;
  const bool result = memo.find(key, call.collect(), selnorm_step_memo_inputs, &found);
  if (value != nullptr) *value = found;
  return result;
}

void store(SelNormNodeMemo &memo, const void *key, Call &call, double value)
{
  memo.store(key, call.collect(), selnorm_step_memo_inputs, value);
}

void report(const char *name, bool ok)
{
  std::printf("%s %d\n", name, ok ? 1 : 0);
}

}  // namespace

int main()
{
  const void *node = reinterpret_cast<const void *>(0x1000);
  {
    SelNormNodeMemo memo;
    Call call;
    report("empty_misses", !hit(memo, node, call));
    store(memo, node, call, 1.5);
    double value = 0.0;
    report("stored_call_hits", hit(memo, node, call, &value) && value == 1.5);
    report("other_node_misses", !hit(memo, reinterpret_cast<const void *>(0x2000), call));
    call.type = 1;
    report("density_type_is_input", !hit(memo, node, call));
    call.type = 0;
    report("type_restored_hits", hit(memo, node, call));
  }
  {
    // Every element of every input, mutated in place, must miss; the design
    // array must not matter. Then a store must carry the mutation, and only
    // the changed array, so the next call with the old value misses.
    SelNormNodeMemo memo;
    Call call;
    store(memo, node, call, 1.5);
    int checked = 0, missed = 0, restored = 0, stored = 0, reverted = 0;
    for (std::vector<double> *array : input_arrays(call)) {
      for (std::size_t element = 0; element < array->size(); ++element) {
        const double original = (*array)[element];
        (*array)[element] = std::nextafter(original, std::numeric_limits<double>::infinity());
        ++checked;
        if (!hit(memo, node, call)) ++missed;
        store(memo, node, call, 2.5);
        double value = 0.0;
        if (hit(memo, node, call, &value) && value == 2.5) ++stored;
        (*array)[element] = original;
        if (!hit(memo, node, call)) ++reverted;
        store(memo, node, call, 1.5);
        if (hit(memo, node, call, &value) && value == 1.5) ++restored;
      }
    }
    std::printf("elements_checked %d\n", checked);
    report("mutated_element_misses", missed == checked);
    report("store_carries_mutation", stored == checked);
    report("store_replaces_old_value", reverted == checked);
    report("restored_call_hits", restored == checked);
    int design_hits = 0;
    for (std::size_t element = 0; element < call.data[14].size(); ++element) {
      const double original = call.data[14][element];
      call.data[14][element] = 0.75;
      if (hit(memo, node, call)) ++design_hits;
      call.data[14][element] = original;
    }
    report("qmc_design_is_not_an_input", design_hits == static_cast<int>(call.data[14].size()));
  }
  {
    // A different length is a different call; the entry follows the new layout.
    SelNormNodeMemo memo;
    Call four(4U), three(3U);
    store(memo, node, four, 1.5);
    report("length_change_misses", !hit(memo, node, three));
    store(memo, node, three, 3.5);
    double value = 0.0;
    report("new_layout_hits", hit(memo, node, three, &value) && value == 3.5);
    report("old_layout_misses", !hit(memo, node, four));
  }
  {
    // The map is bounded: a new key beyond the limit empties it, while
    // updating an existing key at the limit does not.
    SelNormNodeMemo memo;
    Call call;
    const int limit = SelNormNodeMemo::maximum_keys;
    for (int i = 1; i <= limit; ++i) {
      store(memo, reinterpret_cast<const void *>(static_cast<std::size_t>(i)), call, i);
    }
    report("limit_keys_all_hit", hit(memo, reinterpret_cast<const void *>(std::size_t(1)), call) &&
           hit(memo, reinterpret_cast<const void *>(static_cast<std::size_t>(limit)), call));
    store(memo, reinterpret_cast<const void *>(static_cast<std::size_t>(limit)), call, -1.0);
    report("update_at_limit_keeps_keys", hit(memo, reinterpret_cast<const void *>(std::size_t(1)), call));
    store(memo, reinterpret_cast<const void *>(static_cast<std::size_t>(limit + 1)), call, 7.0);
    report("new_key_beyond_limit_clears",
           !hit(memo, reinterpret_cast<const void *>(std::size_t(1)), call) &&
           !hit(memo, reinterpret_cast<const void *>(static_cast<std::size_t>(limit)), call) &&
           hit(memo, reinterpret_cast<const void *>(static_cast<std::size_t>(limit + 1)), call));
  }
  return 0;
}
