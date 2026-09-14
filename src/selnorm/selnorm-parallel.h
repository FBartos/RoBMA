// Row-parallel execution helpers for the R-facing selection batches.
//
// These helpers parallelize independent posterior rows of the post-fit batch
// kernels. Every row keeps its deterministic serial arithmetic: rows only
// write to their own output slots and read call-common inputs, so results are
// identical to the serial evaluation regardless of the thread count. The R API
// is never called from worker threads; interrupts are checked between chunks
// on the main thread and per-row validation happens before the parallel
// region. JAGS-side single-state evaluations never enter this path.
#ifndef ROBMA_SELNORM_PARALLEL_H
#define ROBMA_SELNORM_PARALLEL_H

#include <Rinternals.h>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <string>
#include <utility>

#if defined(_OPENMP)
#include <omp.h>
#endif

// Process-wide row budget shared by the R-facing selection batches. The value
// is owned by the R option 'native_threads' and is applied from the main
// thread; batch code only reads it. Zero means the OpenMP default count.
inline std::atomic<int> &robma_native_threads_value()
{
  static std::atomic<int> value{0};
  return value;
}

inline void robma_set_native_threads(int threads)
{
  robma_native_threads_value().store(threads < 1 ? 0 : threads,
                                     std::memory_order_relaxed);
}

// Thread count for a batch of 'rows' independent rows. Small batches stay
// serial: the per-row work does not amortize the parallel-region overhead.
inline int robma_batch_threads(int rows, int minimum_rows = 16)
{
#if defined(_OPENMP)
  if (rows < minimum_rows) return 1;
  int requested = robma_native_threads_value().load(std::memory_order_relaxed);
  if (requested <= 0) requested = omp_get_max_threads();
  const int available = omp_get_max_threads();
  if (available <= 1 || requested <= 1) return 1;
  return std::max(1, std::min(std::min(requested, available), rows));
#else
  (void)rows;
  (void)minimum_rows;
  return 1;
#endif
}

namespace robma_parallel {

// Applies 'body(s)' to every row in [0, rows). With more than one thread the
// rows are split into chunks processed by independent OpenMP workers; 'body'
// must not call the R API and must not fail (validate inputs beforehand).
// A C++ exception escaping 'body' in any worker is re-raised as an R error on
// the main thread after the loop completes, mirroring the serial behaviour.
template <typename Body>
void for_rows(int rows, int threads, Body &&body)
{
  if (rows <= 0) return;
  bool failed = false;
  std::mutex failure_mutex;
  std::string failure_message;
  const auto guarded_body = [&](int s) {
    try {
      body(s);
    } catch (const std::exception &error) {
      std::lock_guard<std::mutex> lock(failure_mutex);
      if (!failed) {
        failed = true;
        failure_message = error.what();
      }
    } catch (...) {
      std::lock_guard<std::mutex> lock(failure_mutex);
      if (!failed) {
        failed = true;
        failure_message = "unknown native row failure";
      }
    }
  };
  if (threads <= 1) {
    for (int s = 0; s < rows; ++s) {
      if (failed) break;
      guarded_body(s);
      R_CheckUserInterrupt();
    }
  } else {
#if defined(_OPENMP)
    // Chunks keep interrupt responsiveness and bound region-spawn overhead:
    // each chunk carries at least eight rows per worker of useful work.
    const int chunk_rows = std::max(8, rows / (threads * 8));
    for (int start = 0; start < rows && !failed; start += chunk_rows) {
      R_CheckUserInterrupt();
      const int end = std::min(rows, start + chunk_rows);
      #pragma omp parallel num_threads(threads)
      {
        #pragma omp for schedule(dynamic, 1)
        for (int s = start; s < end; ++s) guarded_body(s);
      }
    }
#else
    (void)threads;
    for (int s = 0; s < rows; ++s) {
      if (failed) break;
      guarded_body(s);
    }
#endif
  }
  if (failed) Rf_error("%s", failure_message.c_str());
}

}

#endif
