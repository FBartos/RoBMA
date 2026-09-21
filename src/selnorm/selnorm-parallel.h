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
// thread; batch code only reads it. A fresh process starts serial, matching
// the option's NA default: post-fit calls thread only after the R side
// resolves a budget for them. Zero means the OpenMP default count.
inline std::atomic<int> &robma_native_threads_value()
{
  static std::atomic<int> value{1};
  return value;
}

// Installs a new budget and returns the one it replaces, so a scoped caller
// can restore what it found rather than a fixed value.
inline int robma_set_native_threads(int threads)
{
  return robma_native_threads_value().exchange(threads < 1 ? 0 : threads,
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

// How a batch of independent rows is distributed. 'chunk_rows' is the number of
// rows one parallel region covers; the remaining rows follow in further
// regions, between which the main thread checks for an interrupt.
struct RobmaRowSchedule {
  int threads;
  int chunk_rows;
};

// Row count alone does not say whether threads pay for themselves: a batch of
// many rows that each evaluate the kernel a handful of times spends more on
// opening parallel regions than it saves. Callers therefore state the work one
// row carries in kernel evaluations - observations, times grid points, times
// selection bins where the kernel walks them - and these two constants turn
// that into a thread count and a region size.
//
// ROBMA_WORK_PER_THREAD is the work a thread must be handed before it is worth
// waking; ROBMA_WORK_PER_THREAD_REGION is how much work each thread gets inside
// one parallel region, which bounds both the region-spawn overhead per unit of
// work and how long the main thread waits before its next interrupt check.
// Both are calibrated over the shapes the post-fit batches use; see
// tests/testthat/test-00-selection-kernel-threads.R for the invariance the
// choice must preserve.
#define ROBMA_WORK_PER_THREAD 2.0e4
#define ROBMA_WORK_PER_THREAD_REGION 1.0e6

inline RobmaRowSchedule robma_row_schedule(int rows, double work_per_row)
{
  RobmaRowSchedule schedule;
  schedule.threads = 1;
  schedule.chunk_rows = rows;
#if defined(_OPENMP)
  if (rows < 2 || !(work_per_row > 0.0)) return schedule;
  int requested = robma_native_threads_value().load(std::memory_order_relaxed);
  if (requested <= 0) requested = omp_get_max_threads();
  const int available = omp_get_max_threads();
  if (available <= 1 || requested <= 1) return schedule;
  const int budget = std::min(std::min(requested, available), rows);
  if (budget <= 1) return schedule;
  const double total = static_cast<double>(rows) * work_per_row;
  const double affordable = std::floor(total / ROBMA_WORK_PER_THREAD);
  if (!(affordable >= 2.0)) return schedule;
  schedule.threads = affordable >= static_cast<double>(budget) ?
    budget : static_cast<int>(affordable);
  const double region = ROBMA_WORK_PER_THREAD_REGION *
    static_cast<double>(schedule.threads) / work_per_row;
  const double capped = std::min(static_cast<double>(rows), std::ceil(region));
  schedule.chunk_rows = std::max(schedule.threads,
    capped >= 1.0 ? static_cast<int>(capped) : 1);
#else
  (void)work_per_row;
#endif
  return schedule;
}

// The same gate for a batch whose parallel region the caller opens itself.
inline int robma_work_threads(int rows, double work_per_row)
{
  return robma_row_schedule(rows, work_per_row).threads;
}

namespace robma_parallel {

// Applies 'body(s)' to every row in [0, rows). With more than one thread the
// rows are split into chunks processed by independent OpenMP workers; 'body'
// must not call the R API and must not fail (validate inputs beforehand).
// A C++ exception escaping 'body' in any worker is re-raised as an R error on
// the main thread after the loop completes, mirroring the serial behaviour.
// 'chunk_rows' of zero keeps the row-count-only default chunking.
template <typename Body>
void for_rows(int rows, int threads, int chunk_rows, Body &&body)
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
    // Chunks keep interrupt responsiveness and bound region-spawn overhead.
    // Without a work-aware chunk size each chunk carries at least eight rows
    // per worker of useful work.
    const int step = chunk_rows > 0 ? chunk_rows : std::max(8, rows / (threads * 8));
    for (int start = 0; start < rows && !failed; start += step) {
      R_CheckUserInterrupt();
      const int end = std::min(rows, start + step);
      #pragma omp parallel num_threads(threads)
      {
        #pragma omp for schedule(dynamic, 1)
        for (int s = start; s < end; ++s) guarded_body(s);
      }
    }
#else
    (void)threads;
    (void)chunk_rows;
    for (int s = 0; s < rows; ++s) {
      if (failed) break;
      guarded_body(s);
    }
#endif
  }
  if (failed) Rf_error("%s", failure_message.c_str());
}

template <typename Body>
void for_rows(int rows, int threads, Body &&body)
{
  for_rows(rows, threads, 0, std::forward<Body>(body));
}

// Runs 'body(s)' over the rows with the thread count and region size the
// schedule chose. One thread keeps the plain serial loop, so a serial batch
// gains neither an interrupt check per row nor a changed access order.
template <typename Body>
void for_rows(int rows, const RobmaRowSchedule &schedule, Body &&body)
{
  if (schedule.threads <= 1) {
    for (int s = 0; s < rows; ++s) body(s);
    return;
  }
  for_rows(rows, schedule.threads, schedule.chunk_rows, std::forward<Body>(body));
}

}

#endif
