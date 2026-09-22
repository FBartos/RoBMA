#ifndef ROBMA_R_NATIVE_BOUNDARY_H
#define ROBMA_R_NATIVE_BOUNDARY_H

#include <Rinternals.h>
#include <R_ext/Error.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>
#include <Rmath.h>

#include <csetjmp>
#include <cstdarg>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <exception>
#include <initializer_list>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>
#if defined(_OPENMP)
#include <omp.h>
#endif

// R's documented continuation-token API, with a setjmp trampoline as used by
// cpp11, catches an R longjmp at each R API call, then unwinds the caller's C++
// scopes. Protecting the whole entry body would skip those scopes instead.
// No exceptions cross .Call or OpenMP boundaries, and no dependency is added.
namespace robma_native {

struct Unwind { SEXP token; }; // Do not let local std::exception handlers swallow it.

class Error : public std::exception {
 public:
  char message[8192] = {};
  const char *what() const noexcept override { return message; }
};

[[noreturn]] inline void error(const char *format, ...)
{
  Error failure;
  va_list args;
  va_start(args, format);
  std::vsnprintf(failure.message, sizeof(failure.message), format, args);
  va_end(args);
  throw failure;
}

inline bool in_worker()
{
#if defined(_OPENMP)
  return omp_in_parallel() != 0;
#else
  return false;
#endif
}

struct Buffer { SEXP object; void *data; R_xlen_t length; bool preserved; };
struct Context {
  SEXP token;
  bool rng_open = false;
  std::vector<Buffer> buffers;
  std::initializer_list<SEXP> arguments;
  ~Context() {
    // R_ReleaseObject is a nonallocating removal from the preserve list.
    for (const Buffer &buffer : buffers) {
      if (buffer.preserved) R_ReleaseObject(buffer.object);
    }
  }
};

// Only the main thread changes this pointer, before/after worker regions.
// OpenMP publication/join barriers make the prepared ALTREP views read-only in
// workers. Nested main-thread .Call entries restore their parent's context.
inline Context *active_context = nullptr;

template <class Function>
SEXP protect_call(Function &&function)
{
  if (in_worker() || active_context == nullptr) {
    error("An R API operation was requested outside the native main-thread boundary.");
  }
  SEXP token = active_context->token;
  std::jmp_buf jump_buffer;
  if (setjmp(jump_buffer)) throw Unwind{token};
  using Callback = typename std::remove_reference<Function>::type;
  SEXP result = R_UnwindProtect(
    [](void *data) -> SEXP { return (*static_cast<Callback *>(data))(); },
    &function,
    [](void *data, Rboolean jump) {
      if (jump) std::longjmp(*static_cast<std::jmp_buf *>(data), 1);
    }, &jump_buffer, token);
  SETCAR(token, R_NilValue);
  return result;
}

// Arguments are evaluated BEFORE entering R_UnwindProtect. In particular,
// ScalarReal(native_calculation()) must unwind native_calculation normally.
template <class Function, class... Args>
auto call(Function function, Args... args) -> decltype(function(args...))
{
  using Result = decltype(function(args...));
  if constexpr (std::is_same<Result, SEXP>::value) {
    return protect_call([&]() -> SEXP { return function(args...); });
  } else if constexpr (std::is_void<Result>::value) {
    protect_call([&]() -> SEXP { function(args...); return R_NilValue; });
  } else {
    Result result{};
    protect_call([&]() -> SEXP { result = function(args...); return R_NilValue; });
    return result;
  }
}

inline void check_interrupt() { call(R_CheckUserInterrupt); }

inline void *altrep_data(SEXP object)
{
  if (active_context != nullptr) {
    for (const Buffer &buffer : active_context->buffers) {
      if (buffer.object == object) return buffer.data;
    }
  }
  if (in_worker()) error("ALTREP inputs must be materialized before native workers start.");
  const R_xlen_t length = call(Rf_xlength, object);
  void *data = nullptr;
  switch (TYPEOF(object)) {
  case REALSXP: data = call(REAL, object); break;
  case INTSXP: data = call(INTEGER, object); break;
  case LGLSXP: data = call(LOGICAL, object); break;
  case RAWSXP: data = call(RAW, object); break;
  default: error("Unsupported native ALTREP data type.");
  }
  active_context->buffers.push_back(Buffer{object, data, length, false});
  call(R_PreserveObject, object);
  active_context->buffers.back().preserved = true;
  return data;
}

inline double *real(SEXP object)
{
  if (TYPEOF(object) != REALSXP) return call(REAL, object);
  return ALTREP(object) ? static_cast<double *>(altrep_data(object)) : REAL(object);
}
inline int *integer(SEXP object)
{
  if (TYPEOF(object) != INTSXP) return call(INTEGER, object);
  return ALTREP(object) ? static_cast<int *>(altrep_data(object)) : INTEGER(object);
}
inline int *logical(SEXP object)
{
  if (TYPEOF(object) != LGLSXP) return call(LOGICAL, object);
  return ALTREP(object) ? static_cast<int *>(altrep_data(object)) : LOGICAL(object);
}
inline Rbyte *raw(SEXP object)
{
  if (TYPEOF(object) != RAWSXP) return call(RAW, object);
  return ALTREP(object) ? static_cast<Rbyte *>(altrep_data(object)) : RAW(object);
}
inline void prepare_workers()
{
  if (in_worker() || active_context == nullptr) {
    error("Native workers must be prepared on the main thread.");
  }
  // Delay materialization until input validation has completed. Eagerly
  // realizing an invalid huge ALTREP control before its length check could
  // turn a cheap validation error into an unnecessary large allocation.
  for (SEXP argument : active_context->arguments) {
    const int type = TYPEOF(argument);
    if (ALTREP(argument) && call(Rf_xlength, argument) <= std::numeric_limits<int>::max() &&
        (type == REALSXP || type == INTSXP || type == LGLSXP || type == RAWSXP)) {
      altrep_data(argument);
    }
  }
}
inline SEXP vector_elt(SEXP object, R_xlen_t index)
{
  if (TYPEOF(object) != VECSXP || ALTREP(object) || index < 0 || index >= XLENGTH(object)) {
    return call(VECTOR_ELT, object, index);
  }
  return VECTOR_ELT(object, index);
}
inline SEXP string_elt(SEXP object, R_xlen_t index)
{
  if (TYPEOF(object) != STRSXP || ALTREP(object) || index < 0 || index >= XLENGTH(object)) {
    return call(STRING_ELT, object, index);
  }
  return STRING_ELT(object, index);
}

inline R_xlen_t xlength(SEXP object)
{
  if (!ALTREP(object)) {
    switch (TYPEOF(object)) {
    case NILSXP: case CHARSXP: case LGLSXP: case INTSXP: case REALSXP:
    case CPLXSXP: case STRSXP: case VECSXP: case EXPRSXP: case RAWSXP:
      return XLENGTH(object);
    default: return call(Rf_xlength, object);
    }
  }
  if (active_context != nullptr) {
    for (const Buffer &buffer : active_context->buffers) {
      if (buffer.object == object) return buffer.length;
    }
  }
  return call(Rf_xlength, object);
}
inline int length(SEXP object)
{
  const R_xlen_t size = xlength(object);
  if (size > std::numeric_limits<int>::max()) error("Native input exceeds integer indexing limits.");
  return static_cast<int>(size);
}
inline const char *character(SEXP object)
{
  return TYPEOF(object) == CHARSXP ? CHAR(object) : call(R_CHAR, object);
}

// Ordinary Rmath inputs remain on their original fast path. Domain-warning
// cases can become R errors under options(warn = 2), so protect those calls.
inline double normal_density(double x, double mean, double sd, int log)
{
  if (sd < 0 || (!std::isfinite(x) && x == mean)) {
    return call(Rf_dnorm4, x, mean, sd, log);
  }
  return Rf_dnorm4(x, mean, sd, log);
}
inline double normal_cdf(double x, double mean, double sd, int lower, int log)
{
  if (sd < 0 || (!std::isfinite(x) && x == mean)) {
    return call(Rf_pnorm5, x, mean, sd, lower, log);
  }
  return Rf_pnorm5(x, mean, sd, lower, log);
}

inline void get_rng_state()
{
  call(GetRNGstate);
  active_context->rng_open = true;
}
inline void put_rng_state()
{
  call(PutRNGstate);
  active_context->rng_open = false;
}

template <class Function>
SEXP entry(std::initializer_list<SEXP> arguments, Function &&function)
{
  // These R allocations precede all owned C++ state.
  SEXP token = PROTECT(R_MakeUnwindCont());
  SEXP result = R_NilValue;
  SEXP continuation = R_NilValue;
  char message[8192] = {};
  bool rng_open = false;
  {
    Context context{token, false, {}, arguments};
    Context *previous = active_context;
    active_context = &context;
    try {
      result = function();
      SETCAR(token, result); // Root the result across optional RNG-state saving.
    } catch (const Unwind &failure) {
      continuation = failure.token;
    } catch (const std::exception &failure) {
      std::strncpy(message, failure.what(), sizeof(message) - 1);
      if (!message[0]) std::strcpy(message, "Native C++ evaluation failed.");
    } catch (...) {
      std::strcpy(message, "Native C++ evaluation failed with an unknown exception.");
    }
    rng_open = context.rng_open;
    active_context = previous;
  }
  // No owned C++ state or active exception remains at these R longjmp sites.
  // Save consumed RNG draws even if an error/interrupt left the simulation.
  if (rng_open) PutRNGstate();
  if (continuation != R_NilValue) R_ContinueUnwind(continuation);
  if (message[0]) Rf_error("%s", message);
  UNPROTECT(1);
  return result;
}

} // namespace robma_native

#define ROBMA_NATIVE_BEGIN(...) return robma_native::entry({__VA_ARGS__}, [&]() -> SEXP {
#define ROBMA_NATIVE_END });

#endif
