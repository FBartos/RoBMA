// Include last, after library headers, in R-facing translation units only.
// Never include this remapping in JAGS kernels or public headers.
#ifndef ROBMA_R_NATIVE_API_H
#define ROBMA_R_NATIVE_API_H
#include "r-native-boundary.h"

#define Rf_error(...) robma_native::error(__VA_ARGS__)
#define R_CheckUserInterrupt() robma_native::check_interrupt()
#define Rf_allocVector(...) robma_native::call(&(Rf_allocVector), __VA_ARGS__)
#define Rf_allocMatrix(...) robma_native::call(&(Rf_allocMatrix), __VA_ARGS__)
#define R_alloc(...) robma_native::call(&(R_alloc), __VA_ARGS__)
#define Rf_coerceVector(...) robma_native::call(&(Rf_coerceVector), __VA_ARGS__)
#define Rf_getAttrib(...) robma_native::call(&(Rf_getAttrib), __VA_ARGS__)
#define Rf_setAttrib(...) robma_native::call(&(Rf_setAttrib), __VA_ARGS__)
#define Rf_install(...) robma_native::call(&(Rf_install), __VA_ARGS__)
#define Rf_mkChar(...) robma_native::call(&(Rf_mkChar), __VA_ARGS__)
#define Rf_mkString(...) robma_native::call(&(Rf_mkString), __VA_ARGS__)
#define Rf_ScalarReal(...) robma_native::call(&(Rf_ScalarReal), __VA_ARGS__)
#define Rf_ScalarInteger(...) robma_native::call(&(Rf_ScalarInteger), __VA_ARGS__)
#define Rf_ScalarLogical(...) robma_native::call(&(Rf_ScalarLogical), __VA_ARGS__)
#define Rf_asReal(...) robma_native::call(&(Rf_asReal), __VA_ARGS__)
#define Rf_inherits(...) robma_native::call(&(Rf_inherits), __VA_ARGS__)
#define Rf_dnorm4(...) robma_native::normal_density(__VA_ARGS__)
#define Rf_pnorm5(...) robma_native::normal_cdf(__VA_ARGS__)
#define R_MakeExternalPtr(...) robma_native::call(&(R_MakeExternalPtr), __VA_ARGS__)
#define R_RegisterCFinalizerEx(...) robma_native::call(&(R_RegisterCFinalizerEx), __VA_ARGS__)
#define SET_STRING_ELT(...) robma_native::call(&(SET_STRING_ELT), __VA_ARGS__)
#define SET_VECTOR_ELT(...) robma_native::call(&(SET_VECTOR_ELT), __VA_ARGS__)
#undef PROTECT
#define PROTECT(x) robma_native::call(&(Rf_protect), x)
#define REAL(x) robma_native::real(x)
#define INTEGER(x) robma_native::integer(x)
#define LOGICAL(x) robma_native::logical(x)
#define RAW(x) robma_native::raw(x)
#define VECTOR_ELT(...) robma_native::vector_elt(__VA_ARGS__)
#define STRING_ELT(...) robma_native::string_elt(__VA_ARGS__)
#define Rf_length(x) robma_native::length(x)
#define XLENGTH(x) robma_native::xlength(x)
#undef CHAR
#define CHAR(x) robma_native::character(x)
#define GetRNGstate() robma_native::get_rng_state()
#define PutRNGstate() robma_native::put_rng_state()

#endif
