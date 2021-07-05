/*
    This file is based on file in the runjags package (version 2.0)
    The previous version of the file is Copyright (C) Matthew Denwood, licensed under GPL-2.
*/

#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h> // for SEXP

extern void getjagsversions(int *forced, int *assumed, int *detected, int *used);

static R_NativePrimitiveArgType getjagsversions_t[] = {
    INTSXP, INTSXP, INTSXP, INTSXP
};

static const R_CMethodDef cMethods[] = {
    {"getjagsversions", (DL_FUNC) &getjagsversions, 4, getjagsversions_t},
    {NULL, NULL, 0, NULL}
};

void
R_init_RoBMA(DllInfo *dll)
{
    R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
