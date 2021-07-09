/*
    This file is based on file in the runjags package (version 2.0)
    The previous version of the file is Copyright (C) Matthew Denwood, licensed under GPL-2.
 */

#include <R.h>
#include <Rmath.h>

// Checks the JAGS version and sets necessary macros:
#include "jagsversions.h"

// Trivial function to check which version of JAGS the binary was compiled against:
extern "C" {
void getjagsversions(int *forced, int *assumed, int *detected, int *used) // All pointers as everything from R is a vector
{
	forced[0]   = (int) JAGS_MAJOR_FORCED;
	assumed[0]  = (int) JAGS_MAJOR_ASSUMED;
	detected[0] = (int) JAGS_MAJOR;
	used[0]     = (int) JAGS_MAJOR_USED;
}
}
