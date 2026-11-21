/*
    This file is based on file in the runjags package (version 2.0)
    The previous version of the file is Copyright (C) Matthew Denwood, licensed under GPL-2.
 */

#include <version.h>

// Trivial function to check which version of JAGS the binary was compiled against:
extern "C" {
void getjagsversions(int *forced, int *assumed, int *detected, int *used)
{
    *forced = 0;
    *assumed = 0;
    *detected = JAGS_MAJOR;
    *used = JAGS_MAJOR;
}
}
