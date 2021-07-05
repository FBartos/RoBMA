/*
    This file is based on file in the runjags package (version 2.0)
    The previous version of the file is Copyright (C) Matthew Denwood, licensed under GPL-2.
	
	This header file sorts the necessary macros for compiling against 
	JAGS 4
 */

#ifndef JAGS_VERSIONS_H_
#define JAGS_VERSIONS_H_

// Get the installed version of JAGS - undefined for JAGS <= 3:
#include <version.h>
#ifndef JAGS_MAJOR
#define JAGS_MAJOR 3
#endif // JAGS_MAJOR

// Which version of JAGS to assume:
// JAGS_MAJOR_FORCED will always be set - but may be just 0 or may be set from an environmental variable:
#if JAGS_MAJOR_FORCED > 0
#define JAGS_MAJOR_USED JAGS_MAJOR_FORCED
#else
#define JAGS_MAJOR_USED JAGS_MAJOR
#endif

// Will be defined on both platforms - if >0 (windows only) need to make sure that the version of ljags- passed matches what version.h says:
#if JAGS_MAJOR_ASSUMED > 0
#if JAGS_MAJOR_USED != JAGS_MAJOR_ASSUMED
#error "Error detecting the JAGS version - you need to set the 'JAGS_MAJOR_FORCED' environmental variable to the major version of JAGS you have installed"
#endif // JAGS_MAJOR_USED != JAGS_MAJOR_ASSUMED
#else
#define JAGS_MAJOR_ASSUMED 0
#endif // JAGS_MAJOR_ASSUMED

// Check version of JAGS is OK:
#if JAGS_MAJOR_USED > 4
#warning "Compiling against a later version of JAGS than has been tested for this version of RoBMA ... you should probably update the RoBMA package!"
#endif

#if JAGS_MAJOR_USED < 4
#error "This version of the RoBMA package requires compilation against JAGS version 4 or later"
#endif


#endif // JAGS_VERSIONS_H_
