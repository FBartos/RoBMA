#include <Rinternals.h>
#include <R_ext/Error.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "glmm-binomial-loglik.h"

// Keep the likelihood kernels in one translation unit while separating the
// shared validation, marginal, and clustered entry-point implementations.
#include "r-glmm-common.cc.inc"
#include "r-glmm-marginal.cc.inc"
#include "r-glmm-cluster.cc.inc"
