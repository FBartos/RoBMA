#include <Rinternals.h>
#include <R_ext/Error.h>
#include <R_ext/Random.h>

#include <Rmath.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <limits>
#include <vector>

#include "plot-root.h"
#include "selnorm/selnorm.h"

extern "C" double Rf_dnorm4(double, double, double, int);

// Private implementation fragments share one anonymous namespace and are kept in
// this order so helper symbols stay internal to this translation unit.
namespace {
#include "r-selnorm-common.cc.inc"
#include "r-selnorm-funnel-common.cc.inc"
#include "r-selnorm-cluster.cc.inc"
}

#include "r-selnorm-loglik.cc.inc"
#include "r-selnorm-kernel.cc.inc"
#include "r-selnorm-funnel-zcurve.cc.inc"
