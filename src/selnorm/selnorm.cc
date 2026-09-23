#include "selnorm.h"
#include "selnorm-tail.h"

#include <JRmath.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

extern "C" double Rf_dnorm4(double, double, double, int);

namespace {

// Internal helper fragments deliberately share this anonymous namespace.
#include "selnorm-probability.cc.inc"
#include "selnorm-step.cc.inc"
#include "selnorm-phack.cc.inc"

}

#include "selnorm-boundary.cc.inc"
#include "selnorm-api.cc.inc"
