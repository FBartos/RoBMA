#ifndef ROBMA_LOCATION_BLOCK_REPORT_H
#define ROBMA_LOCATION_BLOCK_REPORT_H

// The per-fit report and the test hook of the location block sampler
// (LocationBlock.h), without JAGS headers, for the R-native entry points.

#include <string>
#include <vector>

namespace jags {
namespace RoBMA {

// The per-fit counters of the location block samplers of this process
// (relaxed atomics), read and cleared when clear is true: values in the order
// of location_block_report_names().
std::vector<char const *> const &location_block_report_names();
std::vector<double> location_block_report(bool clear);

// Test hook, not package interface: a deliberately broken kernel frozen into
// the samplers built from now on (0 none; 1 a shear that differs between a
// move and its reverse; 2 dormant slab draws from a law other than the prior
// that the reverse densities assume, and a collapsed flip without its reverse
// density; 3 no final correction). Returns the previous value.
int location_block_test_mutant(int mutant);
// The outcome of the last factory call that found a declaration: "claimed",
// "claimed without latents" or the first check that declined.
std::string location_block_test_decision();

}
}
#endif
