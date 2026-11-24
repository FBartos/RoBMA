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
// flip and its reverse; 2 flips without the reverse proposal density, in both
// directions; 3 no final correction). Returns the previous value.
int location_block_test_mutant(int mutant);
// The outcome of the last factory call that found a declaration: "claimed",
// "claimed without latents" or the first check that declined.
std::string location_block_test_decision();
// The probe (frozen into samplers built from now on): every sweep first
// dry-runs each collapsed Gaussian move and records, per move, its type
// (0 flip on, 1 flip off, 2 refresh, 3 study effect, 4 swap, 5 joint), the
// label (1 selection), the stage-1 log ratio, the full log ratio from the
// direct Gaussian evaluation of priors, children and proposal densities, the
// selection part log A~(old) - log A~(new), and the collapsed log odds.
// probe < 0 reads the flag; otherwise sets it and returns the previous one.
bool location_block_test_probe(int probe);
std::vector<double> location_block_test_take_probe();
// Test hook: switches frozen into the samplers built from now on; a negative
// argument keeps a setting. refresh: the kernels that refresh the active
// coefficients in selection labels (0 the default, 1 the per-group refreshes
// R_j only, 2 one elliptical slice step E, 3 two E steps, 4 E and the R_j);
// fit_rule: the rule of the Newton-Laplace fits (0 the coarsest rung with the
// acceptance safeguard, 1 the coarsest rung always, 2 the sampler's rule);
// timing: 1 accumulates the wall time of the update's parts in the report.
// Returns the previous three settings.
std::vector<int> location_block_test_config(int refresh, int fit_rule, int timing);

}
}
#endif
