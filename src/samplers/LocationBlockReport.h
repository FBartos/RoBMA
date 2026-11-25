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

// The evidence record of the location gates of this process since the last
// clearing read, over the updates after adaptation: per gate (its declaration
// slot) and chain the rows of location_block_evidence_names() (the updates
// and the gate's transitions in each direction between each update's entry
// and the state it retains; of the inner stage, before the exact correction,
// the flips proposed and accepted in each direction, proposals the free
// screen removed, and the count, sum and maximum of the log collapsed odds of
// the other state against the current one for the proposals that reached
// stage 1).
std::vector<char const *> const &location_block_evidence_names();
std::vector<double> location_block_evidence(bool clear);

// The evidence record of an owned bias label (LocationBlockBias.h) since the
// last clearing read, over the updates after adaptation: per chain the rows
// of location_bias_evidence_names() (proposals and acceptances of the label
// update L0, the Gaussian-family move L1, the slope update S and the swap L2
// with its stage-1 passes; per direction of L2, leaving the selection or the
// Gaussian family, the count, sum and maximum of the stage-1 log ratio of
// every attempt and of the swaps the sweep's exact correction kept), and in
// `labels` the rows (chain, label, sweeps ending in the label).
std::vector<char const *> const &location_bias_evidence_names();
std::vector<double> location_bias_evidence(bool clear, std::vector<double> *labels);

// Test hook, not package interface: a deliberately broken kernel frozen into
// the samplers built from now on (0 none; 1 a shear that differs between a
// flip and its reverse; 2 flips without the reverse proposal density, in both
// directions; 3 no final correction; 4 the elliptical slice step around the
// slab prior means instead of the Gaussian working mean; 5 the elliptical
// slice step without the latent shear; with an owned bias label, of its swap
// L2: 6 the label Hastings factor omitted, 7 the gate share omitted, 8 the
// prior density of the slope that becomes dormant omitted, 9 the nested
// weight node's density omitted from the acceptance ratio (the recorded
// inner-target change keeps it), 10 the within-label contrast instead of the
// absolute evidence of a selection label; 11 the joint refresh's slot taken
// from the entry label; 12 the design column of the last coordinate dropped
// after the factory's probe, an untracked mean change; 13 the slice refresh's
// coordinates in the forward order in reversed sweeps; 14 one error thrown
// after the exact correction's write of the final state, in the first sweep
// after adaptation that changes the state). Returns the previous value.
int location_block_test_mutant(int mutant);
// The outcome of the last factory call that found a declaration: "claimed",
// "claimed without latents" or the first check that declined.
std::string location_block_test_decision();
// The probe (frozen into samplers built from now on): every sweep first
// dry-runs each collapsed Gaussian move and records, per move, its type
// (0 flip on, 1 flip off, 2 refresh, 3 study effect, 4 swap, 5 joint), the
// label (1 selection), the stage-1 log ratio, the full log ratio from the
// direct Gaussian evaluation of priors, children and proposal densities, the
// selection part log A~(old) - log A~(new), and the collapsed log odds. It
// also dry-runs one elliptical slice proposal (type 6) and records log T as
// the block computes it, the separately evaluated change of the inner target
// minus that of the Gaussian working model, the inner target's change and the
// working model's change (the first two must agree). With an owned bias label
// each bias component records at its slot (LocationBlockBiasMoves.h): type 7
// L1 (analytic and computed log ratio, prior-cell flags, label), type 8 S
// (its acceptance log ratio), type 9 L2 (kind 0: stage 1, stage 1 plus the
// reverse stage 1, full, full plus the reverse full; kind 1: differences of
// E_0, of the evidences and label proposals and of the shears rebuilt from
// the proposal, and of the accepted flows of the two directions; kind 2: the
// fitted Gaussian plane's largest error and scale; kind 3: stage 1 and the
// full ratio plus their reverses from the forward move's base quantities;
// kind 4: relative differences of G0 and of the plane rebuilt from the
// proposal, the largest difference of the refitted Laplace means in SDs and
// of their precisions, relative) and type 10 L0 (the
// largest difference of a label's buffered exact score from the graph's and
// the scale).
// probe < 0 reads the flag; otherwise sets it and returns the previous one.
bool location_block_test_probe(int probe);
std::vector<double> location_block_test_take_probe();
// Test hook: switches frozen into the samplers built from now on; a negative
// argument keeps a setting. refresh: the kernels that refresh the active
// coefficients in selection labels (0 the default, 1 the per-group refreshes
// R_j only, 2 one elliptical slice step E, 3 two E steps, 4 E and the R_j,
// 5 the R_j with every group's slice update forced from the first adaptOff());
// fit_rule: the rule of the Newton-Laplace fits (0 the coarsest rung with the
// acceptance safeguard, 1 the coarsest rung always, 2 the sampler's rule);
// timing: 1 accumulates the wall time of the update's parts in the report;
// swap: the schedule of the bias label's cross-family swap L2 (0 the
// default, k >= 1 attempted with probability 1 / k per sweep); trace (bits):
// 1 records every L2 attempt (location_block_test_take_trace), 2 every
// sweep's sequence (location_block_test_take_sequence).
// Returns the previous five settings.
std::vector<int> location_block_test_config(int refresh, int fit_rule, int timing, int swap = -1, int trace = -1);
// The per-attempt record of L2 since the last call (trace switch): per
// attempt chain, source label, target label, stage 1 passed, stage-1 log
// ratio, gate drawn, stage-2 log ratio, accepted, the first latent column
// entry (the random-effect SD) and the target label's first six weights.
std::vector<double> location_block_test_take_trace();
// The sweep sequences since the last call (trace switch, bit 2): per sweep a
// header (-1, chain, reversed, entry label is a selection label, entries) and
// its entries (kind, a, b) in the order they ran: the slots (kind 0 flip,
// 1 refresh, 2 swap, 3 joint, 4 ellipse, 5 latent, 6 bias component; a and b
// the group, pair, latent or bias component) and within a composite slice
// refresh its coordinates (kind 100, group, coordinate).
std::vector<double> location_block_test_take_sequence();

}
}
#endif
