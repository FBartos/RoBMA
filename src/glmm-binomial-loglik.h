#ifndef ROBMA_GLMM_BINOMIAL_LOGLIK_H
#define ROBMA_GLMM_BINOMIAL_LOGLIK_H

// Evaluate count * log(probability), including the exact zero-count limit.
inline double glmm_count_log_probability(int count, double log_probability)
{
  return count == 0 ? 0.0 :
    static_cast<double>(count) * log_probability;
}


// Evaluate the two-sample binomial log likelihood from paired log tails.
inline double glmm_binomial_log_likelihood(
  int events_1, int events_2, int total_1, int total_2,
  double log_coefficient,
  double log_probability_1, double log_complement_1,
  double log_probability_2, double log_complement_2)
{
  return
    log_coefficient +
    glmm_count_log_probability(events_1, log_probability_1) +
    glmm_count_log_probability(total_1 - events_1, log_complement_1) +
    glmm_count_log_probability(events_2, log_probability_2) +
    glmm_count_log_probability(total_2 - events_2, log_complement_2);
}

#endif
