#ifndef ROBMA_GLMM_AGHQ_H
#define ROBMA_GLMM_AGHQ_H

#include <Rinternals.h>

#ifdef __cplusplus
extern "C" {
#endif

SEXP RoBMA_glmm_binom_aghq(SEXP ai, SEXP ci, SEXP n1i, SEXP n2i,
                            SEXP mu_samples, SEXP tau_within, SEXP weights,
                            SEXP alpha, SEXP beta, SEXP nodes,
                            SEXP log_weights, SEXP tolerance,
                            SEXP consecutive, SEXP mode_tolerance);

SEXP RoBMA_glmm_binom_aghq_row_sum(SEXP ai, SEXP ci, SEXP n1i, SEXP n2i,
                                    SEXP mu_samples, SEXP tau_within,
                                    SEXP weights, SEXP alpha, SEXP beta,
                                    SEXP nodes, SEXP log_weights,
                                    SEXP tolerance, SEXP consecutive,
                                    SEXP mode_tolerance);

SEXP RoBMA_glmm_pois_aghq(SEXP x1i, SEXP x2i, SEXP t1i, SEXP t2i,
                           SEXP mu_samples, SEXP tau_within, SEXP weights,
                           SEXP mean, SEXP sd, SEXP nodes,
                           SEXP log_weights, SEXP tolerance,
                           SEXP consecutive, SEXP mode_tolerance);

SEXP RoBMA_glmm_pois_aghq_row_sum(SEXP x1i, SEXP x2i, SEXP t1i, SEXP t2i,
                                   SEXP mu_samples, SEXP tau_within,
                                   SEXP weights, SEXP mean, SEXP sd,
                                   SEXP nodes, SEXP log_weights,
                                   SEXP tolerance, SEXP consecutive,
                                   SEXP mode_tolerance);

#ifdef __cplusplus
}
#endif

#endif
