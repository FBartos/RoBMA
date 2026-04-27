/*
    This file is based on file in the runjags package (version 2.0)
    The previous version of the file is Copyright (C) Matthew Denwood, licensed under GPL-2.
*/

#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h> // for SEXP

extern void getjagsversions(int *forced, int *assumed, int *detected, int *used);
extern SEXP RoBMA_wnorm_mix_logpdf_matrix(SEXP yi, SEXP mean, SEXP sd,
                                          SEXP omega, SEXP crit_yi,
                                          SEXP bias_indicator,
                                          SEXP crit_yi_mapping,
                                          SEXP crit_yi_mapping_max);
extern SEXP RoBMA_wnorm_mix_cdf_matrix(SEXP q, SEXP mean, SEXP sd,
                                       SEXP omega, SEXP crit_yi,
                                       SEXP bias_indicator,
                                       SEXP crit_yi_mapping,
                                       SEXP crit_yi_mapping_max,
                                       SEXP lower_tail, SEXP log_p);
extern SEXP RoBMA_wnorm_mix_moments_matrix(SEXP mean, SEXP sd,
                                           SEXP omega, SEXP crit_yi,
                                           SEXP bias_indicator,
                                           SEXP crit_yi_mapping,
                                           SEXP crit_yi_mapping_max);
extern SEXP RoBMA_wnorm_mix_log_norm_matrix(SEXP mean, SEXP sd,
                                            SEXP omega, SEXP crit_yi,
                                            SEXP bias_indicator,
                                            SEXP crit_yi_mapping,
                                            SEXP crit_yi_mapping_max);
extern SEXP RoBMA_wnorm_mix_logpdf_precomp_matrix(SEXP yi, SEXP mean,
                                                  SEXP sd, SEXP omega,
                                                  SEXP crit_yi,
                                                  SEXP bias_indicator,
                                                  SEXP crit_yi_mapping,
                                                  SEXP crit_yi_mapping_max,
                                                  SEXP log_norm);
extern SEXP RoBMA_wnorm_mix_cdf_precomp_matrix(SEXP q, SEXP mean,
                                               SEXP sd, SEXP omega,
                                               SEXP crit_yi,
                                               SEXP bias_indicator,
                                               SEXP crit_yi_mapping,
                                               SEXP crit_yi_mapping_max,
                                               SEXP log_norm,
                                               SEXP lower_tail,
                                               SEXP log_p);
extern SEXP RoBMA_glmm_binom_marginal_loglik(SEXP ai, SEXP ci,
                                             SEXP n1i, SEXP n2i,
                                             SEXP mu_samples,
                                             SEXP tau_within,
                                             SEXP theta_grid,
                                             SEXP log_theta_weights,
                                             SEXP logit_pi_grid,
                                             SEXP log_pi_weights);
extern SEXP RoBMA_glmm_pois_marginal_loglik(SEXP x1i, SEXP x2i,
                                            SEXP t1i, SEXP t2i,
                                            SEXP mu_samples,
                                            SEXP tau_within,
                                            SEXP theta_grid,
                                            SEXP log_theta_weights,
                                            SEXP log_phi_grid,
                                            SEXP log_phi_weights);
extern SEXP RoBMA_glmm_binom_cluster_loglik(SEXP ai, SEXP ci,
                                            SEXP n1i, SEXP n2i,
                                            SEXP mu_samples,
                                            SEXP tau_within,
                                            SEXP tau_between,
                                            SEXP cluster_index,
                                            SEXP cluster_size,
                                            SEXP weights,
                                            SEXP theta_grid,
                                            SEXP log_theta_weights,
                                            SEXP logit_pi_grid,
                                            SEXP log_pi_weights,
                                            SEXP gamma_grid,
                                            SEXP log_gamma_weights);
extern SEXP RoBMA_glmm_pois_cluster_loglik(SEXP x1i, SEXP x2i,
                                           SEXP t1i, SEXP t2i,
                                           SEXP mu_samples,
                                           SEXP tau_within,
                                           SEXP tau_between,
                                           SEXP cluster_index,
                                           SEXP cluster_size,
                                           SEXP weights,
                                           SEXP theta_grid,
                                           SEXP log_theta_weights,
                                           SEXP log_phi_grid,
                                           SEXP log_phi_weights,
                                           SEXP gamma_grid,
                                           SEXP log_gamma_weights);
extern SEXP RoBMA_regplot_normal_mixture_interval(SEXP mean, SEXP sd,
                                                  SEXP probs);
extern SEXP RoBMA_regplot_weighted_mixture_interval(SEXP mean, SEXP sd,
                                                    SEXP probs,
                                                    SEXP omega,
                                                    SEXP bias_indicator,
                                                    SEXP is_weightfunction,
                                                    SEXP crit_yi,
                                                    SEXP crit_yi_mapping,
                                                    SEXP crit_yi_mapping_max,
                                                    SEXP effect_direction);

static R_NativePrimitiveArgType getjagsversions_t[] = {
    INTSXP, INTSXP, INTSXP, INTSXP
};

static const R_CMethodDef cMethods[] = {
    {"getjagsversions", (DL_FUNC) &getjagsversions, 4, getjagsversions_t},
    {NULL, NULL, 0, NULL}
};

static const R_CallMethodDef callMethods[] = {
    {"RoBMA_wnorm_mix_logpdf_matrix",  (DL_FUNC) &RoBMA_wnorm_mix_logpdf_matrix,  8},
    {"RoBMA_wnorm_mix_cdf_matrix",     (DL_FUNC) &RoBMA_wnorm_mix_cdf_matrix,     10},
    {"RoBMA_wnorm_mix_moments_matrix", (DL_FUNC) &RoBMA_wnorm_mix_moments_matrix, 7},
    {"RoBMA_wnorm_mix_log_norm_matrix", (DL_FUNC) &RoBMA_wnorm_mix_log_norm_matrix, 7},
    {"RoBMA_wnorm_mix_logpdf_precomp_matrix", (DL_FUNC) &RoBMA_wnorm_mix_logpdf_precomp_matrix, 9},
    {"RoBMA_wnorm_mix_cdf_precomp_matrix", (DL_FUNC) &RoBMA_wnorm_mix_cdf_precomp_matrix, 11},
    {"RoBMA_glmm_binom_marginal_loglik", (DL_FUNC) &RoBMA_glmm_binom_marginal_loglik, 10},
    {"RoBMA_glmm_pois_marginal_loglik",  (DL_FUNC) &RoBMA_glmm_pois_marginal_loglik,  10},
    {"RoBMA_glmm_binom_cluster_loglik", (DL_FUNC) &RoBMA_glmm_binom_cluster_loglik, 16},
    {"RoBMA_glmm_pois_cluster_loglik",  (DL_FUNC) &RoBMA_glmm_pois_cluster_loglik,  16},
    {"RoBMA_regplot_normal_mixture_interval", (DL_FUNC) &RoBMA_regplot_normal_mixture_interval, 3},
    {"RoBMA_regplot_weighted_mixture_interval", (DL_FUNC) &RoBMA_regplot_weighted_mixture_interval, 10},
    {NULL, NULL, 0}
};

void
R_init_RoBMA(DllInfo *dll)
{
    R_registerRoutines(dll, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
