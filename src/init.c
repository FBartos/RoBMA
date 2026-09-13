/*
    This file is based on file in the runjags package (version 2.0)
    The previous version of the file is Copyright (C) Matthew Denwood, licensed under GPL-2.
*/

#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h> // for SEXP

#include "glmm-aghq.h"

extern void getjagsversions(int *forced, int *assumed, int *detected, int *used);
extern SEXP RoBMA_selnorm_cache_control(SEXP capacity_bytes, SEXP clear);
extern SEXP RoBMA_selnorm_cache_snapshot(void);
extern SEXP RoBMA_selnorm_cache_restore(SEXP snapshots);
extern SEXP RoBMA_selnorm_sampler_control(SEXP enabled, SEXP settings, SEXP clear);
extern SEXP RoBMA_glmm_binom_marginal_loglik(SEXP ai, SEXP ci,
                                             SEXP n1i, SEXP n2i,
                                             SEXP mu_samples,
                                             SEXP tau_within,
                                             SEXP weights,
                                             SEXP theta_grid,
                                             SEXP log_theta_weights,
                                             SEXP logit_pi_grid,
                                             SEXP log_pi_weights);
extern SEXP RoBMA_glmm_binom_marginal_loglik_row_sum(SEXP ai, SEXP ci,
                                                     SEXP n1i, SEXP n2i,
                                                     SEXP mu_samples,
                                                     SEXP tau_within,
                                                     SEXP weights,
                                                     SEXP theta_grid,
                                                     SEXP log_theta_weights,
                                                     SEXP logit_pi_grid,
                                                     SEXP log_pi_weights);
extern SEXP RoBMA_glmm_pois_marginal_loglik(SEXP x1i, SEXP x2i,
                                            SEXP t1i, SEXP t2i,
                                            SEXP mu_samples,
                                            SEXP tau_within,
                                            SEXP weights,
                                            SEXP theta_grid,
                                            SEXP log_theta_weights,
                                            SEXP log_phi_grid,
                                            SEXP log_phi_weights);
extern SEXP RoBMA_glmm_pois_marginal_loglik_row_sum(SEXP x1i, SEXP x2i,
                                                    SEXP t1i, SEXP t2i,
                                                    SEXP mu_samples,
                                                    SEXP tau_within,
                                                    SEXP weights,
                                                    SEXP theta_grid,
                                                    SEXP log_theta_weights,
                                                    SEXP log_phi_grid,
                                                    SEXP log_phi_weights);
extern SEXP RoBMA_glmm_binom_conditional_loglik_sum(SEXP ai, SEXP ci,
                                                    SEXP n1i, SEXP n2i,
                                                    SEXP mu_samples,
                                                    SEXP logit_baserate,
                                                    SEXP weights);
extern SEXP RoBMA_glmm_pois_conditional_loglik_sum(SEXP x1i, SEXP x2i,
                                                   SEXP t1i, SEXP t2i,
                                                   SEXP mu_samples,
                                                   SEXP log_phi,
                                                   SEXP weights);
extern SEXP RoBMA_plot_normal_mixture_quantiles(SEXP mean, SEXP sd,
                                                SEXP probs, SEXP weights);
extern SEXP RoBMA_plot_selnorm_mixture_quantiles(SEXP mean, SEXP sd,
                                                 SEXP se, SEXP probs,
                                                 SEXP weights, SEXP selected,
                                                 SEXP omega, SEXP alpha,
                                                 SEXP phack_kind,
                                                 SEXP kernel_mode,
                                                 SEXP z_lower, SEXP z_upper,
                                                 SEXP sign, SEXP q,
                                                 SEXP phack_z_source,
                                                 SEXP phack_z_dest,
                                                 SEXP segment_bounds,
                                                 SEXP segment_step_bin,
                                                 SEXP segment_phack_region,
                                                 SEXP telescope_probabilities);
extern SEXP RoBMA_selnorm_kernel_loglik_matrix(SEXP yi, SEXP mu_num,
                                               SEXP sigma_num, SEXP mu_norm,
                                               SEXP sigma_norm, SEXP sei,
                                               SEXP weights, SEXP omega,
                                               SEXP alpha, SEXP phack_kind,
                                               SEXP kernel_mode,
                                               SEXP z_lower, SEXP z_upper,
                                               SEXP obs_bin, SEXP sign,
                                               SEXP q,
                                               SEXP phack_z_source,
                                                   SEXP phack_z_dest,
                                                   SEXP segment_bounds,
                                                   SEXP segment_step_bin,
                                                   SEXP segment_phack_region,
                                                   SEXP telescope_probabilities);
extern SEXP RoBMA_selnorm_mnorm_step_loglik_batch(
    SEXP yi, SEXP means, SEXP covariance_lower, SEXP sei, SEXP omega,
    SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
    SEXP telescope_probabilities, SEXP kernel_mode, SEXP qmc,
    SEXP points, SEXP scrambles, SEXP relative_tolerance,
    SEXP return_normalizer, SEXP vector_rule, SEXP normalizer_quadrature);
extern SEXP RoBMA_selnorm_sampling_conditioned_batch(
    SEXP yi, SEXP means, SEXP sampling_lower, SEXP diagonal, SEXP loading,
    SEXP rank, SEXP sampling_auxiliary, SEXP random_auxiliary, SEXP sei,
    SEXP omega, SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
    SEXP kernel_mode, SEXP telescope_probabilities, SEXP vector_rule,
    SEXP group_index, SEXP qmc, SEXP initial_points, SEXP max_points,
    SEXP scrambles, SEXP relative_tolerance,
    SEXP cluster_nodes, SEXP cluster_log_weights, SEXP cluster_orders,
    SEXP factor_nodes, SEXP factor_log_weights, SEXP factor_orders,
    SEXP factor_rule_counts);
extern SEXP RoBMA_selnorm_sampling_deletion_loglik_batch(
    SEXP yi, SEXP means, SEXP sampling_variances, SEXP integrated_variances,
    SEXP total_variances, SEXP sei, SEXP omega, SEXP z_lower, SEXP z_upper,
    SEXP obs_bin, SEXP kernel_mode, SEXP telescope_probabilities,
    SEXP nodes, SEXP log_weights, SEXP orders, SEXP relative_tolerance);
extern SEXP RoBMA_selnorm_conditioned_normalizer_batch(
    SEXP means, SEXP diagonal, SEXP loading, SEXP rank, SEXP sei, SEXP omega,
    SEXP z_lower, SEXP z_upper, SEXP sign, SEXP kernel_mode,
    SEXP telescope_probabilities, SEXP vector_rule, SEXP group_index, SEXP qmc,
    SEXP initial_points, SEXP max_points, SEXP scrambles, SEXP relative_tolerance,
    SEXP cluster_nodes, SEXP cluster_log_weights, SEXP cluster_orders,
    SEXP factor_nodes, SEXP factor_log_weights, SEXP factor_orders,
    SEXP factor_rule_counts);
extern SEXP RoBMA_selnorm_factor_projection_batch(
    SEXP means, SEXP diagonal, SEXP loading, SEXP rank, SEXP sei, SEXP omega,
    SEXP z_lower, SEXP z_upper, SEXP sign, SEXP kernel_mode,
    SEXP telescope_probabilities, SEXP vector_rule, SEXP group_index, SEXP qmc,
    SEXP initial_points, SEXP max_points, SEXP scrambles, SEXP relative_tolerance,
    SEXP cluster_nodes, SEXP cluster_log_weights, SEXP cluster_orders,
    SEXP factor_nodes, SEXP factor_log_weights, SEXP factor_orders,
    SEXP factor_rule_counts, SEXP z, SEXP kind);
extern SEXP RoBMA_selnorm_mnorm_zplot_batch(
    SEXP yi, SEXP means, SEXP covariance_lower, SEXP sei, SEXP omega,
    SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
    SEXP telescope_probabilities, SEXP kernel_mode, SEXP qmc,
    SEXP points, SEXP scrambles, SEXP relative_tolerance,
    SEXP z, SEXP probability, SEXP factors, SEXP vector_rule);
extern SEXP RoBMA_selnorm_cluster_step_loglik_batch(
    SEXP yi, SEXP means, SEXP residual_sd, SEXP loading, SEXP sei, SEXP omega,
    SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
    SEXP telescope_probabilities, SEXP kernel_mode, SEXP quadrature_nodes,
    SEXP quadrature_log_weights, SEXP quadrature_orders, SEXP qmc,
    SEXP initial_points, SEXP max_points, SEXP scrambles,
    SEXP relative_tolerance, SEXP return_normalizer, SEXP vector_rule);
extern SEXP RoBMA_selnorm_factor_step_loglik_batch(
    SEXP yi, SEXP means, SEXP residual_sd, SEXP loading, SEXP sei, SEXP omega,
    SEXP z_lower, SEXP z_upper, SEXP obs_bin, SEXP sign,
    SEXP telescope_probabilities, SEXP kernel_mode, SEXP quadrature_nodes,
    SEXP quadrature_log_weights, SEXP quadrature_orders,
    SEXP quadrature_rule_counts, SEXP qmc,
    SEXP initial_points, SEXP max_points, SEXP scrambles,
    SEXP relative_tolerance, SEXP return_normalizer, SEXP vector_rule);
extern SEXP RoBMA_selnorm_mnorm_step_rng_batch(
    SEXP means, SEXP covariance, SEXP sei, SEXP omega,
    SEXP z_lower, SEXP z_upper, SEXP sign, SEXP kernel_mode,
    SEXP dependency_blocks, SEXP max_attempts, SEXP vector_rule);
extern SEXP RoBMA_selnorm_gaussian_event_mass_batch(
    SEXP means, SEXP covariance_lower, SEXP sei, SEXP omega,
    SEXP z_lower, SEXP z_upper, SEXP sign, SEXP kernel_mode,
    SEXP vector_rule, SEXP lower, SEXP upper, SEXP qmc,
    SEXP points, SEXP scrambles, SEXP rank_one_loading,
    SEXP normalizer_quadrature, SEXP relative_tolerance);
extern SEXP RoBMA_norm_loglik_row_sum(SEXP yi, SEXP mu_samples,
                                      SEXP tau_within, SEXP sei,
                                      SEXP weights);
extern SEXP RoBMA_known_v_covariance_plan_create(SEXP y,
                                                  SEXP sampling_covariance,
                                                  SEXP random_covariance_factors,
                                                  SEXP block_indices);
extern SEXP RoBMA_known_v_covariance_plan_loglik(SEXP pointer, SEXP mean,
                                                  SEXP random_covariance_factors,
                                                  SEXP extra_variance);
extern SEXP RoBMA_known_v_covariance_plan_loglik_batch(
    SEXP pointer, SEXP means, SEXP random_covariance_states,
    SEXP extra_variances);
extern SEXP RoBMA_known_v_covariance_plan_group_iid_variance_grid_loglik(
    SEXP pointer, SEXP means, SEXP group_variances,
    SEXP diagonal_variances);
extern SEXP RoBMA_known_v_covariance_plan_affine_grid_loglik(
    SEXP pointer, SEXP means, SEXP base_covariances,
    SEXP update_covariances, SEXP reference_coefficient,
    SEXP coefficients);
extern SEXP RoBMA_known_v_covariance_plan_factor_grid_loglik(
    SEXP pointer, SEXP means, SEXP random_covariance_states,
    SEXP extra_variances, SEXP update_grid);
extern SEXP RoBMA_known_v_covariance_plan_location_quadratic_batch(
    SEXP pointer, SEXP means, SEXP bases,
    SEXP random_covariance_states, SEXP extra_variances);
extern SEXP RoBMA_known_v_covariance_plan_conditional_loglik_batch(
    SEXP pointer, SEXP means, SEXP random_covariance_states,
    SEXP extra_variances);
extern SEXP RoBMA_known_v_covariance_plan_conditional_loglik(
    SEXP pointer, SEXP mean, SEXP random_covariance_states,
    SEXP extra_variance);
extern SEXP RoBMA_known_v_covariance_plan_conditional_summary_batch(
    SEXP pointer, SEXP means, SEXP random_covariance_states,
    SEXP extra_variances);
extern SEXP RoBMA_known_v_covariance_plan_precision_residual_batch(
    SEXP pointer, SEXP means, SEXP random_covariance_states,
    SEXP extra_variances);
extern SEXP RoBMA_selnorm_kernel_loglik_row_sum(SEXP yi, SEXP mu_num,
                                                SEXP sigma_num,
                                                SEXP mu_norm,
                                                SEXP sigma_norm,
                                                SEXP sei, SEXP weights,
                                                SEXP omega, SEXP alpha,
                                                SEXP phack_kind,
                                                SEXP kernel_mode,
                                                SEXP z_lower,
                                                SEXP z_upper,
                                                SEXP obs_bin, SEXP sign,
                                                SEXP q,
                                                SEXP phack_z_source,
                                                SEXP phack_z_dest,
                                                SEXP segment_bounds,
                                                SEXP segment_step_bin,
                                                SEXP segment_phack_region,
                                                SEXP telescope_probabilities);
extern SEXP RoBMA_norm_cluster_loglik(SEXP yi, SEXP sei,
                                      SEXP mu_samples,
                                      SEXP tau_within,
                                      SEXP tau_between,
                                      SEXP cluster_index,
                                      SEXP cluster_size, SEXP weights,
                                      SEXP gamma_grid,
                                      SEXP log_gamma_weights);
extern SEXP RoBMA_norm_cluster_analytic_loglik(SEXP yi, SEXP vi,
                                               SEXP mu_samples,
                                               SEXP tau_within,
                                               SEXP tau_between,
                                               SEXP cluster_index,
                                               SEXP cluster_size);
extern SEXP RoBMA_norm_cluster_analytic_loglik_row_sum(SEXP yi, SEXP vi,
                                                       SEXP mu_samples,
                                                       SEXP tau_within,
                                                       SEXP tau_between,
                                                       SEXP cluster_index,
                                                       SEXP cluster_size);
extern SEXP RoBMA_norm_cluster_analytic_rho_grid_loglik(
    SEXP yi, SEXP vi, SEXP mu_samples, SEXP tau_total, SEXP rho,
    SEXP cluster_index, SEXP cluster_size);
extern SEXP RoBMA_norm_cluster_loglik_row_sum(SEXP yi, SEXP sei,
                                             SEXP mu_samples,
                                             SEXP tau_within,
                                             SEXP tau_between,
                                             SEXP cluster_index,
                                             SEXP cluster_size,
                                             SEXP weights,
                                             SEXP gamma_grid,
                                             SEXP log_gamma_weights);
extern SEXP RoBMA_selnorm_cluster_loglik(SEXP yi, SEXP sei,
                                         SEXP mu_samples,
                                         SEXP tau_within,
                                         SEXP tau_between,
                                         SEXP cluster_index,
                                         SEXP cluster_size, SEXP weights,
                                         SEXP gamma_grid,
                                         SEXP log_gamma_weights,
                                         SEXP omega, SEXP alpha,
                                         SEXP phack_kind,
                                         SEXP kernel_mode,
                                         SEXP z_lower, SEXP z_upper,
                                         SEXP obs_bin, SEXP sign, SEXP q,
                                         SEXP phack_z_source,
                                         SEXP phack_z_dest,
                                         SEXP segment_bounds,
                                         SEXP segment_step_bin,
                                         SEXP segment_phack_region,
                                         SEXP telescope_probabilities);
extern SEXP RoBMA_selnorm_cluster_loglik_row_sum(SEXP yi, SEXP sei,
                                                 SEXP mu_samples,
                                                 SEXP tau_within,
                                                 SEXP tau_between,
                                                 SEXP cluster_index,
                                                 SEXP cluster_size,
                                                 SEXP weights,
                                                 SEXP gamma_grid,
                                                 SEXP log_gamma_weights,
                                                 SEXP omega, SEXP alpha,
                                                 SEXP phack_kind,
                                                 SEXP kernel_mode,
                                                 SEXP z_lower,
                                                 SEXP z_upper,
                                                 SEXP obs_bin, SEXP sign,
                                                 SEXP q,
                                                 SEXP phack_z_source,
                                                 SEXP phack_z_dest,
                                                 SEXP segment_bounds,
                                                 SEXP segment_step_bin,
                                                 SEXP segment_phack_region,
                                                 SEXP telescope_probabilities);
extern SEXP RoBMA_selnorm_cluster_location_grid(SEXP yi, SEXP sei,
                                                SEXP mu_samples,
                                                SEXP mu_basis,
                                                SEXP current,
                                                SEXP values,
                                                SEXP tau_within,
                                                SEXP tau_between,
                                                SEXP cluster_index,
                                                SEXP cluster_size,
                                                SEXP weights,
                                                SEXP gamma_grid,
                                                SEXP log_gamma_weights,
                                                SEXP omega, SEXP alpha,
                                                SEXP phack_kind,
                                                SEXP kernel_mode,
                                                SEXP z_lower,
                                                SEXP z_upper,
                                                SEXP obs_bin, SEXP sign,
                                                SEXP q,
                                                SEXP phack_z_source,
                                                SEXP phack_z_dest,
                                                SEXP segment_bounds,
                                                SEXP segment_step_bin,
                                                SEXP segment_phack_region,
                                                SEXP telescope_probabilities);
extern SEXP RoBMA_selnorm_kernel_log_norm_matrix(SEXP mean, SEXP sd,
                                                 SEXP sei, SEXP omega,
                                                 SEXP alpha,
                                                 SEXP phack_kind,
                                                 SEXP kernel_mode,
                                                 SEXP z_lower,
                                                 SEXP z_upper,
                                                 SEXP sign, SEXP q,
                                                 SEXP phack_z_source,
                                                 SEXP phack_z_dest,
                                                 SEXP segment_bounds,
                                                 SEXP segment_step_bin,
                                                 SEXP segment_phack_region,
                                                 SEXP telescope_probabilities);
extern SEXP RoBMA_selnorm_kernel_log_norm_delta_grid(SEXP mean, SEXP sd,
                                                     SEXP basis,
                                                     SEXP current_log_norm,
                                                     SEXP current,
                                                     SEXP values, SEXP sei,
                                                     SEXP weights,
                                                     SEXP omega, SEXP alpha,
                                                     SEXP phack_kind,
                                                     SEXP kernel_mode,
                                                     SEXP z_lower,
                                                     SEXP z_upper,
                                                     SEXP sign, SEXP q,
                                                     SEXP phack_z_source,
                                                     SEXP phack_z_dest,
                                                     SEXP segment_bounds,
                                                     SEXP segment_step_bin,
                                                     SEXP segment_phack_region,
                                                     SEXP telescope_probabilities);
extern SEXP RoBMA_selnorm_kernel_cdf_matrix(SEXP x, SEXP mean, SEXP sd,
                                            SEXP sei, SEXP omega,
                                            SEXP alpha, SEXP phack_kind,
                                            SEXP kernel_mode,
                                            SEXP z_lower,
                                            SEXP z_upper, SEXP sign,
                                            SEXP q, SEXP phack_z_source,
                                            SEXP phack_z_dest,
                                            SEXP segment_bounds,
                                            SEXP segment_step_bin,
                                            SEXP segment_phack_region,
                                            SEXP lower_tail,
                                            SEXP telescope_probabilities);
extern SEXP RoBMA_selnorm_kernel_moments_matrix(SEXP mean, SEXP sd,
                                                SEXP sei, SEXP omega,
                                                SEXP alpha,
                                                SEXP phack_kind,
                                                SEXP kernel_mode,
                                                SEXP z_lower,
                                                SEXP z_upper,
                                                SEXP sign, SEXP q,
                                                SEXP phack_z_source,
                                                SEXP phack_z_dest,
                                                SEXP segment_bounds,
                                                SEXP segment_step_bin,
                                                SEXP segment_phack_region,
                                                SEXP telescope_probabilities);
extern SEXP RoBMA_selnorm_kernel_rng_matrix(SEXP mean, SEXP sd,
                                            SEXP sei, SEXP omega,
                                            SEXP alpha,
                                            SEXP phack_kind,
                                            SEXP kernel_mode,
                                            SEXP z_lower,
                                            SEXP z_upper,
                                            SEXP sign, SEXP q,
                                            SEXP phack_z_source,
                                            SEXP phack_z_dest,
                                            SEXP segment_bounds,
                                            SEXP segment_step_bin,
                                            SEXP segment_phack_region,
                                            SEXP telescope_probabilities);
extern SEXP RoBMA_selnorm_kernel_weighted_summary(SEXP yi, SEXP mean,
                                                  SEXP sd, SEXP sei,
                                                  SEXP psis_weights,
                                                  SEXP omega, SEXP alpha,
                                                  SEXP phack_kind,
                                                  SEXP kernel_mode,
                                                  SEXP z_lower,
                                                  SEXP z_upper,
                                                  SEXP sign, SEXP q,
                                                  SEXP phack_z_source,
                                                  SEXP phack_z_dest,
                                                  SEXP segment_bounds,
                                                  SEXP segment_step_bin,
                                                  SEXP segment_phack_region,
                                                  SEXP telescope_probabilities);
extern SEXP RoBMA_zcurve_normal_density_matrix(SEXP z_sequence,
                                               SEXP mean, SEXP sd,
                                               SEXP sei);
extern SEXP RoBMA_selnorm_zcurve_threshold_summary(SEXP z_threshold,
                                                   SEXP mean, SEXP sd,
                                                   SEXP sei, SEXP omega,
                                                   SEXP alpha,
                                                   SEXP phack_kind,
                                                   SEXP kernel_mode,
                                                   SEXP z_lower,
                                                   SEXP z_upper,
                                                   SEXP sign, SEXP q,
                                                   SEXP phack_z_source,
                                                   SEXP phack_z_dest,
                                                   SEXP segment_bounds,
                                                   SEXP segment_step_bin,
                                                   SEXP segment_phack_region,
                                                   SEXP extrapolate,
                                                   SEXP telescope_probabilities);
extern SEXP RoBMA_selnorm_zcurve_density_matrix(SEXP z_sequence,
                                                SEXP mean, SEXP sd,
                                                SEXP sei, SEXP omega,
                                                SEXP alpha,
                                                SEXP phack_kind,
                                                SEXP kernel_mode,
                                                SEXP z_lower,
                                                SEXP z_upper,
                                                SEXP sign, SEXP q,
                                                SEXP phack_z_source,
                                                SEXP phack_z_dest,
                                                SEXP segment_bounds,
                                                SEXP segment_step_bin,
                                                SEXP segment_phack_region,
                                                SEXP extrapolate,
                                                SEXP telescope_probabilities,
                                                SEXP latent_sd,
                                                SEXP quadrature);

static R_NativePrimitiveArgType getjagsversions_t[] = {
    INTSXP, INTSXP, INTSXP, INTSXP
};

static const R_CMethodDef cMethods[] = {
    {"getjagsversions", (DL_FUNC) &getjagsversions, 4, getjagsversions_t},
    {NULL, NULL, 0, NULL}
};

extern SEXP RoBMA_selnorm_covariance_envelope_components(SEXP, SEXP, SEXP);
extern SEXP RoBMA_selnorm_context_star_zplot(
  SEXP means, SEXP covariance_lower, SEXP selection_se, SEXP omega,
  SEXP lower, SEXP upper, SEXP sign, SEXP telescope, SEXP log_normalizers,
  SEXP context_log_weights, SEXP nodes, SEXP log_weights, SEXP z, SEXP allowance);


static const R_CallMethodDef callMethods[] = {
    {"RoBMA_selnorm_covariance_envelope_components", (DL_FUNC) &RoBMA_selnorm_covariance_envelope_components, 3},
    {"RoBMA_selnorm_context_star_zplot", (DL_FUNC) &RoBMA_selnorm_context_star_zplot, 14},
    {"RoBMA_selnorm_cache_control", (DL_FUNC) &RoBMA_selnorm_cache_control, 2},
    {"RoBMA_selnorm_cache_snapshot", (DL_FUNC) &RoBMA_selnorm_cache_snapshot, 0},
    {"RoBMA_selnorm_cache_restore", (DL_FUNC) &RoBMA_selnorm_cache_restore, 1},
    {"RoBMA_selnorm_sampler_control", (DL_FUNC) &RoBMA_selnorm_sampler_control, 3},
    {"RoBMA_glmm_binom_aghq", (DL_FUNC) &RoBMA_glmm_binom_aghq, 14},
    {"RoBMA_glmm_binom_aghq_row_sum", (DL_FUNC) &RoBMA_glmm_binom_aghq_row_sum, 14},
    {"RoBMA_glmm_pois_aghq", (DL_FUNC) &RoBMA_glmm_pois_aghq, 14},
    {"RoBMA_glmm_pois_aghq_row_sum", (DL_FUNC) &RoBMA_glmm_pois_aghq_row_sum, 14},
    {"RoBMA_glmm_binom_marginal_loglik", (DL_FUNC) &RoBMA_glmm_binom_marginal_loglik, 11},
    {"RoBMA_glmm_pois_marginal_loglik",  (DL_FUNC) &RoBMA_glmm_pois_marginal_loglik,  11},
    {"RoBMA_glmm_binom_marginal_loglik_row_sum", (DL_FUNC) &RoBMA_glmm_binom_marginal_loglik_row_sum, 11},
    {"RoBMA_glmm_pois_marginal_loglik_row_sum",  (DL_FUNC) &RoBMA_glmm_pois_marginal_loglik_row_sum,  11},
    {"RoBMA_glmm_binom_conditional_loglik_sum", (DL_FUNC) &RoBMA_glmm_binom_conditional_loglik_sum, 7},
    {"RoBMA_glmm_pois_conditional_loglik_sum",  (DL_FUNC) &RoBMA_glmm_pois_conditional_loglik_sum,  7},
    {"RoBMA_plot_normal_mixture_quantiles", (DL_FUNC) &RoBMA_plot_normal_mixture_quantiles, 4},
    {"RoBMA_plot_selnorm_mixture_quantiles", (DL_FUNC) &RoBMA_plot_selnorm_mixture_quantiles, 20},
    {"RoBMA_selnorm_kernel_loglik_matrix", (DL_FUNC) &RoBMA_selnorm_kernel_loglik_matrix, 22},
    {"RoBMA_selnorm_mnorm_step_loglik_batch", (DL_FUNC) &RoBMA_selnorm_mnorm_step_loglik_batch, 18},
    {"RoBMA_selnorm_sampling_conditioned_batch", (DL_FUNC) &RoBMA_selnorm_sampling_conditioned_batch, 30},
    {"RoBMA_selnorm_sampling_deletion_loglik_batch", (DL_FUNC) &RoBMA_selnorm_sampling_deletion_loglik_batch, 16},
    {"RoBMA_selnorm_conditioned_normalizer_batch", (DL_FUNC) &RoBMA_selnorm_conditioned_normalizer_batch, 25},
    {"RoBMA_selnorm_factor_projection_batch", (DL_FUNC) &RoBMA_selnorm_factor_projection_batch, 27},
    {"RoBMA_selnorm_mnorm_zplot_batch", (DL_FUNC) &RoBMA_selnorm_mnorm_zplot_batch, 19},
    {"RoBMA_selnorm_cluster_step_loglik_batch", (DL_FUNC) &RoBMA_selnorm_cluster_step_loglik_batch, 22},
    {"RoBMA_selnorm_factor_step_loglik_batch", (DL_FUNC) &RoBMA_selnorm_factor_step_loglik_batch, 23},
    {"RoBMA_selnorm_mnorm_step_rng_batch", (DL_FUNC) &RoBMA_selnorm_mnorm_step_rng_batch, 11},
    {"RoBMA_selnorm_gaussian_event_mass_batch", (DL_FUNC) &RoBMA_selnorm_gaussian_event_mass_batch, 17},
    {"RoBMA_norm_loglik_row_sum", (DL_FUNC) &RoBMA_norm_loglik_row_sum, 5},
    {"RoBMA_known_v_covariance_plan_create", (DL_FUNC) &RoBMA_known_v_covariance_plan_create, 4},
    {"RoBMA_known_v_covariance_plan_loglik", (DL_FUNC) &RoBMA_known_v_covariance_plan_loglik, 4},
    {"RoBMA_known_v_covariance_plan_loglik_batch", (DL_FUNC) &RoBMA_known_v_covariance_plan_loglik_batch, 4},
    {"RoBMA_known_v_covariance_plan_group_iid_variance_grid_loglik", (DL_FUNC) &RoBMA_known_v_covariance_plan_group_iid_variance_grid_loglik, 4},
    {"RoBMA_known_v_covariance_plan_affine_grid_loglik", (DL_FUNC) &RoBMA_known_v_covariance_plan_affine_grid_loglik, 6},
    {"RoBMA_known_v_covariance_plan_factor_grid_loglik", (DL_FUNC) &RoBMA_known_v_covariance_plan_factor_grid_loglik, 5},
    {"RoBMA_known_v_covariance_plan_location_quadratic_batch", (DL_FUNC) &RoBMA_known_v_covariance_plan_location_quadratic_batch, 5},
    {"RoBMA_known_v_covariance_plan_conditional_loglik", (DL_FUNC) &RoBMA_known_v_covariance_plan_conditional_loglik, 4},
    {"RoBMA_known_v_covariance_plan_conditional_loglik_batch", (DL_FUNC) &RoBMA_known_v_covariance_plan_conditional_loglik_batch, 4},
    {"RoBMA_known_v_covariance_plan_conditional_summary_batch", (DL_FUNC) &RoBMA_known_v_covariance_plan_conditional_summary_batch, 4},
    {"RoBMA_known_v_covariance_plan_precision_residual_batch", (DL_FUNC) &RoBMA_known_v_covariance_plan_precision_residual_batch, 4},
    {"RoBMA_selnorm_kernel_loglik_row_sum", (DL_FUNC) &RoBMA_selnorm_kernel_loglik_row_sum, 22},
    {"RoBMA_norm_cluster_loglik", (DL_FUNC) &RoBMA_norm_cluster_loglik, 10},
    {"RoBMA_norm_cluster_loglik_row_sum", (DL_FUNC) &RoBMA_norm_cluster_loglik_row_sum, 10},
    {"RoBMA_norm_cluster_analytic_loglik", (DL_FUNC) &RoBMA_norm_cluster_analytic_loglik, 7},
    {"RoBMA_norm_cluster_analytic_loglik_row_sum", (DL_FUNC) &RoBMA_norm_cluster_analytic_loglik_row_sum, 7},
    {"RoBMA_norm_cluster_analytic_rho_grid_loglik", (DL_FUNC) &RoBMA_norm_cluster_analytic_rho_grid_loglik, 7},
    {"RoBMA_selnorm_cluster_loglik", (DL_FUNC) &RoBMA_selnorm_cluster_loglik, 25},
    {"RoBMA_selnorm_cluster_loglik_row_sum", (DL_FUNC) &RoBMA_selnorm_cluster_loglik_row_sum, 25},
    {"RoBMA_selnorm_cluster_location_grid", (DL_FUNC) &RoBMA_selnorm_cluster_location_grid, 28},
    {"RoBMA_selnorm_kernel_log_norm_matrix", (DL_FUNC) &RoBMA_selnorm_kernel_log_norm_matrix, 17},
    {"RoBMA_selnorm_kernel_log_norm_delta_grid", (DL_FUNC) &RoBMA_selnorm_kernel_log_norm_delta_grid, 22},
    {"RoBMA_selnorm_kernel_cdf_matrix", (DL_FUNC) &RoBMA_selnorm_kernel_cdf_matrix, 19},
    {"RoBMA_selnorm_kernel_moments_matrix", (DL_FUNC) &RoBMA_selnorm_kernel_moments_matrix, 17},
    {"RoBMA_selnorm_kernel_rng_matrix", (DL_FUNC) &RoBMA_selnorm_kernel_rng_matrix, 17},
    {"RoBMA_selnorm_kernel_weighted_summary", (DL_FUNC) &RoBMA_selnorm_kernel_weighted_summary, 19},
    {"RoBMA_zcurve_normal_density_matrix", (DL_FUNC) &RoBMA_zcurve_normal_density_matrix, 4},
    {"RoBMA_selnorm_zcurve_threshold_summary", (DL_FUNC) &RoBMA_selnorm_zcurve_threshold_summary, 19},
    {"RoBMA_selnorm_zcurve_density_matrix", (DL_FUNC) &RoBMA_selnorm_zcurve_density_matrix, 21},
    {NULL, NULL, 0}
};

void
R_init_RoBMA(DllInfo *dll)
{
    R_registerRoutines(dll, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
