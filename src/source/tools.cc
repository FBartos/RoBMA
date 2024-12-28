#include "tools.h"
#include <JRmath.h>
#include <cmath>
#include <algorithm>
#include <numeric> 
#include <vector>
#include "mnorm.h"

// assigns log of weights product corresponding to the outcome x based on a one/two-sided weightfunction
double log_weight_onesided(double const *x, double const *crit_x, double const *omega, const int J)
{
	double w = -68;

	if(*x >= *(crit_x + J - 2)){
	// return the last omega if x > last crit_x
		w = *(omega + J - 1);
	}else if(*x < *crit_x){
	// return the first omega if x < first crit_x
		w = *omega;
	}else{
	// check the remaining cutpoints sequentially
		for(int i = 1; i < J; i++){
			if( (*x >= *(crit_x + i - 1)) && (*x < *(crit_x + i)) ){
				w = *(omega + i);
				break;
			}
		}
	}

	return std::log(w);
}

double log_weight_twosided(double const *x, double const *crit_x, double const *omega, const int J)
{
	double w = -68;
	double abs_x = std::fabs(*x);

	if(abs_x >= *(crit_x + J - 2)){
	// return the last omega if abs(x) > last crit_x
		w = *(omega + J - 1);
	}else if(abs_x < *crit_x){
	// return the first omega if abs(x) < first crit_x
		w = *omega;
	}else{
	// check the remaining cutpoints sequentially
		for(int i = 1; i < J; i++){
			if( (abs_x >= *(crit_x + i - 1)) && (abs_x < *(crit_x + i)) ){
				w = *(omega + i);
				break;
			}
		}
	}
	return std::log(w);
}

// compute the standardizing constant for a one/two-sided weightfunctions
double log_std_constant_onesided(double const *x, double const *const_mu, double const *sigma, double const *crit_x, double const *omega, const int J) {

    // Pre-allocate the denoms vector with size J
    std::vector<double> denoms(J, 0.0);

    // Compute the first denominator
    denoms[0] = pnorm(crit_x[0], *const_mu, *sigma, true, false);
    denoms[0] = (denoms[0] < 0.0) ? 0.0 : denoms[0];
    double denom_sum = denoms[0];

    // Compute the middle denominators
    for (int j = 1; j < J - 1; ++j) {
        denoms[j] = pnorm(crit_x[j], *const_mu, *sigma, true, false) - denom_sum;
        denoms[j] = (denoms[j] < 0.0) ? 0.0 : denoms[j];
        denom_sum += denoms[j];
    }

    // Compute the last denominator
    if (J > 1) { // Ensure there is more than one element
        denoms[J - 1] = 1.0 - denom_sum;
        denoms[J - 1] = (denoms[J - 1] < 0.0) ? 0.0 : denoms[J - 1];
    }

    // Pre-allocate log_terms vector with maximum possible size J
    std::vector<double> log_terms;
    log_terms.reserve(J);

    // Populate log_terms with log(denoms[k]) + log(omega[k]) where denoms[k] and omega[k] > 0
    for (int k = 0; k < J; ++k) {
        if (denoms[k] > 0.0 && omega[k] > 0.0) {
            log_terms.push_back(std::log(denoms[k]) + std::log(omega[k]));
        }
    }

    // If no valid terms, return negative infinity
    if (log_terms.empty()) {
        return -INFINITY;
    }

    // Find the maximum log_term to use in the log-sum-exp trick
    double max_log = *std::max_element(log_terms.begin(), log_terms.end());

    // Compute the sum of exponentials of (log_terms - max_log)
    double sum_exp = 0.0;
    for (const auto& log_val : log_terms) {
        sum_exp += std::exp(log_val - max_log);
    }

    // Compute log_sum_exp
    double log_sum_exp = max_log + std::log(sum_exp);

    return log_sum_exp;
}

double log_std_constant_twosided(double const *x, double const *const_mu, double const *sigma, double const *crit_x, double const *omega, const int J) {

    // Pre-allocate the denoms vector with size J
    std::vector<double> denoms(J, 0.0);

    // Compute the first denominator: P(X < crit_x[0]) - P(X < -crit_x[0])
    denoms[0] = pnorm(crit_x[0], *const_mu, *sigma, true, false) - pnorm(-crit_x[0], *const_mu, *sigma, true, false);
    denoms[0] = (denoms[0] < 0.0) ? 0.0 : denoms[0];
    double denom_sum = denoms[0];

    // Compute the middle denominators: P(X < crit_x[j]) - P(X < -crit_x[j]) - denom_sum
    for (int j = 1; j < J - 1; ++j) {
        denoms[j] = pnorm(crit_x[j], *const_mu, *sigma, true, false) - pnorm(-crit_x[j], *const_mu, *sigma, true, false) - denom_sum;
        denoms[j] = (denoms[j] < 0.0) ? 0.0 : denoms[j];
        denom_sum += denoms[j];
    }

    // Compute the last denominator: 1.0 - denom_sum
    if (J > 1) { // Ensure there is more than one element
        denoms[J - 1] = 1.0 - denom_sum;
        denoms[J - 1] = (denoms[J - 1] < 0.0) ? 0.0 : denoms[J - 1];
    }

    // Pre-allocate log_terms vector with maximum possible size J
    std::vector<double> log_terms;
    log_terms.reserve(J);

    // Populate log_terms with log(denoms[k]) + log(omega[k]) where denoms[k] and omega[k] > 0
    for (int k = 0; k < J; ++k) {
        if (denoms[k] > 0.0 && omega[k] > 0.0) {
            log_terms.push_back(std::log(denoms[k]) + std::log(omega[k]));
        }
    }

    // If no valid terms, return negative infinity
    if (log_terms.empty()) {
        return -INFINITY;
    }

    // Find the maximum log_term to use in the log-sum-exp trick
    double max_log = *std::max_element(log_terms.begin(), log_terms.end());

    // Compute the sum of exponentials of (log_terms[k] - max_log)
    double sum_exp = 0.0;
    for (const auto& log_val : log_terms) {
        sum_exp += std::exp(log_val - max_log);
    }

    // Compute log_sum_exp
    double log_sum_exp = max_log + std::log(sum_exp);

    return log_sum_exp;
}

double log_std_m_constant_onesided(double const *x, double const *const_mu, double const *sigma, double const *crit_x, double const *omega ,const int K, const int J)
{
  double std_constant = 0;

  // create dynamically allocated arrays for covariance matrix decomposition (and non-constat mu)
  double * sigma_stdev;
  double * sigma_corr;
  double * mu;

  sigma_stdev = new double [K];
  sigma_corr  = new double [K * (K - 1) / 2];
  mu          = new double [K];

  for(int k = 0; k < K; k++){
    *(sigma_stdev + k) = std::sqrt( *(sigma+ K * k + k) );
    *(mu + k)          = *(const_mu + k);
  }
  for(int k = 0; k < K; k++){
    for(int j = 0; j < K && j < k; j++){
      *(sigma_corr + k * (k - 1) / 2 + j) = *(sigma + K * k + j) / std::sqrt( *(sigma + K * k + k) * *(sigma + K * j + j));
    }
  }

  // create dynamically allocated arrays for the current lower bounds, upper bounds and  infinity indicator
  double * temp_lower;
  double * temp_upper;
  int    * temp_infin;

  temp_lower = new double [K];
  temp_upper = new double [K];
  temp_infin = new int    [K];

  // current weight holder
  double temp_log_weight;

  // indexes tracker for weights assignment
  int * index_weights;
  index_weights = new int [K];
  for(int k = 0; k < K; k++){
    *(index_weights + k) = 0;
  }

  // iterate over all subspaces created by the cutpoints
  for(int i = 0; i < std::pow(J, K); i++){

  	// reset the current weight
  	temp_log_weight = 0;

    // assign proper upper and lower cutopoints & weight
    for(int k = 0; k < K; k++){

      // weight
      temp_log_weight += std::log(*(omega + *(index_weights + k)));

      // the upper and lower bounds
      if(*(index_weights + k) == 0){
        // the lower bound is - infinity for the first cutpoint
        *(temp_lower + k) = 0;
        *(temp_upper + k) = *(crit_x + k * (J - 1) + *(index_weights + k));
        *(temp_infin + k) = 0; // 0 for lower bound =  -inf
      }else if(*(index_weights + k) == J - 1){
        // the upper bound is + infinity for the last cutpoint
        *(temp_lower + k) = *(crit_x + k * (J - 1) + *(index_weights + k) - 1);
        *(temp_upper + k) = 0;
        *(temp_infin + k) = 1; // 1 for upper bound = +inf
      }else{
        *(temp_lower + k) = *(crit_x + k * (J - 1) + *(index_weights + k) - 1);
        *(temp_upper + k) = *(crit_x + k * (J - 1) + *(index_weights + k));
        *(temp_infin + k) = 2; // 2 for neither of the bound is infinite
      }
    }

  	// get the current weighted probability
  	double temp_prob = cpp_mnorm_cdf(&temp_lower[0], &temp_upper[0], &temp_infin[0], &mu[0], &sigma_stdev[0], &sigma_corr[0], K);
  	// check and skip negative numbers due to numerical impressions for very low probabilities
  	if(temp_prob > 0){
  	  std_constant += std::exp(std::log(temp_prob) + temp_log_weight);
  	}

	  // increase index for the next iteration
	  if((i + 1) < std::pow(J, K)){
	    increase_index( &index_weights[0], K-1, J-1);
	  }
  }

  // clean the memory
  delete[] sigma_stdev;
  delete[] sigma_corr;
  delete[] mu;
  delete[] temp_lower;
  delete[] temp_upper;
  delete[] temp_infin;
  delete[] index_weights;

  return std::log(std_constant);
}

double log_std_m_constant_twosided(double const *x, double const *const_mu, double const *sigma, double const *crit_x, double const *omega ,const int K, const int J)
{
  double log_std_constant = 0;

  // turn the two-sided selection into one-sided selection
  double * omega_onesided;
  double * crit_x_onesided;

  omega_onesided  = new double [2 * J - 1];
  crit_x_onesided = new double [2 * (J - 1) * K];

  for(int i = 0; i < 2 * J - 1; i++){
    if(i < J){
        *(omega_onesided + i) = *(omega + (J - 1) - i);
    }else{
      *(omega_onesided + i) = *(omega + i - (J - 1));
    }
  }

  for(int k = 0; k < K; k++){
    for(int j = 0; j < 2 * (J - 1); j++){
      if(j < (J - 1)){
          *(crit_x_onesided + k * 2 * (J - 1) + j) = - *(crit_x + k * (J - 1) + (J - 2) - j);
        }else{
          *(crit_x_onesided + k * 2 * (J - 1) + j) =   *(crit_x + k * (J - 1) + j - (J - 1));
        }
      }
  }

  log_std_constant = log_std_m_constant_onesided(&x[0], &const_mu[0], &sigma[0], &crit_x_onesided[0], &omega_onesided[0], K, 2 * J - 1);

  // clean the memory
  delete[] omega_onesided;
  delete[] crit_x_onesided;

  return log_std_constant;
}

//// additional helper functions
// increases index carrying the cutoff coordinates by one
void increase_index(int *index, int i, int max_i)
{
  if( *(index + i) == (max_i)){
    *(index + i) = 0;
    increase_index(&index[0], i-1, max_i);
  }else{
    *(index + i) += 1;
  }
}

// parsing vector input into different pieces
double * extract_x_v(const double *x_v, int indx_start, int K)
{
  double * this_x;
  this_x = new double [K];

  // cout << "x = ";
  for(int k = 0; k < K; k++){
    *(this_x + k) = *(x_v  + indx_start + k);
    // cout << *(this_x + k) << "\t";
  }

  return &this_x[0];
}

double * extract_mu_v(const double *mu_v,  int indx_start, int K)
{
  double * this_mu;
  this_mu = new double [K];

  // cout << "mu = ";
  for(int k = 0; k < K; k++){
    *(this_mu + k) = *(mu_v + indx_start + k);
    // cout << *(this_mu + k) << "\t";
  }

  return &this_mu[0];
}

double * extract_sigma_v(const double *se2_v, const double *tau2, double cov, int indx_start, int K)
{
  double * this_sigma;
  this_sigma  = new double [K * K];

  for(int k = 0; k < K; k++){
    // cout << "sigma[" << k << ",] = ";
    for(int j = 0; j < K; j++){
      if(k == j){
        *(this_sigma + k * K + j) = *(se2_v + indx_start + k) + *tau2;
      }else{
        *(this_sigma + k * K + j) = cov;
      }
      // cout << *(this_sigma + k * K + j) << "\t";
    }
  }

  return &this_sigma[0];
}

double * extract_crit_x_v(const double *crit_x_v, int indx_start, int K, const int J)
{
  double * this_crit_x;
  this_crit_x  = new double [K * (J-1)];

  for(int k = 0; k < K; k++){
    // cout << "crit_x[" << k << ",] = ";
    for(int j = 0; j < J - 1 ; j++){
      *(this_crit_x + k * (J - 1) + j) = *(crit_x_v + (indx_start + k) * (J - 1) + j);
      // cout << *(this_crit_x + k * (J-1) + j) << "\t";
    }
  }

  return &this_crit_x[0];
}

double ddirichlet(const std::vector<double>& x, const std::vector<double>& alpha){

  int k = x.size();

  double prod_gamma = 0.0, sum_alpha = 0.0, p_tmp = 0.0, beta_const = 0.0;

  for (int j = 0; j < k; j++) {
    double alpha_val = alpha[j];
    double x_val = x[j];
    sum_alpha += alpha_val;
    prod_gamma += std::lgamma(alpha_val);
    p_tmp += std::log(x_val) * (alpha_val - 1.0);
  }

  beta_const = prod_gamma - std::lgamma(sum_alpha);
  double p = p_tmp - beta_const;

  return p;
}
