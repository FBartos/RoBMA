#include "get_weight.h"
#include <JRmath.h>
#include <cmath>

using std::log;
using std::fabs;

double log_weight_onesided(double const *x, double const *crit_x, double const *omega, const int J)
{
	double w;

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

	return log(w);
}

double log_weight_twosided(double const *x, double const *crit_x, double const *omega, const int J)
{
	double w;
	double abs_x = fabs(*x);

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
	return log(w);
}


