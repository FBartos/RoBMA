#include "transformations_common.h"
#include <util/nainf.h>
#include <JRmath.h>
#include <cmath>

using std::exp;
using std::sqrt;
using std::pow;
using std::log;

double pi(){
  return 3.14159265358979323846;
}

// based on Borenstein, M., Hedges, L. V., Higgins, J. P., & Rothstein, H. R. (2011). Introduction to meta-analysis. John Wiley & Sons.
// main effect size transformation functions
double cpp_d2r        (double d){
  return d / sqrt(pow(d, 2) + 4);
}
double cpp_d2logOR    (double d){
  return d * pi() / sqrt(3);
}
double cpp_r2d        (double r){
  return 2 * r / sqrt( 1 - pow(r, 2) );
}
double cpp_r2z        (double r){
  return 0.5 * log( (1 + r) / (1 - r) );
}
double cpp_logOR2d    (double logOR){
  return logOR * sqrt(3)/pi();
}
double cpp_z2r        (double z){
  return (exp(2 * z) - 1) / (1 + exp(2 * z));
}

// composite functions
double cpp_d2z        (double d){
  return cpp_r2z(cpp_d2r(d));
}
double cpp_r2logOR    (double r){
  return cpp_d2logOR(cpp_r2d(r));
}
double cpp_logOR2z    (double logOR){
  return cpp_d2z(cpp_logOR2d(logOR));
}
double cpp_logOR2r    (double logOR){
  return cpp_d2r(cpp_logOR2d(logOR));
}
double cpp_z2d        (double z){
  return cpp_r2d(cpp_z2r(z));
}
double cpp_z2logOR    (double z){
  return cpp_d2logOR(cpp_z2d(z));
}

// main standard error transformation functions
double cpp_se_d2se_logOR  (double se_d){
  return sqrt(pow(se_d, 2) * pow(pi(), 2) / 3);
}
double cpp_se_d2se_r      (double se_d, double d){
  return sqrt( (pow(4, 2) * pow(se_d, 2)) / pow( pow(d, 2) + 4, 3) );
}
double cpp_se_r2se_d      (double se_r, double r){
  return sqrt( (4 * pow(se_r, 2)) / pow( 1 - pow(r, 2), 3) );
}
double cpp_se_logOR2se_d  (double se_logOR){
  return sqrt( pow(se_logOR, 2) * 3 / pow(pi(), 2) );
}

// composite standard error transformation functions
double cpp_se_d2se_z      (double se_d, double d){
  double n = cpp_n_d(d, se_d);
  return cpp_se_z(n);
}
double cpp_se_r2se_z      (double se_r, double r){
  double n = cpp_n_r(r, se_r);
  return cpp_se_z(n);
}
double cpp_se_r2se_logOR  (double se_r, double r){
  double se_d = cpp_se_r2se_d(se_r, r);
  return cpp_se_d2se_logOR(se_d);
}
double cpp_se_logOR2se_r  (double se_logOR, double logOR){
  double se_d = cpp_se_logOR2se_d(se_logOR);
  double d    = cpp_logOR2d(logOR);
  return cpp_se_d2se_r(se_d, d);
}
double cpp_se_logOR2se_z  (double se_logOR, double logOR){
  double se_d = cpp_se_logOR2se_d(se_logOR);
  double d    = cpp_logOR2d(logOR);
  return cpp_se_d2se_z(se_d, d);
}
double cpp_se_z2se_d      (double se_z, double z){
  double n = cpp_n_z(se_z);
  double d = cpp_z2d(z);
  return cpp_se_d(d, n);
}
double cpp_se_z2se_r      (double se_z, double z){
  double r = cpp_z2r(z);
  double n = cpp_n_z(se_z);
  return cpp_se_r(r, n);
}
double cpp_se_z2se_logOR  (double se_z, double z){
  double se_d = cpp_se_z2se_d(se_z, z);
  return cpp_se_d2se_logOR(se_d);
}

// linear scaling function
double cpp_scale_d2logOR (double d){
  return d * (pi()/sqrt(3));
}
double cpp_scale_d2z     (double d){
  return d / 2;
}
double cpp_scale_logOR2d (double logOR){
  return logOR / (pi()/sqrt(3));
}
double cpp_scale_z2d     (double z){
  return z * 2;
}

// composite standard error transformation functions
double cpp_scale_z2logOR (double z){
  return cpp_scale_d2logOR(cpp_scale_z2d(z));
}
double cpp_scale_logOR2z (double logOR){
  return cpp_scale_d2z(cpp_scale_logOR2d(logOR));
}

// helper functions
double cpp_n_d            (double d, double se){
  return (pow(d,2) + 8) / (2 * pow(se,2));
}
double cpp_n_r            (double r, double se){
  // (r^4 - 2*r^2 + se^2 + 1)/se^2 : according to the Borenstein, however, it is not consistent with the remaining
  double d    = cpp_r2d(r);
  double se_d = cpp_se_r2se_d(se, r);
  return cpp_n_d(d, se_d);
}
double cpp_n_z            (double se){
  return (3 * pow(se,2) + 1) / pow(se,2);
}
double cpp_se_d           (double d, double n){
  return sqrt( 4 / n + pow(d,2) / ( 2 * n));
}
double cpp_se_r           (double r, double n){
  // sqrt((1-r^2)^2/(n-1)) : according to the Borenstein, however, it is not consistent with the remaining transformations
  double d    = cpp_r2d(r);
  double se_d = cpp_se_d(d, n);
  return cpp_se_d2se_r(se_d, d);
}
double cpp_se_z           (double n){
  return sqrt( 1 / (n - 3));
}
