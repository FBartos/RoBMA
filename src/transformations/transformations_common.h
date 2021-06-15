#ifndef transformations_common_H_
#define transformations_common_H_


double pi ();

// main effect size transformation functions
double cpp_d2r     (double d);
double cpp_d2logOR (double d);
double cpp_r2d     (double r);
double cpp_r2z     (double r);
double cpp_logOR2d (double logOR);
double cpp_z2r     (double z);

// composite effect size transformation functions
double cpp_d2z        (double d);
double cpp_r2logOR    (double r);
double cpp_logOR2z    (double logOR);
double cpp_logOR2r    (double logOR);
double cpp_z2d        (double z);
double cpp_z2logOR    (double z);

// main standard error transformation functions
double cpp_se_d2se_logOR  (double se_d);
double cpp_se_d2se_r      (double se_d, double d);
double cpp_se_r2se_d      (double se_r, double r);
double cpp_se_logOR2se_d  (double se_logOR);

// composite standard error transformation functions
double cpp_se_d2se_z      (double se_d, double d);
double cpp_se_r2se_z      (double se_r, double r);
double cpp_se_r2se_logOR  (double se_r, double r);
double cpp_se_logOR2se_r  (double se_logOR, double logOR);
double cpp_se_logOR2se_z  (double se_logOR, double logOR);
double cpp_se_z2se_d      (double se_z, double z);
double cpp_se_z2se_r      (double se_z, double z);
double cpp_se_z2se_logOR  (double se_z, double z);

// linear scaling function (not great, but help as aproximations for setting priors, especially when transformation formula break for too high variances)
double cpp_scale_d2logOR (double d);
double cpp_scale_d2z     (double d);
double cpp_scale_logOR2d (double logOR);
double cpp_scale_z2d     (double z);

// composite standard error transformation functions
double cpp_scale_z2logOR (double z);
double cpp_scale_logOR2z (double logOR);

// helper functions
double cpp_se_d           (double d, double n);
double cpp_se_r           (double r, double n);
double cpp_se_z           (double n);
double cpp_n_d            (double d, double se);
double cpp_n_r            (double r, double se);
double cpp_n_z            (double se);
#endif /* transformations_common_H_ */
