#include <math>

double protonConfLevel_pass(float beta_meas, float beta_true, float p, std::vector<double> fit_beta_mean, std::vector<double> fit_beta_sig ){
  
  double beta_mean = fit_beta_mean[0]/sqrt( (1.0 + pow(fit_beta_mean[1]/p,2) ) );
  double beta_sig = fit_beta_sig[0] + fit_beta_sig[1]/sqrt(p);
  
  double cl = ( 1 - erf( fabs( beta_meas - beta_mean )/sigma/sqrt(2.0) ) );

  return cl;


}
