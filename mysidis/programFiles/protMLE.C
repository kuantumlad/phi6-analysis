#include <vector.h>
#include <math.h>

double protonMLE_pass(float beta_meas, float beta_true, float p, std::vector<double> fit_beta_mean, std::vector<double> fit_beta_sig ){


  double beta_mean = fit_beta_mean[0]/sqrt( (1.0 + pow(fit_beta_mean[1]/p,2) ) );
  double beta_sig = fit_beta_sig[0] + fit_beta_sig[1]/sqrt(p);


  double residual = beta_meas - beta_mean;
  double residual_sig = (beta_meas - beta_mean) / beta_sig;
  double exponential = exp( 0.5 * pow( residual_sig, 2 ) );

  double likelihood = (1.0/sqrt(2.0 * 3.14159265358 * beta_sig*beta_sig) ) * exp( -0.5 * pow( (beta_meas - beta_mean )/beta_sigma, 2 ));


  return likelihood;

}
