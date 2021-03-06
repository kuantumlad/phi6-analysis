#include <vector>
#include <iostream> 
#include <map>
#include <cmath>

std::map<int,int> hadronMLEID(int hypothesis, int gpart, int e_index[], Int_t q[], Float_t p[], UChar_t sc_sect[], UChar_t dc_sect[], Float_t sc_t[], Float_t sc_r[], UChar_t sc_pd[], int pip_vvp_strict, int pip_R1fid_strict, int pip_MXcut_strict, int ExpOrSim, Float_t ec_ei[], UChar_t ec_sect[], UChar_t cc_sect[], UShort_t nphe[], Float_t ec_eo[], Float_t cx[], Float_t cy[], Float_t cz[], Float_t b[], Float_t tl1_x[], Float_t tl1_y[],  TLorentzVector V4_H[], int currentrunno, int pim_vvp_strict, int pim_R1fid_strict, int pim_MXcut_strict, std::map<int, std::vector<double> > pr_mean_fit, std::map<int, std::vector<double> > pip_mean_fit, std::map<int, std::vector<double> > kp_mean_fit, std::map<int, std::vector<double> > prot_sigma_fit, std::map<int, std::vector<double> > pip_sigma_fit, std::map<int, std::vector<double> > kp_sigma_fit, double pr_conf, double pip_conf, double kp_conf, double pr_anticonf, double pip_anticonf, double kp_anticonf )
{
  Float_t speed_of_light = 29.9792458; // cm/ns
  Float_t pip_mass = 0.13957; // GeV
  Float_t pim_mass = 0.13957; // GeV
  Float_t prot_mass = 0.938272; // GeV
  Float_t kp_mass = 0.4937; // GeV
  Float_t km_mass = 0.4937; // GeV

  Float_t pi = 3.14159265359;
  Float_t pi180 = pi/180.0;

  std::vector<int> hadrons;
  double temp_p_pr = -1.0;
  double temp_p_pip = -1.0;
  double temp_p_kp = -1.0;

  int temp_index_pr = -1;
  int temp_index_pip = -1;
  int temp_index_kp = -1;

  //std::cout << " >>GOING TO LOOP OVER GPART IN HADRONS " << std::endl;
  //std::cout << " GPART " << gpart << std::endl;
  for(int k = 0; k < gpart; k++) // loop over particles
    {
      float cand_time = -1.0;
      float beta_measured_pim = -1.0;
      float beta_measured_km = -1.0;
	
      float beta_measured_pr = -1.0;
      float beta_measured_kp = -1.0;
      float beta_measured_pip = -1.0;
	
      float beta_measured = -1.0;

      int sc_sector = sc_sect[k]; //from 0 to 6
      //std::cout << " SECTOR " << sc_sector << std::endl;
      float corrStartTime = e_sctimeCorr(ExpOrSim, sc_t[e_index[1]], dc_sect[e_index[1]], sc_pd[e_index[1]], currentrunno) - (sc_r[e_index[1]]/speed_of_light);
      float ehtcorrMeasuredTime = h_sctimeCorr(ExpOrSim, sc_t[k], dc_sect[k], sc_pd[k], currentrunno) - corrStartTime;
      beta_measured = (sc_r[k]/ehtcorrMeasuredTime)/speed_of_light;
	  	
      if( sc_sector >= 1 && q[k] > 0 && k != e_index[1] ){
	//std::cout << " sector " << sc_sector << std::endl;
	
	float CalcBeta_pip = p[k]/sqrt(pow(p[k],2) + pow(pip_mass,2));
	float CalcBeta_pim = p[k]/sqrt(pow(p[k],2) + pow(pim_mass,2));
	float CalcBeta_kp = p[k]/sqrt(pow(p[k],2) + pow(kp_mass,2));
	float CalcBeta_km = p[k]/sqrt(pow(p[k],2) + pow(km_mass,2));
	float CalcBeta_pr = p[k]/sqrt(pow(p[k],2) + pow(prot_mass,2));

	std::vector<float> calc_beta_pos;
	calc_beta_pos.push_back(CalcBeta_pip);
	calc_beta_pos.push_back(CalcBeta_kp);
	calc_beta_pos.push_back(CalcBeta_pr);
	
	std::vector<float> calc_beta_neg;
	calc_beta_neg.push_back(CalcBeta_pim);
	calc_beta_neg.push_back(CalcBeta_km);

	//mle fit values
	//std::cout << " >> GETTING FIT VALUES " << std::endl;
	//std:: cout << " >> PROTON MEAN MAP SIZE"  << pr_mean_fit.size()  << std::endl;
	std::vector<double> prot_mean_a = pr_mean_fit[0];//->second;      
	std::vector<double> prot_mean_b = pr_mean_fit[1];//)->second;      
	std::vector<double> prot_mean_c = pr_mean_fit[2];//second;      

	//std:: cout << " >> PIP MEAN MAP SIZE"  << pip_mean_fit.size()  << std::endl;	
	std::vector<double> pip_mean_a = pip_mean_fit[0];//->second;      
	std::vector<double> pip_mean_b = pip_mean_fit[1];//->second;      
	std::vector<double> pip_mean_c = pip_mean_fit[2];//->second;      

	std::vector<double> Kp_mean_a = kp_mean_fit[0];//->second;      
	std::vector<double> Kp_mean_b = kp_mean_fit[1];//->second;      
	std::vector<double> Kp_mean_c = kp_mean_fit[2];//->second;      

	//std::cout << " proton sigma size " << prot_sigma_fit.size() << std::endl;
	std::vector<double> prot_sigma_a = prot_sigma_fit[0];//->second;      
	std::vector<double> prot_sigma_b = prot_sigma_fit[1];//.find(1)->second;      
	std::vector<double> prot_sigma_c = prot_sigma_fit[2];//.find(2)->second;      
	
	std::vector<double> pip_sigma_a = pip_sigma_fit[0];//.find(0)->second;      
	std::vector<double> pip_sigma_b = pip_sigma_fit[1];//.find(1)->second;      
	std::vector<double> pip_sigma_c = pip_sigma_fit[2];//.find(2)->second;      

	std::vector<double> Kp_sigma_a = kp_sigma_fit[0];//.find(0)->second;      
	std::vector<double> Kp_sigma_b = kp_sigma_fit[1];//.find(1)->second;      
	std::vector<double> Kp_sigma_c = kp_sigma_fit[2];//.find(2)->second;      

	double mom = p[k];

	//std::cout << " PRINTING SIZE " << std::endl;
	//std::cout << " SIZE " << prot_mean_a.size() << " " << prot_mean_a[sc_sector-1] << " " << prot_mean_b[sc_sector-1] << " " << prot_mean_c[sc_sector-1] << std::endl;

	//	double mean_prot = prot_mean[0] + prot_mean[1]*mom/sqrt(mom*mom+prot_mean[2]) + prot_mean[3]*mom*mom/sqrt(mom*mom*mom*mom+prot_mean[4]);
	//	double sigma_prot = prot_sigma[0] + prot_sigma[1]/sqrt(mom) + prot_sigma[2]*exp(-mom);
	double mean_prot =  CalcBeta_pr + prot_mean_a[sc_sector-1]*pow(mom,2) + prot_mean_b[sc_sector-1]*mom + prot_mean_c[sc_sector-1];
	double sigma_prot =  prot_sigma_a[sc_sector-1]*pow(mom,2) + prot_sigma_b[sc_sector-1]*mom + prot_sigma_c[sc_sector-1];
	double prob_prot = (1/(sigma_prot*sqrt(2*3.14159))) * exp(-0.5 * pow((beta_measured - mean_prot)/sigma_prot, 2));
	double conf_prot = (1.0 - TMath::Erf(fabs(beta_measured - mean_prot)/sigma_prot/sqrt(2.0))); 

	//double mean_pip = pip_mean[0] + pip_mean[1]*mom/sqrt(mom*mom+pip_mean[2]) + pip_mean[3]*mom*mom/sqrt(mom*mom*mom*mom+pip_mean[4]);
	//double sigma_pip = pip_sigma[0] + pip_sigma[1]/sqrt(mom) + pip_sigma[2]*exp(-mom);
	
	double mean_pip =  CalcBeta_pip + pip_mean_a[sc_sector-1]*pow(mom,2) + pip_mean_b[sc_sector-1]*mom + pip_mean_c[sc_sector-1];
	double sigma_pip =  pip_sigma_a[sc_sector-1]*pow(mom,2) + pip_sigma_b[sc_sector-1]*mom + pip_sigma_c[sc_sector-1];
	double prob_pip = (1/(sigma_pip*sqrt(2*3.14159))) * exp(-0.5 * pow((beta_measured - mean_pip)/sigma_pip, 2));
	double conf_pip = (1.0 - TMath::Erf(fabs(beta_measured - mean_pip)/sigma_pip/sqrt(2.0))); 

	//double mean_Kp = Kp_mean[0] + Kp_mean[1]*mom/sqrt(mom*mom+Kp_mean[2]) + Kp_mean[3]*mom*mom/sqrt(mom*mom*mom*mom+Kp_mean[4]);
	//double sigma_Kp = Kp_sigma[0] + Kp_sigma[1]/sqrt(mom) + Kp_sigma[2]*exp(-mom);
	
	double mean_Kp =  CalcBeta_kp + Kp_mean_a[sc_sector-1]*pow(mom,2) + Kp_mean_b[sc_sector-1]*mom + Kp_mean_c[sc_sector-1];
	double sigma_Kp =  Kp_sigma_a[sc_sector-1]*pow(mom,2) + Kp_sigma_b[sc_sector-1]*mom + Kp_sigma_c[sc_sector-1];
	double prob_Kp = (1/(sigma_Kp*sqrt(2*3.14159))) * exp(-0.5 * pow((beta_measured - mean_Kp)/sigma_Kp, 2));
	double conf_Kp = (1.0 - TMath::Erf(fabs(beta_measured - mean_Kp)/sigma_Kp/sqrt(2.0)));        

	//std::cout << "PROB " << prob_prot << " " << prob_pip << " " << prob_Kp << std::endl;
	//std::cout << "CONF " << conf_prot << " " << conf_pip << " " << conf_Kp << std::endl;
	
	if(prob_prot > prob_pip  && prob_prot > prob_Kp  && hypothesis == 2212 && conf_prot > pr_conf  && conf_pip < pr_anticonf  && conf_Kp < pr_anticonf){ 
	  
	  //cout << "proton  -  probability: " << prob_prot << "     confidence: " << conf_prot << endl;//"     ( mom: " << mom << " beta:" << CalcBeta_pip << " mean: " << mean_prot << " sigma: " << sigma_prot <<  " )" << endl;
	  //cout << "pip     -  probability: " << prob_pip  << "     confidence: " << conf_pip  << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip  << " sigma: " << sigma_pip  <<  " )" << endl;
	  //cout << "Kp      -  probability: " << prob_Kp   << "     confidence: " << conf_Kp   << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Kp   << " sigma: " << sigma_Kp   <<  " )" << endl;
	  //cout << endl;
	  

	  //keep only the largest momentum candidates
	  if( p[k] > temp_p_pr ){
	    temp_p_pr = p[k];
	    temp_index_pr = k;
	  }
	  //return true;
        }
        if(prob_pip > prob_prot && prob_pip  > prob_Kp  && hypothesis == 211  && conf_pip > pip_conf && conf_prot < pip_anticonf  && conf_Kp  < pip_anticonf){ 
	  //cout << "proton  -  probability: " << prob_prot << "     confidence: " << conf_prot << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << " sigma: " << sigma_prot <<  " )" << endl;
	  //cout << "pip     -  probability: " << prob_pip  << "     confidence: " << conf_pip  << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip  << " sigma: " << sigma_pip  <<  " )" << endl;
	  //cout << "Kp      -  probability: " << prob_Kp   << "     confidence: " << conf_Kp   << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Kp   << " sigma: " << sigma_Kp   <<  " )" << endl;
	  //cout << endl;

 	  if( p[k] > temp_p_pip ){
	    temp_p_pip = p[k];
	    temp_index_pip = k;
	  }
	  //return true;
        }
        if(prob_Kp > prob_prot && prob_Kp > prob_pip && hypothesis == 321  && conf_Kp > kp_conf  && conf_prot < pr_anticonf  && conf_pip < pr_anticonf ){ 
	  ///cout << "proton  -  probability: " << prob_prot << "     confidence: " << conf_prot << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_prot << " sigma: " << sigma_prot <<  " )" << endl;
	  ///cout << "pip     -  probability: " << prob_pip  << "     confidence: " << conf_pip  << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_pip  << " sigma: " << sigma_pip  <<  " )" << endl;
	  //cout << "Kp      -  probability: " << prob_Kp   << "     confidence: " << conf_Kp   << endl;//"     ( mom: " << mom << " beta:" << beta << " mean: " << mean_Kp   << " sigma: " << sigma_Kp   <<  " )" << endl;
	  //cout << endl;
	  if( p[k] > temp_p_kp ){
	    temp_p_kp = p[k];
	    temp_index_kp = k;
	  }
	  //return true;
	}

      }
    }

  std::map<int, int> hadron_candidates;
  hadron_candidates[2212] = temp_index_pr;
  hadron_candidates[211] = temp_index_pip;
  hadron_candidates[321] = temp_index_kp;

  return hadron_candidates;

}

