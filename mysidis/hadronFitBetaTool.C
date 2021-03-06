#include <fstream>
#include <iostream>
#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "programFiles/functions.C"
#include "programFiles/eID.C"
#include "programFiles/hadronID.C"
#include "programFiles/getGenIndices.C"
#include "MomCorr.C"
#include "CalculatorTools.C"

int hadronFitBetaTool( int filestart = 1, int Nfiles = 1){

  TStopwatch *stopwatch = new TStopwatch();
  
  TFile *fIn = new TFile("hadron_beta_out.root");  

  Float_t e_mass = 0.000511; // GeV
  Float_t prot_mass = 0.938272; // GeV
  Float_t pip_mass = 0.13957; // GeV
  Float_t pim_mass = 0.13957; // GeV
  Float_t kp_mass = 0.4937; // GeV
  Float_t km_mass = 0.4937; // GeV
  Float_t speed_of_light = 29.9792458; // cm/ns
  Float_t Beam_Energy = 5.498; // GeV
  Float_t pi = 3.14159265359;
  Float_t pi180 = pi/180.0;
  Float_t pi180_inv = 180.0/pi;
  int ExpOrSim = 1;


  std::vector<TH2F*> h2_betap_pr;
  std::vector<TH2F*> h2_betap_pip;
  std::vector<TH2F*> h2_betap_kp;

  std::vector<TH2F*> h2_betap_pim;
  std::vector<TH2F*> h2_betap_km;

  //create vectors of beta vs p histograms
  for( int s = 0; s <= 6; s++ ){
    //positive hadrons
    h2_betap_pr.push_back( (TH2F*)fIn->Get(Form("h_beta_all_pr_sector_%d", s)) );
    h2_betap_pip.push_back( (TH2F*)fIn->Get(Form("h_beta_all_pip_sector_%d", s)) );
    h2_betap_kp.push_back( (TH2F*)fIn->Get(Form("h_beta_all_kp_sector_%d", s)) );

    //negative hadrons
    h2_betap_pim.push_back( (TH2F*)fIn->Get(Form("h_beta_all_pim_sector_%d", s)) );
    h2_betap_km.push_back( (TH2F*)fIn->Get(Form("h_beta_all_km_sector_%d", s)) );
  }

  std::vector<TGraphErrors*> h_beta_mean_pr;
  std::vector<TGraphErrors*> h_beta_mean_pip;
  std::vector<TGraphErrors*> h_beta_mean_kp;

  std::vector<TGraphErrors*> h_beta_mean_km;
  std::vector<TGraphErrors*> h_beta_mean_pim;  
  
  std::vector<TGraphErrors*> h_beta_sig_pr;
  std::vector<TGraphErrors*> h_beta_sig_kp;
  std::vector<TGraphErrors*> h_beta_sig_pip;

  //create vector of tgraphs for means and sig to fit 
  //fits will be used for the maximum likelihood estimators

  std::vector<TCanvas*> v_can;

  for( int s =  0; s <= 6; s++ ){
    v_can.push_back( new TCanvas(Form("c_%d",s), Form("c_%d",s), 800, 800 ));
    v_can[s]->Divide(5,5);

  }

  for( int i = 0; i <= 6; i++ ){

    int nbins_x_pr = h2_betap_pr[i]->GetNbinsX();
    int nbins_x_pip = h2_betap_pip[i]->GetNbinsX();
    int nbins_x_kp = h2_betap_kp[i]->GetNbinsX();

    int ndiv_pr = 10;
    int ndiv_pip = 10;
    int ndiv_kp = 10;

    int nfits_pr = nbins_x_pr / ndiv_pr;
    int nfits_pip = nbins_x_pip / ndiv_pip;
    int nfits_kp = nbins_x_kp / ndiv_kp;

    int start_bin_pr = 0;
    int start_bin_pip = 0;
    int start_bin_kp = 0;

    int end_bin_pr = ndiv_pr;
    int end_bin_pip = ndiv_pip;
    int end_bin_kp = ndiv_kp;

    //fill with the fit mean and sigma for the tgraphs

    std::vector<double> temp_center_pr;
    std::vector<double> temp_center_pip;
    std::vector<double> temp_center_kp;

    std::vector<double> temp_center_err_pr;
    std::vector<double> temp_center_err_pip;
    std::vector<double> temp_center_err_kp;

    std::vector<double> temp_fit_means_pr;    
    std::vector<double> temp_fit_means_pip;    
    std::vector<double> temp_fit_means_kp;    

    std::vector<double> temp_fit_sigmas_pr;    
    std::vector<double> temp_fit_sigmas_pip;    
    std::vector<double> temp_fit_sigmas_kp;    

    std::vector<double> temp_fit_means_err_pr;    
    std::vector<double> temp_fit_means_err_pip;    
    std::vector<double> temp_fit_means_err_kp;    

    std::vector<double> temp_fit_sigmas_err_pr;    
    std::vector<double> temp_fit_sigmas_err_pip;    
    std::vector<double> temp_fit_sigmas_err_kp;    

    std::cout << " PERFORMING FITS N TIMES: " << nfits_pr << std::endl;
    for( int j = 0; j < nfits_pr; j++ ){
      std::cout << " >> FITTING BETWEEN BINS " << start_bin_pr << " " << end_bin_pr << std::endl;

      TH1D *h_beta =  h2_betap_pr[i]->ProjectionY(Form("pr_betap_projY_s%d_%d",i,j),start_bin_pr, end_bin_pr );
      TF1 *fit_beta = fitHistogram(h_beta);

      v_can[i]->cd(j+1);
      h_beta->Draw("same");
      fit_beta->Draw("same");

      double min_p = h2_betap_pr[i]->GetXaxis()->GetBinCenter(start_bin_pr) - h2_betap_pr[i]->GetXaxis()->GetBinWidth(start_bin_pr)/2.0;
      double max_p = h2_betap_pr[i]->GetXaxis()->GetBinCenter(end_bin_pr) + h2_betap_pr[i]->GetXaxis()->GetBinWidth(end_bin_pr)/2.0;
      double slice_center = ( min_p + max_p )/2.0;

      std::cout << " MNTM RANGE OF SLICE " << min_p << " MAX " << max_p << " CENTER " << slice_center << std::endl;
            
      double mean = fit_beta->GetParameter(1);
      double mean_err = 0.0;//fit_beta->GetParError(1)/mean;

      double sig = fit_beta->GetParameter(2);
      double sig_error = fit_beta->GetParError(2);
      double sig_err = 0.0;// sig*( mean_err/mean + 3*sig/sig_error);
      
      temp_center_pr.push_back(slice_center);
      temp_center_err_pr.push_back(0.0);

      temp_fit_means_err_pr.push_back(mean_err);
      temp_fit_sigmas_err_pr.push_back(sig_err);
    
      temp_fit_means_pr.push_back(mean);
      temp_fit_sigmas_pr.push_back(sig);
     
      start_bin_pr+=ndiv_pr;
      end_bin_pr+=ndiv_pr;           
    
    }

    // %%%%%%%%%%% FITTING PIP BETA VS P NOW %%%%%%%%%%%%%%%%%
    std::cout << " PERFORMING FITS N TIMES: " << nfits_pip << std::endl;
    for( int j = 0; j < nfits_pip; j++ ){
      std::cout << " >> FITTING BETWEEN BINS " << start_bin_pip << " " << end_bin_pip << std::endl;


      TH1D *h_beta =  h2_betap_pip[i]->ProjectionY(Form("pip_betap_projY_s%d_%d",i,j),start_bin_pip, end_bin_pip );
      TF1 *fit_beta = fitHistogram(h_beta);

      v_can[i]->cd(j+1);
      h_beta->Draw("same");
      fit_beta->Draw("same");

      double min_p = h2_betap_pip[i]->GetXaxis()->GetBinCenter(start_bin_pip) - h2_betap_pip[i]->GetXaxis()->GetBinWidth(start_bin_pip)/2.0;
      double max_p = h2_betap_pip[i]->GetXaxis()->GetBinCenter(end_bin_pip) + h2_betap_pip[i]->GetXaxis()->GetBinWidth(end_bin_pip)/2.0;
      double slice_center = ( min_p + max_p )/2.0;

      std::cout << " MNTM RANGE OF SLICE " << min_p << " MAX " << max_p << " CENTER " << slice_center << std::endl;
            
      double mean = fit_beta->GetParameter(1);
      double mean_err = 0.0;//fit_beta->GetParError(1)/mean;

      double sig = fit_beta->GetParameter(2);
      double sig_error = fit_beta->GetParError(2);
      double sig_err = 0.0;// sig*( mean_err/mean + 3*sig/sig_error);
      
      temp_center_pip.push_back(slice_center);
      temp_center_err_pip.push_back(0.0);

      temp_fit_means_err_pip.push_back(mean_err);
      temp_fit_sigmas_err_pip.push_back(sig_err);
    
      temp_fit_means_pip.push_back(mean);
      temp_fit_sigmas_pip.push_back(sig);
     
      start_bin_pip+=ndiv_pip;
      end_bin_pip+=ndiv_pip;           
    
    }

    // %%%%%%%%%%%%%% FITTING KAON PLUS BETA VS P NOW %%%%%%%%%%%%%%%
    std::cout << " PERFORMING FITS N TIMES: " << nfits_kp << std::endl;
    for( int j = 0; j < nfits_kp; j++ ){
      std::cout << " >> FITTING BETWEEN BINS " << start_bin_kp << " " << end_bin_kp << std::endl;

      TH1D *h_beta =  h2_betap_kp[i]->ProjectionY(Form("kp_betap_projY_s%d_%d",i,j),start_bin_kp, end_bin_kp );
      TF1 *fit_beta = fitHistogram(h_beta);

      v_can[i]->cd(j+1);
      h_beta->Draw("same");
      fit_beta->Draw("same");

      double min_p = h2_betap_kp[i]->GetXaxis()->GetBinCenter(start_bin_kp) - h2_betap_kp[i]->GetXaxis()->GetBinWidth(start_bin_kp)/2.0;
      double max_p = h2_betap_kp[i]->GetXaxis()->GetBinCenter(end_bin_kp) + h2_betap_kp[i]->GetXaxis()->GetBinWidth(end_bin_kp)/2.0;
      double slice_center = ( min_p + max_p )/2.0;

      std::cout << " MNTM RANGE OF SLICE " << min_p << " MAX " << max_p << " CENTER " << slice_center << std::endl;
            
      double mean = fit_beta->GetParameter(1);
      double mean_err = 0.0;//fit_beta->GetParError(1)/mean;

      double sig = fit_beta->GetParameter(2);
      double sig_error = fit_beta->GetParError(2);
      double sig_err = 0.0;// sig*( mean_err/mean + 3*sig/sig_error);
      
      temp_center_kp.push_back(slice_center);
      temp_center_err_kp.push_back(0.0);

      temp_fit_means_err_kp.push_back(mean_err);
      temp_fit_sigmas_err_kp.push_back(sig_err);
    
      temp_fit_means_kp.push_back(mean);
      temp_fit_sigmas_kp.push_back(sig);
     
      start_bin_kp+=ndiv_kp;
      end_bin_kp+=ndiv_kp;           
    
    }

    h_beta_mean_pr.push_back( new TGraphErrors(temp_center_pr.size(), &(temp_center_pr[0]), &(temp_fit_means_pr[0]), &(temp_center_err_pr[0]), &(temp_fit_means_err_pr[0])) );
    h_beta_mean_pip.push_back( new TGraphErrors(temp_center_pip.size(), &(temp_center_pip[0]), &(temp_fit_means_pip[0]), &(temp_center_err_pip[0]), &(temp_fit_means_err_pip[0])) );
    h_beta_mean_kp.push_back( new TGraphErrors(temp_center_kp.size(), &(temp_center_kp[0]), &(temp_fit_means_kp[0]), &(temp_center_err_kp[0]), &(temp_fit_means_err_kp[0])) );

    h_beta_sig_pr.push_back( new TGraphErrors( temp_center_pr.size(), &(temp_center_pr[0]), &(temp_fit_sigmas_pr[0]), &(temp_center_err_pr[0]),&(temp_fit_sigmas_err_pr[0])) );
    h_beta_sig_pip.push_back( new TGraphErrors( temp_center_pip.size(), &(temp_center_pip[0]), &(temp_fit_sigmas_pip[0]), &(temp_center_err_pip[0]),&(temp_fit_sigmas_err_pip[0])) );
    h_beta_sig_kp.push_back( new TGraphErrors( temp_center_kp.size(), &(temp_center_kp[0]), &(temp_fit_sigmas_kp[0]), &(temp_center_err_kp[0]),&(temp_fit_sigmas_err_kp[0])) );

  }
    

  TCanvas *c_temp_pr = new TCanvas("c_final_pr", "c_final_pr", 800, 800);
  c_temp_pr->Divide(4,2);

  TCanvas *c_temp_pr_sig = new TCanvas("c_final_sig_pr", "c_final_sig_pr", 800, 800);
  c_temp_pr_sig->Divide(4,2);

  std::cout << " size " << h_beta_mean_pr.size() << std::endl;
  for(int s = 1; s <= 6; s++ ){
    c_temp_pr->cd(s+1);

    /// %%%%% MEAN FIT %%%%%%%
    TF1 *fit_mean_pr = new TF1(Form("fit_mean_pr_s%d",s),"[0] + [1]*x/sqrt(x*x + [2]) + [3]*x*x/sqrt(x*x*x*x + [4])", 0.1, 4.0);
    fit_mean_pr->SetParameter(0,0.0);
    fit_mean_pr->SetParameter(1,1.0);
    fit_mean_pr->SetParameter(2,prot_mass);
    fit_mean_pr->SetParameter(3,1.0);
    fit_mean_pr->SetParameter(4,1.0);

    h_beta_mean_pr[s]->Fit(Form("fit_mean_pr_s%d",s),"R");


    h_beta_mean_pr[s]->SetTitle(Form("PR BETA MEAN FOR SECTOR %d", s));
    h_beta_mean_pr[s]->SetMarkerColor(kBlack);
    h_beta_mean_pr[s]->SetMarkerStyle(8);
    h_beta_mean_pr[s]->Draw("AP");       

    fit_mean_pr->SetLineColor(kRed);
    fit_mean_pr->Draw("same");

    c_temp_pr_sig->cd(s+1);

    /// %%%%% SIGMA FIT %%%%%%%
    TF1 *fit_sigma_pr = new TF1(Form("fit_sigma_pr_s%d",s),"[0] + [1]/x + [2]*exp(-x)", 0.1, 4.0);
    fit_sigma_pr->SetParameter(0,0.0);
    fit_sigma_pr->SetParameter(1,1.0);
    fit_sigma_pr->SetParameter(2,prot_mass);
    fit_sigma_pr->SetParameter(3,1.0);
    fit_sigma_pr->SetParameter(4,1.0);

    h_beta_sig_pr[s]->Fit(Form("fit_mean_pr_s%d",s),"R");

    h_beta_sig_pr[s]->SetTitle(Form("PR BETA MEAN FOR SECTOR %d", s));
    h_beta_sig_pr[s]->SetMarkerColor(kBlack);
    h_beta_sig_pr[s]->SetMarkerStyle(8);
    h_beta_sig_pr[s]->Draw("AP");       

    fit_sigma_pr->SetLineColor(kRed);
    fit_sigma_pr->Draw("same");
  }

  TCanvas *c_temp_pip = new TCanvas("c_final_pip", "c_final_pip", 800, 800);
  c_temp_pip->Divide(4,2);

  TCanvas *c_temp_pip_sig = new TCanvas("c_final_sig_pip", "c_final_sig_pip", 800, 800);
  c_temp_pip_sig->Divide(4,2);

  std::cout << " size " << h_beta_mean_pip.size() << std::endl;
  for(int s = 1; s <= 6; s++ ){
    c_temp_pip->cd(s+1);

    /// %%%%% MEAN FIT %%%%%%%
    TF1 *fit_mean_pip = new TF1(Form("fit_mean_pip_s%d",s),"[0] + [1]*x/sqrt(x*x + [2]) + [3]*x*x/sqrt(x*x*x*x + [4])", 0.1, 4.0);
    fit_mean_pip->SetParameter(0,0.0);
    fit_mean_pip->SetParameter(1,1.0);
    fit_mean_pip->SetParameter(2,pip_mass);
    fit_mean_pip->SetParameter(3,1.0);
    fit_mean_pip->SetParameter(4,1.0);   

    h_beta_mean_pip[s]->Fit(Form("fit_mean_pip_s%d",s),"R");

    h_beta_mean_pip[s]->SetTitle(Form("PIP BETA MEAN FOR SECTOR %d", s));
    h_beta_mean_pip[s]->SetMarkerColor(kBlack);
    h_beta_mean_pip[s]->SetMarkerStyle(8);
    h_beta_mean_pip[s]->Draw("AP");       

    fit_mean_pip->SetLineColor(kRed);
    fit_mean_pip->Draw("same");

    c_temp_pip_sig->cd(s+1);
    /// %%%%% SIGMA FIT %%%%%%%
    TF1 *fit_sigma_pip = new TF1(Form("fit_sigma_pip_s%d",s),"[0] + [1]/x + [2]*exp(-x)", 0.1, 4.0);
    fit_sigma_pip->SetParameter(0,0.0);
    fit_sigma_pip->SetParameter(1,1.0);
    fit_sigma_pip->SetParameter(2,pip_mass);
    fit_sigma_pip->SetParameter(3,1.0);
    fit_sigma_pip->SetParameter(4,1.0);

    h_beta_sig_pip[s]->Fit(Form("fit_mean_pip_s%d",s),"R");

    h_beta_sig_pip[s]->SetTitle(Form("PIP BETA MEAN FOR SECTOR %d", s));
    h_beta_sig_pip[s]->SetMarkerColor(kBlack);
    h_beta_sig_pip[s]->SetMarkerStyle(8);
    h_beta_sig_pip[s]->Draw("AP");       

    fit_sigma_pip->SetLineColor(kRed);
    fit_sigma_pip->SetLineStyle(3);
    fit_sigma_pip->Draw("same");
  }

  TCanvas *c_temp_kp = new TCanvas("c_final_kp", "c_final_kp", 800, 800);
  c_temp_kp->Divide(4,2);

  TCanvas *c_temp_kp_sig = new TCanvas("c_final_sig_kp", "c_final_sig_kp", 800, 800);
  c_temp_kp_sig->Divide(4,2);

  std::cout << " size " << h_beta_mean_kp.size() << std::endl;
  for(int s = 1; s <= 6; s++ ){
    c_temp_kp->cd(s+1);

    /// %%%%% MEAN FIT %%%%%%%%%
    TF1 *fit_mean_kp = new TF1(Form("fit_mean_kp_s%d",s),"[0] + [1]*x/sqrt(x*x + [2]) + [3]*x*x/sqrt(x*x*x*x + [4])", 0.0, 4.0);
    fit_mean_kp->SetParameter(0,0.50);
    fit_mean_kp->SetParameter(1,0.50);
    fit_mean_kp->SetParameter(2,kp_mass);
    fit_mean_kp->SetParameter(3,0.50);
    fit_mean_kp->SetParameter(4,0.50);   

    h_beta_mean_kp[s]->Fit(Form("fit_mean_kp_s%d",s),"R");

    h_beta_mean_kp[s]->SetTitle(Form("KP BETA MEAN FOR SECTOR %d", s));
    h_beta_mean_kp[s]->SetMarkerColor(kBlack);
    h_beta_mean_kp[s]->SetMarkerStyle(8);
    h_beta_mean_kp[s]->Draw("AP");      

    fit_mean_kp->SetLineColor(kRed);
    fit_mean_kp->Draw("same"); 
    
    /// %%%%% SIGMA FIT %%%%%%%
    TF1 *fit_sigma_kp = new TF1(Form("fit_sigma_kp_s%d",s),"[0] + [1]/x + [2]*exp(-x)", 0.1, 4.0);
    fit_sigma_kp->SetParameter(0,0.0);
    fit_sigma_kp->SetParameter(1,1.0);
    fit_sigma_kp->SetParameter(2,kp_mass);
    fit_sigma_kp->SetParameter(3,1.0);
    fit_sigma_kp->SetParameter(4,1.0);
    
    h_beta_sig_kp[s]->Fit(Form("fit_mean_kp_s%d",s),"R");  

    h_beta_sig_kp[s]->SetTitle(Form("KP BETA MEAN FOR SECTOR %d", s));
    h_beta_sig_kp[s]->SetMarkerColor(kBlack);
    h_beta_sig_kp[s]->SetMarkerStyle(8);
    h_beta_sig_kp[s]->Draw("AP");      
   
    fit_sigma_kp->SetLineColor(kRed);
    fit_sigma_kp->SetLineStyle(3);
    fit_sigma_kp->Draw("same"); 
  }
  
  return 0;
}


//std::map<int, std::vector< std::vector< double > > > > SliceNFitHistograms
