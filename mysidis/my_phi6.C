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
#include "programFiles/hadronMLEID.C"
#include "programFiles/getGenIndices.C"
#include "MomCorr.C"
#include <vector>
#include <map>
#include "loadBetaParameters.C"
#include "getBeta.C"


int my_phi6(int iteration_number = 0,  int filestart = 1, int fileend = 3, int ExpOrSim = 1, bool do_momCorr_e = 0, bool do_momCorr_pions = 1, int cut_to_vary = 0, int cut_strictness = 0 ){ //int e_zvertex_strict = 0, int e_ECsampling_strict = 0, int e_ECoutVin_strict = 0, int e_ECgeometric_strict = 0, int e_CCthetaMatching_strict = 0, int e_R1fid_strict = 0, int e_R3fid_strict = 0, int e_CCphiMatching_strict = 0, int e_CCfiducial_strict = 0, int yCut_strict = 0, int pip_vvp_strict = 0, int pip_R1fid_strict = 0, int pip_MXcut_strict = 0, int pim_vvp_strict = 0, int pim_R1fid_strict = 0, int pim_MXcut_strict = 0)
  std::cout << " >> BEGINNING PHI ANALYSIS " << std::endl;
  TStopwatch *stopwatch = new TStopwatch();
  
  TChain *h22chain = new TChain("h22");  
  string firstfilename = "";
  string lastfilename = "";

  // %%%%%%%% read the input files into the TChain %%%%%%%%
  int NtotalFiles = 590;// 11625 must be before a merging of files
  ifstream filelist;
  filelist.open("programFiles/dataFiles.txt");
  //int kStop = Nfiles + filestart;
  int kStop = fileend;
  //if(kStop > NtotalFiles+1) kStop = NtotalFiles+1;
  int total_files_to_process = fileend - filestart;
  std::string line;
  int line_counter = 0;

  std::cout << " >> ADDING FILES TOGETHER INTO ONE CHAIN " << std::endl;
  for(int k = 0; k <= NtotalFiles; k++){
    string filename;
    filelist>>filename;
    if( k >= filestart && k <= fileend ){
      std::cout << " >> ADDING FILE " << filename.c_str() << std::endl;
      if(k >= filestart) h22chain->Add(filename.c_str());
      if(k == filestart) firstfilename = filename;
      if(k == kStop - 1) lastfilename = filename;
    }
  }
  
  firstfilename = firstfilename.substr(39,10); // trim the string down so it's just the run#.subrun# (e.g. 38458.a00)
  lastfilename = lastfilename.substr(39,10);
  cout<<"first file: "<<firstfilename<<" ... last file: "<<lastfilename<<endl;
  // %%%%% end read the input files into a TChain %%%%%

  MomCorr_e1f *MomCorr = new MomCorr_e1f();

  std::cout << " DONE WITH MOMCORR CLASS " << std::endl;

  Int_t gpart;
  Int_t mcentr;
  Int_t mcid;
  Float_t mctheta;
  Float_t mcphi;
  Float_t mcp;
    
  Float_t p[35];
  Int_t q[35];
  Float_t cx[35];
  Float_t cy[35];
  Float_t cz[35];
  Float_t vz[35];
  Float_t vy[35];
  Float_t vx[35];
  Float_t b[35];
    
  UChar_t dc_sect[35];
  Float_t tl1_x[35];
  Float_t tl1_y[35];
  Float_t tl3_x[35]; // used to be dc_xsc in h10
  Float_t tl3_y[35];
  Float_t tl3_z[35];
  Float_t tl3_cx[35]; // used to be dc_cxsc in h10
  Float_t tl3_cy[35];
  Float_t tl3_cz[35];
    
  UChar_t ec_sect[35];
  Float_t etot[35];
  Float_t ec_ei[35];
  Float_t ec_eo[35];
  Float_t ec_t[35];
  Float_t ech_x[35];
  Float_t ech_y[35];
  Float_t ech_z[35];
    
  UChar_t sc_sect[35];
  Float_t sc_t[35];
  Float_t sc_r[35];
  UChar_t sc_pd[35];
    
  UChar_t cc_sect[35];
  UShort_t cc_segm[35];
  UShort_t nphe[35];
  // %%%%% end define variables %%%%%
    
  // %%%%% set branch addresses %%%%%    
  h22chain->SetBranchAddress("gpart", &gpart);
  h22chain->SetBranchAddress("p", p);
  h22chain->SetBranchAddress("q", q);
  h22chain->SetBranchAddress("cx", cx);
  h22chain->SetBranchAddress("cy", cy);
  h22chain->SetBranchAddress("cz", cz);
  h22chain->SetBranchAddress("vz", vz);
  h22chain->SetBranchAddress("vy", vy);
  h22chain->SetBranchAddress("vx", vx);
  h22chain->SetBranchAddress("b", b);
    
  h22chain->SetBranchAddress("dc_sect", dc_sect);
  h22chain->SetBranchAddress("tl1_x", tl1_x);
  h22chain->SetBranchAddress("tl1_y", tl1_y);
  h22chain->SetBranchAddress("tl3_x", tl3_x);
  h22chain->SetBranchAddress("tl3_y", tl3_y);
  h22chain->SetBranchAddress("tl3_z", tl3_z);
  h22chain->SetBranchAddress("tl3_cx", tl3_cx);
  h22chain->SetBranchAddress("tl3_cy", tl3_cy);
  h22chain->SetBranchAddress("tl3_cz", tl3_cz);
    
  h22chain->SetBranchAddress("ec_sect", ec_sect);
  h22chain->SetBranchAddress("etot", etot);
  h22chain->SetBranchAddress("ec_ei", ec_ei);
  h22chain->SetBranchAddress("ec_eo", ec_eo);
  h22chain->SetBranchAddress("ec_t", ec_t);
  h22chain->SetBranchAddress("ech_x", ech_x);
  h22chain->SetBranchAddress("ech_y", ech_y);
  h22chain->SetBranchAddress("ech_z", ech_z);
    
  h22chain->SetBranchAddress("sc_sect", sc_sect);
  h22chain->SetBranchAddress("sc_t", sc_t);
  h22chain->SetBranchAddress("sc_r", sc_r);
  h22chain->SetBranchAddress("sc_pd", sc_pd);
    
  h22chain->SetBranchAddress("cc_sect", cc_sect);
  h22chain->SetBranchAddress("cc_segm", cc_segm);
  h22chain->SetBranchAddress("nphe", nphe);
  // %%%%% end set branch addresses %%%%%
    
  // %%%%% stuff for output root files %%%%%
  std::string cut_number = std::to_string(cut_to_vary);
  std::string cut_strict_lvl;
  std::string out_dir;
  if( cut_strictness == -1 ){
    cut_strict_lvl = "l";
    out_dir = "/home/bclary/clas6/e1f/retro-sidis/mysidis/phi_loose/";
  }
  else if( cut_strictness == 0 ){
    cut_strict_lvl = "n";
    out_dir = "/home/bclary/clas6/e1f/retro-sidis/mysidis/phi_nominal/";
  }
  else if( cut_strictness == 1 ){
    cut_strict_lvl = "t";
    out_dir = "/home/bclary/clas6/e1f/retro-sidis/mysidis/phi_tight/";
  }

  
  string outfilename = out_dir +  "my_phi6_"+ std::to_string(filestart) + "_" + std::to_string(fileend) + "_clvl_" + cut_number + "_cs_" + cut_strict_lvl + ".root";
  /*TFile *outputfile;
  outputfile = new TFile(outfilename.c_str(),"recreate");
  std::cout << " >> SAVING TO FILE " << outfilename << std::endl;
  TTree *output_tree = new TTree("output_tree","Final Particles");
  TLorentzVector *final_el = new TLorentzVector(0,0,0,0);
  TLorentzVector *final_pr = new TLorentzVector(0,0,0,0);
  TLorentzVector *final_kp = new TLorentzVector(0,0,0,0);
  TLorentzVector *final_km = new TLorentzVector(0,0,0,0);
  output_tree->Branch("final_el","TLorentzVector",&final_el);
  output_tree->Branch("final_pr","TLorentzVector",&final_pr);
  output_tree->Branch("final_kp","TLorentzVector",&final_kp);
  //output_tree->Branch("final_el","TLorentzVector",&final_el);
  */
  std::cout << " GETTING NUMNER OF ENTRIES " << std::endl;
  cout<<"entries: "<<h22chain->GetEntries()/1000<<" thousand"<<endl;

  Float_t e_mass = 0.000511; // GeV
  Float_t prot_mass = 0.938272; // GeV
  Float_t pip_mass = 0.13957; // GeV
  Float_t pim_mass = 0.13957; // GeV
  Float_t kp_mass = 0.4937; // GeV
  Float_t km_mass = 0.4937; // GeV
  
  int proton_id = 2212;
  int kaon_plus_id = 321;
  int pion_plus_id = 211;

  Float_t speed_of_light = 29.9792458; // cm/ns
  Float_t Beam_Energy = 5.498; // GeV
  Float_t pi = 3.14159265359;
  Float_t pi180 = pi/180.0;
  Float_t pi180_inv = 180.0/pi;

  TLorentzVector V4k(0.0, 0.0, Beam_Energy, Beam_Energy); // x, y, z, t
  TLorentzVector V4ISproton(0.0, 0.0, 0.0, prot_mass); // IS = Initial State

  // %%%%%%%%%% cut settings %%%%%%%%%%%
  int e_zvertex_strict = 0;
  int e_ECsampling_strict = 0;
  int e_ECoutVin_strict = 0;
  int e_CCthetaMatching_strict = 0;
  int e_ECgeometric_strict = 0;
  int e_R1fid_strict = 0;
  int e_R3fid_strict = 0;
  int e_CCphiMatching_strict = 0;
  int e_CCfiducial_strict = 0;

  int pip_vvp_strict = 0;
  int pip_R1fid_strict = 0;
  int pip_MXcut_strict = 0; 
  int pim_vvp_strict = 0; 
  int pim_R1fid_strict = 0;
  int pim_MXcut_strict = 0;


  if( cut_to_vary < 0 ){
    std::cout <<" NOMINAL - NO CUTS VARIED " << std::endl;;
  }
  else if( cut_to_vary == 0 ){
    e_zvertex_strict = cut_strictness;
  }
  else if( cut_to_vary == 1 ){
    e_ECsampling_strict = cut_strictness;
  }
  else if( cut_to_vary == 2 ){
    e_ECoutVin_strict = cut_strictness;
  }
  else if( cut_to_vary == 3 ){
    e_CCthetaMatching_strict = cut_strictness;
  }
  else if( cut_to_vary == 4 ){
    e_ECgeometric_strict = cut_strictness;
  }
  else if( cut_to_vary == 5 ){
    e_R1fid_strict = cut_strictness;
  }
  else if( cut_to_vary == 6 ){
    e_R3fid_strict = cut_strictness;
  }
  else if( cut_to_vary == 7 ){
    e_CCphiMatching_strict = cut_strictness;
  }
  else if( cut_to_vary == 8 ){
    e_CCfiducial_strict = cut_strictness;
  }
  else if( cut_to_vary == 9 ){ //changes the pip to pr and the pim to kp as these are the only particles of interest for now.include km later?
    pip_vvp_strict = cut_strictness;
  }
  else if( cut_to_vary == 10 ){
    pip_R1fid_strict = cut_strictness;
  }
  else if( cut_to_vary == 11 ){
    pip_MXcut_strict = cut_strictness;
  }
  else if( cut_to_vary == 12 ){
    pim_vvp_strict = cut_strictness;
  }
  else if( cut_to_vary == 13 ){
    pim_R1fid_strict = cut_strictness;
  }
  else if( cut_to_vary == 14 ){
    pim_MXcut_strict = cut_strictness;
  }

  // %%%%%%%%%% end cut settings %%%%%%%%%

  //int do_momCorr_e = 0;
  int do_momCorr_proton = 0;
  int do_momCorr_kaon = 0;

  double pr_conf = 0.01;// DR
  double pip_conf = 0.01; 
  double kp_conf = 0.01;
  double pr_anticonf = 0.95;
  double pip_anticonf = 0.95;
  double kp_anticonf = 0.95;
  std::cout << " >> CREATING HISTOGRAMS " << std::endl;

  //%%%%%%%%%%%%%%% cut containers for beta MLE parameters %%%%%%%%%%%%%%
  //include run dependence here
  std::map<int, std::vector<double> > pr_mean_fit = loadBetaParameters("protonBetaFitParameters.txt");
  std::map< int , std::vector<double> >::iterator it;
  for( it = pr_mean_fit.begin(); it != pr_mean_fit.end(); it++ ){
    //std::cout << " FIT PARAMETER " << it->first << std::endl;
    //std::cout << "SIZE PF VECTOR " << (it->second).size()  << std::endl;
    for( std::vector<double>::iterator it2 = (it->second).begin(); it2 != (it->second).end(); ++it2){
      //std::cout << " SECTOR VALUES " << *it2 << std::endl;
    }

  }
  std::map<int, std::vector<double> > pip_mean_fit = loadBetaParameters("pipMeanBetaFitParameters.txt");
  std::map<int, std::vector<double> > kp_mean_fit = loadBetaParameters("kpMeanBetaFitParameters.txt");

  std::map<int, std::vector<double> > pr_sigma_fit = loadBetaParameters("protonSigmaBetaFitParameters.txt");
  std::map<int, std::vector<double> > pip_sigma_fit = loadBetaParameters("pipSigmaFitParameters.txt");
  std::map<int, std::vector<double> > kp_sigma_fit = loadBetaParameters("kpSigmaBetaFitParameters.txt");;

  //%%%%%%%%%%%%%%% cut container end %%%%%%%%%%%%%%

  //%%%%%%%%%%%%%%% histograms for electrons, protons, and kaons %%%%%%%%%%%%%%%%%%

  TH1D *h_el_p = new TH1D("h_el_p","h_el_p", 200, 0.0, Beam_Energy);
  TH1D *h_pr_p = new TH1D("h_pr_p","h_pr_p", 200, 0.0, Beam_Energy);
  TH1D *h_kp_p = new TH1D("h_kp_p","h_kp_p", 200, 0.0, Beam_Energy);

  TH1D *h_el_p_pc = new TH1D("h_el_p_pc","h_el_p_pc", 200, 0.0, Beam_Energy);
  TH1D *h_pr_p_pc = new TH1D("h_pr_p_pc","h_pr_p_pc", 200, 0.0, Beam_Energy);
  TH1D *h_kp_p_pc = new TH1D("h_kp_p_pc","h_kp_p_pc", 200, 0.0, Beam_Energy);

  TH2D *h_el_thetap_pc = new TH2D("h_el_thetap_pc","h_el_thetap_pc", 200, 0.0, Beam_Energy, 200, 0.0, 65.0);
  TH2D *h_pr_thetap_pc = new TH2D("h_pr_thetap_pc","h_pr_thetap_pc", 200, 0.0, Beam_Energy, 200, 0.0, 70.0);
  TH2D *h_kp_thetap_pc = new TH2D("h_kp_thetap_pc","h_kp_thetap_pc", 200, 0.0, Beam_Energy, 200, 0.0, 70.0);
  
  TH2D *h_el_thetaphi = new TH2D("h_el_thetaphi","h_el_thetaphi", 200, 0.0, 60.0, 200, -30.0, 30.0);

  TH2D *h_el_vzphi = new TH2D("h_el_vzphi","h_el_vzphi",200, -35.0, -15.0, 200, -180.0, 180.0 );
  TH2D *h_el_vzphi_vc = new TH2D("h_el_vzphi_vc","h_el_vzphi_vc",200, -35.0, -15.0, 200, -180.0, 180.0 );
  TH2D *h_el_vzphi_pc = new TH2D("h_el_vzphi_pc","h_el_vzphi_pc",200, -35.0, -15.0, 200, -180.0, 180.0 );
  TH2D *h_el_vzphi_vcpc = new TH2D("h_el_vzphi_vcpc","h_el_vzphi_vcpc", 200, -35.0, -15.0, 200, -180.0, 180.0 );

  TH2D *h_el_etotp = new TH2D("h_el_etotp","h_el_etotp", 200, 0.0, Beam_Energy, 200, 0.0, 0.50);
  TH2D *h_el_etotp_pc = new TH2D("h_el_etotp_pc","h_el_etotp_pc", 200, 0.0, Beam_Energy, 200, 0.0, 0.50);
  TH1D *h_el_nphe = new TH1D("h_el_nphe","h_el_nphe", 200, 0.0, 200);

  TH1D *h_el_q2 = new TH1D("h_el_q2","h_el_q2", 200, 0.0, Beam_Energy);
  TH1D *h_el_q2_pc = new TH1D("h_el_q2_pc","h_el_q2_pc",200, 0.0, Beam_Energy);
  
  //delta time assuming kaon mass (kaon should be centered about 0)
  TH2D *h_deltime_p_kp = new TH2D("h_deltime_p_kp","h_deltime_p_kp", 200, 0.0, 3.0, 200, -10.0, 10.0 );

  TH2D *h_pr_betap = new TH2D("h_pr_betap","h_pr_betap",200, 0.0, 3.0, 200, 0.0, 1.1 );
  TH2D *h_kp_betap = new TH2D("h_kp_betap","h_kp_betap",200, 0.0, 3.0, 200, 0.0, 1.1 );

  TH2D *h_pr_thetaphi = new TH2D("h_pr_thetaphi","h_pr_thetaphi", 200, 0.0, 70.0, 200, -45.0, 45.0);
  TH2D *h_kp_thetaphi = new TH2D("h_kp_thetaphi","h_kp_thetaphi", 200, 0.0, 70.0, 200, -45.0, 45.0);
  
  //fiducials plots
  TH2D *h_el_ec = new TH2D("h_el_ec","h_el_ec", 1000, -500, 500, 1000, -500, 500);
  TH2D *h_el_dcr1 = new TH2D("h_el_dcr1","h_el_dcr1", 500, -100.0, 100.0, 500, -100, 100 );
  TH2D *h_el_dcr3 = new TH2D("h_el_dcr3","h_el_dcr3", 500, 50.0, 400.0, 500, -200.0, 200.0);

  //all negative plots
  TH2D *h_el_etotp_b = new TH2D("h_el_etotp_b","h_el_etotp_b", 200, 0.0, Beam_Energy, 200, 0.0, 0.50);

  TH2D *h_el_ec_b = new TH2D("h_el_ec_b","h_el_ec_b", 1000, -500, 500, 1000, -500, 500);
  TH2D *h_el_dcr1_b = new TH2D("h_el_dcr1_b","h_el_dcr1_b", 500, -100.0, 100.0, 500, -100, 100 );
  TH2D *h_el_dcr3_b = new TH2D("h_el_dcr3_b","h_el_dcr3_b", 500, 50.0, 400.0, 500, -200.0, 200.0);

  //all positive betas
  TH2D *h_pos_betap = new TH2D("h_pos_betap","h_pos_betap", 200, 0.0, 4.5, 200, 0.0, 1.1);

  //all negative betas
  TH2D *h_neg_betap = new TH2D("h_neg_betap","h_neg_betap", 200, 0.0, 4.5, 200, 0.0, 1.1);

  // mle plots
  TH1D *h_prot_conf = new TH1D("h_prot_conf","h_prot_conf", 200, 0.0, 1.05);
  TH1D *h_kp_conf = new TH1D("h_kp_conf","h_kp_conf", 200, 0.0, 1.05);
  TH1D *h_pip_conf = new TH1D("h_pip_conf","h_pip_conf", 200, 0.0, 1.05);
    
  std::vector<TH1D*> h_el_sect_p;
  std::vector<TH1D*> h_pr_sect_p;
  std::vector<TH1D*> h_kp_sect_p;
  std::vector<TH1D*> h_el_sect_nphe;
  std::vector<TH1D*> h_el_sect_q2;

  std::vector<TH1D*> h_el_sect_p_pc;
  std::vector<TH1D*> h_pr_sect_p_pc;
  std::vector<TH1D*> h_kp_sect_p_pc;
  std::vector<TH1D*> h_el_sect_q2_pc;

  std::vector<TH2D*> h_el_sect_etotp;
  std::vector<TH2D*> h_el_sect_etotp_pc;

  std::vector<TH1D*> h_sect_deltime_p_kp;
  std::vector<TH2D*> h_pr_sect_betap;
  std::vector<TH2D*> h_kp_sect_betap;

  std::vector<TH2D*> h_pr_sect_betap_pc;
  std::vector<TH2D*> h_kp_sect_betap_pc;

  std::vector<TH2D*> h_el_sect_thetaphi;
  std::vector<TH2D*> h_pr_sect_thetaphi;
  std::vector<TH2D*> h_kp_sect_thetaphi;

  std::vector<TH2D*> h_el_sect_cctheta;
  std::vector<TH2D*> h_el_sect_ccphi;

  std::vector<TH2D*> h_el_sect_dcr3;
  
  // all negative charges plotted 
  std::vector<TH2D*> h_el_sect_etotp_b;
  std::vector<TH2D*> h_el_sect_cctheta_b;
  std::vector<TH2D*> h_neg_sect_betap;

  // all positive charges plots
  std::vector<TH2D*> h_pos_sect_betap;

  

  for( int i = 0; i <= 6; i++ ){

    h_el_sect_p.push_back( new TH1D(Form("h_el_s%d_p",i),Form("h_el_s%d_p",i), 200, 0.0, Beam_Energy) );
    h_el_sect_nphe.push_back( new TH1D(Form("h_el_s%d_nphe",i),Form("h_el_s%d_nphe",i), 200, 0.0, 200) );
    h_el_sect_q2.push_back( new TH1D(Form("h_el_s%d_q2",i),Form("h_el_s%d_q2",i), 200, 0.0, Beam_Energy) );
    h_el_sect_p_pc.push_back( new TH1D(Form("h_el_s%d_p_pc",i),Form("h_el_s%d_p_pc",i), 200, 0.0, Beam_Energy) );
    h_el_sect_q2_pc.push_back( new TH1D(Form("h_el_s%d_q2_pc",i),Form("h_el_s%d_q2_pc",i), 200, 0.0, Beam_Energy) );
    h_el_sect_etotp.push_back( new TH2D(Form("h_el_s%d_etotp",i),Form("h_el_s%d_etotp",i), 200, 0.0, Beam_Energy, 200, 0.0, 0.50) );
    h_el_sect_etotp_pc.push_back( new TH2D(Form("h_el_s%d_etotp_pc",i),Form("h_el_s%d_etotp_pc",i), 200, 0.0, Beam_Energy, 200, 0.0, 0.50) );
    h_el_sect_thetaphi.push_back( new TH2D(Form("h_el_s%d_thetaphi",i),Form("h_el_s%d_thetaphi",i), 200, 0.0, 60.0, 200, -30.0, 30.0) );
    h_el_sect_cctheta.push_back( new TH2D(Form("h_el_s%d_cctheta",i), Form("h_el_s%d_cctheta",i), 19, 0.0, 19.0, 100, 0.0, 60.0 ) );
    h_el_sect_ccphi.push_back( new TH2D(Form("h_el_s%d_ccphi",i), Form("h_el_s%d_ccphi",i), 19, 0.0, 19.0, 100, 0.0, 360.0 ) );    
    h_el_sect_dcr3.push_back( new TH2D(Form("h_el_s%d_dcr3",i), Form("h_el_s%d_dcr3",i), 500, 50, 400, 500, -200.0, 200) );
    
    h_pr_sect_p.push_back( new TH1D(Form("h_pr_s%d_p",i),Form("h_pr_s%d_p",i), 200, 0.0, 4.0) );
    h_pr_sect_p_pc.push_back( new TH1D(Form("h_pr_s%d_p_pc",i),Form("h_pr_s%d_p_pc",i), 200, 0.0, 4.0) );
    h_pr_sect_betap.push_back( new TH2D(Form("h_pr_s%d_betap",i),Form("h_pr_s%d_betap",i), 200, 0.0, 4.0, 200, 0.0, 1.1) );
    h_pr_sect_betap_pc.push_back( new TH2D(Form("h_pr_s%d_betap_pc",i),Form("h_pr_s%d_betap_pc",i), 200, 0.0, 4.0, 200, 0.0, 1.1) );
    h_pr_sect_thetaphi.push_back( new TH2D(Form("h_pr_s%d_thetaphi",i), Form("h_pr_s%d_thetaphi",i), 200, 0.0, 65.0, 200, -30.0, 30.0));

    h_kp_sect_p.push_back( new TH1D(Form("h_kp_s%d_p",i), Form("h_kp_s%d_p",i), 200, 0.0, 4.0) );
    h_kp_sect_p_pc.push_back( new TH1D(Form("h_kp_s%d_p_pc",i), Form("h_kp_s%d_p_pc",i), 200, 0.0, 4.0) );
    h_kp_sect_betap.push_back( new TH2D(Form("h_kp_s%d_betap",i), Form("h_kp_s%d_betap",i), 200, 0.0, 4.0, 200, 0.0, 1.1) );
    h_kp_sect_betap_pc.push_back( new TH2D(Form("h_kp_s%d_betap_pc",i), Form("h_kp_s%d_betap_pc",i), 200, 0.0, 4.0, 200, 0.0, 1.1) );
    h_kp_sect_thetaphi.push_back( new TH2D(Form("h_kp_s%d_thetaphi",i), Form("h_kp_s%d_thetaphi",i), 200, 0.0, 65.0, 200, -30.0, 30.0) );

    h_el_sect_etotp_b.push_back( new TH2D(Form("h_el_s%d_etotp_b",i),Form("h_el_s%d_etotp_b",i), 200, 0.0, Beam_Energy, 200, 0.0, 0.50) );
    h_el_sect_cctheta_b.push_back( new TH2D(Form("h_el_s%d_cctheta_b",i), Form("h_el_s%d_cctheta_b",i), 19, 0.0, 19.0, 100, 0.0, 60.0 ) );

    h_pos_sect_betap.push_back( new TH2D(Form("h_pos_s%d_betap",i), Form("h_pos_s%d_betap",i), 200, 0.0, 4.5, 200, 0.0, 1.1) );
    h_neg_sect_betap.push_back( new TH2D(Form("h_neg_s%d_betap",i), Form("h_neg_s%d_betap",i), 200, 0.0, 4.5, 200, 0.0, 1.1) );
   
  }

  int limit = h22chain->GetEntries() / 100;
  std::cout << " LIMIT " << limit << std::endl;
  
  std::cout << " TOTAL ENTRIES TO PROCESS " << h22chain->GetEntries() << std::endl;
  for(int i = 0; i < h22chain->GetEntries(); i++){
    if( i == 0 ) std::cout << " >> EXECUTING MAIN LOOP ENTRY 0 " << std::endl;
    //std::cout << i << std::endl;
    if ( i % limit == 0 ){
      double completed = (double)i / (double)h22chain->GetEntries();
      std::cout << " >> COMPLETED " << ceil(completed * 100) << " % --- ON EVENT " << i << std::endl;
    }
    
    h22chain->GetEntry(i);
    
    string currentrunno_string = "-123";
    if(ExpOrSim == 1)
      {
	currentrunno_string = h22chain->GetCurrentFile()->GetName();
	currentrunno_string = currentrunno_string.substr(36,5);
      }
    int currentrunno = atoi(currentrunno_string.c_str()); // converts the string to an int
    
    //counter number of good candidate particles
    int n_good_el = 0;
    int n_good_pr = 0;
    int n_good_kp = 0;

    TVector3 V3_e[2]; // 2 for gen(0) and rec(1)
    TLorentzVector V4_e[2];
    TLorentzVector V4_e_uncorr[2];
        
    // %%%%% electron ID %%%%%
    int e_index[2] = {-123,-123};
    //gets the gpart index for the electron
    e_index[1] = eID(gpart, q, p, cc_sect, sc_sect, ec_sect, dc_sect, cx, cy, cz, tl1_x, tl1_y, tl3_x, tl3_y, tl3_z, tl3_cx, tl3_cy, tl3_cz, e_zvertex_strict, vz, vy, vx, e_ECsampling_strict, ExpOrSim, etot, e_ECoutVin_strict, ec_ei, ech_x, ech_y, ech_z, e_CCthetaMatching_strict, cc_segm, e_ECgeometric_strict, e_R1fid_strict, e_R3fid_strict, e_CCphiMatching_strict, sc_pd, e_CCfiducial_strict);
    //%%%%% end electron ID %%%%%    
        
    if(e_index[1] > -122){
      V3_e[1].SetXYZ(p[e_index[1]]*cx[e_index[1]], p[e_index[1]]*cy[e_index[1]], p[e_index[1]]*cz[e_index[1]]);
      V4_e[1].SetXYZT(V3_e[1].X(), V3_e[1].Y(), V3_e[1].Z(), sqrt(V3_e[1].Mag2() + pow(e_mass,2)));
      V4_e_uncorr[1].SetXYZT(V3_e[1].X(), V3_e[1].Y(), V3_e[1].Z(), sqrt(V3_e[1].Mag2() + pow(e_mass,2)));
      if(do_momCorr_e && ExpOrSim == 1) V4_e[1] = MomCorr->PcorN(V4_e[1], -1, 11);
      n_good_el++;
      h_el_thetap_pc->Fill(V4_e[1].P(), 180.0/TMath::Pi() * V4_e[1].Theta() );
      //final_el->SetXYZT(V4_e[1].X(), V4_e[1].Y(), V4_e[1].Z(), V4_e[1].T() );
    }
    
    TVector3 V3_H[2];
    TLorentzVector V4_q[2], V4_H[2];    

    for( int i = 0; i < gpart; i++ ){
      if( q[i] < 0 ){

	//plot fiducial regions  and sampling fraction for all negatives 
	double ectot = ec_ei[e_index[1]] + ec_eo[e_index[1]];
	h_el_etotp_b->Fill(p[e_index[1]],ectot/p[e_index[1]]);

      	h_el_ec_b->Fill(ech_x[i], ech_y[i]);
	h_el_dcr1_b->Fill(tl1_x[i], tl1_y[i]);
	h_el_dcr3_b->Fill(tl3_x[i], tl3_y[i]);
	
	int temp_ec_sector = dc_sect[i];
	int temp_cc_sector = cc_sect[i];
	if( temp_ec_sector != 0 ){
	  h_el_sect_etotp_b[temp_ec_sector-1]->Fill(p[i],ectot/p[i]);
	}
	if ( temp_cc_sector != 0 ){
	  int pmt = cc_segm[i]/1000 -1;
	  int segment = cc_segm[i]%1000/10;
	  Float_t thetaCC = 57.2957795*get_thetaCC(tl3_x[i], tl3_y[i], tl3_z[i], tl3_cx[i], tl3_cy[i], tl3_cz[i]);
	  h_el_sect_cctheta_b[temp_cc_sector-1]->Fill(segment,thetaCC);
	}
      }     
      if( e_index[1] >= 0 ){
	if( q[i] > 0 ){
	   double pos_beta = getBeta(i, e_index, sc_sect, dc_sect, sc_t, sc_r, sc_pd, currentrunno, ExpOrSim);
	  h_pos_betap->Fill(p[i], pos_beta);
	  if( dc_sect[i] != 0 ){
	    h_pos_sect_betap[dc_sect[i] - 1]->Fill(p[i], pos_beta);	    
	  }
	}
	else if( q[i] < 0 ){
	  double neg_beta = getBeta(i, e_index, sc_sect, dc_sect, sc_t, sc_r, sc_pd, currentrunno, ExpOrSim);
	  h_neg_betap->Fill(p[i], neg_beta);
	  if( dc_sect[i] != 0 ){
	    h_neg_sect_betap[dc_sect[i] - 1]->Fill(p[i], neg_beta);	    
	  }
	}
      }

    }

    if( e_index[1] >= 0 ){

      //%%%%%%%%%% fill electron plots %%%%%%%%%%%%
      int golden_electron = e_index[1];
      
      if( dc_sect[e_index[1]] != 0 && cc_sect[e_index[1]] != 0 && ec_sect[e_index[1]] != 0 ){
	h_el_p->Fill(p[e_index[1]]);
	h_el_p_pc->Fill(V4_e[1].P());
	h_el_thetaphi->Fill( 180.0/TMath::Pi() * V4_e[1].Theta(), 180.0/TMath::Pi() * V4_e[1].Phi());
	//std::cout << " vz " << vz[e_index[1]] << std::endl;
	h_el_vzphi->Fill( vz[e_index[1]], 180.0/TMath::Pi() * V4_e_uncorr[1].Phi() );
	double vz_corr = getCorrZ(ExpOrSim, vx[e_index[1]], vy[e_index[1]], vz[e_index[1]], p[e_index[1]]*cx[e_index[1]], p[e_index[1]]*cy[e_index[1]], p[e_index[1]]*cz[e_index[1]], dc_sect[e_index[1]]);
	h_el_vzphi_vc->Fill( vz_corr, 180.0/TMath::Pi() * V4_e_uncorr[1].Phi() );
	h_el_vzphi_vcpc->Fill( vz_corr, 180.0/TMath::Pi() * V4_e[1].Phi() );
	h_el_nphe->Fill(nphe[e_index[1]]);
	double ectot = ec_ei[e_index[1]] + ec_eo[e_index[1]];
	h_el_etotp->Fill(p[e_index[1]],ectot/p[e_index[1]]);
	h_el_etotp_pc->Fill( V4_e[1].P(), ectot/V4_e[1].P() );     
	
	//fiducial plots
	h_el_ec->Fill(ech_x[e_index[1]], ech_y[e_index[1]]);
	h_el_dcr1->Fill(tl1_x[e_index[1]], tl1_y[e_index[1]]);
	h_el_dcr3->Fill(tl3_x[golden_electron], tl3_y[golden_electron]);

	Float_t thetaCC = 57.2957795*get_thetaCC(tl3_x[golden_electron], tl3_y[golden_electron], tl3_z[golden_electron], tl3_cx[golden_electron], tl3_cy[golden_electron], tl3_cz[golden_electron]);
	Float_t phi = atan3(cy[golden_electron],cx[golden_electron])*57.2957795;
	      
	//now per sector - shifted from 0 - 6 to -1 to 5, if -1 then no hit in detector sectors
	int dc_sector = dc_sect[e_index[1]]-1;
	int ec_sector = ec_sect[e_index[1]]-1;
	int cc_sector = cc_sect[e_index[1]]-1;
		
	if( dc_sector >= 0 ){
	  h_el_sect_p[dc_sector]->Fill(p[e_index[1]]);
	  h_el_sect_p_pc[dc_sector]->Fill(V4_e[1].P());
	  h_el_sect_thetaphi[dc_sector]->Fill(180.0/TMath::Pi() * V4_e[1].Theta(), 180.0/TMath::Pi() * V4_e[1].Phi());
	}
	if( cc_sector >= 0 ){
	  h_el_sect_nphe[cc_sector]->Fill( nphe[e_index[1]] );
	  int pmt = cc_segm[golden_electron]/1000 -1;
	  int segment = cc_segm[golden_electron]%1000/10;
	  h_el_sect_cctheta[cc_sector]->Fill(segment, thetaCC);
	  //std::cout << " >> phi " << phi << std::endl;
	  h_el_sect_ccphi[cc_sector]->Fill(segment, phi);
	}
	if( ec_sector >= 0 ){
	  h_el_sect_etotp[ec_sector]->Fill(p[e_index[1]],ectot/p[e_index[1]]);
	  h_el_sect_etotp_pc[ec_sector]->Fill( V4_e[1].P(), ectot/V4_e[1].P() );      
	}
      }
      //%%%%%%%%%%%% end electron plots %%%%%%%%%%%%
                 

      //%%%%%%%%%%%% begin hadron analysis %%%%%%%%%%%%

      std::map<int, int> proton_pid = hadronMLEID(proton_id, gpart, e_index, q, p, sc_sect, dc_sect, sc_t, sc_r, sc_pd, pip_vvp_strict, pip_R1fid_strict, pip_MXcut_strict, ExpOrSim, ec_ei, ec_sect, cc_sect, nphe, ec_eo, cx, cy, cz, b, tl1_x, tl1_y, V4_H, currentrunno, pim_vvp_strict, pim_R1fid_strict, pim_MXcut_strict, pr_mean_fit, pip_mean_fit, kp_mean_fit, pr_sigma_fit, pip_sigma_fit, kp_sigma_fit, pr_conf, pip_conf, kp_conf, pr_anticonf, pip_anticonf, kp_anticonf);

      std::map<int, int> kaon_pid = hadronMLEID(kaon_plus_id, gpart, e_index, q, p, sc_sect, dc_sect, sc_t, sc_r, sc_pd, pip_vvp_strict, pip_R1fid_strict, pip_MXcut_strict, ExpOrSim, ec_ei, ec_sect, cc_sect, nphe, ec_eo, cx, cy, cz, b, tl1_x, tl1_y, V4_H, currentrunno, pim_vvp_strict, pim_R1fid_strict, pim_MXcut_strict, pr_mean_fit, pip_mean_fit, kp_mean_fit, pr_sigma_fit, pip_sigma_fit, kp_sigma_fit, pr_conf, pip_conf, kp_conf, pr_anticonf, pip_anticonf, kp_anticonf);
    
      TVector3 V3_kp[2], V3_prot[2]; // 2 for gen(0) and rec(1)
      TLorentzVector V4_kp[2], V4_prot[2];
      int prot_index = proton_pid[2212];
      int kp_index = kaon_pid[321];
      if( prot_index > 0 ){
	V3_prot[1].SetXYZ(p[prot_index]*cx[prot_index], p[prot_index]*cy[prot_index], p[prot_index]*cz[prot_index]);
	V4_prot[1].SetXYZT(V3_prot[1].X(), V3_prot[1].Y(), V3_prot[1].Z(), sqrt(V3_prot[1].Mag2() + pow(prot_mass,2)));
	if(do_momCorr_proton && ExpOrSim == 1) V4_prot[1] = MomCorr->PcorN(V4_prot[1], 1, 2212);
	n_good_pr++;
	//final_pr->SetXYZT(V4_prot[1].X(),V4_prot[1].Y(),V4_prot[1].Z(),V4_prot[1].T()) ;

	double pr_beta = getBeta(prot_index, e_index, sc_sect, dc_sect, sc_t, sc_r, sc_pd, currentrunno, ExpOrSim);	
	h_pr_betap->Fill(V4_prot[1].P(), pr_beta);
	h_pr_thetaphi->Fill(180.0/TMath::Pi() * V4_prot[1].Theta(), 180.0/TMath::Pi() * V4_prot[1].Phi());

	h_pr_thetap_pc->Fill(V4_prot[1].P(), 180.0/TMath::Pi() * V4_prot[1].Theta());
	
	//std::cout << " >> DC SECTOR FOR PROTON " << " "  << " gpart " << gpart << std::endl;       
	if( dc_sect[prot_index] != 0 ){
	  //std::cout << " >> DC SECTOR FOR PROTON " << dc_sect[prot_index] - 1 << std::endl;
	  h_pr_sect_p[dc_sect[prot_index]-1]->Fill(p[prot_index]);
	  h_pr_sect_p_pc[dc_sect[prot_index]-1]->Fill(V4_prot[1].P());
	  h_pr_sect_betap[dc_sect[prot_index]-1]->Fill(p[prot_index], pr_beta);
	  h_pr_sect_betap_pc[dc_sect[prot_index]-1]->Fill(V4_prot[1].P(), pr_beta);
	  h_pr_sect_thetaphi[dc_sect[prot_index]-1]->Fill(180.0/TMath::Pi() * V4_prot[1].Theta(), 180.0/TMath::Pi() * V4_prot[1].Phi());
	}
	
      }
      
      if( kp_index > 0 ){
	V3_kp[1].SetXYZ(p[kp_index]*cx[kp_index], p[kp_index]*cy[kp_index], p[kp_index]*cz[kp_index]);
	V4_kp[1].SetXYZT(V3_kp[1].X(), V3_kp[1].Y(), V3_kp[1].Z(), sqrt(V3_kp[1].Mag2() + pow(kp_mass,2)));
	if(do_momCorr_kaon && ExpOrSim == 1) V4_kp[1] = MomCorr->PcorN(V4_kp[1], 1, 321);
	//final_kp->SetXYZT(V4_kp[1].X(),V4_kp[1].Y(),V4_kp[1].Z(),V4_kp[1].T()) ;


	n_good_kp++;
	double kp_beta = getBeta(kp_index, e_index, sc_sect, dc_sect, sc_t, sc_r, sc_pd, currentrunno, ExpOrSim);
	h_kp_betap->Fill(V4_kp[1].P(), kp_beta);
	h_kp_thetaphi->Fill(180.0/TMath::Pi() * V4_kp[1].Theta(), 180.0/TMath::Pi() * V4_kp[1].Phi());
	
	h_kp_thetap_pc->Fill(V4_kp[1].P(), 180.0/TMath::Pi() * V4_kp[1].Theta());

	if( dc_sect[kp_index] != 0 ){
	  h_kp_sect_p[dc_sect[kp_index]-1]->Fill(p[kp_index]);
	  h_kp_sect_p_pc[dc_sect[kp_index]-1]->Fill(V4_kp[1].P());
	  h_kp_sect_betap[dc_sect[kp_index]-1]->Fill(p[kp_index], kp_beta);
	  h_kp_sect_betap_pc[dc_sect[kp_index]-1]->Fill(V4_kp[1].P(), kp_beta);
	  h_kp_sect_thetaphi[dc_sect[kp_index]-1]->Fill(180.0/TMath::Pi() * V4_kp[1].Theta(), 180.0/TMath::Pi() * V4_kp[1].Phi());
	}	
      }      
    }
    
    //write out good particle lv to a tree and associated conf lvl
    if( n_good_el == 1 && n_good_pr == 1 && n_good_kp == 1 ){
      //output_tree->Fill();
    }
    
  }
  //output_tree->Write();
  //outputfile->Write();
  cout<<endl<<"Done!"<<endl; 
  stopwatch->Print();
  return 0;

}
