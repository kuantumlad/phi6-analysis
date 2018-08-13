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

int hadronBetaTool( int filestart = 1, int Nfiles = 1){

  TStopwatch *stopwatch = new TStopwatch();
  
  TChain *h22chain = new TChain("h22");  
  string firstfilename = "";
  string lastfilename = "";


  // %%%%%%%% read the input files into the TChain %%%%%%%%
  int NtotalFiles = 11625;  
  ifstream filelist;
  filelist.open("programFiles/dataFiles.txt");
  int kStop = Nfiles + filestart;
  if(kStop > NtotalFiles+1) kStop = NtotalFiles+1;
  for(int k = 1; k < kStop; k++){
    string filename;
    filelist>>filename;
    std::cout << " >> ADDING FILE " << filename.c_str() << std::endl;
    if(k >= filestart) h22chain->Add(filename.c_str());
    if(k == filestart) firstfilename = filename;
    if(k == kStop - 1) lastfilename = filename;
  }
  
  firstfilename = firstfilename.substr(36,9); // trim the string down so it's just the run#.subrun# (e.g. 38458.a00)
  lastfilename = lastfilename.substr(36,9);
  cout<<"first file: "<<firstfilename<<" ... last file: "<<lastfilename<<endl;
  // %%%%% end read the input files into a TChain %%%%%

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

  int e_zvertex_strict = 0;
  int e_ECsampling_strict = 0;
  int e_ECoutVin_strict = 0;
  int e_CCthetaMatching_strict = 0;
  int e_ECgeometric_strict = 0;
  int e_R1fid_strict = 0;
  int e_R3fid_strict = 0;
  int e_CCphiMatching_strict = 0;
  int e_CCfiducial_strict = 0;
  int do_momCorr_e = 0;

  MomCorr_e1f *MomCorr = new MomCorr_e1f();

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
  string outfilename;
  TFile *outputfile;
  outputfile = new TFile("hadron_beta_out.root","RECREATE");

  cout<<"entries: "<<h22chain->GetEntries()/1000<<" thousand"<<endl;

  TLorentzVector V4k(0.0, 0.0, Beam_Energy, Beam_Energy); // x, y, z, t
  TLorentzVector V4ISproton(0.0, 0.0, 0.0, prot_mass); // IS = Initial State

  
  int limit = h22chain->GetEntries() / 100;
  std::cout << " LIMIT " << limit << std::endl;
  
  std::cout << " TOTAL ENTRIES TO PROCESS " << h22chain->GetEntries() << std::endl;


  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // HISTOGRAMS
  TH2D *h_el_thetap_pc = new TH2D("h_el_thetap_pc","h_el_thetap_pc", 200, 0.0, Beam_Energy, 200, 0.0, 60.0);
  TH2D *h2_betap_proton = new TH2D("h2_betap_proton","h2_betap_proton",200, 0.01, 5.0, 200, 0.01, 1.1 );
  TH2D *h2_betap_kp = new TH2D("h2_betap_kp","h2_betap_kp",200, 0.01, 5.0, 200, 0.01, 1.1 );
  TH2D *h2_betap_km = new TH2D("h2_betap_km","h2_betap_km",200, 0.01, 5.0, 200, 0.01, 1.1 );

  TH2D *h2_betap_pip = new TH2D("h2_betap_pip","h2_betap_pip",200, 0.01, 5.0, 200, 0.01, 1.1 );
  TH2D *h2_betap_pim = new TH2D("h2_betap_pim","h2_betap_pim",200, 0.01, 5.0, 200, 0.01, 1.1 );

  std::vector<TH2F*> beta_all_pos_sector;
  std::vector<TH2F*> beta_all_neg_sector;
  std::vector<TH2F*> beta_all_pr_sector;
  std::vector<TH2F*> beta_all_pip_sector;
  std::vector<TH2F*> beta_all_kp_sector;
  std::vector<TH2F*> beta_all_pim_sector;
  std::vector<TH2F*> beta_all_km_sector;
  for( int s = 0; s <= 6; s++ ){
    beta_all_pos_sector.push_back( new TH2F(Form("h_beta_all_pos_sector_%d", s),Form("h_beta_all_pos_sector_%d",s), 250, 0.01, 5.0, 250, 0.01, 1.1));
    beta_all_neg_sector.push_back( new TH2F(Form("h_beta_all_neg_sector_%d", s),Form("h_beta_all_neg_sector_%d",s), 250, 0.01, 5.0, 250, 0.01, 1.1));
    beta_all_pr_sector.push_back( new TH2F(Form("h_beta_all_pr_sector_%d", s),Form("h_beta_all_pr_sector_%d",s), 250, 0.01, 5.0, 250, 0.01, 1.1));
    beta_all_pip_sector.push_back( new TH2F(Form("h_beta_all_pip_sector_%d", s),Form("h_beta_all_pip_sector_%d",s), 250, 0.01, 5.0, 250, 0.01, 1.1));
    beta_all_kp_sector.push_back( new TH2F(Form("h_beta_all_kp_sector_%d", s),Form("h_beta_all_kp_sector_%d",s), 250, 0.01, 5.0, 250, 0.01, 1.1));
    beta_all_pim_sector.push_back( new TH2F(Form("h_beta_all_pim_sector_%d", s),Form("h_beta_all_pim_sector_%d",s), 250, 0.01, 5.0, 250, 0.01, 1.1));
    beta_all_km_sector.push_back( new TH2F(Form("h_beta_all_km_sector_%d", s),Form("h_beta_all_km_sector_%d",s), 250, 0.01, 5.0, 250, 0.01, 1.1));
  }

  // %%%%%%%%%%%%% end histograms %%%%%%%%%%%%%%%%%

  for(int i = 0; i < h22chain->GetEntries(); i++){    
    //std::cout << i << std::endl;
    
    h22chain->GetEntry(i);

    if ( i % limit == 0 ){
      double completed = (double)i / (double)h22chain->GetEntries();
      std::cout << " >> COMPLETED " << ceil(completed * 100) << " % --- ON EVENT " << i << std::endl;
    }


    string currentrunno_string = "-123";
    if(ExpOrSim == 1)
      {
	currentrunno_string = h22chain->GetCurrentFile()->GetName();
	currentrunno_string = currentrunno_string.substr(36,5);
    }
    int currentrunno = atoi(currentrunno_string.c_str()); // converts the string to an int
    
    //std::cout << " CURRENT RUN " << currentrunno << std::endl;

    TVector3 V3_e[2]; // 2 for gen(0) and rec(1)
    TLorentzVector V4_e[2];
    
    // %%%%% electron ID %%%%%
    int e_index[2] = {-123,-123};
    //gets the gpart index for the electron
    e_index[1] = eID(gpart, q, p, cc_sect, sc_sect, ec_sect, dc_sect, cx, cy, cz, tl1_x, tl1_y, tl3_x, tl3_y, tl3_z, tl3_cx, tl3_cy, tl3_cz, e_zvertex_strict, vz, vy, vx, e_ECsampling_strict, ExpOrSim, etot, e_ECoutVin_strict, ec_ei, ech_x, ech_y, ech_z, e_CCthetaMatching_strict, cc_segm, e_ECgeometric_strict, e_R1fid_strict, e_R3fid_strict, e_CCphiMatching_strict, sc_pd, e_CCfiducial_strict);
    //%%%%% end electron ID %%%%%    
    
    if(e_index[1] > -122){
      V3_e[1].SetXYZ(p[e_index[1]]*cx[e_index[1]], p[e_index[1]]*cy[e_index[1]], p[e_index[1]]*cz[e_index[1]]);
      V4_e[1].SetXYZT(V3_e[1].X(), V3_e[1].Y(), V3_e[1].Z(), sqrt(V3_e[1].Mag2() + pow(e_mass,2)));
      if(do_momCorr_e && ExpOrSim == 1) V4_e[1] = MomCorr->PcorN(V4_e[1], -1, 11);            
    }

    // %%%%%%% hadron beta vs p for fitting %%%%%
    if( e_index[1] >= 0 ){
      for( int k = 0; k < gpart; k++ ){

	float cand_time = -1.0;
	float beta_measured_pim = -1.0;
	float beta_measured_km = -1.0;
	
	float beta_measured_pr = -1.0;
 	float beta_measured_kp = -1.0;
	float beta_measured_pip = -1.0;
	
	float beta_measured = -1.0;

	int sc_sector = sc_sect[k]; //from 1 to 6

	float corrStartTime = e_sctimeCorr(ExpOrSim, sc_t[e_index[1]], dc_sect[e_index[1]], sc_pd[e_index[1]], currentrunno) - (sc_r[e_index[1]]/speed_of_light);
	float ehtcorrMeasuredTime = h_sctimeCorr(ExpOrSim, sc_t[k], dc_sect[k], sc_pd[k], currentrunno) - corrStartTime;
	beta_measured = (sc_r[k]/ehtcorrMeasuredTime)/speed_of_light;
	  	
	if( sc_sector >= 1 ){
	  //std::cout << " sector " << sc_sector << std::endl;
	  beta_all_pos_sector[sc_sector]->Fill(p[k], beta_measured);	
	}

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
	
	
	int candidate_hadron = -1;       
	float temp_t = 100000.0;
	int temp_c = -1;
	std::map<int, float> beta_map;
	std::vector<int> pid_values;
	pid_values.push_back(211);
	pid_values.push_back(321);
	pid_values.push_back(2212);
	
	for( int c = 0; c < 3; c++ ){
	 
	  float delta_t = fabs(ehtcorrMeasuredTime - (sc_r[k]/(calc_beta_pos[c]*speed_of_light) ));
						   
	  //std::cout << delta_t  << std::endl;
	  beta_map[pid_values[c]] = delta_t;	  
	  
	  if( delta_t < temp_t ){
	    temp_t = delta_t;
	    temp_c = c;
	  }
	}
	
	if ( temp_c >= 0 ){
	//  std::cout << " smallest delta_t " << temp_t << " index " << temp_c << std::endl;
	  candidate_hadron = temp_c;
	}
	
	
	
	//std::cout << " HADRON IS " << candidate_hadron << std::endl;
	
	if( candidate_hadron == 0 ){
	  //std::cout << " PIP  " << std::endl;
	  //std::cout << " sc_sector " << sc_sector << std::endl;
	  beta_all_pip_sector[sc_sector]->Fill(p[k],beta_measured);
	}
	else if( candidate_hadron == 1 ){
	  //std::cout << " KP  " << std::endl;
	  beta_all_kp_sector[sc_sector]->Fill(p[k],beta_measured);
	}
	else if( candidate_hadron == 2 ){
	  //std::cout << " PROTON  " << std::endl;
	  beta_all_pr_sector[sc_sector]->Fill(p[k],beta_measured);	
	}

	beta_map.clear();


      }
    }
  }      


  
  outputfile->Write();
  outputfile->Save();
  
  return 0;
}
