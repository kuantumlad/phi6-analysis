#include <vector>
#include <iostream> 
#include <map>
#include <cmath>


double getBeta(int index, int e_index[], UChar_t sc_sect[], UChar_t dc_sect[], Float_t sc_t[], Float_t sc_r[], UChar_t sc_pd[], int currentrunno, int ExpOrSim ){

  double beta_measured = -1.0;
  Float_t speed_of_light = 29.9792458; // cm/ns

  int k = index;

  int sc_sector = sc_sect[k]; //from 0 to 6
  //std::cout << " SECTOR " << sc_sector << std::endl;
  float corrStartTime = e_sctimeCorr(ExpOrSim, sc_t[e_index[1]], dc_sect[e_index[1]], sc_pd[e_index[1]], currentrunno) - (sc_r[e_index[1]]/speed_of_light);
  float ehtcorrMeasuredTime = h_sctimeCorr(ExpOrSim, sc_t[k], dc_sect[k], sc_pd[k], currentrunno) - corrStartTime;
  beta_measured = (sc_r[k]/ehtcorrMeasuredTime)/speed_of_light;
 
  return beta_measured;

}
