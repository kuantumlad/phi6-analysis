TF1* fitHistogram( TH1D *h_temp ){
      
  //System.out.println(" >> FITTING HISTOGRAM " + h_temp.getName() );
  double xlow, xhigh, histmax;
  int binlow, binhigh, binmax;

  double percentofmax = 0.45;
    
  TF1 *fit = NULL;

  //if( h_temp.getEntries() > 0 ){
  binmax = h_temp->GetMaximumBin();
  histmax = h_temp->GetMaximum();
  binlow = binmax;
  binhigh = binmax;
  
  while( h_temp->GetBinContent(binhigh++) >= percentofmax*histmax && binhigh <= h_temp->GetNbinsX() ){}
  while( h_temp->GetBinContent(binlow--) >= percentofmax*histmax && binlow > 1 ){}
  
  xlow = h_temp->GetXaxis()->GetBinCenter(binlow) - ( h_temp->GetXaxis()->GetBinCenter(binlow) - h_temp->GetXaxis()->GetBinCenter(binlow+1))/2.0; // needs to be low edge, only center now
  xhigh = h_temp->GetXaxis()->GetBinCenter(binhigh+1) - (h_temp->GetXaxis()->GetBinCenter(binhigh+1) - h_temp->GetXaxis()->GetBinCenter(binhigh))/2.0; // needs to be low edge, only center now
												  
  std::cout << " >> values used " << xlow <<  " " << xhigh << " " << histmax << std::endl;
  
  TF1 *fit_temp = new TF1("fit_temp","gaus", xlow, xhigh );
  fit_temp->SetParameter(0, histmax);
  fit_temp->SetParameter(1, h_temp->GetMean() );
  fit_temp->SetParameter(2, h_temp->GetRMS() );
  
  h_temp->Fit(fit_temp, "REQ"); //was only R at first
  fit = fit_temp;  
  
  std::cout << " >> PARAMETER SET " << fit_temp->GetParameter(0) << " " << fit_temp->GetParameter(1) << " " << fit_temp->GetParameter(2) << std::endl;


  return fit;	    
}
