#include "PFlowMonitor.h"

void PFlowMonitor()
{

  std::cout<<"=================="<<std::endl;
  std::cout<<"PFlowAna monitor "<<std::endl;
  std::cout<<"=================="<<std::endl;
  
  //=================================
  // Chose here your config options:
  //=================================
  m_SinglePionLowPerformanceStudies = false;
  m_DijetLowPerformance = false;
  m_DijetSubtraction = true;
  m_Zmumu = false;
  
 
  
  /* WIP  Narrow or Wide eta and pt range */
  if (m_UseNarrowPtRange) _ptRange= {0, 2, 5, 10, 20};
  else _ptRange = {0, 2, 5, 10, 20, 40, 60, 80, 100, 150, 200, 500, 1000};
  if (m_UseNarrowEtaRange) _etaRange= {0, 1, 2, 2.5};
  else _etaRange= {0.0, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.5};
  
  
  m_PrintDebug = false; // Printing message criteria 
  
  if(m_DijetLowPerformance)
    HistFile = new TFile ("");
  else if(m_DijetSubtraction)
    //HistFile = new TFile ("/afs/cern.ch/user/m/moles/pflow/Performance/PFlowAnaPackage/SubTest2/hist-user.moles.mc15_13TeV.361024.JZ4W.test_AOD.95245123.root");
    HistFile = new TFile ("../../../MyDir1/hist-Run.root");
  else if(m_Zmumu)
    HistFile = new TFile ("");
  else std::cout<<" None file found"<<std::endl;

  setStyle();
 
  ChargeShowerSubtraction(15);
  printPS();
}


//////////////////////////////
void ChargeShowerSubtraction(int DeltaR)
{

  std::vector<string>  etaLabelVect = getEtaLabelVect(_etaRange);
  
  std::vector<string> h_EntriesVect = FillHistNameVector("h_Entries_vs_Pull",DeltaR,_etaRange);
  getHist(h_EntriesVect, DeltaR, etaLabelVect);
  std::vector<string> h_CalHitsRemainVect = FillHistNameVector("h_CalHitsRemainingOverPt_vs_Pull",DeltaR,_etaRange);
  getHist(h_CalHitsRemainVect, DeltaR, etaLabelVect);
  std::vector<string> h_CalHitsVect = FillHistNameVector("h_CalHitsOverPt",999,_etaRange);
  getHist(h_CalHitsVect, DeltaR, etaLabelVect);

  Compute_CaloResol(h_CalHitsVect);

  // getFinalHist();

  
  return;
}


std::vector<string> getEtaLabelVect(std::vector<float> EtaRange){

  std::string complete_name = "WrongName";
  std::vector<string> h_LabelVect;
  h_LabelVect.clear();
  
  for (unsigned i_eta = 0; i_eta<EtaRange.size(); i_eta++){
    std::string label; 
    if(i_eta < EtaRange.size()-1) label= (to_string_with_precision(EtaRange.at(i_eta),2)+" #leq |#eta| #leq "+to_string_with_precision(EtaRange.at(i_eta+1),2)).c_str(); 
    else label = ("|#eta| #leq "+to_string_with_precision(EtaRange.at(i_eta),2)).c_str(); 
    h_LabelVect.push_back(label);
  }
  return h_LabelVect; 
}

std::vector<string> FillHistNameVector(std::string name, int DeltaR, std::vector<float> EtaRange){

  std::vector<string> h_NameVect;
  h_NameVect.clear();
  std::string complete_name = "WrongName";
  
  for (unsigned i_eta = 0; i_eta<EtaRange.size(); i_eta++){
    std::string pull_name = ("Pull"+std::to_string(DeltaR)).c_str();
    
    if(DeltaR == 999){
      if (i_eta != EtaRange.size()-1) 
	complete_name = (name + "_eta" + std::to_string((int) ((10 * EtaRange.at(i_eta)))) + "_" + std::to_string((int) ((10 * EtaRange.at(i_eta+1))))).c_str();
      else 
	complete_name = (name + "_eta" +
			 std::to_string((int) ((10 * EtaRange.at(i_eta))))).c_str();
    }
    else{
      if (i_eta != EtaRange.size()-1) {
	complete_name = (name + std::to_string(DeltaR) + "_eta" +
			 std::to_string((int) ((10 * EtaRange.at(i_eta)))) + "_" + std::to_string((int) ((10 * EtaRange.at(i_eta+1))))).c_str();
      }
      else {
	complete_name = (name + std::to_string(DeltaR) + "_eta" +
			 std::to_string((int) ((10 * EtaRange.at(i_eta))))).c_str();
      }
    }
    h_NameVect.push_back(complete_name);
  }
  return h_NameVect; 
}



void getHist(vector<std::string> name, int DeltaR, vector<std::string> etaLabel ){
  
  TCanvas* Can =new TCanvas((name.at(0)).c_str(),(name.at(0)).c_str(),950, 750);
  Can->Divide(3,3);
    
  TH2D* hist;
  for (unsigned i_hist = 0; i_hist<name.size(); i_hist++){

    std::cout<<(name.at(i_hist)).c_str()<<std::endl;
    
    hist = (TH2D*)HistFile->Get((name.at(i_hist)).c_str());
    Can->cd(i_hist+1);
    
    //Different axis label depeding of the histogram
    std::string pull_name = ("Pull"+std::to_string(DeltaR)).c_str();
    string xAxis_name = "CalHitsOverPt";
    if(strstr((name.at(i_hist)).c_str(),xAxis_name.c_str()))
      Draw2DHist(hist,(name.at(i_hist)).c_str(), "p_{T}^{true}[MeV]", "E^{true}_{cl}/p^{true}_{track}");
    else
      Draw2DHist(hist,(name.at(i_hist)).c_str(), "p_{T}^{true}[MeV]", pull_name.c_str());
    
    hist->Draw("COLZ");
    DrawLatex(0.75, 0.80, (etaLabel.at(i_hist)).c_str());
  }
    
  m_CanVect.push_back(Can);
  return; 
}



void Compute_CaloResol(vector<std::string> name)
{
  
  gStyle->SetPadGridX(false);
  gStyle->SetPadGridY(false);
  gStyle->SetOptLogx(false);
  
  TH2D* hist;
  for (unsigned i_hist = 0; i_hist<name.size(); i_hist++){
    hist = (TH2D*)HistFile->Get((name.at(i_hist)).c_str());

    TCanvas* Can = new TCanvas(("CanCaloResol_"+std::to_string(i_hist)).c_str(),("CanCaloResol_"+std::to_string(i_hist)).c_str(),950, 750);
    Can->Divide(5,5);
    int NBin_X = hist->GetNbinsX();
    std::vector<float> EffectiveCaloResolution; // store mu and sigma
    EffectiveCaloResolution.resize(NBin_X);

    //Gausian fit over y for each x bin
    for (int i=0;i< NBin_X;i++){
      Can->cd(i+1);
      TH1D* ProjectionY = hist->ProjectionY("hist_px",i,i+1,"");
      if (i<11)ProjectionY ->Rebin(2);
      ProjectionY ->DrawCopy("hist");
      
      float mean = FitGaussHist(ProjectionY , 1, 0);
      float sigma = FitGaussHist(ProjectionY , 0, 1);
      
      //Fill the vector 
      EffectiveCaloResolution.at(i) = sigma/mean;
  }
    m_CanVect.push_back(Can);
  }
  return;
}


  
void  Compute_EafterSub(TH2D *h_EAfter, TH2D *h_Entries){
  // Calculate the ratio between both histograms


  
  // //Calculate E_true_after 
  // TProfile2D *hist2=hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_04->Clone();
  // hist2->Sumw2();
  // hist2->Divide(hTurnOff_Entries_vs_Pull_eta_00_04);

  
  // TCanvas* Can = new TCanvas(hist->GetName(),hist->GetName(),950, 750);
  // Can->Divide(5,5);
  // int NBin_X = hist->GetNbinsX();
  // std::vector<float> EffectiveCaloResolution; // store mu and sigma
  // EffectiveCaloResolution.resize(NBin_X);
  
  // for (int i=0;i< NBin_X;i++){
  //   Can->cd(i+1);
  //   TH1D* ProjectionY = hist->ProjectionY("hist_px",i,i+1,"");
  //   ProjectionY ->Rebin(2);
  //   ProjectionY ->DrawCopy("hist");
    
  //   float mean = FitGaussHist(ProjectionY , 1, 0);
  //   float sigma = FitGaussHist(ProjectionY , 0, 1);
    
  //   //Fill the vector 
  //   EffectiveCaloResolution.at(i) = sigma/mean;
  // }
  // m_CanVect.push_back(Can);
  return;
}


void getFinalHist(){
  
  // int n_binX =25;
  // double ptTrack_trueRange[26] = {500., 650.,850.,1100., 1400., 1800., 2300., 2950., 3800., 4800., 6100., 7700., 9300., 11300., 13800., 16800., 20800., 25800., 31800.,
  // 				  39800., 50000., 60000., 70000., 80000., 90000., 100000.,};
  // int n_binY =25;
  // double PullRange[26] = {-5.0,-4.0,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0,10.,12.,14.,16.,20.,25.};
  // hFinal= new TH2D("hFinal","hFinal",n_binX,ptTrack_trueRange,n_binY, PullRange);
  
  // for (int x = 0; x< hFinal->GetNbinsX();x++){
  //   for (int y = 0; y< hFinal->GetNbinsY();y++){
  //     double E_after = hist2->GetBinContent(x+1,y+1);
  //     hFinal->SetBinContent(x+1,y+1,E_after/EffectiveCaloResolution[x]);
  //   }
  // }

  // //for different eta ranges!
  // TCanvas* Can_Subtraction = new TCanvas("Can_Subtraction","Can_Subtraction",900, 800);
  // Can_Subtraction->Divide(1,1);
  // Can_Subtraction->cd(1);
  // gPad->SetLogx();
  // hFinal->Draw("COLZ");
  

  
  return; 
}

//////////////////////////////
void Draw2DHist(TH2D *hist, const char* title, const char* Xaxis, const char* Yaxis)
{   
  hist ->SetStats(kFALSE);
  hist ->SetTitle(title);
  hist ->GetXaxis()->SetTitle(Xaxis);
  hist ->GetYaxis()->SetTitle(Yaxis);
  
  return;
}
//////////////////////////////
void DrawLatex( float x, float y, const char* label){
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextFont(46);
  // t->SetTextSizePixels(20);
  t->DrawLatex(x,y,label);
  
  return;
} 

//////////////////////////////
float FitGaussHist( TH1D* hist, bool GetMean, bool GetSigma){

  
  float Peak = hist->GetBinCenter(hist->GetMaximumBin()); 
  float mean = hist->GetMean();
  float RMS = hist->GetRMS();
  TF1 *fitg = new TF1("fitg","gaus",mean-1.5*RMS,mean+1.5*RMS);
  fitg->SetLineWidth(1);
  fitg->SetLineColor(kGreen);
  hist->Fit(fitg,"R+");
  hist->DrawCopy("hist");
  fitg->Draw("same");
  
  double p2 = fitg->GetParameter(2); double e_p2 = fitg->GetParError(2);   
  double p1 = fitg->GetParameter(1); double e_p1 = fitg->GetParError(1);
  float chi2 = fitg->GetChisquare();
  float DeltaPar =  fabs(p1-Peak); 
  
  float f_p2 = p2; 
  float f_p1 = p1;
  float f_chi2 = chi2;
  int IterNum = 0;
  
  while (0.01 < DeltaPar && IterNum <= 20 ){
    cout<<"Initial DeltaPar"<<DeltaPar<<endl; 
    f_p2 = p2;
    f_p1 = p1;
    f_chi2 = chi2;
    cout<<"Parameters p1="<<p1<<"  p2="<<p2<<endl;
    TF1 *fitg2 = new TF1("fitg2","gaus",p1-2*p2,p1+2*p2);
    fitg2->SetLineWidth(1);
    fitg2->SetLineColor(kGreen);
    hist->Fit(fitg2,"R");
    hist->DrawCopy("hist");
    fitg2->Draw("same");
    
    float new_p2 = fitg2->GetParameter(2);   
    float new_p1 = fitg2->GetParameter(1);
    float new_ep2 = fitg2->GetParError(2);   
    float new_ep1 = fitg2->GetParError(1);
    float new_chi2 =  fitg2->GetChisquare();
    DeltaPar = fabs(p1-new_p1);
    cout<<"Final DeltaPar"<<DeltaPar<<endl; 
    cout<<"new_p1="<<new_p1<<"  new_p2="<<new_p2<<endl;
    p2 = new_p2;
    p1 = new_p1;
    e_p2 = new_ep2;
    e_p1 = new_ep1;
    chi2 = new_chi2;
    IterNum++;
  } 

 

  TPaveText *pav = new TPaveText(0.12, 0.60, 0.42,0.85,"NDC");
  pav->SetBorderSize(0);  pav->SetFillColorAlpha(0,0.);
  char texto1[100]; 
  if(IterNum == 20) pav->AddText("Fit failed");
  else{
    sprintf(texto1, "  #mu = %2.3f #pm %2.3f", f_p1,e_p1);
    pav->AddText(texto1);
    sprintf(texto1, " #sigma = %2.3f #pm %2.3f", f_p2,e_p2);
    pav->AddText(texto1);
    sprintf(texto1, " #chi^{2} = %2.3f", f_chi2);
    pav->AddText(texto1);
  }
  pav->Draw();

  if(GetMean) return p1;
  if (GetSigma) return p2;
  return 999;
}

////////////////////////////
void printPS() {
  
  cout << "\nStoring the plots in a ps file..." << endl;
  Char_t filename[] = "Subtraction_plots.ps";
  Char_t command[] = "";
  Char_t name[] = "";
  
  TCanvas c0;
  sprintf(command,"%s[",filename);
  c0.Print(command);
  
  // Add canvas
  for (int i=0; i<m_CanVect.size();i++){
    m_CanVect.at(i)->Print(filename);
  }
  
  sprintf(command,"%s]",filename);
  c0.Print(command);
  
  // Compress it!
  sprintf(name,".!gzip -f %s",filename);
  gROOT->ProcessLine(name);
  
  cout << " - Plots stored successfully!" << endl;
  
  return;
}

/////////////////////
void setStyle() {
  
  // Extra options to the ATLAS Style
  Float_t Labelsize = 0.04;
  Float_t Titlesize = 0.05;
  
  gStyle->SetLabelSize(Labelsize,"x");
  gStyle->SetLabelSize(Labelsize,"y");
  gStyle->SetLabelSize(Labelsize-0.01,"z");
  
  gStyle->SetTitleSize(Titlesize,"x");
  gStyle->SetTitleOffset(0.95,"x");
  gStyle->SetTitleSize(Titlesize,"y");
  gStyle->SetTitleOffset(0.75,"y");
  gStyle->SetTitleSize(Titlesize,"z");
  gStyle->SetTextFont(42);
  
  gStyle->SetStripDecimals(false);    
  TGaxis::SetMaxDigits(4);
  gStyle->SetPalette(1);

  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);
  gStyle->SetGridColor(10);

  gStyle->SetOptLogx(true);
  
  gROOT->ForceStyle();
  return; 
}


std::string to_string_with_precision(float value, int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << value;
    return out.str();
}

  

  // /*Construct final histogram*/
  // /*
  // double TurnOffBins[20]={500.,631.,794.,1000.,1259.,1585.,1995.,2512.,3162.,3981.,5012.,6310.,7943.,10000.,12589.,15849.,19953.,25119.,31623.,40000.};
  // double TurnOffBins2[26]={-5.0,-4.0,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0,10.,12.,14.,16.,20.,25.};
  // hFinal= new TH2D("hFinal","hFinal",19,TurnOffBins,25,TurnOffBins2);
  
  // for (int x = 0; x< hFinal->GetNbinsX();x++){
  //   for (int y = 0; y< hFinal->GetNbinsY();y++){
  //     double E_after = hist2->GetBinContent(x+1,y+1);
  //     hFinal->SetBinContent(x+1,y+1,E_after/EffectiveCaloResolution[x]);
  //   }
  // }

  // //hFinal
  // TCanvas* Can_Subtraction = new TCanvas("Can_Subtraction","Can_Subtraction",900, 800);
  // Can_Subtraction->Divide(1,1);
  // Can_Subtraction->cd(1);
  // gPad->SetLogx();
  // hFinal->Draw("COLZ");
  // */


  
  // /*
  // //For different eta ranges eta<1; 1<eta<2; 2<eta<2.5
  
  // TCanvas* Can_hTurnOf_etaRegion = new TCanvas("Can_hTurnOff_etaRegion","Can_hTurnOff_etaRegion",950, 750);
  // Can_hTurnOff_etaRegion->Divide(3,3);
  
  // Can_hTurnOff_etaRegion->cd(1); gPad->SetLogx();
  // TH2D* hTurnOff_CalHitsOverPt_eta_00_10 = (TH2D*)HistFile->Get("hTurnOff_CalHitsOverPt_eta_00_10_hist");
  // Draw2DHist(hTurnOff_CalHitsOverPt_eta_00_10,"CalHitsOverPt", "p_{T}^{true}", "E^{true}_{cl}/p^{true}_{track}");
  // hTurnOff_CalHitsOverPt_eta_00_10->Draw("COLZ");
  // DrawLatex("|#eta|<1.0"); 
  
  // Can_hTurnOf_etaRegion->cd(2); gPad->SetLogx();
  // TProfile2D* hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_10 = (TProfile2D*)HistFile->Get("hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_10_hist");
  // Draw2DHist(hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_10,"CalHitsRemainingOverPt", "p_{T}^{true}", "E^{true}_{cl}(after)/p^{true}_{track}");
  // hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_10 ->DrawCopy("COLZ");
  // DrawLatex("|#eta|<1.0");
  
  // Can_hTurnOf_etaRegion->cd(3);gPad->SetLogx();
  // TH2D* hTurnOff_Entries_vs_Pull_eta_00_10 = (TH2D*)HistFile->Get("hTurnOff_Entries_vs_Pull_eta_00_10_hist");
  // Draw2DHist(hTurnOff_Entries_vs_Pull_eta_00_10,"Entries_vs_Pull15", "p_{T}^{true}", "pull15");
  // hTurnOff_Entries_vs_Pull_eta_00_10->Draw("COLZ");
  // DrawLatex("|#eta|<1.0");
  // */
  

//   return;
// }



 



// void PrintHistosCode()
// //////////////////////////////
// {
//   cout << " ** PrintHistosCode ** # elements " << HistogramList.size() << endl;

//   // print histogram color code
//   float Legend_lx = 0.60;
//   float Legend_ux = 0.88;
//   float Legend_ly = 0.66;
//   float Legend_uy = 0.76;

//   TLegend *Rotllo = new TLegend(Legend_lx, Legend_ly, Legend_ux, Legend_uy);
//   Rotllo->SetFillStyle(3005);
//   Rotllo->SetFillColor(10);
//   char QueDir[80];
//   for (int i=0; i < HistogramList.size(); i++) {
//     sprintf(QueDir,"hola %d", i);
//     if (i==0) sprintf(QueDir,"Truth #pi^{+}");
//     if (i==1) sprintf(QueDir,"Truth #pi^{+} (p_{T}>400MeV &  |#eta<=2.5|");
//     if (i==2) sprintf(QueDir,"PFO");
//     Rotllo->AddEntry(HistogramList.at(i), QueDir);
//   }
//   Rotllo->Draw();

//   return;
// }





