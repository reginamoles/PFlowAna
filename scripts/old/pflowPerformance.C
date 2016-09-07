//=======================================================================================================
// This macro was originally created to draw histograms for PFO performance study (from HISTOS ntuple)
//=======================================================================================================

#include <vector>
#include <cmath>
using namespace std;
//
void InitStyle();
void colorins();

int color_histos = kAzure-5;
int truth_color = kBlue+2;
int truthCuts_color = kGreen+2;
int pfo_color = kYellow;
void PrintHistosCode();

void Draw2DHist(TH2D* hist, const char*, const char*, const char*);
void DrawLatex(const char*);
float FitGaussHist(TH1D*,bool,bool);
  
void kinematic();


TFile *HistFile;
//TFile *HistFileBug; TFile *HistFileNoBugFixed;


std::vector<TH1F*> HistogramList;


void pflowPerformance(TFile* MyFile=NULL)
{

  InitStyle();
  colorins();
  
  // load histograms file:
  if (MyFile != NULL) HistFile = MyFile;
  else {
    //Di-jet sample (calibration Hit fix)
    //HistFile = new TFile ("/Users/moles/cern/testarea/BarcodeBugTest/CalibHitFixed/user.moles.CalibHitFixed.root");
    HistFileBug = new TFile ("/../../../MyDir1/hist-Run.root");
  }
  
  ChargeShowerSubtraction();
  
  
  //EfficiencyLead();
  //PurityLead();
  //NCluster09Lead();
  //DeltaE_09_07();
  // kinematic();
  // TrackExtrapolator();
  // Purity();
  // Efficiency();
  // NCluster09();
  // DeltaR_EM2();
    
}



//////////////////////////////
void ChargeShowerSubtraction()
{
  TCanvas* Can_hTurnOf_eta04 = new TCanvas("Can_hTurnOff_eta04","Can_hTurnOff_eta04",950, 350);
  Can_hTurnOff_eta04->Divide(3,1);
  
  Can_hTurnOff_eta04->cd(1);
  gPad->SetLogx();
  TH2D* hTurnOff_CalHitsOverPt_eta_00_04 = (TH2D*)HistFile->Get("hTurnOff_CalHitsOverPt_eta_00_04_hist");
  Draw2DHist(hTurnOff_CalHitsOverPt_eta_00_04,"CalHitsOverPt", "p_{T}^{true}", "E^{true}_{cl}/p^{true}_{track}");
  hTurnOff_CalHitsOverPt_eta_00_04->Draw("COLZ");
  DrawLatex("|#eta|<0.4"); 
  
  //test barcode bug
  //TH2D* hTurnOff_CalHitsOverPt_eta_00_04_Bug = (TH2D*)HistFileBug->Get("hTurnOff_CalHitsOverPt_eta_00_04_hist");
  //TH2D* hTurnOff_CalHitsOverPt_eta_00_04_NoBugFixed = (TH2D*)HistFileNoBugFixed->Get("hTurnOff_CalHitsOverPt_eta_00_04_hist");

  Can_hTurnOf_eta04->cd(2);
  gPad->SetLogx();
  TProfile2D* hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_04 = (TProfile2D*)HistFile->Get("hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_04_hist");
  Draw2DHist(hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_04,"CalHitsRemainingOverPt", "p_{T}^{true}", "E^{true}_{cl}(after)/p^{true}_{track}");
  hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_04 ->Draw("COLZ");
  DrawLatex("|#eta|<0.4");
  
  Can_hTurnOf_eta04->cd(3);
  gPad->SetLogx();
  TH2D* hTurnOff_Entries_vs_Pull_eta_00_04 = (TH2D*)HistFile->Get("hTurnOff_Entries_vs_Pull_eta_00_04_hist");
  Draw2DHist(hTurnOff_Entries_vs_Pull_eta_00_04,"Entries_vs_Pull15", "p_{T}^{true}", "pull15");
  hTurnOff_Entries_vs_Pull_eta_00_04->Draw("COLZ");
  DrawLatex("|#eta|<0.4");
  
  /*Gausian fit over y for each x bin*/
  TCanvas* CanFit = new TCanvas("CanFit","CanFit",950, 900);
  CanFit->Divide(5,4);
  double EffectiveCaloResolution[20]; // store mu and sigma
  int NBin_X = hTurnOff_CalHitsOverPt_eta_00_04->GetNbinsX();
  
  for (int i=0;i< NBin_X;i++){
    EffectiveCaloResolution[i] = 0;
  
    CanFit->cd(i+1);
    TH1D* hist = hTurnOff_CalHitsOverPt_eta_00_04->ProjectionY("hist_px",i,i+1,"");
    hist->DrawCopy("hist");

    float mean = FitGaussHist(hist, 1, 0);
    float sigma = FitGaussHist(hist, 0, 1);
    
    //Fill a vector with sigma/mu:
    EffectiveCaloResolution[i] = sigma/mean;
  }

   /*
  //Calculate E_true_after 
  //PROBLEM IN THE NORMALIZATION of E_after
  TProfile2D *hist2=hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_04->Clone();
  hist2->Sumw2();
  hist2->Divide(hTurnOff_Entries_vs_Pull_eta_00_04);
  */
  

  /*Construct final histogram*/
  /*
  double TurnOffBins[20]={500.,631.,794.,1000.,1259.,1585.,1995.,2512.,3162.,3981.,5012.,6310.,7943.,10000.,12589.,15849.,19953.,25119.,31623.,40000.};
  double TurnOffBins2[26]={-5.0,-4.0,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0,10.,12.,14.,16.,20.,25.};
  hFinal= new TH2D("hFinal","hFinal",19,TurnOffBins,25,TurnOffBins2);
  
  for (int x = 0; x< hFinal->GetNbinsX();x++){
    for (int y = 0; y< hFinal->GetNbinsY();y++){
      double E_after = hist2->GetBinContent(x+1,y+1);
      hFinal->SetBinContent(x+1,y+1,E_after/EffectiveCaloResolution[x]);
    }
  }

  //hFinal
  TCanvas* Can_Subtraction = new TCanvas("Can_Subtraction","Can_Subtraction",900, 800);
  Can_Subtraction->Divide(1,1);
  Can_Subtraction->cd(1);
  gPad->SetLogx();
  hFinal->Draw("COLZ");
  */


  
  /*
  //For different eta ranges eta<1; 1<eta<2; 2<eta<2.5
  
  TCanvas* Can_hTurnOf_etaRegion = new TCanvas("Can_hTurnOff_etaRegion","Can_hTurnOff_etaRegion",950, 750);
  Can_hTurnOff_etaRegion->Divide(3,3);
  
  Can_hTurnOff_etaRegion->cd(1); gPad->SetLogx();
  TH2D* hTurnOff_CalHitsOverPt_eta_00_10 = (TH2D*)HistFile->Get("hTurnOff_CalHitsOverPt_eta_00_10_hist");
  Draw2DHist(hTurnOff_CalHitsOverPt_eta_00_10,"CalHitsOverPt", "p_{T}^{true}", "E^{true}_{cl}/p^{true}_{track}");
  hTurnOff_CalHitsOverPt_eta_00_10->Draw("COLZ");
  DrawLatex("|#eta|<1.0"); 
  
  Can_hTurnOf_etaRegion->cd(2); gPad->SetLogx();
  TProfile2D* hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_10 = (TProfile2D*)HistFile->Get("hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_10_hist");
  Draw2DHist(hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_10,"CalHitsRemainingOverPt", "p_{T}^{true}", "E^{true}_{cl}(after)/p^{true}_{track}");
  hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_10 ->DrawCopy("COLZ");
  DrawLatex("|#eta|<1.0");
  
  Can_hTurnOf_etaRegion->cd(3);gPad->SetLogx();
  TH2D* hTurnOff_Entries_vs_Pull_eta_00_10 = (TH2D*)HistFile->Get("hTurnOff_Entries_vs_Pull_eta_00_10_hist");
  Draw2DHist(hTurnOff_Entries_vs_Pull_eta_00_10,"Entries_vs_Pull15", "p_{T}^{true}", "pull15");
  hTurnOff_Entries_vs_Pull_eta_00_10->Draw("COLZ");
  DrawLatex("|#eta|<1.0");
  */
  

  return;
}






//////////////////////////////
void kinematic()
{
  //////////////////////////////
  // Testing variables
  //////////////////////////////
  TCanvas* Can_Kinematics = new TCanvas("kinematics","kinematics",900, 450);
  Can_Kinematics->Divide(3,1);

  TLegend *MyLegend = new TLegend (0.45, 0.70, 0.88, 0.85);
  MyLegend->SetFillStyle(1001);
  MyLegend->SetFillColor(10); 

  TH1 *h_mcPt =  (TH1F*)HistFile->Get("h_mcPt");
  h_mcPt->SetTitle("Transverse momentum #pi^{+}");
  h_mcPt ->GetXaxis()->SetTitle("p_{T} [GeV]");
  h_mcPt->SetFillColor(truth_color);
  h_mcPt->SetFillStyle(3002);
  
  TH1 *h_mcPtCut =  (TH1F*)HistFile->Get("h_mcPtCut");
  h_mcPtCut->SetFillColor(truthCuts_color);
 
  TH1 *h_cpfoPt =  (TH1F*)HistFile->Get("h_cpfoPt");
  h_cpfoPt->SetFillColor(pfo_color);
  h_cpfoPt->SetFillStyle(1001);
  //h_cpfoPt->SetFillStyle(3002);
   
  TH1 *h_mcEta =  (TH1F*)HistFile->Get("h_mcEta");
  h_mcEta->SetTitle(" #eta #pi^{+}");
  h_mcEta ->GetXaxis()->SetTitle("#eta");
  //h_mcEta->GetYaxis()->SetRangeUser(0.,50.);
  h_mcEta->SetFillColor(truth_color);
  h_mcEta->SetFillStyle(3002);
  
  TH1 *h_mcEtaCut =  (TH1F*)HistFile->Get("h_mcEtaCut");
  h_mcEtaCut->SetFillColor(truthCuts_color);
  
  TH1 *h_cpfoEta =  (TH1F*)HistFile->Get("h_cpfoEta");
  h_cpfoEta->SetFillColor(pfo_color);
  h_cpfoEta->SetFillStyle(1001);
  //h_cpfoEta->SetFillStyle(3002);

  //draw histograms
  Can_Kinematics->cd(1);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  h_mcPt->DrawCopy("hist");
  h_cpfoPt->DrawCopy("hist,same");
  h_mcPtCut->DrawCopy("hist,same");

  MyLegend->AddEntry( h_mcPt, "Truth #pi^{+}");
  MyLegend->AddEntry( h_mcPtCut, "Truth #pi^{+} (0.4GeV<p_{T}<40GeV & |#eta|<2.5");
  MyLegend->AddEntry( h_cpfoPt, "PFO","f");
  MyLegend->Draw();
  
  //draw histograms
  Can_Kinematics->cd(2);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  //h_mcPt->GetXaxis()->SetRangeUser(0.,100.);
  //h_mcPt->GetYaxis()->SetRangeUser(0.,100.);
  h_mcPt->SetTitle("Transverse momentum #pi^{+} (ZOOM)");
  h_mcPt->Draw("hist");
  h_cpfoPt->Draw("hist,same");
  h_mcPtCut->Draw("hist,same");
  MyLegend->Draw();
  
  Can_Kinematics->cd(3);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  h_mcEta->Draw("hist");
  h_cpfoEta->Draw("hist,same");
  h_mcEtaCut->Draw("hist,same");
  MyLegend->Draw();


  TCanvas* Can_HadrInter = new TCanvas("Can_HadrInter","Can_HadrInter",900, 450);
  Can_HadrInter->Divide(2,1);
  
  TH1 *h_VtxR =  (TH1F*)HistFile->Get("h_VtxR_PIX");
  h_VtxR->SetFillColor(kOrange+7);
  h_VtxR->SetTitle("");
  Can_HadrInter->cd(1);
  gPad->SetGridx(1);
  h_VtxR ->GetXaxis()->SetTitle("Vertex Radius [mm]");
  h_VtxR->Draw("hist");
  
  
  TH2 *h_VtxXY =  (TH2F*)HistFile->Get("h_VtxXY_PIX");
  h_VtxXY->SetFillColor(truthCuts_color);
  Can_HadrInter->cd(2);
  h_VtxXY->Draw("COLZ");
  
  //Can_Kinematics->Print("kinematics.png");

}

//===========================================
// Track Extrapolation
// (eta,phi parameters at different layers)
//===========================================

void TrackExtrapolator(){

  //Eta parameter

  TLegend *MyEtaLegend = new TLegend (0.70, 0.50, 0.88, 0.85);
  MyEtaLegend->SetFillStyle(1001);
  MyEtaLegend->SetFillColor(10); 
  
  TCanvas* Can_Eta = new TCanvas("Can_Eta","Can_Eta",900, 800);
  Can_Eta->Divide(2,2);
  
  Can_Eta->cd(1);
  gPad->SetGridx(1); gPad->SetGridy(1);
  TH1 *h_Eta_EMB1_OnlyBa =  (TH1F*)HistFile->Get("h_Eta_EMB1_OnlyBa");
  h_Eta_EMB1_OnlyBa->SetStats(kFALSE);
  h_Eta_EMB1_OnlyBa->SetLineColor(kBlue);
  h_Eta_EMB1_OnlyBa->SetTitle("Only Barrel Tracks (EMB)");
  h_Eta_EMB1_OnlyBa->GetXaxis()->SetTitle("#eta");
  //h_Eta_EMB1_OnlyBa->GetYaxis()->SetRangeUser(0.,80.);
  h_Eta_EMB1_OnlyBa->Draw("hist");
  TH1 *h_Eta_EMB2_OnlyBa =  (TH1F*)HistFile->Get("h_Eta_EMB2_OnlyBa");
  h_Eta_EMB2_OnlyBa->SetLineColor(kGreen+2);
  h_Eta_EMB2_OnlyBa->Draw("histsame");
  TH1 *h_Eta_EMB3_OnlyBa =  (TH1F*)HistFile->Get("h_Eta_EMB3_OnlyBa");
  h_Eta_EMB3_OnlyBa->SetLineColor(kRed);
  h_Eta_EMB3_OnlyBa->Draw("histsame");

  char QueDir[80];
  sprintf(QueDir,"%s %d", "EMB1:", h_Eta_EMB1_OnlyBa->Integral());
  MyEtaLegend->AddEntry( h_Eta_EMB1_OnlyBa, QueDir, "l");
  sprintf(QueDir,"%s %d", "EMB2:", h_Eta_EMB2_OnlyBa->Integral());
  MyEtaLegend->AddEntry( h_Eta_EMB2_OnlyBa, QueDir, "l");
  sprintf(QueDir,"%s %d", "EMB3:", h_Eta_EMB3_OnlyBa->Integral());
  MyEtaLegend->AddEntry( h_Eta_EMB3_OnlyBa, QueDir, "l");
  MyEtaLegend->Draw();
  

  TLegend *MyEtaECLegend = new TLegend (0.70, 0.50, 0.88, 0.85);
  MyEtaECLegend->SetFillStyle(1001);
  MyEtaECLegend->SetFillColor(10); 
  
  Can_Eta->cd(2);
  gPad->SetGridx(1); gPad->SetGridy(1);
  TH1 *h_Eta_EME1_OnlyEC =  (TH1F*)HistFile->Get("h_Eta_EME1_OnlyEC");
  h_Eta_EME1_OnlyEC->SetStats(kFALSE);
  h_Eta_EME1_OnlyEC->SetLineColor(kBlue);
  h_Eta_EME1_OnlyEC->SetTitle("Only EC Tracks (EME)");
  h_Eta_EME1_OnlyEC->GetXaxis()->SetTitle("#eta");
  //h_Eta_EME1_OnlyEC->GetYaxis()->SetRangeUser(0,80);
  h_Eta_EME1_OnlyEC->Draw("hist");
  TH1 *h_Eta_EME2_OnlyEC =  (TH1F*)HistFile->Get("h_Eta_EME2_OnlyEC");
  h_Eta_EME2_OnlyEC->SetLineColor(kGreen+2);
  h_Eta_EME2_OnlyEC->Draw("histsame");
  TH1 *h_Eta_EME3_OnlyEC =  (TH1F*)HistFile->Get("h_Eta_EME3_OnlyEC");
  h_Eta_EME3_OnlyEC->SetLineColor(kRed);
  h_Eta_EME3_OnlyEC->Draw("histsame");

  sprintf(QueDir,"%s %d", "EME1:", h_Eta_EME1_OnlyEC->Integral());
  MyEtaECLegend->AddEntry( h_Eta_EME1_OnlyEC, QueDir, "l");
  sprintf(QueDir,"%s %d", "EME2:", h_Eta_EME2_OnlyEC->Integral());
  MyEtaECLegend->AddEntry( h_Eta_EME2_OnlyEC, QueDir, "l");
  sprintf(QueDir,"%s %d", "EME3:", h_Eta_EME3_OnlyEC->Integral());
  MyEtaECLegend->AddEntry( h_Eta_EME3_OnlyEC, QueDir, "l");
  MyEtaECLegend->Draw();

  Can_Eta->cd(3);
  gPad->SetGridx(1); gPad->SetGridy(1);
  TH1 *h_Eta_Tile1_OnlyBa =  (TH1F*)HistFile->Get("h_Eta_Tile1_OnlyBa");
  h_Eta_Tile1_OnlyBa->SetStats(kFALSE);
  h_Eta_Tile1_OnlyBa->SetLineColor(kBlue);
  h_Eta_Tile1_OnlyBa->SetTitle("Only Barrel Tracks (Tile)");
  h_Eta_Tile1_OnlyBa->GetXaxis()->SetTitle("#eta");
  // h_Eta_Tile1_OnlyBa->GetYaxis()->SetRangeUser(0,80);
  h_Eta_Tile1_OnlyBa->Draw("hist");
  TH1 *h_Eta_Tile2_OnlyBa =  (TH1F*)HistFile->Get("h_Eta_Tile2_OnlyBa");
  h_Eta_Tile2_OnlyBa->SetLineColor(kGreen+2);
  h_Eta_Tile2_OnlyBa->Draw("histsame");
  TH1 *h_Eta_Tile3_OnlyBa =  (TH1F*)HistFile->Get("h_Eta_Tile3_OnlyBa");
  h_Eta_Tile3_OnlyBa->SetLineColor(kRed);
  h_Eta_Tile3_OnlyBa->Draw("histsame");

  TLegend *MyEtaTileLegend = new TLegend (0.70, 0.50, 0.88, 0.85);
  MyEtaTileLegend->SetFillStyle(1001);
  MyEtaTileLegend->SetFillColor(10); 
  sprintf(QueDir,"%s %d", "Tile1:", h_Eta_Tile1_OnlyBa->Integral());
  MyEtaTileLegend->AddEntry( h_Eta_Tile1_OnlyBa, QueDir, "l");
  sprintf(QueDir,"%s %d", "Tile2:", h_Eta_Tile2_OnlyBa->Integral());
  MyEtaTileLegend->AddEntry( h_Eta_Tile2_OnlyBa, QueDir, "l");
  sprintf(QueDir,"%s %d", "Tile3:",  h_Eta_Tile3_OnlyBa->Integral());
  MyEtaTileLegend->AddEntry( h_Eta_Tile3_OnlyBa, QueDir, "l");
  MyEtaTileLegend->Draw();

  
  Can_Eta->cd(4);
  gPad->SetGridx(1); gPad->SetGridy(1);
  TH1 *h_Eta_HEC1_OnlyEC =  (TH1F*)HistFile->Get("h_Eta_HEC1_OnlyEC");
  h_Eta_HEC1_OnlyEC->SetStats(kFALSE);
  h_Eta_HEC1_OnlyEC->SetLineColor(kBlue);
  h_Eta_HEC1_OnlyEC->SetTitle("Only EC Tracks (HEC)");
  h_Eta_HEC1_OnlyEC->GetXaxis()->SetTitle("#eta");
  //h_Eta_HEC1_OnlyEC->GetYaxis()->SetRangeUser(0,80);
  h_Eta_HEC1_OnlyEC->Draw("hist");
  TH1 *h_Eta_HEC2_OnlyEC =  (TH1F*)HistFile->Get("h_Eta_HEC2_OnlyEC");
  h_Eta_HEC2_OnlyEC->SetLineColor(kGreen+2);
  h_Eta_HEC2_OnlyEC->Draw("histsame");
  TH1 *h_Eta_HEC3_OnlyEC =  (TH1F*)HistFile->Get("h_Eta_HEC3_OnlyEC");
  h_Eta_HEC3_OnlyEC->SetLineColor(kOrange);
  h_Eta_HEC3_OnlyEC->Draw("histsame");
  TH1 *h_Eta_HEC4_OnlyEC =  (TH1F*)HistFile->Get("h_Eta_HEC4_OnlyEC");
  h_Eta_HEC4_OnlyEC->SetLineColor(kRed);
  h_Eta_HEC4_OnlyEC->Draw("histsame");

  TLegend *MyEtaECTileLegend = new TLegend (0.70, 0.50, 0.88, 0.85);
  MyEtaECTileLegend->SetFillStyle(1001);
  MyEtaECTileLegend->SetFillColor(10); 
  
  sprintf(QueDir,"%s %d", "HEC1:", h_Eta_HEC1_OnlyEC->Integral());
  MyEtaECTileLegend->AddEntry( h_Eta_HEC1_OnlyEC, QueDir, "l");
  sprintf(QueDir,"%s %d", "HEC2:", h_Eta_HEC2_OnlyEC->Integral());
  MyEtaECTileLegend->AddEntry( h_Eta_HEC2_OnlyEC, QueDir, "l");
  sprintf(QueDir,"%s %d", "HEC3:",  h_Eta_HEC3_OnlyEC->Integral());
  MyEtaECTileLegend->AddEntry( h_Eta_HEC3_OnlyEC, QueDir, "l");
  sprintf(QueDir,"%s %d", "HEC4:",  h_Eta_HEC4_OnlyEC->Integral());
  MyEtaECTileLegend->AddEntry( h_Eta_HEC4_OnlyEC, QueDir, "l");
  MyEtaECTileLegend->Draw();
  

  // //Phi parameter
  // TCanvas* Can_Phi = new TCanvas("Can_PhiBarrel","Can_PhiBarrel",900, 450);
  // Can_Phi->Divide(2,3);

  // Can_Phi->cd(1);
  // gPad->SetGridx(1); gPad->SetGridy(1);
  // TH1 *h_Phi_EMB1_OnlyBa =  (TH1F*)HistFile->Get("h_Phi_EMB1_OnlyBa");
  // h_Phi_EMB1_OnlyBa->SetLineColor(kBlue+2);
  // h_Phi_EMB1_OnlyBa->SetTitle("Barrel Tracks");
  // h_Phi_EMB1_OnlyBa->GetXaxis()->SetTitle("phi");
  // h_Phi_EMB1_OnlyBa->Draw("hist");
  // TH1 *h_Phi_EMB2_OnlyBa =  (TH1F*)HistFile->Get("h_Phi_EMB2_OnlyBa");
  // h_Phi_EMB2_OnlyBa->SetLineColor(kBlue-2);
  // h_Phi_EMB2_OnlyBa->Draw("histsame");
  // TH1 *h_Phi_EMB3_OnlyBa =  (TH1F*)HistFile->Get("h_Phi_EMB3_OnlyBa");
  // h_Phi_EMB3_OnlyBa->SetLineColor(kCyan-7);
  // h_Phi_EMB3_OnlyBa->Draw("histsame");
  // TH1 *h_Phi_Tile1_OnlyBa =  (TH1F*)HistFile->Get("h_Phi_Tile1_OnlyBa");
  // h_Phi_Tile1_OnlyBa->SetLineColor(kGreen+10);
  // h_Phi_Tile1_OnlyBa->Draw("histsame");
  // TH1 *h_Phi_Tile2_OnlyBa =  (TH1F*)HistFile->Get("h_Phi_Tile2_OnlyBa");
  // h_Phi_Tile2_OnlyBa->SetLineColor(kGreen+2);
  // h_Phi_Tile2_OnlyBa->Draw("histsame");
  // TH1 *h_Phi_Tile3_OnlyBa =  (TH1F*)HistFile->Get("h_Phi_Tile3_OnlyBa");
  // h_Phi_Tile3_OnlyBa->SetLineColor(kOrange+9);
  // h_Phi_Tile3_OnlyBa->Draw("histsame");

  
  // TH1 *h_EtaOnlyBarrel =  (TH1F*)HistFile->Get("h_EtaOnlyBarrel");
  // h_EtaOnlyBarrel->SetLineColor(kOrange);
  // h_EtaOnlyBarrel->SetLineWidth(2);
  // h_EtaOnlyBarrel->SetTitle("");
  // h_EtaOnlyBarrel->GetXaxis()->SetTitle("#eta");
  // h_EtaOnlyBarrel->Draw("hist");



  
  TCanvas* Can_TracExtrapolation = new TCanvas("TracExtrapolation","TracExtrapolation",900, 450);
  Can_TracExtrapolation->Divide(3,1);
  Can_TracExtrapolation->cd(1);
  
  TH1 *h_ExtrapolationQuality_EMB1 =  (TH1F*)HistFile->Get("h_ExtrapolationQuality_EMB1");
  h_ExtrapolationQuality_EMB1->SetLineColor(kOrange);
  h_ExtrapolationQuality_EMB1->SetLineWidth(2);
  h_ExtrapolationQuality_EMB1->SetTitle("Barrel (EMB1(Orange),EMB2(Green),EMB3(Blue))");
  gPad->SetGridx(1); gPad->SetGridy(1);
  // h_ExtrapolationQuality_EMB1->GetYaxis()->SetRangeUser(0.,800.);
  h_ExtrapolationQuality_EMB1->GetXaxis()->SetTitle("Track extrapolation");
  h_ExtrapolationQuality_EMB1->Draw("hist");

  TH1 *h_ExtrapolationQuality_EMB2 =  (TH1F*)HistFile->Get("h_ExtrapolationQuality_EMB2");
  h_ExtrapolationQuality_EMB2->SetLineColor(kGreen);
  h_ExtrapolationQuality_EMB2->SetLineWidth(2);
  h_ExtrapolationQuality_EMB2->Draw("histsame");
  
  TH1 *h_ExtrapolationQuality_EMB3 =  (TH1F*)HistFile->Get("h_ExtrapolationQuality_EMB3");
  h_ExtrapolationQuality_EMB3->SetLineColor(kBlue);
  h_ExtrapolationQuality_EMB3->SetLineWidth(2);
  h_ExtrapolationQuality_EMB3->Draw("histsame");

  Can_TracExtrapolation->cd(2);
  TH1 *h_ExtrapolationQuality_EME1 =  (TH1F*)HistFile->Get("h_ExtrapolationQuality_EME1");
  h_ExtrapolationQuality_EME1->SetLineColor(kOrange);
  h_ExtrapolationQuality_EME1->SetLineWidth(2);
  h_ExtrapolationQuality_EME1->SetTitle("End-cap (EME1(Orange),EME2(Green),EME3(Blue))");
  gPad->SetGridx(1); gPad->SetGridy(1);
3  //h_ExtrapolationQuality_EME1->GetYaxis()->SetRangeUser(0.,800.);
  h_ExtrapolationQuality_EME1->GetXaxis()->SetTitle("Track extrapolation");
  h_ExtrapolationQuality_EME1->Draw("hist");

  TH1 *h_ExtrapolationQuality_EME2 =  (TH1F*)HistFile->Get("h_ExtrapolationQuality_EME2");
  h_ExtrapolationQuality_EME2->SetLineColor(kGreen);
  h_ExtrapolationQuality_EME2->SetLineWidth(2);
  h_ExtrapolationQuality_EME2->Draw("histsame");
  
  TH1 *h_ExtrapolationQuality_EME3 =  (TH1F*)HistFile->Get("h_ExtrapolationQuality_EME3");
  h_ExtrapolationQuality_EME3->SetLineColor(kBlue);
  h_ExtrapolationQuality_EME3->SetLineWidth(2);
  h_ExtrapolationQuality_EME3->Draw("histsame");

  Can_TracExtrapolation->cd(3);
  TH1 *h_Eta_Tile1 =  (TH1F*)HistFile->Get("h_Eta_Tile1");
  h_Eta_Tile1->SetLineColor(kOrange+2);
  h_Eta_Tile1->SetLineWidth(2);
  h_Eta_Tile1->SetTitle("Tile1(Orange) & Tile3(Blue)");
  gPad->SetGridx(1); gPad->SetGridy(1);
  // h_ExtrapolationQuality_EME1->GetYaxis()->SetRangeUser(0.,800.);
  h_Eta_Tile1->GetXaxis()->SetTitle("#eta");
  h_Eta_Tile1->Draw("hist");

  TH1 *h_Eta_Tile3 =  (TH1F*)HistFile->Get("h_Eta_Tile3");
  h_Eta_Tile3->SetLineWidth(2);
  h_Eta_Tile3->SetLineColor(kBlue);
  h_Eta_Tile3->Draw("histsame");


  //Can Extrapolator Barrel vs EC (EM
  gStyle->SetOptStat(10);
  TCanvas* Can_TracExtrapolation_EM = new TCanvas("TracExtrapolation_EM","TracExtrapolation_EM",900, 450);
  Can_TracExtrapolation_EM->Divide(3,1);
  
  Can_TracExtrapolation_EM->cd(1);
  TH2 *h_ExtrapolationQuality_EM1 =  (TH2F*)HistFile->Get("h_ExtrapolationQuality_EM1");
  h_ExtrapolationQuality_EM1->GetZaxis()->SetRangeUser(0,h_ExtrapolationQuality_EM1->GetEntries());
  h_ExtrapolationQuality_EM1->GetXaxis()->SetTitle("#eta_{EMB1}");
  h_ExtrapolationQuality_EM1->GetYaxis()->SetTitle("#eta_{EME1}");
  h_ExtrapolationQuality_EM1->Draw("TEXTCOLZ");
  

  Can_TracExtrapolation_EM->cd(2);
  TH2 *h_ExtrapolationQuality_EM2 =  (TH2F*)HistFile->Get("h_ExtrapolationQuality_EM2");
  h_ExtrapolationQuality_EM2->GetZaxis()->SetRangeUser(0,h_ExtrapolationQuality_EM2->GetEntries());
  h_ExtrapolationQuality_EM2->GetXaxis()->SetTitle("#eta_{EMB2}");
  h_ExtrapolationQuality_EM2->GetYaxis()->SetTitle("#eta_{EME2}");
  h_ExtrapolationQuality_EM2->Draw("TEXTCOLZ");

  Can_TracExtrapolation_EM->cd(3);
  TH2 *h_ExtrapolationQuality_EM3 =  (TH2F*)HistFile->Get("h_ExtrapolationQuality_EM3");
  h_ExtrapolationQuality_EM3->GetZaxis()->SetRangeUser(0,h_ExtrapolationQuality_EM3->GetEntries());
  h_ExtrapolationQuality_EM3->GetXaxis()->SetTitle("#eta_{EMB3}");
  h_ExtrapolationQuality_EM3->GetYaxis()->SetTitle("#eta_{EME3}");
  h_ExtrapolationQuality_EM3->Draw("TEXTCOLZ");

  
  return;
}
  

//////////////////////////////
void Efficiency()
{
  
  TCanvas* Can_Efficiency = new TCanvas("Efficiency","Efficiency",900, 800);
  Can_Efficiency->Divide(3,2);

  Can_Efficiency->cd(1);
  TH1 *h_Eff_0_2GeV_eta1 =  (TH1F*)HistFile->Get("EfficiencyPlot_EM2_0_2GeV_eta1");
  h_Eff_0_2GeV_eta1->SetStats(kFALSE);
  h_Eff_0_2GeV_eta1->SetTitle("matched cluster efficiency in |#eta_{EM2}| < 1.0");
  h_Eff_0_2GeV_eta1 ->GetXaxis()->SetTitle("#varepsilon_{clu}");
  h_Eff_0_2GeV_eta1 ->GetYaxis()->SetTitle("Fraction of particles");
  h_Eff_0_2GeV_eta1->SetLineColor(kBlue);
  h_Eff_0_2GeV_eta1->Scale(1.0/h_Eff_0_2GeV_eta1->Integral());
  h_Eff_0_2GeV_eta1->Rebin(5);
  h_Eff_0_2GeV_eta1->GetYaxis()->SetRangeUser(0,1);
  h_Eff_0_2GeV_eta1->Draw();
  TH1 *h_Eff_2_5GeV_eta1 =  (TH1F*)HistFile->Get("EfficiencyPlot_EM2_2_5GeV_eta1");
  h_Eff_2_5GeV_eta1->SetLineColor(kRed);
  h_Eff_2_5GeV_eta1->Scale(1.0/h_Eff_2_5GeV_eta1->Integral());
  h_Eff_2_5GeV_eta1->Rebin(5);
  h_Eff_2_5GeV_eta1->Draw("same");
  TH1 *h_Eff_5GeV_eta1 =  (TH1F*)HistFile->Get("EfficiencyPlot_EM2_5GeV_eta1");
  h_Eff_5GeV_eta1->SetLineColor(kMagenta);
  h_Eff_5GeV_eta1->Scale(1.0/h_Eff_5GeV_eta1->Integral());
  h_Eff_5GeV_eta1->Rebin(5);
  h_Eff_5GeV_eta1->Draw("same");

  TLegend *MyEffLegend = new TLegend (0.30, 0.70, 0.6, 0.85);
  MyEffLegend->SetFillStyle(1001);
  MyEffLegend->SetFillColor(10);
  MyEffLegend->SetLineColor(0);
  MyEffLegend->AddEntry(  h_Eff_0_2GeV_eta1, "0< p_{track} < 2 GeV", "l");
  MyEffLegend->AddEntry(  h_Eff_2_5GeV_eta1, "2 < p_{track} < 5 GeV", "l");
  MyEffLegend->AddEntry(  h_Eff_5GeV_eta1, "p_{track} > 5 GeV", "l");
  MyEffLegend->Draw();

  
  Can_Efficiency->cd(2);
  TH1 *h_Eff_0_2GeV_eta2 =  (TH1F*)HistFile->Get("EfficiencyPlot_EM2_0_2GeV_eta2");
  h_Eff_0_2GeV_eta2->SetStats(kFALSE);
  h_Eff_0_2GeV_eta2->SetTitle("matched cluster efficiency in 1.0 < |#eta_{EM2}| < 2.0");
  h_Eff_0_2GeV_eta2 ->GetXaxis()->SetTitle("#varepsilon_{clu}");
  h_Eff_0_2GeV_eta2->SetLineColor(kBlue);
  h_Eff_0_2GeV_eta2->Scale(1.0/h_Eff_0_2GeV_eta2->Integral());
  h_Eff_0_2GeV_eta2->Rebin(5);
  h_Eff_0_2GeV_eta2->GetYaxis()->SetRangeUser(0,1);
  h_Eff_0_2GeV_eta2->Draw();
  TH1 *h_Eff_2_5GeV_eta2 =  (TH1F*)HistFile->Get("EfficiencyPlot_EM2_2_5GeV_eta2");
  h_Eff_2_5GeV_eta2->SetLineColor(kRed);
  h_Eff_2_5GeV_eta2->Scale(1.0/h_Eff_2_5GeV_eta2->Integral());
  h_Eff_2_5GeV_eta2->Rebin(5);
  h_Eff_2_5GeV_eta2->Draw("same");
  TH1 *h_Eff_5GeV_eta2 =  (TH1F*)HistFile->Get("EfficiencyPlot_EM2_5GeV_eta2");
  h_Eff_5GeV_eta2->SetLineColor(kPink);
  h_Eff_5GeV_eta2->Scale(1.0/h_Eff_5GeV_eta2->Integral());
  h_Eff_5GeV_eta2->Rebin(5);
  h_Eff_5GeV_eta2->Draw("same");
  MyEffLegend->Draw();


  Can_Efficiency->cd(3);
  TH1 *h_Eff_0_2GeV_eta25 =  (TH1F*)HistFile->Get("EfficiencyPlot_EM2_0_2GeV_eta25");
  h_Eff_0_2GeV_eta25->SetStats(kFALSE);
  h_Eff_0_2GeV_eta25->SetTitle("matched cluster efficiency in 2.0 < |#eta_{EM2}| < 2.5");
  h_Eff_0_2GeV_eta25 ->GetXaxis()->SetTitle("#varepsilon_{clu}");
  h_Eff_0_2GeV_eta25->SetLineColor(kBlue);
  h_Eff_0_2GeV_eta25->Scale(1.0/h_Eff_0_2GeV_eta25->Integral());
  h_Eff_0_2GeV_eta25->Rebin(5);
  h_Eff_0_2GeV_eta25->GetYaxis()->SetRangeUser(0,1);
  h_Eff_0_2GeV_eta25->Draw();
  TH1 *h_Eff_2_5GeV_eta25 =  (TH1F*)HistFile->Get("EfficiencyPlot_EM2_2_5GeV_eta25");
  h_Eff_2_5GeV_eta25->SetLineColor(kRed);
  h_Eff_2_5GeV_eta25->Scale(1.0/h_Eff_2_5GeV_eta25->Integral());
  h_Eff_2_5GeV_eta25->Rebin(5);
  h_Eff_2_5GeV_eta25->Draw("same");
  TH1 *h_Eff_5GeV_eta25 =  (TH1F*)HistFile->Get("EfficiencyPlot_EM2_5GeV_eta25");
  h_Eff_5GeV_eta25->SetLineColor(kPink);
  h_Eff_5GeV_eta25->Scale(1.0/h_Eff_5GeV_eta25->Integral());
  h_Eff_5GeV_eta25->Rebin(5);
  h_Eff_5GeV_eta25->Draw("same");
  MyEffLegend->Draw();
  
  return;
}


//////////////////////////////
void Purity()
{
  
  TCanvas* Can_Purity = new TCanvas("Purity","Purity",900, 800);
  Can_Purity->Divide(3,2);

  Can_Purity->cd(1);
  TH1 *h_Pur_0_2GeV_eta1 =  (TH1F*)HistFile->Get("PurityPlot_EM2_0_2GeV_eta1");
  h_Pur_0_2GeV_eta1->SetStats(kFALSE);
  h_Pur_0_2GeV_eta1->SetTitle("matched cluster purity in |#eta_{EM2}| < 1.0");
  h_Pur_0_2GeV_eta1 ->GetXaxis()->SetTitle("#rho_{clu}");
  h_Pur_0_2GeV_eta1 ->GetYaxis()->SetTitle("Fraction of particles");
  h_Pur_0_2GeV_eta1->SetLineColor(kBlue);
  h_Pur_0_2GeV_eta1->Scale(1.0/h_Pur_0_2GeV_eta1->Integral());
  h_Pur_0_2GeV_eta1->Rebin(5);
  h_Pur_0_2GeV_eta1->GetYaxis()->SetRangeUser(0,1);
  h_Pur_0_2GeV_eta1->Draw();
  TH1 *h_Pur_2_5GeV_eta1 =  (TH1F*)HistFile->Get("PurityPlot_EM2_2_5GeV_eta1");
  h_Pur_2_5GeV_eta1->SetLineColor(kRed);
  h_Pur_2_5GeV_eta1->Scale(1.0/h_Pur_2_5GeV_eta1->Integral());
  h_Pur_2_5GeV_eta1->Rebin(5);
  h_Pur_2_5GeV_eta1->Draw("same");
  TH1 *h_Pur_5GeV_eta1 =  (TH1F*)HistFile->Get("PurityPlot_EM2_5GeV_eta1");
  h_Pur_5GeV_eta1->SetLineColor(kGreen+3);
  h_Pur_5GeV_eta1->Scale(1.0/h_Pur_5GeV_eta1->Integral());
  h_Pur_5GeV_eta1->Rebin(5);
  h_Pur_5GeV_eta1->Draw("same");

  TLegend *MyPurLegend = new TLegend (0.30, 0.70, 0.6, 0.85);
  MyPurLegend->SetFillStyle(1001);
  MyPurLegend->SetFillColor(10);
  MyPurLegend->SetLineColor(0);
  MyPurLegend->AddEntry(  h_Pur_0_2GeV_eta1, "p_{track} < 2 GeV", "l");
  MyPurLegend->AddEntry(  h_Pur_2_5GeV_eta1, "2 < p_{track} < 5 GeV", "l");
  MyPurLegend->AddEntry(  h_Pur_5GeV_eta1,  "p_{track} > 5 GeV", "l");
  MyPurLegend->Draw();

  
  Can_Purity->cd(2);
  TH1 *h_Pur_0_2GeV_eta2 =  (TH1F*)HistFile->Get("PurityPlot_EM2_0_2GeV_eta2");
  h_Pur_0_2GeV_eta2->SetStats(kFALSE);
  h_Pur_0_2GeV_eta2->SetTitle("matched cluster purity in 1.0 < |#eta_{EM2}| < 2.0");
  h_Pur_0_2GeV_eta2 ->GetXaxis()->SetTitle("#rho_{clu}");
  h_Pur_0_2GeV_eta2->SetLineColor(kBlue);
  h_Pur_0_2GeV_eta2->Scale(1.0/h_Pur_0_2GeV_eta2->Integral());
  h_Pur_0_2GeV_eta2->Rebin(5);
  h_Pur_0_2GeV_eta2->GetYaxis()->SetRangeUser(0,1);
  h_Pur_0_2GeV_eta2->Draw();
  TH1 *h_Pur_2_5GeV_eta2 =  (TH1F*)HistFile->Get("PurityPlot_EM2_2_5GeV_eta2");
  h_Pur_2_5GeV_eta2->SetLineColor(kRed);
  h_Pur_2_5GeV_eta2->Scale(1.0/h_Pur_2_5GeV_eta2->Integral());
  h_Pur_2_5GeV_eta2->Rebin(5);
  h_Pur_2_5GeV_eta2->Draw("same");
  TH1 *h_Pur_5GeV_eta2 =  (TH1F*)HistFile->Get("PurityPlot_EM2_5GeV_eta2");
  h_Pur_5GeV_eta2->SetLineColor(kGreen+3);
  h_Pur_5GeV_eta2->Scale(1.0/h_Pur_5GeV_eta2->Integral());
  h_Pur_5GeV_eta2->Rebin(5);
  h_Pur_5GeV_eta2->Draw("same");
  MyPurLegend->Draw();

  Can_Purity->cd(3);
  TH1 *h_Pur_0_2GeV_eta25 =  (TH1F*)HistFile->Get("PurityPlot_EM2_0_2GeV_eta25");
  h_Pur_0_2GeV_eta25->SetStats(kFALSE);
  h_Pur_0_2GeV_eta25->SetTitle("matched cluster purity in 2.0 < |#eta_{EM2}| < 2.5");
  h_Pur_0_2GeV_eta25 ->GetXaxis()->SetTitle("#rho_{clu}");
  h_Pur_0_2GeV_eta25->SetLineColor(kBlue);
  h_Pur_0_2GeV_eta25->Rebin(5);
  h_Pur_0_2GeV_eta25->Scale(1.0/h_Pur_0_2GeV_eta25->Integral());
  h_Pur_0_2GeV_eta25->GetYaxis()->SetRangeUser(0,1);
  h_Pur_0_2GeV_eta25->Draw();
  TH1 *h_Pur_2_5GeV_eta25 =  (TH1F*)HistFile->Get("PurityPlot_EM2_2_5GeV_eta25");
  h_Pur_2_5GeV_eta25->SetLineColor(kRed);
  h_Pur_2_5GeV_eta25->Rebin(5);
  h_Pur_2_5GeV_eta25->Scale(1.0/h_Pur_2_5GeV_eta25->Integral());
  h_Pur_2_5GeV_eta25->Draw("same");
  TH1 *h_Pur_5GeV_eta25 =  (TH1F*)HistFile->Get("PurityPlot_EM2_5GeV_eta25");
  h_Pur_5GeV_eta25->SetLineColor(kGreen+3);
  h_Pur_5GeV_eta25->Rebin(5);
  h_Pur_5GeV_eta25->Scale(1.0/h_Pur_5GeV_eta25->Integral());
  h_Pur_5GeV_eta25->Draw("same");
  MyPurLegend->Draw();
  
  return;
}


//////////////////////////////
void NCluster09()
{
  TCanvas* Can_NCluster09 = new TCanvas("NCluster09","NCluster09",900, 800);
  Can_NCluster09->Divide(3,2);

  Can_NCluster09->cd(1);
  TH1 *h_nClus09_0_2GeV_eta1 = (TH1F*)HistFile->Get("NCluster_09_0_2GeV_eta1");
  h_nClus09_0_2GeV_eta1->SetStats(kFALSE);
  h_nClus09_0_2GeV_eta1->SetTitle("|#eta_{EM2}| < 1.0");
  h_nClus09_0_2GeV_eta1 ->GetXaxis()->SetTitle("N_{clu}(#Sigma E^{true}>90%)");
  h_nClus09_0_2GeV_eta1 ->GetYaxis()->SetTitle("Fraction of particles");
  h_nClus09_0_2GeV_eta1 ->SetLineColor(kBlue);
  h_nClus09_0_2GeV_eta1->Scale(1.0/h_nClus09_0_2GeV_eta1->Integral());
  h_nClus09_0_2GeV_eta1->GetYaxis()->SetRangeUser(0,1);
  h_nClus09_0_2GeV_eta1->Draw();
  TH1 *h_nClus09_2_5GeV_eta1 =  (TH1F*)HistFile->Get("NCluster_09_2_5GeV_eta1");
  h_nClus09_2_5GeV_eta1->SetLineColor(kRed);
  h_nClus09_2_5GeV_eta1->Scale(1.0/h_nClus09_2_5GeV_eta1->Integral());
  h_nClus09_2_5GeV_eta1->Draw("same");
  TH1 *h_nClus09_5GeV_eta1 =  (TH1F*)HistFile->Get("NCluster_09_5GeV_eta1");
  h_nClus09_5GeV_eta1->SetLineColor(kGreen+3);
  h_nClus09_5GeV_eta1->Scale(1.0/h_nClus09_5GeV_eta1->Integral());
  h_nClus09_5GeV_eta1->Draw("same");

  TLegend *MyClus09 = new TLegend (0.30, 0.70, 0.6, 0.85);
  MyClus09->SetFillStyle(1001);
  MyClus09->SetFillColor(10);
  MyClus09->SetLineColor(0);
  MyClus09->AddEntry(  h_nClus09_0_2GeV_eta1, "p^{track} < 2 GeV", "l");
  MyClus09->AddEntry(  h_nClus09_2_5GeV_eta1, "2 < p^{track} < 5 GeV", "l");
  MyClus09->AddEntry(  h_nClus09_5GeV_eta1,  "p^{track} > 5 GeV", "l");
  MyClus09->Draw();

  Can_NCluster09->cd(2);
  TH1 *h_nClus09_0_2GeV_eta2 =  (TH1F*)HistFile->Get("NCluster_09_0_2GeV_eta2");
  h_nClus09_0_2GeV_eta2->SetStats(kFALSE);
  h_nClus09_0_2GeV_eta2->SetTitle("1.0 <|#eta_{EM2}| < 2.0");
  h_nClus09_0_2GeV_eta2 ->GetXaxis()->SetTitle("N_{clu}(#Sigma E^{true}>90%)");
  h_nClus09_0_2GeV_eta2->SetLineColor(kBlue);
  h_nClus09_0_2GeV_eta2->Scale(1.0/h_nClus09_0_2GeV_eta2->Integral());
  h_nClus09_0_2GeV_eta2->GetYaxis()->SetRangeUser(0,1);
  h_nClus09_0_2GeV_eta2->Draw();
  TH1 *h_nClus09_2_5GeV_eta2 =  (TH1F*)HistFile->Get("NCluster_09_2_5GeV_eta2");
  h_nClus09_2_5GeV_eta2->SetLineColor(kRed);
  h_nClus09_2_5GeV_eta2->Scale(1.0/h_nClus09_2_5GeV_eta2->Integral());
  h_nClus09_2_5GeV_eta2->Draw("same");
  TH1 *h_nClus09_5GeV_eta2 =  (TH1F*)HistFile->Get("NCluster_09_5GeV_eta2");
  h_nClus09_5GeV_eta2->SetLineColor(kGreen+3);
  h_nClus09_5GeV_eta2->Scale(1.0/h_nClus09_5GeV_eta2->Integral());
  h_nClus09_5GeV_eta2->Draw("same");
  MyClus09->Draw();

  Can_NCluster09->cd(3);
  TH1 *h_nClus09_0_2GeV_eta25 =  (TH1F*)HistFile->Get("NCluster_09_0_2GeV_eta25");
  h_nClus09_0_2GeV_eta25->SetStats(kFALSE);
  h_nClus09_0_2GeV_eta25->SetTitle("|#eta_{EM2}| < 2.5");
  h_nClus09_0_2GeV_eta25 ->GetXaxis()->SetTitle("N_{clu}(#Sigma E^{true}>90%)");
  h_nClus09_0_2GeV_eta25->SetLineColor(kBlue);
  h_nClus09_0_2GeV_eta25->Scale(1.0/h_nClus09_0_2GeV_eta25->Integral());
  h_nClus09_0_2GeV_eta25->GetYaxis()->SetRangeUser(0,1);
  h_nClus09_0_2GeV_eta25->Draw();
  TH1 *h_nClus09_2_5GeV_eta25 =  (TH1F*)HistFile->Get("NCluster_09_2_5GeV_eta25");
  h_nClus09_2_5GeV_eta25->SetLineColor(kRed);
  h_nClus09_2_5GeV_eta25->Scale(1.0/h_nClus09_2_5GeV_eta25->Integral());
  h_nClus09_2_5GeV_eta25->Draw("same");
  TH1 *h_nClus09_5GeV_eta25 =  (TH1F*)HistFile->Get("NCluster_09_5GeV_eta25");
  h_nClus09_5GeV_eta25->SetLineColor(kGreen+3);
  h_nClus09_5GeV_eta25->Scale(1.0/h_nClus09_5GeV_eta25->Integral());
  h_nClus09_5GeV_eta25->Draw("same");
  MyClus09->Draw();
  
  return;
}

//////////////////////////////
void DeltaR_EM2()
{
  TCanvas* Can_DeltaR = new TCanvas("DeltaR","DeltaR",900, 800);
  Can_DeltaR->Divide(3,2);

  Can_DeltaR->cd(1);
  
  TH1 *h_DeltaREM2_0_2GeV = (TH1F*)HistFile->Get("dR_EM2_0_2GeV");
  h_DeltaREM2_0_2GeV ->SetStats(kFALSE);
  h_DeltaREM2_0_2GeV ->SetTitle("#DeltaR");
  h_DeltaREM2_0_2GeV ->GetXaxis()->SetTitle("DeltaR");
  h_DeltaREM2_0_2GeV ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaREM2_0_2GeV ->SetLineColor(kBlue);
  h_DeltaREM2_0_2GeV->Scale(1.0/h_DeltaREM2_0_2GeV->Integral());
  h_DeltaREM2_0_2GeV->GetYaxis()->SetRangeUser(0,0.6);
  h_DeltaREM2_0_2GeV->Draw();
  TH1 *h_DeltaREM2_2_5GeV =  (TH1F*)HistFile->Get("dR_EM2_2_5GeV");
  h_DeltaREM2_2_5GeV->SetLineColor(kRed);
  h_DeltaREM2_2_5GeV->Scale(1.0/h_DeltaREM2_2_5GeV->Integral());
  h_DeltaREM2_2_5GeV->Draw("same");
  TH1 *h_DeltaREM2_5GeV =  (TH1F*)HistFile->Get("dR_EM2_5GeV");
  h_DeltaREM2_5GeV->SetLineColor(kGreen+3);
  h_DeltaREM2_5GeV->Scale(1.0/h_DeltaREM2_5GeV->Integral());
  h_DeltaREM2_5GeV->Draw("same");

  TLegend *MyClus09 = new TLegend (0.30, 0.70, 0.6, 0.85);
  MyClus09->SetFillStyle(1001);
  MyClus09->SetFillColor(10);
  MyClus09->SetLineColor(0);
  MyClus09->AddEntry(  h_DeltaREM2_0_2GeV, "p^{track} < 2 GeV", "l");
  MyClus09->AddEntry(  h_DeltaREM2_2_5GeV, "2 < p^{track} < 5 GeV", "l");
  MyClus09->AddEntry(  h_DeltaREM2_5GeV,  "p^{track} > 5 GeV", "l");
  MyClus09->Draw();

 
  return;
}


//////////////////////////////
void EfficiencyLead()
{
  
  TCanvas* Can_EfficiencyLead = new TCanvas("EfficiencyLead","EfficiencyLead",900, 800);
  Can_EfficiencyLead->Divide(3,2);

  Can_EfficiencyLead->cd(1);
  TH1 *h_EffLead_1_2GeV_eta1 =  (TH1F*)HistFile->Get("Eff_Lead_1_2GeV_eta1");
  h_EffLead_1_2GeV_eta1->SetStats(kFALSE);
  h_EffLead_1_2GeV_eta1->SetTitle("matched cluster efficiency in |#eta_{EM2}| < 1.0");
  h_EffLead_1_2GeV_eta1 ->GetXaxis()->SetTitle("#varepsilon_{clu}^{lead}");
  h_EffLead_1_2GeV_eta1 ->GetYaxis()->SetTitle("Fraction of particles");
  h_EffLead_1_2GeV_eta1->SetLineColor(kBlue);
  h_EffLead_1_2GeV_eta1 ->SetFillColor(kBlue);h_EffLead_1_2GeV_eta1 ->SetFillStyle(3003);
  if(h_EffLead_1_2GeV_eta1)h_EffLead_1_2GeV_eta1->Scale(1.0/h_EffLead_1_2GeV_eta1->Integral());
  h_EffLead_1_2GeV_eta1->Rebin(2);
  h_EffLead_1_2GeV_eta1->GetYaxis()->SetRangeUser(0,0.7);
  h_EffLead_1_2GeV_eta1->Draw();
  TH1 *h_EffLead_2_5GeV_eta1 =  (TH1F*)HistFile->Get("Eff_Lead_2_5GeV_eta1");
  h_EffLead_2_5GeV_eta1->SetLineColor(kRed);
  h_EffLead_2_5GeV_eta1 ->SetFillColor(kRed);h_EffLead_2_5GeV_eta1 ->SetFillStyle(3003);
  if(h_EffLead_2_5GeV_eta1)h_EffLead_2_5GeV_eta1->Scale(1.0/h_EffLead_2_5GeV_eta1->Integral());
  h_EffLead_2_5GeV_eta1->Rebin(2);
  h_EffLead_2_5GeV_eta1->Draw("same");
  TH1 *h_EffLead_5_10GeV_eta1 =  (TH1F*)HistFile->Get("Eff_Lead_5_10GeV_eta1");
  h_EffLead_5_10GeV_eta1->SetLineColor(kGreen+3);
  h_EffLead_5_10GeV_eta1 ->SetFillColor(kGreen);h_EffLead_5_10GeV_eta1 ->SetFillStyle(3003);
  if(h_EffLead_5_10GeV_eta1)h_EffLead_5_10GeV_eta1->Scale(1.0/h_EffLead_5_10GeV_eta1->Integral());
  h_EffLead_5_10GeV_eta1->Rebin(2);
  h_EffLead_5_10GeV_eta1->Draw("same");
  TH1 *h_EffLead_10GeV_eta1 =  (TH1F*)HistFile->Get("Eff_Lead_10GeV_eta1");
  h_EffLead_10GeV_eta1->SetLineColor(kMagenta);
  h_EffLead_10GeV_eta1 ->SetFillColor(kMagenta);h_EffLead_10GeV_eta1 ->SetFillStyle(3003);
  if(h_EffLead_10GeV_eta1)h_EffLead_10GeV_eta1->Scale(1.0/h_EffLead_10GeV_eta1->Integral());
  h_EffLead_10GeV_eta1->Rebin(2);
  h_EffLead_10GeV_eta1->Draw("same");
  
  TPaveText *EtaRange = new TPaveText(0.20, 0.57, 0.40, 0.67,"NDC");
  EtaRange->AddText("|#eta| < 1.0");
  EtaRange->SetFillColor(0);
  EtaRange->SetLineColor(0);
  EtaRange->SetShadowColor(0);
  EtaRange->Draw();

  TLegend *MyEffLegend = new TLegend (0.15, 0.70, 0.45, 0.85);
  MyEffLegend->SetFillStyle(1001);
  MyEffLegend->SetFillColor(10);
  MyEffLegend->SetLineColor(0);
  MyEffLegend->AddEntry(  h_EffLead_1_2GeV_eta1, "1 < p_{T}^{track} < 2 GeV", "l");
  MyEffLegend->AddEntry(  h_EffLead_2_5GeV_eta1, "2 < p_{T}^{track} < 5 GeV", "l");
  MyEffLegend->AddEntry(  h_EffLead_5_10GeV_eta1,"5 < p_{T}^{track} < 10 GeV", "l");
  MyEffLegend->AddEntry(  h_EffLead_10GeV_eta1,  "10 < p_{T}^{track} < 20 GeV", "l");
  MyEffLegend->Draw();

  
  Can_EfficiencyLead->cd(2);
  TH1 *h_EffLead_1_2GeV_eta2 =  (TH1F*)HistFile->Get("Eff_Lead_1_2GeV_eta2");
  h_EffLead_1_2GeV_eta2->SetStats(kFALSE);
  h_EffLead_1_2GeV_eta2->SetTitle("matched cluster efficiency in 1.0 < |#eta_{EM2}| < 2.0");
  h_EffLead_1_2GeV_eta2 ->GetXaxis()->SetTitle("#varepsilon_{clu}^{lead}");
  h_EffLead_1_2GeV_eta2->SetLineColor(kBlue);
  h_EffLead_1_2GeV_eta2 ->SetFillColor(kBlue);h_EffLead_1_2GeV_eta2 ->SetFillStyle(3003);
  if(h_EffLead_1_2GeV_eta2) h_EffLead_1_2GeV_eta2->Scale(1.0/h_EffLead_1_2GeV_eta2->Integral());
  h_EffLead_1_2GeV_eta2->Rebin(2);
  h_EffLead_1_2GeV_eta2->GetYaxis()->SetRangeUser(0,0.7);
  h_EffLead_1_2GeV_eta2->Draw();
  TH1 *h_EffLead_2_5GeV_eta2 =  (TH1F*)HistFile->Get("Eff_Lead_2_5GeV_eta2");
  h_EffLead_2_5GeV_eta2->SetLineColor(kRed);
  h_EffLead_2_5GeV_eta2 ->SetFillColor(kRed);h_EffLead_2_5GeV_eta2 ->SetFillStyle(3003);
  if( h_EffLead_2_5GeV_eta2)h_EffLead_2_5GeV_eta2->Scale(1.0/h_EffLead_2_5GeV_eta2->Integral());
  h_EffLead_2_5GeV_eta2->Rebin(2);
  h_EffLead_2_5GeV_eta2->Draw("same");
  TH1 *h_EffLead_5_10GeV_eta2 =  (TH1F*)HistFile->Get("Eff_Lead_5_10GeV_eta2");
  h_EffLead_5_10GeV_eta2->SetLineColor(kGreen+3);
  h_EffLead_5_10GeV_eta2 ->SetFillColor(kGreen);h_EffLead_5_10GeV_eta2 ->SetFillStyle(3003);
  if(h_EffLead_5_10GeV_eta2)h_EffLead_5_10GeV_eta2->Scale(1.0/h_EffLead_5_10GeV_eta2->Integral());
  h_EffLead_5_10GeV_eta2->Rebin(2);
  h_EffLead_5_10GeV_eta2->Draw("same");
  TH1 *h_EffLead_10GeV_eta2 =  (TH1F*)HistFile->Get("Eff_Lead_10GeV_eta2");
  h_EffLead_10GeV_eta2->SetLineColor(kMagenta);
  h_EffLead_10GeV_eta2 ->SetFillColor(kMagenta);h_EffLead_10GeV_eta2 ->SetFillStyle(3003);
  if(h_EffLead_10GeV_eta2)h_EffLead_10GeV_eta2->Scale(1.0/h_EffLead_10GeV_eta2->Integral());
  h_EffLead_10GeV_eta2->Rebin(2);
  h_EffLead_10GeV_eta2->Draw("same");
  MyEffLegend->Draw();
  TPaveText *EtaRange2 = new TPaveText(0.15, 0.57, 0.45, 0.67,"NDC");
  EtaRange2->AddText("1.0 <|#eta| < 2.0");
  EtaRange2->SetFillColor(0);
  EtaRange2->SetLineColor(0);
  EtaRange2->SetShadowColor(0);
  EtaRange2->Draw();
  
  Can_EfficiencyLead->cd(3);
  TH1 *h_EffLead_1_2GeV_eta25 =  (TH1F*)HistFile->Get("Eff_Lead_1_2GeV_eta25");
  h_EffLead_1_2GeV_eta25->SetStats(kFALSE);
  h_EffLead_1_2GeV_eta25->SetTitle("matched cluster efficiency in 2.0 < |#eta_{EM2}| < 2.5");
  h_EffLead_1_2GeV_eta25 ->GetXaxis()->SetTitle("#varepsilon_{clu}^{lead}");
  h_EffLead_1_2GeV_eta25->SetLineColor(kBlue);
  h_EffLead_1_2GeV_eta25 ->SetFillColor(kBlue);h_EffLead_1_2GeV_eta25 ->SetFillStyle(3003);
  if(h_EffLead_1_2GeV_eta25)h_EffLead_1_2GeV_eta25->Scale(1.0/h_EffLead_1_2GeV_eta25->Integral());
  h_EffLead_1_2GeV_eta25->Rebin(2);
  h_EffLead_1_2GeV_eta25->GetYaxis()->SetRangeUser(0,0.7);
  h_EffLead_1_2GeV_eta25->Draw();
  TH1 *h_EffLead_2_5GeV_eta25 =  (TH1F*)HistFile->Get("Eff_Lead_2_5GeV_eta25");
  h_EffLead_2_5GeV_eta25->SetLineColor(kRed);
  h_EffLead_2_5GeV_eta25 ->SetFillColor(kRed);h_EffLead_2_5GeV_eta25 ->SetFillStyle(3003);
  if(h_EffLead_2_5GeV_eta25)h_EffLead_2_5GeV_eta25->Scale(1.0/h_EffLead_2_5GeV_eta25->Integral());
  h_EffLead_2_5GeV_eta25->Rebin(2);
  h_EffLead_2_5GeV_eta25->Draw("same");
  TH1 *h_EffLead_5_10GeV_eta25 =  (TH1F*)HistFile->Get("Eff_Lead_5_10GeV_eta25");
  h_EffLead_5_10GeV_eta25->SetLineColor(kGreen+3);
  h_EffLead_5_10GeV_eta25 ->SetFillColor(kGreen);h_EffLead_5_10GeV_eta25 ->SetFillStyle(3003);
  if(h_EffLead_5_10GeV_eta25)h_EffLead_5_10GeV_eta25->Scale(1.0/h_EffLead_5_10GeV_eta25->Integral());
  h_EffLead_5_10GeV_eta25->Rebin(2);
  h_EffLead_5_10GeV_eta25->Draw("same");
  TH1 *h_EffLead_10GeV_eta25 =  (TH1F*)HistFile->Get("Eff_Lead_10GeV_eta25");
  h_EffLead_10GeV_eta25->SetLineColor(kMagenta);
  h_EffLead_10GeV_eta25 ->SetFillColor(kMagenta);h_EffLead_10GeV_eta25 ->SetFillStyle(3003);
  if(h_EffLead_10GeV_eta25)h_EffLead_10GeV_eta25->Scale(1.0/h_EffLead_10GeV_eta25->Integral());
  h_EffLead_10GeV_eta25->Rebin(2);
  h_EffLead_10GeV_eta25->Draw("same");
  MyEffLegend->Draw();
  TPaveText *EtaRange3 = new TPaveText(0.50, 0.57, 0.80, 0.67,"NDC");
  EtaRange3->AddText("2.0 < |#eta| < 2.5");
  EtaRange3->SetFillColor(0);
  EtaRange3->SetLineColor(0);
  EtaRange3->SetShadowColor(0);
  EtaRange3->Draw();
  
  return;
}


//////////////////////////////
void PurityLead()
{
  
  TCanvas* Can_PurityLead = new TCanvas("PurityLead","PurityLead",900, 800);
  Can_PurityLead->Divide(3,2);

  Can_PurityLead->cd(1);
  TH1 *h_PurLead_1_2GeV_eta1 =  (TH1F*)HistFile->Get("Pur_Lead_1_2GeV_eta1");
  h_PurLead_1_2GeV_eta1->SetStats(kFALSE);
  h_PurLead_1_2GeV_eta1->SetTitle("matched cluster purity in |#eta_{EM2}| < 1.0");
  h_PurLead_1_2GeV_eta1 ->GetXaxis()->SetTitle("#rho_{clu}");
  h_PurLead_1_2GeV_eta1 ->GetYaxis()->SetTitle("Fraction of particles");
  h_PurLead_1_2GeV_eta1->SetLineColor(kBlue);
  h_PurLead_1_2GeV_eta1->Scale(1.0/h_PurLead_1_2GeV_eta1->Integral());
  h_PurLead_1_2GeV_eta1->Rebin(2);
  h_PurLead_1_2GeV_eta1->GetYaxis()->SetRangeUser(0,1.2);
  h_PurLead_1_2GeV_eta1->Draw();
  TH1 *h_PurLead_2_5GeV_eta1 =  (TH1F*)HistFile->Get("Pur_Lead_2_5GeV_eta1");
  h_PurLead_2_5GeV_eta1->SetLineColor(kRed);
  h_PurLead_2_5GeV_eta1->Scale(1.0/h_PurLead_2_5GeV_eta1->Integral());
  h_PurLead_2_5GeV_eta1->Rebin(2);
  h_PurLead_2_5GeV_eta1->Draw("same");
  TH1 *h_PurLead_5_10GeV_eta1 =  (TH1F*)HistFile->Get("Pur_Lead_5_10GeV_eta1");
  h_PurLead_5_10GeV_eta1->SetLineColor(kGreen+3);
  h_PurLead_5_10GeV_eta1->Scale(1.0/h_PurLead_5_10GeV_eta1->Integral());
  h_PurLead_5_10GeV_eta1->Rebin(2);
  h_PurLead_5_10GeV_eta1->Draw("same");
  TH1 *h_PurLead_10GeV_eta1 =  (TH1F*)HistFile->Get("Pur_Lead_10GeV_eta1");
  h_PurLead_10GeV_eta1->SetLineColor(kMagenta);
  h_PurLead_10GeV_eta1->Scale(1.0/h_PurLead_10GeV_eta1->Integral());
  h_PurLead_10GeV_eta1->Rebin(2);
  h_PurLead_10GeV_eta1->Draw("same");
  
  TPaveText *EtaRange = new TPaveText(0.20, 0.57, 0.40, 0.67,"NDC");
  EtaRange->AddText("|#eta| < 1.0");
  EtaRange->SetFillColor(0);
  EtaRange->SetLineColor(0);
  EtaRange->SetShadowColor(0);
  EtaRange->Draw();
  
  TLegend *MyPurLegend = new TLegend (0.15, 0.70, 0.45, 0.85);
  MyPurLegend->SetFillStyle(1001);
  MyPurLegend->SetFillColor(10);
  MyPurLegend->SetLineColor(0);
  MyPurLegend->AddEntry(  h_PurLead_1_2GeV_eta1, "1 < p_{T}^{track} < 2 GeV", "l");
  MyPurLegend->AddEntry(  h_PurLead_2_5GeV_eta1, "2 < p_{T}^{track} < 5 GeV", "l");
  MyPurLegend->AddEntry(  h_PurLead_5_10GeV_eta1,  "5 < p_{T}^{track} < 10 GeV", "l");
  MyPurLegend->AddEntry(  h_PurLead_10GeV_eta1,  " 10 < p_{T}^{track} < 20 GeV", "l");
  MyPurLegend->Draw();
  
  Can_PurityLead->cd(2);
  TH1 *h_PurLead_1_2GeV_eta2 =  (TH1F*)HistFile->Get("Pur_Lead_1_2GeV_eta2");
  h_PurLead_1_2GeV_eta2->SetStats(kFALSE);
  h_PurLead_1_2GeV_eta2->SetTitle("matched cluster purity in 1.0 < |#eta_{EM2}| < 2.0");
  h_PurLead_1_2GeV_eta2 ->GetXaxis()->SetTitle("#rho_{clu}");
  h_PurLead_1_2GeV_eta2->SetLineColor(kBlue);
  h_PurLead_1_2GeV_eta2->Scale(1.0/h_PurLead_1_2GeV_eta2->Integral());
  h_PurLead_1_2GeV_eta2->Rebin(2);
  h_PurLead_1_2GeV_eta2->GetYaxis()->SetRangeUser(0,1.2);
  h_PurLead_1_2GeV_eta2->Draw();
  TH1 *h_PurLead_2_5GeV_eta2 =  (TH1F*)HistFile->Get("Pur_Lead_2_5GeV_eta2");
  h_PurLead_2_5GeV_eta2->SetLineColor(kRed);
  h_PurLead_2_5GeV_eta2->Scale(1.0/h_PurLead_2_5GeV_eta2->Integral());
  h_PurLead_2_5GeV_eta2->Rebin(2);
  h_PurLead_2_5GeV_eta2->Draw("same");
  TH1 *h_PurLead_5_10GeV_eta2 =  (TH1F*)HistFile->Get("Pur_Lead_5_10GeV_eta2");
  h_PurLead_5_10GeV_eta2->SetLineColor(kGreen+3);
  h_PurLead_5_10GeV_eta2->Scale(1.0/h_PurLead_5_10GeV_eta2->Integral());
  h_PurLead_5_10GeV_eta2->Rebin(2);
  h_PurLead_5_10GeV_eta2->Draw("same");
  TH1 *h_PurLead_10GeV_eta2 =  (TH1F*)HistFile->Get("Pur_Lead_10GeV_eta2");
  h_PurLead_10GeV_eta2->SetLineColor(kMagenta);
  h_PurLead_10GeV_eta2->Scale(1.0/h_PurLead_10GeV_eta2->Integral());
  h_PurLead_10GeV_eta2->Rebin(2);
  h_PurLead_10GeV_eta2->Draw("same");
  MyPurLegend->Draw();
  TPaveText *EtaRange2 = new TPaveText(0.15, 0.57, 0.45, 0.67,"NDC");
  EtaRange2->AddText("1.0 <|#eta| < 2.0");
  EtaRange2->SetFillColor(0);
  EtaRange2->SetLineColor(0);
  EtaRange2->SetShadowColor(0);
  EtaRange2->Draw();

  
  Can_PurityLead->cd(3);
  TH1 *h_PurLead_1_2GeV_eta25 =  (TH1F*)HistFile->Get("Pur_Lead_1_2GeV_eta25");
  h_PurLead_1_2GeV_eta25->SetStats(kFALSE);
  h_PurLead_1_2GeV_eta25->SetTitle("matched cluster purity in 2.0 < |#eta_{EM2}| < 2.5");
  h_PurLead_1_2GeV_eta25 ->GetXaxis()->SetTitle("#rho_{clu}");
  h_PurLead_1_2GeV_eta25->SetLineColor(kBlue);
  h_PurLead_1_2GeV_eta25->Rebin(2);
  h_PurLead_1_2GeV_eta25->Scale(1.0/h_PurLead_1_2GeV_eta25->Integral());
  h_PurLead_1_2GeV_eta25->GetYaxis()->SetRangeUser(0,1.2);
  h_PurLead_1_2GeV_eta25->Draw();
  TH1 *h_PurLead_2_5GeV_eta25 =  (TH1F*)HistFile->Get("Pur_Lead_2_5GeV_eta25");
  h_PurLead_2_5GeV_eta25->SetLineColor(kRed);
  h_PurLead_2_5GeV_eta25->Rebin(2);
  h_PurLead_2_5GeV_eta25->Scale(1.0/h_PurLead_2_5GeV_eta25->Integral());
  h_PurLead_2_5GeV_eta25->Draw("same");
  TH1 *h_PurLead_5_10GeV_eta25 =  (TH1F*)HistFile->Get("Pur_Lead_5_10GeV_eta25");
  h_PurLead_5_10GeV_eta25->SetLineColor(kGreen+3);
  h_PurLead_5_10GeV_eta25->Rebin(2);
  h_PurLead_5_10GeV_eta25->Scale(1.0/h_PurLead_5_10GeV_eta25->Integral());
  h_PurLead_5_10GeV_eta25->Draw("same");
  TH1 *h_PurLead_10GeV_eta25 =  (TH1F*)HistFile->Get("Pur_Lead_10GeV_eta25");
  h_PurLead_10GeV_eta25->SetLineColor(kMagenta);
  h_PurLead_10GeV_eta25->Rebin(2);
  h_PurLead_10GeV_eta25->Scale(1.0/h_PurLead_10GeV_eta25->Integral());
  h_PurLead_10GeV_eta25->Draw("same");
  MyPurLegend->Draw();
  TPaveText *EtaRange3 = new TPaveText(0.20, 0.57, 0.40, 0.67,"NDC");
  EtaRange3->AddText("2.0 < |#eta| < 2.5");
  EtaRange3->SetFillColor(0);
  EtaRange3->SetLineColor(0);
  EtaRange3->SetShadowColor(0);
  EtaRange3->Draw();

  
  return;
}


//////////////////////////////
void NCluster09Lead()
{
  TCanvas* Can_NCluster09Lead = new TCanvas("NCluster09Lead","NCluster09Lead",900, 800);
  Can_NCluster09Lead->Divide(3,2);

  Can_NCluster09Lead->cd(1);
  TH1 *h_nClus09_1_2GeV_eta1 = (TH1F*)HistFile->Get("NClus_09_1_2GeV_eta1"); 
  h_nClus09_1_2GeV_eta1->SetStats(kFALSE);
  h_nClus09_1_2GeV_eta1->SetTitle(" ");
  h_nClus09_1_2GeV_eta1 ->GetXaxis()->SetTitle("N_{clu}(#Sigma E^{true}>90%)");
  h_nClus09_1_2GeV_eta1 ->GetYaxis()->SetTitle("Fraction of particles");
  h_nClus09_1_2GeV_eta1 ->SetLineColor(kBlue);
  h_nClus09_1_2GeV_eta1 ->SetFillColor(kBlue);h_nClus09_1_2GeV_eta1 ->SetFillStyle(3003);
  h_nClus09_1_2GeV_eta1->Scale(1.0/h_nClus09_1_2GeV_eta1->Integral());
  h_nClus09_1_2GeV_eta1->GetYaxis()->SetRangeUser(0,0.7);
  h_nClus09_1_2GeV_eta1->GetXaxis()->SetTitleSize(0.05);
  h_nClus09_1_2GeV_eta1->Draw();
  TH1 *h_nClus09_2_5GeV_eta1 =  (TH1F*)HistFile->Get("NClus_09_2_5GeV_eta1");
  h_nClus09_2_5GeV_eta1->SetLineColor(kRed);
  h_nClus09_2_5GeV_eta1 ->SetFillColor(kRed);h_nClus09_2_5GeV_eta1 ->SetFillStyle(3003);
  h_nClus09_2_5GeV_eta1->Scale(1.0/h_nClus09_2_5GeV_eta1->Integral());
  h_nClus09_2_5GeV_eta1->Draw("same");
  TH1 *h_nClus09_5_10GeV_eta1 =  (TH1F*)HistFile->Get("NClus_09_5_10GeV_eta1");
  h_nClus09_5_10GeV_eta1->SetLineColor(kGreen+3);
  h_nClus09_5_10GeV_eta1 ->SetFillColor(kGreen);h_nClus09_5_10GeV_eta1 ->SetFillStyle(3003);
  h_nClus09_5_10GeV_eta1->Scale(1.0/h_nClus09_5_10GeV_eta1->Integral());
  h_nClus09_5_10GeV_eta1->Draw("same");
  TH1 *h_nClus09_10GeV_eta1 =  (TH1F*)HistFile->Get("NClus_09_10GeV_eta1");
  h_nClus09_10GeV_eta1->SetLineColor(kMagenta);
  h_nClus09_10GeV_eta1 ->SetFillColor(kMagenta);h_nClus09_10GeV_eta1 ->SetFillStyle(3003);
  if(h_nClus09_10GeV_eta1)h_nClus09_10GeV_eta1->Scale(1.0/h_nClus09_10GeV_eta1->Integral());
  h_nClus09_10GeV_eta1->Draw("same");
  
  TLegend *MyClus09 = new TLegend (0.50, 0.70, 0.87, 0.85);
  MyClus09->SetFillStyle(1001);
  MyClus09->SetFillColor(10);
  MyClus09->SetLineColor(0);
  MyClus09->AddEntry(  h_nClus09_1_2GeV_eta1, "1 < p_{T}^{track} < 2 GeV", "l");
  MyClus09->AddEntry(  h_nClus09_2_5GeV_eta1, "2 < p_{T}^{track} < 5 GeV", "l");
  MyClus09->AddEntry(  h_nClus09_5_10GeV_eta1,  " 5 < p_{T}^{track} < 10 GeV", "l");
  MyClus09->AddEntry(  h_nClus09_10GeV_eta1,  " p_{T}^{track} > 10 GeV", "l");
  MyClus09->Draw();

  TPaveText *EtaRange = new TPaveText(0.55, 0.57, 0.75, 0.67,"NDC");
  EtaRange->AddText("|#eta| < 1.0");
  EtaRange->SetFillColor(0);
  EtaRange->SetLineColor(0);
  EtaRange->SetShadowColor(0);
  EtaRange->Draw();
    
  Can_NCluster09Lead->cd(2);
  TH1 *h_nClus09_1_2GeV_eta2 =  (TH1F*)HistFile->Get("NClus_09_1_2GeV_eta2");
  h_nClus09_1_2GeV_eta2->SetStats(kFALSE);
  h_nClus09_1_2GeV_eta2->SetTitle(" ");
  h_nClus09_1_2GeV_eta2 ->GetXaxis()->SetTitle("N_{clu}(#Sigma E^{true}>90%)");
  h_nClus09_1_2GeV_eta2->SetLineColor(kBlue);
  h_nClus09_1_2GeV_eta2 ->SetFillColor(kBlue);h_nClus09_1_2GeV_eta2 ->SetFillStyle(3003);
  if(h_nClus09_1_2GeV_eta2)h_nClus09_1_2GeV_eta2->Scale(1.0/h_nClus09_1_2GeV_eta2->Integral());
  h_nClus09_1_2GeV_eta2->GetYaxis()->SetRangeUser(0,0.7);
  h_nClus09_1_2GeV_eta2->GetXaxis()->SetTitleSize(0.05);
  h_nClus09_1_2GeV_eta2->Draw();
  TH1 *h_nClus09_2_5GeV_eta2 =  (TH1F*)HistFile->Get("NClus_09_2_5GeV_eta2");
  h_nClus09_2_5GeV_eta2->SetLineColor(kRed);
  h_nClus09_2_5GeV_eta2 ->SetFillColor(kRed);h_nClus09_2_5GeV_eta2 ->SetFillStyle(3003);
  if(h_nClus09_2_5GeV_eta2)h_nClus09_2_5GeV_eta2->Scale(1.0/h_nClus09_2_5GeV_eta2->Integral());
  h_nClus09_2_5GeV_eta2->Draw("same");
  TH1 *h_nClus09_5_10GeV_eta2 =  (TH1F*)HistFile->Get("NClus_09_5_10GeV_eta2");
  h_nClus09_5_10GeV_eta2->SetLineColor(kGreen+3);
  h_nClus09_5_10GeV_eta2 ->SetFillColor(kGreen);h_nClus09_5_10GeV_eta2 ->SetFillStyle(3003);
  if(h_nClus09_5_10GeV_eta2)h_nClus09_5_10GeV_eta2->Scale(1.0/h_nClus09_5_10GeV_eta2->Integral());
  h_nClus09_5_10GeV_eta2->Draw("same");
  TH1 *h_nClus09_10GeV_eta2 =  (TH1F*)HistFile->Get("NClus_09_10GeV_eta2");
  h_nClus09_10GeV_eta2->SetLineColor(kMagenta);
  h_nClus09_10GeV_eta2 ->SetFillColor(kMagenta);h_nClus09_10GeV_eta2 ->SetFillStyle(3003);
  if(h_nClus09_10GeV_eta2)h_nClus09_10GeV_eta2->Scale(1.0/h_nClus09_10GeV_eta2->Integral());
  h_nClus09_10GeV_eta2->Draw("same");
  MyClus09->Draw();

  TPaveText *EtaRange2 = new TPaveText(0.50, 0.57, 0.80, 0.67,"NDC");
  EtaRange2->AddText("1.0 <|#eta| < 2.0");
  EtaRange2->SetFillColor(0);
  EtaRange2->SetLineColor(0);
  EtaRange2->SetShadowColor(0);
  EtaRange2->Draw();

  Can_NCluster09Lead->cd(3);
  TH1 *h_nClus09_1_2GeV_eta25 =  (TH1F*)HistFile->Get("NClus_09_1_2GeV_eta25");
  h_nClus09_1_2GeV_eta25->SetStats(kFALSE);
  h_nClus09_1_2GeV_eta25->SetTitle(" ");
  h_nClus09_1_2GeV_eta25 ->GetXaxis()->SetTitle("N_{clu}(#Sigma E^{true}>90%)");
  h_nClus09_1_2GeV_eta25->SetLineColor(kBlue);
  h_nClus09_1_2GeV_eta25 ->SetFillColor(kBlue);h_nClus09_1_2GeV_eta25 ->SetFillStyle(3003);
  h_nClus09_1_2GeV_eta25->Scale(1.0/h_nClus09_1_2GeV_eta25->Integral());
  h_nClus09_1_2GeV_eta25->GetYaxis()->SetRangeUser(0,0.7);
  h_nClus09_1_2GeV_eta25->GetXaxis()->SetTitleSize(0.05);
  h_nClus09_1_2GeV_eta25->Draw();
  TH1 *h_nClus09_2_5GeV_eta25 =  (TH1F*)HistFile->Get("NClus_09_2_5GeV_eta25");
  h_nClus09_2_5GeV_eta25->SetLineColor(kRed);
  h_nClus09_2_5GeV_eta25 ->SetFillColor(kRed);h_nClus09_2_5GeV_eta25 ->SetFillStyle(3003);
  h_nClus09_2_5GeV_eta25->Scale(1.0/h_nClus09_2_5GeV_eta25->Integral());
  h_nClus09_2_5GeV_eta25->Draw("same");
  TH1 *h_nClus09_5_10GeV_eta25 =  (TH1F*)HistFile->Get("NClus_09_5_10GeV_eta25");
  h_nClus09_5_10GeV_eta25->SetLineColor(kGreen+3);
  h_nClus09_5_10GeV_eta25 ->SetFillColor(kGreen);h_nClus09_5_10GeV_eta25 ->SetFillStyle(3003);
  h_nClus09_5_10GeV_eta25->Scale(1.0/h_nClus09_5_10GeV_eta25->Integral());
  h_nClus09_5_10GeV_eta25->Draw("same");
  TH1 *h_nClus09_10GeV_eta25 =  (TH1F*)HistFile->Get("NClus_09_10GeV_eta25");
  h_nClus09_10GeV_eta25->SetLineColor(kMagenta);
  h_nClus09_10GeV_eta25 ->SetFillColor(kMagenta);h_nClus09_10GeV_eta25 ->SetFillStyle(3003);
  h_nClus09_10GeV_eta25->Scale(1.0/h_nClus09_10GeV_eta25->Integral());
  h_nClus09_10GeV_eta25->Draw("same");
  MyClus09->Draw();
  
  TPaveText *EtaRange3 = new TPaveText(0.50, 0.57, 0.80, 0.67,"NDC");
  EtaRange3->AddText("2.0 < |#eta| < 2.5");
  EtaRange3->SetFillColor(0);
  EtaRange3->SetLineColor(0);
  EtaRange3->SetShadowColor(0);
  EtaRange3->Draw();

  
  return;
}




void DeltaE_09_07()
{
  TCanvas* Can_DeltaE_09_07 = new TCanvas("Can_DeltaE_09_07","Can_DeltaE_09_07",900, 950);
  Can_DeltaE_09_07->Divide(3,4);
  
  Can_DeltaE_09_07->cd(1);
  TH1 *h_DeltaE09_1_2GeV_eta1 = (TH1F*)HistFile->Get("DeltaE_1_2GeV_eta1");
  h_DeltaE09_1_2GeV_eta1->SetStats(kFALSE);
  h_DeltaE09_1_2GeV_eta1->SetTitle(" ");
  h_DeltaE09_1_2GeV_eta1 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_1_2GeV_eta1 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_1_2GeV_eta1 ->SetLineColor(kRed);
  h_DeltaE09_1_2GeV_eta1 ->SetFillColor(kRed);
  h_DeltaE09_1_2GeV_eta1 ->SetFillStyle(3004);
  h_DeltaE09_1_2GeV_eta1->Scale(1.0/h_DeltaE09_1_2GeV_eta1->Integral());
  h_DeltaE09_1_2GeV_eta1->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_1_2GeV_eta1->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_1_2GeV_eta1->Draw();
  TH1 *h_DeltaE07_1_2GeV_eta1 = (TH1F*)HistFile->Get("DeltaE07_1_2GeV_eta1");
  h_DeltaE07_1_2GeV_eta1 ->SetLineColor(kBlue);
  h_DeltaE07_1_2GeV_eta1 ->SetFillColor(kBlue);
  h_DeltaE07_1_2GeV_eta1 ->SetFillStyle(3004);
  h_DeltaE07_1_2GeV_eta1->Scale(1.0/h_DeltaE07_1_2GeV_eta1->Integral());
  h_DeltaE07_1_2GeV_eta1->Draw("same");
  
  TLegend *MyClus09 = new TLegend (0.70, 0.70, 0.87, 0.85);
  MyClus09->SetFillStyle(1001);
  MyClus09->SetFillColor(10);
  MyClus09->SetLineColor(0);
  MyClus09->AddEntry(  h_DeltaE09_1_2GeV_eta1, "#varepsilon_{clu}>90%", "l");
  MyClus09->AddEntry(  h_DeltaE07_1_2GeV_eta1, "#varepsilon_{clu}<70%", "l");
  MyClus09->Draw();
  
  TPaveText *EtaRange = new TPaveText(0.15, 0.7, 0.55, 0.87,"NDC");
  EtaRange->AddText("1 < p_{T}^{track} < 2 GeV   |#eta| < 1.0");
  EtaRange->SetFillColor(0);
  EtaRange->SetLineColor(0);
  EtaRange->SetShadowColor(0);
  EtaRange->Draw();
  
  Can_DeltaE_09_07->cd(2);
  TH1 *h_DeltaE09_1_2GeV_eta2 = (TH1F*)HistFile->Get("DeltaE_1_2GeV_eta2");
  h_DeltaE09_1_2GeV_eta2->SetStats(kFALSE);
  h_DeltaE09_1_2GeV_eta2->SetTitle(" ");
  h_DeltaE09_1_2GeV_eta2 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_1_2GeV_eta2 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_1_2GeV_eta2 ->SetLineColor(kRed);
  h_DeltaE09_1_2GeV_eta2 ->SetFillColor(kRed);
  h_DeltaE09_1_2GeV_eta2 ->SetFillStyle(3005);
  h_DeltaE09_1_2GeV_eta2->Scale(1.0/h_DeltaE09_1_2GeV_eta2->Integral());
  h_DeltaE09_1_2GeV_eta2->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_1_2GeV_eta2->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_1_2GeV_eta2->Draw();
  TH1 *h_DeltaE07_1_2GeV_eta2 = (TH1F*)HistFile->Get("DeltaE07_1_2GeV_eta2");
  h_DeltaE07_1_2GeV_eta2 ->SetLineColor(kBlue);
  h_DeltaE07_1_2GeV_eta2 ->SetFillColor(kBlue);
  h_DeltaE07_1_2GeV_eta2 ->SetFillStyle(3005);
  h_DeltaE07_1_2GeV_eta2->Scale(1.0/h_DeltaE07_1_2GeV_eta2->Integral());
  h_DeltaE07_1_2GeV_eta2->Draw("same");
  MyClus09->Draw();
  
  TPaveText *EtaRange2 = new TPaveText(0.15, 0.7, 0.65, 0.87,"NDC");
  EtaRange2->AddText("1 < p_{T}^{track} < 2 GeV  1.0 < |#eta| < 2.0");
  EtaRange2->SetFillColor(0);
  EtaRange2->SetLineColor(0);
  EtaRange2->SetShadowColor(0);
  EtaRange2->Draw();
  
  Can_DeltaE_09_07->cd(3);
  TH1 *h_DeltaE09_1_2GeV_eta25 = (TH1F*)HistFile->Get("DeltaE_1_2GeV_eta25");
  h_DeltaE09_1_2GeV_eta25->SetStats(kFALSE);
  h_DeltaE09_1_2GeV_eta25->SetTitle(" ");
  h_DeltaE09_1_2GeV_eta25 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_1_2GeV_eta25 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_1_2GeV_eta25 ->SetLineColor(kRed);
  h_DeltaE09_1_2GeV_eta25 ->SetFillColor(kRed);
  h_DeltaE09_1_2GeV_eta25 ->SetFillStyle(3005);
  h_DeltaE09_1_2GeV_eta25->Scale(1.0/h_DeltaE09_1_2GeV_eta25->Integral());
  h_DeltaE09_1_2GeV_eta25->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_1_2GeV_eta25->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_1_2GeV_eta25->Draw();
  TH1 *h_DeltaE07_1_2GeV_eta25 = (TH1F*)HistFile->Get("DeltaE07_1_2GeV_eta25");
  h_DeltaE07_1_2GeV_eta25 ->SetLineColor(kBlue);
  h_DeltaE07_1_2GeV_eta25 ->SetFillColor(kBlue);
  h_DeltaE07_1_2GeV_eta25 ->SetFillStyle(3005);
  h_DeltaE07_1_2GeV_eta25->Scale(1.0/h_DeltaE07_1_2GeV_eta25->Integral());
  h_DeltaE07_1_2GeV_eta25->Draw("same");
  MyClus09->Draw();
  
  TPaveText *EtaRange3 = new TPaveText(0.15, 0.7, 0.65, 0.87,"NDC");
  EtaRange3->AddText("1 < p_{T}^{track} < 2 GeV  2.0 < |#eta| < 2.5");
  EtaRange3->SetFillColor(0);
  EtaRange3->SetLineColor(0);
  EtaRange3->SetShadowColor(0);
  EtaRange3->Draw();
  
  //Other pt range
  Can_DeltaE_09_07->cd(4);
  TH1 *h_DeltaE09_2_5GeV_eta1 = (TH1F*)HistFile->Get("DeltaE_2_5GeV_eta1");
  h_DeltaE09_2_5GeV_eta1->SetStats(kFALSE);
  h_DeltaE09_2_5GeV_eta1->SetTitle(" ");
  h_DeltaE09_2_5GeV_eta1 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_2_5GeV_eta1 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_2_5GeV_eta1 ->SetLineColor(kRed);
  h_DeltaE09_2_5GeV_eta1 ->SetFillColor(kRed);
  h_DeltaE09_2_5GeV_eta1 ->SetFillStyle(3004);
  h_DeltaE09_2_5GeV_eta1->Scale(1.0/h_DeltaE09_2_5GeV_eta1->Integral());
  h_DeltaE09_2_5GeV_eta1->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_2_5GeV_eta1->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_2_5GeV_eta1->Draw();
  TH1 *h_DeltaE07_2_5GeV_eta1 = (TH1F*)HistFile->Get("DeltaE07_2_5GeV_eta1");
  h_DeltaE07_2_5GeV_eta1 ->SetLineColor(kBlue);
  h_DeltaE07_2_5GeV_eta1 ->SetFillColor(kBlue);
  h_DeltaE07_2_5GeV_eta1 ->SetFillStyle(3005);
  h_DeltaE07_2_5GeV_eta1->Scale(1.0/h_DeltaE07_2_5GeV_eta1->Integral());
  h_DeltaE07_2_5GeV_eta1->Draw("same");
  MyClus09->Draw();
  
  TPaveText *EtaRange4 = new TPaveText(0.15, 0.7, 0.55, 0.87,"NDC");
  EtaRange4->AddText("2 < p_{T}^{track} < 5 GeV   |#eta| < 1.0");
  EtaRange4->SetFillColor(0);
  EtaRange4->SetLineColor(0);
  EtaRange4->SetShadowColor(0);
  EtaRange4->Draw();
 
  Can_DeltaE_09_07->cd(5);
  TH1 *h_DeltaE09_2_5GeV_eta2 = (TH1F*)HistFile->Get("DeltaE_2_5GeV_eta2");
  h_DeltaE09_2_5GeV_eta2->SetStats(kFALSE);
  h_DeltaE09_2_5GeV_eta2->SetTitle(" ");
  h_DeltaE09_2_5GeV_eta2 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_2_5GeV_eta2 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_2_5GeV_eta2 ->SetLineColor(kRed);
  h_DeltaE09_2_5GeV_eta2 ->SetFillColor(kRed);
  h_DeltaE09_2_5GeV_eta2 ->SetFillStyle(3005);
  h_DeltaE09_2_5GeV_eta2->Scale(1.0/h_DeltaE09_2_5GeV_eta2->Integral());
  h_DeltaE09_2_5GeV_eta2->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_2_5GeV_eta2->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_2_5GeV_eta2->Draw();
  TH1 *h_DeltaE07_2_5GeV_eta2 = (TH1F*)HistFile->Get("DeltaE07_2_5GeV_eta2");
  h_DeltaE07_2_5GeV_eta2 ->SetLineColor(kBlue);
  h_DeltaE07_2_5GeV_eta2 ->SetFillColor(kBlue);
  h_DeltaE07_2_5GeV_eta2 ->SetFillStyle(3005);
  h_DeltaE07_2_5GeV_eta2->Scale(1.0/h_DeltaE07_2_5GeV_eta2->Integral());
  h_DeltaE07_2_5GeV_eta2->Draw("same");
  MyClus09->Draw();
  
  TPaveText *EtaRange5 = new TPaveText(0.15, 0.7, 0.65, 0.87,"NDC");
  EtaRange5->AddText("2 < p_{T}^{track} < 5 GeV  1.0 < |#eta| < 2.0");
  EtaRange5->SetFillColor(0);
  EtaRange5->SetLineColor(0);
  EtaRange5->SetShadowColor(0);
  EtaRange5->Draw();
  
  Can_DeltaE_09_07->cd(6);
  TH1 *h_DeltaE09_2_5GeV_eta25 = (TH1F*)HistFile->Get("DeltaE_2_5GeV_eta25");
  h_DeltaE09_2_5GeV_eta25->SetStats(kFALSE);
  h_DeltaE09_2_5GeV_eta25->SetTitle(" ");
  h_DeltaE09_2_5GeV_eta25 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_2_5GeV_eta25 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_2_5GeV_eta25 ->SetLineColor(kRed);
  h_DeltaE09_2_5GeV_eta25 ->SetFillColor(kRed);
  h_DeltaE09_2_5GeV_eta25 ->SetFillStyle(3005);
  h_DeltaE09_2_5GeV_eta25->Scale(1.0/h_DeltaE09_2_5GeV_eta25->Integral());
  h_DeltaE09_2_5GeV_eta25->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_2_5GeV_eta25->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_2_5GeV_eta25->Draw();
  TH1 *h_DeltaE07_2_5GeV_eta25 = (TH1F*)HistFile->Get("DeltaE07_2_5GeV_eta25");
  h_DeltaE07_2_5GeV_eta25 ->SetLineColor(kBlue);
  h_DeltaE07_2_5GeV_eta25 ->SetFillColor(kBlue);
  h_DeltaE07_2_5GeV_eta25 ->SetFillStyle(3005);
  h_DeltaE07_2_5GeV_eta25->Scale(1.0/h_DeltaE07_2_5GeV_eta25->Integral());
  h_DeltaE07_2_5GeV_eta25->Draw("same");
  MyClus09->Draw();
  
  TPaveText *EtaRange6 = new TPaveText(0.15, 0.7, 0.65, 0.87,"NDC");
  EtaRange6->AddText("2 < p_{T}^{track} < 5 GeV  2.0 < |#eta| < 2.5");
  EtaRange6->SetFillColor(0);
  EtaRange6->SetLineColor(0);
  EtaRange6->SetShadowColor(0);
  EtaRange6->Draw();
  
  //Other pt range
  Can_DeltaE_09_07->cd(7);
  TH1 *h_DeltaE09_5_10GeV_eta1 = (TH1F*)HistFile->Get("DeltaE_5_10GeV_eta1");
  h_DeltaE09_5_10GeV_eta1->SetStats(kFALSE);
  h_DeltaE09_5_10GeV_eta1->SetTitle(" ");
  h_DeltaE09_5_10GeV_eta1 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_5_10GeV_eta1 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_5_10GeV_eta1 ->SetLineColor(kRed);
  h_DeltaE09_5_10GeV_eta1 ->SetFillColor(kRed);
  h_DeltaE09_5_10GeV_eta1 ->SetFillStyle(3004);
  h_DeltaE09_5_10GeV_eta1->Scale(1.0/h_DeltaE09_5_10GeV_eta1->Integral());
  h_DeltaE09_5_10GeV_eta1->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_5_10GeV_eta1->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_5_10GeV_eta1->Draw();
  TH1 *h_DeltaE07_5_10GeV_eta1 = (TH1F*)HistFile->Get("DeltaE07_5_10GeV_eta1");
  h_DeltaE07_5_10GeV_eta1 ->SetLineColor(kBlue);
  h_DeltaE07_5_10GeV_eta1 ->SetFillColor(kBlue);
  h_DeltaE07_5_10GeV_eta1 ->SetFillStyle(3005);
  h_DeltaE07_5_10GeV_eta1->Scale(1.0/h_DeltaE07_5_10GeV_eta1->Integral());
  h_DeltaE07_5_10GeV_eta1->Draw("same");
  MyClus09->Draw();
  
  TPaveText *EtaRange7 = new TPaveText(0.15, 0.7, 0.55, 0.87,"NDC");
  EtaRange7->AddText("5 < p_{T}^{track} < 10 GeV   |#eta| < 1.0");
  EtaRange7->SetFillColor(0);
  EtaRange7->SetLineColor(0);
  EtaRange7->SetShadowColor(0);
  EtaRange7->Draw();
  
  Can_DeltaE_09_07->cd(8);
  TH1 *h_DeltaE09_5_10GeV_eta2 = (TH1F*)HistFile->Get("DeltaE_5_10GeV_eta2");
  h_DeltaE09_5_10GeV_eta2->SetStats(kFALSE);
  h_DeltaE09_5_10GeV_eta2->SetTitle(" ");
  h_DeltaE09_5_10GeV_eta2 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_5_10GeV_eta2 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_5_10GeV_eta2 ->SetLineColor(kRed);
  h_DeltaE09_5_10GeV_eta2 ->SetFillColor(kRed);
  h_DeltaE09_5_10GeV_eta2 ->SetFillStyle(3005);
  h_DeltaE09_5_10GeV_eta2->Scale(1.0/h_DeltaE09_5_10GeV_eta2->Integral());
  h_DeltaE09_5_10GeV_eta2->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_5_10GeV_eta2->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_5_10GeV_eta2->Draw();
  TH1 *h_DeltaE07_5_10GeV_eta2 = (TH1F*)HistFile->Get("DeltaE07_5_10GeV_eta2");
  h_DeltaE07_5_10GeV_eta2 ->SetLineColor(kBlue);
  h_DeltaE07_5_10GeV_eta2 ->SetFillColor(kBlue);
  h_DeltaE07_5_10GeV_eta2 ->SetFillStyle(3005);
  h_DeltaE07_5_10GeV_eta2->Scale(1.0/h_DeltaE07_5_10GeV_eta2->Integral());
  h_DeltaE07_5_10GeV_eta2->Draw("same");
  MyClus09->Draw();
  
  TPaveText *EtaRange8 = new TPaveText(0.15, 0.7, 0.65, 0.87,"NDC");
  EtaRange8->AddText("5 < p_{T}^{track} < 10 GeV  1.0 < |#eta| < 2.0");
  EtaRange8->SetFillColor(0);
  EtaRange8->SetLineColor(0);
  EtaRange8->SetShadowColor(0);
  EtaRange8->Draw();

  Can_DeltaE_09_07->cd(9);
  h_DeltaE09_1_2GeV_eta25->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_1_2GeV_eta25->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_1_2GeV_eta25->Draw();
  TH1 *h_DeltaE07_1_2GeV_eta25 = (TH1F*)HistFile->Get("DeltaE07_1_2GeV_eta25");
  h_DeltaE07_1_2GeV_eta25 ->SetLineColor(kBlue);
  h_DeltaE07_1_2GeV_eta25 ->SetFillColor(kBlue);
  h_DeltaE07_1_2GeV_eta25 ->SetFillStyle(3005);
  h_DeltaE07_1_2GeV_eta25->Scale(1.0/h_DeltaE07_1_2GeV_eta25->Integral());
  h_DeltaE07_1_2GeV_eta25->Draw("same");
  MyClus09->Draw();
  
  TPaveText *EtaRange3 = new TPaveText(0.15, 0.7, 0.65, 0.87,"NDC");
  EtaRange3->AddText("1 < p_{T}^{track} < 2 GeV  2.0 < |#eta| < 2.5");
  EtaRange3->SetFillColor(0);
  EtaRange3->SetLineColor(0);
  EtaRange3->SetShadowColor(0);
  EtaRange3->Draw();

  //Other pt range
  Can_DeltaE_09_07->cd(4);
  TH1 *h_DeltaE09_2_5GeV_eta1 = (TH1F*)HistFile->Get("DeltaE_2_5GeV_eta1");
  h_DeltaE09_2_5GeV_eta1->SetStats(kFALSE);
  h_DeltaE09_2_5GeV_eta1->SetTitle(" ");
  h_DeltaE09_2_5GeV_eta1 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_2_5GeV_eta1 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_2_5GeV_eta1 ->SetLineColor(kRed);
  h_DeltaE09_2_5GeV_eta1 ->SetFillColor(kRed);
  h_DeltaE09_2_5GeV_eta1 ->SetFillStyle(3004);
  h_DeltaE09_2_5GeV_eta1->Scale(1.0/h_DeltaE09_2_5GeV_eta1->Integral());
  h_DeltaE09_2_5GeV_eta1->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_2_5GeV_eta1->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_2_5GeV_eta1->Draw();
  TH1 *h_DeltaE07_2_5GeV_eta1 = (TH1F*)HistFile->Get("DeltaE07_2_5GeV_eta1");
  h_DeltaE07_2_5GeV_eta1 ->SetLineColor(kBlue);
  h_DeltaE07_2_5GeV_eta1 ->SetFillColor(kBlue);
  h_DeltaE07_2_5GeV_eta1 ->SetFillStyle(3005);
  h_DeltaE07_2_5GeV_eta1->Scale(1.0/h_DeltaE07_2_5GeV_eta1->Integral());
  h_DeltaE07_2_5GeV_eta1->Draw("same");
  MyClus09->Draw();

  TPaveText *EtaRange4 = new TPaveText(0.15, 0.7, 0.55, 0.87,"NDC");
  EtaRange4->AddText("2 < p_{T}^{track} < 5 GeV   |#eta| < 1.0");
  EtaRange4->SetFillColor(0);
  EtaRange4->SetLineColor(0);
  EtaRange4->SetShadowColor(0);
  EtaRange4->Draw();
  
  Can_DeltaE_09_07->cd(5);
  TH1 *h_DeltaE09_2_5GeV_eta2 = (TH1F*)HistFile->Get("DeltaE_2_5GeV_eta2");
  h_DeltaE09_2_5GeV_eta2->SetStats(kFALSE);
  h_DeltaE09_2_5GeV_eta2->SetTitle(" ");
  h_DeltaE09_2_5GeV_eta2 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_2_5GeV_eta2 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_2_5GeV_eta2 ->SetLineColor(kRed);
  h_DeltaE09_2_5GeV_eta2 ->SetFillColor(kRed);
  h_DeltaE09_2_5GeV_eta2 ->SetFillStyle(3005);
  h_DeltaE09_2_5GeV_eta2->Scale(1.0/h_DeltaE09_2_5GeV_eta2->Integral());
  h_DeltaE09_2_5GeV_eta2->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_2_5GeV_eta2->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_2_5GeV_eta2->Draw();
  TH1 *h_DeltaE07_2_5GeV_eta2 = (TH1F*)HistFile->Get("DeltaE07_2_5GeV_eta2");
  h_DeltaE07_2_5GeV_eta2 ->SetLineColor(kBlue);
  h_DeltaE07_2_5GeV_eta2 ->SetFillColor(kBlue);
  h_DeltaE07_2_5GeV_eta2 ->SetFillStyle(3005);
  h_DeltaE07_2_5GeV_eta2->Scale(1.0/h_DeltaE07_2_5GeV_eta2->Integral());
  h_DeltaE07_2_5GeV_eta2->Draw("same");
  MyClus09->Draw();
  
  TPaveText *EtaRange5 = new TPaveText(0.15, 0.7, 0.65, 0.87,"NDC");
  EtaRange5->AddText("2 < p_{T}^{track} < 5 GeV  1.0 < |#eta| < 2.0");
  EtaRange5->SetFillColor(0);
  EtaRange5->SetLineColor(0);
  EtaRange5->SetShadowColor(0);
  EtaRange5->Draw();
  
  Can_DeltaE_09_07->cd(6);
  TH1 *h_DeltaE09_2_5GeV_eta25 = (TH1F*)HistFile->Get("DeltaE_2_5GeV_eta25");
  h_DeltaE09_2_5GeV_eta25->SetStats(kFALSE);
  h_DeltaE09_2_5GeV_eta25->SetTitle(" ");
  h_DeltaE09_2_5GeV_eta25 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_2_5GeV_eta25 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_2_5GeV_eta25 ->SetLineColor(kRed);
  h_DeltaE09_2_5GeV_eta25 ->SetFillColor(kRed);
  h_DeltaE09_2_5GeV_eta25 ->SetFillStyle(3005);
  h_DeltaE09_2_5GeV_eta25->Scale(1.0/h_DeltaE09_2_5GeV_eta25->Integral());
  h_DeltaE09_2_5GeV_eta25->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_2_5GeV_eta25->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_2_5GeV_eta25->Draw();
  TH1 *h_DeltaE07_2_5GeV_eta25 = (TH1F*)HistFile->Get("DeltaE07_2_5GeV_eta25");
  h_DeltaE07_2_5GeV_eta25 ->SetLineColor(kBlue);
  h_DeltaE07_2_5GeV_eta25 ->SetFillColor(kBlue);
  h_DeltaE07_2_5GeV_eta25 ->SetFillStyle(3005);
  h_DeltaE07_2_5GeV_eta25->Scale(1.0/h_DeltaE07_2_5GeV_eta25->Integral());
  h_DeltaE07_2_5GeV_eta25->Draw("same");
  MyClus09->Draw();
  
  TPaveText *EtaRange6 = new TPaveText(0.15, 0.7, 0.65, 0.87,"NDC");
  EtaRange6->AddText("2 < p_{T}^{track} < 5 GeV  2.0 < |#eta| < 2.5");
  EtaRange6->SetFillColor(0);
  EtaRange6->SetLineColor(0);
  EtaRange6->SetShadowColor(0);
  EtaRange6->Draw();
  
  //Other pt range
  Can_DeltaE_09_07->cd(7);
  TH1 *h_DeltaE09_5_10GeV_eta1 = (TH1F*)HistFile->Get("DeltaE_5_10GeV_eta1");
  h_DeltaE09_5_10GeV_eta1->SetStats(kFALSE);
  h_DeltaE09_5_10GeV_eta1->SetTitle(" ");
  h_DeltaE09_5_10GeV_eta1 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_5_10GeV_eta1 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_5_10GeV_eta1 ->SetLineColor(kRed);
  h_DeltaE09_5_10GeV_eta1 ->SetFillColor(kRed);
  h_DeltaE09_5_10GeV_eta1 ->SetFillStyle(3004);
  h_DeltaE09_5_10GeV_eta1->Scale(1.0/h_DeltaE09_5_10GeV_eta1->Integral());
  h_DeltaE09_5_10GeV_eta1->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_5_10GeV_eta1->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_5_10GeV_eta1->Draw();
  TH1 *h_DeltaE07_5_10GeV_eta1 = (TH1F*)HistFile->Get("DeltaE07_5_10GeV_eta1");
  h_DeltaE07_5_10GeV_eta1 ->SetLineColor(kBlue);
  h_DeltaE07_5_10GeV_eta1 ->SetFillColor(kBlue);
  h_DeltaE07_5_10GeV_eta1 ->SetFillStyle(3005);
  h_DeltaE07_5_10GeV_eta1->Scale(1.0/h_DeltaE07_5_10GeV_eta1->Integral());
  h_DeltaE07_5_10GeV_eta1->Draw("same");
  MyClus09->Draw();

  TPaveText *EtaRange7 = new TPaveText(0.15, 0.7, 0.55, 0.87,"NDC");
  EtaRange7->AddText("5 < p_{T}^{track} < 10 GeV   |#eta| < 1.0");
  EtaRange7->SetFillColor(0);
  EtaRange7->SetLineColor(0);
  EtaRange7->SetShadowColor(0);
  EtaRange7->Draw();
  
  Can_DeltaE_09_07->cd(8);
  TH1 *h_DeltaE09_5_10GeV_eta2 = (TH1F*)HistFile->Get("DeltaE_5_10GeV_eta2");
  h_DeltaE09_5_10GeV_eta2->SetStats(kFALSE);
  h_DeltaE09_5_10GeV_eta2->SetTitle(" ");
  h_DeltaE09_5_10GeV_eta2 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_5_10GeV_eta2 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_5_10GeV_eta2 ->SetLineColor(kRed);
  h_DeltaE09_5_10GeV_eta2 ->SetFillColor(kRed);
  h_DeltaE09_5_10GeV_eta2 ->SetFillStyle(3005);
  h_DeltaE09_5_10GeV_eta2->Scale(1.0/h_DeltaE09_5_10GeV_eta2->Integral());
  h_DeltaE09_5_10GeV_eta2->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_5_10GeV_eta2->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_5_10GeV_eta2->Draw();
  TH1 *h_DeltaE07_5_10GeV_eta2 = (TH1F*)HistFile->Get("DeltaE07_5_10GeV_eta2");
  h_DeltaE07_5_10GeV_eta2 ->SetLineColor(kBlue);
  h_DeltaE07_5_10GeV_eta2 ->SetFillColor(kBlue);
  h_DeltaE07_5_10GeV_eta2 ->SetFillStyle(3005);
  h_DeltaE07_5_10GeV_eta2->Scale(1.0/h_DeltaE07_5_10GeV_eta2->Integral());
  h_DeltaE07_5_10GeV_eta2->Draw("same");
  MyClus09->Draw();
  
  TPaveText *EtaRange8 = new TPaveText(0.15, 0.7, 0.65, 0.87,"NDC");
  EtaRange8->AddText("5 < p_{T}^{track} < 10 GeV  1.0 < |#eta| < 2.0");
  EtaRange8->SetFillColor(0);
  EtaRange8->SetLineColor(0);
  EtaRange8->SetShadowColor(0);
  EtaRange8->Draw();
  
  Can_DeltaE_09_07->cd(9);
  TH1 *h_DeltaE09_5_10GeV_eta25 = (TH1F*)HistFile->Get("DeltaE_5_10GeV_eta25");
  h_DeltaE09_5_10GeV_eta25->SetStats(kFALSE);
  h_DeltaE09_5_10GeV_eta25->SetTitle(" ");
  h_DeltaE09_5_10GeV_eta25 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_5_10GeV_eta25 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_5_10GeV_eta25 ->SetLineColor(kRed);
  h_DeltaE09_5_10GeV_eta25 ->SetFillColor(kRed);
  h_DeltaE09_5_10GeV_eta25 ->SetFillStyle(3005);
  h_DeltaE09_5_10GeV_eta25->Scale(1.0/h_DeltaE09_5_10GeV_eta25->Integral());
  h_DeltaE09_5_10GeV_eta25->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_5_10GeV_eta25->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_5_10GeV_eta25->Draw();
  TH1 *h_DeltaE07_5_10GeV_eta25 = (TH1F*)HistFile->Get("DeltaE07_5_10GeV_eta25");
  h_DeltaE07_5_10GeV_eta25 ->SetLineColor(kBlue);
  h_DeltaE07_5_10GeV_eta25 ->SetFillColor(kBlue);
  h_DeltaE07_5_10GeV_eta25 ->SetFillStyle(3005);
  h_DeltaE07_5_10GeV_eta25->Scale(1.0/h_DeltaE07_5_10GeV_eta25->Integral());
  h_DeltaE07_5_10GeV_eta25->Draw("same");
  MyClus09->Draw();
  
  TPaveText *EtaRange9 =  new TPaveText(0.15, 0.7, 0.65, 0.87,"NDC");
  TH1 *h_DeltaE09_5_10GeV_eta25 = (TH1F*)HistFile->Get("DeltaE_5_10GeV_eta25");
  h_DeltaE09_5_10GeV_eta25->SetStats(kFALSE);
  h_DeltaE09_5_10GeV_eta25->SetTitle(" ");
  h_DeltaE09_5_10GeV_eta25 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_5_10GeV_eta25 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_5_10GeV_eta25 ->SetLineColor(kRed);
  h_DeltaE09_5_10GeV_eta25 ->SetFillColor(kRed);
  h_DeltaE09_5_10GeV_eta25 ->SetFillStyle(3005);
  h_DeltaE09_5_10GeV_eta25->Scale(1.0/h_DeltaE09_5_10GeV_eta25->Integral());
  h_DeltaE09_5_10GeV_eta25->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_5_10GeV_eta25->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_5_10GeV_eta25->Draw();
  TH1 *h_DeltaE07_5_10GeV_eta25 = (TH1F*)HistFile->Get("DeltaE07_5_10GeV_eta25");
  h_DeltaE07_5_10GeV_eta25 ->SetLineColor(kBlue);
  h_DeltaE07_5_10GeV_eta25 ->SetFillColor(kBlue);
  h_DeltaE07_5_10GeV_eta25 ->SetFillStyle(3005);
  h_DeltaE07_5_10GeV_eta25->Scale(1.0/h_DeltaE07_5_10GeV_eta25->Integral());
  h_DeltaE07_5_10GeV_eta25->Draw("same");
  MyClus09->Draw();
  
  TPaveText *EtaRange9 =  new TPaveText(0.15, 0.7, 0.65, 0.87,"NDC");
  EtaRange9->AddText("5 < p_{T}^{track} < 10 GeV  2.0 < |#eta| < 2.5");
  EtaRange9->SetFillColor(0);
  EtaRange9->SetLineColor(0);
  EtaRange9->SetShadowColor(0);
  EtaRange9->Draw();
  //Other pt range
  Can_DeltaE_09_07->cd(10);
  TH1 *h_DeltaE09_10GeV_eta1 = (TH1F*)HistFile->Get("DeltaE_10GeV_eta1");
  h_DeltaE09_10GeV_eta1->SetStats(kFALSE);
  h_DeltaE09_10GeV_eta1->SetTitle(" ");
  h_DeltaE09_10GeV_eta1 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_10GeV_eta1 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_10GeV_eta1 ->SetLineColor(kRed);
  h_DeltaE09_10GeV_eta1 ->SetFillColor(kRed);
  h_DeltaE09_10GeV_eta1 ->SetFillStyle(3004);
  h_DeltaE09_10GeV_eta1->Scale(1.0/h_DeltaE09_10GeV_eta1->Integral());
  h_DeltaE09_10GeV_eta1->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_10GeV_eta1->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_10GeV_eta1->Draw();
  TH1 *h_DeltaE07_10GeV_eta1 = (TH1F*)HistFile->Get("DeltaE07_10GeV_eta1");
  h_DeltaE07_10GeV_eta1 ->SetLineColor(kBlue);
  h_DeltaE07_10GeV_eta1 ->SetFillColor(kBlue);
  h_DeltaE07_10GeV_eta1 ->SetFillStyle(3005);
  h_DeltaE07_10GeV_eta1->Scale(1.0/h_DeltaE07_10GeV_eta1->Integral());
  h_DeltaE07_10GeV_eta1->Draw("same");
  MyClus09->Draw();
  
  TPaveText *EtaRange10 = new TPaveText(0.15, 0.7, 0.55, 0.87,"NDC");
  EtaRange10->AddText("10 < p_{T}^{track} < 20 GeV   |#eta| < 1.0");
  EtaRange10->SetFillColor(0);
  EtaRange10->SetLineColor(0);
  EtaRange10->SetShadowColor(0);
  EtaRange10->Draw();
  
  Can_DeltaE_09_07->cd(11);
  TH1 *h_DeltaE09_10GeV_eta2 = (TH1F*)HistFile->Get("DeltaE_10GeV_eta2");
  h_DeltaE09_10GeV_eta2->SetStats(kFALSE);
  h_DeltaE09_10GeV_eta2->SetTitle(" ");
  h_DeltaE09_10GeV_eta2 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_10GeV_eta2 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_10GeV_eta2 ->SetLineColor(kRed);
  h_DeltaE09_10GeV_eta2 ->SetFillColor(kRed);
  h_DeltaE09_10GeV_eta2 ->SetFillStyle(3005);
  h_DeltaE09_10GeV_eta2->Scale(1.0/h_DeltaE09_10GeV_eta2->Integral());
  h_DeltaE09_10GeV_eta2->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_10GeV_eta2->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_10GeV_eta2->Draw();
  TH1 *h_DeltaE07_10GeV_eta2 = (TH1F*)HistFile->Get("DeltaE07_10GeV_eta2");
  h_DeltaE07_10GeV_eta2 ->SetLineColor(kBlue);
  h_DeltaE07_10GeV_eta2 ->SetFillColor(kBlue);
  h_DeltaE07_10GeV_eta2 ->SetFillStyle(3005);
  h_DeltaE07_10GeV_eta2->Scale(1.0/h_DeltaE07_10GeV_eta2->Integral());
  h_DeltaE07_10GeV_eta2->Draw("same");
  MyClus09->Draw();
  
  TPaveText *EtaRange11 = new TPaveText(0.15, 0.7, 0.65, 0.87,"NDC");
  EtaRange11->AddText("10 < p_{T}^{track} < 20 GeV  1.0 < |#eta| < 2.0");
  EtaRange11->SetFillColor(0);
  EtaRange11->SetLineColor(0);
  EtaRange11->SetShadowColor(0);
  EtaRange11->Draw();
  
  Can_DeltaE_09_07->cd(12);
  TH1 *h_DeltaE09_10GeV_eta25 = (TH1F*)HistFile->Get("DeltaE_10GeV_eta25");
  h_DeltaE09_10GeV_eta25->SetStats(kFALSE);
  h_DeltaE09_10GeV_eta25->SetTitle(" ");
  h_DeltaE09_10GeV_eta25 ->GetXaxis()->SetTitle("E_{clu}-E_{exp}/#Sigma(E^{exp})");
  h_DeltaE09_10GeV_eta25 ->GetYaxis()->SetTitle("Fraction of particles");
  h_DeltaE09_10GeV_eta25 ->SetLineColor(kRed);
  h_DeltaE09_10GeV_eta25 ->SetFillColor(kRed);
  h_DeltaE09_10GeV_eta25 ->SetFillStyle(3005);
  h_DeltaE09_10GeV_eta25->Scale(1.0/h_DeltaE09_10GeV_eta25->Integral());
  h_DeltaE09_10GeV_eta25->GetYaxis()->SetRangeUser(0,0.25);
  h_DeltaE09_10GeV_eta25->GetXaxis()->SetTitleSize(0.05);
  h_DeltaE09_10GeV_eta25->Draw();
  TH1 *h_DeltaE07_10GeV_eta25 = (TH1F*)HistFile->Get("DeltaE07_10GeV_eta25");
  h_DeltaE07_10GeV_eta25 ->SetLineColor(kBlue);
  h_DeltaE07_10GeV_eta25 ->SetFillColor(kBlue);
  h_DeltaE07_10GeV_eta25 ->SetFillStyle(3005);
  h_DeltaE07_10GeV_eta25->Scale(1.0/h_DeltaE07_10GeV_eta25->Integral());
  h_DeltaE07_10GeV_eta25->Draw("same");
  MyClus09->Draw();
  
  TPaveText *EtaRange12 = new TPaveText(0.15, 0.7, 0.65, 0.87,"NDC");
  EtaRange12->AddText("10 < p_{T}^{track} < 20 GeV  2.0 < |#eta| < 2.5");
  EtaRange12->SetFillColor(0);
  EtaRange12->SetLineColor(0);
  EtaRange12->SetShadowColor(0);
  EtaRange12->Draw();
  return;
}












//////////////////////////////
void InitStyle()
{
  gROOT->ForceStyle();

  //gStyle->SetHistLineWidth(2);
  //gStyle->SetHistFillStyle(3003);
  //gStyle->SetHistFillColor(color_histos);
  //gStyle->SetHistLineColor(gStyle->GetHistFillColor());

  gStyle->SetLabelSize(0.035,"X");
  gStyle->SetLabelSize(0.035,"Y");
  gStyle->SetLabelSize(0.03,"Z");

  gStyle->SetTitleOffset(1.,"X");
  gStyle->SetTitleOffset(1.2,"Y");

  // gStyle->SetStatColor(kRed);
  gStyle->SetStatFormat("5.3g");

  gStyle->SetOptStat("emr");

  return;
}



void PrintHistosCode()
//////////////////////////////
{
  cout << " ** PrintHistosCode ** # elements " << HistogramList.size() << endl;

  // print histogram color code
  float Legend_lx = 0.60;
  float Legend_ux = 0.88;
  float Legend_ly = 0.66;
  float Legend_uy = 0.76;

  TLegend *Rotllo = new TLegend(Legend_lx, Legend_ly, Legend_ux, Legend_uy);
  Rotllo->SetFillStyle(3005);
  Rotllo->SetFillColor(10);
  char QueDir[80];
  for (int i=0; i < HistogramList.size(); i++) {
    sprintf(QueDir,"hola %d", i);
    if (i==0) sprintf(QueDir,"Truth #pi^{+}");
    if (i==1) sprintf(QueDir,"Truth #pi^{+} (p_{T}>400MeV &  |#eta<=2.5|");
    if (i==2) sprintf(QueDir,"PFO");
    Rotllo->AddEntry(HistogramList.at(i), QueDir);
  }
  Rotllo->Draw();

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


void DrawLatex( const char* label){
  TLatex *t = new TLatex();
  t->SetNDC();
  t->SetTextAlign(22);
  t->SetTextFont(42);
  t->SetTextSizePixels(20);
  t->DrawLatex(0.25,0.85,label);
  
  return;
} 

float FitGaussHist( TH1D* hist, bool GetMean, bool GetSigma){

  float Peak = hist->GetBinCenter(hist->GetMaximumBin()); 
  float mean = hist->GetMean();
  float RMS = hist->GetRMS();
  TF1 *fitg = new TF1("fitg","gaus",mean-1.5*RMS,mean+1.5*RMS);
  hist->Fit(fitg,"R+");
  hist->DrawCopy();

  double p2 = fitg->GetParameter(2); double e_p2 = fitg->GetParError(2);   
  double p1 = fitg->GetParameter(1); double e_p1 = fitg->GetParError(1);    
  float DeltaPar =  fabs(p1-Peak); 
  
  float f_p2 = p2; 
  float f_p1 = p1; 
  int IterNum = 0;
  
  while (0.1 < DeltaPar && IterNum < 10 ){
    cout<<"Initial DeltaPar"<<DeltaPar<<endl; 
    f_p2 = p2;
    f_p1 = p1;
    cout<<"Parameters p1="<<p1<<"  p2="<<p2<<endl;
    TF1 *fitg2 = new TF1("fitg2","gaus",p1-2*p2,p1+2*p2);
    hist->Fit(fitg2,"R");
    hist->DrawCopy();
    float new_p2 = fitg2->GetParameter(2);   
    float new_p1 = fitg2->GetParameter(1);
    float new_ep2 = fitg2->GetParError(2);   
    float new_ep1 = fitg2->GetParError(1);
    DeltaPar = fabs(p1-new_p1);
    cout<<"Final DeltaPar"<<DeltaPar<<endl; 
    cout<<"new_p1="<<new_p1<<"  new_p2="<<new_p2<<endl;
    p2 = new_p2;
    p1 = new_p1;
    e_p2 = new_ep2;
    e_p1 = new_ep1;
    IterNum++;
  } 

  TPaveText *pav = new TPaveText(0.05, 0.80, 0.50,0.95,"NDC");
  pav->SetBorderSize(1);  pav->SetFillColor(0);
  char texto1[100]; 
  sprintf(texto1, "  #mu(E^{true}) = %2.3f #pm %2.3f", f_p1,e_p1);
  pav->AddText(texto1);
  sprintf(texto1, " #sigma(E^{true}) = %2.3f #pm %2.3f", f_p2,e_p2);
  pav->AddText(texto1);
  pav->Draw();
   
  if(GetMean) return p1;
  if (GetSigma) return p2;
}
  
  




void colorins()
{
  if (false) {
    // escala de ferro
    const Int_t NRGBs = 4;
    const Int_t NCont = 50; 
    
    Double_t stops[NRGBs] = { 0.20, 0.40, 0.80, 1.00};
    Double_t red[NRGBs]   = { 0.20, 1.00, 0.95, 0.95};
    Double_t green[NRGBs] = { 0.00, 0.10, 0.95, 0.95};
    Double_t blue[NRGBs]  = { 0.00, 0.10, 0.00, 1.00};
    
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }

  if (true) {
    const Int_t NRGBs = 5;
    const Int_t NCont = 99; // 255 originaly but the DrawPannel complaints
    
    Double_t stops[NRGBs] = { 0.20, 0.40, 0.60, 0.80, 1.00 };
    Double_t red[NRGBs]   = { 0.95, 0.40, 0.50, 1.00, 0.00 };
    Double_t green[NRGBs] = { 0.95, 0.40, 0.90, 0.10, 0.00 };
    Double_t blue[NRGBs]  = { 0.95, 0.95, 0.20, 0.02, 0.00 };
    
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }

  if (false) {
    // escala de rojos
    const Int_t NRGBs = 2;
    const Int_t NCont = 25; 
    
    Double_t stops[NRGBs] = { 0.10, 1.00};
    Double_t red[NRGBs]   = { 1.00, 1.00};
    Double_t green[NRGBs] = { 1.00, 0.05 };
    Double_t blue[NRGBs]  = { 1.00, 0.05 };
    
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }
  if (false) {
    // escala de blavets
    const Int_t NRGBs = 2;
    const Int_t NCont = 50; 
    
    Double_t stops[NRGBs] = { 0.10, 1.00};
    Double_t red[NRGBs]   = { 0.05, 0.10};
    Double_t green[NRGBs] = { 0.00, 0.90 };
    Double_t blue[NRGBs]  = { 0.05, 1.00 };
    
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
  }


  return;
}
