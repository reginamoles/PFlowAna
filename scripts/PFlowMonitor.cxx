#include "PFlowMonitor.h"
#include <TLegend.h>

void PFlowMonitor()
{

  std::cout<<"=================="<<std::endl;
  std::cout<<"PFlowAna monitor "<<std::endl;
  std::cout<<"=================="<<std::endl;
  
  //=================================
  // Chose here your config options:
  //=================================
  m_1to2matching = true;
  m_UseNarrowPtRange = true;
  m_UseNarrowEtaRange = true;
  
  /* WIP  Narrow or Wide eta and pt range */
  if (m_UseNarrowPtRange) _ptRange= {0, 2, 5, 10, 20};
  else _ptRange = {0, 2, 5, 10, 20, 40, 60, 80, 100, 150, 200, 500, 1000};
  if (m_UseNarrowEtaRange) _etaRange= {0, 1, 2, 2.5};
  else _etaRange= {0.0, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.5};
  
  HistFile = new TFile ("/afs/cern.ch/work/z/zhangr/eflowRec/PFlowAnaPackage/MyDir/hist-Run.root");

  setStyle();
  Efficiency();
//  printPS();
}

////
//////////////////////////////
void Efficiency()
{

  TCanvas* Can_Efficiency = new TCanvas("Efficiency", "Efficiency of 1st matched cluster", 900, 800);
  Can_Efficiency->Divide(2, 2);
  TH1F* h_Eff1[5];

  int tcolor[5] = {4, 2, 3, 8, 6};
  std::string catagory = "EffMatch1";


  int etabin(0);

  Can_Efficiency->cd(etabin+1);
  double xpos(0.25), ypos(0.84);
  TLegend* Legend = new TLegend(xpos, ypos - 0.08 * 3, xpos + 0.3, ypos);
  Legend->SetFillStyle(0);
  Legend->SetBorderSize(0);
  Legend->SetTextFont(43);
  Legend->SetTextSize(20);

  for (int ipt = 0; ipt < _ptRange.size(); ++ipt) {
    std::pair<std::string, std::string> names = histName(ipt, 0, catagory, "", _ptRange, _etaRange);

    h_Eff1[ipt] = (TH1F*) HistFile->Get(names.first.c_str());
    if (!h_Eff1[ipt]) {
      std::cerr << "[ERROR]\t Histogram " << names.first << " not exist!" << std::endl;
    }
    h_Eff1[ipt]->SetLineWidth(2);
    h_Eff1[ipt]->Print();
    h_Eff1[ipt]->SetStats(kFALSE);
    h_Eff1[ipt]->SetLineColor(tcolor[ipt]);
    h_Eff1[ipt]->GetXaxis()->SetTitle("#varepsilon_{1st cluster}");
    h_Eff1[ipt]->GetYaxis()->SetTitle("Fraction of particles");
    h_Eff1[ipt]->SetTitle(Form("%1.1f < |#eta_{EM2}| < %1.1f", _etaRange[etabin], _etaRange[etabin+1]));
    h_Eff1[ipt]->Scale(1./h_Eff1[ipt]->Integral());
    h_Eff1[ipt]->GetYaxis()->SetRangeUser(0, 1.2);
    if (ipt == 0) {
      h_Eff1[ipt]->Draw("hist");
    } else {
      h_Eff1[ipt]->Draw("samehist");
    }
    std::string lable = names.second;
    Legend->AddEntry(h_Eff1[ipt], lable.c_str(), "le");
  }
  Legend->Draw();
  Can_Efficiency->SaveAs(Form("plots/%s_%d.eps", catagory.c_str(), etabin));


  return;







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






////////////////////////////
void printPS() {
  
  cout << "\nStoring the plots in a ps file..." << endl;
  system("mkdir -vp plots/");
  Char_t filename[] = "plots/Subtraction_plots.eps";
  Char_t command[] = "";
  Char_t name[] = "";
  
  TCanvas c0;
  sprintf(command,"%s",filename);
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
  
  gStyle->SetTextFont(42);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetPadLeftMargin(0.17);
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetPadRightMargin(0.10);
  gStyle->SetPadTopMargin(0.15);
  // axis' scale
  gStyle->SetLabelSize(0.08, "xyz");
  gStyle->SetLabelOffset(0.01, "xyz");
  gStyle->SetNdivisions(505, "xyz");
  gStyle->SetLineWidth(2);

  // axis' name
  gStyle->SetTitleColor(1, "xyz");
  gStyle->SetTitleSize(0.08, "xyz");
  gStyle->SetTitleOffset(1.0, "xyz");

  gStyle->SetStripDecimals(false);
  TGaxis::SetMaxDigits(4);
  gStyle->SetPalette(1);

  gStyle->SetPadGridX(true);
  gStyle->SetPadGridY(true);
  gStyle->SetGridColor(10);

  gROOT->ForceStyle();
  return; 
}


std::string to_string_with_precision(float value, int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << value;
    return out.str();
}

std::pair<std::string, std::string> histName(unsigned i_pt, unsigned i_eta, const std::string& name, const std::string& matchScheme, std::vector<float>& PtRange,
                                   std::vector<float>& EtaRange) {

  std::string complete_name = "WrongName", legend_name = "WrongName";
  if (i_pt != PtRange.size()-1 && i_eta != EtaRange.size()-1) {
    complete_name = (name + matchScheme + "_" + std::to_string((int) (PtRange.at(i_pt))) + "_" + std::to_string((int) (PtRange.at(i_pt+1))) + "GeV__eta"
                     + std::to_string((int) ((10 * EtaRange.at(i_eta)))) + "_" + std::to_string((int) ((10 * EtaRange.at(i_eta+1))))).c_str();
    legend_name = std::to_string((int) (PtRange.at(i_pt))) + "-" + std::to_string((int) (PtRange.at(i_pt+1))) + "GeV";
  }
  else {
    complete_name = (name + matchScheme + "_" + std::to_string((int) (PtRange.at(i_pt))) + "GeV__eta"
                      + std::to_string((int) ((10 * EtaRange.at(i_eta)))) ).c_str();
    legend_name = ">" + std::to_string((int) (PtRange.at(i_pt))) + "GeV";
  }

  return std::make_pair(complete_name, legend_name);
}
