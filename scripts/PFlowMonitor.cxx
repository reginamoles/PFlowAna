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
  if (m_UseNarrowPtRange) m_ptRange= {0, 2, 5, 10, 20};
  else m_ptRange = {0, 2, 5, 10, 20, 40, 60, 80, 100, 150, 200, 500, 1000};
  if (m_UseNarrowEtaRange) m_etaRange= {0, 1, 2, 2.5};
  else m_etaRange= {0.0, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.5};
  
  HistFile = new TFile ("/afs/cern.ch/work/z/zhangr/eflowRec/PFlowAnaPackage/MyDir/hist-Run.root");

  setStyle();
  Efficiency();
//  printPS();
}

////
//////////////////////////////
void Efficiency()
{
  int tcolor[5] = {4, 2, 8, 6, 28};

  std::vector<std::string> catagory = {"EffMatch1", "EffMatchboth"};
  std::vector<std::string> xTitle = {"#varepsilon_{1st cluster}", "#varepsilon_{both clusters}"};

  for (int icat = 0; icat < catagory.size(); ++icat) {
    TCanvas* Can_Efficiency = new TCanvas(catagory[icat].c_str(), catagory[icat].c_str(), 900, 800);
    Can_Efficiency->Divide(2, 2);
    TH1F* h_pTs[5];

    for (int ieta = 0; ieta < m_etaRange.size(); ++ieta) {
      Can_Efficiency->cd(ieta + 1);
      double xpos(0.25), ypos(0.88);
      TLegend* Legend = new TLegend(xpos, ypos - 0.08 * 3, xpos + 0.3, ypos);
      Legend->SetFillStyle(0);
      Legend->SetBorderSize(0);
      Legend->SetTextFont(43);
      Legend->SetTextSize(20);

      int first(0);
      for (int ipt = 0; ipt < m_ptRange.size(); ++ipt) {
        std::pair<std::string, std::string> names = histName(ipt, ieta, catagory[icat], "", m_ptRange, m_etaRange);

        h_pTs[ipt] = (TH1F*) HistFile->Get(names.first.c_str());
        if (!h_pTs[ipt]) {
          std::cerr << "[ERROR]\t Histogram " << names.first << " not exist!" << std::endl;
        }
        h_pTs[ipt]->SetLineWidth(2);
        h_pTs[ipt]->Print();
        h_pTs[ipt]->SetStats(kFALSE);
        h_pTs[ipt]->SetLineColor(tcolor[ipt]);
        h_pTs[ipt]->GetXaxis()->SetTitle(xTitle[icat].c_str());
        h_pTs[ipt]->GetYaxis()->SetTitle("Fraction of particles");
        h_pTs[ipt]->SetTitle(Form("%1.1f < |#eta_{EM2}| < %1.1f", m_etaRange[ieta], m_etaRange[ieta + 1]));
        double entries = h_pTs[ipt]->Integral();
        if (entries > 0) {
          h_pTs[ipt]->Scale(1. / entries);
          first++;
        }
        h_pTs[ipt]->GetYaxis()->SetRangeUser(0, 1.2);
        std::string lable = names.second;
        Legend->AddEntry(h_pTs[ipt], lable.c_str(), "le");
        if (entries == 0) {
          continue;
        }
        if (first == 1) {
          h_pTs[ipt]->Draw("histE");
        } else {
          h_pTs[ipt]->Draw("samehistE");
        }
      }
      if (!first) {
        std::cout << "[WARNING]\t Histograms are empty in eta bin " << ieta << std::endl;
        continue;
      }
      Legend->Draw();
      h_pTs[0]->Draw("axissame");
      Can_Efficiency->GetPad(ieta + 1)->SaveAs(Form("plots/%s_%d.eps", catagory[icat].c_str(), ieta));
    }
  }
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
  gStyle->SetPadTopMargin(0.1);
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
