#include "PFlowMonitor.h"
#include <TLegend.h>
#include <sstream>
#include <string>

PFlowMonitor::PFlowMonitor() {
  m_1to2matching = true;
  m_UseNarrowPtRange = true;
  m_UseNarrowEtaRange = true;
}

PFlowMonitor::~PFlowMonitor() {};

void PFlowMonitor::run(char* inputs)
{

  std::cout<<"=================="<<std::endl;
  std::cout<<"PFlowAna monitor "<<std::endl;
  std::cout<<"=================="<<std::endl;
  
  //=================================
  // Chose here your config options:
  //=================================
  
  /* WIP  Narrow or Wide eta and pt range */
  if (m_UseNarrowPtRange) {
    m_ptRange.push_back(0);
    m_ptRange.push_back(2);
    m_ptRange.push_back(5);
    m_ptRange.push_back(10);
    m_ptRange.push_back(20);
  }
  else {
    m_ptRange.push_back(0);
    m_ptRange.push_back(2);
    m_ptRange.push_back(5);
    m_ptRange.push_back(10);
    m_ptRange.push_back(20);
    m_ptRange.push_back(40);
    m_ptRange.push_back(60);
    m_ptRange.push_back(80);
    m_ptRange.push_back(100);
    m_ptRange.push_back(150);
    m_ptRange.push_back(200);
    m_ptRange.push_back(500);
    m_ptRange.push_back(1000);
  }
  if (m_UseNarrowEtaRange) {
    m_etaRange.push_back(0);
    m_etaRange.push_back(1);
    m_etaRange.push_back(2);
    m_etaRange.push_back(2.5);
  } else {
    m_etaRange.push_back(0.0);
    m_etaRange.push_back(0.6);
    m_etaRange.push_back(0.8);
    m_etaRange.push_back(1.0);
    m_etaRange.push_back(1.2);
    m_etaRange.push_back(1.4);
    m_etaRange.push_back(1.6);
    m_etaRange.push_back(2.0);
    m_etaRange.push_back(2.5);
  }

  ifstream inputfiles(inputs);
  std::string inputfile;
  while (!inputfiles.eof()) {
    inputfiles >> inputfile;
    /* Enable to comment in configVariable */
    if ((inputfile.compare(0, 1, "#") == 0)) continue;
    if (!inputfiles.good()) break;
    HistFile.push_back(new TFile(inputfile.c_str()));
  }

  setStyle();
  Efficiency();
}

////
//////////////////////////////
void PFlowMonitor::Efficiency()
{
  int tcolor[5] = {4, 2, 8, 6, 28};

  std::vector<std::string> catagory;
  catagory.push_back("EffMatch1");
  catagory.push_back("EffMatchboth");
  catagory.push_back("PurMatch1");
  catagory.push_back("PurMatch2");
  catagory.push_back("SubtractStatus");
  std::vector<std::string> xTitle;
  xTitle.push_back("#varepsilon_{1st cluster}");
  xTitle.push_back("#varepsilon_{both clusters}");
  xTitle.push_back("P_{1st cluster}");
  xTitle.push_back("P_{2nd cluster}");
  xTitle.push_back("Stage");

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

        for(unsigned int ifile = 0; ifile < HistFile.size(); ++ifile) {
          if (ifile == 0) {
            h_pTs[ipt] = (TH1F*) HistFile[ifile]->Get(names.first.c_str());
          } else {
            h_pTs[ipt]->Add((TH1F*) HistFile[ifile]->Get(names.first.c_str()));
          }
        }
        if (!h_pTs[ipt]) {
          std::cerr << "[ERROR]\t Histogram " << names.first << " not exist!" << std::endl;
        }
        h_pTs[ipt]->SetLineWidth(2);
        h_pTs[ipt]->Print();
        h_pTs[ipt]->SetStats(kFALSE);
        h_pTs[ipt]->SetLineColor(tcolor[ipt]);
        h_pTs[ipt]->GetXaxis()->SetTitle(xTitle[icat].c_str());
        h_pTs[ipt]->GetYaxis()->SetTitle("Fraction of particles");
        if(ieta == m_etaRange.size()-1) {
          h_pTs[ipt]->SetTitle(Form("%1.1f < |#eta_{EM2}|", m_etaRange[ieta]));
        } else {
          h_pTs[ipt]->SetTitle(Form("%1.1f < |#eta_{EM2}| < %1.1f", m_etaRange[ieta], m_etaRange[ieta + 1]));
        }
        double entries = h_pTs[ipt]->Integral();
        std::cout<<ipt<<" "<<entries<<std::endl;
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
          h_pTs[ipt]->Draw("hist");
        } else {
          h_pTs[ipt]->Draw("samehist");
        }
      }
      if (!first) {
        std::cout << "[WARNING]\t Histograms are empty in eta bin " << ieta << std::endl;
        continue;
      }
      Legend->Draw();
      h_pTs[0]->Draw("axissame");
      system("mkdir -vp plots/");
      Can_Efficiency->GetPad(ieta + 1)->SaveAs(Form("plots/%s_%d.eps", catagory[icat].c_str(), ieta));
    }
  }
  return;

}



/////////////////////
void PFlowMonitor::setStyle() {
  
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



std::pair<std::string, std::string> PFlowMonitor::histName(unsigned i_pt, unsigned i_eta, const std::string& name, const std::string& matchScheme, std::vector<float>& PtRange,
                                   std::vector<float>& EtaRange) {

  std::string complete_name = "WrongName", legend_name = "WrongName";
  if (i_pt != PtRange.size()-1 && i_eta != EtaRange.size()-1) {
    complete_name = (name + matchScheme + "_" + Int_to_String((int) (PtRange.at(i_pt))) + "_" + Int_to_String((int) (PtRange.at(i_pt+1))) + "GeV__eta"
                     + Int_to_String((int) ((10 * EtaRange.at(i_eta)))) + "_" + Int_to_String((int) ((10 * EtaRange.at(i_eta+1))))).c_str();
    legend_name = Int_to_String((int) (PtRange.at(i_pt))) + "-" + Int_to_String((int) (PtRange.at(i_pt+1))) + "GeV";
  }
  else {
    complete_name = (name + matchScheme + "_" + Int_to_String((int) (PtRange.at(i_pt))) + "GeV__eta"
                      + Int_to_String((int) ((10 * EtaRange.at(i_eta)))) ).c_str();
    legend_name = ">" + Int_to_String((int) (PtRange.at(i_pt))) + "GeV";
  }

  return std::make_pair(complete_name, legend_name);
}

std::string PFlowMonitor::Int_to_String(int n) {
  std::ostringstream stream;
  stream << n;
  return stream.str();
}
