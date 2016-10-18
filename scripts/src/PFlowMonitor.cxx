#include "PFlowMonitor.h"
#include <TLegend.h>
#include <sstream>
#include <string>

PFlowMonitor::PFlowMonitor() {
  m_1to2matching = true;
  m_UseNarrowPtRange = true;
  m_UseNarrowEtaRange = true;
  m_debug = false;
}

PFlowMonitor::~PFlowMonitor() {};

void PFlowMonitor::run(char* inputs, char* outfolder)
{

  std::cout<<"=================="<<std::endl;
  std::cout<<"PFlowAna monitor "<<std::endl;
  std::cout<<"=================="<<std::endl;
  
  //=================================
  // Chose here your config options:
  //=================================
  
  /* WIP  Narrow or Wide eta and pt range */
  if (m_UseNarrowPtRange) {
//    m_ptRange.push_back(0);
//    m_ptRange.push_back(100);
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
//    m_etaRange.push_back(0);
//    m_etaRange.push_back(100);
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

  std::ifstream inputfiles(inputs);
  std::string inputfile;
  while (!inputfiles.eof()) {
    inputfiles >> inputfile;
    /* Enable to comment in configVariable */
    if ((inputfile.compare(0, 1, "#") == 0)) continue;
    if (!inputfiles.good()) break;
    HistFile.push_back(new TFile(inputfile.c_str()));
  }

  system(Form("mkdir -vp plots/%s/", outfolder));
  setStyle();
//  Plot(Form("plots/%s/", outfolder));
//  Efficiency(Form("plots/%s/", outfolder));
//  eflowdRp(Form("plots/%s/", outfolder), 0);
//  eflowdRp(Form("plots/%s/", outfolder), 1);

  AverageEfficiencyPerPt(Form("plots/%s/", outfolder));

}

////
//////////////////////////////
void PFlowMonitor::Plot(const char* outfolder) {
  bool c_showNumber(!false);
  bool c_overlay(true);
  int tcolor[5] = {4, 2, 8, 6, 28};

  std::vector<std::string> catagory;
  catagory.push_back("EffMatch1");
  catagory.push_back("EffMatchboth");
  catagory.push_back("PurMatch1");
  catagory.push_back("PurMatch2");
  catagory.push_back("dRp_leading");
  catagory.push_back("dRp_1st");
  catagory.push_back("dRp_2nd");
  catagory.push_back("SubtractStatus");
  catagory.push_back("NClus_09");

  catagory.push_back("EffLeading");

  std::vector<std::string> xTitle;
  xTitle.push_back("#varepsilon_{1st cluster}");
  xTitle.push_back("#varepsilon_{both clusters}");
  xTitle.push_back("P_{1st cluster}");
  xTitle.push_back("P_{2nd cluster}");
  xTitle.push_back("#Delta R'_{leading cluster}^2");
  xTitle.push_back("#Delta R'_{1st cluster}^2");
  xTitle.push_back("#Delta R'_{2nd cluster}^2");
  xTitle.push_back("Stage");
  xTitle.push_back("N_{cluster}(#Sigma E^{true}>90%)");
  xTitle.push_back("#varepsilon_{leading cluster}");

  for (unsigned int icat = 0; icat < (m_debug ? 1 : catagory.size()); ++icat) {
    if (icat > 8) c_overlay = false;
    TCanvas* Can_Efficiency = new TCanvas(catagory[icat].c_str(), catagory[icat].c_str(), 900, 800);
    Can_Efficiency->Divide(2, 2);
    TH1F* h_pTs[5];

    for (unsigned int ieta = 0; ieta < (m_debug ? 1 : m_etaRange.size()); ++ieta) {
      Can_Efficiency->cd(ieta + 1);
      double xpos(0.25), ypos(0.88);
      TLegend* Legend = new TLegend(xpos, ypos - 0.08 * 3, xpos + 0.3, ypos);
      Legend->SetFillStyle(0);
      Legend->SetBorderSize(0);
      Legend->SetTextFont(43);
      Legend->SetTextSize(20);

      int first(0);
      for (unsigned int ipt = 0; ipt < m_ptRange.size(); ++ipt) {
        std::pair<std::string, std::string> names = histName(ipt, ieta, catagory[icat], "", m_ptRange, m_etaRange);

        for(unsigned int ifile = 0; ifile < HistFile.size(); ++ifile) {
          if (m_debug) std::cout<<HistFile[ifile]->GetName()<<std::endl;
          if (ifile == 0) {
            h_pTs[ipt] = (TH1F*) HistFile[ifile]->Get(names.first.c_str());
          } else {
            h_pTs[ipt]->Add((TH1F*) HistFile[ifile]->Get(names.first.c_str()));
            std::cout<<__LINE__<<" Add: "<<names.first<<std::endl;
          }
        }
        if (!h_pTs[ipt]) {
          std::cerr << "[ERROR]\t Histogram " << names.first << " not exist!" << std::endl;
        }
        h_pTs[ipt]->SetLineWidth(2);
//        h_pTs[ipt]->Print();
        h_pTs[ipt]->SetStats(kFALSE);
        h_pTs[ipt]->SetLineColor(tcolor[ipt]);
        h_pTs[ipt]->GetXaxis()->SetTitle(xTitle[icat].c_str());
        if (c_overlay) {
          h_pTs[ipt]->GetYaxis()->SetTitle("Fraction of particles");
        } else {
          h_pTs[ipt]->GetYaxis()->SetTitle("Number of particles");
        }
        if(ieta == m_etaRange.size()-1) {
          h_pTs[ipt]->SetTitle(Form("%1.1f < |#eta_{EM2}|", m_etaRange[ieta]));
        } else {
          h_pTs[ipt]->SetTitle(Form("%1.1f < |#eta_{EM2}| < %1.1f", m_etaRange[ieta], m_etaRange[ieta + 1]));
        }
        double entries = h_pTs[ipt]->Integral();
        if (c_overlay && entries > 0) {
          h_pTs[ipt]->Scale(1. / entries);
          first++;
        }
        double maxY = h_pTs[ipt]->GetBinCenter(h_pTs[ipt]->GetMaximumBin());
        if (c_overlay) {
          h_pTs[ipt]->GetYaxis()->SetRangeUser(0, (maxY > 0.8 ? maxY * 1.2 : 0.8));
        } else {
//          h_pTs[ipt]->GetYaxis()->SetRangeUser(0, maxY);
        }
        std::string lable;
        if (c_showNumber) {
          lable = Form("%s (%.0f)", names.second.c_str(), entries);
        } else {
          lable = Form("%s", names.second.c_str());
        }

        Legend->AddEntry(h_pTs[ipt], lable.c_str(), "le");
        if (entries == 0) {
          if (!c_overlay) {
            Legend->Clear();
            Can_Efficiency->GetPad(ieta + 1)->Update();
          }
          continue;
        }
        if (c_overlay) {
          if (first == 1) {
            h_pTs[ipt]->Draw("hist");
          } else {
            h_pTs[ipt]->Draw("samehist");
          }
        } else {
          Can_Efficiency->cd(ieta + 1);
          h_pTs[ipt]->SetLineColor(1);
          h_pTs[ipt]->Draw("hist");
          Legend->Draw();
          Can_Efficiency->GetPad(ieta + 1)->SaveAs(Form("%s%s_%d_%d.eps", outfolder, catagory[icat].c_str(), ieta, ipt));
          Legend->Clear();
          Can_Efficiency->GetPad(ieta + 1)->Update();
        }
      }
      if (!first) {
        std::cout << "[WARNING]\t Histograms are empty in eta bin " << ieta << std::endl;
        continue;
      }

      Legend->Draw();
      h_pTs[0]->Draw("axissame");
      if (c_overlay) {
        Can_Efficiency->GetPad(ieta + 1)->SaveAs(Form("%s%s_%d.eps", outfolder, catagory[icat].c_str(), ieta));
      }
    }
  }
  return;

}

void PFlowMonitor::Efficiency(const char* outfolder) {
  bool c_showNumber(!false);
  int tcolor[4] = {1, 4, 8, 2};

  std::vector<std::string> catagory;
  catagory.push_back("EffClusterboth_total");
  catagory.push_back("EffClusterboth_CLS1");
  catagory.push_back("EffClusterboth_CLS2");
  catagory.push_back("EffClusterboth_RSS");

  std::vector<std::string> display;
  display.push_back("Total");
  display.push_back("Matched 1 cluster");
  display.push_back("Matched 2 clusters");
  display.push_back("Recover split shower");

  std::string xTitle = "#varepsilon_{both clusters}";

  for (unsigned int ieta = 0; ieta < (m_debug ? 1 : m_etaRange.size()); ++ieta) {

    for (unsigned int ipt = 0; ipt < m_ptRange.size(); ++ipt) {
      TCanvas* Can_Efficiency = new TCanvas("Efficiency", "Efficiency", 450, 400);
      TH1F* h_cats[4];

      double xpos(0.25), ypos(0.88);
      TLegend* Legend = new TLegend(xpos, ypos - 0.08 * 3, xpos + 0.3, ypos);
      Legend->SetFillStyle(0);
      Legend->SetBorderSize(0);
      Legend->SetTextFont(43);
      Legend->SetTextSize(20);

      for (unsigned int icat = 0; icat < (m_debug ? 1 : catagory.size()); ++icat) {

        std::pair<std::string, std::string> names = histName(ipt, ieta, catagory[icat], "", m_ptRange, m_etaRange);

        for (unsigned int ifile = 0; ifile < HistFile.size(); ++ifile) {
          if (m_debug) std::cout << HistFile[ifile]->GetName() << std::endl;
          if (ifile == 0) {
            h_cats[icat] = (TH1F*) HistFile[ifile]->Get(names.first.c_str());
          } else {
            h_cats[icat]->Add((TH1F*) HistFile[ifile]->Get(names.first.c_str()));
            std::cout<<__LINE__<<" Add: "<<names.first<<std::endl;
          }
        }

        if (!h_cats[icat]) {
          std::cerr << "[ERROR]\t Histogram " << names.first << " not exist!" << std::endl;
        }

        h_cats[icat]->SetLineWidth(2);
//        h_cats[icat]->Print();
        h_cats[icat]->SetStats(kFALSE);
        h_cats[icat]->GetXaxis()->SetTitle(xTitle.c_str());
        h_cats[icat]->GetYaxis()->SetTitle("Number of particles");
        if (ieta == m_etaRange.size() - 1 && ipt == m_ptRange.size()) {
          h_cats[icat]->SetTitle(Form("%1.1f < |#eta_{EM2}| && %.0f GeV < p_{T}", m_etaRange[ieta], m_ptRange[ipt]));
        } else if (ieta == m_etaRange.size() - 1 && ipt != m_ptRange.size()) {
          h_cats[icat]->SetTitle(Form("%1.1f < |#eta_{EM2}| && %.0f < p_{T} < %.0f GeV", m_etaRange[ieta], m_ptRange[ipt], m_ptRange[ipt + 1]));
        } else if (ieta != m_etaRange.size() - 1 && ipt == m_ptRange.size() - 1) {
          h_cats[icat]->SetTitle(Form("%1.1f < |#eta_{EM2}| < %1.1f && %.0f GeV < p_{T}", m_etaRange[ieta], m_etaRange[ieta + 1], m_ptRange[ipt]));
        } else {
          h_cats[icat]->SetTitle(Form("%1.1f < |#eta_{EM2}| < %1.1f && %.0f < p_{T} < %.0f GeV", m_etaRange[ieta], m_etaRange[ieta + 1], m_ptRange[ipt], m_ptRange[ipt + 1]));
        }
        h_cats[icat]->SetLineColor(tcolor[icat]);
        std::string lable;
        double entries = h_cats[icat]->Integral();
        if (c_showNumber) {
          lable = Form("%s (%.0f)", display[icat].c_str(), entries);
        } else {
          lable = Form("%s", display[icat].c_str());
        }

        Legend->AddEntry(h_cats[icat], lable.c_str(), "l");
        if (entries == 0) {
          continue;
        }
        if (icat == 0) {
          h_cats[icat]->Draw("hist");
        } else {
          h_cats[icat]->Draw("samehist");
        }

        Legend->Draw();
      }
      h_cats[0]->Draw("axissame");
      Can_Efficiency->SaveAs(Form("%s%s_%d_%d.eps", outfolder, catagory[0].c_str(), ieta, ipt));
      delete Can_Efficiency;
    }
  }
  return;
}

void PFlowMonitor::eflowdRp(const char* outfolder, const int mode) {
  bool c_showNumber(!false);
  int tcolor[3] = {1, 4, 2};

  std::vector<std::string> catagory; catagory.clear();
  catagory.push_back("eflowdR1");
  if (mode == 0) {
      catagory.push_back("eflowdR1_CLS");
      catagory.push_back("eflowdR1_RSS");
  } else if (mode == 1) {
      catagory.push_back("eflowdR1_correct");
      catagory.push_back("eflowdR1_wrong");
  }

  std::vector<std::string> display;
  display.push_back("Total");
  if (mode == 0) {
      display.push_back("dR' matching");
      display.push_back("Recover split shower");
  } else if (mode == 1) {
      display.push_back("correct matching");
      display.push_back("wrong matching");
  }

  std::string xTitle = "eflowRec dR'(1st cluster)^2";

  for (unsigned int ieta = 0; ieta < (m_debug ? 1 : m_etaRange.size()); ++ieta) {

    for (unsigned int ipt = 0; ipt < m_ptRange.size(); ++ipt) {
      TCanvas* Can_Efficiency = new TCanvas("dR1", "dR1", 450, 400);
      TH1F* h_cats[4];

      double xpos(0.25), ypos(0.88);
      TLegend* Legend = new TLegend(xpos, ypos - 0.08 * 3, xpos + 0.3, ypos);
      Legend->SetFillStyle(0);
      Legend->SetBorderSize(0);
      Legend->SetTextFont(43);
      Legend->SetTextSize(20);

      for (unsigned int icat = 0; icat < (m_debug ? 1 : catagory.size()); ++icat) {

        std::pair<std::string, std::string> names = histName(ipt, ieta, catagory[icat], "", m_ptRange, m_etaRange);

        for (unsigned int ifile = 0; ifile < HistFile.size(); ++ifile) {
          if (m_debug) std::cout << HistFile[ifile]->GetName() << std::endl;
          if (ifile == 0) {
            h_cats[icat] = (TH1F*) HistFile[ifile]->Get(names.first.c_str());
            std::cout<<__LINE__<<" Get: "<<names.first<<" "<<h_cats[icat]->Integral()<<std::endl;
          } else {
            h_cats[icat]->Add((TH1F*) HistFile[ifile]->Get(names.first.c_str()));
            std::cout<<__LINE__<<" Add: "<<names.first<<" "<<h_cats[icat]->Integral()<<std::endl;
          }
        }

        if (!h_cats[icat]) {
          std::cerr << "[ERROR]\t Histogram " << names.first << " not exist!" << std::endl;
        }
        h_cats[icat]->SetLineWidth(2);
//        h_cats[icat]->Print();
        h_cats[icat]->SetStats(kFALSE);
        h_cats[icat]->GetXaxis()->SetTitle(xTitle.c_str());
        h_cats[icat]->GetYaxis()->SetTitle("Number of particles");
        if (ieta == m_etaRange.size() - 1 && ipt == m_ptRange.size()) {
          h_cats[icat]->SetTitle(Form("%1.1f < |#eta_{EM2}| && %.0f GeV < p_{T}", m_etaRange[ieta], m_ptRange[ipt]));
        } else if (ieta == m_etaRange.size() - 1 && ipt != m_ptRange.size()) {
          h_cats[icat]->SetTitle(Form("%1.1f < |#eta_{EM2}| && %.0f < p_{T} < %.0f GeV", m_etaRange[ieta], m_ptRange[ipt], m_ptRange[ipt + 1]));
        } else if (ieta != m_etaRange.size() - 1 && ipt == m_ptRange.size() - 1) {
          h_cats[icat]->SetTitle(Form("%1.1f < |#eta_{EM2}| < %1.1f && %.0f GeV < p_{T}", m_etaRange[ieta], m_etaRange[ieta + 1], m_ptRange[ipt]));
        } else {
          h_cats[icat]->SetTitle(Form("%1.1f < |#eta_{EM2}| < %1.1f && %.0f < p_{T} < %.0f GeV", m_etaRange[ieta], m_etaRange[ieta + 1], m_ptRange[ipt], m_ptRange[ipt + 1]));
        }
        h_cats[icat]->SetLineColor(tcolor[icat]);
        std::string lable;
        double entries = h_cats[icat]->Integral();
        if (c_showNumber) {
          lable = Form("%s (%.0f)", display[icat].c_str(), entries);
        } else {
          lable = Form("%s", display[icat].c_str());
        }

        Legend->AddEntry(h_cats[icat], lable.c_str(), "l");
        if (entries == 0) {
          continue;
        }
        if (icat == 0) {
          h_cats[icat]->Draw("hist");
        } else {
          h_cats[icat]->Draw("samehist");
        }

        Legend->Draw();
      }
      h_cats[0]->Draw("axissame");
      Can_Efficiency->SaveAs(Form("%s%s%d_%d_%d.eps", outfolder, catagory[0].c_str(), mode, ieta, ipt));
      delete Can_Efficiency;
    }
  }
  return;
}

void PFlowMonitor::AverageEfficiencyPerPt(const char* outfolder) {

  bool c_showNumber(!false);
  int tcolor[3] = { 4, 2, 8 };

  std::vector<std::string> catagory;
  catagory.clear();
  TH1F* h_cats[4], *h_result[2];
  catagory.push_back("h_efficiency9_pt");
  catagory.push_back("h_ntracks9_pt");
  catagory.push_back("h_efficiency5_pt");
  catagory.push_back("h_ntracks5_pt");

  double xpos(0.25), ypos(0.88);
  TLegend* Legend = new TLegend(xpos, ypos - 0.08 * 3, xpos + 0.3, ypos);
  Legend->SetFillStyle(0);
  Legend->SetBorderSize(0);
  Legend->SetTextFont(43);
  Legend->SetTextSize(20);

  std::string xTitle = "p^{true}_{T} [GeV]";

  for (unsigned int icat = 0; icat < (m_debug ? 1 : catagory.size() / 2); ++icat) {
    bool empty(true);
    TCanvas* Can_AverEfficiency = new TCanvas("AverEfficiency", "AverEfficiency", 450, 400);
    for (unsigned int ieta = 0; ieta < (m_debug ? 1 : m_etaRange.size()); ++ieta) {

      std::string names = histName(ieta, catagory[2*icat], m_etaRange);
      std::string avers = histName(ieta, catagory[2*icat+1], m_etaRange);
      for (unsigned int ifile = 0; ifile < HistFile.size(); ++ifile) {
        if (m_debug) std::cout << HistFile[ifile]->GetName() << std::endl;
        if (ifile == 0) {
          h_cats[2*icat] = (TH1F*) HistFile[ifile]->Get(names.c_str());
          h_cats[2*icat+1] = (TH1F*) HistFile[ifile]->Get(avers.c_str());
          std::cout << __LINE__ << " Add: " << HistFile[ifile]->GetName() << " " << names << std::endl;
        } else {
          h_cats[2*icat]->Add((TH1F*) HistFile[ifile]->Get(names.c_str()));
          h_cats[2*icat+1] = (TH1F*) HistFile[ifile]->Get(avers.c_str());
          std::cout << __LINE__ << " Add: " << HistFile[ifile]->GetName() << " "<< names << std::endl;
        }
      }

      if (!h_cats[2*icat] || !h_cats[2*icat+1]) {
        std::cerr << "[ERROR]\t Histogram " << names <<" or " << avers << " not exist!" << std::endl;
      }

      h_result[icat] = (TH1F*) h_cats[2 * icat]->Clone();
      if (h_result[icat]) continue;

      h_result[icat]->Divide(h_cats[2 * icat + 1]);

      h_result[icat]->SetLineWidth(2);
      h_result[icat]->SetStats(kFALSE);
      h_result[icat]->GetXaxis()->SetTitle(xTitle.c_str());
      h_result[icat]->GetYaxis()->SetTitle("Matching efficiency");
      h_result[icat]->SetLineColor(tcolor[ieta]);
      std::string lable;
      double entries = h_cats[2 * icat + 1]->Integral(); // number of tracks

      std::string display;
      if (ieta == m_etaRange.size() - 1) {
        display = (Form("%1.1f < |#eta_{EM2}|", m_etaRange[ieta]));
      } else {
        display = (Form("%1.1f < |#eta_{EM2}| < %1.1f", m_etaRange[ieta], m_etaRange[ieta + 1]));
      }

      if (c_showNumber) {
        lable = Form("%s (%.0f)", display.c_str(), entries);
      } else {
        lable = Form("%s", display.c_str());
      }

      Legend->AddEntry(h_result[icat], lable.c_str(), "l");

      if (entries == 0) {
        continue;
      }
      if (ieta == 0) {
        h_result[icat]->GetYaxis()->SetRangeUser(0., 1);
        h_result[icat]->Draw("hist");
      } else {
        h_result[icat]->Draw("samehist");
      }
      empty = false;
    }
    if (!empty) {
      Legend->Draw();
      h_result[0]->Draw("axissame");
      Can_AverEfficiency->SaveAs(Form("%sAverageEfficiency_%d.eps", outfolder, icat));
      delete Legend;
    }
      delete Can_AverEfficiency;
  }
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
  if (i_pt != PtRange.size() - 1 && i_eta != EtaRange.size() - 1) {
    complete_name = (name + matchScheme + "_" + Int_to_String((int) (PtRange.at(i_pt))) + "_" + Int_to_String((int) (PtRange.at(i_pt + 1))) + "GeV__eta"
                     + Int_to_String((int) ((10 * EtaRange.at(i_eta)))) + "_" + Int_to_String((int) ((10 * EtaRange.at(i_eta + 1))))).c_str();
    legend_name = Int_to_String((int) (PtRange.at(i_pt))) + "-" + Int_to_String((int) (PtRange.at(i_pt + 1))) + "GeV";
  } else {
    complete_name = (name + matchScheme + "_" + Int_to_String((int) (PtRange.at(i_pt))) + "GeV__eta" + Int_to_String((int) ((10 * EtaRange.at(i_eta))))).c_str();
    if (i_pt != PtRange.size() - 1 && i_eta == EtaRange.size() - 1) {
      legend_name = Int_to_String((int) (PtRange.at(i_pt))) + "-" + Int_to_String((int) (PtRange.at(i_pt + 1))) + "GeV";
    } else {
      legend_name = ">" + Int_to_String((int) (PtRange.at(i_pt))) + "GeV";
    }
  }

  return std::make_pair(complete_name, legend_name);
}

std::string PFlowMonitor::histName(unsigned i_eta, const std::string& name, std::vector<float>& EtaRange) {

  std::string complete_name = "WrongName";
  if (i_eta != EtaRange.size() - 1) {
    complete_name = (name + "_eta" + std::to_string((int) ((10 * EtaRange.at(i_eta)))) + "_" + std::to_string((int) ((10 * EtaRange.at(i_eta + 1))))).c_str();
  } else {
    complete_name = (name + "_eta" + std::to_string((int) ((10 * EtaRange.at(i_eta))))).c_str();
  }

  return complete_name;
}

std::string PFlowMonitor::Int_to_String(int n) {
  std::ostringstream stream;
  stream << n;
  return stream.str();
}
