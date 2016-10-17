#ifndef PFlowMonitor_H
#define PFlowMonitor_H

//=======================================================================================================
// This macro was originally created to draw histograms for PFO performance study (from HISTOS ntuple)
//=======================================================================================================

// ROOT include(s)
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TGaxis.h>


//STL includes
#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstring>
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>

class PFlowMonitor {

public:
  PFlowMonitor();
  ~PFlowMonitor();

  void run(char* inputs, char* outfolder);

  void Plot(const char* outfolder);
  void Efficiency(const char* outfolder);
  void eflowdRp(const char* outfolder, const int mode);
  void AverageEfficiencyPerPt(const char* outfolder);
  void setStyle();
  std::pair<std::string, std::string> histName(unsigned i_pt, unsigned i_eta, const std::string& name, const std::string& matchScheme, std::vector<float>& PtRange,
                                               std::vector<float>& EtaRange);
  std::string histName(unsigned i_eta, const std::string& name, std::vector<float>& EtaRange);

  std::string Int_to_String(int n);

private:
  bool m_UseNarrowPtRange;
  bool m_UseNarrowEtaRange;
  bool m_1to2matching;
  bool m_debug;

  std::vector<float> m_ptRange;
  std::vector<float> m_etaRange;

  std::vector<TFile*> HistFile;
};

#endif
