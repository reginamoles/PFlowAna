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

bool m_UseNarrowPtRange;
bool m_UseNarrowEtaRange;
bool m_1to2matching;

std::vector<float> _ptRange;
std::vector<float> _etaRange;

TFile *HistFile;

//Store plots in a ps file
std::vector<TCanvas*> m_CanVect;
void Efficiency();
void printPS();
void setStyle();
std::pair<std::string, std::string> histName(unsigned i_pt, unsigned i_eta, const std::string& name, const std::string& matchScheme, std::vector<float>& PtRange, std::vector<float>& EtaRange);


#endif
