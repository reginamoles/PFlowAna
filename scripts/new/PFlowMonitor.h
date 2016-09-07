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
#include <TProfile.h>
#include <TProfile2D.h>
#include <TCanvas.h>


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



bool m_SinglePionLowPerformanceStudies;
bool m_DijetLowPerformance;
bool m_DijetSubtraction;
bool m_Zmumu;
bool m_UseNarrowPtRange;
bool m_UseNarrowEtaRange;
bool m_PrintDebug;


std::vector<float> _ptRange;
std::vector<float> _etaRange;

TFile *HistFile;

//Subtraction
void ChargeShowerSubtraction(int);

//Histograms plotting
std::vector<string> FillHistNameVector(std::string, int, std::vector<float>);
void getHist(std::vector<std::string>, int, vector<std::string>);
//void getHist(std::string, int, std::vector<float>);
void Draw2DHist(TH2D* hist, const char*, const char*, const char*); 
void DrawLatex(float, float, const char*);
void Compute_CaloResol(std::vector<std::string>);
//void Compute_CaloResol(TH2D*);
void Compute_EafterSub(TH2D*, TH2D*);
void getFinalHist();

//Tools
std::vector<string> getEtaLabelVect(std::vector<float>);
float FitGaussHist(TH1D*,bool,bool); 
std::string to_string_with_precision(float, int);

//Store plots in a ps file
std::vector<TCanvas*> m_CanVect;
void printPS();
void setStyle();

/* int color_histos = kAzure-5; */
/* int truth_color = kBlue+2; */
/* int truthCuts_color = kGreen+2; */
  /* int pfo_color = kYellow; */
  


/* //Histograms */
/* std::vector<TH1F*> HistogramList; */

/* // Style functions */
/* void InitStyle(); */
/* void colorins(); */

/* //Plot functions */
/* void PrintHistosCode(); */ 



#endif
