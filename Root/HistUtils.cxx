#include <PFlowAna/xAODPFlowAna.h>
#include <EventLoop/Worker.h>

//////////////////////////
// Performance Histograms
//////////////////////////
void xAODPFlowAna :: PerformanceHistos(){
    return;
}


//////////////////////////
// Histograms booking
//////////////////////////

void xAODPFlowAna :: bookH1DHistogram(std::string name, int n_bins, float x_low, float x_up){
  Info("bookH1DHistogram () ", "bookH1DHistograms...");
  
  TH1D* h1 = new TH1D(name.c_str(),name.c_str(),n_bins, x_low, x_up);
  h1->Sumw2();
  m_H1Dict[name] = h1;
  wk()->addOutput (m_H1Dict[name]);
  return; 
}

std::string xAODPFlowAna::histName(unsigned i_pt, unsigned i_eta, const std::string& name, const std::string& matchScheme, std::vector<float>& PtRange,
                                   std::vector<float>& EtaRange) {

  std::string complete_name = "WrongName";
  if (i_pt != PtRange.size()-1 && i_eta != EtaRange.size()-1) {
    complete_name = (name + matchScheme + "_" + std::to_string((int) (PtRange.at(i_pt))) + "_" + std::to_string((int) (PtRange.at(i_pt+1))) + "GeV__eta"
                     + std::to_string((int) ((10 * EtaRange.at(i_eta)))) + "_" + std::to_string((int) ((10 * EtaRange.at(i_eta+1))))).c_str();
  }
  else {
    complete_name = (name + matchScheme + "_" + std::to_string((int) (PtRange.at(i_pt))) + "GeV__eta"
                      + std::to_string((int) ((10 * EtaRange.at(i_eta)))) ).c_str();
  }
  return complete_name;
}


void xAODPFlowAna :: bookH1DPerformanceHistogram(std::string name, std::string matchScheme, std::vector<float> PtRange, std::vector<float> EtaRange, int n_bins, float x_low, float x_up)
{
  Info("bookH1DPerformanceHistogram () ", "bookH1DPerformanceHistograms...");

   for (unsigned i_pt = 0; i_pt<PtRange.size(); i_pt++){
    for (unsigned i_eta =0; i_eta<EtaRange.size(); i_eta++){
      
      std::string complete_name = histName(i_pt, i_eta, name, matchScheme, PtRange, EtaRange);

      TH1D* h1 = new TH1D(complete_name.c_str(), complete_name.c_str(),n_bins, x_low, x_up);
      h1->Sumw2();
      m_H1Dict[complete_name] = h1;
      wk()->addOutput (m_H1Dict[complete_name]);
    }
  }
  return ;
} 


std::string xAODPFlowAna::histSubName(unsigned i_R, unsigned i_eta, const std::string& name, std::vector<int>& DeltaR, std::vector<float>& EtaRange) {

  std::string complete_name = "WrongName";
  if (i_eta != EtaRange.size()-1) {
    complete_name = (name + std::to_string((int) (DeltaR.at(i_R))) + "_eta" +
		     std::to_string((int) ((10 * EtaRange.at(i_eta)))) + "_" + std::to_string((int) ((10 * EtaRange.at(i_eta+1))))).c_str();
  }
  else {
    complete_name = (name + std::to_string((int) (DeltaR.at(i_R))) + "_eta" +
		     std::to_string((int) ((10 * EtaRange.at(i_eta))))).c_str();
  }
  
  return complete_name;
}



void xAODPFlowAna :: bookSubHistogram(std::string name,  std::vector<int> DeltaR, std::vector<float> EtaRange, int n_binsX, const Double_t *Xaxis,  int n_binsY,  const Double_t *Yaxis, std::string histType){
  
  for (unsigned i_R = 0; i_R<DeltaR.size(); i_R++){
    for (unsigned i_eta = 0; i_eta<EtaRange.size(); i_eta++){
      
      std::string complete_name = histSubName(i_R, i_eta, name, DeltaR, EtaRange);
      
      if(histType.compare("TH2D") == 0){
	TH2D* h2 = new TH2D(complete_name.c_str(), complete_name.c_str(), n_binsX, Xaxis, n_binsY, Yaxis);
       	h2->Sumw2();
       	m_H2Dict[complete_name] = h2;
       	wk()->addOutput (m_H2Dict[complete_name]);
      }
      
      else if(histType.compare("TProfile") == 0){
       	TProfile2D* TProfile = new TProfile2D(complete_name.c_str(), complete_name.c_str(), n_binsX, Xaxis, n_binsY, Yaxis);
      	TProfile->Sumw2();
       	m_TProfDict[complete_name] = TProfile;
       	wk()->addOutput (m_TProfDict[complete_name]);
      }
      else  Info("bookSubHistogram ()", "ERROR: no histogram type valid!");
    }
  }
  
  return; 
}



std::string xAODPFlowAna::histSubName2(unsigned i_eta, const std::string& name,std::vector<float>& EtaRange) {
  
  std::string complete_name = "WrongName";
  if (i_eta != EtaRange.size()-1) {
    complete_name = (name + "_eta" +
		     std::to_string((int) ((10 * EtaRange.at(i_eta)))) + "_" + std::to_string((int) ((10 * EtaRange.at(i_eta+1))))).c_str();
  }
  else {
    complete_name = (name + "_eta" +
		     std::to_string((int) ((10 * EtaRange.at(i_eta))))).c_str();
  }
  return complete_name;
}




void xAODPFlowAna :: bookSubHistogram2(std::string name,  std::vector<float> EtaRange, int n_binsX, const Double_t *Xaxis,  int n_binsY,   float y_low, float y_up){
  
 
  for (unsigned i_eta = 0; i_eta<EtaRange.size(); i_eta++){

    std::string complete_name = histSubName2(i_eta, name, EtaRange);
    TH2D* h2 = new TH2D(complete_name.c_str(), complete_name.c_str(), n_binsX, Xaxis, n_binsY, y_low, y_up);
    h2->Sumw2();
    m_H2Dict[complete_name] = h2;
    wk()->addOutput (m_H2Dict[complete_name]);
  }
  
  return; 
}

