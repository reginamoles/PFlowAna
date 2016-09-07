	#include <PFlowAna/xAODPFlowAna.h>

//////////////////////////
// Zmumu Selection
//////////////////////////

bool xAODPFlowAna :: ZmumuSelection(const xAOD::ElectronContainer* goodElectrons,const xAOD::MuonContainer* goodMuons){
  
  Info("", "------------------- ");
  Info("", "   Zmumu selection   ");
  Info("", "------------------- ");
  
  if (goodElectrons->size()!=0) return false;
  std::cout<<"no electrons"<<std::endl;
  if (goodMuons->size()!=2) return false;
  std::cout<<"two muons"<<std::endl;
  if (goodMuons->at(0)->charge() * (goodMuons->at(1)->charge())!=-1) return false;
  std::cout<<"two muons with opposite charge"<<std::endl;
  if ((goodMuons->at(0)->pt()/GEV) < 25 || (goodMuons->at(1)->pt()/GEV)<25) return false;
  std::cout<<"two muons with pt>25 GeV"<<std::endl;
  if (goodMuons->at(0)->eta() > 2.4 || (goodMuons->at(1)->eta()) > 2.4) return false;
  std::cout<<"two muons central eta"<<std::endl;
  
  // Z TLorentzVector
  TLorentzVector Z = (goodMuons->at(0)->p4()+goodMuons->at(1)->p4());
  Info("Z boson candidate", "Z E = %.2f GeV  pt  = %.2f GeV eta = %.2f  phi =  %.2f",Z.E(), Z.Pt(), Z.Eta(), Z.Phi());

  if (Z.M()/GEV<82 || Z.M()/GEV>102) return false; //Mass widow - to be checked by Christian
  if (Z.Pt()/GEV < 30) return false;

  return true;
}



void xAODPFlowAna :: FillZmumuHistograms(const xAOD::MuonContainer* goodMuons){
  
  if (goodMuons->size()!=2)  Info("", "  ERROR: Something is wrong in the Zmumu selection!   ");
  
  for(int t = 0; t < 2; t++){
    m_H1Dict["h_muonPt"]->Fill((goodMuons->at(t)->pt())/GEV);
    m_H1Dict["h_muonE"]->Fill((goodMuons->at(t)->e())/GEV);
    m_H1Dict["h_muonM"]->Fill((goodMuons->at(t)->m())/GEV);
    m_H1Dict["h_muonEta"]->Fill(goodMuons->at(t)->eta());
    m_H1Dict["h_muonPhi"]->Fill(goodMuons->at(t)->phi());
  }
  
  m_H1Dict["h_ZPt"]->Fill(((goodMuons->at(0)->p4())+(goodMuons->at(1)->p4())).Pt()/GEV);
  m_H1Dict["h_ZE"]->Fill(((goodMuons->at(0)->p4())+(goodMuons->at(1)->p4())).E()/GEV);
  m_H1Dict["h_ZM"]->Fill(((goodMuons->at(0)->p4())+(goodMuons->at(1)->p4())).M()/GEV);
  m_H1Dict["h_ZEta"]->Fill(((goodMuons->at(0)->p4())+(goodMuons->at(1)->p4())).Eta());
  m_H1Dict["h_ZPhi"]->Fill(((goodMuons->at(0)->p4())+(goodMuons->at(1)->p4())).Phi());
  
  return; 
}
  

void xAODPFlowAna :: JetRecoil_Zmumu(const xAOD::MuonContainer* goodMuons, const xAOD::JetContainer* goodJets){
  
  TLorentzVector Z = goodMuons->at(0)->p4()+goodMuons->at(1)->p4();
  //Loop over jets
  int n_RecoilingJets = 0;
  float SumPt_RecoilingJets = 0;

  xAOD::JetContainer::const_iterator jet_itr = goodJets->begin();
  xAOD::JetContainer::const_iterator jet_end = goodJets->end();
  for( ; jet_itr != jet_end; ++jet_itr ) {
    if( (*jet_itr)->pt()/GEV < 20 ) continue;
    if( fabs(deltaPhi((*jet_itr)->phi(), Z.Phi())) > (M_PI - 0.4)) continue;
    n_RecoilingJets++;
    SumPt_RecoilingJets = SumPt_RecoilingJets + (*jet_itr)->pt();
    // *WIP* These distributions has to be crosscheck (odd results for n_RecoilingJets == 1)
    if (n_RecoilingJets == 1)  m_H1Dict["h_ZPt_to_JetPt"]->Fill( (*jet_itr)->pt()/Z.Pt() );
    else if (n_RecoilingJets > 1 )  m_H1Dict["h_ZPt_to_JetPt_sum"]->Fill(SumPt_RecoilingJets/Z.Pt());
  }
  
  return; 
}

