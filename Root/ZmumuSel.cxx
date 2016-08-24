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

  if (Z.M()/GEV<55 || Z.M()/GEV>135) return false; //Mass widow - to be checked by Christian
  if (Z.Pt()/GEV < 30) return false;

  // *WIP* Move the histograms to another function: FillZmumuHistograms();
  
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

  return true;
}

void xAODPFlowAna :: JetRecoil_Zmumu(const xAOD::ElectronContainer* goodElectrons,const xAOD::MuonContainer* goodMuons, const xAOD::JetContainer* goodJets){
  
 //  // Loop for jets
//   // Cut on pt > 40GeV
//   // Find the leading jet, is it recoiling the Zmumu?  DeltaPhi > pi - 0.4
//   // Is it 
  
//   TLorentzVector Z = goodMuons->at(0)->p4()+goodMuons->at(1)->p4();
//   //Loop over jets
//   for( ; jet_itr != jet_end; ++jet_itr ) {
//     //Are the jets ordered by pt? 
//     if(((jet->pt())*/GEV)<20)continue;
//     if(fabs(deltaPhi((jet->phi()), ((goodMuons->at(0)->p4())+(goodMuons->at(1)->p4())).Phi())) > (M_PI - 0.4))continue;
    
// //     if((goodJets->size()!=0) && (numBadJets==0)){
// //       int choice = 0;
// //       double jetptsum = 0;
// //       if((goodJets->size())>1){
// // 	int m_size = 0;
// // 	m_size = goodJets->size();
// // 	for(int counter = 1; counter < m_size; counter ++){
// // 	  jetptsum += goodJets->at(counter)->pt();
// // 	  if ( ( goodJets->at(counter)->pt() ) > ( goodJets->at(choice)->pt()))choice = counter;
// // 	}
// //       }
      
// //       else jetptsum = goodJets->at(0)->pt(); 
// //       Info("execute()", "Jet  %i has passed all the filters", choice + 1);
// //       m_select++;



// //       h_jetPtcorr->Fill( ( goodJets->at(choice)->pt()*0.001));
// //       //Info("execute()", " jet eta = %.2f ", ((jet)->eta()));
// //       h_jetEtacorr->Fill( ( goodJets->at(choice)->eta()));
// //       //Info("execute()", " jet phi = %.2f ", ((jet)->phi()));
// //       h_jetPhicorr->Fill( ( goodJets->at(choice)->phi()));
// //       //Info("execute()", " jet m = %.2f ", ((jet)->m()));
// //       h_jetMcorr->Fill( ( goodJets->at(choice)->m()));
// //       //Info("execute()", "jet energy corrected = %.2f ", ((jet)->e()*0.001));
// //     h_jetEcorr->Fill( ( goodJets->at(choice)->e()*0.001));
    

    
// //     h_Zpt_to_Jetpt->Fill(( goodJets->at(choice)->pt())/(((goodMuons->at(0)->p4())+(goodMuons->at(1)->p4())).Pt()));
// //     h_Zpt_to_Jetpt_sum->Fill(jetptsum/(((goodMuons->at(0)->p4())+(goodMuons->at(1)->p4())).Pt()));

  return; 
    
}


 
