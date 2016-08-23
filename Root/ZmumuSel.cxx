#include <PFlowAna/xAODPFlowAna.h>

//////////////////////////
// Zmumu Selection
//////////////////////////

bool xAODPFlowAna :: ZmumuSelection(const xAOD::ElectronContainer* goodElectrons,const xAOD::MuonContainer* goodMuons){
  
  // Info("", "------------------- ");
  // Info("", "   Zmumu selection   ");
  // Info("", "------------------- ");
  
  // if (goodElectrons->size()!=0) return false;
  // if (goodMuons->size()!=2) return false;
  // if (goodMuons->at(0)->charge() * (goodMuons->at(1)->charge())!=-1) return false;
  // if ((goodMuons->at(0)->pt()/GEV) < 25 || (goodMuons->at(1)->pt()/GEV)<25) return false;
  // if (goodMuons->at(0)->eta() < 2.4 || (goodMuons->at(1)->eta())<2.4) return false;

  // // Z TLorentzVector
  // TLorentzVector Z = goodMuons->at(0)->p4()+goodMuons->at(1)->p4();
  // if (Z.m()<55 || Z.m()>135) return false; //Mass widow - to be checked by Christian
  // if (Z.pt/GeV < 30) return false;

  // //Fill histograms
  // for(int t = 0; t < 2; t++){
  //   h_muonpt->Fill((goodMuons->at(t)->pt())*0.001);
  //   h_muoneta->Fill(goodMuons->at(t)->eta());
  //   h_muonphi->Fill(goodMuons->at(t)->phi());
  //   h_muone->Fill((goodMuons->at(t)->e())*0.001);
  //   h_muonm->Fill((goodMuons->at(t)->m())*0.001);
  // }
  
  // h_Zpt->Fill(((goodMuons->at(0)->p4())+(goodMuons->at(1)->p4())).Pt()*0.001);
  // h_Zeta->Fill(((goodMuons->at(0)->p4())+(goodMuons->at(1)->p4())).Eta());
  // h_Zphi->Fill(((goodMuons->at(0)->p4())+(goodMuons->at(1)->p4())).Phi());
  // h_Ze->Fill(((goodMuons->at(0)->p4())+(goodMuons->at(1)->p4())).E()*0.001);
  // h_Zm->Fill(((goodMuons->at(0)->p4())+(goodMuons->at(1)->p4())).M()*0.001);
    
  
  return true;
}

//((goodMuons->at(0)->pt()/GEV + goodMuons->at(1)->pt()/GEV))>=30) &&
//((((goodMuons->at(0)->p4()/GEV) + (goodMuons->at(1)->p4())).M()/GEV)>=(65)) &&
//((((goodMuons->at(0)->p4()) + (goodMuons->at(1)->p4())).M())<=(115000.)) ){					


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


 
