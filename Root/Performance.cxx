#include <PFlowAna/xAODPFlowAna.h>


//////////////////////////
// Performance Histograms
//////////////////////////
void xAODPFlowAna :: PerformanceHistos(){  
  return;
}


/////////////////////////////
// Resize TruthParticle vectors
/////////////////////////////
void xAODPFlowAna :: resize_tpVectors(const xAOD::TruthParticleContainer* TruthParticles){  
  
  //Vectors for storing true information
  _mc_hasEflowTrack.resize(TruthParticles->size()); 
  _mc_hasEflowTrackIndex.resize(TruthParticles->size());
  _mc_trueE.resize(TruthParticles->size());
  _mc_trueEafter.resize(TruthParticles->size());

  return;
}

/////////////////////////////
// Resize PFO vectors 
/////////////////////////////
void xAODPFlowAna :: resize_PFOVectors(const xAOD::PFOContainer* JetETMissChargedParticleFlowObjects){  

 //Vector for storing PFO variables
  _pfo_Pt.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_iniEoPexp.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_inisigmaEoPexp.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_LFI.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_hasClusterMatched.resize(JetETMissChargedParticleFlowObjects->size());
  _clMatchedEflow.resize(JetETMissChargedParticleFlowObjects->size());  
  _clMatchedEflowEcone15.resize(JetETMissChargedParticleFlowObjects->size());
  
  return;
}

/////////////////////////////
// Fill PFO vectors 
/////////////////////////////
void xAODPFlowAna :: fill_PFOVectors(const xAOD::PFOContainer* JetETMissChargedParticleFlowObjects){  
  xAOD::PFOContainer::const_iterator cpfo_itr = JetETMissChargedParticleFlowObjects->begin();
  xAOD::PFOContainer::const_iterator cpfo_end = JetETMissChargedParticleFlowObjects->end();

  for( ; cpfo_itr != cpfo_end; ++cpfo_itr ) {
    
    int cpfo_index = std::distance(JetETMissChargedParticleFlowObjects->begin(),cpfo_itr);
    _pfo_Pt.at(cpfo_index) = (*cpfo_itr)->pt();
    _pfo_iniEoPexp.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("EExpect");
    _pfo_inisigmaEoPexp.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("varEExpect");
    _pfo_LFI.at(cpfo_index) = (*cpfo_itr)->auxdata< int >("FirstIntLayer");

    //Not clear if this is the correct way to match clusters!
    const xAOD::CaloCluster* matchedCluster = (*cpfo_itr)->cluster(0);
    float matchedCluster_E = -999;   float matchedCluster_eta = -999;   float matchedCluster_phi = -999;
    if (matchedCluster) {matchedCluster_E = matchedCluster->rawE(); matchedCluster_eta = matchedCluster->eta(); matchedCluster_phi = matchedCluster->phi();}
    Info("fill_PFOVectors()", " matchedCluster_E  = %.3f, eta = %.3f, phi = %.3f ",matchedCluster_E, matchedCluster_eta, matchedCluster_phi);
    
    if(matchedCluster) _pfo_hasClusterMatched.at(cpfo_index) = 1;
    else _pfo_hasClusterMatched.at(cpfo_index) = 0;
  }
  return;
}


/////////////////////////////
// Truth particle SELECTION
/////////////////////////////
//Step 1: Calculate the link for charge Pflow objects and mc particles
void xAODPFlowAna :: tp_Selection(const xAOD::TruthParticleContainer* TruthParticles,const xAOD::PFOContainer* JetETMissChargedParticleFlowObjects){  
  
  xAOD::TruthParticleContainer::const_iterator tp_itr = TruthParticles->begin();
  xAOD::TruthParticleContainer::const_iterator tp_end = TruthParticles->end();
  int tp_index = 0;
  for( ; tp_itr != tp_end; ++tp_itr ) {
    tp_index = std::distance(TruthParticles->begin(),tp_itr);

    //Associate truth particle to PFO if it exists
    if ((*tp_itr)->status()==1){ //Stable final-state particle status = 1
      if (fabs((*tp_itr)->pdgId())==211      //pi+=211;
	  || fabs((*tp_itr)->pdgId())==321   //K+=321; 
	  || fabs((*tp_itr)->pdgId())==2212  //p+=2212; 
	  || fabs((*tp_itr)->pdgId())==3222  //Sigma+=3222; 
	  || fabs((*tp_itr)->pdgId())==3112  //Sigma-=3112; 
	  || fabs((*tp_itr)->pdgId())==3312  //Xi-=3312; ; 
	  || fabs((*tp_itr)->pdgId())==3334  //Omega-=3334; 
	  || fabs((*tp_itr)->pdgId())==11    //e-=11; 
	  || fabs((*tp_itr)->pdgId())==13){  //mu-=13
	
	const xAOD::TruthVertex* prodVtx =  (*tp_itr)->prodVtx();
	float tp_z0 = 999;
	if(prodVtx) tp_z0 = prodVtx->z();

	xAOD::PFOContainer::const_iterator cpfo_itr = JetETMissChargedParticleFlowObjects->begin();
	xAOD::PFOContainer::const_iterator cpfo_end = JetETMissChargedParticleFlowObjects->end();
	
	int cpfo_index = 0;
	for( ; cpfo_itr != cpfo_end; ++cpfo_itr ) {
	  cpfo_index = std::distance(JetETMissChargedParticleFlowObjects->begin(),cpfo_itr);
	  Info("tp_Selection()", " ChargedParticleFlowObjects charge = %f   E = %.2f GeV  pt = %.2f GeV  eta = %.2f  phi = %.2f ",
	       (*cpfo_itr)->charge(), ((*cpfo_itr)->e()/GEV), ((*cpfo_itr)->pt()/GEV), (*cpfo_itr)->eta(), (*cpfo_itr)->phi());

	  const xAOD::TrackParticle* ptrk = (*cpfo_itr)->track(0);
	  float z0 = 999; 
	  if (ptrk) { z0 = ptrk->z0();}
	  
	  //Check the criteria, why we don´t use DeltaR. Are there truth particles are associated to the same cpfo?
	  
	  if((*cpfo_itr)->charge()!=0 && (*cpfo_itr)->pt()!=0 
	     && fabs((*cpfo_itr)->eta()) < 2.5 
	     && fabs((*cpfo_itr)->eta()-(*tp_itr)->eta()) < 0.02
	     && acos(cos((*cpfo_itr)->phi())*cos((*tp_itr)->phi())+sin((*cpfo_itr)->phi())*sin((*tp_itr)->phi()))<0.02
	     && (fabs(z0-tp_z0)*sin(2*atan(exp(-1.0*(*cpfo_itr)->eta())))) < 2.){
	    
	    _mc_hasEflowTrack.at(tp_index) = 1; //=1 indicates that we have a eflowTrack matched to the mc particle
	    _mc_hasEflowTrackIndex.at(tp_index) = cpfo_index; //say us which eflowObject corresponds for each mc particle (not association = 0)
	    //std::cout<<" Has eflow track associated "<<std::endl;
	  }
	}
      } 
    }
  }
  return;
}
 
 
/////////////////////////////////////
// Calibration Hits for mc particle 
/////////////////////////////////////
//Step 2: Calculate the true energy per particle
//The sum of all calibration hits in topoclusters for the ith mc particle

void xAODPFlowAna :: ComputeCalibHitsPerParticle(const xAOD::CalCellInfoContainer* CalCellInfo_TopoCluster,const xAOD::CalCellInfoContainer* CalCellInfo,
						 const xAOD::TruthParticleContainer* TruthParticles){  
   
  xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_itr = CalCellInfo_TopoCluster->begin();
  xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_end = CalCellInfo_TopoCluster->end();
  for( ; CalCellInfoTopoCl_itr != CalCellInfoTopoCl_end; ++CalCellInfoTopoCl_itr ) {

    xAOD::TruthParticleContainer::const_iterator tp_itr = TruthParticles->begin();
    xAOD::TruthParticleContainer::const_iterator tp_end = TruthParticles->end();
    for( ; tp_itr != tp_end; ++tp_itr ) {
      int i_mcPart = std::distance( TruthParticles->begin(),tp_itr);
      //checking the barcode of the particle is the same associated to the cell and fill information related with mc particles
      if((*CalCellInfoTopoCl_itr)->barcode()==(*tp_itr)->barcode()){
	_mc_trueE.at(i_mcPart) += (*CalCellInfoTopoCl_itr)->EMEnergy()*(*CalCellInfoTopoCl_itr)->cellWeight();
	_mc_trueE.at(i_mcPart)  += (*CalCellInfoTopoCl_itr)->nonEMEnergy()*(*CalCellInfoTopoCl_itr)->cellWeight();
      }	   
    }//end loop over mc particles
  }//end loop over CalCellInfo_TopoCluster
  
  
   //The sum of all calibration hits in neutral eflow objects for the ith mc particle
  xAOD::CalCellInfoContainer::const_iterator CalCellInfoPFO_itr = CalCellInfo->begin();
  xAOD::CalCellInfoContainer::const_iterator CalCellInfoPFO_end = CalCellInfo->end();
  for( ; CalCellInfoPFO_itr != CalCellInfoPFO_end; ++CalCellInfoPFO_itr ) {
    xAOD::TruthParticleContainer::const_iterator tp_itr = TruthParticles->begin();
    xAOD::TruthParticleContainer::const_iterator tp_end = TruthParticles->end();
    for( ; tp_itr != tp_end; ++tp_itr ) {
      int i_mcPart = std::distance( TruthParticles->begin(),tp_itr);
      if((*CalCellInfoPFO_itr)->barcode()==(*tp_itr)->barcode()){
 	_mc_trueEafter.at(i_mcPart) += (*CalCellInfoPFO_itr)->EMEnergy()*(*CalCellInfoPFO_itr)->cellWeight();
	_mc_trueEafter.at(i_mcPart)  += (*CalCellInfoPFO_itr)->nonEMEnergy()*(*CalCellInfoPFO_itr)->cellWeight();
      }
    }//end loop over mc particles
  }//end loop over CalCellPFO
  
  
  //printing information 
  for (unsigned int i=0;i<_mc_trueEafter.size();i++){
    if(_mc_hasEflowTrack.at(i)==1){
      Info("ComputeCalibHitsPerParticle()", " The MC particle with index %i has a PFO associated with trueE = %.2f GeV and trueEafter = %.2f GeV ",
	   _mc_hasEflowTrackIndex[i], _mc_trueE[i]/GEV , _mc_trueEafter[i]/GEV);
    }
  }
  
  return;
}


/////////////////////////////////////
// Subtraction algorithm
/////////////////////////////////////
// Calculate the energy in a cone 0.15 around the CPFO axis

void xAODPFlowAna :: SubtractionPerf(const xAOD::PFOContainer* JetETMissChargedParticleFlowObjects,const xAOD::CaloClusterContainer* topocluster, const xAOD::TruthParticleContainer* TruthParticles){ 
  
  xAOD::PFOContainer::const_iterator cpfo_itr = JetETMissChargedParticleFlowObjects->begin();
  xAOD::PFOContainer::const_iterator cpfo_end = JetETMissChargedParticleFlowObjects->end();
  int cpfo_index = 0;
  for( ; cpfo_itr != cpfo_end; ++cpfo_itr ) {
    cpfo_index = std::distance(JetETMissChargedParticleFlowObjects->begin(),cpfo_itr);
    double Econe15 = 0.0;

    xAOD::CaloClusterContainer::const_iterator CaloCluster_itr = topocluster->begin();
    xAOD::CaloClusterContainer::const_iterator CaloCluster_end = topocluster->end();
    for( ; CaloCluster_itr != CaloCluster_end; ++CaloCluster_itr  ) {
      double dEta = ((*cpfo_itr)->eta()-(*CaloCluster_itr)->eta());
      double dPhi = fabs((*cpfo_itr)->phi()-(*CaloCluster_itr)->phi());
      if(dPhi > M_PI) dPhi = 2*M_PI - dPhi;
      if(sqrt(dEta*dEta + dPhi*dPhi)<=0.15){
	Econe15 += (*CaloCluster_itr)->rawE();
      }
    }
    _clMatchedEflowEcone15.at(cpfo_index) = Econe15;
  }

  //printing information 
  for (unsigned int i=0;i<_clMatchedEflowEcone15.size();i++) Info("SubtractionPerf()", " _clMatchedEflowEcone15 E  = %.2f GeV ", _clMatchedEflowEcone15[i]/GEV);
  
  //Calculate the variables
  xAOD::TruthParticleContainer::const_iterator tp_itr = TruthParticles->begin();
  xAOD::TruthParticleContainer::const_iterator tp_end = TruthParticles->end();
  for( ; tp_itr != tp_end; ++tp_itr ) {
    int i = std::distance(TruthParticles->begin(),tp_itr);
    
    if (_mc_hasEflowTrack[i]==1){
      //if (_mc_hasEflowTrack[i]==1 && (_mc_trueE[i]/(*tp_itr)->pt()*cosh((*tp_itr)->eta()))>0.1){
      // std::cout<<" Has eflow track associated "<<std::endl;
            
      if (std::fabs((_pfo_Pt.at(_mc_hasEflowTrackIndex[i])-(*tp_itr)->pt())/(*tp_itr)->pt())<0.05){
	//std::cout<<" (_pfo_Pt.at(_mc_hasEflowTrackIndex[i])-(*tp_itr)->pt())/(*tp_itr)->pt())<0.05 "<<std::endl;
	
	if (std::fabs((*tp_itr)->eta())>0.0 && std::fabs((*tp_itr)->eta())<0.4){
	  // std::cout<<" std::fabs((*tp_itr)->eta())>0.0 && std::fabs((*tp_itr)->eta())<0.4 "<<std::endl;
	  // std::cout<<"_mc_trueE[i]="<<_mc_trueE[i]<<" ((*tp_itr)->pt()*cosh((*tp_itr)->eta()))="<<((*tp_itr)->pt()*cosh((*tp_itr)->eta()))<< "Ratio="<<_mc_trueE[i]/((*tp_itr)->pt()*cosh((*tp_itr)->eta()))<<std::endl;
	  
	  hTurnOff_CalHitsOverPt_eta_00_04_hist->Fill((*tp_itr)->pt(),_mc_trueE[i]/((*tp_itr)->pt()*cosh((*tp_itr)->eta())),m_EvtWeight);
	  
	  if (_pfo_hasClusterMatched[_mc_hasEflowTrackIndex[i]]==1 && _pfo_LFI[_mc_hasEflowTrackIndex[i]]!=999){
	    
	    float pull15=(_clMatchedEflowEcone15[_mc_hasEflowTrackIndex[i]]-_pfo_iniEoPexp[_mc_hasEflowTrackIndex[i]])/_pfo_inisigmaEoPexp[_mc_hasEflowTrackIndex[i]];
	    Info("SubtractionPerf", " pull15  = %.2f GeV ", pull15);
	    
	    hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_04_hist->Fill((*tp_itr)->pt(),pull15,_mc_trueEafter[i]/((*tp_itr)->pt()*cosh((*tp_itr)->eta())), m_EvtWeight);
	    hTurnOff_Entries_vs_Pull_eta_00_04_hist->Fill((*tp_itr)->pt(),pull15,m_EvtWeight);
	  }
	}
	
      }
    }
  }

  //  -------------------
  //  _mc_LinkedToTruthJets.resize(TruthParticles->size());
  //if (mc_LinkedToTruthJets[i]>=0){
  //if (MyJetsAreaTruth[mc_LinkedToTruthJets[i]].Pt()>20.*1000. && MyJetsAreaTruth[mc_LinkedToTruthJets[i]].Pt()<30.*1000.)
  //hTurnOff_Entries_vs_Pull_pT_20_30_eta_00_04_hist->Fill(mc_pt->at(i),pull,EvtWeight);
  //hTurnOff_CalHitsRemainingOverPt_vs_Pull_pT_20_30_eta_00_04_hist->Fill(mc_pt->at(i),pull,mc_trueEafter[i]/(mc_pt->at(i)*cosh(mc_eta->at(i))),m_EvtWeight);
  
  return;
}




///////////////////////
// Clear vectors
////////////////////////
void xAODPFlowAna :: clear_PerformanceVectors(){

  //tpVectors
  _mc_hasEflowTrack.clear();
  _mc_hasEflowTrackIndex.clear();
  _mc_trueE.clear();
  _mc_trueEafter.clear();

  //PFOVectors
  _pfo_Pt.clear();
  _pfo_iniEoPexp.clear();
  _pfo_inisigmaEoPexp.clear();
  _pfo_LFI.clear();
  _pfo_hasClusterMatched.clear();
  _clMatchedEflow.clear();
  _clMatchedEflowEcone15.clear();

  //_mc_LinkedToTruthJets.clear();
  return;
}




// Some variable definitions;
// mc_hasEflowTrack[i] - this is 1 if there is a charged eflow object associated with that mc particle
// mc_hasEflowTrackIndex[i] - the index of the eflow charged object associated with the ith mc particle.
// mc_trueE[i] - the sum of all calibration hits in topoclusters for the ith mc particle
// mc_trueEafter[i] - the sum of all calibration hits in neutral eflow objects for the ith mc particle
// clMatchedEflow[p] - the index of the cluster matched to the pth charged eflow object
// eflow_eflow_iniEoPexp[p] - the expected E/p value of the pth charged object
// eflow_eflow_inisigmaEoPexp[p] - the expected width of E/p value of the pth charged object
// clMatchedEflowEcone10[q] - the sum of cluster energies in a cone of radius 0.1 around the qth charged eflow object
// mc_LinkedToTruthJets[i] - the index of the truth jet the mc particle belongs to

//  //turnOff TProfiles
//  for (unsigned int i=0; i<mc_pt->size(); i++){

//    if (mc_hasEflowTrack[i]==1){

//      if (std::fabs((eflow_eflow_pt->at(mc_hasEflowTrackIndex[i])-mc_pt->at(i))/mc_pt->at(i))<0.05){

//        if (std::fabs(mc_eta->at(i))>0.0 && std::fabs(mc_eta->at(i))<0.4){

//        hTurnOff_CalHitsOverPt_eta_00_04_hist->Fill(mc_pt->at(i),mc_trueE[i]/(mc_pt->at(i)*cosh(mc_eta->at(i))),EvtWeight);

//          if (clMatchedEflow[mc_hasEflowTrackIndex[i]]>=0 && eflow_eflow_iniLFIindex->at(mc_hasEflowTrackIndex[i])!=8){

//            float pull=(cl_em_E->at(clMatchedEflow[mc_hasEflowTrackIndex[i]])-eflow_eflow_iniEoPexp->at(mc_hasEflowTrackIndex[i]))/eflow_eflow_inisigmaEoPexp->at(mc_hasEflowTrackIndex[i]);
//            float pull10=(clMatchedEflowEcone10[mc_hasEflowTrackIndex[i]]-eflow_eflow_iniEoPexp->at(mc_hasEflowTrackIndex[i]))/eflow_eflow_inisigmaEoPexp->at(mc_hasEflowTrackIndex[i]);
//            float pull15=(clMatchedEflowEcone15[mc_hasEflowTrackIndex[i]]-eflow_eflow_iniEoPexp->at(mc_hasEflowTrackIndex[i]))/eflow_eflow_inisigmaEoPexp->at(mc_hasEflowTrackIndex[i]);
//            float pull20=(clMatchedEflowEcone20[mc_hasEflowTrackIndex[i]]-eflow_eflow_iniEoPexp->at(mc_hasEflowTrackIndex[i]))/eflow_eflow_inisigmaEoPexp->at(mc_hasEflowTrackIndex[i]);
// hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_04_hist->Fill(mc_pt->at(i),pull,mc_trueEafter[i]/(mc_pt->at(i)*cosh(mc_eta->at(i))),EvtWeight);
// hTurnOff_Entries_vs_Pull_eta_00_04_hist->Fill(mc_pt->at(i),pull,EvtWeight);

// -------------------

//            if (mc_LinkedToTruthJets[i]>=0){
//              if (MyJetsAreaTruth[mc_LinkedToTruthJets[i]].Pt()>20.*1000. && MyJetsAreaTruth[mc_LinkedToTruthJets[i]].Pt()<30.*1000.) hTurnOff_Entries_vs_Pull_pT_20_30_eta_00_04_hist->Fill(mc_pt->at(i),pull,EvtWeight);
//              if (MyJetsAreaTruth[mc_LinkedToTruthJets[i]].Pt()>20.*1000. && MyJetsAreaTruth[mc_LinkedToTruthJets[i]].Pt()<30.*1000.) hTurnOff_CalHitsRemainingOverPt_vs_Pull_pT_20_30_eta_00_04_hist->Fill(mc_pt->at(i),pull,mc_trueEafter[i]/(mc_pt->at(i)*cosh(mc_eta->at(i))),EvtWeight);
