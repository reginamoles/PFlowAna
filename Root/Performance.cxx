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

  Info("resize_tpVectors () ", "In resize_tpVectors...");
  //Vectors for storing true information
  _mc_hasEflowTrack.resize(TruthParticles->size());
  _mc_hasEflowTrackIndex.resize(TruthParticles->size());
  _mc_hasEflowTrackP.resize(TruthParticles->size());
  _mc_hasEflowTrackPt.resize(TruthParticles->size());
  _mc_hasEflowTrackEtaAtLayer.resize(TruthParticles->size());
  _mc_matchedClusterHash.resize(TruthParticles->size());
  return;
}

/////////////////////////////
// Resize PFO vectors 
/////////////////////////////
void xAODPFlowAna :: resize_PFOVectors(const xAOD::PFOContainer* JetETMissChargedParticleFlowObjects){
  
  Info("resize_PFOVectors () ", "In resize_PFOVectors...");
 //Vector for storing PFO variables
  _pfo_Pt.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_iniEoPexp.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_inisigmaEoPexp.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_LFI.resize(JetETMissChargedParticleFlowObjects->size());
  if(m_1to2matching) {
  _pfo_SubtractStatus.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_EtaEMB1.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_PhiEMB1.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_EtaEME1.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_PhiEME1.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_EtaEMB2.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_PhiEMB2.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_EtaEME2.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_PhiEME2.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_EtaEMB3.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_PhiEMB3.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_EtaEME3.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_PhiEME3.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_EtaHEC1.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_PhiHEC1.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_EtaHEC2.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_PhiHEC2.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_EtaHEC3.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_PhiHEC3.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_EtaHEC4.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_PhiHEC4.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_EtaTile1.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_PhiTile1.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_EtaTile2.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_PhiTile2.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_EtaTile3.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_PhiTile3.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_EOP1.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_hashCluster1.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_hashCluster2.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_EOPTotal.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_NMatchedClusterInCellLevelSubtraction.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_eMatchedCluster1.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_eMatchedCluster2.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_RpMatchedCluster1.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_RpMatchedCluster2.resize(JetETMissChargedParticleFlowObjects->size());
  }

  _clMatchedEflow.resize(JetETMissChargedParticleFlowObjects->size());  
  _clMatchedEflowEcone15.resize(JetETMissChargedParticleFlowObjects->size());
  
  _pfo_hasClusterMatched.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_hasClusterMatched_Eta.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_hasClusterMatched_Phi.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_hasClusterMatched_E.resize(JetETMissChargedParticleFlowObjects->size());
  _pfo_hasClusterMatched_Index.resize(JetETMissChargedParticleFlowObjects->size());

  
  return;
}

/////////////////////////////
// Resize PFO vectors 
/////////////////////////////
void xAODPFlowAna :: initialise_PFOVectors(int n_mcParticles, int n_clusters, int n_cPFO){
  
  for (int i_mcPart=0; i_mcPart< n_mcParticles; i_mcPart++) {
    _CalHitEPerPar.push_back(0);
    _CalHitEPerParAfterSubtraction.push_back(0);
    for (int i_clus=0; i_clus < n_clusters; i_clus++){_CalHitEPerClusFromOnePart.push_back(0);}
  }
  
  for (int i_clus=0; i_clus < n_clusters; i_clus++){_CalHitEPerClusFromAllPart.push_back(0);}
  
  for (int i_cpfo=0; i_cpfo<n_cPFO; i_cpfo++) {_pfo_hasClusterMatched_E.at(i_cpfo) = -1;}

  return; 
} 

/////////////////////////////
// Fill PFO vectors 
/////////////////////////////
void xAODPFlowAna :: fill_PFOVectors(const xAOD::PFOContainer* JetETMissChargedParticleFlowObjects){

  Info("fill_PFOVectors () ", "In fill_PFOVectors...");
  xAOD::PFOContainer::const_iterator cpfo_itr = JetETMissChargedParticleFlowObjects->begin();
  xAOD::PFOContainer::const_iterator cpfo_end = JetETMissChargedParticleFlowObjects->end();
  int cpfo_index = 0;
  for( ; cpfo_itr != cpfo_end; ++cpfo_itr ) {
    
    cpfo_index = std::distance(JetETMissChargedParticleFlowObjects->begin(),cpfo_itr);
    _pfo_Pt.at(cpfo_index) = (*cpfo_itr)->pt();
    _pfo_iniEoPexp.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("EExpect"); //eflowRec_tracksExpectedEnergyDeposit
    _pfo_inisigmaEoPexp.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("varEExpect"); //eflowRec_tracksExpectedEnergyDepositVariance
    _pfo_LFI.at(cpfo_index) = (*cpfo_itr)->auxdata< int >("FirstIntLayer"); //eflowRec_FirstIntLayer
    if(m_1to2matching) {
    _pfo_SubtractStatus.at(cpfo_index) = (*cpfo_itr)->auxdata< int >("SubtractStatus");
    _pfo_EtaEMB1.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("EtaEMB1");
    _pfo_PhiEMB1.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("PhiEMB1");
    _pfo_EtaEME1.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("EtaEME1");
    _pfo_PhiEME1.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("PhiEME1");
    _pfo_EtaEMB2.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("EtaEMB2");
    _pfo_PhiEMB2.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("PhiEMB2");
    _pfo_EtaEME2.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("EtaEME2");
    _pfo_PhiEME2.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("PhiEME2");
    _pfo_EtaEMB3.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("EtaEMB3");
    _pfo_PhiEMB3.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("PhiEMB3");
    _pfo_EtaEME3.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("EtaEME3");
    _pfo_PhiEME3.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("PhiEME3");
    _pfo_EtaHEC1.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("EtaHEC1");
    _pfo_PhiHEC1.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("PhiHEC1");
    _pfo_EtaHEC2.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("EtaHEC2");
    _pfo_PhiHEC2.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("PhiHEC2");
    _pfo_EtaHEC3.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("EtaHEC3");
    _pfo_PhiHEC3.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("PhiHEC3");
    _pfo_EtaHEC4.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("EtaHEC4");
    _pfo_PhiHEC4.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("PhiHEC4");
    _pfo_EtaTile1.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("EtaTile1");
    _pfo_PhiTile1.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("PhiTile1");
    _pfo_EtaTile2.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("EtaTile2");
    _pfo_PhiTile2.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("PhiTile2");
    _pfo_EtaTile3.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("EtaTile3");
    _pfo_PhiTile3.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("PhiTile3");
    _pfo_hashCluster1.at(cpfo_index) = (*cpfo_itr)->auxdata< long int >("hashCluster1");
    std::cout<<"zhangrui _pfo_hashCluster1.at("<<cpfo_index<<")="<<_pfo_hashCluster1.at(cpfo_index)<<std::endl;
    _pfo_hashCluster2.at(cpfo_index) = (*cpfo_itr)->auxdata< long int >("hashCluster2");
    std::cout<<"zhangrui _pfo_hashCluster2.at("<<cpfo_index<<")="<<_pfo_hashCluster2.at(cpfo_index)<<std::endl;
    _pfo_EOP1.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("EOP1");
    _pfo_EOPTotal.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("EOPTotal");
    _pfo_NMatchedClusterInCellLevelSubtraction.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("NMatchedClusterInCellLevelSubtraction");
    _pfo_eMatchedCluster1.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("eMatchedCluster1");
    _pfo_eMatchedCluster2.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("eMatchedCluster2");
    _pfo_RpMatchedCluster1.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("RpMatchedCluster1");
    _pfo_RpMatchedCluster2.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("RpMatchedCluster2");
    }

    //Not clear if this is the correct way to match clusters!
    const xAOD::CaloCluster* matchedCluster = (*cpfo_itr)->cluster(0);
    if(matchedCluster) {
      _pfo_hasClusterMatched.at(cpfo_index) = 1;
      _pfo_hasClusterMatched_Eta.at(cpfo_index) = matchedCluster->eta();
      _pfo_hasClusterMatched_Phi.at(cpfo_index) = matchedCluster->phi();
      _pfo_hasClusterMatched_E.at(cpfo_index) = matchedCluster->e();}
    else {
      _pfo_hasClusterMatched.at(cpfo_index) = 0;
      _pfo_hasClusterMatched_Eta.at(cpfo_index) = 999;
      _pfo_hasClusterMatched_Phi.at(cpfo_index) = 999;
      _pfo_hasClusterMatched_E.at(cpfo_index) = -999;
    }

    Info("execute()", "Cluster mathced: E  = %.3f, eta = %.3f, phi = %.3f ",
	 _pfo_hasClusterMatched_E.at(cpfo_index),
	 _pfo_hasClusterMatched_Eta.at(cpfo_index),
	 _pfo_hasClusterMatched_Phi.at(cpfo_index));
  }
  
  std::cout<<"Finish fill_PFOVectors function."<<std::endl;
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
	
	  const xAOD::TrackParticle* ptrk = (*cpfo_itr)->track(0);
	  float z0 = 999; 
	  if (ptrk) { z0 = ptrk->z0();}
	  
	  //Check the criteria, why we don´t use DeltaR. Are there truth particles are associated to the same cpfo?
	  
	  if((*cpfo_itr)->charge()!=0 && (*cpfo_itr)->pt()!=0 
	     && fabs((*cpfo_itr)->eta()) < 2.5 
	     && fabs((*cpfo_itr)->eta()-(*tp_itr)->eta()) < 0.02
	     && acos(cos((*cpfo_itr)->phi())*cos((*tp_itr)->phi())+sin((*cpfo_itr)->phi())*sin((*tp_itr)->phi()))<0.02
	     && (fabs(z0-tp_z0)*sin(2*atan(exp(-1.0*(*cpfo_itr)->eta())))) < 2.){

	    std::cout<<"zhangrui in tp_Selection"<<std::endl;
	    _mc_hasEflowTrack.at(tp_index) = 1; //=1 indicates that we have a eflowTrack matched to the mc particle
	    _mc_hasEflowTrackIndex.at(tp_index) = cpfo_index; //say us which eflowObject corresponds for each mc particle (not association = 0)
      _mc_hasEflowTrackP.at(tp_index) = fabs(1. / ptrk->qOverP());
      _mc_hasEflowTrackPt.at(tp_index) =  (*cpfo_itr)->pt();
      if(m_1to2matching) {
      std::cout<<"zhangrui before "<<cpfo_index<<", "<<_pfo_hashCluster1.size()<<","<<_pfo_hashCluster2.size()<<std::endl;
      std::cout<<"zhangrui before _mc_matchedClusterHash="<<_pfo_hashCluster1.at(cpfo_index)<<","<<_pfo_hashCluster2.at(cpfo_index)<<std::endl;
      _mc_matchedClusterHash.at(tp_index) = std::make_pair(_pfo_hashCluster1.at(cpfo_index), _pfo_hashCluster2.at(cpfo_index));
      std::cout<<"zhangrui after _mc_matchedClusterHash="<<_pfo_hashCluster1.at(cpfo_index)<<","<<_pfo_hashCluster2.at(cpfo_index)<<std::endl;
      }
	  }
	}
      } 
    }
  }
  return;
}
 
 
///////////////////////////////////////////
// Calculate the true energy per particle
//////////////////////////////////////////
//Step 2: Calculate the true energy per particle
//The sum of all calibration hits in topoclusters for the ith mc particle

void xAODPFlowAna :: ComputeCalibHitsPerParticle(const xAOD::CalCellInfoContainer* CalCellInfo_TopoCluster,const xAOD::CalCellInfoContainer* CalCellInfo,
						 const xAOD::TruthParticleContainer* TruthParticles){  
  
  xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_itr = CalCellInfo_TopoCluster->begin();
  xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_end = CalCellInfo_TopoCluster->end();
  for( ; CalCellInfoTopoCl_itr != CalCellInfoTopoCl_end; ++CalCellInfoTopoCl_itr ) {
    int i_CalCell = std::distance(CalCellInfo_TopoCluster->begin(),CalCellInfoTopoCl_itr);
    _CalCellInfo_index.push_back(-1);
    xAOD::TruthParticleContainer::const_iterator tp_itr = TruthParticles->begin();
    xAOD::TruthParticleContainer::const_iterator tp_end = TruthParticles->end();
    for( ; tp_itr != tp_end; ++tp_itr ) {
      int i_mcPart = std::distance(TruthParticles->begin(),tp_itr);
      //checking the barcode of the particle is the same associated to the cell and fill information related with mc particles
      if((*CalCellInfoTopoCl_itr)->barcode()==(*tp_itr)->barcode()){
	_CalCellInfo_index.at(i_CalCell)=i_mcPart;  // tell us cell info associated to which mc particle
	_CalHitEPerPar.at(i_mcPart) += (*CalCellInfoTopoCl_itr)->EMEnergy()*(*CalCellInfoTopoCl_itr)->cellWeight();
	_CalHitEPerPar.at(i_mcPart)  += (*CalCellInfoTopoCl_itr)->nonEMEnergy()*(*CalCellInfoTopoCl_itr)->cellWeight();
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
	_CalHitEPerParAfterSubtraction.at(i_mcPart) += (*CalCellInfoPFO_itr)->EMEnergy()*(*CalCellInfoPFO_itr)->cellWeight();
	_CalHitEPerParAfterSubtraction.at(i_mcPart)  += (*CalCellInfoPFO_itr)->nonEMEnergy()*(*CalCellInfoPFO_itr)->cellWeight();
      }
    }//end loop over mc particles
  }//end loop over CalCellPFO
  
   //printing information 
  for (unsigned int i=0;i<_CalHitEPerParAfterSubtraction.size();i++){
    if(_mc_hasEflowTrack.at(i)==1){
      Info("ComputeCalibHitsPerParticle()", " The MC particle with index %i has a PFO associated with CalHitEPerPar = %.2f GeV and CalHitEPerParAfterSubtraction = %.2f GeV ",
   	   _mc_hasEflowTrackIndex[i], _CalHitEPerPar[i]/GEV , _CalHitEPerParAfterSubtraction[i]/GEV);
    }
  }
  
  return;
}

//void xAODPFlowAna :: FindMatchedClusterIndex(const xAOD::PFOContainer* JetETMissChargedParticleFlowObjects, const xAOD::CaloClusterContainer* topocluster) {
//
//  std::map<int, std::vector<long int, long int>> pfoToClusters;
//  xAOD::PFOContainer::const_iterator cpfo_itr = JetETMissChargedParticleFlowObjects->begin();
//  xAOD::PFOContainer::const_iterator cpfo_end = JetETMissChargedParticleFlowObjects->end();
//  int cpfo_index = 0;
//  for( ; cpfo_itr != cpfo_end; ++cpfo_itr ) {
//    cpfo_index = std::distance(JetETMissChargedParticleFlowObjects->begin(),cpfo_itr);
//
//    std::vector<long int, long int> index; index.clear();
//    index.push_back( _pfo_hashCluster1.at(cpfo_index));
//    index.push_back( _pfo_hashCluster2.at(cpfo_index));
//
//    xAOD::CaloClusterContainer::const_iterator CaloCluster_itr = topocluster->begin();
//    xAOD::CaloClusterContainer::const_iterator CaloCluster_end = topocluster->end();
//    for( ; CaloCluster_itr != CaloCluster_end; ++CaloCluster_itr  ) {
//
//      double dEta = ((*cpfo_itr)->eta()-(*CaloCluster_itr)->eta());
//      double dPhi = fabs((*cpfo_itr)->phi()-(*CaloCluster_itr)->phi());
//      if(dPhi > M_PI) dPhi = 2*M_PI - dPhi;
//
//      // matched calo cluster with pflow cluster
//      if(sqrt(dEta*dEta + dPhi*dPhi)<=0.15){
//
//
//        Econe15 += (*CaloCluster_itr)->rawE();
//      }
//    }
//    _clMatchedEflowEcone15.at(cpfo_index) = Econe15;
//  }
//
//
//
//
//
//
//
//
//  std::vector<int> index; index.clear();
//  //For determining the energy before subtraction loop over CalCellInfo_TopoCluster
//  xAOD::CaloClusterContainer::const_iterator CaloCluster_itr = topocluster->begin();
//  xAOD::CaloClusterContainer::const_iterator CaloCluster_end = topocluster->end();
//  for( ; CaloCluster_itr != CaloCluster_end; ++CaloCluster_itr  ) {
//    int i_clus = std::distance(topocluster->begin(),CaloCluster_itr);
//    long int hash =  (*CaloCluster_itr)->pt()*1000 * (*CaloCluster_itr)->eta() * 10 * (*CaloCluster_itr)->phi() * 10;
//
//    if (hash == )
//
//  }
//}

void xAODPFlowAna :: ComputeCalibHitsPerCluster(const xAOD::CalCellInfoContainer* CalCellInfo_TopoCluster, const xAOD::CaloClusterContainer* topocluster, int _n_mcParticles){

  xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_itr = CalCellInfo_TopoCluster->begin();
  xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_end = CalCellInfo_TopoCluster->end();
  for( ; CalCellInfoTopoCl_itr != CalCellInfoTopoCl_end; ++CalCellInfoTopoCl_itr ) {
    int i_CalCell = std::distance(CalCellInfo_TopoCluster->begin(),CalCellInfoTopoCl_itr);
    
    //For determining the energy before subtraction loop over CalCellInfo_TopoCluster
    xAOD::CaloClusterContainer::const_iterator CaloCluster_itr = topocluster->begin();
    xAOD::CaloClusterContainer::const_iterator CaloCluster_end = topocluster->end();
    for( ; CaloCluster_itr != CaloCluster_end; ++CaloCluster_itr  ) {
      int i_clus = std::distance(topocluster->begin(),CaloCluster_itr);
      
//      long int hash =  (*CaloCluster_itr)->pt()*1000 * (*CaloCluster_itr)->eta() * 10 * (*CaloCluster_itr)->phi() * 10;
//      std::cout<<"zhangrui calculate hash="<<hash<<std::endl;

      if (fabs((*CalCellInfoTopoCl_itr)->cellEta() - (*CaloCluster_itr)->rawEta()) < 0.4){
	if (fabs((*CalCellInfoTopoCl_itr)->clusterRecoEnergy()-(*CaloCluster_itr)->rawE())/fabs((*CalCellInfoTopoCl_itr)->clusterRecoEnergy())<0.00001){
	  //energy resolution? Condition is not clear to me!!  	  
  	  if( (_CalCellInfo_index.at(i_CalCell)) != -1){
  	    //Calibration hits E from ONE particle in one cluster
  	    _CalHitEPerClusFromOnePart.at(i_clus*_n_mcParticles+ _CalCellInfo_index.at(i_CalCell)) += (*CalCellInfoTopoCl_itr)->EMEnergy()*(*CalCellInfoTopoCl_itr)->cellWeight();
  	    _CalHitEPerClusFromOnePart.at(i_clus*_n_mcParticles+ _CalCellInfo_index.at(i_CalCell)) += (*CalCellInfoTopoCl_itr)->nonEMEnergy()*(*CalCellInfoTopoCl_itr)->cellWeight();
  	    //Calibration hits E from ALL particle in one cluster
  	    _CalHitEPerClusFromAllPart.at(i_clus)+= (*CalCellInfoTopoCl_itr)->EMEnergy()*(*CalCellInfoTopoCl_itr)->cellWeight();
  	    _CalHitEPerClusFromAllPart.at(i_clus)+= (*CalCellInfoTopoCl_itr)->nonEMEnergy()*(*CalCellInfoTopoCl_itr)->cellWeight();
  	  }
  	}
      }
    }
  }//end loop over cells
  
  return;
}

//===================
// Fill histograms
//====================
void xAODPFlowAna :: fill_RPlus_R0(const xAOD::TruthParticleContainer* TruthParticles){
  xAOD::TruthParticleContainer::const_iterator tp_itr = TruthParticles->begin();
  xAOD::TruthParticleContainer::const_iterator tp_end = TruthParticles->end();
  for( ; tp_itr != tp_end; ++tp_itr ) {
    int i_mcPart = std::distance( TruthParticles->begin(),tp_itr);
    if(_mc_hasEflowTrack.at(i_mcPart)==1 && _CalHitEPerPar.at(i_mcPart)!=0){
      _1MinusChargedR->Fill(1-(_CalHitEPerParAfterSubtraction.at(i_mcPart)/_CalHitEPerPar.at(i_mcPart)));
    }
  }
  return;
}

void xAODPFlowAna::FillEffPurHisto(int i_mcPart, xAOD::TruthParticleContainer::const_iterator tp_itr, const std::vector<double>& v_Efficiency,
                                   const std::vector<double>& v_Purity) {
  // Fill efficiency & purity for cluster matched to histograms
  double max_eff = *max_element(v_Efficiency.begin(), v_Efficiency.end());
  double i_max_eff = distance(v_Efficiency.begin(), max_element(v_Efficiency.begin(), v_Efficiency.end()));


  _mc_hasEflowTrackEtaAtLayer.at(i_mcPart) = (*tp_itr)->eta();
  for (unsigned iptbin = 0; iptbin < _ptRange.size(); ++iptbin) {
    for (unsigned ietabin = 0; ietabin < _etaRange.size(); ++ietabin) {
      bool inRegion[2] = { false, false };
      if (iptbin == _ptRange.size() - 1 && _mc_hasEflowTrackPt.at(i_mcPart) / GEV > _ptRange.at(iptbin)) {
        inRegion[0] = true;
      } else if (_mc_hasEflowTrackPt.at(i_mcPart) / GEV > _ptRange.at(iptbin) && _mc_hasEflowTrackPt.at(i_mcPart) / GEV <= _ptRange.at(iptbin + 1)) {
        inRegion[0] = true;
      }

      if (fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart)) > _etaRange.at(ietabin) && ietabin == _etaRange.size() - 1) {
        inRegion[1] = true;
      } else if (fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart)) > _etaRange.at(ietabin) && fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart)) <= _etaRange.at(ietabin + 1)) {
        inRegion[1] = true;
      }

      if (!(inRegion[0] && inRegion[1])) continue;

      std::string complete_name = histName(iptbin, ietabin, "Eff", "", _ptRange, _etaRange);
      m_H1Dict[complete_name]->Fill(max_eff);
      if (v_Efficiency.at(i_max_eff) > 0.5) {
        std::string complete_name = histName(iptbin, ietabin, "Pur", "", _ptRange, _etaRange);
        m_H1Dict[complete_name]->Fill(v_Purity.at(i_max_eff));
      } //Purity for those clusters with eff>50%
    }
  }
}

void xAODPFlowAna::fillEffPurVectorDefault(
                                                                                 const xAOD::CaloClusterContainer* topocluster, int i_mcPart,
                                                                                 const xAOD::TruthParticleContainer* TruthParticles, std::vector<double>& v_Efficiency,
                                                                                 std::vector<double>& v_Purity) {
  xAOD::CaloClusterContainer::const_iterator CaloCluster_itr = topocluster->begin();
  xAOD::CaloClusterContainer::const_iterator CaloCluster_end = topocluster->end();

  for (; CaloCluster_itr != CaloCluster_end; ++CaloCluster_itr) {
    int i_clus = std::distance(topocluster->begin(), CaloCluster_itr);
    if (_CalHitEPerPar.at(i_mcPart) != 0) {
      float Eff = _CalHitEPerClusFromOnePart.at(i_clus * TruthParticles->size() + i_mcPart) / _CalHitEPerPar.at(i_mcPart);
      v_Efficiency.at(i_clus) = Eff;
    } else {
      v_Efficiency.at(i_clus) = -999;
    }
    Info("fill_Eff_LeadCluster", " v_Efficicency.at(%i)  = %.3f ", i_clus, v_Efficiency.at(i_clus));
    if (_CalHitEPerClusFromAllPart.at(i_clus) != 0) {
      float Pur = _CalHitEPerClusFromOnePart.at(i_clus * TruthParticles->size() + i_mcPart) / _CalHitEPerClusFromAllPart.at(i_clus);
      v_Purity.at(i_clus) = Pur;
    } else {
      v_Purity.at(i_clus) = -999;
    }
    Info("fill_Eff_LeadCluster", " v_PURITY.at(%i)  = %.3f ", i_clus, v_Purity.at(i_clus));
  } //end loop over cluster
  return ;
}

// Efficiency & putity for cluster matched
void xAODPFlowAna::Calculate_Efficiency_Purity(const xAOD::TruthParticleContainer* TruthParticles, int _n_clusters, const xAOD::CaloClusterContainer* topocluster) {
  xAOD::TruthParticleContainer::const_iterator tp_itr = TruthParticles->begin();
  xAOD::TruthParticleContainer::const_iterator tp_end = TruthParticles->end();

  for (; tp_itr != tp_end; ++tp_itr) {
    int i_mcPart = std::distance(TruthParticles->begin(), tp_itr);

    //We rerequire at least 15% of the energy of the true particle
    if (!(_mc_hasEflowTrack.at(i_mcPart) == 1 && (_CalHitEPerPar.at(i_mcPart) / ((*tp_itr)->pt() * cosh((*tp_itr)->eta()))) > 0.15)) continue;
    std::vector<double> v_Efficiency(_n_clusters);
    std::vector<double> v_Purity(_n_clusters);

    if(m_1to2matching) {
    long int clusterHash1 = _mc_matchedClusterHash.at(i_mcPart).first;
    long int clusterHash2 = _mc_matchedClusterHash.at(i_mcPart).second;
    std::cout<<"zhangrui hash="<<clusterHash1<<","<<clusterHash2<<std::endl;

    if(clusterHash2==-1) {
      fillEffPurVectorDefault(topocluster, i_mcPart, TruthParticles, v_Efficiency, v_Purity);
    } else {
      xAOD::CaloClusterContainer::const_iterator CaloCluster_itr = topocluster->begin();
      xAOD::CaloClusterContainer::const_iterator CaloCluster_end = topocluster->end();
      int pos1(-1), pos2(-1);
      std::cout<<"zhangrui pos1="<<pos1<<", pos2="<<pos2<<std::endl;
      for (; CaloCluster_itr != CaloCluster_end; ++CaloCluster_itr) {
        int i_clus = std::distance(topocluster->begin(), CaloCluster_itr);
        long int hash =  (*CaloCluster_itr)->rawE()*1000 * (*CaloCluster_itr)->rawEta() * 10 * (*CaloCluster_itr)->rawPhi() * 10;
        std::cout<<"zhangrui iclus="<<i_clus<<"/"<<topocluster->size()<<", hash="<<hash<<" hash info "<<(*CaloCluster_itr)->rawE()<<", "<<(*CaloCluster_itr)->rawEta()<<", "<<(*CaloCluster_itr)->rawPhi()<<std::endl;
        if (clusterHash1 == hash) {
          pos1 = i_clus;
        }
        if (clusterHash2 == hash) {
          pos2 = i_clus;
        }
      }

      std::cout<<"pos1="<<pos1<<", pos2="<<pos2<<std::endl;
      CaloCluster_itr = topocluster->begin();

      if (_CalHitEPerPar.at(i_mcPart) != 0) {
        float Eff = (_CalHitEPerClusFromOnePart.at(pos1 * TruthParticles->size() + i_mcPart) + _CalHitEPerClusFromOnePart.at(pos2 * TruthParticles->size() + i_mcPart))
            / (_CalHitEPerPar.at(i_mcPart) + _CalHitEPerPar.at(i_mcPart));
        std::cout << "Eff=" << Eff << std::endl;
        v_Efficiency.at(pos1) = Eff;
        v_Efficiency.at(pos2) = Eff;
        std::cout << "zhangrui in eff " << Eff << " pos1=" << pos1 << " pos2=" << pos2 << std::endl;
      }

      Info("fill_Eff_TwoCluster", " v_Efficicency.at(%i) and v_Efficicency.at(%i) = %.3f ", pos1, pos2, v_Efficiency.at(pos1));


      if (_CalHitEPerClusFromAllPart.at(pos1) != 0) {
        float Pur = (_CalHitEPerClusFromOnePart.at(pos1 * TruthParticles->size() + i_mcPart) + _CalHitEPerClusFromOnePart.at(pos2 * TruthParticles->size() + i_mcPart))
            / (_CalHitEPerClusFromAllPart.at(pos1) + _CalHitEPerClusFromAllPart.at(pos2));
        v_Purity.at(pos1) = Pur;
        v_Purity.at(pos2) = Pur;
      }

      Info("fill_Eff_TwoCluster", " v_PURITY.at(%i) and v_PURITY.at(%i)  = %.3f ", pos1, pos2, v_Purity.at(pos1));

    }
    } else {
      fillEffPurVectorDefault(topocluster, i_mcPart, TruthParticles, v_Efficiency, v_Purity);
    }

    FillEffPurHisto(i_mcPart, tp_itr, v_Efficiency, v_Purity);
  }
  return;
}


/*
      
      //Fill Eff & Purity histograms here
      double max_eff = *max_element(v_Efficiency.begin(), v_Efficiency.end());
      double i_max_eff = distance(v_Efficiency.begin(), max_element(v_Efficiency.begin(), v_Efficiency.end()));

    
      if(_mc_hasEflowTrackPt.at(i_mcPart)>1000 && _mc_hasEflowTrackPt.at(i_mcPart)<=2000){
	if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=1){
	  _eff_Lead_1_2GeV_eta1->Fill(max_eff);
	  if(v_Efficiency.at(i_max_eff)>0.5){_pur_Lead_1_2GeV_eta1->Fill(v_Purity.at(i_max_eff));}  //Purity for those clusters with eff>50%
	}
	else if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))>1 && fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=2.0){
	  _eff_Lead_1_2GeV_eta2->Fill(max_eff);
	  if(v_Efficiency.at(i_max_eff)>0.5){_pur_Lead_1_2GeV_eta2->Fill(v_Purity.at(i_max_eff));}
	}
	else{
	  _eff_Lead_1_2GeV_eta25->Fill(max_eff);
	  if(v_Efficiency.at(i_max_eff)>0.5){_pur_Lead_1_2GeV_eta25->Fill(v_Purity.at(i_max_eff));}
	}
      }
      else if(_mc_hasEflowTrackPt.at(i_mcPart)>2000 && _mc_hasEflowTrackPt.at(i_mcPart)<=5000){
	if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=1){
	  _eff_Lead_2_5GeV_eta1->Fill(max_eff);
	  if(v_Efficiency.at(i_max_eff)>0.5){_pur_Lead_2_5GeV_eta1->Fill(v_Purity.at(i_max_eff));} 
	}
	else if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))>1 && fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=2.0){
	  _eff_Lead_2_5GeV_eta2->Fill(max_eff);
	  if(v_Efficiency.at(i_max_eff)>0.5){_pur_Lead_2_5GeV_eta2->Fill(v_Purity.at(i_max_eff));}
	}
	else{
	  _eff_Lead_2_5GeV_eta25->Fill(max_eff);
	  if(v_Efficiency.at(i_max_eff)>0.5){_pur_Lead_2_5GeV_eta25->Fill(v_Purity.at(i_max_eff));}
	}	
      }
      else if(_mc_hasEflowTrackPt.at(i_mcPart)>5000 && _mc_hasEflowTrackPt.at(i_mcPart)<=10000){
	if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=1){
	  _eff_Lead_5_10GeV_eta1->Fill(max_eff);
	  if(v_Efficiency.at(i_max_eff)>0.5){_pur_Lead_5_10GeV_eta1->Fill(v_Purity.at(i_max_eff));} 
	}
	else if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))>1 && fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=2.0){
	  _eff_Lead_5_10GeV_eta2->Fill(max_eff);
	  if(v_Efficiency.at(i_max_eff)>0.5){_pur_Lead_5_10GeV_eta2->Fill(v_Purity.at(i_max_eff));}
	}
	else{
	  _eff_Lead_5_10GeV_eta25->Fill(max_eff);
	  if(v_Efficiency.at(i_max_eff)>0.5){_pur_Lead_5_10GeV_eta25->Fill(v_Purity.at(i_max_eff));}
	}	
      }	
      if(_mc_hasEflowTrackPt.at(i_mcPart)>10000 && _mc_hasEflowTrackPt.at(i_mcPart)<=20000){
	if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=1){
	  _eff_Lead_10GeV_eta1->Fill(max_eff);
	  if(v_Efficiency.at(i_max_eff)>0.5){_pur_Lead_10GeV_eta1->Fill(v_Purity.at(i_max_eff));} 
	}
	else if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))>1 && fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=2.0){
	  _eff_Lead_10GeV_eta2->Fill(max_eff);
	  if(v_Efficiency.at(i_max_eff)>0.5){_pur_Lead_10GeV_eta2->Fill(v_Purity.at(i_max_eff));}
	}
	else{
	  _eff_Lead_10GeV_eta25->Fill(max_eff);
	  if(v_Efficiency.at(i_max_eff)>0.5){_pur_Lead_10GeV_eta25->Fill(v_Purity.at(i_max_eff));}
	}	
      }
      
      //Fill distributions E_exp-E_cl
      cpfo_itr = JetETMissChargedParticleFlowObjects->begin();
      cpfo_end = JetETMissChargedParticleFlowObjects->end();
      for( ; cpfo_itr != cpfo_end; ++cpfo_itr ) {
       	int cpfo_index = std::distance(JetETMissChargedParticleFlowObjects->begin(),cpfo_itr);
	
	if( _pfo_hasClusterMatched_Index.at(cpfo_index) >=0 && v_Efficiency.at(_pfo_hasClusterMatched_Index.at(cpfo_index))>=0.9){
       	  std::cout<<"Pt="<< _mc_hasEflowTrackPt.at(i_mcPart)  <<"  Eta="<<fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))
       		   << "  Ratio="<<((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index))<<std::endl;

	  if(_mc_hasEflowTrackPt.at(i_mcPart)>1000 && _mc_hasEflowTrackPt.at(i_mcPart)<=2000){
      	    if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=1)
      	      _DeltaE_1_2GeV_eta1->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	    if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))>1 && fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=2.0)
      	      _DeltaE_1_2GeV_eta2->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	    else _DeltaE_1_2GeV_eta25->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	  }
	  
      	  if(_mc_hasEflowTrackPt.at(i_mcPart)>2000 && _mc_hasEflowTrackPt.at(i_mcPart)<=5000){
      	    if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=1)
      	      _DeltaE_2_5GeV_eta1->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	    if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))>1 && fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=2.0)
      	      _DeltaE_2_5GeV_eta2->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	    else _DeltaE_2_5GeV_eta25->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	  }

	  if(_mc_hasEflowTrackPt.at(i_mcPart)>5000 && _mc_hasEflowTrackPt.at(i_mcPart)<=10000){
      	    if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=1)
      	      _DeltaE_5_10GeV_eta1->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	    if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))>1 && fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=2.0)
      	      _DeltaE_5_10GeV_eta2->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	    else _DeltaE_5_10GeV_eta25->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	  }

	  if(_mc_hasEflowTrackPt.at(i_mcPart)>10000 && _mc_hasEflowTrackPt.at(i_mcPart)<=20000){
      	    if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=1)
      	      _DeltaE_10GeV_eta1->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	    if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))>1 && fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=2.0)
      	      _DeltaE_10GeV_eta2->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	    else _DeltaE_10GeV_eta25->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	  }
	}
	
      	if( _pfo_hasClusterMatched_Index.at(cpfo_index) >=0 && v_Efficiency.at(_pfo_hasClusterMatched_Index.at(cpfo_index))<=0.7){

	  if(_mc_hasEflowTrackPt.at(i_mcPart)>1000 && _mc_hasEflowTrackPt.at(i_mcPart)<=2000){
      	    if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=1)
      	      _DeltaE07_1_2GeV_eta1->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	    if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))>1 && fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=2.0)
      	      _DeltaE07_1_2GeV_eta2->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	    else _DeltaE07_1_2GeV_eta25->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	  }
	  if(_mc_hasEflowTrackPt.at(i_mcPart)>2000 && _mc_hasEflowTrackPt.at(i_mcPart)<=5000){
      	    if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=1)
      	      _DeltaE07_2_5GeV_eta1->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	    if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))>1 && fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=2.0)
      	      _DeltaE07_2_5GeV_eta2->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	    else _DeltaE07_2_5GeV_eta25->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	  }

	  if(_mc_hasEflowTrackPt.at(i_mcPart)>5000 && _mc_hasEflowTrackPt.at(i_mcPart)<=10000){
      	    if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=1)
      	      _DeltaE07_5_10GeV_eta1->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	    if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))>1 && fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=2.0)
      	      _DeltaE07_5_10GeV_eta2->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	    else _DeltaE07_5_10GeV_eta25->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	  }
	  if(_mc_hasEflowTrackPt.at(i_mcPart)>10000 && _mc_hasEflowTrackPt.at(i_mcPart)<=20000){
      	    if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=1)
      	      _DeltaE07_10GeV_eta1->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	    if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))>1 && fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=2.0)
      	      _DeltaE07_10GeV_eta2->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	    else _DeltaE07_10GeV_eta25->Fill((_pfo_hasClusterMatched_E.at(cpfo_index)-_pfo_iniEoPexp.at(cpfo_index))/_pfo_inisigmaEoPexp.at(cpfo_index));
      	  }
      	}
      }
	
      //link with cluster index destroyed!
      std::sort (v_Efficiency.begin(), v_Efficiency.end());
      
      int NClusters_09 = 0.0;
      double Eff_09 = v_Efficiency[v_Efficiency.size()-1];
      if(Eff_09>0.9) NClusters_09 = 1;
      else{
	Eff_09 += v_Efficiency[v_Efficiency.size()-2];
	if(Eff_09>0.9) NClusters_09 = 2;
	else{
	  Eff_09 += v_Efficiency[v_Efficiency.size()-3];
	  if(Eff_09>0.9) NClusters_09 = 3;
	  else{
	    Eff_09 += v_Efficiency[v_Efficiency.size()-4];
	    if(Eff_09>0.9) NClusters_09 = 4;
	    else{
	      Eff_09 += v_Efficiency[v_Efficiency.size()-5];
	      if(Eff_09>0.9) NClusters_09 = 5;
	      else{
		Eff_09 += v_Efficiency[v_Efficiency.size()-6];
		if(Eff_09>0.9) NClusters_09 = 6;
		else{
		  Eff_09 += v_Efficiency[v_Efficiency.size()-7];
		  if(Eff_09>0.9) NClusters_09 = 7;
		  else{
		    NClusters_09 = 8;
		  }}}}}}}
      _nCluster09->Fill(NClusters_09);
      
      
      //for plotting
      
      if(_mc_hasEflowTrackPt.at(i_mcPart)>1000 && _mc_hasEflowTrackPt.at(i_mcPart)<=2000) {
       	if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=1){_nClus09_1_2GeV_eta1->Fill(NClusters_09);}
      	else if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))>1 && fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=2)
    	  {_nClus09_1_2GeV_eta2->Fill(NClusters_09);}
  	else{_nClus09_1_2GeV_eta25->Fill(NClusters_09);    }
      }
      else if (_mc_hasEflowTrackPt.at(i_mcPart)>2000 && _mc_hasEflowTrackPt.at(i_mcPart)<=5000) {
	if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=1){_nClus09_2_5GeV_eta1->Fill(NClusters_09);}
	else if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))>1 && fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=2)
	  {_nClus09_2_5GeV_eta2->Fill(NClusters_09);}   
	else{_nClus09_2_5GeV_eta25->Fill(NClusters_09);}
      } 
      else if (_mc_hasEflowTrackPt.at(i_mcPart)>5000 && _mc_hasEflowTrackPt.at(i_mcPart)<=10000) {
	if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=1){_nClus09_5_10GeV_eta1->Fill(NClusters_09);}
	else if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))>1 && fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=2)
	  {_nClus09_5_10GeV_eta2->Fill(NClusters_09);}
	else{_nClus09_5_10GeV_eta25->Fill(NClusters_09);}
      }
       else if (_mc_hasEflowTrackPt.at(i_mcPart)>10000 && _mc_hasEflowTrackPt.at(i_mcPart)<=20000) {
	if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=1){_nClus09_10GeV_eta1->Fill(NClusters_09);}
      	else if(fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))>1 && fabs(_mc_hasEflowTrackEtaAtLayer.at(i_mcPart))<=2)
	  {_nClus09_10GeV_eta2->Fill(NClusters_09);}
	else{_nClus09_10GeV_eta25->Fill(NClusters_09);}
      }
    }                   
  }                   
  
  
*/



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
	  hTurnOff_CalHitsOverPt_eta_00_04_hist->Fill((*tp_itr)->pt(), _CalHitEPerPar[i]/((*tp_itr)->pt()*cosh((*tp_itr)->eta())),m_EvtWeight);
	  
	  if (_pfo_hasClusterMatched[_mc_hasEflowTrackIndex[i]]==1 && _pfo_LFI[_mc_hasEflowTrackIndex[i]]!=999){
	    float pull15=(_clMatchedEflowEcone15[_mc_hasEflowTrackIndex[i]]-_pfo_iniEoPexp[_mc_hasEflowTrackIndex[i]])/_pfo_inisigmaEoPexp[_mc_hasEflowTrackIndex[i]];
	    Info("SubtractionPerf", " pull15  = %.2f GeV ", pull15);
	    hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_04_hist->Fill((*tp_itr)->pt(),pull15,_CalHitEPerParAfterSubtraction[i]/((*tp_itr)->pt()*cosh((*tp_itr)->eta())), m_EvtWeight);
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
  _mc_hasEflowTrackP.clear();
  _mc_hasEflowTrackPt.clear();
  _mc_hasEflowTrackEtaAtLayer.clear();
  _mc_matchedClusterHash.clear();

  //PFOVectors
  _pfo_Pt.clear();
  _pfo_iniEoPexp.clear();
  _pfo_inisigmaEoPexp.clear();
  _pfo_LFI.clear();
  if(m_1to2matching) {
  _pfo_SubtractStatus.clear();
  _pfo_EtaEMB1.clear();
  _pfo_PhiEMB1.clear();
  _pfo_EtaEME1.clear();
  _pfo_PhiEME1.clear();
  _pfo_EtaEMB2.clear();
  _pfo_PhiEMB2.clear();
  _pfo_EtaEME2.clear();
  _pfo_PhiEME2.clear();
  _pfo_EtaEMB3.clear();
  _pfo_PhiEMB3.clear();
  _pfo_EtaEME3.clear();
  _pfo_PhiEME3.clear();
  _pfo_EtaHEC1.clear();
  _pfo_PhiHEC1.clear();
  _pfo_EtaHEC2.clear();
  _pfo_PhiHEC2.clear();
  _pfo_EtaHEC3.clear();
  _pfo_PhiHEC3.clear();
  _pfo_EtaHEC4.clear();
  _pfo_PhiHEC4.clear();
  _pfo_EtaTile1.clear();
  _pfo_PhiTile1.clear();
  _pfo_EtaTile2.clear();
  _pfo_PhiTile2.clear();
  _pfo_EtaTile3.clear();
  _pfo_PhiTile3.clear();
  _pfo_hashCluster1.clear();
  _pfo_hashCluster2.clear();
  _pfo_EOP1.clear();
  _pfo_EOPTotal.clear();
  _pfo_NMatchedClusterInCellLevelSubtraction.clear();
  _pfo_eMatchedCluster1.clear();
  _pfo_eMatchedCluster2.clear();
  _pfo_RpMatchedCluster1.clear();
  _pfo_RpMatchedCluster2.clear();
}
  _pfo_hasClusterMatched.clear();
  _pfo_hasClusterMatched_Index.clear();
  _pfo_hasClusterMatched_Eta.clear();
  _pfo_hasClusterMatched_Phi.clear();
  _pfo_hasClusterMatched_E.clear();
  _clMatchedEflow.clear();
  _clMatchedEflowEcone15.clear();


  //CalibrationHits
  _CalCellInfo_index.clear();
  _CalHitEPerPar.clear(); 
  _CalHitEPerParAfterSubtraction.clear();
  _CalHitEPerClusFromOnePart.clear();
  _CalHitEPerClusFromAllPart.clear();
  
   
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
