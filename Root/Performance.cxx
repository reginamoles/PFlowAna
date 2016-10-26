#include <PFlowAna/xAODPFlowAna.h>
#include "xAODCaloEvent/CaloCluster.h"

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

//  Info("resize_tpVectors () ", "In resize_tpVectors...");
  //Vectors for storing true information
  _mc_hasEflowTrack.resize(TruthParticles->size());
  _mc_LFI.resize(TruthParticles->size());
  _mc_hasEflowTrackIndex.resize(TruthParticles->size());
  for(int i=0; i<_mc_hasEflowTrackIndex.size();++i) _mc_hasEflowTrackIndex.at(i) = -1;
  _mc_hasEflowTrackP.resize(TruthParticles->size());
  _mc_hasEflowTrackPt.resize(TruthParticles->size());
  _mc_hasEflowTrackEtaAtLayer.resize(TruthParticles->size());
  _mc_matchedClusterHash.resize(TruthParticles->size());
  _mc_subtractStatus.resize(TruthParticles->size());
  _mc_RpMatchedCluster1.resize(TruthParticles->size());
  _mc_RpMatchedCluster2.resize(TruthParticles->size());
  _mc_etaTrkEM1.resize(TruthParticles->size());
  _mc_phiTrkEM1.resize(TruthParticles->size());
  _mc_etaTrkEM2.resize(TruthParticles->size());
  _mc_phiTrkEM2.resize(TruthParticles->size());
  _mc_etaTrkEM3.resize(TruthParticles->size());
  _mc_phiTrkEM3.resize(TruthParticles->size());
  _mc_pos1.resize(TruthParticles->size());
  _mc_imax.resize(TruthParticles->size());
  _mc_dRp_componets.resize(12*TruthParticles->size());

  return;
}

/////////////////////////////
// Resize PFO vectors 
/////////////////////////////
void xAODPFlowAna :: resize_PFOVectors(const xAOD::PFOContainer* JetETMissChargedParticleFlowObjects){
  
//  Info("resize_PFOVectors () ", "In resize_PFOVectors...");
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
// Resize CaloFromPFO vectors
/////////////////////////////
void xAODPFlowAna::resize_CaloFromPFO(const xAOD::CaloClusterContainer* m_JetETMissCaloClusterObjects) {
  _calo_hash.resize(m_JetETMissCaloClusterObjects->size());
  _calo_EtaVariance.resize(m_JetETMissCaloClusterObjects->size());
  _calo_PhiVariance.resize(m_JetETMissCaloClusterObjects->size());
  _calo_MeanEta.resize(m_JetETMissCaloClusterObjects->size());
  _calo_MeanPhi.resize(m_JetETMissCaloClusterObjects->size());
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
  for (int i_clus=0; i_clus < n_clusters; i_clus++){_CalClusEta.push_back(0);}
  for (int i_clus=0; i_clus < n_clusters; i_clus++){_CalClusPhi.push_back(0);}
  
  for (int i_cpfo=0; i_cpfo<n_cPFO; i_cpfo++) {_pfo_hasClusterMatched_E.at(i_cpfo) = -1;}

  return; 
} 

/////////////////////////////
// Fill PFO vectors 
/////////////////////////////
void xAODPFlowAna :: fill_PFOVectors(const xAOD::PFOContainer* JetETMissChargedParticleFlowObjects){

//  Info("fill_PFOVectors () ", "In fill_PFOVectors...");
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
    _pfo_hashCluster2.at(cpfo_index) = (*cpfo_itr)->auxdata< long int >("hashCluster2");
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

//    Info("execute()", "Cluster mathced: E  = %.3f, eta = %.3f, phi = %.3f ", _pfo_hasClusterMatched_E.at(cpfo_index), _pfo_hasClusterMatched_Eta.at(cpfo_index), _pfo_hasClusterMatched_Phi.at(cpfo_index));
  }
  
  return;
}

void xAODPFlowAna::fill_CaloFromPFO(const xAOD::CaloClusterContainer* JetETMissCaloClusterObjects) {

  xAOD::CaloClusterContainer::const_iterator calopfo_itr = JetETMissCaloClusterObjects->begin();
  xAOD::CaloClusterContainer::const_iterator cpfo_end = JetETMissCaloClusterObjects->end();
  int cpfo_index = 0;
  for( ; calopfo_itr != cpfo_end; ++calopfo_itr ) {

    cpfo_index = std::distance(JetETMissCaloClusterObjects->begin(),calopfo_itr);
    _calo_hash.at(cpfo_index) = (*calopfo_itr)->auxdata< long int >("AllClusterHash");
    _calo_EtaVariance.at(cpfo_index) = (*calopfo_itr)->auxdata< float >("MatchedClusterEtaVariance");
    _calo_PhiVariance.at(cpfo_index) = (*calopfo_itr)->auxdata< float >("MatchedClusterPhiVariance");
    _calo_MeanEta.at(cpfo_index) = (*calopfo_itr)->auxdata< float >("MeanEta");
    _calo_MeanPhi.at(cpfo_index) = (*calopfo_itr)->auxdata< float >("MeanPhi");

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
	
	  const xAOD::TrackParticle* ptrk = (*cpfo_itr)->track(0);
	  float z0 = 999; 
	  if (ptrk) { z0 = ptrk->z0();}
	  
	  //Check the criteria, why we don´t use DeltaR. Are there truth particles are associated to the same cpfo?

	  if((*cpfo_itr)->charge()!=0 && (*cpfo_itr)->pt()!=0 
	     && fabs((*cpfo_itr)->eta()) < 2.5 
	     && AreBothTracksMatched(tp_index,cpfo_index))
	     //&& fabs((*cpfo_itr)->eta()-(*tp_itr)->eta()) < 0.02
	     //&& acos(cos((*cpfo_itr)->phi())*cos((*tp_itr)->phi())+sin((*cpfo_itr)->phi())*sin((*tp_itr)->phi()))<0.02
	     //&& (fabs(z0-tp_z0)*sin(2*atan(exp(-1.0*(*cpfo_itr)->eta())))) < 2.)
      {
	    _mc_hasEflowTrack.at(tp_index) = 1; //=1 indicates that we have a eflowTrack matched to the mc particle
	    _mc_LFI.at(tp_index) = _pfo_LFI.at(cpfo_index);
	    _mc_hasEflowTrackIndex.at(tp_index) = cpfo_index; //say us which eflowObject corresponds for each mc particle (not association = 0)
      _mc_hasEflowTrackP.at(tp_index) = fabs(1. / ptrk->qOverP());
      _mc_hasEflowTrackPt.at(tp_index) =  (*cpfo_itr)->pt();
      _mc_etaTrkEM1.at(tp_index) = _pfo_EtaEMB1.at(cpfo_index) > -990 ? _pfo_EtaEMB1.at(cpfo_index) : _pfo_EtaEME1.at(cpfo_index);
      _mc_phiTrkEM1.at(tp_index) = _pfo_PhiEMB1.at(cpfo_index) > -990 ? _pfo_PhiEMB1.at(cpfo_index) : _pfo_PhiEME1.at(cpfo_index);
      _mc_etaTrkEM2.at(tp_index) = _pfo_EtaEMB2.at(cpfo_index) > -990 ? _pfo_EtaEMB2.at(cpfo_index) : _pfo_EtaEME2.at(cpfo_index);
      _mc_phiTrkEM2.at(tp_index) = _pfo_PhiEMB2.at(cpfo_index) > -990 ? _pfo_PhiEMB2.at(cpfo_index) : _pfo_PhiEME2.at(cpfo_index);
      _mc_etaTrkEM3.at(tp_index) = _pfo_EtaEMB3.at(cpfo_index) > -990 ? _pfo_EtaEMB3.at(cpfo_index) : _pfo_EtaEME3.at(cpfo_index);
      _mc_phiTrkEM3.at(tp_index) = _pfo_PhiEMB3.at(cpfo_index) > -990 ? _pfo_PhiEMB3.at(cpfo_index) : _pfo_PhiEME3.at(cpfo_index);

      if(m_1to2matching) {
      _mc_matchedClusterHash.at(tp_index) = std::make_pair(_pfo_hashCluster1.at(cpfo_index), _pfo_hashCluster2.at(cpfo_index));
      _mc_subtractStatus.at(tp_index) = _pfo_SubtractStatus.at(cpfo_index);
      _mc_RpMatchedCluster1.at(tp_index) = _pfo_RpMatchedCluster1.at(cpfo_index) * _pfo_RpMatchedCluster1.at(cpfo_index);
      _mc_RpMatchedCluster2.at(tp_index) = _pfo_RpMatchedCluster2.at(cpfo_index) * _pfo_RpMatchedCluster2.at(cpfo_index);
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
//      Info("ComputeCalibHitsPerParticle()", " The MC particle with index %i has a PFO associated with CalHitEPerPar = %.2f GeV and CalHitEPerParAfterSubtraction = %.2f GeV ", _mc_hasEflowTrackIndex[i], _CalHitEPerPar[i]/GEV , _CalHitEPerParAfterSubtraction[i]/GEV);
    }
  }
  
  return;
}


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

void xAODPFlowAna :: FillCaloClusterR(const xAOD::CaloClusterContainer* topocluster){
  xAOD::CaloClusterContainer::const_iterator CaloCluster_itr = topocluster->begin();
   xAOD::CaloClusterContainer::const_iterator CaloCluster_end = topocluster->end();
   for( ; CaloCluster_itr != CaloCluster_end; ++CaloCluster_itr  ) {
     int i_clus = std::distance(topocluster->begin(),CaloCluster_itr);
     _CalClusEta.at(i_clus) = (*CaloCluster_itr)->rawEta();
     _CalClusPhi.at(i_clus) = (*CaloCluster_itr)->rawPhi();
//     std::cout<<"input "<<(*CaloCluster_itr)->rawEta()<<","<< (*CaloCluster_itr)->rawPhi()<<std::endl;
     }
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


void xAODPFlowAna::fillEffPurVectorDefault(const xAOD::CaloClusterContainer* topocluster, int i_mcPart, const xAOD::TruthParticleContainer* TruthParticles,
                                           std::vector<double>& full_Efficiency, std::vector<double>& full_Purity, double tketa, double tkphi) {
  xAOD::CaloClusterContainer::const_iterator CaloCluster_itr = topocluster->begin();
  xAOD::CaloClusterContainer::const_iterator CaloCluster_end = topocluster->end();

  for (; CaloCluster_itr != CaloCluster_end; ++CaloCluster_itr) {
    int i_clus = std::distance(topocluster->begin(), CaloCluster_itr);
    if (_CalHitEPerPar.at(i_mcPart) != 0) {
      float Eff = _CalHitEPerClusFromOnePart.at(i_clus * TruthParticles->size() + i_mcPart) / _CalHitEPerPar.at(i_mcPart);
      full_Efficiency.at(i_clus) = Eff;
    } else {
      full_Efficiency.at(i_clus) = -999;
    }
//    Info("LeadCluster", " v_Efficicency.at(%i)  = %.3f ", i_clus, full_Efficiency.at(i_clus));
//    std::cout<<"fillEffPurVectorDefault energy "<<(*CaloCluster_itr)->rawPhi()<<","<<(*CaloCluster_itr)->rawEta()<<","<<(*CaloCluster_itr)->rawE()<<std::endl;
    if (_CalHitEPerClusFromAllPart.at(i_clus) != 0) {
      float Pur = _CalHitEPerClusFromOnePart.at(i_clus * TruthParticles->size() + i_mcPart) / _CalHitEPerClusFromAllPart.at(i_clus);
      full_Purity.at(i_clus) = Pur;
    } else {
      full_Purity.at(i_clus) = -999;
    }
//    Info("LeadCluster", " v_PURITY.at(%i)  = %.3f ", i_clus, full_Purity.at(i_clus));
//    std::cout<<"LeadCluster "<<i_clus<<" eta="<<_CalClusEta.at(i_clus)<<" phi="<<_CalClusPhi.at(i_clus)<<" track "<<tketa<<","<<tkphi<<" R2="<<(tketa-_CalClusEta.at(i_clus))*(tketa-_CalClusEta.at(i_clus)) + (tkphi-_CalClusPhi.at(i_clus))*(tkphi-_CalClusPhi.at(i_clus))<<std::endl;
  } //end loop over cluster
  return ;
}

// Efficiency & putity for cluster matched
void xAODPFlowAna::Calculate_Efficiency_Purity(const xAOD::TruthParticleContainer* TruthParticles, int _n_clusters, const xAOD::CaloClusterContainer* topocluster, const xAOD::CalCellInfoContainer* CalCellInfo_TopoCluster) {
  xAOD::TruthParticleContainer::const_iterator tp_itr = TruthParticles->begin();
  xAOD::TruthParticleContainer::const_iterator tp_end = TruthParticles->end();

  for (; tp_itr != tp_end; ++tp_itr) {
    int i_mcPart = std::distance(TruthParticles->begin(), tp_itr);
    _mc_pos1.at(i_mcPart) = -1;
    _mc_imax.at(i_mcPart) = -1;

    //We rerequire at least 20% of the energy of the true particle
    if(trackQuality(i_mcPart, tp_itr)) continue;

    bool twoClusters(false);
    int pos1(-1), pos2(-1);

    if (m_1to2matching) {

      long int clusterHash1 = _mc_matchedClusterHash.at(i_mcPart).first;
      long int clusterHash2 = _mc_matchedClusterHash.at(i_mcPart).second;

      twoClusters = ((clusterHash2 == -1) ? false : true);
      _full_Efficiency_max.resize(TruthParticles->size());
      _v_Efficiency.resize(3*TruthParticles->size()); // cluster1, cluster2, cluster1+2
      _v_Purity.resize(3); // cluster1, cluster2, cluster1+2
      _v_dRp.resize(4*3); // cluster1, cluster2, leading
      xAOD::CaloClusterContainer::const_iterator CaloCluster_itr = topocluster->begin();
      xAOD::CaloClusterContainer::const_iterator CaloCluster_end = topocluster->end();
      for (; CaloCluster_itr != CaloCluster_end; ++CaloCluster_itr) {
        int i_clus = std::distance(topocluster->begin(), CaloCluster_itr);

        long int hash = (*CaloCluster_itr)->rawE() * 1000 * (*CaloCluster_itr)->rawEta() * 10 * (*CaloCluster_itr)->rawPhi() * 10;
        if (clusterHash1 == hash) {
          pos1 = i_clus;
        }
        if (clusterHash2 != -1 && clusterHash2 == hash) {
          pos2 = i_clus;
        }
      }

      double etaVar(-1), phiVar(-1);
      double etaMean(-1), phiMean(-1);
      CaloCluster_itr = topocluster->begin() + pos1;
      _v_dRp.at(0) = (pos1 == -1 ? -1 : distanceRprime(_mc_etaTrkEM2[i_mcPart], _mc_phiTrkEM2[i_mcPart], CaloCluster_itr, CalCellInfo_TopoCluster, etaVar, phiVar, etaMean, phiMean));
      _mc_dRp_componets[12*i_mcPart+0] = (pos1 == -1 ? -1 : etaMean);
      _mc_dRp_componets[12*i_mcPart+1] = (pos1 == -1 ? -1 : phiMean);
      _mc_dRp_componets[12*i_mcPart+2] = (pos1 == -1 ? -1 : etaVar);
      _mc_dRp_componets[12*i_mcPart+3] = (pos1 == -1 ? -1 : phiVar);
      CaloCluster_itr = topocluster->begin() + pos2;
      _v_dRp.at(1) = (pos2 == -1 ? -1 : distanceRprime(_mc_etaTrkEM2[i_mcPart], _mc_phiTrkEM2[i_mcPart], CaloCluster_itr, CalCellInfo_TopoCluster, etaVar, phiVar, etaMean, phiMean));
      _mc_dRp_componets[12*i_mcPart+4] = (pos2 == -1 ? -1 : etaMean);
      _mc_dRp_componets[12*i_mcPart+5] = (pos2 == -1 ? -1 : phiMean);
      _mc_dRp_componets[12*i_mcPart+6] = (pos2 == -1 ? -1 : etaVar);
      _mc_dRp_componets[12*i_mcPart+7] = (pos2 == -1 ? -1 : phiVar);

      if (_CalHitEPerPar.at(i_mcPart) != 0) {
        float Eff1 = (pos1 == -1 ? -1 : _CalHitEPerClusFromOnePart.at(pos1 * TruthParticles->size() + i_mcPart) / _CalHitEPerPar.at(i_mcPart));
        float Eff2 = (pos2 == -1 ? (pos1 == -1 ? -1 : 0) : _CalHitEPerClusFromOnePart.at(pos2 * TruthParticles->size() + i_mcPart) / _CalHitEPerPar.at(i_mcPart));
        float Effboth =
            (pos2 == -1 ?
                Eff1 :
                (_CalHitEPerClusFromOnePart.at(pos1 * TruthParticles->size() + i_mcPart) + _CalHitEPerClusFromOnePart.at(pos2 * TruthParticles->size() + i_mcPart)) / _CalHitEPerPar.at(
                    i_mcPart));
        _v_Efficiency.at(3*i_mcPart + 0) = (_v_dRp.at(0) < 0 ? -1 : Eff1);
        _v_Efficiency.at(3*i_mcPart + 1) = (_v_dRp.at(1) < 0 ? -1 : Eff2);
        _v_Efficiency.at(3*i_mcPart + 2) = (_v_dRp.at(0) > 0 || _v_dRp.at(1) > 0 ? Effboth : -1);
      }
//        Info("TwoCluster", " v_Efficicency both = %d, %d, %.3f, %.3f, %.3f ", pos1, pos2, v_Efficiency.at(0), v_Efficiency.at(1), v_Efficiency.at(2));

      if (pos1 != -1 && _CalHitEPerClusFromAllPart.at(pos1) != 0) {
        float Pur1 = (pos1 == -1 ? -1 : _CalHitEPerClusFromOnePart.at(pos1 * TruthParticles->size() + i_mcPart) / _CalHitEPerClusFromAllPart.at(pos1));
        float Pur2 = (pos2 == -1 ? (pos1 == -1 ? -1 : 0) : _CalHitEPerClusFromOnePart.at(pos2 * TruthParticles->size() + i_mcPart) / _CalHitEPerClusFromAllPart.at(pos2));
        float Purboth =
            (pos2 == -1 ?
                Pur1 :
                (_CalHitEPerClusFromOnePart.at(pos1 * TruthParticles->size() + i_mcPart) + _CalHitEPerClusFromOnePart.at(pos2 * TruthParticles->size() + i_mcPart)) / (_CalHitEPerClusFromAllPart.at(
                    pos1)
                                                                                                                                                                       + _CalHitEPerClusFromAllPart.at(
                                                                                                                                                                           pos2)));

        _v_Purity.at(0) = Pur1;
        _v_Purity.at(1) = Pur2;
        _v_Purity.at(2) = Purboth;
      }
    }
    _full_Efficiency.resize(_n_clusters);
    _full_Purity.resize(_n_clusters);
    fillEffPurVectorDefault(topocluster, i_mcPart, TruthParticles, _full_Efficiency, _full_Purity, (*tp_itr)->eta(), (*tp_itr)->phi());

    int imax = fillEffPurHistoDefault(i_mcPart, tp_itr, _full_Efficiency, _full_Purity);
    if (m_1to2matching) {
      xAOD::CaloClusterContainer::const_iterator CaloCluster_itr = topocluster->begin() + imax;
      double etaVar(-1), phiVar(-1);
      double etaMean(-1), phiMean(-1);
      _v_dRp.at(2) = distanceRprime(_mc_etaTrkEM2[i_mcPart], _mc_phiTrkEM2[i_mcPart], CaloCluster_itr, CalCellInfo_TopoCluster, etaVar, phiVar, etaMean, phiMean);
      _mc_dRp_componets[12*i_mcPart+8] = etaMean;
      _mc_dRp_componets[12*i_mcPart+9] = phiMean;
      _mc_dRp_componets[12*i_mcPart+10] = etaVar;
      _mc_dRp_componets[12*i_mcPart+11] = phiVar;

      filldRpHisto(i_mcPart, tp_itr, _v_dRp);
    }

    if (m_1to2matching) {
      const double max_eff = *max_element(_full_Efficiency.begin(), _full_Efficiency.end());
      fillEffPurHistoMatch(i_mcPart, tp_itr, _v_Efficiency, _v_Purity, twoClusters, (pos1 == imax), max_eff);
      Info("CheckMatch", "1.(truth particle - %d, pt: imax == pos1, effmax, effpos1 [LFI]): %.2f: %d == %d, %.3f, %.3f [%d] , pflow index = %d", i_mcPart, _mc_hasEflowTrackPt.at(i_mcPart) / GEV, imax, pos1, _full_Efficiency.at(imax),
           _v_Efficiency.at(3*i_mcPart + 0), _mc_LFI.at(i_mcPart), _mc_hasEflowTrackIndex.at(i_mcPart));

      _mc_pos1.at(i_mcPart) = pos1;
      _mc_imax.at(i_mcPart) = imax;

      double pos1rawEtmp;
      if (pos1 != -1) {
        xAOD::CaloClusterContainer::const_iterator pos1cluster = topocluster->begin() + pos1;
        pos1rawEtmp = (*pos1cluster)->rawE() / GEV;
      } else {
        pos1rawEtmp = -99;
      }

      xAOD::CaloClusterContainer::const_iterator maxcluster = topocluster->begin() + imax;
      Info(
          "CheckMatch",
          "2.trackEtaPhiPtE: (%.3f, %.3f; %.1f, %.1f); leadingEtaPhiE,dRp: (%.3f/%.3f, %.3f/%.3f; %.1f): %.3f; pos1EtaPhiE,dRp_cal_eflow: (%.3f/%.3f, %.3f/%.3f; %.1f): %.3f, %.3f;",
          _mc_etaTrkEM2[i_mcPart], _mc_phiTrkEM2[i_mcPart], _mc_hasEflowTrackPt[i_mcPart] / GEV, _mc_hasEflowTrackP[i_mcPart] / GEV, _mc_dRp_componets.at(12 * i_mcPart + 8),
          _mc_dRp_componets.at(12 * i_mcPart + 10), _mc_dRp_componets.at(12 * i_mcPart + 9), _mc_dRp_componets.at(12 * i_mcPart + 11), (*maxcluster)->rawE() / GEV, _v_dRp.at(2),
          _mc_dRp_componets.at(12 * i_mcPart + 0), _mc_dRp_componets.at(12 * i_mcPart + 2), _mc_dRp_componets.at(12 * i_mcPart + 1), _mc_dRp_componets.at(12 * i_mcPart + 3),
          pos1rawEtmp, _v_dRp.at(0), _mc_RpMatchedCluster1[i_mcPart]);

    }



    unsigned etabin;
    for (unsigned ietabin = 0; ietabin < _etaRange.size(); ++ietabin) {
      if (ietabin == _etaRange.size() - 1 && fabs(_mc_etaTrkEM2[i_mcPart]) > _etaRange.at(ietabin)) {
        etabin = ietabin;
      } else if (fabs(_mc_etaTrkEM2[i_mcPart]) > _etaRange.at(ietabin) && fabs(_mc_etaTrkEM2[i_mcPart]) < _etaRange.at(ietabin + 1)) {
        etabin = ietabin;
      }
    }

    const double max_eff = *max_element(_full_Efficiency.begin(), _full_Efficiency.end());
    int ptbin(_mc_hasEflowTrackPt[i_mcPart] / GEV * 10);
    if (pos1 != -1 && _mc_etaTrkEM2[i_mcPart] > -990) {
      std::string complete_name;
      if (max_eff > 0.5) {
        complete_name = histName(etabin, "h_efficiency5_pt", _etaRange);
        m_H1Dict[complete_name]->Fill(ptbin, _v_Efficiency.at(3 * i_mcPart + 0));
        complete_name = histName(etabin, "h_ntracks5_pt", _etaRange);
        m_H1Dict[complete_name]->Fill(ptbin, 1);
        if (max_eff > 0.9) {
          complete_name = histName(etabin, "h_efficiency9_pt", _etaRange);
          m_H1Dict[complete_name]->Fill(ptbin, _v_Efficiency.at(3 * i_mcPart + 0));
          complete_name = histName(etabin, "h_ntracks9_pt", _etaRange);
          m_H1Dict[complete_name]->Fill(ptbin, 1);
        }
      }
    }

    // Fill NClusters reach 90% efficiency: need full_Efficiency
    int NClusters_09 = getNClustersFor90Eff(i_mcPart, _full_Efficiency);
    fillNClustersFor90Eff(i_mcPart, tp_itr, NClusters_09);
  }


  return;
}

void xAODPFlowAna::fillEffPurHistoMatch(int i_mcPart, xAOD::TruthParticleContainer::const_iterator tp_itr, const std::vector<double>& v_Efficiency,
                                        const std::vector<double>& v_Purity, bool twoClusters, const bool correctMatch, const double max_eff) {
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

      if (max_eff < 0.5) continue;

      std::string complete_name = histName(iptbin, ietabin, "EffClusterboth_total", "", _ptRange, _etaRange);
      m_H1Dict[complete_name]->Fill(v_Efficiency[3*i_mcPart + 2]);
      complete_name = histName(iptbin, ietabin, "SubtractStatus", "", _ptRange, _etaRange);
      m_H1Dict[complete_name]->Fill(_mc_subtractStatus[i_mcPart]);
      complete_name = histName(iptbin, ietabin, "LHED", "", _ptRange, _etaRange);
      m_H1Dict[complete_name]->Fill(_mc_LFI[i_mcPart]);
      complete_name = histName(iptbin, ietabin, "eflowdR1", "", _ptRange, _etaRange);
      m_H1Dict[complete_name]->Fill(_mc_RpMatchedCluster1[i_mcPart]);
      complete_name = histName(iptbin, ietabin, "eflowdR2", "", _ptRange, _etaRange);
      m_H1Dict[complete_name]->Fill(_mc_RpMatchedCluster2[i_mcPart]);
      //      Info("EffClusterboth_total", " v_Efficicency both fill = %.3f ", v_Efficiency[2]);

      // CellLevelMatching
      if (_mc_subtractStatus[i_mcPart] == 1) {
        std::string complete_name = histName(iptbin, ietabin, "EffMatch1", "", _ptRange, _etaRange);
        m_H1Dict[complete_name]->Fill(v_Efficiency[3*i_mcPart + 0]);
        complete_name = histName(iptbin, ietabin, "PurMatch1", "", _ptRange, _etaRange);
        m_H1Dict[complete_name]->Fill(v_Purity[0]);
        complete_name = histName(iptbin, ietabin, "EffMatchboth", "", _ptRange, _etaRange);
        m_H1Dict[complete_name]->Fill(v_Efficiency[3*i_mcPart + 2]);
        complete_name = histName(iptbin, ietabin, "PurMatch2", "", _ptRange, _etaRange);
        m_H1Dict[complete_name]->Fill(v_Purity[1]);

        if (twoClusters) {
          complete_name = histName(iptbin, ietabin, "EffClusterboth_CLS2", "", _ptRange, _etaRange);
          m_H1Dict[complete_name]->Fill(v_Efficiency[3*i_mcPart + 2]);
          m_H1Dict[complete_name]->SetFillColor(4);
        } else {
          complete_name = histName(iptbin, ietabin, "EffClusterboth_CLS1", "", _ptRange, _etaRange);
          m_H1Dict[complete_name]->Fill(v_Efficiency[3*i_mcPart + 2]);
          m_H1Dict[complete_name]->SetFillColor(9);
        }
        complete_name = histName(iptbin, ietabin, "eflowdR1_CLS", "", _ptRange, _etaRange);
        m_H1Dict[complete_name]->Fill(_mc_RpMatchedCluster1[i_mcPart]);
        m_H1Dict[complete_name]->SetFillColor(9);
      } else {
        complete_name = histName(iptbin, ietabin, "EffClusterboth_RSS", "", _ptRange, _etaRange);
        m_H1Dict[complete_name]->Fill(v_Efficiency[3*i_mcPart + 2]);
        m_H1Dict[complete_name]->SetFillColor(2);
        complete_name = histName(iptbin, ietabin, "eflowdR1_RSS", "", _ptRange, _etaRange);
        m_H1Dict[complete_name]->Fill(_mc_RpMatchedCluster1[i_mcPart]);
        m_H1Dict[complete_name]->SetFillColor(9);
      }

      if (correctMatch) {
        complete_name = histName(iptbin, ietabin, "eflowdR1_correct", "", _ptRange, _etaRange);
        m_H1Dict[complete_name]->Fill(_mc_RpMatchedCluster1[i_mcPart]);
        m_H1Dict[complete_name]->SetFillColor(4);
      } else {
        complete_name = histName(iptbin, ietabin, "eflowdR1_wrong", "", _ptRange, _etaRange);
        m_H1Dict[complete_name]->Fill(_mc_RpMatchedCluster1[i_mcPart]);
        m_H1Dict[complete_name]->SetFillColor(2);
      }
    }
  }

}

int xAODPFlowAna::fillEffPurHistoDefault(int i_mcPart, xAOD::TruthParticleContainer::const_iterator tp_itr, const std::vector<double>& full_Efficiency,
                                   const std::vector<double>& full_Purity) {
  // Fill efficiency & purity for cluster matched to histograms
  double max_eff = *max_element(full_Efficiency.begin(), full_Efficiency.end());
  _full_Efficiency_max[i_mcPart] = max_eff;
  int i_max_eff = distance(full_Efficiency.begin(), max_element(full_Efficiency.begin(), full_Efficiency.end()));
//  std::cout<<"fillDefault begin "<<max_eff<<" "<<i_max_eff<<std::endl;

  _mc_hasEflowTrackEtaAtLayer.at(i_mcPart) = (*tp_itr)->eta();
  for (unsigned iptbin = 0; iptbin < _ptRange.size(); ++iptbin) {
    for (unsigned ietabin = 0; ietabin < _etaRange.size(); ++ietabin) {
//      std::cout<<"iptbin="<<iptbin<<" ietabin="<<ietabin<<std::endl;
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

      std::string complete_name = histName(iptbin, ietabin, "EffLeading", "", _ptRange, _etaRange);
//      std::cout<<"before fill max_eff="<<max_eff<<std::endl;
      m_H1Dict[complete_name]->Fill(max_eff);
//      std::cout<<"after fill max_eff="<<max_eff<<std::endl;
      if (full_Efficiency.at(i_max_eff) > 0.5) {
        std::string complete_name = histName(iptbin, ietabin, "PurLeading", "", _ptRange, _etaRange);
        m_H1Dict[complete_name]->Fill(full_Purity.at(i_max_eff));
      } //Purity for those clusters with eff>50%
    }
  }
//  std::cout<<"fillDefault end "<<max_eff<<" "<<i_max_eff<<std::endl;
  return i_max_eff;
}

void xAODPFlowAna::filldRpHisto(int i_mcPart, xAOD::TruthParticleContainer::const_iterator tp_itr, std::vector<double>& v_dRp) {
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

      std::string complete_name = histName(iptbin, ietabin, "dRp_leading", "", _ptRange, _etaRange);
      m_H1Dict[complete_name]->Fill(v_dRp[2]);
      complete_name = histName(iptbin, ietabin, "dRp_1st", "", _ptRange, _etaRange);
      m_H1Dict[complete_name]->Fill(v_dRp[0]);
      complete_name = histName(iptbin, ietabin, "dRp_2nd", "", _ptRange, _etaRange);
      m_H1Dict[complete_name]->Fill(v_dRp[1]);
    }
  }
}



int xAODPFlowAna::getNClustersFor90Eff(int i_mcPart, std::vector<double>& full_Efficiency) {
  //link with cluster index destroyed!
  std::sort(full_Efficiency.begin(), full_Efficiency.end());
  int NClusters_09 = 0;
  unsigned int i(1); double Eff_09(0);
  while (i <= full_Efficiency.size()) {
    Eff_09 += full_Efficiency[full_Efficiency.size() - i];
    if (Eff_09 > 0.9) {
      NClusters_09 = i;
      break;
    } else {
      i++;
    }
  }
  return NClusters_09;
}

void xAODPFlowAna::fillNClustersFor90Eff(int i_mcPart, xAOD::TruthParticleContainer::const_iterator tp_itr, int NClusters_09) {

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

      std::string complete_name = histName(iptbin, ietabin, "NClus_09", "", _ptRange, _etaRange);
      m_H1Dict[complete_name]->Fill(NClusters_09);
    }
  }
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
//  for (unsigned int i=0;i<_clMatchedEflowEcone15.size();i++) Info("SubtractionPerf()", " _clMatchedEflowEcone15 E  = %.2f GeV ", _clMatchedEflowEcone15[i]/GEV);
  
  //Calculate the variables
  xAOD::TruthParticleContainer::const_iterator tp_itr = TruthParticles->begin();
  xAOD::TruthParticleContainer::const_iterator tp_end = TruthParticles->end();
  for( ; tp_itr != tp_end; ++tp_itr ) {
    int i = std::distance(TruthParticles->begin(),tp_itr);
       
    if (_mc_hasEflowTrack[i]==1){
      //if (_mc_hasEflowTrack[i]==1 && (_mc_trueE[i]/(*tp_itr)->pt()*cosh((*tp_itr)->eta()))>0.1){
      if (std::fabs((_pfo_Pt.at(_mc_hasEflowTrackIndex[i])-(*tp_itr)->pt())/(*tp_itr)->pt())<0.05){
	if (std::fabs((*tp_itr)->eta())>0.0 && std::fabs((*tp_itr)->eta())<0.4){
	  hTurnOff_CalHitsOverPt_eta_00_04_hist->Fill((*tp_itr)->pt(), _CalHitEPerPar[i]/((*tp_itr)->pt()*cosh((*tp_itr)->eta())),m_EvtWeight);
	  
	  if (_pfo_hasClusterMatched[_mc_hasEflowTrackIndex[i]]==1 && _pfo_LFI[_mc_hasEflowTrackIndex[i]]!=999){
	    float pull15=(_clMatchedEflowEcone15[_mc_hasEflowTrackIndex[i]]-_pfo_iniEoPexp[_mc_hasEflowTrackIndex[i]])/_pfo_inisigmaEoPexp[_mc_hasEflowTrackIndex[i]];
//	    Info("SubtractionPerf", " pull15  = %.2f GeV ", pull15);
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
  _mc_LFI.clear();
  _mc_hasEflowTrackIndex.clear();
  _mc_hasEflowTrackP.clear();
  _mc_hasEflowTrackPt.clear();
  _mc_hasEflowTrackEtaAtLayer.clear();
  _mc_matchedClusterHash.clear();
  _mc_subtractStatus.clear();
  _mc_RpMatchedCluster1.clear();
  _mc_RpMatchedCluster2.clear();
  _mc_etaTrkEM2.clear();
  _mc_phiTrkEM2.clear();
  _mc_pos1.clear();
  _mc_imax.clear();

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
  _full_Efficiency_max.clear();
  _v_Efficiency.clear();
  _v_Purity.clear();
  _full_Efficiency.clear();
  _full_Purity.clear();
  _v_dRp.clear();
  _mc_dRp_componets.clear();
}
  _mc_MinDeltaREflowTrackPair.clear();
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
  _CalClusEta.clear();
  _CalClusPhi.clear();
  
   
  //_mc_LinkedToTruthJets.clear();
  return;
}

double xAODPFlowAna::distanceRprime(const double tr_eta, const double tr_phi, xAOD::CaloClusterContainer::const_iterator& icluster, const xAOD::CalCellInfoContainer* CalCellInfo_TopoCluster, double& etaVar, double& phiVar, double& etaMean, double& phiMean) {
  if(tr_eta<-100) return -2;

  getClusterVariance(icluster, etaMean, phiMean, etaVar, phiVar, CalCellInfo_TopoCluster);
  double dEta = tr_eta - etaMean;
  double dPhi = fabs(tr_phi - phiMean);
  dPhi = dPhi <= M_PI ? dPhi : 2*M_PI - dPhi;

  double dRp = dEta * dEta / etaVar + dPhi * dPhi / phiVar;
//  Info("TwoCluster", " calculate Rp both = (%.5f-%.5f) / %.5f, ((%.5f-%.5f)=%.5f) / %.5f = %.5f", tr_eta, etaMean, etaVar, tr_phi, phiMean, dPhi, phiVar, dRp);
//  std::cout<<"zhangrui "<<getNCells(icluster, CalCellInfo_TopoCluster)<<" cells; eta "<<(*icluster)->rawEta()<<"-"<<(*icluster)->eta()<<"-"<<etaMean<<" phi "<< (*icluster)->rawPhi()<<"-"<< (*icluster)->phi()<<"-"<<phiMean<<" track tr_eta="<<tr_eta<<" tr_phi="<<tr_phi<<std::endl;

  return dRp;
}

void xAODPFlowAna::getClusterVariance(xAOD::CaloClusterContainer::const_iterator icluster, double& etaMean, double& phiMean, double& etaVar, double& phiVar, const xAOD::CalCellInfoContainer* CalCellInfo_TopoCluster) {
  double m_etaPhiLowerLimit(0.0025);

  /* Sum eta, eta^2, phi and phi^2 of all cells */
  double sumeta = 0;
  double sumeta2 = 0;
  double sumphi = 0;
  double sumphi2 = 0;
  double thisCellPhi;
  int nCells = getNCells(icluster, CalCellInfo_TopoCluster);

    /* Catch empty clusters */
    if (nCells == 0 || nCells == 1) {
//      std::cout << "setCluster()\tWARNING\tEmpty cluster passed!" << std::endl;
      etaVar = m_etaPhiLowerLimit;
      phiVar = m_etaPhiLowerLimit;
      return;
    }
  assert(nCells > 0);

  const std::vector<float> cellEta = getCellEta(icluster, CalCellInfo_TopoCluster, nCells);
  const std::vector<float> cellPhi = getCellPhi(icluster, CalCellInfo_TopoCluster, nCells);

  for (unsigned int icell=0; icell < nCells; ++icell) {
    sumeta += cellEta[icell];
    sumeta2 += cellEta[icell] * cellEta[icell];
    eflowAzimuth tmp;
    tmp.m_value = cellPhi[icell];
    thisCellPhi = tmp.cycle(cellPhi[icell]);
    sumphi += thisCellPhi;
    sumphi2 += thisCellPhi * thisCellPhi;
  }

  /* Calculate mean eta and phi */
  etaMean = sumeta / ((double) nCells);
  phiMean = sumphi / ((double) nCells);

  /* Calculate variance of eta and phi (but don't let them go below the lower limit) */
  double varianceCorrection = (double) nCells / (double) (nCells - 1);
  etaVar = std::max(m_etaPhiLowerLimit, (sumeta2 / (double) nCells - etaMean * etaMean) * varianceCorrection);
  phiVar = std::max(m_etaPhiLowerLimit, (sumphi2 / (double) nCells - phiMean * phiMean) * varianceCorrection);
//  etaVar = m_etaPhiLowerLimit;
//  phiVar = m_etaPhiLowerLimit;
  //  std::cout<<"zhangrui etaVar="<<etaVar<<" =>"<<sumeta2<<"/("<<nCells<<"-"<<etaMean<<"*"<<etaMean<<")  *  "<<varianceCorrection<<std::endl;
//  std::cout<<"zhangrui phiVar="<<phiVar<<" =>"<<sumphi2<<"/("<<nCells<<"-"<<phiMean<<"*"<<phiMean<<")  *  "<<varianceCorrection<<std::endl;

  return;
}

unsigned int xAODPFlowAna::getNCells(xAOD::CaloClusterContainer::const_iterator icluster, const xAOD::CalCellInfoContainer* CalCellInfoTopoCluster) const {
  unsigned int nCells = 0;
  xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_itr =  CalCellInfoTopoCluster->begin();
  xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_end =  CalCellInfoTopoCluster->end();
  for( ; CalCellInfoTopoCl_itr != CalCellInfoTopoCl_end; ++CalCellInfoTopoCl_itr ) {
    if((*icluster)->rawE()==(*CalCellInfoTopoCl_itr)->clusterRecoEnergy()){
//      std::cout<<"getNCells "<<(*icluster)->rawE()<<" == "<<(*CalCellInfoTopoCl_itr)->clusterRecoEnergy()<<" nCells="<<nCells<<std::endl;
  nCells++;
    }
  }
  return nCells;
}

const std::vector<float> xAODPFlowAna::getCellEta(xAOD::CaloClusterContainer::const_iterator icluster, const xAOD::CalCellInfoContainer* _CalCellInfoTopoCluster, const unsigned int nCells) const{
  std::vector<float> _cellEta(nCells);
  xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_itr =  _CalCellInfoTopoCluster->begin();
  xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_end =  _CalCellInfoTopoCluster->end();
  int index = 0;
  for( ; CalCellInfoTopoCl_itr != CalCellInfoTopoCl_end; ++CalCellInfoTopoCl_itr ) {
    if((*icluster)->rawE()==(*CalCellInfoTopoCl_itr)->clusterRecoEnergy()){
  _cellEta.at(index)=(*CalCellInfoTopoCl_itr)->cellEta();
  index++;
    }
  }
  /* std::cout << "_cellEta size ="<< _cellEta.size()<<std::endl; */
  /* for (unsigned int i=0;i<_cellEta.size();i++) */
  /* std::cout << " *** CellEta["<<i<<"]" << _cellEta[i]; */
  /* std::cout << '\n'; */
  return _cellEta;
}

const std::vector<float>  xAODPFlowAna::getCellPhi(xAOD::CaloClusterContainer::const_iterator icluster, const xAOD::CalCellInfoContainer* _CalCellInfoTopoCluster, const unsigned int nCells) const{
  std::vector<float> _cellPhi(nCells);
  xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_itr =  _CalCellInfoTopoCluster->begin();
  xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_end =  _CalCellInfoTopoCluster->end();
  int index = 0;
  for( ; CalCellInfoTopoCl_itr != CalCellInfoTopoCl_end; ++CalCellInfoTopoCl_itr ) {
    if((*icluster)->rawE()==(*CalCellInfoTopoCl_itr)->clusterRecoEnergy()){
  _cellPhi[index]=(*CalCellInfoTopoCl_itr)->cellPhi();
  index++;
    }
  }
  return _cellPhi;
}
