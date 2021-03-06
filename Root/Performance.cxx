#include <PFlowAna/xAODPFlowAna.h>

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
  _mc_hasEflowTrackEtaAtLayer.resize(TruthParticles->size());
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
  _clMatchedEflow.resize(JetETMissChargedParticleFlowObjects->size());  
  _clMatchedEflowEcone10.resize(JetETMissChargedParticleFlowObjects->size());
  _clMatchedEflowEcone15.resize(JetETMissChargedParticleFlowObjects->size());
  _clMatchedEflowEcone20.resize(JetETMissChargedParticleFlowObjects->size());
  
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
    _pfo_iniEoPexp.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("eflowRec_tracksExpectedEnergyDeposit"); //EExpect
    _pfo_inisigmaEoPexp.at(cpfo_index) = (*cpfo_itr)->auxdata< float >("eflowRec_tracksExpectedEnergyDepositVariance"); //varEExpect
    _pfo_LFI.at(cpfo_index) = (*cpfo_itr)->auxdata< int >("eflowRec_FirstIntLayer"); //FirstIntLayer

    // Info("fill_PFOVectors()", "pfo information: EoPExp  = %.3f, SigmaEoPt = %.3f, LFI = %d ",
    // 	 _pfo_iniEoPexp.at(cpfo_index),
    // 	 _pfo_inisigmaEoPexp.at(cpfo_index),
    // 	  _pfo_LFI.at(cpfo_index));
    
    
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
    
    // Info("execute()", "Cluster matched: E  = %.3f, eta = %.3f, phi = %.3f ",
	 // _pfo_hasClusterMatched_E.at(cpfo_index),
	 // _pfo_hasClusterMatched_Eta.at(cpfo_index),
	 // _pfo_hasClusterMatched_Phi.at(cpfo_index));
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
	  
	  // Cuts and matching
	  // Charged particle with pt > 0 in the central part of the detector
	  // Match tp and cPFO based on MinDeltaR (filled in CalculateMatrix_MinDeltaR)
	  // WIP: Vertex requirements?  
	  if((*cpfo_itr)->charge()!=0 && (*cpfo_itr)->pt()!=0 
	     && fabs((*cpfo_itr)->eta()) < 2.5 
	     && AreBothTracksMatched(tp_index,cpfo_index))
	    // && acos(cos((*cpfo_itr)->phi())*cos((*tp_itr)->phi())+sin((*cpfo_itr)->phi())*sin((*tp_itr)->phi()))<0.02
	    // && (fabs(z0-tp_z0)*sin((*cpfo_itr)->theta())) < 0.2)
	    {
	      
	      _mc_hasEflowTrack.at(tp_index) = 1; //=1 indicates that we have a eflowTrack matched to the mc particle
	      _mc_hasEflowTrackIndex.at(tp_index) = cpfo_index; //say us which eflowObject corresponds for each mc particle (not association = 0)
	      _mc_hasEflowTrackP.at(tp_index) = fabs(1. / ptrk->qOverP());
	      _mc_hasEflowTrackPt.at(tp_index) =  (*cpfo_itr)->pt();
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
      //Info("ComputeCalibHitsPerParticle()", " The MC particle with index %i has a PFO associated with CalHitEPerPar = %.2f GeV and CalHitEPerParAfterSubtraction = %.2f GeV ", _mc_hasEflowTrackIndex[i], _CalHitEPerPar[i]/GEV , _CalHitEPerParAfterSubtraction[i]/GEV);
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

// Efficiency & putity for cluster matched
void xAODPFlowAna::Calculate_Efficiency_Purity(const xAOD::TruthParticleContainer* TruthParticles, int _n_clusters, const xAOD::CaloClusterContainer* topocluster) {
  xAOD::TruthParticleContainer::const_iterator tp_itr = TruthParticles->begin();
  xAOD::TruthParticleContainer::const_iterator tp_end = TruthParticles->end();

  for (; tp_itr != tp_end; ++tp_itr) {
    int i_mcPart = std::distance(TruthParticles->begin(), tp_itr);

    //We rerequire at least 15% of the energy of the true particle
    if(_mc_hasEflowTrack.at(i_mcPart) == 1 && (_CalHitEPerPar.at(i_mcPart)/((*tp_itr)->pt()*cosh((*tp_itr)->eta())))>0.15) {
      std::vector<double> v_Efficiency(_n_clusters);
      std::vector<double> v_Purity(_n_clusters);

      xAOD::CaloClusterContainer::const_iterator CaloCluster_itr = topocluster->begin();
      xAOD::CaloClusterContainer::const_iterator CaloCluster_end = topocluster->end();
      for(; CaloCluster_itr != CaloCluster_end; ++CaloCluster_itr ) {
        int i_clus = std::distance(topocluster->begin(),CaloCluster_itr);

        if (_CalHitEPerPar.at(i_mcPart) != 0) {
          float Eff = _CalHitEPerClusFromOnePart.at(i_clus*TruthParticles->size()+i_mcPart)/_CalHitEPerPar.at(i_mcPart);
          v_Efficiency.at(i_clus) = Eff;
        }
        else {v_Efficiency.at(i_clus) = -999;}

	//       Info("fill_Eff_LeadCluster", " v_Efficicency.at(%i)  = %.3f ",i_clus,v_Efficiency.at(i_clus));

        if (_CalHitEPerClusFromAllPart.at(i_clus) !=0) {
          float Pur = _CalHitEPerClusFromOnePart.at(i_clus*TruthParticles->size()+i_mcPart)/_CalHitEPerClusFromAllPart.at(i_clus);
          v_Purity.at(i_clus) = Pur;
        }
        else {v_Purity.at(i_clus) =-999;}

        //Info("fill_Eff_LeadCluster", " v_PURITY.at(%i)  = %.3f ",i_clus,v_Purity.at(i_clus));
      }  //end loop over cluster

      // Fill efficiency & putity for cluster matched to histograms
      double max_eff = *max_element(v_Efficiency.begin(), v_Efficiency.end());
      double i_max_eff = distance(v_Efficiency.begin(), max_element(v_Efficiency.begin(), v_Efficiency.end()));

      _mc_hasEflowTrackEtaAtLayer.at(i_mcPart) = (*tp_itr)->eta();

      for (unsigned iptbin = 0; iptbin < _ptRange.size(); ++iptbin) {
        for (unsigned ietabin = 0; ietabin < _etaRange.size(); ++ietabin) {

          bool inRegion[2] = {false, false};

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

          if (inRegion[0] && inRegion[1]) {
            std::string complete_name = histName(iptbin, ietabin, "Eff", "", _ptRange, _etaRange);
            m_H1Dict[complete_name]->Fill(max_eff);
            if (v_Efficiency.at(i_max_eff) > 0.5) {
              std::string complete_name = histName(iptbin, ietabin, "Pur", "", _ptRange, _etaRange);
              m_H1Dict[complete_name]->Fill(v_Purity.at(i_max_eff));
            }  //Purity for those clusters with eff>50%

          }

        }
      }
    }
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
// Calculate the energy in a cone R around the CPFO axis

void xAODPFlowAna :: SumClusterE_ConeR(const xAOD::PFOContainer* JetETMissChargedParticleFlowObjects,const xAOD::CaloClusterContainer* topocluster, float DeltaR){ 
  
  xAOD::PFOContainer::const_iterator cpfo_itr = JetETMissChargedParticleFlowObjects->begin();
  xAOD::PFOContainer::const_iterator cpfo_end = JetETMissChargedParticleFlowObjects->end();
  int cpfo_index = 0;
  
  for( ; cpfo_itr != cpfo_end; ++cpfo_itr ) {
    cpfo_index = std::distance(JetETMissChargedParticleFlowObjects->begin(),cpfo_itr);
    
    //WIP: To be remove after subtraction debugging
    // Info("PrintPFOInfo", "Charged PFO %d E = %.2f GeV  pt = %.2f GeV  eta = %.2f  phi = %.2f",
    // 	 cpfo_index, (*cpfo_itr)->e()/GEV,
    // 	 (*cpfo_itr)->pt()/GEV,
    // 	 (*cpfo_itr)->eta(),
    // 	 (*cpfo_itr)->phi());     
    
    
    double EconeR = 0.0;
    
    xAOD::CaloClusterContainer::const_iterator CaloCluster_itr = topocluster->begin();
    xAOD::CaloClusterContainer::const_iterator CaloCluster_end = topocluster->end();
    for( ; CaloCluster_itr != CaloCluster_end; ++CaloCluster_itr  ) {
 
      double dEta = ((*cpfo_itr)->eta()-(*CaloCluster_itr)->eta());
      double dPhi = fabs((*cpfo_itr)->phi()-(*CaloCluster_itr)->phi());
      if(dPhi > M_PI) dPhi = 2*M_PI - dPhi;
      if(sqrt(dEta*dEta + dPhi*dPhi)<=DeltaR){
	EconeR += (*CaloCluster_itr)->rawE();

	//WIP: To be remove after subtraction debugging
	// std::cout<<" DeltaR: "<< sqrt(dEta*dEta + dPhi*dPhi)<<std::endl;
	// Info("PrintPFOInfo", "Associated cluster E = %.2f eta = %.2f  phi = %.2f",
	//      (*CaloCluster_itr)->rawE()/GEV,
	//      (*CaloCluster_itr)->eta(),
	//      (*CaloCluster_itr)->phi()); 
      }
    }
    
    if(AreTheSame(DeltaR, 0.10)) _clMatchedEflowEcone10.at(cpfo_index) = EconeR;
    else if(AreTheSame(DeltaR, 0.15)) _clMatchedEflowEcone15.at(cpfo_index) = EconeR;
    else if(AreTheSame(DeltaR, 0.20)) _clMatchedEflowEcone20.at(cpfo_index) = EconeR;
    else Info("SumClusterE_ConeR() ", "ERROR: DeltaR cone not defined!");

    //WIP: To be remove after subtraction debugging
    //std::cout<<" DeltaR: "<< DeltaR<<" ECone = "<<EconeR/GEV<<std::endl;
  }
  
  return;
}
  
void xAODPFlowAna :: SubtractionPerf(const xAOD::TruthParticleContainer* TruthParticles){ 

  std::vector<int> PullDeltaR= {10, 15, 20}; //**WIP: maybe move to the header and define only at the beginning
  
  int testIndex = 0;
  
  //Calculate the variables
  xAOD::TruthParticleContainer::const_iterator tp_itr = TruthParticles->begin();
  xAOD::TruthParticleContainer::const_iterator tp_end = TruthParticles->end();
  for( ; tp_itr != tp_end; ++tp_itr ) {
    int i = std::distance(TruthParticles->begin(),tp_itr);

    if (_mc_hasEflowTrack[i]!=1) continue;

    //WIP: Copmments hve to be removed after subtraction debugging
    // std::cout<<"Counter = "<< testIndex<<"  tp = "<<i<<" _mc_hasEflowTrack[i]:"<<_mc_hasEflowTrack[i]<<"  mc_hasEflowTrackIndex[i] = "<<_mc_hasEflowTrackIndex[i]<<std::endl;
    //testIndex++;
    // std::cout<<"_CalHitEPerPar[i]="<< _CalHitEPerPar[i]<<" _CalHitEPerPar[i]/(*tp_itr)->pt()*cosh((*tp_itr)->eta()):"<<_CalHitEPerPar[i]/(*tp_itr)->pt()*cosh((*tp_itr)->eta())<<std::endl;
    // std::cout<<"std::fabs((_pfo_Pt.at(_mc_hasEflowTrackIndex[i])-(*tp_itr)->pt())/(*tp_itr)->pt()) =" <<std::fabs((_pfo_Pt.at(_mc_hasEflowTrackIndex[i])-(*tp_itr)->pt())/(*tp_itr)->pt()) <<std::endl;
    // std::cout<< "_pfo_hasClusterMatched[_mc_hasEflowTrackIndex[i]]"<<_pfo_hasClusterMatched[_mc_hasEflowTrackIndex[i]]<<std::endl;
    // std::cout<< "_pfo_LFI[_mc_hasEflowTrackIndex[i]]"<<_pfo_LFI[_mc_hasEflowTrackIndex[i]]<<std::endl;
    
    //Is this cut needed ? 
    if (_CalHitEPerPar[i]/(*tp_itr)->pt()*cosh((*tp_itr)->eta()) < 0.05) continue;
    //std::cout<< "Edeposited cut"<<std::endl;
    
    if (std::fabs((_pfo_Pt.at(_mc_hasEflowTrackIndex[i])-(*tp_itr)->pt())/(*tp_itr)->pt()) > 0.05) continue; 
    // std::cout<< "DeltaPt cut"<<std::endl;

    //shoul be here or after filling h_CalHitsOverPt?
    if (_pfo_hasClusterMatched[_mc_hasEflowTrackIndex[i]]!=1) continue;
    //std::cout<< "hasClusterMatched"<<std::endl;
    if (_pfo_LFI[_mc_hasEflowTrackIndex[i]]==999) continue;
    //std::cout<< "LFI"<<std::endl;
    
    //filled for different DeltaR and eta values
    for (unsigned i_eta = 0; i_eta<_etaRange.size()-1; i_eta++){
      if (std::fabs((*tp_itr)->eta()) > _etaRange.at(i_eta) && std::fabs((*tp_itr)->eta()) <= _etaRange.at(i_eta+1)){
	
	std::string complete_name = histSubName2(i_eta, "h_CalHitsOverPt", _etaRange);
	m_H2Dict[complete_name]->Fill((*tp_itr)->pt(), _CalHitEPerPar[i]/((*tp_itr)->pt()*cosh((*tp_itr)->eta())),m_EvtWeight);
	
	for (unsigned i_R = 0; i_R<PullDeltaR.size(); i_R++){
	  float pull = 0;
	  if (i_R == 0) pull=(_clMatchedEflowEcone10[_mc_hasEflowTrackIndex[i]]-_pfo_iniEoPexp[_mc_hasEflowTrackIndex[i]])/_pfo_inisigmaEoPexp[_mc_hasEflowTrackIndex[i]];
	  else if (i_R == 1) pull=(_clMatchedEflowEcone15[_mc_hasEflowTrackIndex[i]]-_pfo_iniEoPexp[_mc_hasEflowTrackIndex[i]])/_pfo_inisigmaEoPexp[_mc_hasEflowTrackIndex[i]];
	  else if (i_R == 2) pull=(_clMatchedEflowEcone20[_mc_hasEflowTrackIndex[i]]-_pfo_iniEoPexp[_mc_hasEflowTrackIndex[i]])/_pfo_inisigmaEoPexp[_mc_hasEflowTrackIndex[i]];
	  else Info("SubtractionPerf() ", "ERROR: DeltaR cone not defined!");
	  
	  //std::cout<< " Eta = "<<std::fabs((*tp_itr)->eta()) <<" DeltaR = "<< i_R <<"  pull = "<<pull << std::endl;
	  
	  complete_name = histSubName(i_R, i_eta, "h_Entries_vs_Pull", PullDeltaR, _etaRange);
	  m_H2Dict[complete_name]->Fill((*tp_itr)->pt(),pull,m_EvtWeight);
	  
	  complete_name = histSubName(i_R, i_eta, "h_CalHitsRemainingOverPt_vs_Pull", PullDeltaR, _etaRange);
	  m_TProfDict[complete_name]->Fill((*tp_itr)->pt(), pull, _CalHitEPerParAfterSubtraction[i]/((*tp_itr)->pt()*cosh((*tp_itr)->eta())),m_EvtWeight);
	}
      }
    }
    // inclusive
    // WIP: Refactoring needed! 
    std::string complete_name = histSubName2(_etaRange.size()-1, "h_CalHitsOverPt", _etaRange);
    m_H2Dict[complete_name]->Fill((*tp_itr)->pt(), _CalHitEPerPar[i]/((*tp_itr)->pt()*cosh((*tp_itr)->eta())),m_EvtWeight);
    for (unsigned i_R = 0; i_R<PullDeltaR.size(); i_R++){
      float pull = 0;
      if (i_R == 0) pull=(_clMatchedEflowEcone10[_mc_hasEflowTrackIndex[i]]-_pfo_iniEoPexp[_mc_hasEflowTrackIndex[i]])/_pfo_inisigmaEoPexp[_mc_hasEflowTrackIndex[i]];
      else if (i_R == 1) pull=(_clMatchedEflowEcone15[_mc_hasEflowTrackIndex[i]]-_pfo_iniEoPexp[_mc_hasEflowTrackIndex[i]])/_pfo_inisigmaEoPexp[_mc_hasEflowTrackIndex[i]];
      else if (i_R == 2) pull=(_clMatchedEflowEcone20[_mc_hasEflowTrackIndex[i]]-_pfo_iniEoPexp[_mc_hasEflowTrackIndex[i]])/_pfo_inisigmaEoPexp[_mc_hasEflowTrackIndex[i]];
      else Info("SubtractionPerf() ", "ERROR: DeltaR cone not defined!");
      
      complete_name = histSubName(i_R, _etaRange.size()-1, "h_Entries_vs_Pull", PullDeltaR, _etaRange);
      m_H2Dict[complete_name]->Fill((*tp_itr)->pt(),pull,m_EvtWeight);
      
      complete_name = histSubName(i_R, _etaRange.size()-1, "h_CalHitsRemainingOverPt_vs_Pull", PullDeltaR, _etaRange);
      m_TProfDict[complete_name]->Fill((*tp_itr)->pt(), pull, _CalHitEPerParAfterSubtraction[i]/((*tp_itr)->pt()*cosh((*tp_itr)->eta())),m_EvtWeight);
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
  _mc_MinDeltaREflowTrackPair.clear();

  //PFOVectors
  _pfo_Pt.clear();
  _pfo_iniEoPexp.clear();
  _pfo_inisigmaEoPexp.clear();
  _pfo_LFI.clear();
  _pfo_hasClusterMatched.clear();
  _pfo_hasClusterMatched_Index.clear();
  _pfo_hasClusterMatched_Eta.clear();
  _pfo_hasClusterMatched_Phi.clear();
  _pfo_hasClusterMatched_E.clear();
  _clMatchedEflow.clear();
  _clMatchedEflowEcone10.clear();
  _clMatchedEflowEcone15.clear();
  _clMatchedEflowEcone20.clear();


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
