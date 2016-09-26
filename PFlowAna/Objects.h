#ifndef Objects_h
#define Objects_h

#include <math.h>
#include <iostream>
#include <vector>

#include "PFMatchInterfaces.h"
#include "eflowUtil.h"
#include "xAODTruth/TruthParticle.h"
#include "xAODCaloEvent/CaloCluster.h"
#include "xAODCalCellInfo/CalCellInfoContainer.h"

class Track : public PFMatch::ITrack {

 public:
  Track(TrackLayerHelper *LayerHelper, int trackIndex);
  virtual ~Track();
  
  double LayerEta(TrackLayer::EflowCaloLayer layer);
  double LayerPhi(TrackLayer::EflowCaloLayer layer);
  eflowEtaPhiPosition etaPhiInLayer(TrackLayer::EflowCaloLayer layer) const;
  void Print();
 private:
  TrackLayerHelper* _LayerHelper;
};




class Cluster : public PFMatch::ICluster {

public:
  Cluster(TrackLayerHelper *LayerHelper,const xAOD::CaloCluster& CaloCalTopoClusters, const xAOD::CalCellInfoContainer* CalCellInfoTopoCluster, int clIndex);
  virtual ~Cluster();

  double eta() const;
  double phi() const;
  double energy() const;

  unsigned int nCells() const; 
  const std::vector<float> cellPhi() const; 
  const std::vector<float> cellEta() const; 
  
  double clEnInLayer(ClusterLayer::CaloSample layer) const;
  
  int getClIndex() const {
    return _clIndex;
  }
  void Print();


 protected:
  
  unsigned int getNCells(const xAOD::CaloCluster& _CaloCalTopoClusters, const xAOD::CalCellInfoContainer* _CalCellInfoTopoCluster) const {
    unsigned int nCells = 0;
    xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_itr =  _CalCellInfoTopoCluster->begin();
    xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_end =  _CalCellInfoTopoCluster->end();
    for( ; CalCellInfoTopoCl_itr != CalCellInfoTopoCl_end; ++CalCellInfoTopoCl_itr ) {
      if(_CaloCalTopoClusters.rawE()==(*CalCellInfoTopoCl_itr)->clusterRecoEnergy()){
	nCells++;
      }
    }
    //std::cout<<"NCells = "<< nCells << std::endl;
    return nCells;
  }
  
  const std::vector<float> getCellEta(const xAOD::CaloCluster& _CaloCalTopoClusters,const xAOD::CalCellInfoContainer* _CalCellInfoTopoCluster) const{
    std::vector<float> _cellEta(getNCells(_CaloCalTopoClusters,_CalCellInfoTopoCluster));
    xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_itr =  _CalCellInfoTopoCluster->begin();
    xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_end =  _CalCellInfoTopoCluster->end();
    int index = 0;
    for( ; CalCellInfoTopoCl_itr != CalCellInfoTopoCl_end; ++CalCellInfoTopoCl_itr ) {
      if(_CaloCalTopoClusters.rawE()==(*CalCellInfoTopoCl_itr)->clusterRecoEnergy()){
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
  
  const std::vector<float>  getCellPhi(const xAOD::CaloCluster& _CaloCalTopoClusters,const xAOD::CalCellInfoContainer* _CalCellInfoTopoCluster) const{
    std::vector<float> _cellPhi(getNCells(_CaloCalTopoClusters,_CalCellInfoTopoCluster));
    xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_itr =  _CalCellInfoTopoCluster->begin();
    xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_end =  _CalCellInfoTopoCluster->end();
    int index = 0;
    for( ; CalCellInfoTopoCl_itr != CalCellInfoTopoCl_end; ++CalCellInfoTopoCl_itr ) {
      if(_CaloCalTopoClusters.rawE()==(*CalCellInfoTopoCl_itr)->clusterRecoEnergy()){
  	_cellPhi[index]=(*CalCellInfoTopoCl_itr)->cellPhi();
  	index++;
      }
    }
    return _cellPhi;
  }
  

 private:
  TrackLayerHelper* _LayerHelper;
  const xAOD::CaloCluster& _CaloCalTopoClusters;
  const xAOD::CalCellInfoContainer* _CalCellInfoTopoCluster; 
  int _clIndex;
  //const std::vector<float> _cellEta;
  //const std::vector<float> _cellPhi;
    
};

#endif
