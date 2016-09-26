/*
 * CollectionTreeHelper.h
 *
 *  Created on: Apr 1, 2014
 *      Author: mahan
 */

#ifndef TrackLayerHelper_h_
#define TrackLayerHelper_h_

#include <cassert>
#include <vector>

#include "xAODPFlow/PFOContainer.h"
#include "xAODCaloEvent/CaloClusterContainer.h"

namespace TrackLayer
{
enum EflowCaloLayer {EMB1, EMB2, EMB3, EME1, EME2, EME3, HEC1, HEC2, HEC3, HEC4, Tile1, Tile2, Tile3};
}
namespace ClusterLayer
{
enum CaloSample {preSamplerB, EMB1, EMB2, EMB3, preSamplerE, EME1, EME2, EME3, HEC0, HEC1, HEC2, HEC3, TileBar0, TileBar1,
   TileBar2, TileGap1,TileGap2, TileGap3, TileExt0, TileExt1, TileExt2};
}

class TrackLayerHelper {
public:
  TrackLayerHelper(const xAOD::PFOContainer* JetETMissNeutralParticleFlowObjects,const xAOD::CaloClusterContainer* topocluster); 
  //TrackLayerHelper(const xAOD::PFO& pfo);
  virtual ~TrackLayerHelper();

  void FillTrackParamLayer(const xAOD::PFOContainer* JetETMissChargedParticleFlowObjects);
  void FillClusterEnergyLayer(const xAOD::CaloClusterContainer* topocluster);
  
  float trackLayerEta(TrackLayer::EflowCaloLayer layer) {
    assert(layer < (int) _trackEta.size());
    return _trackEta[layer];
  }
  float trackLayerPhi(TrackLayer::EflowCaloLayer layer) {
    assert(layer < (int) _trackPhi.size());
    return _trackPhi[layer];
  }
  
  float clusterLayerE(ClusterLayer::CaloSample layer) {
     assert(layer < (int)_clusterE.size());
     return _clusterE[layer];
   }

  void Print();
  
 private:
  std::vector<float> _trackEta;
  std::vector<float> _trackPhi;
  std::vector<float> _clusterE;
  
};

#endif /* COLLECTIONTREEHELPER_H_ */
