/*
 * PFMatchPositions.h
 *
 *  Created on: 25.03.2014
 *      Author: tlodd
 */

#ifndef PFMATCHPOSITION_H_
#define PFMATCHPOSITION_H_

#include <vector>

#include "eflowUtil.h"
#include "PFMatchInterfaces.h"

namespace PFMatch {

/* Position classes */

class EtaPhi: public IPosition, public eflowEtaPhiPosition {
public:
  EtaPhi(const eflowEtaPhiPosition& etaphi): eflowEtaPhiPosition(etaphi) { }
  virtual ~EtaPhi() { }
};

class AllLayersEtaPhi: public IPosition {
public:
  AllLayersEtaPhi() { }
  virtual ~AllLayersEtaPhi() {
    unsigned int nLay = _etaphiInLayer.size();
    for(unsigned int iLay = 0; iLay < nLay; ++iLay) { delete _etaphiInLayer[iLay]; }
  }

  std::vector<EtaPhi*> _etaphiInLayer;
};

class EtaPhiWithVariance: public EtaPhi {
public:
  EtaPhiWithVariance(eflowEtaPhiPosition etaphi, double etaVar, double phiVar):
    EtaPhi(etaphi), _etaVariance(etaVar), _phiVariance(phiVar) { }
  virtual ~EtaPhiWithVariance() { }

  double _etaVariance;
  double _phiVariance;
};


/* Cluster position provider classes */

class ClusterPlainEtaPhiProvider: public IClusterPositionProvider {
public:
  ClusterPlainEtaPhiProvider() { }
  virtual ~ClusterPlainEtaPhiProvider() { }

  virtual IPosition* getPosition(const ICluster* cluster);
};

class ClusterGeometricalCenterProvider: public IClusterPositionProvider {
public:
  ClusterGeometricalCenterProvider() { }
  virtual ~ClusterGeometricalCenterProvider() { }

private:
  virtual IPosition* getPosition(const ICluster* cluster);

  static const double m_etaPhiLowerLimit;
};


/* Track position provider classes */

class TrackEtaPhiInFixedLayersProvider: public ITrackPositionProvider {
public:
  TrackEtaPhiInFixedLayersProvider(LayerType barrelLayer, LayerType endcapLayer):
    _barrelLayer(barrelLayer), _endcapLayer(endcapLayer) { }
  virtual ~TrackEtaPhiInFixedLayersProvider() { }

  virtual IPosition* getPosition(const ITrack* track);

private:
  LayerType _barrelLayer;
  LayerType _endcapLayer;
};

}

#endif /* PFMATCHPOSITION_H_ */
