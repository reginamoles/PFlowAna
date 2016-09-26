/*
 * PFMatchInterfaces.h
 *
 *  Created on: 03.04.2014
 *      Author: tlodd
 */

#ifndef PFMATCHINTERFACES_H_
#define PFMATCHINTERFACES_H_

#include <vector>
#include <cassert>

#include "eflowUtil.h"
#include "TrackLayerHelper.h"

namespace PFMatch {

 typedef TrackLayer::EflowCaloLayer LayerType;
// #define TrackLayer eflowCalo

/* The track and cluster abstract interface classes */

class ITrack {
protected:
  ITrack() { }
public:
  virtual ~ITrack() { }

  virtual eflowEtaPhiPosition etaPhiInLayer(LayerType layer) const = 0;
};

class ICluster {
protected:
  ICluster() { }
public:
  virtual ~ICluster () { }

  virtual double eta() const = 0;
  virtual double phi() const = 0;

  virtual unsigned int nCells() const = 0;
  virtual const std::vector<float> cellPhi() const = 0;
  virtual const std::vector<float> cellEta() const = 0;
};


/* Position base class */

class IPosition {
protected:
  IPosition() { }
public:
  virtual ~IPosition() { }
};


/* PositionProvider interface */

template <class ObjectType>
class IPositionProvider {
protected:
  IPositionProvider() { }
public:
  virtual ~IPositionProvider() { }

  virtual IPosition* getPosition(const ObjectType* cluster) = 0;
};

typedef IPositionProvider<ICluster> IClusterPositionProvider;
typedef IPositionProvider<ITrack>     ITrackPositionProvider;


/* Distance calculator interface */

class IDistanceCalculator {
protected:
  IDistanceCalculator() { }
public:
  virtual ~IDistanceCalculator() { }

  virtual double distanceBetween(IPosition* position1, IPosition* position2) = 0;
  template<class PositionType>
  PositionType* convertPosition(IPosition* position);
};

template<class PositionType>
inline PositionType* IDistanceCalculator::convertPosition(IPosition* position) {
  PositionType* converted = dynamic_cast<PositionType*>(position);
  assert(converted);
  return converted;
}


/* The distance provider */

class DistanceProvider {
public:
  DistanceProvider(ITrackPositionProvider*     trackPosition,
                   IClusterPositionProvider* clusterPosition,
                   IDistanceCalculator*   distanceCalculator) :
      _trackPosition(trackPosition),
      _clusterPosition(clusterPosition),
      _distanceCalculator(distanceCalculator) { }
  virtual ~DistanceProvider() {
    delete _trackPosition;
    delete _clusterPosition;
    delete _distanceCalculator;
  }

  double distanceBetween(const ITrack* track, const ICluster* cluster) {
    std::cout<<"distanceBetween"<<std::endl; 
    return _distanceCalculator->distanceBetween(_trackPosition->getPosition(track),
                                                _clusterPosition->getPosition(cluster));
  }

private:
  ITrackPositionProvider*     _trackPosition;
  IClusterPositionProvider* _clusterPosition;
  IDistanceCalculator*   _distanceCalculator;
};

}

#endif /* PFMATCHINTERFACES_H_ */
