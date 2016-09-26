/*
 * PFMatchDistance.h
 *
 *  Created on: 28.03.2014
 *      Author: tlodd
 */

#ifndef PFMATCHDISTANCE_H_
#define PFMATCHDISTANCE_H_

#include "PFMatchInterfaces.h"

namespace PFMatch {

/* Distance calculator classes */

class EtaPhiSqDistanceCalculator: public IDistanceCalculator {
public:
  EtaPhiSqDistanceCalculator() { }
  virtual ~EtaPhiSqDistanceCalculator() { }

  virtual double distanceBetween(IPosition* position1, IPosition* position2);
};

class EtaPhiSqSignificanceCalculator: public IDistanceCalculator {
public:
  EtaPhiSqSignificanceCalculator() { }
  virtual ~EtaPhiSqSignificanceCalculator() { }

  virtual double distanceBetween(IPosition* position1, IPosition* position2);
};


}

#endif /* PFMATCHDISTANCE_H_ */
