/*
 * PFMatchDistance.cxx
 *
 *  Created on: 28.03.2014
 *      Author: tlodd
 */

#include <PFlowAna/PFMatchDistance.h>
#include <PFlowAna/PFMatchPositions.h>

namespace PFMatch {

double EtaPhiSqDistanceCalculator::distanceBetween(IPosition* position1,
                                                   IPosition* position2) {
  EtaPhi* etaphi1 = convertPosition<EtaPhi>(position1);
  EtaPhi* etaphi2 = convertPosition<EtaPhi>(position2);

  return etaphi1->dRSq(*etaphi2);
}

double EtaPhiSqSignificanceCalculator::distanceBetween(IPosition* position1,
                                                       IPosition* position2) {
  std::cout<<"EtaPhiSqSignificanceCalculator::distanceBetween"<<std::endl;
  
  EtaPhi* etaphi1 = convertPosition<EtaPhi>(position1);
  EtaPhiWithVariance* etaphi2 = convertPosition<EtaPhiWithVariance>(position2);
  
  std::cout<< "Track: eta="<<etaphi1->getEta()<<std::endl;
  //std::cout<< "Track: phi="<<etaphi1->getPhi()<<std::endl;
  std::cout<< "Cluster:eta ="<<etaphi2->getEta()<<std::endl;
  //std::cout<< "phi="<<etaphi2->getPhi()<<std::endl;
  std::cout<< "Variance:eta ="<<etaphi2->_etaVariance<< "  phi="<<etaphi2->_phiVariance<<std::endl;

  //Has to be checked (should be in some part of the code)
  double EtaVariance = etaphi2->_etaVariance;
  double PhiVariance = etaphi2->_phiVariance;
  if(EtaVariance<0.05) EtaVariance = 0.05;
  if(PhiVariance<0.05) PhiVariance = 0.05;
  
  std::cout<< "Variance:eta ="<<EtaVariance<< "  phi="<<PhiVariance<<std::endl;
  
  double dEta = etaphi1->getEta() - etaphi2->getEta();
  double dPhi = etaphi1->getPhi().getAbsDifference(etaphi2->getPhi());
  
  return dEta*dEta/EtaVariance + dPhi*dPhi/PhiVariance;
}

}
