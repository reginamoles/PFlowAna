/*
 * PFMatchPositions.cxx
 *
 *  Created on: 25.03.2014
 *      Author: tlodd
 */

#include <cassert>
#include <iostream>

#include <PFlowAna/PFMatchPositions.h>

namespace PFMatch {

/* Track position providers */

IPosition* TrackEtaPhiInFixedLayersProvider::getPosition(const ITrack* track) {
  eflowEtaPhiPosition etaphi = track->etaPhiInLayer(_barrelLayer);
  if (etaphi.getEta() == -999.){
    etaphi = track->etaPhiInLayer(_endcapLayer);
  }
  return new EtaPhi(etaphi);
}


/* Cluster position providers */

IPosition* ClusterPlainEtaPhiProvider::getPosition(const ICluster* cluster) {
  return new EtaPhi(eflowEtaPhiPosition(cluster->eta(), cluster->phi()));
}

const double ClusterGeometricalCenterProvider::m_etaPhiLowerLimit(6.25e-4);

IPosition* ClusterGeometricalCenterProvider::getPosition(const ICluster* cluster) {

  /* Remainder cluster (TODO: What is this about and can this ever happen?!?) */
  if (cluster->eta() == -10) {
    std::cout << "eflowGeometricalClusterCenter::getPosition()\tWARNING\tBad cluster passed!" << std::endl;
    return new EtaPhiWithVariance(eflowEtaPhiPosition(cluster->eta(), cluster->phi()), 0., 0.);
  }

  //Necessitem calcular nCell
  unsigned int nCells = cluster->nCells();
    
  /* Catch empty clusters */
  if (nCells == 0){
    std::cout << "eflowGeometricalClusterCenter::setCluster()\tWARNING\tEmpty cluster passed!" << std::endl;
    return new EtaPhiWithVariance(eflowEtaPhiPosition(cluster->eta(), cluster->phi()), m_etaPhiLowerLimit, m_etaPhiLowerLimit);
  }

  assert(nCells > 0);

  /* Sum eta, eta^2, phi and phi^2 of all cells */
  double sumeta = 0;
  double sumeta2 = 0;
  double sumphi = 0;
  double sumphi2 = 0;
  //Necessitem calcular cellEta & cellPhi
  const std::vector<float> cellEta = cluster->cellEta();
  const std::vector<float> cellPhi = cluster->cellPhi();

  for(unsigned int iCell = 0; iCell < nCells; ++iCell) {
    //std::cout << " CellEta["<<iCell<<"]" << (cluster->cellEta())[iCell]<<std::cout << '\n';
    sumeta += cellEta[iCell];
    sumeta2 += cellEta[iCell]*cellEta[iCell];
    double thisCellPhi =  eflowAzimuth(cellPhi[iCell]).cycle(cluster->phi());
    sumphi += thisCellPhi;
    sumphi2 += thisCellPhi*thisCellPhi;
  }

    
  /* Calculate mean eta and phi */
  double etaMean = sumeta/((double)nCells);
  double phiMean = sumphi/((double)nCells);
  
  /* Calculate variance of eta and phi (but don't let them go below the lower limit) */
  double varianceCorrection = (double)nCells / (double)(nCells-1);
  double etaVariance = std::max(m_etaPhiLowerLimit, (sumeta2/(double)nCells - etaMean*etaMean) * varianceCorrection);
  double phiVariance = std::max(m_etaPhiLowerLimit, (sumphi2/(double)nCells - phiMean*phiMean) * varianceCorrection);
  return new EtaPhiWithVariance(eflowEtaPhiPosition(etaMean, phiMean), etaVariance, phiVariance);
}
  
}
