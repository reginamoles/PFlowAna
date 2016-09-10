#include <PFlowAna/xAODPFlowAna.h>
#include <TSystem.h>

//////////////////////////
// Utils
//////////////////////////

// Apply data event/lumi-block cleaning
bool xAODPFlowAna :: isGoodDataEvent (const xAOD::EventInfo* eventInfo, GoodRunsListSelectionTool *m_grl){
  
  // assume the event is good, then do checks and change to false if fails any selections 
  bool goodEvt = true;

  // check for detector imperfections
  if(   (eventInfo->errorState(xAOD::EventInfo::LAr)==xAOD::EventInfo::Error )  ||
	(eventInfo->errorState(xAOD::EventInfo::Tile)==xAOD::EventInfo::Error ) ||
	(eventInfo->errorState(xAOD::EventInfo::SCT)==xAOD::EventInfo::Error )  ||
	(eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18) ) ){
    goodEvt = false;
  }
  
  if(!m_grl->passRunLB(*eventInfo)) goodEvt = false;

  return goodEvt;
}


void xAODPFlowAna :: EnsurePhiInMinusPiToPi(double& phi) {
  phi = fmod(phi, (2*M_PI));
  if (phi < -M_PI) phi += 2*M_PI;
  if (phi > M_PI)  phi -= 2*M_PI;
  return;
}

double xAODPFlowAna :: deltaPhi(double phi1, double phi2) {
  EnsurePhiInMinusPiToPi(phi1);
  EnsurePhiInMinusPiToPi(phi2);
  double dPhi=phi1-phi2;
  if (dPhi>M_PI) dPhi=2*M_PI-dPhi;
  else if(dPhi<-M_PI) dPhi=2*M_PI+dPhi;
  return dPhi;
}


bool xAODPFlowAna ::AreTheSame(float a, float b)
{
  float EPSILON = 0.0001;
  return fabs(a - b) < EPSILON;
  
}

// double xAODPFlowAna :: DeltaR ( xAOD::TruthParticle cPFO, xAOD::CaloCluster cluster){

//   double dEta = ((*cpfo_itr)->eta()-(*CaloCluster_itr)->eta());
//   double dPhi = fabs(part1->phi()- part2->phi());
//   if(dPhi > M_PI) dPhi = 2*M_PI - dPhi;
  
//   return sqrt(dEta*dEta + dPhi*dPhi);
// }
