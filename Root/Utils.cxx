#include <PFlowAna/xAODPFlowAna.h>

//////////////////////////
// Utils
//////////////////////////

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
