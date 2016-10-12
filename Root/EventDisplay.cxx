/*
 * EventDisplay.cxx
 *
 *  Created on: Oct 11, 2016
 *      Author: zhangrui
 */

#include <PFlowAna/xAODPFlowAna.h>
#include "TCanvas.h"
#include "TMarker.h"
#include "TLegend.h"
#include "TLegend.h"
#include "TLegend.h"


void xAODPFlowAna::eventDisplay(const xAOD::CaloClusterContainer* topocluster, int EventNumber, int pflowNo) {

  xAOD::CaloClusterContainer::const_iterator CaloCluster_itr = topocluster->begin();
  xAOD::CaloClusterContainer::const_iterator CaloCluster_end = topocluster->end();

  TCanvas* EDCanEtaPhi = new TCanvas ("EDCanEtaPhi","Event Display: Eta-Phi view", 0, 0, 600, 600);
  // declare eta-phi histogram
  char htitle[100];
  sprintf(htitle, "Event #%d, pflow #%d : #eta-#phi view",EventNumber, pflowNo);
  TH2F * h_ED_EtaPhi = new TH2F ("h_ED_EtaPhi", htitle, 2, -2.5, 2.5, 2, -TMath::Pi(), TMath::Pi());
  h_ED_EtaPhi->SetXTitle("#eta");
  h_ED_EtaPhi->SetYTitle("#phi");
  h_ED_EtaPhi->SetStats(kFALSE);
  h_ED_EtaPhi->Draw();


  // draw point
  TMarker* track_mark = new TMarker();

  for (int i_mcPart = 0; i_mcPart < _mc_pos1.size(); ++i_mcPart) {
    std::cout<<"bbbb i_mcPart="<<i_mcPart<<" _mc_hasEflowTrackIndex.at(i_mcPart)="<<_mc_hasEflowTrackIndex.at(i_mcPart)<<std::endl;
    if (_mc_hasEflowTrackIndex.at(i_mcPart) == pflowNo) {
      std::cout << "if " << i_mcPart <<" _mc_hasEflowTrackIndex.at(i_mcPart)==pflowNo "<<_mc_hasEflowTrackIndex.at(i_mcPart)<<"=="<<pflowNo<< " pos=" << _mc_pos1[i_mcPart] << " max=" << _mc_imax[i_mcPart] << " track " << _mc_etaExtra[i_mcPart] << ", " << _mc_phiExtra[i_mcPart]
                << std::endl;
      track_mark->SetMarkerStyle(20);
      track_mark->SetMarkerSize(2);
      track_mark->SetMarkerColor(kRed);
      track_mark->DrawMarker(_mc_etaExtra[i_mcPart], _mc_phiExtra[i_mcPart]);
    } else if (_mc_pos1[i_mcPart] != -1) {
      std::cout << "else if " << i_mcPart << " pos=" << _mc_pos1[i_mcPart] << " max=" << _mc_imax[i_mcPart] << " track " << _mc_etaExtra[i_mcPart] << ", " << _mc_phiExtra[i_mcPart]
                << std::endl;
      track_mark->SetMarkerStyle(29);
      track_mark->SetMarkerSize(1);
      track_mark->SetMarkerColor(6);
      track_mark->DrawMarker(_mc_etaExtra[i_mcPart], _mc_phiExtra[i_mcPart]);
    } else if (_mc_imax[i_mcPart] != -1) {
      std::cout << "else " << i_mcPart << " pos=" << _mc_pos1[i_mcPart] << " max=" << _mc_imax[i_mcPart] << " track " << _mc_etaExtra[i_mcPart] << ", " << _mc_phiExtra[i_mcPart]
                << std::endl;
      track_mark->SetMarkerStyle(29);
      track_mark->SetMarkerSize(1);
      track_mark->SetMarkerColor(7);
      track_mark->DrawMarker(_mc_etaExtra[i_mcPart], _mc_phiExtra[i_mcPart]);
    }

  }

  EDCanEtaPhi->SaveAs("a.eps");
  delete h_ED_EtaPhi;
  delete EDCanEtaPhi;


}
