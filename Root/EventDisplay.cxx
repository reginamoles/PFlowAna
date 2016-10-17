/*
 * EventDisplay.cxx
 *
 *  Created on: Oct 11, 2016
 *      Author: zhangrui
 */

#include <PFlowAna/xAODPFlowAna.h>
#include "TCanvas.h"
#include "TMarker.h"
#include "TEllipse.h"
#include "TText.h"
#include "TLegend.h"

void xAODPFlowAna::pflowDisplay(int EventNumber, int pflowNo, double etalow, double etahi, double philow, double phihi, const xAOD::CaloClusterContainer* topocluster) {
  TCanvas* EDCanEtaPhi = new TCanvas("EDCanEtaPhi", "Event Display: Eta-Phi view", 0, 0, 600, 600);
  /* declare eta-phi histogram */
  char htitle[100];
  sprintf(htitle, "Event #%d, pflow #%d : #eta-#phi view", EventNumber, pflowNo);
  TH2F* h_ED_EtaPhi = new TH2F("h_ED_EtaPhi", htitle, 2, etalow, etahi, 2, philow, phihi);
  h_ED_EtaPhi->SetXTitle("#eta");
  h_ED_EtaPhi->SetYTitle("#phi");
  h_ED_EtaPhi->SetStats(kFALSE);
  h_ED_EtaPhi->Draw();
  /* draw cluster marker */
  TEllipse* cluster_ellipse = new TEllipse();
  TEllipse* max_ellipse = 0;
  TEllipse* pos1_ellipse = 0;
  int color = cluster_ellipse->GetFillColor();
  int style = cluster_ellipse->GetFillStyle();
  int ind(-1);
  for (int i_mcPart = 0; i_mcPart < _mc_pos1.size(); ++i_mcPart) {
    if (_mc_hasEflowTrackIndex.at(i_mcPart) == pflowNo) {
      ind = i_mcPart;
      break;
    }
  }


  xAOD::CaloClusterContainer::const_iterator CaloCluster_itr = topocluster->begin();
  xAOD::CaloClusterContainer::const_iterator CaloCluster_end = topocluster->end();
  for (; CaloCluster_itr != CaloCluster_end; ++CaloCluster_itr) {
    int i_clus = std::distance(topocluster->begin(), CaloCluster_itr);
    double eta = (*CaloCluster_itr)->rawEta();
    double phi = (*CaloCluster_itr)->rawPhi();
    double etavar(0.05), phivar(0.05);
    if (eta < etalow + 0.05 || eta > etahi - 0.05 || phi < philow + 0.05 || phi > phihi - 0.05) continue;

    cluster_ellipse->SetX1(eta);
    cluster_ellipse->SetY1(phi);
    cluster_ellipse->SetR1(etavar);
    cluster_ellipse->SetR2(phivar);
    cluster_ellipse->SetFillColor(color);
    cluster_ellipse->SetFillStyle(style);
    cluster_ellipse->DrawClone();
//    Info("EventDisplay", "Other clusters: (%.3f, %.3f)", eta, phi);
    if (i_clus == _mc_imax.at(ind)) {
      eta = _mc_dRp_componets.at(12 * ind + 8);
      phi = _mc_dRp_componets.at(12 * ind + 9);
      etavar = sqrt(_mc_dRp_componets.at(12 * ind + 10));
      phivar = sqrt(_mc_dRp_componets.at(12 * ind + 11));
      cluster_ellipse->SetX1(eta);
      cluster_ellipse->SetY1(phi);
      cluster_ellipse->SetR1(etavar);
      cluster_ellipse->SetR2(phivar);
      cluster_ellipse->SetFillStyle(3544);
      cluster_ellipse->SetFillColor(4);
      max_ellipse = (TEllipse*) (cluster_ellipse->Clone());
      Info("EventDisplay", "Leading cluster: (%.3f, %.3f) (%.3f, %.3f)", eta, phi, etavar, phivar);
    } else if (i_clus == _mc_pos1.at(ind)) {
      eta = _mc_dRp_componets.at(12 * ind + 0);
      phi = _mc_dRp_componets.at(12 * ind + 1);
      etavar = sqrt(_mc_dRp_componets.at(12 * ind + 2));
      phivar = sqrt(_mc_dRp_componets.at(12 * ind + 3));
      cluster_ellipse->SetX1(eta);
      cluster_ellipse->SetY1(phi);
      cluster_ellipse->SetR1(etavar);
      cluster_ellipse->SetR2(phivar);
      cluster_ellipse->SetFillStyle(3352);
      cluster_ellipse->SetFillColor(5);
      pos1_ellipse = (TEllipse*) (cluster_ellipse->Clone());
      Info("EventDisplay", "Matched cluster: (%.3f, %.3f) (%.3f, %.3f)", eta, phi, etavar, phivar);
    }
  }
  if (max_ellipse)max_ellipse->DrawClone(); //kBlue
  if (pos1_ellipse)pos1_ellipse->DrawClone(); // kYellow



  /* draw track marker */
  TMarker* track_mark = new TMarker();
  TText* momentum = 0;
  for (int i_mcPart = 0; i_mcPart < _mc_pos1.size(); ++i_mcPart) {

    if (_mc_etaExtra[i_mcPart] < etalow + 0.05 || _mc_etaExtra[i_mcPart] > etahi - 0.05 || _mc_phiExtra[i_mcPart] < philow + 0.05 || _mc_phiExtra[i_mcPart] > phihi - 0.05)
      continue;

    if (_mc_hasEflowTrackIndex.at(i_mcPart) == pflowNo) {
      track_mark->SetMarkerStyle(20);
      track_mark->SetMarkerSize(2);
      track_mark->SetMarkerColor(kRed);
      track_mark->DrawMarker(_mc_etaExtra[i_mcPart], _mc_phiExtra[i_mcPart]);
      momentum = new TText(_mc_etaExtra[i_mcPart], _mc_phiExtra[i_mcPart], Form("%.1f", _mc_hasEflowTrackPt.at(ind) / GEV));
      momentum->SetTextColor(kRed+2);
      momentum->SetTextSize(0.03);
      momentum->Draw();
      Info("EventDisplay", "Extrapolated studied track: (%.3f, %.3f)", _mc_etaExtra[i_mcPart], _mc_phiExtra[i_mcPart]);

      m_H1Dict["h_Extratrack_eta"]->Fill(_mc_etaExtra[i_mcPart]);
      m_H1Dict["h_Extratrack_phi"]->Fill(_mc_etaExtra[i_mcPart]);

    } else if (_mc_pos1[i_mcPart] != -1) {
      track_mark->SetMarkerStyle(29);
      track_mark->SetMarkerSize(1);
      track_mark->SetMarkerColor(6);
      track_mark->DrawMarker(_mc_etaExtra[i_mcPart], _mc_phiExtra[i_mcPart]);
      Info("EventDisplay", "Track w/ matched cluster: (%.3f, %.3f)", _mc_etaExtra[i_mcPart], _mc_phiExtra[i_mcPart]);
    } else if (_mc_imax[i_mcPart] != -1) {
      track_mark->SetMarkerStyle(29);
      track_mark->SetMarkerSize(1);
      track_mark->SetMarkerColor(7);
      track_mark->DrawMarker(_mc_etaExtra[i_mcPart], _mc_phiExtra[i_mcPart]);
      Info("EventDisplay", "Track w/o matched cluster: (%.3f, %.3f)", _mc_etaExtra[i_mcPart], _mc_phiExtra[i_mcPart]);
    } else {
      continue;
    }

  }
  if (max_ellipse) {
    TText* eff_imax = new TText(max_ellipse->GetX1(), max_ellipse->GetY1(), Form("%.2f", _full_Efficiency_max[ind]));
    eff_imax->SetTextSize(0.03);
    eff_imax->Draw();
  }
  if (pos1_ellipse) {
    TText* eff_pos1 = new TText(pos1_ellipse->GetX1(), pos1_ellipse->GetY1(), Form("%.2f", _v_Efficiency[3 * ind + 0]));
    eff_pos1->SetTextSize(0.03);
    eff_pos1->Draw();
  }
  if(momentum) momentum->Draw();


  h_ED_EtaPhi->Draw("axis same");
  system(Form("mkdir -v -p EvtDisplay_%s/", m_folder.c_str()));
  EDCanEtaPhi->SaveAs(Form("EvtDisplay_%s/Evt%d_pflow%d_pt%d.eps", m_folder.c_str(), EventNumber, pflowNo, int(_mc_hasEflowTrackPt.at(ind) / GEV)));
  delete h_ED_EtaPhi;
  delete EDCanEtaPhi;
}

void xAODPFlowAna::eventDisplay(const xAOD::CaloClusterContainer* topocluster, int EventNumber) {

  for (int i_mcPart = 0; i_mcPart < _mc_pos1.size(); ++i_mcPart) {
    if (_mc_pos1[i_mcPart] == _mc_imax[i_mcPart]) continue;
    if (_mc_etaExtra[i_mcPart] < -990 || _mc_phiExtra[i_mcPart] < -990) continue;

//    if(_mc_hasEflowTrackIndex.at(i_mcPart)!=36) continue;

    int etalow = _mc_etaExtra[i_mcPart] * 10 - 5;
    int etahi = _mc_etaExtra[i_mcPart] * 10 + 5;
    int philow = _mc_phiExtra[i_mcPart] * 10 - 5;
    int phihi = _mc_phiExtra[i_mcPart] * 10 + 5;
    pflowDisplay(m_eventCounter-1, _mc_hasEflowTrackIndex.at(i_mcPart), etalow*0.1, etahi*0.1, philow*0.1, phihi*0.1, topocluster);
  }
}
