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
#include "TArrow.h"


void xAODPFlowAna::pflowDisplay(int EventNumber, int pflowNo, double etalow, double etahi, double philow, double phihi, const xAOD::CaloClusterContainer* topocluster, const xAOD::CaloClusterContainer* JetETMissCaloClusterObjects) {
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
  int imc(-1);
  for (int i_mcPart = 0; i_mcPart < _mc_pos1.size(); ++i_mcPart) {
    if (_mc_hasEflowTrackIndex.at(i_mcPart) == pflowNo) {
      imc = i_mcPart;
      break;
    }
  }

  /* draw cluster cloud */
  xAOD::CaloClusterContainer::const_iterator CaloCluster_itr = topocluster->begin();
  xAOD::CaloClusterContainer::const_iterator CaloCluster_end = topocluster->end();
  for (; CaloCluster_itr != CaloCluster_end; ++CaloCluster_itr) {
    int i_clus = std::distance(topocluster->begin(), CaloCluster_itr);
    double eta = (*CaloCluster_itr)->rawEta();
    double phi = (*CaloCluster_itr)->rawPhi();
//    double etavar(0.05), phivar(0.05);
    double etavar = sqrt(_calo_EtaVariance[i_clus]);
    double phivar = sqrt(_calo_PhiVariance[i_clus]);
    if (etavar < 0.00001) continue; // too little energy in the cluster
    if (eta < etalow + 0.05 || eta > etahi - 0.05 || phi < philow + 0.05 || phi > phihi - 0.05) continue;
    long int calohash = _calo_hash[i_clus];

    cluster_ellipse->SetX1(eta);
    cluster_ellipse->SetY1(phi);
    cluster_ellipse->SetR1(etavar);
    cluster_ellipse->SetR2(phivar);
    cluster_ellipse->SetFillColor(color);
    cluster_ellipse->SetFillStyle(style);
    cluster_ellipse->DrawClone();
    double dRp = sqrt(pow((eta-_mc_etaExtra[imc])/etavar, 2.) + pow((phi - _mc_phiExtra[imc])/phivar, 2));
    Info("EventDisplay", "Other clusters: (%.3f, %.3f) (%.3f, %.3f) : %.3f", eta, phi, etavar, phivar, dRp);


    if (i_clus == _mc_imax.at(imc)) {
      eta = _mc_dRp_componets.at(12 * imc + 8);
      phi = _mc_dRp_componets.at(12 * imc + 9);
      etavar = sqrt(_mc_dRp_componets.at(12 * imc + 10));
      phivar = sqrt(_mc_dRp_componets.at(12 * imc + 11));
      cluster_ellipse->SetX1(eta);
      cluster_ellipse->SetY1(phi);
      cluster_ellipse->SetR1(etavar);
      cluster_ellipse->SetR2(phivar);
      cluster_ellipse->SetFillStyle(3544);
      cluster_ellipse->SetFillColor(4);
      max_ellipse = (TEllipse*) (cluster_ellipse->Clone());
      double dRp = sqrt(pow((eta-_mc_etaExtra[imc])/etavar, 2.) + pow((phi - _mc_phiExtra[imc])/phivar, 2));
      Info("EventDisplay", "Leading cluster: (%.3f, %.3f) (%.3f, %.3f) : %.3f", eta, phi, etavar, phivar, dRp);
    } else if (i_clus == _mc_pos1.at(imc)) {
      assert(calohash ==_mc_matchedClusterHash.at(imc).first);
      eta = _mc_dRp_componets.at(12 * imc + 0);
      phi = _mc_dRp_componets.at(12 * imc + 1);
      etavar = sqrt(_mc_dRp_componets.at(12 * imc + 2));
      phivar = sqrt(_mc_dRp_componets.at(12 * imc + 3));
      cluster_ellipse->SetX1(eta);
      cluster_ellipse->SetY1(phi);
      cluster_ellipse->SetR1(etavar);
      cluster_ellipse->SetR2(phivar);
      cluster_ellipse->SetFillStyle(3352);
      cluster_ellipse->SetFillColor(5);
      pos1_ellipse = (TEllipse*) (cluster_ellipse->Clone());
      double dRp = sqrt(pow((eta-_mc_etaExtra[imc])/etavar, 2.) + pow((phi - _mc_phiExtra[imc])/phivar, 2));
      Info("EventDisplay", "Matched cluster: (%.3f, %.3f) (%.3f, %.3f) : %.3f", eta, phi, etavar, phivar, dRp);
    }
  }
  if (max_ellipse) {
    max_ellipse->DrawClone(); //kBlue
    TArrow* maxArror = 0;
    if (max_ellipse->GetX1() < etalow) {
      maxArror = new TArrow(0,0.1,max_ellipse->GetY1(), max_ellipse->GetY1());
    }
    if (max_ellipse->GetX1() > etahi) {
      maxArror = new TArrow(0.9,1,max_ellipse->GetY1(), max_ellipse->GetY1());
    }
    if (max_ellipse->GetY1() > phihi) {
      maxArror = new TArrow(max_ellipse->GetX1(), max_ellipse->GetX1(), 0.9, 1);
    }
    if (max_ellipse->GetY1() < philow) {
      maxArror = new TArrow(max_ellipse->GetX1(), max_ellipse->GetX1(), 0., 0.1);
    }
    if (maxArror) {
      maxArror->Draw();
    }
  }
  if (pos1_ellipse)pos1_ellipse->DrawClone(); // kYellow

  TLegend* Legend = new TLegend(0.6, 0.89, 0.9, 0.89);

  /* draw track marker */
  TMarker* track_mark = new TMarker(), *track, *otrack = 0, *wotrack = 0;
  TText* momentum = 0;
  for (int i_mcPart = 0; i_mcPart < _mc_pos1.size(); ++i_mcPart) {

    if (_mc_etaExtra[i_mcPart] < etalow + 0.05 || _mc_etaExtra[i_mcPart] > etahi - 0.05 || _mc_phiExtra[i_mcPart] < philow + 0.05 || _mc_phiExtra[i_mcPart] > phihi - 0.05)
      continue;

    if (_mc_hasEflowTrackIndex.at(i_mcPart) == pflowNo) {

      track_mark->SetMarkerStyle(20);
      track_mark->SetMarkerSize(2);
      track_mark->SetMarkerColor(kRed);
      track_mark->DrawMarker(_mc_etaExtra[i_mcPart], _mc_phiExtra[i_mcPart]);
      momentum = new TText(_mc_etaExtra[i_mcPart], _mc_phiExtra[i_mcPart], Form("%.1f[%d]", _mc_hasEflowTrackPt.at(imc) / GEV, _mc_LFI.at(imc)));
      momentum->SetTextColor(kRed+2);
      momentum->SetTextSize(0.03);
      momentum->Draw();
      Info("EventDisplay", "Extrapolated studied track: (%.3f, %.3f), LFI = %d", _mc_etaExtra[i_mcPart], _mc_phiExtra[i_mcPart], _mc_LFI.at(imc));
      track = (TMarker*)track_mark->Clone();
      Legend->SetY1(Legend->GetY1() - 0.04);
      Legend->AddEntry(track, "Extrapolated track", "p");

      m_H1Dict["h_Extratrack_eta"]->Fill(_mc_etaExtra[i_mcPart]);
      m_H1Dict["h_Extratrack_phi"]->Fill(_mc_etaExtra[i_mcPart]);

    } else if (_mc_pos1[i_mcPart] != -1) {
      track_mark->SetMarkerStyle(29);
      track_mark->SetMarkerSize(1);
      track_mark->SetMarkerColor(6);
      track_mark->DrawMarker(_mc_etaExtra[i_mcPart], _mc_phiExtra[i_mcPart]);
      Info("EventDisplay", "Track w/ matched cluster: (%.3f, %.3f)", _mc_etaExtra[i_mcPart], _mc_phiExtra[i_mcPart]);
      otrack = (TMarker*)track_mark->Clone();
    } else if (_mc_imax[i_mcPart] != -1) {
      track_mark->SetMarkerStyle(29);
      track_mark->SetMarkerSize(1);
      track_mark->SetMarkerColor(7);
      track_mark->DrawMarker(_mc_etaExtra[i_mcPart], _mc_phiExtra[i_mcPart]);
      Info("EventDisplay", "Track w/o matched cluster: (%.3f, %.3f)", _mc_etaExtra[i_mcPart], _mc_phiExtra[i_mcPart]);
      wotrack = (TMarker*)track_mark->Clone();
    } else {
      continue;
    }

  }
  if (max_ellipse) {
    TText* eff_imax = new TText(max_ellipse->GetX1(), max_ellipse->GetY1(), Form("%.2f", _full_Efficiency_max[imc]));
    eff_imax->SetTextSize(0.03);
    eff_imax->Draw();
  }
  if (pos1_ellipse) {
    TText* eff_pos1 = new TText(pos1_ellipse->GetX1(), pos1_ellipse->GetY1(), Form("%.2f", _v_Efficiency[3 * imc + 0]));
    eff_pos1->SetTextSize(0.03);
    eff_pos1->Draw();
  }
  if(momentum) {
    momentum->Draw();
  }

  h_ED_EtaPhi->Draw("axis same");

  Legend->SetFillStyle(0);
  Legend->SetBorderSize(0);
  Legend->SetTextFont(43);
  Legend->SetTextSize(20);
  if(otrack) {
    Legend->SetY1(Legend->GetY1() - 0.04);
    Legend->AddEntry(otrack, "Other tracks", "p");
  }
  if(wotrack) {
    Legend->SetY1(Legend->GetY1() - 0.04);
    Legend->AddEntry(wotrack, "Tracks failed matching", "p");
  }
  if(max_ellipse) {
    Legend->SetY1(Legend->GetY1() - 0.04);
    Legend->AddEntry(max_ellipse, "Leading cluster", "f");
  }
  if(pos1_ellipse) {
    Legend->SetY1(Legend->GetY1() - 0.04);
    Legend->AddEntry(pos1_ellipse, "Matched cluster", "f");
  }
  Legend->Draw();

  system(Form("mkdir -v -p EvtDisplay_%s/", m_folder.c_str()));
  EDCanEtaPhi->SaveAs(Form("EvtDisplay_%s/Evt%d_pflow%d_pt%d.eps", m_folder.c_str(), EventNumber, pflowNo, int(_mc_hasEflowTrackPt.at(imc) / GEV)));
  delete h_ED_EtaPhi;
  delete EDCanEtaPhi;
}

void xAODPFlowAna::eventDisplay(const xAOD::CaloClusterContainer* JetETMissCaloClusterObjects, const xAOD::CaloClusterContainer* topocluster, int EventNumber) {
  assert(JetETMissCaloClusterObjects->size() == topocluster->size());

  for (int i_mcPart = 0; i_mcPart < _mc_pos1.size(); ++i_mcPart) {
    if (_mc_pos1[i_mcPart] == _mc_imax[i_mcPart]) continue;
    if (_mc_etaExtra[i_mcPart] < -990 || _mc_phiExtra[i_mcPart] < -990) continue;

//    if(_mc_hasEflowTrackIndex.at(i_mcPart)!=0) continue;

    int etalow = _mc_etaExtra[i_mcPart] * 10 - 5;
    int etahi = _mc_etaExtra[i_mcPart] * 10 + 5;
    int philow = _mc_phiExtra[i_mcPart] * 10 - 5;
    int phihi = _mc_phiExtra[i_mcPart] * 10 + 5;
    pflowDisplay(m_eventCounter-1, _mc_hasEflowTrackIndex.at(i_mcPart), etalow*0.1, etahi*0.1, philow*0.1, phihi*0.1, topocluster, JetETMissCaloClusterObjects);
  }
}
