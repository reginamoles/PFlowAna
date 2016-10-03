#ifndef PFlowAna_xAODPFlowAna_H
#define PFlowAna_xAODPFlowAna_H

#include <EventLoop/Algorithm.h>

#include "xAODRootAccess/TEvent.h"

#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertexContainer.h"
#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/TrackParticleContainer.h"
#include "xAODTracking/TrackParticlexAODHelpers.h"
#include "xAODJet/JetContainer.h"
#include "xAODCaloEvent/CaloCluster.h"
#include "xAODCaloEvent/CaloClusterContainer.h"
#include "xAODPFlow/PFOContainer.h"
#include "xAODPFlow/PFO.h"
#include "xAODCalCellInfo/CalCellInfo.h"
#include "xAODCalCellInfo/CalCellInfoContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODMuon/MuonContainer.h"


// ROOT include(s)
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <TProfile2D.h>


//STL includes
#include <string>
#include <fstream>
#include <iostream>
#include <cstring>
#include <vector>
#include <map>
#include <vector>



class xAODPFlowAna : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
 public:
  // float cutValue;
  
  
  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
  
  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();
  xAODPFlowAna ();
  xAODPFlowAna (bool SinglePionLowPerformanceStudies, bool DijetLowPerformance, bool DijetSubtraction, bool Zmumu, bool matching);
  
 private:
  
  float GEV; //!
  
  //Ranges define in histInitialize
  std::vector<float> _ptRange;//!
  std::vector<float> _etaRange;//!
  
  /* WIP: create a config file to select of this options */
  //Example: https://svnweb.cern.ch/trac/atlasperf/browser/CombPerf/JetETMiss/Run2/Jet/Calibration/JetCalibrationTools/MC/DeriveGSC/trunk/Root/GSC_analysis.cxx
  bool m_SinglePionLowPerformanceStudies; 
  bool m_DijetLowPerformance;
  bool m_DijetSubtraction;
  bool m_Zmumu;
  bool m_1to2matching;

  xAOD::TEvent *m_event;//!
  int m_eventCounter; //!

  double m_EvtWeight; //!
  // Tree *myTree; //!
  // TH1 *myHist; //!


  //CutFlow variables (worik in progress)
  int m_select; //!
  int m_trigger; //!
  int m_number; //!
  int m_jvt; //!
  int m_angle; //!
  int m_nojetsafterfilter; //!
  int m_jetpt; //!


  
  //----------------------------------
  //  Printing varibles and functions
  //----------------------------------
  float PrintDebug;//!
  
  void PrintTruthInfo(const xAOD::TruthParticleContainer*,const xAOD::TruthVertexContainer*, bool);
  void PrintTrackInfo (const xAOD::TrackParticleContainer*, bool);
  void PrintPFOInfo(const xAOD::PFOContainer*, const xAOD::PFOContainer*, bool);
  void PrintClusterInfo(const xAOD::CaloClusterContainer*, const xAOD::CaloClusterContainer*, bool);
  void PrintCalCellInfo(const xAOD::CalCellInfoContainer* , const xAOD::CalCellInfoContainer*, bool);
  void PrintJetCollections(const xAOD::JetContainer*, const xAOD::JetContainer*, bool);


  //----------------------------------
  //  Utils functions
  //----------------------------------
  void EnsurePhiInMinusPiToPi(double& );
  double deltaPhi(double , double );
    
  //----------------------------
  // Performance studies functions
  //----------------------------
  void resize_tpVectors(const xAOD::TruthParticleContainer*);
  void resize_PFOVectors(const xAOD::PFOContainer*);
  void initialise_PFOVectors(int, int, int);
  void fill_PFOVectors(const xAOD::PFOContainer*);  
  void tp_Selection(const xAOD::TruthParticleContainer* ,const xAOD::PFOContainer*);
  void ComputeCalibHitsPerParticle(const xAOD::CalCellInfoContainer* ,const xAOD::CalCellInfoContainer*, const xAOD::TruthParticleContainer*);
  void ComputeCalibHitsPerCluster(const xAOD::CalCellInfoContainer*, const xAOD::CaloClusterContainer*, int);
  void FillCaloClusterR(const xAOD::CaloClusterContainer* topocluster);
  void Calculate_Efficiency_Purity(const xAOD::TruthParticleContainer*,int,const xAOD::CaloClusterContainer*,const xAOD::CalCellInfoContainer*);
  void FindMatchedClusterIndex(const xAOD::PFOContainer* JetETMissChargedParticleFlowObjects, const xAOD::CaloClusterContainer* topocluster);
  void SubtractionPerf(const xAOD::PFOContainer*,const xAOD::CaloClusterContainer*, const xAOD::TruthParticleContainer*);
  int getNClustersFor90Eff(int i_mcPart, std::vector<double>& full_Efficiency);
  void fillNClustersFor90Eff(int i_mcPart, xAOD::TruthParticleContainer::const_iterator tp_itr, int NClusters_09);

  void clear_PerformanceVectors();

  void PerformanceHistos(); //SinglePions performance --> WIP
  
  //---------------------------------
  //  Create and fill Histograms
  //--------------------------------
  //Subtraction studies
  TH2D *hTurnOff_CalHitsOverPt_eta_00_04_hist;//!
  TProfile2D *hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_04_hist;//!
  TH2D *hTurnOff_Entries_vs_Pull_eta_00_04_hist;//!

  void bookH1DHistogram(std::string, int, float, float);
  
  //Low performance studies
  std::map<std::string, TH1D*> m_H1Dict;//!
  void bookH1DPerformanceHistogram(std::string, std::string, std::vector<float>, std::vector<float>, int, float, float);
 
  TH1F *_R0;//!
  TH1F *_1MinusChargedR;//!
  void fill_RPlus_R0(const xAOD::TruthParticleContainer*);
                     
  
  //------------------------------------------------------
  // Performance studies: vectors to store the selection
  //-------------------------------------------------------
  std::vector<float> _pfo_Pt;//!
  //Layer of first interaction
  std::vector<int> _pfo_LFI;//!
  //Extrapolate track eta
  std::vector<float> _pfo_etaExtra;
  //Extrapolate track phi
  std::vector<float> _pfo_phiExtra;
  //The expected E/p value of the pth charged object
  std::vector<float> _pfo_iniEoPexp;//!
  //The expected width of E/p value of the pth charged object
  std::vector<float> _pfo_inisigmaEoPexp;//!
  //The subtract status of the pth charged object
  std::vector<int> _pfo_SubtractStatus;//!
  //The EtaEMB1 of extract of the pth charged object
  std::vector<float> _pfo_EtaEMB1;//!
  //The PhiEMB1 of extract of the pth charged object
  std::vector<float> _pfo_PhiEMB1;//!
  //The EtaEME1 of extract of the pth charged object
  std::vector<float> _pfo_EtaEME1;//!
  //The PhiEME1 of extract of the pth charged object
  std::vector<float> _pfo_PhiEME1;//!
  //The EtaEMB2 of extract of the pth charged object
  std::vector<float> _pfo_EtaEMB2;//!
  //The PhiEMB2 of extract of the pth charged object
  std::vector<float> _pfo_PhiEMB2;//!
  //The EtaEME2 of extract of the pth charged object
  std::vector<float> _pfo_EtaEME2;//!
  //The PhiEME2 of extract of the pth charged object
  std::vector<float> _pfo_PhiEME2;//!
  //The EtaEMB3 of extract of the pth charged object
  std::vector<float> _pfo_EtaEMB3;//!
  //The PhiEMB3 of extract of the pth charged object
  std::vector<float> _pfo_PhiEMB3;//!
  //The EtaEME3 of extract of the pth charged object
  std::vector<float> _pfo_EtaEME3;//!
  //The PhiEME3 of extract of the pth charged object
  std::vector<float> _pfo_PhiEME3;//!
  //The EtaHEC1 of extract of the pth charged object
  std::vector<float> _pfo_EtaHEC1;//!
  //The PhiHEC1 of extract of the pth charged object
  std::vector<float> _pfo_PhiHEC1;//!
  //The EtaHEC2 of extract of the pth charged object
  std::vector<float> _pfo_EtaHEC2;//!
  //The PhiHEC2 of extract of the pth charged object
  std::vector<float> _pfo_PhiHEC2;//!
  //The EtaHEC3 of extract of the pth charged object
  std::vector<float> _pfo_EtaHEC3;//!
  //The PhiHEC3 of extract of the pth charged object
  std::vector<float> _pfo_PhiHEC3;//!
  //The EtaHEC4 of extract of the pth charged object
  std::vector<float> _pfo_EtaHEC4;//!
  //The PhiHEC4 of extract of the pth charged object
  std::vector<float> _pfo_PhiHEC4;//!
  //The EtaTile1 of extract of the pth charged object
  std::vector<float> _pfo_EtaTile1;//!
  //The PhiTile1 of extract of the pth charged object
  std::vector<float> _pfo_PhiTile1;//!
  //The EtaTile2 of extract of the pth charged object
  std::vector<float> _pfo_EtaTile2;//!
  //The PhiTile2 of extract of the pth charged object
  std::vector<float> _pfo_PhiTile2;//!
  //The EtaTile3 of extract of the pth charged object
  std::vector<float> _pfo_EtaTile3;//!
  //The PhiTile3 of extract of the pth charged object
  std::vector<float> _pfo_PhiTile3;//!
  //The first cluster hash code
  std::vector<long int> _pfo_hashCluster1;//!
  //The second cluster hash code
  std::vector<long int> _pfo_hashCluster2;//!
  //The first cluster EOP of the pth charged object
  std::vector<float> _pfo_EOP1;//!
  //The all clusters EOP of the pth charged object
  std::vector<float> _pfo_EOPTotal;//!
  //The number of CellLevel mathing of the pth charged object
  std::vector<float> _pfo_NMatchedClusterInCellLevelSubtraction;//!
  //The energy of first cluster
  std::vector<float> _pfo_eMatchedCluster1;//!
  //The energy of second cluster
  std::vector<float> _pfo_eMatchedCluster2;//!
  //The Rprime of the first cluster
  std::vector<float> _pfo_RpMatchedCluster1;//!
  //The Rprime of the second cluster
  std::vector<float> _pfo_RpMatchedCluster2;//!

  //Topocluster energy at em scale
  std::vector<float> _topocluster_em_E;//!
  
  //Charge shower subtraction vectors (from Chris code):
  //This is 1 if there is a charged eflow object associated with that mc particle
  std::vector<int> _mc_hasEflowTrack;//! 
  //The index of the eflow charged object associated with the ith mc particle.
   std::vector<int> _mc_hasEflowTrackIndex;//!
   std::vector<float> _mc_hasEflowTrackP;    //!
   std::vector<float> _mc_hasEflowTrackPt;    //!
   std::vector<float> _mc_hasEflowTrackEtaAtLayer;    //!
   std::vector<std::pair<long int, long int>> _mc_matchedClusterHash; //!
   std::vector<int> _mc_subtractStatus;    //!
   std::vector<double> _mc_RpMatchedCluster1;    //!
   std::vector<double> _mc_RpMatchedCluster2;    //!
   std::vector<double> _mc_etaExtra;    //!
   std::vector<double> _mc_phiExtra;    //!


   std::vector<double> _CalHitEPerClusFromOnePart; //!   //calibration energy per cluster from a certain particle
   std::vector<double> _CalHitEPerClusFromAllPart; //!   //calibration energy per cluster from all particles
   std::vector<double> _CalClusEta; //!
   std::vector<double> _CalClusPhi; //!
   std::vector<double> _CalClusEtaVar; //!
   std::vector<double> _CalClusPhiVar; //!

   std::vector< std::pair<int,int> > _mc_MinDeltaREflowTrackPair;//! //indices for tp and cpfo with MinDeltaR
  
   //This is 1 if there is a cluster matched to the CPFO
   std::vector<int> _pfo_hasClusterMatched;//!   
   std::vector<int> _pfo_hasClusterMatched_Index;//!   
   std::vector<float> _pfo_hasClusterMatched_Eta;//!   
   std::vector<float> _pfo_hasClusterMatched_Phi;//!   
   std::vector<float> _pfo_hasClusterMatched_E;//!   
   //The index of the cluster matched to the pth charged eflow object
   std::vector<int> _clMatchedEflow;//!

   std::vector<int> _CalCellInfo_index;//!                   // tell us cell info associated to which mc particle 
   std::vector<double> _CalHitEPerPar; //!                  //calibration hit energy per particle
   std::vector<double> _CalHitEPerParAfterSubtraction; //! //calibration hit energy per particle


   //The sum of all calibration hits in topoclusters for the ith mc particle
   //std::vector<float> _mc_trueE;//! 
   //The sum of all calibration hits in neutral eflow objects for the ith mc particle
   //std::vector<float> _mc_trueEafter;//!  
   
   
   //The sum of cluster energies in a cone of radius 0.1 around the qth charged eflow object
   std::vector<float> _clMatchedEflowEcone10;//!
   //The sum of cluster energies in a cone of radius 0.15 around the qth charged eflow object
   std::vector<float> _clMatchedEflowEcone15;//! 
   //The index of the truth jet the mc particle belongs to
   std::vector<int> _mc_LinkedToTruthJets;//!  

  
   // DeltaR calculation
   void CalculateMatrix_MinDeltaR (const xAOD::TruthParticleContainer*, const xAOD::PFOContainer*, float);
   bool AreBothTracksMatched (int, int);

   /* WIP: Jet matching tool */
  
  //BadJetsScan
  void BadJetsScan(const xAOD::Jet&);
  void MatchJetCollections(const xAOD::JetContainer*, const xAOD::JetContainer*);
  bool HasPFlowJetMatched(const xAOD::Jet&); //return a true is has been matched
  int  WhichPFlowJetMatched(const xAOD::Jet&); //return the index of the PFlowJet matched
  //Zmumu
  bool ZmumuSelection(const xAOD::ElectronContainer*,const xAOD::MuonContainer*); //return a trueif event pass the selection
  void JetRecoil_Zmumu(const xAOD::ElectronContainer*, const xAOD::MuonContainer*, const xAOD::JetContainer*);
  std::string histName(unsigned i_pt, unsigned i_eta, const std::string& name, const std::string& matchScheme, std::vector<float>& PtRange, std::vector<float>& EtaRange);
  void fillEffPurVectorDefault(const xAOD::CaloClusterContainer* topocluster, int i_mcPart, const xAOD::TruthParticleContainer* TruthParticles, std::vector<double>& v_Efficiency,
                               std::vector<double>& v_Purity, double tketa, double tkphi);
  void fillEffPurHistoMatch(int i_mcPart, xAOD::TruthParticleContainer::const_iterator tp_itr, const std::vector<double>& v_Efficiency, const std::vector<double>& v_Purity, bool twoClusters);
  int fillEffPurHistoDefault(int i_mcPart, xAOD::TruthParticleContainer::const_iterator tp_itr, const std::vector<double>& v_Efficiency, const std::vector<double>& v_Purity);
  void filldRpHistoLeading(xAOD::TruthParticleContainer::const_iterator tp_itr, const xAOD::CaloClusterContainer* topocluster, const std::vector<double>& full_Efficiency);

  double distanceRprime(double tr_eta, double tr_phi, xAOD::CaloClusterContainer::const_iterator& cluster, const xAOD::CalCellInfoContainer* CalCellInfo_TopoCluster);
  void getClusterVariance(xAOD::CaloClusterContainer::const_iterator icluster, double& etaVar, double& phiVar, const xAOD::CalCellInfoContainer* CalCellInfo_TopoCluster);
  unsigned int getNCells(xAOD::CaloClusterContainer::const_iterator icluster, const xAOD::CalCellInfoContainer* _CalCellInfoTopoCluster) const;
  const std::vector<float> getCellEta(xAOD::CaloClusterContainer::const_iterator icluster, const xAOD::CalCellInfoContainer* _CalCellInfoTopoCluster, const unsigned int nCells) const;
  const std::vector<float>  getCellPhi(xAOD::CaloClusterContainer::const_iterator icluster, const xAOD::CalCellInfoContainer* _CalCellInfoTopoCluster, const unsigned int nCells) const;
  void filldRpHisto(int i_mcPart, xAOD::TruthParticleContainer::const_iterator tp_itr, std::vector<double>& v_dRp);




public:

  // this is needed to distribute the algorithm to the workers
  ClassDef(xAODPFlowAna, 1);
};

struct eflowAzimuth {
  double m_value;

  double cycle(double phi) {
    double plainDifference = phi - m_value;
    if (plainDifference > M_PI) {
      return m_value + 2.0 * M_PI;
    } else if (plainDifference < -M_PI) {
      return m_value - 2.0 * M_PI;
    } else {
      return m_value;
    }
  }

  void adjustRange() {
    if (m_value <= -M_PI) {
      m_value += (2 * M_PI * floor(-(m_value - M_PI) / (2 * M_PI)));
    } else if (m_value > M_PI) {
      m_value -= (2 * M_PI * floor((m_value + M_PI) / (2 * M_PI)));
    }
  }
};

#endif
