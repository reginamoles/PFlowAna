#ifndef PFlowAna_xAODPFlowAna_H
#define PFlowAna_xAODPFlowAna_H

#include <EventLoop/Algorithm.h>

#include "xAODRootAccess/TEvent.h"

#include "xAODEventInfo/EventInfo.h"
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

// Duplicated events & GRL
#include "GoodRunsLists/GoodRunsListSelectionTool.h"
//#include "EventLoopAlgs/DuplicateChecker.h"

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
  xAODPFlowAna (bool SinglePionLowPerformanceStudies, bool DijetLowPerformance, bool DijetSubtraction, bool Zmumu, std::string matchScheme, bool UseNarrowPtRange, bool UseNarrowEtaRange, bool PrintDebug);


  
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
  std::string m_matchScheme;
  bool m_UseNarrowPtRange;
  bool m_UseNarrowEtaRange;
  bool m_PrintDebug; 
  

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
  void PrintClusterInfo(const xAOD::CaloClusterContainer*, bool);
  void PrintPFOClusterInfo(const xAOD::CaloClusterContainer*, bool);
  void PrintCalCellInfo(const xAOD::CalCellInfoContainer* , const xAOD::CalCellInfoContainer*, bool);
  void PrintJetCollectionInfo(const xAOD::JetContainer*, const xAOD::JetContainer*, bool);
  void PrintElectronInfo(const xAOD::ElectronContainer*, bool);
  void PrintMuonInfo(const xAOD::MuonContainer*, bool);

  //----------------------------------
  //  Utils functions
  //----------------------------------
  bool isGoodDataEvent (const xAOD::EventInfo*, GoodRunsListSelectionTool*);
  void EnsurePhiInMinusPiToPi(double& );
  double deltaPhi(double , double );
  bool AreTheSame(float , float);
    
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
  void Calculate_Efficiency_Purity(const xAOD::TruthParticleContainer*,int,const xAOD::CaloClusterContainer*);
  void SumClusterE_ConeR(const xAOD::PFOContainer*,const xAOD::CaloClusterContainer*, float);
  void SubtractionPerf(const xAOD::TruthParticleContainer*);
    
  void clear_PerformanceVectors();

  void PerformanceHistos(); //SinglePions performance --> WIP

  //----------------------------
  // Zmumu studies
  //----------------------------
  bool ZmumuSelection(const xAOD::ElectronContainer*,const xAOD::MuonContainer*); //return a true if event pass the selection
  void JetRecoil_Zmumu(const xAOD::MuonContainer*, const xAOD::JetContainer*);
  void FillZmumuHistograms(const xAOD::MuonContainer*);
  
  //---------------------------------
  //  Create and fill Histograms
  //--------------------------------
  //Low performance studies
  std::map<std::string, TH1D*> m_H1Dict;//!
  std::map<std::string, TH2D*> m_H2Dict;//!
  void bookH1DPerformanceHistogram(std::string, std::string, std::vector<float>, std::vector<float>, int, float, float);
  std::string histName(unsigned i_pt, unsigned i_eta, const std::string& name, const std::string& matchScheme, std::vector<float>& PtRange, std::vector<float>& EtaRange);
  
  //Subtraction studies
  std::map<std::string, TProfile2D*> m_TProfDict;//!
  void bookSubHistogram (std::string, std::vector<int>, std::vector<float>, int,  const Double_t *, int,  const Double_t *, std::string);
  std::string histSubName(unsigned i_R, unsigned i_eta, const std::string& name, std::vector<int>& DeltaR, std::vector<float>& EtaRange);
  void bookSubHistogram2(std::string, std::vector<float>, int, const Double_t * ,  int,   float, float);
  std::string  histSubName2(unsigned i_eta, const std::string& name,std::vector<float>& EtaRange);
  
  //Zmumu studies
  void bookH1DHistogram(std::string, int, float, float);
  
  //R0 R+ studies
  TH1F *_R0;//!
  TH1F *_1MinusChargedR;//!
  void fill_RPlus_R0(const xAOD::TruthParticleContainer*);
  
  //------------------------------------------------------
  // Performance studies: vectors to store the selection
  //-------------------------------------------------------
  std::vector<float> _pfo_Pt;//!
  //Layer of first interaction
  std::vector<int> _pfo_LFI;//!  
  //The expected E/p value of the pth charged object
  std::vector<float> _pfo_iniEoPexp;//!
  //The expected width of E/p value of the pth charged object
  std::vector<float> _pfo_inisigmaEoPexp;//!
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
   std::vector<double> _CalHitEPerClusFromOnePart; //!   //calibration energy per cluster from a certain particle
   std::vector<double> _CalHitEPerClusFromAllPart; //!   //calibration energy per cluster from all particles
  
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
   
   
   //The sum of cluster energies in a cone of radius 0.1, 0.15 and 0.2 around the qth charged eflow object
   std::vector<float> _clMatchedEflowEcone10;//!
   std::vector<float> _clMatchedEflowEcone15;//!
   std::vector<float> _clMatchedEflowEcone20;//! 
   //The index of the truth jet the mc particle belongs to
   std::vector<int> _mc_LinkedToTruthJets;//!  

  

   /* WIP: Jet matching tool */
  
  //BadJetsScan
  void BadJetsScan(const xAOD::Jet&);
  void MatchJetCollections(const xAOD::JetContainer*, const xAOD::JetContainer*);
  bool HasPFlowJetMatched(const xAOD::Jet&); //return a true is has been matched
  int  WhichPFlowJetMatched(const xAOD::Jet&); //return the index of the PFlowJet matched
 
public:

  // this is needed to distribute the algorithm to the workers
  ClassDef(xAODPFlowAna, 1);
};

#endif
