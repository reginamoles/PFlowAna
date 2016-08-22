// ASG status code check
#include <AsgTools/MessageCheck.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <PFlowAna/xAODPFlowAnaEDM.h>
#include <PFlowAna/xAODPFlowAna.h>

#include <iostream>
#include <vector>
#include <utility> // The pair template is defined in the standard header <utility>
#include <algorithm> //min_element


// this is needed to distribute the algorithm to the workers
ClassImp(xAODPFlowAna)



xAODPFlowAna :: xAODPFlowAna ()
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}


xAODPFlowAna :: xAODPFlowAna (bool SinglePionLowPerformanceStudies, bool DijetLowPerformance, bool DijetSubtraction, bool Zmumu)
{
  m_SinglePionLowPerformanceStudies = SinglePionLowPerformanceStudies;
  m_DijetLowPerformance = DijetLowPerformance;
  m_DijetSubtraction = DijetSubtraction;
  m_Zmumu = Zmumu;

}


EL::StatusCode xAODPFlowAna :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  //In order to tell EventLoop that we want to use the xAODRootAccess we initialize the algorithm to use the xAODRootAccess package
  job.useXAOD ();

  ANA_CHECK_SET_TYPE (EL::StatusCode); // set type of return code you are expecting (add to top of each function once)
  ANA_CHECK(xAOD::Init());

    
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode xAODPFlowAna :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.This method gets called before any input files are
  // connected.

  // ** Subtraction Histograms **
  std::cout<<"m_DijetSubtraction = "<<m_DijetSubtraction<<std::endl;
  if(m_DijetSubtraction){
    double TurnOffBins[20]={500.,631.,794.,1000.,1259.,1585.,1995.,2512.,3162.,3981.,5012.,6310.,7943.,10000.,12589.,15849.,19953.,25119.,31623.,40000.};
    double TurnOffBins2[26]={-5.0,-4.0,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0,10.,12.,14.,16.,20.,25.};
    hTurnOff_Entries_vs_Pull_eta_00_04_hist = new TH2D("hTurnOff_Entries_vs_Pull_eta_00_04_hist","hTurnOff_Entries_vs_Pull_eta_00_04_hist",19,TurnOffBins,25,TurnOffBins2);
    wk()->addOutput (hTurnOff_Entries_vs_Pull_eta_00_04_hist);
    hTurnOff_CalHitsOverPt_eta_00_04_hist = new TH2D ("hTurnOff_CalHitsOverPt_eta_00_04_hist","hTurnOff_CalHitsOverPt_eta_00_04_hist",19,TurnOffBins,40,-0.5,1.5);
    wk()->addOutput (hTurnOff_CalHitsOverPt_eta_00_04_hist);
    hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_04_hist =
      new TProfile2D("hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_04_hist","hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_04_hist",19,TurnOffBins,25,TurnOffBins2,"S");
    wk()->addOutput (hTurnOff_CalHitsRemainingOverPt_vs_Pull_eta_00_04_hist);
  }
  


  
  
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode xAODPFlowAna :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
 
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode xAODPFlowAna :: changeInput (bool firstFile)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode xAODPFlowAna :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.
 
  //Called once, before the first event is executed 

  ANA_CHECK_SET_TYPE (EL::StatusCode);

  m_event = wk()->xaodEvent();
  Info("initialize()", "Number of events = %lli", m_event->getEntries() );

  m_store = new xAOD::TStore();

  // Variable initialization
  m_eventCounter = 0; //Count number of events
  PrintDebug = false; //Printing message criteria -->  Should be chosen from the ATestRun

  // Ranges from histograms


  
  // Conversion factors
  GEV = 1000.; //Units
  
  //----------
  // Tools
  //----------
  //Jet cleaning Tool initialized and configured
  m_jetCleaning = new JetCleaningTool("JetCleaning");
  m_jetCleaning->msg().setLevel( MSG::DEBUG ); 
  ANA_CHECK(m_jetCleaning->setProperty( "CutLevel", "LooseBad"));
  ANA_CHECK(m_jetCleaning->setProperty("DoUgly", false));
  ANA_CHECK(m_jetCleaning->initialize());

  
  // JetCalibration tool for EMTopo jets 
  const std::string name = "akt4EMTopoCalibrationTool";
  TString jetAlgo = "AntiKt4EMTopo";  
  TString config = "JES_MC15cRecommendation_May2016.config"; 
  TString calibSeq = "JetArea_Residual_Origin_EtaJES_GSC_Insitu"; 
  bool isData = true; 
  
  m_akt4EMTopoCalibrationTool = new JetCalibrationTool(name);
  ANA_CHECK(m_akt4EMTopoCalibrationTool->setProperty("JetCollection",jetAlgo.Data()));
  ANA_CHECK(m_akt4EMTopoCalibrationTool->setProperty("ConfigFile",config.Data()));
  ANA_CHECK(m_akt4EMTopoCalibrationTool->setProperty("CalibSequence",calibSeq.Data()));
  ANA_CHECK(m_akt4EMTopoCalibrationTool->setProperty("IsData",isData));
  // Initialize the tool
  ANA_CHECK( m_akt4EMTopoCalibrationTool->initializeTool(name));
  
  /*
  // JetCalibration tool for PFlow jets 
  const std::string name = "JetCalibration_PFlow"; //string describing the current thread, for logging
  TString jetAlgo_PFlow = AntiKt4EMTopo;  //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo (see below)
  TString config_PFlow = ; //Path to global config used to initialize the tool (see below)
  TString calibSeq_PFlow = ; //String describing the calibration sequence to apply (see below)
  bool isData = ; //bool describing if the events are data or from simulation
  
  m_jetCalibration_PFlow = new JetCalibrationTool(name);
  ANA_CHECK(m_jetCalibration_PFlow->setProperty("JetCollection",jetAlgo.Data()));
  ANA_CHECK(m_jetCalibration_PFlow->setProperty("ConfigFile",config.Data()));
  ANA_CHECK(m_jetCalibration_PFlow->setProperty("CalibSequence",calibSeq.Data()));
  ANA_CHECK(m_jetCalibration_PFlow->setProperty("IsData",isData));
  // Initialize the tool
  ANA_CHECK(m_jetCalibration->initializeTool(name));
  */

  
  // Configure the JERTool.
  m_JERTool = new JERTool("JERTool");
  //jerTool.msg()->setLevel(MSG::DEBUG);
  ANA_CHECK( m_JERTool->setProperty("PlotFileName", "JetResolution/Prerec2015_xCalib_2012JER_ReducedTo9NP_Plots_v2.root") );
  ANA_CHECK( m_JERTool->setProperty("CollectionName", "AntiKt4EMTopoJets") );
  ANA_CHECK( m_JERTool->initialize() );

  
  // Configure the JERSmearingTool
  m_SmearTool = new  JERSmearingTool("JERSmearingTool");
  //smearTool.msg()->setLevel(MSG::DEBUG);
  ToolHandle<IJERTool> jerHandle(m_JERTool->name());
  ANA_CHECK( m_SmearTool->setProperty("JERTool", jerHandle) );
  ANA_CHECK( m_SmearTool->setProperty("ApplyNominalSmearing", false) );
  ANA_CHECK( m_SmearTool->setProperty("isMC", true) );
  ANA_CHECK( m_SmearTool->setProperty("SystematicMode", "Full") );
  ANA_CHECK( m_SmearTool->initialize() );
    
 return EL::StatusCode::SUCCESS;
}



EL::StatusCode xAODPFlowAna :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  
  //called once, before the first event is executed 
  ANA_CHECK_SET_TYPE (EL::StatusCode);
    
  if( (m_eventCounter % 1) == 0 ){
    Info("execute()", "----------------" );
    Info("execute()", "   Event %i   ", m_eventCounter );
    Info("execute()", "----------------" );
  }
  m_eventCounter++;

 
  //----------------------------
  // Event information
  //--------------------------- 
  const xAOD::EventInfo* eventInfo = 0;
  ANA_CHECK(m_event->retrieve( eventInfo, "EventInfo"));  
  
  // check if the event is data or MC
  bool isMC = false;
  m_EvtWeight = 1.0;
  
  if( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) )
    isMC = true; 
   
  if( isMC ) { m_EvtWeight = eventInfo->mcEventWeight();}
  Info("execute()", "Event number = %llu  Run Number =  %d  Event weight = %.2f  isMC = %s",eventInfo->eventNumber(), eventInfo->runNumber(), m_EvtWeight, (isMC ? "true" : "false"));
  

  if( isMC ){
    //---------------------------
    // Truth particles & vertices
    //---------------------------
    m_TruthParticles = 0;
    ANA_CHECK(m_event->retrieve( m_TruthParticles,"TruthParticles"));
    m_TruthVertices = 0;
    ANA_CHECK(m_event->retrieve(m_TruthVertices,"TruthVertices"));
    PrintTruthInfo(m_TruthParticles, m_TruthVertices, PrintDebug);
    
    //---------------------------
    // CalCellInfo_TopoCluster
    //---------------------------
    m_CalCellInfo_TopoCluster = 0;
    ANA_CHECK(m_event->retrieve(m_CalCellInfo_TopoCluster, "CalCellInfo_TopoCluster"));
    m_CalCellInfo = 0; //CalCellInfo PFO
    ANA_CHECK(m_event->retrieve(m_CalCellInfo, "CalCellInfo"));
    PrintCalCellInfo(m_CalCellInfo_TopoCluster,m_CalCellInfo, PrintDebug);
  }

  //---------------------------
  // Track Collection
  //---------------------------
  m_InDetTrackParticles  = 0;
  ANA_CHECK(m_event->retrieve( m_InDetTrackParticles ,"InDetTrackParticles"));
  PrintTrackInfo(m_InDetTrackParticles,PrintDebug);
  
  //---------------------------
  // cPFO and nPFO
  //---------------------------
  m_JetETMissChargedParticleFlowObjects  = 0;
  ANA_CHECK(m_event->retrieve( m_JetETMissChargedParticleFlowObjects ,"JetETMissChargedParticleFlowObjects"));
  m_JetETMissNeutralParticleFlowObjects = 0;
  ANA_CHECK(m_event->retrieve(m_JetETMissNeutralParticleFlowObjects,"JetETMissNeutralParticleFlowObjects"));
  PrintPFOInfo( m_JetETMissChargedParticleFlowObjects,m_JetETMissNeutralParticleFlowObjects, PrintDebug);

  //---------------------------
  // EMTopoCluster and PFO cluster
  //---------------------------
  m_topocluster = 0;
  ANA_CHECK(m_event->retrieve( m_topocluster, "CaloCalTopoClusters"));
  m_PFOcluster = 0;
//  ANA_CHECK(m_event->retrieve( m_PFOcluster, "PFOClusters_JetETMiss"));
//  PrintClusterInfo(m_topocluster,m_PFOcluster, PrintDebug);
  
  //----------------------------
  // Jet information
  //--------------------------- 
  m_Jets = 0;
  m_PFlowJets = 0;
  ANA_CHECK(m_event->retrieve( m_Jets, "AntiKt4EMTopoJets" ));
  Info("execute()", "  number of jets = %lu", m_Jets->size());
  ANA_CHECK(m_event->retrieve( m_PFlowJets, "AntiKt4EMPFlowJets" ));
  Info("execute()", "  number of PFlow jets = %lu", m_PFlowJets->size());
  PrintJetCollections(m_Jets,m_PFlowJets, true);


  //---------------------------
  // GRL 
  //--------------------------- 
  //Code to be added

  //------------------------------
  // Trigger (Efficiency and SF)
  //------------------------------
  
  //--------------------
  // Pileup Reweighting 
  //--------------------
  //Code to be added

  //-------------------------
  // Tool for muons
  //-------------------------

  //---------------------------
  // Tools for Jets
  //--------------------------- 
  //the order for the tools has to be checked!

  int numGoodJets = 0;

  // Create the new container and its auxiliary store to store the good and calibrated jets .
  xAOD::JetContainer* m_akt4CalibEMTopo = new xAOD::JetContainer();
  xAOD::AuxContainerBase* m_akt4CalibEMTopoAux = new xAOD::AuxContainerBase();
  m_akt4CalibEMTopo->setStore( m_akt4CalibEMTopoAux ); //< Connect the two
  

  xAOD::JetContainer::const_iterator jet_itr =  m_Jets->begin();
  xAOD::JetContainer::const_iterator jet_end =  m_Jets->end();
  for( ; jet_itr != jet_end; ++jet_itr ) {

    Info("Execute () ", "Jet before calibration E = %.2f GeV  pt = %.2f GeV eta = %.2f  phi =  %.2f",
	 (*jet_itr)->e()/GEV,(*jet_itr)->pt()/GEV, (*jet_itr)->eta(), (*jet_itr)->phi());
    
    //Cleaning TOOL
    //Should we remove the whole event or only the jet? Top analyses remove the whole event.
    if( !m_jetCleaning->accept( **jet_itr )) continue; //only keep good clean jets
    numGoodJets++;

    // //Calibration Tool (Jet Calibration)
    xAOD::Jet* jet = new xAOD::Jet();
    m_akt4EMTopoCalibrationTool->calibratedCopy(**jet_itr,jet); //make a calibrated copy, assuming a copy hasn't been made already, alternative is:
    m_akt4CalibEMTopo->push_back(jet); // jet acquires the m_akt4CalibEMTopo auxstore

    Info("Execute () ", "Jet after Calibration E = %.10f GeV  pt = %.10f GeV eta = %.2f  phi =  %.2f",
     	 jet->e()/GEV, jet->pt()/GEV, jet->eta(), jet->phi());
    
    //JER Tool (Jet Energy Resolution)
    double resMC = m_JERTool->getRelResolutionMC(jet);
    double resData = m_JERTool->getRelResolutionData(jet);

    Info("Execute () ","resMC = %.10f  resData = %.10f ", resMC, resData);
	 
    ANA_CHECK(m_SmearTool->applyCorrection(*jet));
    //virtual CP::CorrectionCode applyCorrection(xAOD::Jet& jet);

    Info("Execute () ", "Jet after Smearing E = %.10f GeV  pt = %.10f GeV eta = %.2f  phi =  %.2f",
	 jet->e()/GEV, jet->pt()/GEV, jet->eta(), jet->phi());
       
  }
  

  // Zmumu selection


  //---------------------
  // Performance studies
  //----------------------
  
  if(m_SinglePionLowPerformanceStudies || m_DijetLowPerformance || m_DijetSubtraction){
    resize_tpVectors(m_TruthParticles);
    resize_PFOVectors(m_JetETMissChargedParticleFlowObjects);
    fill_PFOVectors(m_JetETMissChargedParticleFlowObjects);
    //truth particle selection
    tp_Selection(m_TruthParticles,m_JetETMissChargedParticleFlowObjects);
    //associate calibration hits to truth particles
    ComputeCalibHitsPerParticle(m_CalCellInfo_TopoCluster,m_CalCellInfo,m_TruthParticles);
    //subtraction code
    SubtractionPerf(m_JetETMissChargedParticleFlowObjects,m_topocluster, m_TruthParticles); 
  }
  

  /*
    MatchJetCollections(m_Jets, m_PFlowJets); 
    //int numGoodJets = 0;
    //int numBadJets = 0;
    // loop over the jets in the container
    xAOD::JetContainer::const_iterator jet_itr = m_Jets->begin();
    xAOD::JetContainer::const_iterator jet_end = m_Jets->end();
    for( ; jet_itr != jet_end; ++jet_itr ) {
    
    // Bad Jets
    if(!m_jetCleaning->accept( **jet_itr )) {
    BadJetsScan(**jet_itr);
    numBadJets++;
    }
    
    // Good jets 
    if( !m_jetCleaning->accept( **jet_itr )) continue; //only keep good clean jets
    numGoodJets++;
    Info("execute()", " GOOD jet pt = %.2f GeV", ((*jet_itr)->pt()/GEV)); // just to print out something
    }
    
    Info("execute()", "  number of jets = %lu numGoodJets = %i  numBadJets = %i", m_Jets->size(), numGoodJets, numBadJets);
  */


  
  
  
  
  ///////////////////////
  // Clear copy containers
  ////////////////////////
  // Deep copies. Clearing containers deletes contents including AuxStore.
  //if(m_akt4CalibEMTopo) m_akt4CalibEMTopo->clear();
  //m_store->clear();


  ///////////////////////
  // Clear vectors
  ////////////////////////
  clear_PerformanceVectors();


  return EL::StatusCode::SUCCESS;
}



EL::StatusCode xAODPFlowAna :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode xAODPFlowAna :: finalize ()
{
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.
 
  // finalize(): called once, after the final event has completed 
  
  ANA_CHECK_SET_TYPE(EL::StatusCode);
  
  //Remove the jet tool if it has been created
  if( m_jetCleaning ) {
    delete m_jetCleaning;
    m_jetCleaning = 0;
  }

  if( m_akt4EMTopoCalibrationTool ) {
    delete m_akt4EMTopoCalibrationTool;
    m_akt4EMTopoCalibrationTool = 0;
  }

  if( m_JERTool && m_SmearTool ) {
    delete m_JERTool;
    m_JERTool = 0;
    delete m_SmearTool;
    m_SmearTool = 0;
    
  }
  
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode xAODPFlowAna :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  return EL::StatusCode::SUCCESS;
}





void xAODPFlowAna :: BadJetsScan (const xAOD::Jet& jet) {
  
  Info("", "--- Bad Jets Scanning ---");
  Info("BadJetsScan", "jet E = %.2f GeV  pt  = %.2f GeV eta = %.2f  phi =  %.2f",
       jet.e()/GEV, jet.pt()/GEV, jet.phi(), jet.eta());
  
  // if(HasPFlowJetMatched(jet)){
  //   Info("BadJetsScan", "PFlow jet matched  %d E = %.2f GeV  pt  = %.2f GeV eta = %.2f  phi =  %.2f",
  // 	 jet.e()/GEV, jet.pt()/GEV, jet.phi(), jet.eta());
  // }

   
  return;
}



