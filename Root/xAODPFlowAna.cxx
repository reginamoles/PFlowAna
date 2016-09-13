// ASG status code check
#include <AsgTools/MessageCheck.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <PFlowAna/xAODPFlowAnaEDM.h>
#include <PFlowAna/xAODPFlowAna.h>
//#include <PFlowAna/Utils.h>

#include <iostream>
#include <vector>
#include <TH1.h>
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


xAODPFlowAna :: xAODPFlowAna (bool data, bool SinglePionLowPerformanceStudies, bool DijetLowPerformance, bool DijetSubtraction, bool Zmumu, std::string matchScheme, bool UseNarrowPtRange, bool UseNarrowEtaRange, bool PrintDebug)
{
	
  m_data = data;
  m_SinglePionLowPerformanceStudies = SinglePionLowPerformanceStudies;
  m_DijetLowPerformance = DijetLowPerformance;
  m_DijetSubtraction = DijetSubtraction;
  m_Zmumu = Zmumu;
  
  m_matchScheme = matchScheme;
  m_UseNarrowPtRange = UseNarrowPtRange;
  m_UseNarrowEtaRange = UseNarrowEtaRange;
  m_PrintDebug = PrintDebug; 
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
  
  // Binning selected (pt range in GeV)
  if (m_UseNarrowPtRange) _ptRange= {0, 2, 5, 10, 20};
  else _ptRange = {0, 2, 5, 10, 20, 40, 60, 80, 100, 150, 200, 500, 1000};
  
  if (m_UseNarrowEtaRange) _etaRange= {0, 1, 2, 2.5};
  else _etaRange= {0.0, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.5};
  
  
  /* WIP: a directory structure has to be created to store the histograms */
 
  if(m_DijetSubtraction){
    
    //===================================================
    // Track-jet multiplicity and leading jet-track pt 
    //===================================================
    int n_bins = 100; float x_low = 0.; float x_up = 100.;
    bookH1DPerformanceHistogram("h_TrackJetMultiplicity","", _ptRange, _etaRange, n_bins, x_low, x_up);
    bookH1DPerformanceHistogram("h_TrackJetLeadPt","", _ptRange, _etaRange, n_bins, x_low, x_up);
 
    //====================================
    // Subtraction criteria
    //====================================
    // Which is the criteria that should be use for the binning?
    // from Chris: {500.,600.,750.,950.,1200.,1500.,1850.,2250.,2700.,3200.,3750.,4350.,5000.,5700.,6450.,15849.,19953.,25119.,31623.,40000.};
    int n_binX =25;
    double ptTrack_trueRange[26] = {500., 650.,850.,1100., 1400., 1800., 2300., 2950., 3800., 4800., 6100., 7700., 9300., 11300., 13800., 16800., 20800., 25800., 31800.,
					    39800., 50000., 60000., 70000., 80000., 90000., 100000.,};
    int n_binY =25;
    double PullRange[26] = {-5.0,-4.0,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5,3.0,4.0,5.0,6.0,7.0,8.0,10.,12.,14.,16.,20.,25.};
    std::vector<int> PullDeltaR= {10, 15, 20};
    
    bookSubHistogram("h_Entries_vs_Pull", PullDeltaR, _etaRange, n_binX, ptTrack_trueRange,  n_binY, PullRange, "TH2D");
    bookSubHistogram("h_CalHitsRemainingOverPt_vs_Pull", PullDeltaR, _etaRange, n_binX, ptTrack_trueRange,  n_binY, PullRange, "TProfile");
    bookSubHistogram2("h_CalHitsOverPt", _etaRange, n_binX, ptTrack_trueRange, 40, -0.5, 1.5);
         
    //====================================
    // Track-jet density
    //====================================

    
    //====================================
    // Jet pt resolution
    //====================================
    
  }

 
  if(m_SinglePionLowPerformanceStudies || m_DijetLowPerformance || m_DijetSubtraction){
    //====================================
    // Track cluster matching histograms
    //====================================
    int n_bins = 100; float x_low = 0.; float x_up = 10.;
    bookH1DPerformanceHistogram("dR",m_matchScheme, _ptRange, _etaRange, n_bins, x_low, x_up);
    
    //====================================
    // Eficiency and purity
    //====================================
    int n_effbins = 20; float eff_low = 0; float eff_up = 1.0001;
    bookH1DPerformanceHistogram("Eff","",_ptRange, _etaRange, n_effbins, eff_low, eff_up);
    bookH1DPerformanceHistogram("Pur","",_ptRange, _etaRange, n_effbins, eff_low, eff_up);
    
    //====================================
    // Cluster with 90% of energy
    //====================================
    int nClusBin = 8; float nClus_low = 0.5; float nClus_up = 8.5;
    bookH1DPerformanceHistogram("NClus_09","",_ptRange, _etaRange, nClusBin, nClus_low, nClus_up);
    
    //====================================
    // Difference Ecl-Eexp/sigma(Eexp)                                                                                         
    //====================================
    int EResol_bin = 8; float EResol_low = -5; float EResol_up = 5;
    bookH1DPerformanceHistogram("DeltaE","",_ptRange, _etaRange, EResol_bin, EResol_low, EResol_up);
    bookH1DPerformanceHistogram("DeltaE_07","",_ptRange, _etaRange, EResol_bin, EResol_low, EResol_up); 
    
    //====================================
    // R0 and R+
    //====================================
    _R0= new TH1F("_R0","_R0", 100, 0, 1); wk()->addOutput(_R0);
    _1MinusChargedR= new TH1F("_1MinusChargedR","_1MinusChargedR", 100, 0, 1); wk()->addOutput(_1MinusChargedR);
    
    //====================================
    // DeltaR (using PFOCluster)
    //====================================
    // Not clear if it is preperly done
    int ClusMatch_bin = 8; float ClusMatch_low = -5; float ClusMatch_up = 5;
    bookH1DPerformanceHistogram("MatchedClus","",_ptRange, _etaRange, ClusMatch_bin, ClusMatch_low, ClusMatch_up);
    
    //==============================================
    // HistsQuality histograms (from the old code)
    //==============================================
    //ALSO 1-->2 studies
    //==============================================
    // Extrapolator Histograms (from the old code)
    //==============================================
  }

  if(m_Zmumu){
    //==============================================
    // Zmumu code
    //==============================================
    //Binning and edges to be checked with the PFlow paper
    int pt_bin = 20; float pt_low = 0; float pt_up = 400;
    int E_bin = 70; float E_low = -100; float E_up = 600;
    int eta_bin = 50; float eta_low = -5; float eta_up = 5;
    int phi_bin = 40; float phi_low = -4; float phi_up = 4;
    
    //histograms for testing (to be removed after checks!)
    bookH1DHistogram("h_jetPtdirty", pt_bin, pt_low, pt_up);
    bookH1DHistogram("h_jetEdirty", E_bin, E_low, E_up);
    bookH1DHistogram("h_jetEtadirty", eta_bin, eta_low, eta_up);
    bookH1DHistogram("h_jetPhidirty", phi_bin, phi_low, phi_up);
    
    bookH1DHistogram("h_jetPtcorr", pt_bin, pt_low, pt_up);
    bookH1DHistogram("h_jetEcorr", E_bin, E_low, E_up);
    bookH1DHistogram("h_jetEtacorr", eta_bin, eta_low, eta_up);
    bookH1DHistogram("h_jetPhicorr", phi_bin, phi_low, phi_up);
    
    //Counters
    bookH1DHistogram("h_jetcount", 20,0,20);
    bookH1DHistogram("h_jetcount", 3,1,4);
    
    //Final histograms
    
    //Jet histograms
    bookH1DHistogram("h_jetPt", pt_bin, pt_low, pt_up);
    bookH1DHistogram("h_jetE", E_bin, E_low, E_up);
    bookH1DHistogram("h_jetM", 20, 0, 100); //Is this a proper range? 
    bookH1DHistogram("h_jetEta", eta_bin, eta_low, eta_up);
    bookH1DHistogram("h_jetPhi", phi_bin, phi_low, phi_up);
    
    //Muon histograms
    bookH1DHistogram("h_muonPt", pt_bin, pt_low, pt_up);
    bookH1DHistogram("h_muonE", E_bin, E_low, E_up);
    bookH1DHistogram("h_muonM", 20, 0, 100); //Is this a proper range? 
    bookH1DHistogram("h_muonEta", eta_bin, eta_low, eta_up);
    bookH1DHistogram("h_muonPhi", phi_bin, phi_low, phi_up);
    
    //Z distributions
    bookH1DHistogram("h_ZPt", pt_bin, pt_low, pt_up);
    bookH1DHistogram("h_ZE", E_bin, E_low, E_up);
    bookH1DHistogram("h_ZM", 10, 70, 120); //Is this a proper range? 
    bookH1DHistogram("h_ZEta", eta_bin, eta_low, eta_up);
    bookH1DHistogram("h_ZPhi", phi_bin, phi_low, phi_up);
    //Z+jet system
    bookH1DHistogram("h_ZPt_to_JetPt", 20, 0, 5);
    bookH1DHistogram("h_ZPt_to_JetPt_sum", 20, 0, 5);
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
  
  // Conversion factors
  GEV = 1000.; //Units
  
  // *WIP* Counters for Data-MC CutFlow (Christian)
  m_select = 0;
  m_trigger = 0;
  m_number = 0;
  m_jvt = 0;
  m_angle = 0;
  m_jetpt = 0;
  m_nojetsafterfilter = 0;

  
  /* WIP: Should isData be initialized here for the tool? */ 
  
  //----------
  // Tools
  //----------
  //GRL *WIP* check which file has to be used
  m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
  //const char* GRLFilePath = "$ALRB_TutorialData/data15_13TeV.periodAllYear_DetStatus-v73-pro19-08_DQDefects-00-01-02_PHYS_StandardGRL_All_Good_25ns.xml";
  const char* GRLFilePath = "/afs/cern.ch/user/a/atlasdqm/grlgen/All_Good/data16_13TeV.periodAllYear_DetStatus-v82-pro20-12_DQDefects-00-02-04_PHYS_StandardGRL_All_Good_25ns.xml";
  const char* fullGRLFilePath = gSystem->ExpandPathName (GRLFilePath);
  std::vector<std::string> vecStringGRL;
  vecStringGRL.push_back(fullGRLFilePath);
  ANA_CHECK(m_grl->setProperty( "GoodRunsListVec", vecStringGRL));
  ANA_CHECK(m_grl->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
  ANA_CHECK(m_grl->initialize());

  //Jet cleaning Tool initialized and configured
  m_jetCleaning = new JetCleaningTool("JetCleaning");
  m_jetCleaning->msg().setLevel( MSG::DEBUG ); 
  ANA_CHECK(m_jetCleaning->setProperty( "CutLevel", "LooseBad"));
  ANA_CHECK(m_jetCleaning->setProperty("DoUgly", false));
  ANA_CHECK(m_jetCleaning->initialize());
	/*
  // JetCalibration tool for EMTopo jets 
  const std::string name = "akt4EMTopoCalibrationTool";
  TString jetAlgo = "AntiKt4EMTopo";  
  TString config = "JES_MC15cRecommendation_May2016.config"; 
  TString calibSeq = "JetArea_Residual_Origin_EtaJES_GSC_Insitu";
  bool isData = true; 
  //***Two differences with respect to Christian configuration 1-JetArea_Residual_Origin_EtaJES_GSC and 2-isData = false
  
  m_akt4EMTopoCalibrationTool = new JetCalibrationTool(name);
  ANA_CHECK(m_akt4EMTopoCalibrationTool->setProperty("JetCollection",jetAlgo.Data()));
  ANA_CHECK(m_akt4EMTopoCalibrationTool->setProperty("ConfigFile",config.Data()));
  ANA_CHECK(m_akt4EMTopoCalibrationTool->setProperty("CalibSequence",calibSeq.Data()));
  ANA_CHECK(m_akt4EMTopoCalibrationTool->setProperty("IsData",isData));
  // Initialize the tool
  ANA_CHECK( m_akt4EMTopoCalibrationTool->initializeTool(name));
  */
  
  // JetCalibration tool for PFlow jets 
  const std::string name = "JetCalibration_PFlow"; //string describing the current thread, for logging
  TString jetAlgo = "AntiKt4EMPFlow";  //String describing your jet collection, for example AntiKt4EMTopo or AntiKt4LCTopo (see below)
  TString config = "JES_MC15cRecommendation_PFlow_Aug2016.config"; //Path to global config used to initialize the tool (see below)
  TString calibSeq = "JetArea_Residual_EtaJES"; //String describing the calibration sequence to apply (see below)
  bool isData = false; //bool describing if the events are data or from simulation
  
  m_akt4EMPFlowCalibrationTool = new JetCalibrationTool(name);
  ANA_CHECK(m_akt4EMPFlowCalibrationTool->setProperty("JetCollection",jetAlgo.Data()));
  ANA_CHECK(m_akt4EMPFlowCalibrationTool->setProperty("ConfigFile",config.Data()));
  ANA_CHECK(m_akt4EMPFlowCalibrationTool->setProperty("CalibSequence",calibSeq.Data()));
  ANA_CHECK(m_akt4EMPFlowCalibrationTool->setProperty("IsData",isData));
  // Initialize the tool
  ANA_CHECK(m_akt4EMPFlowCalibrationTool->initializeTool(name));
  
  
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

  
  //***Configuration of the tools below have to be checked!!!!
  // Muon calibration and smearing tool
  m_muonCalibrationAndSmearingTool = new CP::MuonCalibrationAndSmearingTool( "MuonCorrectionTool" );
  //m_muonCalibrationAndSmearingTool->msg().setLevel( MSG::DEBUG );
  ANA_CHECK(m_muonCalibrationAndSmearingTool->initialize());
  //Get systematic uncertainties
  const CP::SystematicRegistry& registry = CP::SystematicRegistry::getInstance();
  const CP::SystematicSet& recommendedSystematics = registry.recommendedSystematics(); // get list of recommended systematics
  m_sysList = CP::make_systematics_vector(recommendedSystematics); 
  
  //Electron calibration	
  m_electronCalibrationAndSmearingTool = new CP::EgammaCalibrationAndSmearingTool("ElectronCorrectionTool"); 
  asg::setProperty(m_electronCalibrationAndSmearingTool, "ESModel", "es2016PRE");  
  asg::setProperty(m_electronCalibrationAndSmearingTool, "decorrelationModel", "FULL_v1");  
  ANA_CHECK(m_electronCalibrationAndSmearingTool->initialize());
 
  
  //Isolation Tool
  m_iso = new CP::IsolationSelectionTool( "m_iso" );
  ANA_CHECK( m_iso->setProperty("MuonWP","Gradient") );
  ANA_CHECK( m_iso->setProperty("ElectronWP","Tight") );
  //ANA_CHECK( m_iso->setProperty("PhotonWP","Cone40") );
  ANA_CHECK( m_iso->initialize() );
  
  //Electron ID Tool
  //m_MediumLH = new AsgElectronLikelihoodTool ("MediumLH");
  // // select the working point
  //ANA_CHECK(m_MediumLH->setProperty("WorkingPoint", "MediumLHElectron"));
  //ANA_CHECK(m_MediumLH->initialize());
  
  m_VeryLooseLHElectron = new AsgElectronLikelihoodTool ("VeryLooseLHElectron");
  // select the working point
  ANA_CHECK(m_VeryLooseLHElectron->setProperty("WorkingPoint", "VeryLooseLHElectron"));
  ANA_CHECK(m_VeryLooseLHElectron->initialize());

  
  //Trigger tools
  m_trigConfigTool = new TrigConf::xAODConfigTool("xAODConfigTool"); // gives us access to the meta-data
  ANA_CHECK( m_trigConfigTool->initialize() );
  ToolHandle< TrigConf::ITrigConfigTool > trigConfigHandle( m_trigConfigTool );
  m_trigDecisionTool = new Trig::TrigDecisionTool("TrigDecisionTool");
  ANA_CHECK(m_trigDecisionTool->setProperty( "ConfigTool", trigConfigHandle ) ); // connect the TrigDecisionTool to the ConfigTool
  ANA_CHECK(m_trigDecisionTool->setProperty( "TrigDecisionKey", "xTrigDecision" ) );
  ANA_CHECK(m_trigDecisionTool->initialize() );
  
  //JVT Tool
  m_jetsf = new CP::JetJvtEfficiency("m_jetsf");
  ANA_CHECK(m_jetsf->setProperty("WorkingPoint","Tight") );
  //ANA_CHECK(m_jetsf->setProperty("SFFile","JetJvtEfficiency/JvtSFFile.root") );
  ANA_CHECK(m_jetsf->initialize() ); 

  
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
  m_EventInfo = 0;
  ANA_CHECK(m_event->retrieve( m_EventInfo, "EventInfo"));  
  
  // check if the event is data or MC
  bool isMC = false;
  m_EvtWeight = 1.0;
  
  if( m_EventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ) isMC = true; 
  //info mismatch
  if( isMC ) { m_EvtWeight = m_EventInfo->mcEventWeight();}
  Info("execute()", "Event number = %llu  Run Number =  %d  Event weight = %.2f  isMC = %s",m_EventInfo->eventNumber(), m_EventInfo->runNumber(), m_EvtWeight, (isMC ? "true" : "false"));


  //--------------------------------------------------------------------------------------------------
  // Event cleaning to be applied on data:
  // Cuts defined to remove problematic luminosity blocks (~1 minute of data taking) based on the GRL
  // and individual events that suffer from detector-level or reconstruction problems. 
  //---------------------------------------------------------------------------------------------------

  if(!isMC){ // it's data!
    bool dataEventPasses = isGoodDataEvent (m_EventInfo, m_grl);
    if(!dataEventPasses){ 
		   //Info("execute()", "something wrong in utils");
      return EL::StatusCode::SUCCESS; // go to the next event
    } 
  }
  
  
  //trigger tools: here the trigger chain is chosen
  auto chainGroup = m_trigDecisionTool->getChainGroup("HLT_mu20_iloose_L1MU15, HLT_mu50");
  std::map<std::string,int> triggerCounts;
  int trigger = 0;
  for(auto &trig : chainGroup->getListOfTriggers()) {
    auto cg = m_trigDecisionTool->getChainGroup(trig);
    std::string thisTrig = trig;
    Info( "execute()", "%30s chain passed(1)/failed(0): %d total chain prescale (L1*HLT): %.1f", thisTrig.c_str(), cg->isPassed(), cg->getPrescale() );
    if(cg->isPassed())trigger = 1;
    Info("execute()", "Trigger + %i", trigger);
  } 

  //***Later we will require pass the trigger (not added yet)
  
  
  if( isMC ){
    //---------------------------
    // Truth particles & vertices
    //---------------------------
    m_TruthParticles = 0;
    ANA_CHECK(m_event->retrieve( m_TruthParticles,"TruthParticles"));
    m_TruthVertices = 0;
    ANA_CHECK(m_event->retrieve(m_TruthVertices,"TruthVertices"));
    PrintTruthInfo(m_TruthParticles, m_TruthVertices, m_PrintDebug);

    
    //---------------------------
    // CalCellInfo_TopoCluster
    //---------------------------
    m_CalCellInfo_TopoCluster = 0;
    if(m_SinglePionLowPerformanceStudies || m_DijetLowPerformance || m_DijetSubtraction){
      ANA_CHECK(m_event->retrieve(m_CalCellInfo_TopoCluster, "CalCellInfo_TopoCluster"));
      m_CalCellInfo = 0; //CalCellInfo PFO
      ANA_CHECK(m_event->retrieve(m_CalCellInfo, "CalCellInfo"));
      PrintCalCellInfo(m_CalCellInfo_TopoCluster,m_CalCellInfo, m_PrintDebug);
    }
  }
  
  //---------------------------
  // Track Collection
  //---------------------------
  m_InDetTrackParticles  = 0;
  ANA_CHECK(m_event->retrieve( m_InDetTrackParticles ,"InDetTrackParticles"));
  PrintTrackInfo(m_InDetTrackParticles,m_PrintDebug);
  
  //---------------------------
  // cPFO and nPFO
  //---------------------------
  m_JetETMissChargedParticleFlowObjects  = 0;
  ANA_CHECK(m_event->retrieve( m_JetETMissChargedParticleFlowObjects ,"JetETMissChargedParticleFlowObjects"));
  m_JetETMissNeutralParticleFlowObjects = 0;
  ANA_CHECK(m_event->retrieve(m_JetETMissNeutralParticleFlowObjects,"JetETMissNeutralParticleFlowObjects"));
  PrintPFOInfo( m_JetETMissChargedParticleFlowObjects,m_JetETMissNeutralParticleFlowObjects, m_PrintDebug);

  //---------------------------
  // EMTopoCluster and PFO cluster
  //---------------------------
  m_topocluster = 0;
  ANA_CHECK(m_event->retrieve( m_topocluster, "CaloCalTopoClusters"));
  PrintClusterInfo(m_topocluster, m_PrintDebug);
  m_PFOcluster = 0;
  if(m_SinglePionLowPerformanceStudies || m_DijetLowPerformance || m_DijetSubtraction) {
    ANA_CHECK(m_event->retrieve( m_PFOcluster, "PFOClusters_JetETMiss"));
    PrintPFOClusterInfo(m_PFOcluster, m_PrintDebug);
  }
  

  //----------------------------
  // Jet information
  //--------------------------- 
  m_Jets = 0;
  m_PFlowJets = 0;
  ANA_CHECK(m_event->retrieve( m_Jets, "AntiKt4EMTopoJets" ));
  Info("execute()", "  number of jets = %lu", m_Jets->size());
  ANA_CHECK(m_event->retrieve( m_PFlowJets, "AntiKt4EMPFlowJets" ));
  Info("execute()", "  number of PFlow jets = %lu", m_PFlowJets->size());
  PrintJetCollectionInfo(m_Jets,m_PFlowJets, m_PrintDebug);

  //----------------------------
  // Electrons 
  //--------------------------- 
  m_Electrons = 0;
  ANA_CHECK(m_event->retrieve(m_Electrons, "Electrons") );
  Info("execute()", "  number of electrons = %lu", m_Electrons->size());
  PrintElectronInfo(m_Electrons, true);
  
  //***Do we need the forward electrons?

  
  //---------------------------
  // Muons 
  //--------------------------- 
  m_Muons = 0;
  ANA_CHECK(m_event->retrieve( m_Muons, "Muons" ));
  Info("execute()", "  number of muons = %lu", m_Muons->size());
  PrintMuonInfo(m_Muons, true);


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


  
  if(m_Zmumu && !trigger){ // I added the trigger here. It seems to be the wrong chain
    int numGoodJets = 0;
    
    
    //first loop over the cleaning because it is not using PFlow
    
    xAOD::JetContainer::const_iterator jetEM_itr =  m_Jets->begin();
    xAOD::JetContainer::const_iterator jetEM_end =  m_Jets->end();
    for( ; jetEM_itr != jetEM_end; ++jetEM_itr ) {
		
		if( !m_jetCleaning->accept( **jetEM_itr )){
		   Info("execute()", "We failed");
		   return EL::StatusCode::SUCCESS; //only keep good clean EVENTS due to missing JetCleaning for PFLOW
		   
	   }
      numGoodJets++;
		
		
	}

    // Create the new container and its auxiliary store to store the good and calibrated jets .
    // AuxContainerBase is used because if we were using a derivation as input some of the original auxillary variables
    // may have been slimmed away (removed to make the container smaller), so if we were to do a deep-copy of the full JetAuxContainer then we would make our container larger than necessary.
    xAOD::JetContainer* goodPFlowJets = new xAOD::JetContainer();
    xAOD::AuxContainerBase* goodPFlowJetsAux = new xAOD::AuxContainerBase();
    goodPFlowJets->setStore( goodPFlowJetsAux );
    
    xAOD::JetContainer::const_iterator jet_itr =  m_PFlowJets->begin();
    xAOD::JetContainer::const_iterator jet_end =  m_PFlowJets->end();
    for( ; jet_itr != jet_end; ++jet_itr ) {
      
      Info("Execute () ", "Jet before calibration E = %.2f GeV  pt = %.2f GeV eta = %.2f  phi =  %.2f",
	   (*jet_itr)->e()/GEV,(*jet_itr)->pt()/GEV, (*jet_itr)->eta(), (*jet_itr)->phi());
      
      //Cleaning TOOL
      
      
      //Calibration Tool (Jet Calibration)
      xAOD::Jet* jet = 0;//new xAOD::Jet();
      m_akt4EMPFlowCalibrationTool->calibratedCopy(**jet_itr,jet); //make a calibrated copy, assuming a copy hasn't been made already
      
      Info("Execute () ", "Jet after Calibration E = %.10f GeV  pt = %.10f GeV eta = %.2f  phi =  %.2f",
	   jet->e()/GEV, jet->pt()/GEV, jet->eta(), jet->phi());
      goodPFlowJets->push_back(jet); // jet acquires the m_akt4CalibEMTopo auxstore
      //JER Tool (Jet Energy Resolution)
      double resMC = m_JERTool->getRelResolutionMC(jet);
      double resData = m_JERTool->getRelResolutionData(jet);
      
      Info("Execute () ","resMC = %.10f  resData = %.10f ", resMC, resData);
      
      ANA_CHECK(m_SmearTool->applyCorrection(*jet));
      //virtual CP::CorrectionCode applyCorrection(xAOD::Jet& jet);
      
      if(!(m_jetsf->passesJvtCut(*jet))) continue; //*** Does it have to be applied to all jets or only those with pt < 40GeV ?
      
      Info("Execute () ", "Jet after Smearing E = %.10f GeV  pt = %.10f GeV eta = %.2f  phi =  %.2f",
	   jet->e()/GEV, jet->pt()/GEV, jet->eta(), jet->phi());
	   
	   
     //*jet= **jet_itr; // copies auxdata from one auxstore to the other
      
      
    }
    
    //Just to check that GoodEMTopoJets has been store properly
    jet_itr =  goodPFlowJets->begin();
    jet_end =  goodPFlowJets->end();
    for( ; jet_itr != jet_end; ++jet_itr ) {
      Info("Execute () ", "goodEMTopoJets: E = %.2f GeV  pt = %.2f GeV eta = %.2f  phi =  %.2f",
	   (*jet_itr)->e()/GEV,(*jet_itr)->pt()/GEV, (*jet_itr)->eta(), (*jet_itr)->phi());
    }
    
    Info("execute()", "Number of good Jets: %i", numGoodJets);
    
    
    //** It has to be revisited. CorrectCopy doing nothing, not sure if the Medium is OK implemented, d0&z0 cuts still missing
    //---------------------------
    // Tools for Electrons
    //--------------------------- 
    // Create the new container and its auxiliary store to store the good and calibrated jets.
    xAOD::ElectronContainer* goodElectrons = new xAOD::ElectronContainer();
    xAOD::AuxContainerBase* goodElectronsAux = new xAOD::AuxContainerBase();
    goodElectrons->setStore( goodElectronsAux );
    
    // //Loop over vertices and look for good primary vertex (needed for z0 electron cut)
    // for (xAOD::VertexContainer::const_iterator vxIter = vxContainer->begin(); vxIter != vxContainer->end(); ++vxIter) {
    //   // Select good primary vertex
    //   if ((*vxIter)->vertexType() == xAOD::VxType::PriVtx) {
    //     //This is the primary vertex
    //   }
    // }
    
    xAOD::ElectronContainer::const_iterator el_itr = m_Electrons->begin();
    xAOD::ElectronContainer::const_iterator el_end = m_Electrons->end();
    for( ; el_itr != el_end; ++el_itr ) {
      
      Info("Execute () ", "Electron before slecetion/calibration E = %.2f GeV  pt = %.2f GeV eta = %.2f  phi =  %.2f",
	   (*el_itr)->e()/GEV,(*el_itr)->pt()/GEV, (*el_itr)->eta(), (*el_itr)->phi());
      
      //Reject bad electrons
      if( !(*el_itr)->isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON) ) continue;
      std::cout<<"Pass isGoodOQ"<<std::endl;
      // Formally unnecessary because all electrons in the container have these authors by construction
      if ( !(*el_itr)->author(xAOD::EgammaParameters::AuthorElectron) && !(*el_itr)->author(xAOD::EgammaParameters::AuthorAmbiguous) ) continue;
      std::cout<<"Pass Author"<<std::endl;
      // Calorimeter crack excluded
      if  ( fabs( (*el_itr)->caloCluster()->etaBE(2) ) >1.37 &&  fabs( (*el_itr)->caloCluster()->etaBE(2) ) <1.52) continue;
      std::cout<<"Is not in the crack region"<<std::endl;
      // Reject electrons outside the kinematic acceptance
      if (( (*el_itr)->pt() < 7000) || (fabs( (*el_itr)->caloCluster()->etaBE(2) ) > 2.47 )) continue;
      
      std::cout<<"Kinematic region selected"<<std::endl;
      // // Electron d0 and z0 cut: |d0BL significance |<5 and |Δz0BL*sinθ|<0.5 mm (***still in progress)
      // xAOD::TrackParticle *tp_el = (*el_itr)->trackParticle() ; 
      // double d0sig = xAOD::TrackingHelpers::d0significance( tp_el, m_event->beamPosSigmaX(), m_event->beamPosSigmaY(), m_event->beamPosSigmaXY() );
      // xAOD::Vertex vtx; //(***) Primary vertex has to be computed
      // //PV has to be chosen properly
      // float delta_z0 = fabs(tp_el->z0() + tp_el->vz() - vtx->z());
      // if(fabs(d0sig) < 5) >= continue;
      // if(fabs(delta_z0*TMath::sin(tp_el->theta)) >= 0.5) continue;
      // Electron Identification tool
      // std::cout<<"MediumLH = "<<m_MediumLH->accept(*el_itr)<<std::endl;
      std::cout<<"VeryLooseLHElectron = "<<m_VeryLooseLHElectron->accept(*el_itr)<<std::endl;
      if(!m_VeryLooseLHElectron->accept(*el_itr)) continue; 
      
      //    if(!m_MediumLH->accept(*el_itr)) continue; 
      //std::cout<<"is LHMediumSel"<<std::endl;
      //Isolation
      //if (!m_iso->accept( **el_itr )) continue; //***Has to be applied after or before calibration?  
      //std::cout<<"is Isolated"<<std::endl;
      
      //Calibration Tool (Jet Calibration)
      xAOD::Electron* el = NULL;
      m_electronCalibrationAndSmearingTool->correctedCopy(**el_itr,el); //make a calibrated copy, assuming a copy hasn't been made already
      goodElectrons->push_back(el); // jet acquires the m_akt4CalibEMTopo auxstore
      *el= **el_itr; // copies auxdata from one auxstore to the other
      
      //Info("execute()", "corrected electron pt = %.2f GeV", ((*el_itr)->pt() * 0.001));  
      //if (!m_iso->accept( **el_itr )) continue; //***Has to be applied after or before calibration? 
      // m_electronCalibrationAndSmearingTool->applyCorrection(**el_itr);
      // xAOD::Electron* electrondc = new xAOD::Electron();
      //goodElectrons->push_back( electrondc ); // jet acquires the goodJets auxstore
      // *electrondc= **el_itr; // copies auxdata from one auxstore to the other
    }
    
    Info("execute()", "  number of good electrons = %lu", goodElectrons->size());
    //Just to check that GoodElectrons has been store properly
    el_itr =  goodElectrons->begin();
    el_end =  goodElectrons->end();
    for( ; el_itr != el_end; ++el_itr ) {
      Info("Execute () ", "goodElectrons: E = %.2f GeV  pt = %.2f GeV eta = %.2f  phi =  %.2f",
	   (*el_itr)->e()/GEV,(*el_itr)->pt()/GEV, (*el_itr)->eta(), (*el_itr)->phi());
    }
    
    
    //*** Problem with corrected copy ; not stored properly! Have a look into the push_back
    //---------------------------
    // Tools for Muons
    //---------------------------
    // Create the new container and its auxiliary store to store the good and calibrated jets .
    xAOD::MuonContainer* goodMuons = new xAOD::MuonContainer();
    xAOD::AuxContainerBase* goodMuonsAux = new xAOD::AuxContainerBase();
    goodMuons->setStore( goodMuonsAux ); //< Connect the two
    
    xAOD::MuonContainer::const_iterator mu_itr = m_Muons->begin();
    xAOD::MuonContainer::const_iterator mu_end = m_Muons->end();
    for( ; mu_itr != mu_end; ++mu_itr ) {
      
      Info("Execute () ", "Muons before calibration: E = %.2f GeV  pt = %.2f GeV eta = %.2f  phi =  %.2f",
	   (*mu_itr)->e()/GEV,(*mu_itr)->pt()/GEV, (*mu_itr)->eta(), (*mu_itr)->phi());
      
      xAOD::Muon* mu = 0;
      m_muonCalibrationAndSmearingTool->correctedCopy(**mu_itr, mu);
      
      if (!m_iso->accept( **mu_itr )) continue;
      if ((((*mu_itr)->pt()/GEV)<25) || (fabs((*mu_itr)->eta())>2.4))continue;
      
      Info("execute()", "corrected muon pt = %.2f GeV", ((*mu_itr)->pt()/GEV));
      Info("execute()", "corrected muon pt (from copy) = %.2f GeV", (mu->pt()/GEV)); 
      
      goodMuons->push_back( mu ); // jet acquires the goodJets auxstore
      *mu= **mu_itr; // copies auxdata from one auxstore to the other
      
    } // end for loop over shallow copied muons
    
    
    Info("execute()", "  number of good and calibrated muons = %lu", goodMuons->size());
    //Just to check that GoodMuons has been store properly
    mu_itr =  goodMuons->begin();
    mu_end =  goodMuons->end();
    for( ; mu_itr != mu_end; ++mu_itr ) {
      Info("Execute () ", "GoodMuons: E = %.2f GeV  pt = %.2f GeV eta = %.2f  phi =  %.2f",
	   (*mu_itr)->e()/GEV,(*mu_itr)->pt()/GEV, (*mu_itr)->eta(), (*mu_itr)->phi());
    }
    
    //---------------------------
    // Zmumu selection
    //---------------------------
    if (ZmumuSelection(goodElectrons, goodMuons)){
      FillZmumuHistograms(goodMuons);
      JetRecoil_Zmumu(goodMuons, goodPFlowJets);
    }
  }
  
  
  
  //---------------------
  // Performance studies
  //----------------------
  
  if(m_SinglePionLowPerformanceStudies || m_DijetLowPerformance || m_DijetSubtraction){
    resize_tpVectors(m_TruthParticles);
    resize_PFOVectors(m_JetETMissChargedParticleFlowObjects);
    initialise_PFOVectors((int)m_TruthParticles->size(), (int)m_topocluster->size(), (int)m_JetETMissChargedParticleFlowObjects->size());
    fill_PFOVectors(m_JetETMissChargedParticleFlowObjects);
    //truth particle selection
    tp_Selection(m_TruthParticles,m_JetETMissChargedParticleFlowObjects);
    //associate calibration hits to truth particles
    ComputeCalibHitsPerParticle(m_CalCellInfo_TopoCluster,m_CalCellInfo,m_TruthParticles);
    //associate calibration hits per cluster
    ComputeCalibHitsPerCluster(m_CalCellInfo_TopoCluster,m_topocluster, int(m_TruthParticles->size()));
    //calculate efficieny and purity
    Calculate_Efficiency_Purity(m_TruthParticles, (int)m_topocluster->size(), m_topocluster);    
    //subtraction code
    if(m_DijetSubtraction){
      std::cout<< "In subtraction"<<std::endl;
      SumClusterE_ConeR(m_JetETMissChargedParticleFlowObjects,m_topocluster,0.10);
      SumClusterE_ConeR(m_JetETMissChargedParticleFlowObjects,m_topocluster,0.15);
      SumClusterE_ConeR(m_JetETMissChargedParticleFlowObjects,m_topocluster,0.20);
      SubtractionPerf(m_TruthParticles);
    }
    //fill histograms
    fill_RPlus_R0(m_TruthParticles);
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
  if (m_grl) {
    delete m_grl;
    m_grl = 0;
  }
  
  if( m_jetCleaning ) {
    delete m_jetCleaning;
    m_jetCleaning = 0;
  }
  
  if( m_akt4EMTopoCalibrationTool  ) {
    delete m_akt4EMTopoCalibrationTool;
    m_akt4EMTopoCalibrationTool = 0;
  }
  
  if(m_akt4EMPFlowCalibrationTool){
    delete m_akt4EMPFlowCalibrationTool;
    m_akt4EMPFlowCalibrationTool = 0;
  }
  
  if( m_JERTool && m_SmearTool ) {
    delete m_JERTool;
    m_JERTool = 0;
    delete m_SmearTool;
    m_SmearTool = 0;
  }
  
  if(m_muonCalibrationAndSmearingTool){
    delete m_muonCalibrationAndSmearingTool;
    m_muonCalibrationAndSmearingTool = 0;
  }
  
  if(m_electronCalibrationAndSmearingTool){
    delete m_electronCalibrationAndSmearingTool;
    m_electronCalibrationAndSmearingTool = 0;
  }
  
  if( m_iso){
      delete m_iso;
      m_iso=0;
  }
  
  if( m_MediumLH){
    delete m_MediumLH;
    m_MediumLH = 0;
  }
  
  if( m_trigConfigTool ) {
    delete m_trigConfigTool;
    m_trigConfigTool = 0;
  }
  
  if( m_trigDecisionTool ) {
    delete m_trigDecisionTool;
    m_trigDecisionTool = 0;
  }

  if(m_jetsf) {
    delete m_jetsf;
    m_jetsf = 0;
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



