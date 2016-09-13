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


xAODPFlowAna :: xAODPFlowAna (bool SinglePionLowPerformanceStudies, bool DijetLowPerformance, bool DijetSubtraction, bool Zmumu, bool matching)
{
  m_SinglePionLowPerformanceStudies = SinglePionLowPerformanceStudies;
  m_DijetLowPerformance = DijetLowPerformance;
  m_DijetSubtraction = DijetSubtraction;
  m_Zmumu = Zmumu;
  m_1to2matching = matching;

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


///////////////////////////
// Histograms booking
//////////////////////////
void xAODPFlowAna :: bookH1DHistogram(std::string name, int n_bins, float x_low, float x_up){
//  Info("bookH1DHistogram () ", "bookH1DHistograms...");
  
  TH1D* h1 = new TH1D(name.c_str(),name.c_str(),n_bins, x_low, x_up);
  h1->Sumw2();
  m_H1Dict[name] = h1;
  wk()->addOutput (m_H1Dict[name]);
  return; 
}

std::string xAODPFlowAna::histName(unsigned i_pt, unsigned i_eta, const std::string& name, const std::string& matchScheme, std::vector<float>& PtRange,
                                   std::vector<float>& EtaRange) {

  std::string complete_name = "WrongName";
  if (i_pt != PtRange.size()-1 && i_eta != EtaRange.size()-1) {
    complete_name = (name + matchScheme + "_" + std::to_string((int) (PtRange.at(i_pt))) + "_" + std::to_string((int) (PtRange.at(i_pt+1))) + "GeV__eta"
                     + std::to_string((int) ((10 * EtaRange.at(i_eta)))) + "_" + std::to_string((int) ((10 * EtaRange.at(i_eta+1))))).c_str();
  }
  else {
    complete_name = (name + matchScheme + "_" + std::to_string((int) (PtRange.at(i_pt))) + "GeV__eta"
                      + std::to_string((int) ((10 * EtaRange.at(i_eta)))) ).c_str();
  }

  return complete_name;
}

void xAODPFlowAna :: bookH1DPerformanceHistogram(std::string name, std::string matchScheme, std::vector<float> PtRange, std::vector<float> EtaRange, int n_bins, float x_low, float x_up)
{
//  Info("bookH1DPerformanceHistogram () ", "bookH1DPerformanceHistograms...");

   for (unsigned i_pt = 0; i_pt<PtRange.size(); i_pt++){
    for (unsigned i_eta =0; i_eta<EtaRange.size(); i_eta++){
      
      std::string complete_name = histName(i_pt, i_eta, name, matchScheme, PtRange, EtaRange);

      TH1D* h1 = new TH1D(complete_name.c_str(), complete_name.c_str(),n_bins, x_low, x_up);
      h1->Sumw2();
      m_H1Dict[complete_name] = h1;
      wk()->addOutput (m_H1Dict[complete_name]);
    }
  }
  return ;
} 


EL::StatusCode xAODPFlowAna :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.This method gets called before any input files are
  // connected.


  //Options form the histograms - to the config file
  std::string _matchScheme = (std::string) "_EM2";
  bool m_UseNarrowPtRange = true;
  bool m_UseNarrowEtaRange = true;
  
  if (m_UseNarrowPtRange) _ptRange= {0, 2, 5, 10, 20};
  else _ptRange = {0, 2, 5, 10, 20, 40, 60, 80, 100, 150, 200, 500, 1000};
  
  if (m_UseNarrowEtaRange) _etaRange= {0, 1, 2, 2.5};
  else _etaRange= {0.0, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 2.0, 2.5};

  
  // _ptRange= {"_0_2GeV","_2_5GeV","_5_10GeV","_10_20GeV","_20_40GeV","_40_60GeV","_60_2GeV","_2_5GeV","_5_10GeV","_0_2GeV","_2_5GeV","_5_10GeV" };
  //  std::vector<std::string> _etaRange= {"","_eta06","_eta08","_eta1","_eta12","_eta14","_eta16","_eta2","_eta25"};
  
  //====================================
  // DiJet subtraction
  //====================================
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

  /* WIP: a directory structure has to be created to store the histograms */
  
  //====================================
  // Track cluster matching histograms
  //====================================
  int n_bins = 100; float x_low = 0.; float x_up = 10.;
  bookH1DPerformanceHistogram("dR",_matchScheme, _ptRange, _etaRange, n_bins, x_low, x_up);
  
  //====================================
  // Efficiency and purity
  //====================================
  int n_effbins = 20; float eff_low = 0; float eff_up = 1.0001;
  if (m_1to2matching) {
    bookH1DPerformanceHistogram("EffMatch1","",_ptRange, _etaRange, n_effbins, eff_low, eff_up);
    bookH1DPerformanceHistogram("PurMatch1","",_ptRange, _etaRange, n_effbins, eff_low, eff_up);
    bookH1DPerformanceHistogram("EffMatchboth","",_ptRange, _etaRange, n_effbins, eff_low, eff_up);
    bookH1DPerformanceHistogram("PurMatch2","",_ptRange, _etaRange, n_effbins, eff_low, eff_up);
  } else {
    bookH1DPerformanceHistogram("Eff", "", _ptRange, _etaRange, n_effbins, eff_low, eff_up);
    bookH1DPerformanceHistogram("Pur","",_ptRange, _etaRange, n_effbins, eff_low, eff_up);
  }

  //====================================
  // 1->2 Matching
  //====================================
  bookH1DPerformanceHistogram("SubtractStatus","",_ptRange, _etaRange, 5, 0, 5);
  bookH1DPerformanceHistogram("EOP","",_ptRange, _etaRange, 5, 0, 1.5);
  bookH1DPerformanceHistogram("EOPtotal","",_ptRange, _etaRange, 5, 0, 1.5);
  bookH1DPerformanceHistogram("Energy1st","",_ptRange, _etaRange, 20, 0, 20);
  bookH1DPerformanceHistogram("Energy2rd","",_ptRange, _etaRange, 20, 0, 20);

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

  //FInal histograms
  //Muon histograms
  bookH1DHistogram("h_muonPt", pt_bin, pt_low, pt_up);
  bookH1DHistogram("h_muonE", E_bin, E_low, E_up);
  bookH1DHistogram("h_muonM", 20, 0, 100); //Is this a proper range? 
  bookH1DHistogram("h_muonEta", eta_bin, eta_low, eta_up);
  bookH1DHistogram("h_muonPhi", phi_bin, phi_low, phi_up);

  //Z distributions
  bookH1DHistogram("h_ZPt", pt_bin, pt_low, pt_up);
  bookH1DHistogram("h_ZE", E_bin, E_low, E_up);
  bookH1DHistogram("h_ZM", 20, 30, 180); //Is this a proper range? 
  bookH1DHistogram("h_ZEta", eta_bin, eta_low, eta_up);
  bookH1DHistogram("h_ZPhi", phi_bin, phi_low, phi_up);
  //Z+jet system
  bookH1DHistogram("h_ZPt_to_JetPt", 20, 0, 5);
  bookH1DHistogram("h_ZPt_to_JetPt_sum", 20, 0, 5);
    
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
//  Info("initialize()", "Number of events = %lli", m_event->getEntries() );

  m_store = new xAOD::TStore();

  // Variable initialization
  m_eventCounter = 0; //Count number of events
  // Printing
  PrintDebug = false; //Printing message criteria -->  Should be chosen from the ATestRun  
  
  // Conversion factors
  GEV = 1000.; //Units
  
  //Counters for Data-MC CutFlow (under construction)
  m_select = 0;
  m_trigger = 0;
  m_number = 0;
  m_jvt = 0;
  m_angle = 0;
  m_jetpt = 0;
  m_nojetsafterfilter = 0;

  
  /* WIP: Should isData be initialized here foe the tool? */ 

  
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
  //***Two differences with respect to Christian configuration 1-JetArea_Residual_Origin_EtaJES_GSC and 2-isData = false
  
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
  const xAOD::EventInfo* eventInfo = 0;
  ANA_CHECK(m_event->retrieve( eventInfo, "EventInfo"));  
  
  // check if the event is data or MC
  bool isMC = false;
  m_EvtWeight = 1.0;
  
  if( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) )
    isMC = true; 
   
  if( isMC ) { m_EvtWeight = eventInfo->mcEventWeight();}
//  Info("execute()", "Event number = %llu  Run Number =  %d  Event weight = %.2f  isMC = %s",eventInfo->eventNumber(), eventInfo->runNumber(), m_EvtWeight, (isMC ? "true" : "false"));


  //trigger tools: here the trigger chain is chosen
  auto chainGroup = m_trigDecisionTool->getChainGroup("HLT_mu20_iloose_L1MU15, HLT_mu50");
  std::map<std::string,int> triggerCounts;
  int trigger = 0;
  for(auto &trig : chainGroup->getListOfTriggers()) {
    auto cg = m_trigDecisionTool->getChainGroup(trig);
    std::string thisTrig = trig;
    //Info( "execute()", "%30s chain passed(1)/failed(0): %d total chain prescale (L1*HLT): %.1f", thisTrig.c_str(), cg->isPassed(), cg->getPrescale() );
    if(cg->isPassed())trigger = 1;
    //Info("execute()", "Trigger + %i", trigger);
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
//    PrintTruthInfo(m_TruthParticles, m_TruthVertices, PrintDebug);
    
    //---------------------------
    // CalCellInfo_TopoCluster
    //---------------------------
    m_CalCellInfo_TopoCluster = 0;
    ANA_CHECK(m_event->retrieve(m_CalCellInfo_TopoCluster, "CalCellInfo_TopoCluster"));
    m_CalCellInfo = 0; //CalCellInfo PFO
    ANA_CHECK(m_event->retrieve(m_CalCellInfo, "CalCellInfo"));
//    PrintCalCellInfo(m_CalCellInfo_TopoCluster,m_CalCellInfo, PrintDebug);
  }

  //---------------------------
  // Track Collection
  //---------------------------
  m_InDetTrackParticles  = 0;
  ANA_CHECK(m_event->retrieve( m_InDetTrackParticles ,"InDetTrackParticles"));
//  PrintTrackInfo(m_InDetTrackParticles,PrintDebug);
  
  //---------------------------
  // cPFO and nPFO
  //---------------------------
  m_JetETMissChargedParticleFlowObjects  = 0;
  ANA_CHECK(m_event->retrieve( m_JetETMissChargedParticleFlowObjects ,"JetETMissChargedParticleFlowObjects"));
  m_JetETMissNeutralParticleFlowObjects = 0;
  ANA_CHECK(m_event->retrieve(m_JetETMissNeutralParticleFlowObjects,"JetETMissNeutralParticleFlowObjects"));
//  PrintPFOInfo( m_JetETMissChargedParticleFlowObjects,m_JetETMissNeutralParticleFlowObjects, PrintDebug);

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
//  Info("execute()", "  number of jets = %lu", m_Jets->size());
  ANA_CHECK(m_event->retrieve( m_PFlowJets, "AntiKt4EMPFlowJets" ));
//  Info("execute()", "  number of PFlow jets = %lu", m_PFlowJets->size());
//  PrintJetCollections(m_Jets,m_PFlowJets, true);

  //----------------------------
  // Electrons (**Print options)
  //--------------------------- 
  m_Electrons = 0;
  ANA_CHECK(m_event->retrieve(m_Electrons, "Electrons") );
//  Info("execute()", "  number of electrons = %lu", m_Electrons->size());
  
  //***Do we need the forward electrons?

  
  //---------------------------
  // Muons (**Print options)
  //--------------------------- 
  m_Muons = 0;
  ANA_CHECK(m_event->retrieve( m_Muons, "Muons" ));
//  Info("execute()", "  number of muons = %lu", m_Muons->size());
  

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
  // AuxContainerBase is used because if we were using a derivation as input some of the original auxillary variables
  // may have been slimmed away (removed to make the container smaller), so if we were to do a deep-copy of the full JetAuxContainer then we would make our container larger than necessary.
  xAOD::JetContainer* goodEMTopoJets = new xAOD::JetContainer();
  xAOD::AuxContainerBase* goodEMTopoJetsAux = new xAOD::AuxContainerBase();
  goodEMTopoJets->setStore( goodEMTopoJetsAux );
  
  xAOD::JetContainer::const_iterator jet_itr =  m_Jets->begin();
  xAOD::JetContainer::const_iterator jet_end =  m_Jets->end();
  for( ; jet_itr != jet_end; ++jet_itr ) {
    
//    Info("Execute () ", "Jet before calibration E = %.2f GeV  pt = %.2f GeV eta = %.2f  phi =  %.2f", (*jet_itr)->e()/GEV,(*jet_itr)->pt()/GEV, (*jet_itr)->eta(), (*jet_itr)->phi());
    
    //Cleaning TOOL
    //Should we remove the whole event or only the jet? Top analyses remove the whole event.
    if( !m_jetCleaning->accept( **jet_itr )) continue; //only keep good clean jets
    numGoodJets++;
    
    //Calibration Tool (Jet Calibration)
    xAOD::Jet* jet = new xAOD::Jet();
    m_akt4EMTopoCalibrationTool->calibratedCopy(**jet_itr,jet); //make a calibrated copy, assuming a copy hasn't been made already
    
//    Info("Execute () ", "Jet after Calibration E = %.10f GeV  pt = %.10f GeV eta = %.2f  phi =  %.2f", jet->e()/GEV, jet->pt()/GEV, jet->eta(), jet->phi());
    
    //JER Tool (Jet Energy Resolution)
    double resMC = m_JERTool->getRelResolutionMC(jet);
    double resData = m_JERTool->getRelResolutionData(jet);
    
//    Info("Execute () ","resMC = %.10f  resData = %.10f ", resMC, resData);
    
    ANA_CHECK(m_SmearTool->applyCorrection(*jet));
    //virtual CP::CorrectionCode applyCorrection(xAOD::Jet& jet);

    if(!(m_jetsf->passesJvtCut(*jet))) continue; //*** Does it have to be applied to all jets or only those with pt < 40GeV ?
    
    goodEMTopoJets->push_back(jet); // jet acquires the m_akt4CalibEMTopo auxstore
    *jet= **jet_itr; // copies auxdata from one auxstore to the other
    
//    Info("Execute () ", "Jet after Smearing E = %.10f GeV  pt = %.10f GeV eta = %.2f  phi =  %.2f", jet->e()/GEV, jet->pt()/GEV, jet->eta(), jet->phi());
  }

  //Just to check that GoodEMTopoJets has been store properly
  jet_itr =  goodEMTopoJets->begin();
  jet_end =  goodEMTopoJets->end();
  for( ; jet_itr != jet_end; ++jet_itr ) {
//    Info("Execute () ", "goodEMTopoJets: E = %.2f GeV  pt = %.2f GeV eta = %.2f  phi =  %.2f", (*jet_itr)->e()/GEV,(*jet_itr)->pt()/GEV, (*jet_itr)->eta(), (*jet_itr)->phi());
  }
  

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
    
//    Info("Execute () ", "Electron before slecetion/calibration E = %.2f GeV  pt = %.2f GeV eta = %.2f  phi =  %.2f", (*el_itr)->e()/GEV,(*el_itr)->pt()/GEV, (*el_itr)->eta(), (*el_itr)->phi());
    
    //Reject bad electrons
    if( !(*el_itr)->isGoodOQ(xAOD::EgammaParameters::BADCLUSELECTRON) ) continue;
//    std::cout<<"Pass isGoodOQ"<<std::endl;
    // Formally unnecessary because all electrons in the container have these authors by construction
    if ( !(*el_itr)->author(xAOD::EgammaParameters::AuthorElectron) && !(*el_itr)->author(xAOD::EgammaParameters::AuthorAmbiguous) ) continue;
//    std::cout<<"Pass Author"<<std::endl;
    // Calorimeter crack excluded
    if  ( fabs( (*el_itr)->caloCluster()->etaBE(2) ) >1.37 &&  fabs( (*el_itr)->caloCluster()->etaBE(2) ) <1.52) continue;
//    std::cout<<"Is not in the crack region"<<std::endl;
    // Reject electrons outside the kinematic acceptance
    if (( (*el_itr)->pt() < 7000) || (fabs( (*el_itr)->caloCluster()->etaBE(2) ) > 2.47 )) continue;
    
 //   std::cout<<"Kinematic region selected"<<std::endl;
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
 //   std::cout<<"VeryLooseLHElectron = "<<m_VeryLooseLHElectron->accept(*el_itr)<<std::endl;
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
  
//  Info("execute()", "  number of good electrons = %lu", goodElectrons->size());
  //Just to check that GoodElectrons has been store properly
  el_itr =  goodElectrons->begin();
  el_end =  goodElectrons->end();
//  for( ; el_itr != el_end; ++el_itr ) {
//    Info("Execute () ", "goodElectrons: E = %.2f GeV  pt = %.2f GeV eta = %.2f  phi =  %.2f", (*el_itr)->e()/GEV,(*el_itr)->pt()/GEV, (*el_itr)->eta(), (*el_itr)->phi());
//  }
  
  
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
    
//    Info("Execute () ", "Muons before calibration: E = %.2f GeV  pt = %.2f GeV eta = %.2f  phi =  %.2f", (*mu_itr)->e()/GEV,(*mu_itr)->pt()/GEV, (*mu_itr)->eta(), (*mu_itr)->phi());
    
    xAOD::Muon* mu = 0;
    m_muonCalibrationAndSmearingTool->correctedCopy(**mu_itr, mu);
    
    if (!m_iso->accept( **mu_itr )) continue;
    if ((((*mu_itr)->pt()/GEV)<25) || (fabs((*mu_itr)->eta())>2.4))continue;

//    Info("execute()", "corrected muon pt = %.2f GeV", ((*mu_itr)->pt()/GEV));
//    Info("execute()", "corrected muon pt (from copy) = %.2f GeV", (mu->pt()/GEV)); 
     
    goodMuons->push_back( mu ); // jet acquires the goodJets auxstore
    *mu= **mu_itr; // copies auxdata from one auxstore to the other

  } // end for loop over shallow copied muons

  
//  Info("execute()", "  number of good and calibrated muons = %lu", goodMuons->size());
  //Just to check that GoodMuons has been store properly
  mu_itr =  goodMuons->begin();
  mu_end =  goodMuons->end();
//  for( ; mu_itr != mu_end; ++mu_itr ) {
//    Info("Execute () ", "GoodMuons: E = %.2f GeV  pt = %.2f GeV eta = %.2f  phi =  %.2f", (*mu_itr)->e()/GEV,(*mu_itr)->pt()/GEV, (*mu_itr)->eta(), (*mu_itr)->phi());
//  }
  

  
  
  //---------------------------
  // Zmumu selection
  //---------------------------
  //if (ZmumuSelection()){
    //Perform the analysis
  //}
  

  
  //---------------------
  // Performance studies
  //----------------------
  
  if(m_SinglePionLowPerformanceStudies || m_DijetLowPerformance || m_DijetSubtraction){
    resize_tpVectors(m_TruthParticles);
    resize_PFOVectors(m_JetETMissChargedParticleFlowObjects);
    initialise_PFOVectors((int)m_TruthParticles->size(), (int)m_topocluster->size(), (int)m_JetETMissChargedParticleFlowObjects->size());
    fill_PFOVectors(m_JetETMissChargedParticleFlowObjects);
    // Calculate the MinDeltaR between cPFO and truth particles and fill the vector _mc_MinDeltaREflowTrackPair
    CalculateMatrix_MinDeltaR (m_TruthParticles, m_JetETMissChargedParticleFlowObjects, 0.03);
    //truth particle selection
    tp_Selection(m_TruthParticles,m_JetETMissChargedParticleFlowObjects);
    //associate calibration hits to truth particles
    ComputeCalibHitsPerParticle(m_CalCellInfo_TopoCluster,m_CalCellInfo,m_TruthParticles);
    //associate calibration hits per cluster
    ComputeCalibHitsPerCluster(m_CalCellInfo_TopoCluster,m_topocluster, int(m_TruthParticles->size()));
    //calculate efficieny and purity
    Calculate_Efficiency_Purity(m_TruthParticles, (int)m_topocluster->size(), m_topocluster);    
    //subtraction code
    if(m_DijetSubtraction)SubtractionPerf(m_JetETMissChargedParticleFlowObjects,m_topocluster, m_TruthParticles);
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
  
//  Info("", "--- Bad Jets Scanning ---");
//  Info("BadJetsScan", "jet E = %.2f GeV  pt  = %.2f GeV eta = %.2f  phi =  %.2f", jet.e()/GEV, jet.pt()/GEV, jet.phi(), jet.eta());
  
  // if(HasPFlowJetMatched(jet)){
  //   Info("BadJetsScan", "PFlow jet matched  %d E = %.2f GeV  pt  = %.2f GeV eta = %.2f  phi =  %.2f",
  // 	 jet.e()/GEV, jet.pt()/GEV, jet.phi(), jet.eta());
  // }

   
  return;
}



