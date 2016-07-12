// ASG status code check
#include <AsgTools/MessageCheck.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"

#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <PFlowAna/xAODPFlowAnaEDM.h>
#include <PFlowAna/xAODPFlowAna.h>
#include <PFlowAna/PrintInfo.h>

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
  // trees.  This method gets called before any input files are
  // connected.  
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
  
  m_eventCounter = 0; //Count number of events
  PrintDebug = false; //Printing message criteria
  
  GEV = 1000.; //Units
  
  //Initialize and configure the jet cleaning tool
  m_jetCleaning = new JetCleaningTool("JetCleaning");
  m_jetCleaning->msg().setLevel( MSG::DEBUG ); 
  ANA_CHECK(m_jetCleaning->setProperty( "CutLevel", "LooseBad"));
  ANA_CHECK(m_jetCleaning->setProperty("DoUgly", false));
  ANA_CHECK(m_jetCleaning->initialize());

  
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
  double EvtWeight = 1.0;
  
  if( eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) )
    isMC = true; 
   
  if( isMC ) { EvtWeight = eventInfo->mcEventWeight();}
  Info("execute()", "Event number = %llu  Event weight = %.2f  isMC = %s",eventInfo->eventNumber(), EvtWeight, (isMC ? "true" : "false"));
  
  //---------------------------
  // Truth particles & vertices
  //---------------------------
  m_TruthParticles = 0;
  ANA_CHECK(m_event->retrieve( m_TruthParticles,"TruthParticles"));
  m_TruthVertices = 0;
  ANA_CHECK(m_event->retrieve(m_TruthVertices,"TruthVertices"));
  PrintTruthInfo(m_TruthParticles, m_TruthVertices, PrintDebug);

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
  ANA_CHECK(m_event->retrieve( m_PFOcluster, "PFOClusters_JetETMiss"));
  PrintClusterInfo(m_topocluster,m_PFOcluster, PrintDebug);
  
  //---------------------------
  // CalCellInfo_TopoCluster
  //---------------------------
  m_CalCellInfo_TopoCluster = 0;
  ANA_CHECK(m_event->retrieve(m_CalCellInfo_TopoCluster, "CalCellInfo_TopoCluster"));
  m_CalCellInfo = 0; //CalCellInfo PFO
  ANA_CHECK(m_event->retrieve(m_CalCellInfo, "CalCellInfo"));
  PrintCalCellInfo(m_CalCellInfo_TopoCluster,m_CalCellInfo, PrintDebug);

  //----------------------------
  // Jet information
  //--------------------------- 
  m_Jets = 0;
  m_PFlowJets = 0;
  ANA_CHECK(m_event->retrieve( m_Jets, "AntiKt4EMTopoJets" ));
  Info("execute()", "  number of jets = %lu", m_Jets->size());
  ANA_CHECK(m_event->retrieve( m_PFlowJets, "AntiKt4EMPFlowJets" ));
  Info("execute()", "  number of PFlow jets = %lu", m_PFlowJets->size());
  PrintJetCollections(m_Jets,m_PFlowJets, PrintDebug);
  
   //Jet matching
  MatchJetCollections(m_Jets, m_PFlowJets);

 
  /*
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
  
  ANA_CHECK_SET_TYPE (EL::StatusCode);



  if( m_jetCleaning ) {
    delete m_jetCleaning;
    m_jetCleaning = 0;
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



