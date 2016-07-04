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
  // Truth particles
  //---------------------------
  m_TruthParticles = 0;
  ANA_CHECK(m_event->retrieve( m_TruthParticles,"TruthParticles"));

  //---------------------------
  // Truth vertices
  //---------------------------
  m_TruthVertices = 0;
  ANA_CHECK(m_event->retrieve(m_TruthVertices,"TruthVertices"));
  PrintTruthInfo();

  //---------------------------
  // Track Collection
  //---------------------------
  m_InDetTrackParticles  = 0;
  ANA_CHECK(m_event->retrieve( m_InDetTrackParticles ,"InDetTrackParticles"));
  PrintTrackInfo();
  
  //---------------------------
  // cPFO and nPFO
  //---------------------------
  m_JetETMissChargedParticleFlowObjects  = 0;
  ANA_CHECK(m_event->retrieve( m_JetETMissChargedParticleFlowObjects ,"JetETMissChargedParticleFlowObjects"));
  m_JetETMissNeutralParticleFlowObjects = 0;
  ANA_CHECK(m_event->retrieve(m_JetETMissNeutralParticleFlowObjects,"JetETMissNeutralParticleFlowObjects"));
  PrintPFOInfo();

  //---------------------------
  // EMTopoCluster and PFO cluster
  //---------------------------
  m_topocluster = 0;
  ANA_CHECK(m_event->retrieve( m_topocluster, "CaloCalTopoClusters"));
  m_PFOcluster = 0;
  ANA_CHECK(m_event->retrieve( m_PFOcluster, "PFOClusters_JetETMiss"));
  PrintClusterInfo();

  //---------------------------
  // CalCellInfo_TopoCluster
  //---------------------------
  m_CalCellInfo_TopoCluster = 0;
  ANA_CHECK(m_event->retrieve(m_CalCellInfo_TopoCluster, "CalCellInfo_TopoCluster"));
  m_CalCellInfo = 0; //CalCellInfo PFO
  ANA_CHECK(m_event->retrieve(m_CalCellInfo, "CalCellInfo"));
  PrintCalCellInfo();

  //----------------------------
  // Jet information
  //--------------------------- 
  m_Jets = 0;
  m_PFlowJets = 0;

  //Cleaning
  int numGoodJets = 0;
  int numBadJets = 0;
  
  ANA_CHECK(m_event->retrieve( m_Jets, "AntiKt4EMTopoJets" ));
  Info("execute()", "  number of jets = %lu", m_Jets->size());
  ANA_CHECK(m_event->retrieve( m_PFlowJets, "AntiKt4EMPFlowJets" ));
  Info("execute()", "  number of PFlow jets = %lu", m_PFlowJets->size());
  PrintJetCollections();
  
  /*
  //For cleaning Study
  MatchJetCollections(m_Jets, m_PFlowJets);
  
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


 void xAODPFlowAna :: PrintTruthInfo (){
   
   Info("", "------------------- ");
   Info("", "   TruthParticles   ");
   Info("", "------------------- ");
      
   Info("PrintTruthInfo", "Number of truth particles = %lu", m_TruthParticles->size());
   if(PrintDebug){
     xAOD::TruthParticleContainer::const_iterator tp_itr = m_TruthParticles->begin();
     xAOD::TruthParticleContainer::const_iterator tp_end = m_TruthParticles->end();
     for( ; tp_itr != tp_end; ++tp_itr ) {
       int tp_index = std::distance(m_TruthParticles->begin(),tp_itr);
       Info("PrintTruthParticlesInfo () ", "Truth Particle index = %d barcode = %i E = %.2f GeV  pt = %.2f GeV eta = %.2f  phi =  %.2f",
	    tp_index,
	    (*tp_itr)->barcode(),
	    (*tp_itr)->e()/GEV,
	    (*tp_itr)->pt()/GEV,
	    (*tp_itr)->eta(),
	    (*tp_itr)->phi());
       
       // const xAOD::TruthVertex*  prodVtx =  (*tp_itr)->prodVtx();
       // if(prodVtx) Info("PrintTruthParticlesInfo", "Vtx (x,y,z) = (%.2f,%.2f,%.2f)",
       // 			prodVtx->x(),
       // 			prodVtx->y(),
       // 			prodVtx->z());
       // else Info("PrintTruthParticlesInfo", "No Vtx for mc particle");
     }
   }

   Info("", "------------------- ");
   Info("", "   TruthVertex      ");
   Info("", "------------------- ");
  
   Info("PrintTruthInfo", "Truth PV Vertex size = %lu ",m_TruthVertices->size());
   if(PrintDebug){
   Info("PrintTruthInfo", " Truth PV Vertex information (x,y,z) = (%.2lf,%.2lf,%.2lf)",
	m_TruthVertices->at(0)->x(),
	m_TruthVertices->at(0)->y(),
	m_TruthVertices->at(0)->z());
   }
  return;
}
 
 void xAODPFlowAna :: PrintTrackInfo (){

   Info("", "-------------------- ");
   Info("", " InDetTrackParticles ");
   Info("", "-------------------- ");
   
   Info("PrintTrackInfo", "Number of InDetTrackParticles = %lu",m_InDetTrackParticles->size());
   if(PrintDebug){
     xAOD::TrackParticleContainer::const_iterator idtrk_itr = m_InDetTrackParticles->begin();
     xAOD::TrackParticleContainer::const_iterator idtrk_end = m_InDetTrackParticles->end();
     for( ; idtrk_itr != idtrk_end; ++idtrk_itr ) {
       Info("PrintTrackInfo", "InDetTrackParticles charge = %f  E  = %.2f GeV  pt = %.2f GeV  eta = %.2f  phi = %.2f ",
	    (*idtrk_itr)->charge(),
	    (*idtrk_itr)->e()/GEV,
	    (*idtrk_itr)->pt()/GEV,
	    (*idtrk_itr)->eta(),
	    (*idtrk_itr)->phi());
     }
   }
   return;
 }


void xAODPFlowAna :: PrintPFOInfo (){

   Info("", "----------------- ");
   Info("", " Charged PFO      ");
   Info("", "----------------- ");
   
  Info("PrintPFOInfo", "Number of ChargedParticleFlowObjects = %lu",m_JetETMissChargedParticleFlowObjects->size());
  if(PrintDebug){
    xAOD::PFOContainer::const_iterator cpfo_itr = m_JetETMissChargedParticleFlowObjects->begin();
    xAOD::PFOContainer::const_iterator cpfo_end = m_JetETMissChargedParticleFlowObjects->end();
    for( ; cpfo_itr != cpfo_end; ++cpfo_itr ) {
      int cpfo_index = std::distance(m_JetETMissChargedParticleFlowObjects->begin(),cpfo_itr);
      Info("PrintPFOInfo", "Charged PFO %d E = %.2f GeV  pt = %.2f GeV  eta = %.2f  phi = %.2f",
	   cpfo_index, (*cpfo_itr)->e()/GEV,
	   (*cpfo_itr)->pt()/GEV,
	   (*cpfo_itr)->eta(),
	   (*cpfo_itr)->phi());     
      
      //Associated cluster
      const xAOD::CaloCluster* matchedCluster = (*cpfo_itr)->cluster(0);
      // raw = "UNCALIBRATED” – electromagnetic energy scale, cluster energy is cell energy sum including possible topological weights from cluster splitting
      if(matchedCluster)
	Info("PrintPFOInfo", "MatchedCluster_E  = %.3f, eta = %.3f, phi = %.3f ",matchedCluster->rawE(), matchedCluster->eta(), matchedCluster->phi());
      else Info("PrintPFOInfo", "No cluster matched to the cPFO");
    }
  }

  Info("", "----------------- ");
  Info("", " Neutral PFO      ");
  Info("", "----------------- ");
  
  Info("PrintPFOInfo", "Number of NeutralParticleFlowObjects = %lu", m_JetETMissNeutralParticleFlowObjects->size());
  if(PrintDebug){
    xAOD::PFOContainer::const_iterator npfo_itr = m_JetETMissNeutralParticleFlowObjects->begin();
    xAOD::PFOContainer::const_iterator npfo_end = m_JetETMissNeutralParticleFlowObjects->end();
    for( ; npfo_itr != npfo_end; ++npfo_itr ) {
      int npfo_index = std::distance(m_JetETMissNeutralParticleFlowObjects->begin(),npfo_itr);
      Info("PrintPFOInfo", "Neutral PFO %d E = %.2f GeV  pt  = %.2f GeV eta = %.2f  phi =  %.2f",
	   npfo_index,
	   (*npfo_itr)->e()/GEV,
	   (*npfo_itr)->pt()/GEV,
	   (*npfo_itr)->eta(),
	   (*npfo_itr)->phi());  
      
      // Associated clusters
      const xAOD::CaloCluster* matchedCluster = (*npfo_itr)->cluster(0);
      if (matchedCluster)
	Info("PrintPFOInfo", "MatchedCluster_E  = %.3f, eta = %.3f, phi = %.3f ",matchedCluster->rawE(), matchedCluster->eta(), matchedCluster->phi());
      else Info("PrintPFOInfo", "No cluster matched to the nPFO");
      //The energy at different layers can be gotten using clusterN->eSample(xAOD::CaloCluster::CaloSample::EMB1);
    }
  }
  return;
}

void xAODPFlowAna :: PrintClusterInfo (){
  
  Info("", "----------------- ");
  Info("", " TopoClusters     ");
  Info("", "----------------- ");

 Info("execute()", "Number of TopoClusters = %lu", m_topocluster->size());
 if(PrintDebug){
   xAOD::CaloClusterContainer::const_iterator CaloCluster_itr = m_topocluster->begin();
   xAOD::CaloClusterContainer::const_iterator CaloCluster_end = m_topocluster->end();
   for( ; CaloCluster_itr != CaloCluster_end; ++CaloCluster_itr  ) {
     int CaloCluster_index = std::distance(m_topocluster->begin(),CaloCluster_itr);
     Info("PrintTopoClusterInfo", "CaloClusterContainer %d E_em = %.2f GeV E_cal = %.2f GeV  pt  = %.2f GeV eta = %.2f  phi =  %.2f",
	  CaloCluster_index,
	  (*CaloCluster_itr)->rawE()/GEV,
	 (*CaloCluster_itr)->calE()/GEV,
	  (*CaloCluster_itr)->pt()/GEV,
	  (*CaloCluster_itr)->eta(),
	  (*CaloCluster_itr)->phi()); 
   }
 }
 Info("", "----------------- ");
  Info("", "  PFOCluster      ");
  Info("", "----------------- ");
  
  Info("PrintTopoClusterInfo", "Number of PFOClusters = %lu", m_PFOcluster->size());
  if(PrintDebug){
    xAOD::CaloClusterContainer::const_iterator pfo_cl_itr = m_PFOcluster->begin();
    xAOD::CaloClusterContainer::const_iterator pfo_cl_end = m_PFOcluster->end();
    for( ; pfo_cl_itr != pfo_cl_end; ++pfo_cl_itr ) {
      int pfo_cl_index = std::distance(m_PFOcluster->begin(),pfo_cl_itr);
      Info("PrintTopoClusterInfo", "PFO cluster %d E = %.2f GeV  pt  = %.2f GeV eta = %.2f  phi =  %.2f",
	   pfo_cl_index,
	   (*pfo_cl_itr)->e()/GEV,
	   (*pfo_cl_itr)->pt()/GEV,
	   (*pfo_cl_itr)->eta(),
	   (*pfo_cl_itr)->phi());
    }
  }
  return;
}

void xAODPFlowAna :: PrintCalCellInfo () {
  
  Info("", "--------------------------- ");
  Info("", "  CalCellInfoTopoCluster    ");
  Info("", "--------------------------- ");
   
  Info("PrintCalCellInfo", "Number of CalCellInfo_TopoCluster = %lu", m_CalCellInfo_TopoCluster->size());
  if(PrintDebug){
    xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_itr = m_CalCellInfo_TopoCluster->begin();
    xAOD::CalCellInfoContainer::const_iterator CalCellInfoTopoCl_end = m_CalCellInfo_TopoCluster->end();
    for( ; CalCellInfoTopoCl_itr != CalCellInfoTopoCl_end; ++CalCellInfoTopoCl_itr ) {
      int index = std::distance(m_CalCellInfo_TopoCluster->begin(),CalCellInfoTopoCl_itr);
      Info("PrintCalCellInfo", "CalCellInfo TopoCluster %d  barcode = %i  particleID =  %i cl_E  = %.2f GeV  cell_eta =  %.2f  cell_phi = %.2f  EMEnergy = %.2f GeV   nonEMEnergy = %.2f GeV",
	   index,
	   (*CalCellInfoTopoCl_itr)->barcode(),
	   (*CalCellInfoTopoCl_itr)->particleID(),
	   (*CalCellInfoTopoCl_itr)->clusterRecoEnergy()/GEV,
	   (*CalCellInfoTopoCl_itr)->cellEta(),
	   (*CalCellInfoTopoCl_itr)->cellPhi(),
	   (*CalCellInfoTopoCl_itr)->EMEnergy()/GEV,
	   (*CalCellInfoTopoCl_itr)->nonEMEnergy()/GEV);
    }
  }

  Info("PrintCalCellInfo", "Number of CalCellInfo = %lu", m_CalCellInfo->size());
  
  return;
}


void xAODPFlowAna :: PrintJetCollections () {
  
  Info("", "--- CalCellInfoTopoCluster ---");
  
   
  Info("PrintJetCollections", "Number of TopoJets = %lu", m_Jets->size());
  if(true){
    xAOD::JetContainer::const_iterator Jets_itr = m_Jets->begin();
    xAOD::JetContainer::const_iterator Jets_end = m_Jets->end();
    for( ; Jets_itr != Jets_end; ++Jets_itr ) {
      int index = std::distance(m_Jets->begin(),Jets_itr);
      Info("PrintJetCollections", "TopoJets E  = %.2f GeV  pt =  %.2f  eta = %.2f  phi = %.2f GeV",
	   (*Jets_itr)->e()/GEV,
	   (*Jets_itr)->pt()/GEV,
	   (*Jets_itr)->eta(),
	   (*Jets_itr)->phi());
    }
  }

  Info("PrintJetCollections", "Number of PFlowJets = %lu", m_PFlowJets->size());
  if(true){
    xAOD::JetContainer::const_iterator PFlowJets_itr = m_PFlowJets->begin();
    xAOD::JetContainer::const_iterator PFlowJets_end = m_PFlowJets->end();
    for( ; PFlowJets_itr != PFlowJets_end; ++PFlowJets_itr ) {
	int index = std::distance(m_PFlowJets->begin(),PFlowJets_itr);
	Info("PrintJetCollections", "PFlowJets E  = %.2f GeV  pt =  %.2f  eta = %.2f  phi = %.2f GeV",
	     (*PFlowJets_itr)->e()/GEV,
	     (*PFlowJets_itr)->pt()/GEV,
	     (*PFlowJets_itr)->eta(),
	     (*PFlowJets_itr)->phi());
    }
  }
  
  return;
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


//Could we return vector< pair<int,int> >
void xAODPFlowAna :: MatchJetCollections (const xAOD::JetContainer* TopoJet, const xAOD::JetContainer* PFlowJet) {
  
  Info("", "--- MatchJetCollections ---");
  
  // Matrix to store all DeltaR values

  std::vector<std::vector<float> > matrix_DeltaR;
  matrix_DeltaR.resize(TopoJet->size()); 
  
  xAOD::JetContainer::const_iterator topojet_itr = TopoJet->begin();
  xAOD::JetContainer::const_iterator topojet_end = TopoJet->end();
    
  for(; topojet_itr != topojet_end; topojet_itr++){
    int topojet_index = std::distance(TopoJet->begin(),topojet_itr);
    matrix_DeltaR[topojet_index].resize(PFlowJet->size());

    xAOD::JetContainer::const_iterator pflowjet_itr = PFlowJet->begin();
    xAOD::JetContainer::const_iterator pflowjet_end = PFlowJet->end();

    for(; pflowjet_itr != pflowjet_end; pflowjet_itr++){
      int pflowjet_index = std::distance(PFlowJet->begin(),pflowjet_itr);
      matrix_DeltaR[topojet_index][pflowjet_index] = ((*topojet_itr)->p4()).DeltaR((*pflowjet_itr)->p4());
    }
  }
  
  // Read Matrix
  for(int i=0;i<TopoJet->size();i++){for(int j=0;j<PFlowJet->size();j++){std::cout<<"i: "<<i<<"  j: "<<j<<"  deltaR: "<<matrix_DeltaR[i][j]<<std::endl;}}

  
  double DeltaRCut = 0.3; 
  std::pair<int,int> MatchedPair_pair;
  std::vector< std::pair<int,int> > MatchedPair_vector;

  for(int i = 0; i< TopoJet->size(); i++){
    int topojet =999;
    int pflowjet   =999;
    float DeltaRMin = 999;
    
    for(int j=0; j<PFlowJet->size(); j++){
      if( matrix_DeltaR[i][j] < DeltaRCut ){
	DeltaRMin = matrix_DeltaR[i][j];
	topojet = i;
	pflowjet   = j;
      }
    }
    
    std::cout<<"DeltaRMin: "<<DeltaRMin<<"i: "<<topojet<<" j:"<<pflowjet<<std::endl;
    
    MatchedPair_pair.first = topojet;
    MatchedPair_pair.second = pflowjet ;
    MatchedPair_vector.push_back(MatchedPair_pair);
    
    //eliminate the colum  
    for(int j=0; j<PFlowJet->size(); j++){
      if(i == topojet) matrix_DeltaR[topojet][j] = 999; 
      if(j == pflowjet) matrix_DeltaR[i][pflowjet] = 999; 
    }
    
    for(int i=0;i<TopoJet->size();i++){for(int j=0;j<PFlowJet->size();j++){std::cout<<"i: "<<i<<"  j: "<<j<<"  deltaR: "<<matrix_DeltaR[i][j]<<std::endl;}}
  }

  return;
}

  
  
 
bool xAODPFlowAna :: HasPFlowJetMatched (const xAOD::Jet& jet){

  return true;
  
}

int xAODPFlowAna :: WhichPFlowJetMatched(const xAOD::Jet& jet){

  return 0;
}




