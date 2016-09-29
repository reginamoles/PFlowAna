#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ScanDir.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <TSystem.h>
#include "EventLoopAlgs/DuplicateChecker.h"

#include "PFlowAna/xAODPFlowAna.h"

int main( int argc, char* argv[] ) {

  // Take the submit directory from the input if provided:
  std::string submitDir = "submitDir";
  if( argc > 1 ) submitDir = argv[ 1 ];

  // Set up the job for xAOD access:
  xAOD::Init().ignore();

  //=================================
  // Chose here your config options:
  //=================================
  bool SinglePionLowPerformanceStudies = false;
  bool DijetLowPerformance = false;
  bool DijetSubtraction = false;
  bool Zmumu = true;
  //data or MC
  bool data = false;
  bool MC = !data;
  std::string matchScheme = (std::string)"_EM2";
  bool UseNarrowPtRange = false;
  bool UseNarrowEtaRange = false;
  bool PrintDebug = false; // Printing message criteria in PFlowAna 


  // Construct the samples to run on:
  // use SampleHandler to scan all of the subdirectories of a directory for particular MC single file:
  SH::SampleHandler sh;

  if(DijetLowPerformance || DijetSubtraction){
    const char* inputFilePath = gSystem->ExpandPathName ("$DataFiles/JZ3/user.moles.mc15_13TeV.361023.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3W.recon.AOD.e3668_s2832_r8014_AOD.92387207/");
    SH::ScanDir().filePattern("user.moles.9200087.AOD._008942.pool.root").scan(sh,inputFilePath); //One indiviudual file

    //test File
    //const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/user/m/moles/pflow/Performance/20.7.5/run/AOD/");
    //SH::ScanDir().filePattern("AOD.pool.root").scan(sh,inputFilePath); //One indiviudual file
    
  }
    
  if(SinglePionLowPerformanceStudies){
    const char* inputFilePath = gSystem->ExpandPathName ("$DataFiles/MyProduction_Subtraction_SingleChargedPions_WholeSample/");
    SH::ScanDir().filePattern("user.moles.6424713.AOD._000153.pool.root").scan(sh,inputFilePath); //One indiviudual file
    //const char* inputFilePath = gSystem->ExpandPathName ("$DataFiles/SinglePions/mc15c_piplus/user.moles.mc15_13TeV.428001.ParticleGun_single_piplus_logE0p2to2000.recon.AOD.e3501_s2832_r8014_AOD/");
    //SH::ScanDir().filePattern("user.moles.8671065.AOD._000411.pool.root").scan(sh,inputFilePath); //One indiviudual file
  }
  
  if(Zmumu){
	if(MC){	
		const char* inputFilePath = gSystem->ExpandPathName ("~/Desktop/");
		SH::ScanDir().filePattern("DAOD_JETM3.08619382._000015.pool.root.1").scan(sh,inputFilePath); //One indiviudual file
	}
	if(data){
		const char* inputFilePath = gSystem->ExpandPathName ("~/Desktop/data/");
		SH::ScanDir().filePattern("DAOD_JETM3.08654688._000082.pool.root.1").scan(sh,inputFilePath); //One indiviudual file
	}
  }
  
  
  //const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/user/z/zhangr/work/eflowRec/r19.05-53/Run");
  //SH::ScanDir().filePattern("AOD.pool.root").scan(sh,inputFilePath); //One indiviudual file
  //SH::ScanDir().scan(sh,inputFilePath); //All files in the directory
  
  //SH::DiskListLocal list("./../../../test/run/AOD3/");
  //SH::scanFiles(sh, list);

  
  // Set the name of the input TTree. It's always "CollectionTree" for xAOD files.
  sh.setMetaString( "nc_tree", "CollectionTree" );

  // Print what we found:
  sh.print();

  // Create an EventLoop job:
  EL::Job job;
  job.sampleHandler( sh );
  job.options()->setDouble (EL::Job::optSkipEvents, 0);
  job.options()->setDouble (EL::Job::optMaxEvents, 10000);

  //-------------------------
  // Check duplicates
  //-------------------------
  EL::DuplicateChecker* duplicates = new EL::DuplicateChecker();
  duplicates->setOutputTreeName ("duplicate_info");
  job.algsAdd (duplicates);
  
  
  xAODPFlowAna* alg = new xAODPFlowAna(data, SinglePionLowPerformanceStudies, DijetLowPerformance, DijetSubtraction, Zmumu, matchScheme, UseNarrowPtRange, UseNarrowEtaRange, PrintDebug);
  job.algsAdd( alg );

  // Run the job using the local/direct driver:
  EL::DirectDriver driver;
  driver.submit( job, submitDir );

  return 0;
}
