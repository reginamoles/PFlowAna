#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ScanDir.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <TSystem.h>


#include "PFlowAna/xAODPFlowAna.h"

int main( int argc, char* argv[] ) {

  // Take the submit directory from the input if provided:
  std::string submitDir = "submitDir";
  if( argc > 1 ) submitDir = argv[ 1 ];

  // Set up the job for xAOD access:
  xAOD::Init().ignore();

  // Construct the samples to run on:
  // use SampleHandler to scan all of the subdirectories of a directory for particular MC single file:
  SH::SampleHandler sh;

  //DiJetSample
  //const char* inputFilePath = gSystem->ExpandPathName ("$DataFiles/JZ3/user.moles.147913.Pythia8_AU2CT10_jetjet_JZ3W.recon.AOD.e3099_s2832_r7617_AOD.67276209/");
  //SH::ScanDir().filePattern("user.moles.7760083.AOD._000102.pool.root").scan(sh,inputFilePath); //One indiviudual file
  //SH::ScanDir().scan(sh,inputFilePath); //All files in the directory

  //SinglePions
  //const char* inputFilePath = gSystem->ExpandPathName ("$DataFiles/SinglePions/mc15c_piplus/user.moles.mc15_13TeV.428001.ParticleGun_single_piplus_logE0p2to2000.recon.AOD.e3501_s2832_r8014_AOD/");
  //SH::ScanDir().filePattern("user.moles.8671065.AOD._000411.pool.root").scan(sh,inputFilePath); //One indiviudual file

  //Zmumu sample
  const char* inputFilePath = gSystem->ExpandPathName ("$DataFiles/Zmumu/mc15_13TeV.361107.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zmumu.merge.DAOD_JETM3.e3601_s2576_s2132_r7725_r7676_p2666_tid08619356_00/");
  SH::ScanDir().filePattern("DAOD_JETM3.08619356._000027.pool.root.1").scan(sh,inputFilePath); //One indiviudual file

  
  //const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/user/z/zhangr/work/eflowRec/r19.05-53/Run");
  //SH::ScanDir().filePattern("AOD.pool.root").scan(sh,inputFilePath); //One indiviudual file

  //SH::DiskListLocal list("./../../../test/run/AOD3/");
  //SH::scanFiles(sh, list);

  
  // Set the name of the input TTree. It's always "CollectionTree"
  // for xAOD files.
  sh.setMetaString( "nc_tree", "CollectionTree" );

  // Print what we found:
  sh.print();

  // Create an EventLoop job:
  EL::Job job;
  job.sampleHandler( sh );
  job.options()->setDouble (EL::Job::optSkipEvents, 0);
  job.options()->setDouble (EL::Job::optMaxEvents, -1);
  
  // Chose here your config options:
  bool SinglePionLowPerformanceStudies = false;
  bool DijetLowPerformance = false;
  bool DijetSubtraction = false;
  bool Zmumu = true;
  std::string matchScheme = (std::string)"_EM2";
  bool UseNarrowPtRange = true;
  bool UseNarrowEtaRange = true;
  bool PrintDebug = false; // Printing message criteria in PFlowAna 

  
  xAODPFlowAna* alg = new xAODPFlowAna(SinglePionLowPerformanceStudies, DijetLowPerformance, DijetSubtraction, Zmumu, matchScheme, UseNarrowPtRange, UseNarrowEtaRange, PrintDebug);
  job.algsAdd( alg );

  // Run the job using the local/direct driver:
  EL::DirectDriver driver;
  driver.submit( job, submitDir );

  return 0;
}
