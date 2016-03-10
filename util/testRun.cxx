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
  SH::SampleHandler sh;

  // use SampleHandler to scan all of the subdirectories of a directory for particular MC single file:
  //const char* inputFilePath = gSystem->ExpandPathName ("$ALRB_TutorialData/p2452/");
  //SH::ScanDir().filePattern("DAOD_SUSY1.07053638._000013.pool.root.1").scan(sh,inputFilePath);

  SH::DiskListLocal list("./../../../../../workspace/MC/AOD/JZ3/user.moles.147913.Pythia8_AU2CT10_jetjet_JZ3W.recon.AOD.e3099_s2832_r7617_AOD.67276\
209/");
  SH::scanFiles(sh, list);



  // Set the name of the input TTree. It's always "CollectionTree"
  // for xAOD files.
  sh.setMetaString( "nc_tree", "CollectionTree" );

  // Print what we found:
  sh.print();

  // Create an EventLoop job:
  EL::Job job;
  job.sampleHandler( sh );
  job.options()->setDouble (EL::Job::optMaxEvents, 500);

  // Add our analysis to the job:
  xAODPFlowAna* alg = new xAODPFlowAna();
  job.algsAdd( alg );

  // Run the job using the local/direct driver:
  EL::DirectDriver driver;
  driver.submit( job, submitDir );

  return 0;
}
