#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ScanDir.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoopGrid/PrunDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <TSystem.h>

#include "PFlowAna/xAODPFlowAna.h"

int main( int argc, char* argv[] ) {

  // Take the submit directory from the input if provided:
  std::string submitDir = "submitDir";
  std::string path = "";
  std::string output = "";
  if( argc > 1 ) submitDir = argv[ 1 ];
  if( argc > 2 ) path = argv[ 2 ];
  if( argc > 3 ) output = argv[ 3 ];

  // Set up the job for xAOD access:
  xAOD::Init().ignore();

  // use SampleHandler to scan all of the subdirectories of a directory for particular MC single file:
  SH::SampleHandler sh;

  //SH::scanRucio (sh, "user.zhangr.mc15_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.recon.AOD.e3569_s28r8014_v4_AOD/");
  SH::scanRucio (sh, path);

  
  // Set the name of the input TTree. It's always "CollectionTree"
  // for xAOD files.
  sh.setMetaString( "nc_tree", "CollectionTree" );

  // Print what we found:
  sh.print();

  // Create an EventLoop job:
  EL::Job job;
  job.sampleHandler( sh );
  job.options()->setDouble (EL::Job::optSkipEvents, 0);
  job.options()->setDouble (EL::Job::optMaxEvents, 2000);

  // Add our analysis to the job:
  // SinglePionLowPerformanceStudies, DijetLowPerformance, DijetSubtraction, Zmumu
  xAODPFlowAna* alg = new xAODPFlowAna(false, false, true, false, true);
  job.algsAdd( alg );

  // Run the job using the local/direct driver:
  //EL::DirectDriver driver;
  EL::PrunDriver driver;
  driver.submit( job, submitDir );
  //driver.options()->setString("nc_outputSampleName", "user.zhangr.PflowAna.%in:name[2]%.%in:name[6]%");
  driver.options()->setString("nc_outputSampleName", output);


  return 0;
}
