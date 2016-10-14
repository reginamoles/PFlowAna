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
  std::string path = "";
  std::string input = "";
  if( argc > 1 ) submitDir = argv[ 1 ];
  if( argc > 2 ) path = argv[ 2 ];
  if( argc > 3 ) input = argv[ 3 ];

  // Set up the job for xAOD access:
  xAOD::Init().ignore();

  // Construct the samples to run on:
  // use SampleHandler to scan all of the subdirectories of a directory for particular MC single file:
  SH::SampleHandler sh;
  //DiJetSample
  //const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/user/z/zhangr/work/eflowRec/data/JZ1W/user.zhangr/");
  const char* inputFilePath = gSystem->ExpandPathName (path.c_str());
  //SH::ScanDir().filePattern("user.zhangr.9367862.AOD._000097.pool.root").scan(sh,inputFilePath);
  SH::ScanDir().filePattern(input).scan(sh,inputFilePath);


  // Set the name of the input TTree. It's always "CollectionTree"
  // for xAOD files.
  sh.setMetaString( "nc_tree", "CollectionTree" );

  // Print what we found:
  sh.print();

  // Create an EventLoop job:
  EL::Job job;
  job.sampleHandler( sh );
  job.options()->setDouble (EL::Job::optSkipEvents, 0);
  job.options()->setDouble (EL::Job::optMaxEvents, 100);

  // Add our analysis to the job:
  // SinglePionLowPerformanceStudies, DijetLowPerformance, DijetSubtraction, Zmumu
  xAODPFlowAna* alg = new xAODPFlowAna(false, false, true, false, true, submitDir);
  std::cout<<"submitDir="<<submitDir<<std::endl;
  //xAODPFlowAna* alg = new xAODPFlowAna(true, false, false, false);
  job.algsAdd( alg );

  // Run the job using the local/direct driver:
  EL::DirectDriver driver;
  driver.submit( job, submitDir );

  return 0;
}
