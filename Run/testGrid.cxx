#include "xAODRootAccess/Init.h"
#include "SampleHandler/SampleHandler.h"
#include "SampleHandler/ScanDir.h"
#include "SampleHandler/ToolsDiscovery.h"
#include "EventLoop/Job.h"
#include "EventLoop/DirectDriver.h"
#include "SampleHandler/DiskListLocal.h"
#include <TSystem.h>


#include "PFlowAna/xAODPFlowAna.h"

void testGrid(const std::string& submitDir) {


  // Set up the job for xAOD access:
  xAOD::Init().ignore();

  // use SampleHandler to scan all of the subdirectories of a directory for particular MC single file:
  SH::SampleHandler sh;

  SH::scanRucio (sh, "mc15_13TeV.370900.MadGraphPythia8EvtGen_A14NNPDF23LO_GG_direct_200_0.merge.DAOD_SUSY1.e4008_a766_a821_r7676_p2666/");

  
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
  driver.options()->setString("nc_outputSampleName", "user.zhangr.PflowAna.%in:name[2]%.%in:name[6]%");


  return;
}
