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
  //const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/work/z/zhangr/eflowRec/data/matchingStudy/user.zhangr.mc15_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.recon.AOD.e3569_s2832_r8014_AOD.94973736/");
  //SH::ScanDir().filePattern("user.zhangr.9335218.AOD._000001.pool.root").scan(sh,inputFilePath); //One indiviudual file
  const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/work/z/zhangr/eflowRec/r19.05-53/Run/");
  SH::ScanDir().filePattern("AOD.pool1.root").scan(sh,inputFilePath); //One indiviudual file

  //SinglePions
  //const char* inputFilePath = gSystem->ExpandPathName ("$DataFiles/SinglePions/mc15c_piplus/user.moles.mc15_13TeV.428001.ParticleGun_single_piplus_logE0p2to2000.recon.AOD.e3501_s2832_r8014_AOD/");
  //SH::ScanDir().filePattern("user.moles.8671065.AOD._000411.pool.root").scan(sh,inputFilePath); //One indiviudual file
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
  job.options()->setDouble (EL::Job::optMaxEvents, 1);

  // Add our analysis to the job:
  // SinglePionLowPerformanceStudies, DijetLowPerformance, DijetSubtraction, Zmumu
  xAODPFlowAna* alg = new xAODPFlowAna(false, false, true, false, true);
  //xAODPFlowAna* alg = new xAODPFlowAna(true, false, false, false);
  job.algsAdd( alg );

  // Run the job using the local/direct driver:
  EL::DirectDriver driver;
  driver.submit( job, submitDir );

  return 0;
}
