#ifndef PFlowAna_xAODPFlowAna_H
#define PFlowAna_xAODPFlowAna_H

#include <EventLoop/Algorithm.h>

#include "xAODRootAccess/TEvent.h"
#include "JetSelectorTools/JetCleaningTool.h"
#include "xAODJet/JetContainer.h"

class xAODPFlowAna : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
 public:
  // float cutValue;
  
  
  
  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
  
  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();
  xAODPFlowAna ();


 
  
 private:
  
  float GEV; //!
  float PrintDebug;//!
  
  xAOD::TEvent *m_event;//!
  int m_eventCounter; //!
  
  JetCleaningTool *m_jetCleaning; //!  

  
  // Tree *myTree; //!
  // TH1 *myHist; //!

  //Printing functions 
  void PrintTruthInfo();
  void PrintTrackInfo ();
  void PrintPFOInfo();
  void PrintClusterInfo();
  void PrintCalCellInfo();
  void PrintJetCollections();
  
  //BadJetsScan
  void BadJetsScan(const xAOD::Jet& jet);
  void MatchJetCollections(const xAOD::JetContainer*, const xAOD::JetContainer* );
  bool HasPFlowJetMatched(const xAOD::Jet& jet); //return a true is has been matched
  int WhichPFlowJetMatched(const xAOD::Jet& jet); //return the index of the PFlowJet matched
  
public:

  // this is needed to distribute the algorithm to the workers
  ClassDef(xAODPFlowAna, 1);
};

#endif
