
#include <PFlowAna/xAODPFlowAna.h>

//////////////////////////////////////////////////////////////////////////
// Functions related with the jet matching between different collections
/////////////////////////////////////////////////////////////////////////


//Could we return vector< pair<int,int> >
void xAODPFlowAna :: MatchJetCollections (const xAOD::JetContainer* TopoJet, const xAOD::JetContainer* PFlowJet) {
  
  Info("", "--- MatchJetCollections ---");

  //Just to not have warning messages. It will be removed asap
  Info("PrintJetCollections", "Number of TopoJets = %lu", TopoJet->size());
  Info("PrintJetCollections", "Number of PFlowJets = %lu", PFlowJet->size());
  
  /*
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
  */
  return;
}

 
bool xAODPFlowAna :: HasPFlowJetMatched (const xAOD::Jet& jet){

  //just for avoiding the warning messages
  Info("PrintJetCollections", "jet pt = %f", jet.pt());

  return true;
  
}

int xAODPFlowAna :: WhichPFlowJetMatched(const xAOD::Jet& jet){
  //just for avoiding the warning messages
  Info("PrintJetCollections", "jet pt = %f", jet.pt());
  return 0;
}




