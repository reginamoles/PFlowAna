#include <PFlowAna/xAODPFlowAna.h>


////////////////////////////////////////////////
// Calculate the objects wich minumum DeltaR
////////////////////////////////////////////////
void xAODPFlowAna :: CalculateMatrix_MinDeltaR ( const xAOD::TruthParticleContainer* TruthParticles, const xAOD::PFOContainer* JetETMissChargedParticleFlowObjects, float DeltaRCut){
  
  // Matrix
  std::vector<std::vector<float> > matrix_DeltaR;
  matrix_DeltaR.resize(TruthParticles->size());
  
  xAOD::TruthParticleContainer::const_iterator tp_itr = TruthParticles->begin();
  xAOD::TruthParticleContainer::const_iterator tp_end = TruthParticles->end();
  int tp_index = 0;
  for( ; tp_itr != tp_end; ++tp_itr ) {
    tp_index = std::distance(TruthParticles->begin(),tp_itr);
    matrix_DeltaR[tp_index].resize(JetETMissChargedParticleFlowObjects->size());  
    xAOD::PFOContainer::const_iterator cpfo_itr = JetETMissChargedParticleFlowObjects->begin();
    xAOD::PFOContainer::const_iterator cpfo_end = JetETMissChargedParticleFlowObjects->end();
    int cpfo_index = 0;
    for( ; cpfo_itr != cpfo_end; ++cpfo_itr ) {
      cpfo_index = std::distance(JetETMissChargedParticleFlowObjects->begin(),cpfo_itr);
      if((*tp_itr)->status()!=1) matrix_DeltaR[tp_index][cpfo_index] = 999;
      else matrix_DeltaR[tp_index][cpfo_index] = ((*tp_itr)->p4()).DeltaR((*cpfo_itr)->p4());
    }
  }
  
  //Read Matrix
  if(false){
    std::cout<<"=================="<<std::endl;
    std::cout<<"   DeltaR Matrix  "<<std::endl; 
    std::cout<<"=================="<<std::endl; 
    for(int i = 0 ; i<(int)TruthParticles->size();i++){
      for(int j = 0; j < (int)JetETMissChargedParticleFlowObjects->size(); j++){
	std::cout<<"i: "<<i<<"  j: "<<j<<"  deltaR: "<<matrix_DeltaR[i][j]<<std::endl;
      }
    }
  }
  
  std::pair<int,int> PairMatched;
  // The cPFO is ordered by pt, then we will loop first on it because we want to associate the higher pt tracks first.
  // WIP: This is not true, they are not associated by pt --> Have a look into it!
  for(int j=0; j<(int)JetETMissChargedParticleFlowObjects->size(); j++){
    int tp_min = 999999999;
    int cpfo_min = 999999999; 
    float DeltaRMin = DeltaRCut;
    
    for(int i = 0; i< (int)TruthParticles->size(); i++){
      if( matrix_DeltaR[i][j] < DeltaRMin ){
	DeltaRMin = matrix_DeltaR[i][j];
   	tp_min = i;
	cpfo_min = j;
      }
    }
    
    PairMatched.first = tp_min;
    PairMatched.second = cpfo_min ;
    _mc_MinDeltaREflowTrackPair.push_back(PairMatched);
    
    // You should eliminate the column to not match twice the same track
    for(int i = 0; i< (int)TruthParticles->size(); i++){
      if(i == tp_min) matrix_DeltaR[tp_min][j] = 999; 
      if(j == cpfo_min) matrix_DeltaR[i][cpfo_min] = 999; 
    }
  }
  
  //Read PairVector
  if(false){
    
    std::cout<<"=========================="<<std::endl;
    std::cout<<"   Pair Vector            "<<std::endl; 
    std::cout<<"=========================="<<std::endl; 
    
    for(int j = 0; j < (int)JetETMissChargedParticleFlowObjects->size(); j++){
      std::cout<<"MinDeltaRPair: cPFO = "<<  (_mc_MinDeltaREflowTrackPair.at(j)).second<<"  tp: "<< (_mc_MinDeltaREflowTrackPair.at(j)).first<<std::endl;
    }
  }
  
  return ;
}

bool xAODPFlowAna :: AreBothTracksMatched (int tp_index,int cpfo_index){
  bool AreMatched = false;
  
  for(int j = 0; j <(int)_mc_MinDeltaREflowTrackPair.size(); j++){
    if (((_mc_MinDeltaREflowTrackPair.at(j)).first) == tp_index && ((_mc_MinDeltaREflowTrackPair.at(j)).second) == cpfo_index) {
      AreMatched = true;
    }
  }
  return AreMatched;
}



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




