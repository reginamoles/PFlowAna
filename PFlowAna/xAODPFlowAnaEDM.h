///////////////////////// -*- C++ -*- /////////////////////////////
// This file has EDM includes, EDM collections and CP tools for xExample.
///////////////////////////////////////////////////////////////////
#ifndef PFlowAna_xAODPFlowAnaEDM_H
#define PFlowAna_xAODPFlowAnaEDM_H

//----------------
// EDM includes:
//----------------
#include "xAODEventInfo/EventInfo.h"
//Truth includes
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthVertexContainer.h"
//TrackParticle Includes
#include "xAODTracking/TrackParticle.h"
#include "xAODTracking/TrackParticleContainer.h"
//Jets and clusters
#include "xAODJet/JetContainer.h"
#include "xAODCaloEvent/CaloCluster.h"
#include "xAODCaloEvent/CaloClusterContainer.h"
//PFO
#include "xAODPFlow/PFOContainer.h"
#include "xAODPFlow/PFO.h"
//CalCellInfo include
//#include "xAODCalCellInfo/CalCellInfo.h"
//#include "xAODCalCellInfo/CalCellInfoContainer.h"


//----------------
// EDM containers
//----------------
const xAOD::TruthParticleContainer* TruthParticles;
const xAOD::TruthVertexContainer* TruthVertices;
const xAOD::TrackParticleContainer* InDetTrackParticles;

const xAOD::CaloClusterContainer* topocluster;
const xAOD::PFOContainer* JetETMissChargedParticleFlowObjects;
const xAOD::PFOContainer* JetETMissNeutralParticleFlowObjects;

const xAOD::CaloClusterContainer* PFOcluster;

//const xAOD::CalCellInfoContainer* CalCellInfo_TopoCluster;
//const xAOD::CalCellInfoContainer* CalCellInfo;

const xAOD::JetContainer* m_jets;

#endif
