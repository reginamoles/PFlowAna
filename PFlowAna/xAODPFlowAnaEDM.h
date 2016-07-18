///////////////////////// -*- C++ -*- /////////////////////////////////////
// This file has EDM includes, EDM collections and CP tools for PFlowAna.
//////////////////////////////////////////////////////////////////////////
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
#include "xAODCalCellInfo/CalCellInfo.h"
#include "xAODCalCellInfo/CalCellInfoContainer.h"

//Create deep and shallow copies
#include "xAODJet/JetAuxContainer.h"
#include "xAODCore/ShallowCopy.h"

//--------------------
// CP Tools includes
//--------------------
#include "JetSelectorTools/JetCleaningTool.h"
#include "JetCalibTools/JetCalibrationTool.h"
#include "JetResolution/JERTool.h"
#include "JetResolution/JERSmearingTool.h"

//----------------
// EDM containers
//----------------
const xAOD::TruthParticleContainer* m_TruthParticles;
const xAOD::TruthVertexContainer* m_TruthVertices;
const xAOD::TrackParticleContainer* m_InDetTrackParticles;

const xAOD::CaloClusterContainer* m_topocluster;
const xAOD::PFOContainer* m_JetETMissChargedParticleFlowObjects;
const xAOD::PFOContainer* m_JetETMissNeutralParticleFlowObjects;

const xAOD::CaloClusterContainer* m_PFOcluster;

const xAOD::CalCellInfoContainer* m_CalCellInfo_TopoCluster;
const xAOD::CalCellInfoContainer* m_CalCellInfo;

const xAOD::JetContainer* m_Jets;
const xAOD::JetContainer* m_PFlowJets;

// Transient object store. Needed for the CP tools.
xAOD::TStore* m_store;



//----------------
// CP Tools
//----------------
JetCleaningTool *m_jetCleaning; 
JetCalibrationTool* m_akt4EMTopoCalibrationTool;
JetCalibrationTool* m_akt4EMPFlowCalibrationTool;

JERTool *m_JERTool; 
JERSmearingTool *m_SmearTool; 




#endif
