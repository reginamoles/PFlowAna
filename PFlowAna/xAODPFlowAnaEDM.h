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
//Electrons
#include "xAODEgamma/ElectronContainer.h"
//Muons
#include "xAODMuon/MuonContainer.h"

//Create deep and shallow copies
#include "xAODJet/JetAuxContainer.h"
#include "xAODCore/ShallowCopy.h"
#include "xAODCore/ShallowAuxContainer.h"
#include "xAODCore/AuxContainerBase.h"

//--------------------
// CP Tools includes
//--------------------
#include "JetSelectorTools/JetCleaningTool.h"
#include "JetCalibTools/JetCalibrationTool.h"
#include "JetResolution/JERTool.h"
#include "JetResolution/JERSmearingTool.h"
#include "MuonMomentumCorrections/MuonCalibrationAndSmearingTool.h"
#include <ElectronPhotonFourMomentumCorrection/EgammaCalibrationAndSmearingTool.h>
#include <IsolationSelection/IsolationSelectionTool.h>
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
#include "TrigConfxAOD/xAODConfigTool.h"
#include "TrigDecisionTool/TrigDecisionTool.h"
#include "JetJvtEfficiency/JetJvtEfficiency.h"
#include "PATInterfaces/CorrectionCode.h"
#include "PATInterfaces/SystematicRegistry.h"
#include "PATInterfaces/SystematicVariation.h" 
#include "PATInterfaces/SystematicsUtil.h"


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

const xAOD::ElectronContainer* m_Electrons;
const xAOD::MuonContainer* m_Muons;

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

//Tool for Muon-Calibration (+systematics)
CP::MuonCalibrationAndSmearingTool *m_muonCalibrationAndSmearingTool; 
std::vector<CP::SystematicSet> m_sysList; 

//Tool for electron calibration
CP::IEgammaCalibrationAndSmearingTool *m_electronCalibrationAndSmearingTool;
//Isolation Tool
CP::IsolationSelectionTool *m_iso; 
//AsgElectronLikelihoodTool, name refers to property chosen
AsgElectronLikelihoodTool *m_MediumLH;
AsgElectronLikelihoodTool *m_VeryLooseLHElectron;
	
//trigger tools
Trig::TrigDecisionTool *m_trigDecisionTool; 
TrigConf::xAODConfigTool *m_trigConfigTool; 

//JVT
CP::JetJvtEfficiency *m_jetsf; 


#endif

