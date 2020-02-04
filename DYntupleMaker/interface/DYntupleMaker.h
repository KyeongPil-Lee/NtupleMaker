#ifndef DYntupleMaker_H
#define DYntupleMaker_H

////////////////////////////////
// -- system include files -- //
////////////////////////////////
#include <memory>
#include <iostream>

//////////////////////
// -- FrameWorks -- //
//////////////////////
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

/////////////////////
// -- For Muons -- //
/////////////////////
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h"

/////////////////////////
// -- For Electrons -- //
/////////////////////////
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"

///////////////////
// -- For MET -- //
///////////////////
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/PFMET.h"

////////////////////
// -- For Jets -- //
////////////////////
#include "DataFormats/PatCandidates/interface/Jet.h"

////////////////////////////
// -- For GenParticles -- //
////////////////////////////
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

//////////////////////////
// -- Track & Vertex -- //
//////////////////////////
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

////////////////////
// -- Triggers -- //
////////////////////
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

////////////////
// -- Else -- //
////////////////
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/GeometrySurface/interface/GloballyPositioned.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "DataFormats/Common/interface/RefToBaseVector.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/NamedCompositeCandidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Common/interface/View.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFrame.h>
#include <TMath.h>

using namespace std;
using namespace pat;
namespace reco { class CandCommonVertexFitterBase; class VertexCompositeCandidate; class CandCommonVertexFitter; }

class DYntupleMaker : public edm::EDAnalyzer
{
public:
	explicit DYntupleMaker(const edm::ParameterSet&);
	~DYntupleMaker();

private:
	virtual void beginJob() ;
	virtual void analyze(const edm::Event&, const edm::EventSetup&);
	virtual void endJob() ;
	virtual void beginRun(const edm::Run &, const edm::EventSetup & );
	virtual void endRun(const edm::Run &, const edm::EventSetup & );

	virtual void fillPrimaryVertex(const edm::Event &iEvent);  // fill primary vertex information
	virtual void fillMET(const edm::Event &iEvent);            // fill MET information
	virtual void fillPhotons(const edm::Event &iEvent);
	virtual void fillMuons(const edm::Event &iEvent, const edm::EventSetup& iSetup);
	virtual void fillElectrons(const edm::Event &iEvent, const edm::EventSetup& iSetup);
	virtual void fillJet(const edm::Event &iEvent);            // fill jet and b-tagging information
	virtual void hltReport(const edm::Event &iEvent);          // fill list of triggers fired in an event
	virtual void fillLHEInfo(const edm::Event &iEvent);
	virtual void fillGENInfo(const edm::Event &iEvent);            // fill MET information
	virtual void fillGenOthersInfo(const edm::Event &iEvent);
	virtual void fillTT(const edm::Event&);

	bool reorder(double &a, double &b)
	{
		return a > b;
	}

	int theDebugLevel;                   // 0 no prints, 1 some, 2 lots
	std::string processName;
	std::string theElectronID;

	HLTConfigProvider hltConfig_;

	// -- Tokens (for 76X) -- //
	edm::EDGetTokenT< std::vector<pat::Muon> > 						MuonToken;
	edm::EDGetTokenT< edm::View<reco::GsfElectron> > 					ElectronToken;
	edm::EDGetTokenT< edm::View<reco::GsfElectron> > 					UnCorrElectronToken;
	edm::EDGetTokenT< edm::View<reco::Photon> > 						PhotonToken;
	edm::EDGetTokenT< std::vector<pat::Jet> > 						JetToken;
	edm::EDGetTokenT< std::vector<pat::MET> > 						MetToken;
	edm::EDGetTokenT< LHEEventProduct > 							LHEEventProductToken;
	edm::EDGetTokenT< LHERunInfoProduct >							LHERunInfoProductToken;
	edm::EDGetTokenT< std::vector<reco::GenParticle> > 					GenParticleToken;

	edm::EDGetTokenT< double > 								RhoToken;
	edm::EDGetTokenT< edm::ValueMap<bool> > 						eleVetoIdMapToken;
	edm::EDGetTokenT< edm::ValueMap<bool> > 						eleLooseIdMapToken;
	edm::EDGetTokenT< edm::ValueMap<bool> > 						eleMediumIdMapToken;
	edm::EDGetTokenT< edm::ValueMap<bool> > 						eleTightIdMapToken;
	edm::EDGetTokenT< edm::ValueMap<bool> > 						eleMVAIdWP80MapToken;
	edm::EDGetTokenT< edm::ValueMap<bool> > 						eleMVAIdWP90MapToken;
	edm::EDGetTokenT< edm::ValueMap<bool> > 						eleHEEPIdMapToken;
	edm::EDGetTokenT< std::vector<reco::Conversion> > 					ConversionsToken;
	edm::EDGetTokenT< std::vector< reco::GsfTrack > > 					GsfTrackToken;

	//edm::EDGetTokenT< edm::ValueMap<float> > 						full5x5SigmaIEtaIEtaMapToken;
	edm::EDGetTokenT< edm::ValueMap<float> > 						phoChargedIsolationToken;
	edm::EDGetTokenT< edm::ValueMap<float> > 						phoNeutralHadronIsolationToken;
	edm::EDGetTokenT< edm::ValueMap<float> > 						phoPhotonIsolationToken;
	edm::EDGetTokenT< edm::ValueMap<bool> > 						phoMediumIdMapToken;

	edm::EDGetTokenT< edm::TriggerResults > 						TriggerToken;
	edm::EDGetTokenT< edm::TriggerResults > 						TriggerTokenPAT;
	edm::EDGetTokenT< std::vector<pat::TriggerObjectStandAlone> > 				TriggerObjectToken;
	//edm::EDGetTokenT< trigger::TriggerEvent > 						TriggerSummaryToken;

	edm::EDGetTokenT< GenEventInfoProduct > 						GenEventInfoToken;
	edm::EDGetTokenT< reco::BeamSpot > 							BeamSpotToken;
	edm::EDGetTokenT< reco::VertexCollection > 						PrimaryVertexToken;
	edm::EDGetTokenT< edm::View<reco::Track> > 						TrackToken;
	edm::EDGetTokenT< std::vector< PileupSummaryInfo > > 					PileUpInfoToken;

	// -- for Level 1 ECAL prefiring -- //
	edm::EDGetTokenT< double > 								prefweight_token;
	edm::EDGetTokenT< double > 								prefweightup_token;
	edm::EDGetTokenT< double > 								prefweightdown_token;



	// // -- For Electrons -- //
	// edm::InputTag theElectronLabel;
	// edm::InputTag rho;
	// edm::InputTag eleVetoIdMap;
	// edm::InputTag eleLooseIdMap;
	// edm::InputTag eleMediumIdMap;
	// edm::InputTag eleTightIdMap;

	// // -- Photon information -- //
	// edm::InputTag thePhotonLabel;
	// edm::InputTag full5x5SigmaIEtaIEtaMapLabel; 
	// edm::InputTag phoChargedIsolationLabel; 
	// edm::InputTag phoNeutralHadronIsolationLabel; 
	// edm::InputTag phoPhotonIsolationLabel;     
	edm::FileInPath effAreaChHadronsFile;
	edm::FileInPath effAreaNeuHadronsFile;
	edm::FileInPath effAreaPhotonsFile;  

	// edm::InputTag theMuonLabel;
	// edm::InputTag theMETLabel;
	// edm::InputTag pfMetCollection_;
	// edm::InputTag theJetLabel;
	// edm::InputTag theTimeTags;
	// edm::InputTag theBeamSpot;
	// edm::InputTag theTrackLabel;
	// edm::InputTag theGenParticleLabel;
	// edm::InputTag thePrimaryVertexLabel;
	// edm::InputTag thePileUpInfoLabel;

	// edm::InputTag theConversionsInputTag;

	// -- Store flags -- // 
	bool theStorePriVtxFlag;                // Yes or No to store primary vertex
	bool theStoreJetFlag;                // Yes or No to store Jet
	bool theStoreMETFlag;                // Yes or No to store MET 
	bool theStoreHLTReportFlag;             // Yes or No to store HLT reuslts (list of triggers fired)
	bool theStoreMuonFlag;
	bool theStoreElectronFlag;
	bool theStoreLHEFlag;
	bool theStoreGENFlag;
	bool theStoreGenOthersFlag;
	bool theStorePhotonFlag;
	bool theStoreTTFlag;
	bool isMC;  //turn gen on and off
	bool theApplyFilter;
	int theFilterType;

	// double theCrossSection;
	// double theFilterEfficiency;
	// double theTotalNevents;
	double theBDiscriminant;
	// double theIntegratedLumi;
	double theBDiscriminant_alg1;
	double theBDiscriminant_alg2;
	double theBDiscriminant_alg3;
	int theNHighQualLeptons;
	
	edm::ESHandle<TransientTrackBuilder> theTTBuilder;

	vector< std::string> vec_inputHLTList_;

	std::vector<std::string > MuonHLT;
	std::vector<int > MuonHLTPS;
	std::vector<std::string > trigModuleNames;
	std::vector<std::string > trigModuleNames_preFil;

	// Pile-up Reweight
	// edm::LumiReWeighting LumiWeights_;
	// reweight::PoissonMeanShifter PShiftUp_;
	// reweight::PoissonMeanShifter PShiftDown_;
	// edm::LumiReWeighting LumiWeightsMuonPhys_;
	// reweight::PoissonMeanShifter PShiftUpMuonPhys_;
	// reweight::PoissonMeanShifter PShiftDownMuonPhys_;

	// std::vector<double> PileUpRD_;
	// std::vector<double> PileUpRDMuonPhys_;
	// std::vector<double> PileUpMC_;

	unsigned int nPileUp;
	double pileUpReweightIn;
	double pileUpReweight;
	double pileUpReweightPlus;
	double pileUpReweightMinus;

	double pileUpReweightInMuonPhys;
	double pileUpReweightMuonPhys;
	double pileUpReweightPlusMuonPhys;
	double pileUpReweightMinusMuonPhys;

	// -- Level 1 ECAL prefiring -- //
	double _prefiringweight;
	double _prefiringweightup;
	double _prefiringweightdown;

	TTree *DYTree;

	int nEvt;

	//FILL a tree
	static const int MPSIZE = 2000;

	// Invariant Mass distribution of SS(OS) di-muon events GG (2 global)
	// GlbTree
	int runNum;
	unsigned long long evtNum;
	int lumiBlock;
	int nMuon;
	double PUweight;
	double sumEt;
	double photonEt;
	double chargedHadronEt;
	double neutralHadronEt;

	// double MET_sumEt;
	// double MET_pt;
	// double MET_px;
	// double MET_py;
	// double MET_phi;
	// double pfMET_sumEt;
	// double pfMET_pt;
	// double pfMET_px;
	// double pfMET_py;
	// double pfMET_phi;
	int Njets;
	int Nelectrons;
	int Nmuons;
	int Nbtagged;
	int NbtaggedCloseMuon;

	// -- PDf weights -- //
	std::vector< double > PDFWeights;
	
	// -- Flags in re-miniAOD -- //
	bool Flag_duplicateMuons;
	bool Flag_badMuons;
	bool Flag_noBadMuons;

	// PV
	int nVertices;
	int PVtrackSize;
	double PVchi2;
	double PVndof;
	double PVnormalizedChi2;
	double PVx;
	double PVy;
	double PVz;
	double PVprob;

	// trigger object
	int _HLT_ntrig;
	int _HLT_trigType[MPSIZE];
	int _HLT_trigFired[MPSIZE];
	std::vector<std::string> _HLT_trigName;
	std::vector<int> _HLT_trigPS;
	double _HLT_trigPt[MPSIZE];
	double _HLT_trigEta[MPSIZE];
	double _HLT_trigPhi[MPSIZE];

	// Jet
	double JETbDiscriminant[MPSIZE];
	double JETcharge[MPSIZE];
	int JETflavour[MPSIZE];
	int JETntracks[MPSIZE];
	double JETpt[MPSIZE];
	double JETeta[MPSIZE];
	double JETphi[MPSIZE];

	int Nbtagged_alg1;
	int NbtaggedCloseMuon_alg1;
	int Nbtagged_alg2;
	int NbtaggedCloseMuon_alg2;
	int Nbtagged_alg3;
	int NbtaggedCloseMuon_alg3;
	double JETbDiscriminant_alg1[MPSIZE];
	double JETbDiscriminant_alg2[MPSIZE];
	double JETbDiscriminant_alg3[MPSIZE];

	double Jet_pT[MPSIZE]; 
	double Jet_eta[MPSIZE]; 
	double Jet_phi[MPSIZE]; 
	double Jet_Charge[MPSIZE]; 
	int Jet_Flavor[MPSIZE]; 

	double Jet_bTag[MPSIZE]; 
	double Jet_CHfrac[MPSIZE]; 
	double Jet_NHfrac[MPSIZE]; 
	double Jet_NHEMfrac[MPSIZE]; 
	double Jet_CHEMfrac[MPSIZE]; 
	int Jet_CHmulti[MPSIZE]; 
	int Jet_NHmulti[MPSIZE];

	
	// PAT Electron
	double Electron_et[MPSIZE];
	double Electron_caloEnergy[MPSIZE];
	double Electron_Energy[MPSIZE];
	double Electron_pT[MPSIZE];
	double Electron_Px[MPSIZE];
	double Electron_Py[MPSIZE];
	double Electron_Pz[MPSIZE];
	double Electron_eta[MPSIZE];
	double Electron_phi[MPSIZE];
	int Electron_charge[MPSIZE];
	double Electron_gsfpT[MPSIZE];
	double Electron_gsfPx[MPSIZE];
	double Electron_gsfPy[MPSIZE];
	double Electron_gsfPz[MPSIZE];
	double Electron_gsfEta[MPSIZE];
	double Electron_gsfPhi[MPSIZE];
	int Electron_gsfCharge[MPSIZE];
	double Electron_etaSC[MPSIZE];
	double Electron_phiSC[MPSIZE];
	double Electron_etaWidth[MPSIZE];
	double Electron_phiWidth[MPSIZE];
	double Electron_dEtaIn[MPSIZE];
	double Electron_dEtaInSeed[MPSIZE];
	double Electron_dPhiIn[MPSIZE];
	double Electron_sigmaIEtaIEta[MPSIZE];
	double Electron_Full5x5_SigmaIEtaIEta[MPSIZE];
	double Electron_HoverE[MPSIZE];
	double Electron_fbrem[MPSIZE];
	double Electron_eOverP[MPSIZE];
	double Electron_energyEC[MPSIZE];
	double Electron_Pnorm[MPSIZE];
	double Electron_InvEminusInvP[MPSIZE];
	double Electron_dxyVTX[MPSIZE];
	double Electron_dzVTX[MPSIZE];
	double Electron_dxy[MPSIZE];
	double Electron_dz[MPSIZE];
	double Electron_dxyBS[MPSIZE];
	double Electron_dzBS[MPSIZE];
	double Electron_AEff03[MPSIZE];
	double Electron_chIso03[MPSIZE];
	double Electron_chIso04[MPSIZE];
	double Electron_nhIso03[MPSIZE];
	double Electron_nhIso04[MPSIZE];
	double Electron_phIso03[MPSIZE];
	double Electron_phIso04[MPSIZE];
	double Electron_pcIso03[MPSIZE];
	double Electron_pcIso04[MPSIZE];
	double Electron_relIsoCom03[MPSIZE];
	double Electron_relIsoCom04[MPSIZE];
	double Electron_relIsoBeta03[MPSIZE];
	double Electron_relIsoBeta04[MPSIZE];
	double Electron_relIsoRho03[MPSIZE];
	bool Electron_hasConversion[MPSIZE];
	int Electron_mHits[MPSIZE];

	int Electron_crack[MPSIZE];
	int Electron_ecalDriven[MPSIZE];
	double Electron_isoEMHADDepth1[MPSIZE];
	double Electron_25over55[MPSIZE];
	double Electron_15over55[MPSIZE];
	double Electron_isoHADDepth2[MPSIZE];
	double Electron_isoPtTrks[MPSIZE];
	double Electron_modIsoEMHADDepth1[MPSIZE];
	double Electron_modIsoPtTrks[MPSIZE];
	double Electron_modIsoEMHADDepth1Orig[MPSIZE];
	double Electron_modIsoPtTrksOrig[MPSIZE];
	double Electron_ambGsf0Pt[MPSIZE];
	double Electron_ambGsf0Eta[MPSIZE];
	double Electron_ambGsf0Phi[MPSIZE];
	double Electron_ambGsf0Charge[MPSIZE];
	double Electron_ambGsf1Pt[MPSIZE];
	double Electron_ambGsf1Eta[MPSIZE];
	double Electron_ambGsf1Phi[MPSIZE];
	double Electron_ambGsf1Charge[MPSIZE];
	double Electron_ambGsf2Pt[MPSIZE];
	double Electron_ambGsf2Eta[MPSIZE];
	double Electron_ambGsf2Phi[MPSIZE];
	double Electron_ambGsf2Charge[MPSIZE];
	double Electron_ambGsf3Pt[MPSIZE];
	double Electron_ambGsf3Eta[MPSIZE];
	double Electron_ambGsf3Phi[MPSIZE];
	double Electron_ambGsf3Charge[MPSIZE];

	double Electron_r9[MPSIZE];
	double Electron_EnergySC[MPSIZE];
	double Electron_preEnergySC[MPSIZE];
	double Electron_rawEnergySC[MPSIZE];
	double Electron_etSC[MPSIZE];
	double Electron_E15[MPSIZE];
	double Electron_E25[MPSIZE];
	double Electron_E55[MPSIZE];
	double Electron_ChIso03FromPU[MPSIZE];
	double Electron_RelPFIso_dBeta[MPSIZE];
	double Electron_RelPFIso_Rho[MPSIZE];
	bool Electron_passConvVeto[MPSIZE];
	bool Electron_passVetoID[MPSIZE];
	bool Electron_passLooseID[MPSIZE];
	bool Electron_passMediumID[MPSIZE];
	bool Electron_passTightID[MPSIZE];
	bool Electron_passMVAID_WP80[MPSIZE];
	bool Electron_passMVAID_WP90[MPSIZE];
	bool Electron_passHEEPID[MPSIZE];

	std::vector<double> vtxTrkDiEChi2;
	std::vector<double> vtxTrkDiEProb;
	std::vector<double> vtxTrkDiENdof;
	std::vector<double> vtxTrkDiE1Pt;
	std::vector<double> vtxTrkDiE2Pt;

	// -- emu vertex -- //
	std::vector<double> vtxTrkEMuChi2;
	std::vector<double> vtxTrkEMuProb;
	std::vector<double> vtxTrkEMuNdof;
	std::vector<double> vtxTrkEMu1Pt;
	std::vector<double> vtxTrkEMu2Pt;
	std::vector<double> vtxTrkEMuChi2_TuneP;
	std::vector<double> vtxTrkEMuProb_TuneP;
	std::vector<double> vtxTrkEMuNdof_TuneP;
	std::vector<double> vtxTrkEMu1Pt_TuneP;
	std::vector<double> vtxTrkEMu2Pt_TuneP;

	// -- Un-corrected electrons -- //
	int nUnCorrElectron;
	double Electron_pTUnCorr[MPSIZE];
	double Electron_etaUnCorr[MPSIZE];
	double Electron_phiUnCorr[MPSIZE];
	double Electron_PxUnCorr[MPSIZE];
	double Electron_PyUnCorr[MPSIZE];
	double Electron_PzUnCorr[MPSIZE];
	double Electron_EnergyUnCorr[MPSIZE];
	double Electron_EnergySCUnCorr[MPSIZE];
	double Electron_etaSCUnCorr[MPSIZE];
	double Electron_phiSCUnCorr[MPSIZE];
	double Electron_etSCUnCorr[MPSIZE];

	// Pat Muon
	//pf isolations
	double Muon_PfChargedHadronIsoR05[MPSIZE];
	double Muon_PfNeutralHadronIsoR05[MPSIZE];
	double Muon_PfGammaIsoR05[MPSIZE];
	double Muon_PfChargedHadronIsoR04[MPSIZE];
	double Muon_PfNeutralHadronIsoR04[MPSIZE];
	double Muon_PfGammaIsoR04[MPSIZE];
	double Muon_PFSumPUIsoR04[MPSIZE];
	double Muon_PfChargedHadronIsoR03[MPSIZE];
	double Muon_PfNeutralHadronIsoR03[MPSIZE];
	double Muon_PfGammaIsoR03[MPSIZE];
	double Muon_PFSumPUIsoR03[MPSIZE];

	int isPFmuon[MPSIZE];
	int isGLBmuon[MPSIZE];
	int isTRKmuon[MPSIZE];
	int isSTAmuon[MPSIZE];

	int Muon_muonType[MPSIZE];
	int Muon_nTrig[MPSIZE];
	int Muon_triggerObjectType[MPSIZE];
	int Muon_filterName[MPSIZE];
	double Muon_dB[MPSIZE];
	double Muon_phi[MPSIZE];
	double Muon_eta[MPSIZE];
	double Muon_pT[MPSIZE];
	double Muon_cktpT[MPSIZE];
	double Muon_cktPx[MPSIZE];
	double Muon_cktPy[MPSIZE];
	double Muon_cktPz[MPSIZE];
	double Muon_cktpTError[MPSIZE];
	double Muon_Px[MPSIZE];
	double Muon_Py[MPSIZE];
	double Muon_Pz[MPSIZE];
	double Muon_sumtrkpt[MPSIZE];
	double Muon_trkiso[MPSIZE];
	double Muon_hcaliso[MPSIZE];
	double Muon_ecaliso[MPSIZE];
	double Muon_trkisoR05[MPSIZE];
	double Muon_hcalisoR05[MPSIZE];
	double Muon_ecalisoR05[MPSIZE];
	int Muon_charge[MPSIZE];
	int Muon_nChambers[MPSIZE];
	int Muon_nMatches[MPSIZE];
	int Muon_nMatchesRPCLayers[MPSIZE];
	int Muon_stationMask[MPSIZE];
	int Muon_nSegments[MPSIZE];
	double Muon_chi2dof[MPSIZE];
	int Muon_nhits[MPSIZE];
	int Muon_trackerHits[MPSIZE];
	int Muon_pixelHits[MPSIZE];
	int Muon_muonHits[MPSIZE];
	int Muon_trackerLayers[MPSIZE];
	int Muon_trackerHitsGLB[MPSIZE];
	int Muon_trackerLayersGLB[MPSIZE];
	int Muon_pixelHitsGLB[MPSIZE];
	double Muon_qoverp[MPSIZE];
	double Muon_theta[MPSIZE];
	double Muon_lambda[MPSIZE];
	double Muon_dxy[MPSIZE];
	double Muon_d0[MPSIZE];
	double Muon_dsz[MPSIZE];
	double Muon_dz[MPSIZE];
	double Muon_dxyBS[MPSIZE];
	double Muon_dzBS[MPSIZE];
	double Muon_dszBS[MPSIZE];
	double Muon_dxyVTX[MPSIZE];
	double Muon_dzVTX[MPSIZE];
	double Muon_dszVTX[MPSIZE];
	double Muon_dxycktVTX[MPSIZE];
	double Muon_dzcktVTX[MPSIZE];
	double Muon_dszcktVTX[MPSIZE];
	double Muon_vx[MPSIZE];
	double Muon_vy[MPSIZE];
	double Muon_vz[MPSIZE];
	std::vector<double> CosAngle;
	std::vector<double> vtxTrkChi2;
	std::vector<double> vtxTrkProb;
	std::vector<double> vtxTrkNdof;
	std::vector<double> vtxTrkCkt1Pt;
	std::vector<double> vtxTrkCkt2Pt;

	std::vector<double> CosAngle_TuneP;
	std::vector<double> vtxTrk1Pt_TuneP;
	std::vector<double> vtxTrk2Pt_TuneP;
	std::vector<double> vtxTrkChi2_TuneP;
	std::vector<double> vtxTrkNdof_TuneP;
	std::vector<double> vtxTrkProb_TuneP;

	//Various track informations
	//MuonBestTrack
	double Muon_Best_pT[MPSIZE];
	double Muon_Best_pTError[MPSIZE];
	double Muon_Best_Px[MPSIZE];
	double Muon_Best_Py[MPSIZE];
	double Muon_Best_Pz[MPSIZE];
	double Muon_Best_eta[MPSIZE];
	double Muon_Best_phi[MPSIZE];
	//Inner Track
	double Muon_Inner_pT[MPSIZE];
	double Muon_Inner_pTError[MPSIZE];
	double Muon_Inner_Px[MPSIZE];
	double Muon_Inner_Py[MPSIZE];
	double Muon_Inner_Pz[MPSIZE];
	double Muon_Inner_eta[MPSIZE];
	double Muon_Inner_phi[MPSIZE];
	//Outer Track
	double Muon_Outer_pT[MPSIZE];
	double Muon_Outer_pTError[MPSIZE];
	double Muon_Outer_Px[MPSIZE];
	double Muon_Outer_Py[MPSIZE];
	double Muon_Outer_Pz[MPSIZE];
	double Muon_Outer_eta[MPSIZE];
	double Muon_Outer_phi[MPSIZE];
	//Global Track
	double Muon_GLB_pT[MPSIZE];
	double Muon_GLB_pTError[MPSIZE];
	double Muon_GLB_Px[MPSIZE];
	double Muon_GLB_Py[MPSIZE];
	double Muon_GLB_Pz[MPSIZE];
	double Muon_GLB_eta[MPSIZE];
	double Muon_GLB_phi[MPSIZE];

	//tuneP MuonBestTrack
	double Muon_TuneP_pT[MPSIZE];
	double Muon_TuneP_pTError[MPSIZE];
	double Muon_TuneP_Px[MPSIZE];
	double Muon_TuneP_Py[MPSIZE];
	double Muon_TuneP_Pz[MPSIZE];
	double Muon_TuneP_eta[MPSIZE];
	double Muon_TuneP_phi[MPSIZE];

	//Medium ID
	double Muon_segmentCompatibility[MPSIZE];
	double Muon_chi2LocalPosition[MPSIZE];
	double Muon_trkKink[MPSIZE];
	double Muon_Inner_validFraction[MPSIZE];

	//Muon IDs
	bool Muon_passLooseID[MPSIZE];
	bool Muon_passMediumID[MPSIZE];
	bool Muon_passMediumID_from_var[MPSIZE];
	bool Muon_passTightID[MPSIZE];
	bool Muon_passHighPtID[MPSIZE];


	// LHE
	int nLHEParticle;
	double LHEParticle_Px[MPSIZE];
	double LHEParticle_Py[MPSIZE];
	double LHEParticle_Pz[MPSIZE];
	double LHEParticle_E[MPSIZE];
	int LHEParticle_ID[MPSIZE];
	int LHEParticle_status[MPSIZE];

	// GEN
	int GENnPair;
	double GENLepton_phi[MPSIZE];
	double GENLepton_eta[MPSIZE];
	double GENLepton_pT[MPSIZE];
	double GENLepton_Px[MPSIZE];
	double GENLepton_Py[MPSIZE];
	double GENLepton_Pz[MPSIZE];
	double GENLepton_E[MPSIZE];
	double GENLepton_mother[MPSIZE];
	double GENLepton_mother_pT[MPSIZE];
	int GENLepton_charge[MPSIZE];
	int GENLepton_status[MPSIZE];
	int GENLepton_ID[MPSIZE];
	int GENLepton_isPrompt[MPSIZE];
	int GENLepton_isPromptFinalState[MPSIZE];
	int GENLepton_isTauDecayProduct[MPSIZE];
	int GENLepton_isPromptTauDecayProduct[MPSIZE];
	int GENLepton_isDirectPromptTauDecayProductFinalState[MPSIZE];
	int GENLepton_isHardProcess[MPSIZE];
	int GENLepton_isLastCopy[MPSIZE];
	int GENLepton_isLastCopyBeforeFSR[MPSIZE];
	int GENLepton_isPromptDecayed[MPSIZE];
	int GENLepton_isDecayedLeptonHadron[MPSIZE];
	int GENLepton_fromHardProcessBeforeFSR[MPSIZE];
	int GENLepton_fromHardProcessDecayed[MPSIZE];
	int GENLepton_fromHardProcessFinalState[MPSIZE];
	int GENLepton_isMostlyLikePythia6Status3[MPSIZE];
	double GENEvt_weight;
	double GENEvt_QScale;
	double GENEvt_x1;
	double GENEvt_x2;
	double GENEvt_alphaQCD;
	double GENEvt_alphaQED;

	int nGenOthers;
	double GenOthers_phi[MPSIZE];
	double GenOthers_eta[MPSIZE];
	double GenOthers_pT[MPSIZE];
	double GenOthers_Px[MPSIZE];
	double GenOthers_Py[MPSIZE];
	double GenOthers_Pz[MPSIZE];
	double GenOthers_E[MPSIZE];
	double GenOthers_mother[MPSIZE];
	int GenOthers_charge[MPSIZE];
	int GenOthers_status[MPSIZE];
	int GenOthers_ID[MPSIZE];
	int GenOthers_isPrompt[MPSIZE];
	int GenOthers_isPromptFinalState[MPSIZE];
	int GenOthers_isTauDecayProduct[MPSIZE];
	int GenOthers_isPromptTauDecayProduct[MPSIZE];
	int GenOthers_isDirectPromptTauDecayProductFinalState[MPSIZE];
	int GenOthers_isHardProcess[MPSIZE];
	int GenOthers_isLastCopy[MPSIZE];
	int GenOthers_isLastCopyBeforeFSR[MPSIZE];
	int GenOthers_isPromptDecayed[MPSIZE];
	int GenOthers_isDecayedLeptonHadron[MPSIZE];
	int GenOthers_fromHardProcessBeforeFSR[MPSIZE];
	int GenOthers_fromHardProcessDecayed[MPSIZE];
	int GenOthers_fromHardProcessFinalState[MPSIZE];
	int GenOthers_isMostlyLikePythia6Status3[MPSIZE];

	// -- Photon information -- //
	int nPhotons;
	double Photon_pT[MPSIZE];
	double Photon_eta[MPSIZE];
	double Photon_phi[MPSIZE];
	double Photon_etaSC[MPSIZE];
	double Photon_phiSC[MPSIZE];
	double Photon_HoverE[MPSIZE];
	int Photon_hasPixelSeed[MPSIZE];
	double Photon_Full5x5_SigmaIEtaIEta[MPSIZE];
	double Photon_ChIso[MPSIZE];
	double Photon_NhIso[MPSIZE];
	double Photon_PhIso[MPSIZE];
	double Photon_ChIsoWithEA[MPSIZE];
	double Photon_NhIsoWithEA[MPSIZE];
	double Photon_PhIsoWithEA[MPSIZE];
	// Effective area constants for all isolation types
	// EffectiveAreas effAreaChHadrons_;
	// EffectiveAreas effAreaNeuHadrons_;
	// EffectiveAreas effAreaPhotons_;
	bool Photon_passMediumID[MPSIZE];

	int NTT;
	double TTrack_dxy[MPSIZE];
	double TTrack_dxyErr[MPSIZE];
	double TTrack_d0[MPSIZE];
	double TTrack_d0Err[MPSIZE];
	double TTrack_dsz[MPSIZE];
	double TTrack_dszErr[MPSIZE];
	double TTrack_dz[MPSIZE];
	double TTrack_dzErr[MPSIZE];
	double TTrack_dxyBS[MPSIZE];
	double TTrack_dszBS[MPSIZE];
	double TTrack_dzBS[MPSIZE];
	double TTrack_pT[MPSIZE];
	double TTrack_Px[MPSIZE];
	double TTrack_Py[MPSIZE];
	double TTrack_Pz[MPSIZE];
	double TTrack_eta[MPSIZE];
	double TTrack_phi[MPSIZE];
	double TTrack_charge[MPSIZE];

	// -- MET -- //
	double pfMET_pT;
	double pfMET_phi; 
	double pfMET_Px; 
	double pfMET_Py; 
	double pfMET_SumEt;

	double pfMET_Type1_pT;
	double pfMET_Type1_phi; 
	double pfMET_Type1_Px; 
	double pfMET_Type1_Py; 
	double pfMET_Type1_SumEt;

	double pfMET_Type1_PhiCor_pT;
	double pfMET_Type1_PhiCor_phi; 
	double pfMET_Type1_PhiCor_Px; 
	double pfMET_Type1_PhiCor_Py; 
	double pfMET_Type1_PhiCor_SumEt;

};
#endif
