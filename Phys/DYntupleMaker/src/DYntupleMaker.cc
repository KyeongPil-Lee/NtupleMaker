//-------------------------------------------------
//
//   Class: DYntupleMaker
//
//   Description: Ntuple maker for DY analysis
//
//
//   Author:
//   H.D.Yoo           Purdue University
//   K.P.Lee           Seoul National University
//
//--------------------------------------------------

#include "Phys/DYntupleMaker/interface/DYntupleMaker.h"

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
#include "FWCore/ParameterSet/interface/ParameterSetfwd.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RegexMatch.h"


/////////////////////
// -- For Muons -- //
/////////////////////
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/L1GlobalMuonTrigger/interface/L1MuGMTCand.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"

/////////////////////////
// -- For Electrons -- //
/////////////////////////
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

///////////////////
// -- For MET -- //
///////////////////
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"

////////////////////
// -- For Jets -- //
////////////////////
#include "DataFormats/PatCandidates/interface/Jet.h"

////////////////////////////
// -- For GenParticles -- //
////////////////////////////
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

//////////////////////////
// -- Track & Vertex -- //
//////////////////////////
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/PatCandidates/interface/Vertexing.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/InvariantMassFromVertex.h"
#include "PhysicsTools/RecoUtils/interface/CandCommonVertexFitter.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

////////////////////
// -- Triggers -- //
////////////////////
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
//#include "DataFormats/PatCandidates/interface/TriggerPrimitive.h"

////////////////
// -- Else -- //
////////////////
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/GeometrySurface/interface/GloballyPositioned.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositVetos.h"
#include "PhysicsTools/IsolationAlgos/interface/IsoDepositVetoFactory.h"
#include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h"
#include "DataFormats/PatCandidates/interface/Isolation.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "MagneticField/Engine/interface/MagneticField.h"

////////////////
// -- ROOT -- //
////////////////
#include <TTree.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TFrame.h>
#include <TMath.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <boost/foreach.hpp>

//
// class decleration
//

using namespace std;
using namespace reco;
using namespace edm;
using namespace pat;
using namespace isodeposit;

// -- Constructor -- //
DYntupleMaker::DYntupleMaker(const edm::ParameterSet& iConfig):
// -- object tokens -- //
MuonToken						( consumes< std::vector<pat::Muon> > 				(iConfig.getUntrackedParameter<edm::InputTag>("Muon")) ),
ElectronToken					( consumes< edm::View<reco::GsfElectron> >			(iConfig.getUntrackedParameter<edm::InputTag>("Electron")) ),
PhotonToken 					( consumes< edm::View<reco::Photon> >				(iConfig.getUntrackedParameter<edm::InputTag>("Photon")) ),
JetToken 						( consumes< std::vector<pat::Jet> >					(iConfig.getUntrackedParameter<edm::InputTag>("Jet")) ),
MetToken 						( consumes< std::vector<pat::MET> >					(iConfig.getUntrackedParameter<edm::InputTag>("MET")) ),
LHEEventProductToken			( consumes< LHEEventProduct >  						(iConfig.getUntrackedParameter<edm::InputTag>("LHEEventProduct")) ),
GenParticleToken 				( consumes< std::vector<reco::GenParticle> >		(iConfig.getUntrackedParameter<edm::InputTag>("GenParticle")) ),
// -- Electron tokens -- //
RhoToken 						( consumes< double >								(iConfig.getUntrackedParameter<edm::InputTag>("rho")) ),
eleVetoIdMapToken 				( consumes< edm::ValueMap<bool> >					(iConfig.getUntrackedParameter<edm::InputTag>("eleVetoIdMap")) ),
eleLooseIdMapToken 				( consumes< edm::ValueMap<bool> >					(iConfig.getUntrackedParameter<edm::InputTag>("eleLooseIdMap")) ),
eleMediumIdMapToken 			( consumes< edm::ValueMap<bool> >					(iConfig.getUntrackedParameter<edm::InputTag>("eleMediumIdMap")) ),
eleTightIdMapToken 				( consumes< edm::ValueMap<bool> >					(iConfig.getUntrackedParameter<edm::InputTag>("eleTightIdMap")) ),
eleMVAIdWP80MapToken 			( consumes< edm::ValueMap<bool> >					(iConfig.getUntrackedParameter<edm::InputTag>("eleMVAIdWP80Map")) ),
eleMVAIdWP90MapToken 			( consumes< edm::ValueMap<bool> >					(iConfig.getUntrackedParameter<edm::InputTag>("eleMVAIdWP90Map")) ),
eleHEEPIdMapToken 				( consumes< edm::ValueMap<bool> >					(iConfig.getUntrackedParameter<edm::InputTag>("eleHEEPIdMap")) ),
ConversionsToken 				( consumes< std::vector<reco::Conversion> > 		(iConfig.getUntrackedParameter<edm::InputTag>("conversionsInputTag")) ),
GsfTrackToken					( consumes< std::vector< reco::GsfTrack > > 		(iConfig.getUntrackedParameter<edm::InputTag>("GsfTrack")) ),
// -- Photon tokens -- //
full5x5SigmaIEtaIEtaMapToken	( consumes< edm::ValueMap<float> > 					(iConfig.getUntrackedParameter<edm::InputTag>("full5x5SigmaIEtaIEtaMap")) ),
phoChargedIsolationToken 		( consumes< edm::ValueMap<float> > 					(iConfig.getUntrackedParameter<edm::InputTag>("phoChargedIsolation")) ),
phoNeutralHadronIsolationToken 	( consumes< edm::ValueMap<float> > 					(iConfig.getUntrackedParameter<edm::InputTag>("phoNeutralHadronIsolation")) ),
phoPhotonIsolationToken 		( consumes< edm::ValueMap<float> > 					(iConfig.getUntrackedParameter<edm::InputTag>("phoPhotonIsolation")) ),
// -- Trigger Token -- //
TriggerToken 					( consumes< edm::TriggerResults >  					(iConfig.getUntrackedParameter<edm::InputTag>("TriggerResults")) ),
TriggerTokenPAT 				( consumes< edm::TriggerResults >  					(iConfig.getUntrackedParameter<edm::InputTag>("TriggerResultsPAT")) ),
TriggerObjectToken 				( consumes< std::vector<pat::TriggerObjectStandAlone> >  	(iConfig.getUntrackedParameter<edm::InputTag>("TriggerObject")) ),
// -- Else -- //
GenEventInfoToken 				( consumes< GenEventInfoProduct >  					(iConfig.getUntrackedParameter<edm::InputTag>("GenEventInfo")) ),
BeamSpotToken					( consumes< reco::BeamSpot > 						(iConfig.getUntrackedParameter<edm::InputTag>("BeamSpot")) ),
PrimaryVertexToken 				( consumes< reco::VertexCollection > 				(iConfig.getUntrackedParameter<edm::InputTag>("PrimaryVertex")) ),
TrackToken 						( consumes< edm::View<reco::Track> >  				(iConfig.getUntrackedParameter<edm::InputTag>("Track")) ),
PileUpInfoToken 				( consumes< std::vector< PileupSummaryInfo > >  	(iConfig.getUntrackedParameter<edm::InputTag>("PileUpInfo")) )
{
	nEvt = 0;

	isMC                              = iConfig.getUntrackedParameter<bool>("isMC");
	processName                       = iConfig.getUntrackedParameter<string>("processName", "HLT");
	theDebugLevel                     = iConfig.getUntrackedParameter<int>("DebugLevel", 0);

	// -- Consumes -- //



	// -- Object Labels -- //
	// theMuonLabel                      = iConfig.getUntrackedParameter<edm::InputTag>("Muon", edm::InputTag("selectedPatMuonsPFlow"));
	// theElectronLabel                  = iConfig.getUntrackedParameter<edm::InputTag>("Electron", edm::InputTag("gedGsfElectrons"));
	// thePhotonLabel                    = iConfig.getUntrackedParameter<edm::InputTag>("Photon", edm::InputTag("gedPhotons"));
	// theJetLabel                       = iConfig.getUntrackedParameter<edm::InputTag>("Jet", edm::InputTag("selectedPatJets"));
	// theMETLabel                       = iConfig.getUntrackedParameter<edm::InputTag>("MET", edm::InputTag("patMETs"));
	// theGenParticleLabel				  = iConfig.getUntrackedParameter<edm::InputTag>("GenParticle", edm::InputTag("genParticles"));
	
	// -- for Electrons -- //
	// rho                               = iConfig.getUntrackedParameter<edm::InputTag>("rho", edm::InputTag("fixedGridRhoFastjetAll"));
	// eleVetoIdMap                      = iConfig.getUntrackedParameter<edm::InputTag>("eleVetoIdMap", edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"));
	// eleLooseIdMap                     = iConfig.getUntrackedParameter<edm::InputTag>("eleLooseIdMap", edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"));
	// eleMediumIdMap                    = iConfig.getUntrackedParameter<edm::InputTag>("eleMediumIdMap", edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"));
	// eleTightIdMap                     = iConfig.getUntrackedParameter<edm::InputTag>("eleTightIdMap", edm::InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"));
	// theConversionsInputTag            = iConfig.getUntrackedParameter<edm::InputTag>("conversionsInputTag", edm::InputTag("allConversions"));

	// -- for Photons -- //
	// full5x5SigmaIEtaIEtaMapLabel      = iConfig.getUntrackedParameter<edm::InputTag>("full5x5SigmaIEtaIEtaMap", edm::InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta")); 
	// phoChargedIsolationLabel          = iConfig.getUntrackedParameter<edm::InputTag>("phoChargedIsolation", edm::InputTag("photonIDValueMapProducer:phoChargedIsolation"));
	// phoNeutralHadronIsolationLabel    = iConfig.getUntrackedParameter<edm::InputTag>("phoNeutralHadronIsolation", edm::InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"));
	// phoPhotonIsolationLabel           = iConfig.getUntrackedParameter<edm::InputTag>("phoPhotonIsolation", edm::InputTag("photonIDValueMapProducer:phoPhotonIsolation"));
	effAreaChHadronsFile              = iConfig.getUntrackedParameter<edm::FileInPath>( "effAreaChHadFile", edm::FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfChargedHadrons_V2.txt") );
	effAreaNeuHadronsFile             = iConfig.getUntrackedParameter<edm::FileInPath>( "effAreaNeuHadFile", edm::FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfNeutralHadrons_V2.txt") );
	effAreaPhotonsFile                = iConfig.getUntrackedParameter<edm::FileInPath>( "effAreaPhoFile", edm::FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfPhotons_V2.txt") );

 	// -- for Jet -- //
	// theBDiscriminant_alg1             = iConfig.getUntrackedParameter<double>("BDiscriminant_tcheff", 0.7);
	// theBDiscriminant_alg2             = iConfig.getUntrackedParameter<double>("BDiscriminant_tchpur", 0.7);
	// theBDiscriminant_alg3             = iConfig.getUntrackedParameter<double>("BDiscriminant_ssv", 2.05);

	// -- for MET -- //
	// pfMetCollection_                  = iConfig.getUntrackedParameter<edm::InputTag>("pfMET", edm::InputTag("pfMet"));

	// -- Else -- //
	// theBeamSpot                       = iConfig.getUntrackedParameter<edm::InputTag>("BeamSpot", edm::InputTag("offlineBeamSpot"));
	// thePrimaryVertexLabel			  = iConfig.getUntrackedParameter<edm::InputTag>("PrimaryVertex", edm::InputTag("offlinePrimaryVerticesWithBS"));
	// theTrackLabel                     = iConfig.getUntrackedParameter<edm::InputTag>("Track", edm::InputTag("generalTracks"));
	// thePileUpInfoLabel				  = iConfig.getUntrackedParameter<edm::InputTag>("PileUpInfo", edm::InputTag("addPileupInfo"));
	

	// -- Store Flags -- //
	theStorePriVtxFlag                = iConfig.getUntrackedParameter<bool>("StorePriVtxFlag", true);
	theStoreJetFlag                   = iConfig.getUntrackedParameter<bool>("StoreJetFlag", false);
	theStoreMETFlag                   = iConfig.getUntrackedParameter<bool>("StoreMETFlag", false);
	theStoreHLTReportFlag             = iConfig.getUntrackedParameter<bool>("StoreHLTReportFlag", true);
	
	theStoreMuonFlag              	  = iConfig.getUntrackedParameter<bool>("StoreMuonFlag", true);
	theStoreElectronFlag              = iConfig.getUntrackedParameter<bool>("StoreElectronFlag", true);
	theStoreLHEFlag                   = iConfig.getUntrackedParameter<bool>("StoreLHEFlag", false);
	theStoreGENFlag                   = iConfig.getUntrackedParameter<bool>("StoreGENFlag", true);
	theStoreGenOthersFlag             = iConfig.getUntrackedParameter<bool>("StoreGenOthersFlag", false);
	theStoreTTFlag                    = iConfig.getUntrackedParameter<bool>("StoreTTFlag", false);
	theStorePhotonFlag                = iConfig.getUntrackedParameter<bool>("StorePhotonFlag", false);

	// -- Filters -- //
	theApplyFilter                    = iConfig.getUntrackedParameter<bool>("ApplyFilter", false);
	theFilterType                     = iConfig.getUntrackedParameter<int>("FilterType", 0); // 0: dilepton (default), 1: single muon (for fake rate), 2: single photon (fake rate)

	// if( isMC )
	// {
	// 	PileUpRD_ = iConfig.getParameter< std::vector<double> >("PileUpRD");
	// 	PileUpRDMuonPhys_ = iConfig.getParameter< std::vector<double> >("PileUpRDMuonPhys");
	// 	PileUpMC_ = iConfig.getParameter< std::vector<double> >("PileUpMC");
	// }

}

DYntupleMaker::~DYntupleMaker() { }

//
// member functions
//

// ------------ method called to for each event  ------------ //
void DYntupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theTTBuilder);
	 
	///////////////////////////////////////////
	// -- initialize for ntuple variables -- //
	///////////////////////////////////////////
	nMuon = -1; 
	// MET_sumEt = -999;
	// MET_pt = -999;
	// MET_px = -999;
	// MET_py = -999;
	// MET_phi = -999;
	// pfMET_sumEt = -999;
	// pfMET_pt = -999;
	// pfMET_px = -999;
	// pfMET_py = -999;
	// pfMET_phi = -999;

	Flag_duplicateMuons = false;
	Flag_badMuons = false;
	Flag_noBadMuons = false;

	nVertices = -1;
	Nmuons = Njets = Nbtagged = NbtaggedCloseMuon = -1;
	Nelectrons = -1;
	PVtrackSize = -1;
	PVchi2 = -1;
	PVndof = -1;
	PVnormalizedChi2 = -1;
	PVx = -1000;
	PVy = -1000;
	PVz = -1000;
	PVprob = -1;

	nLHEParticle = -1;
	GENnPair = -1;

	// -- PF iso deposits -- // 
	sumEt = 0;
	photonEt = 0;
	chargedHadronEt = 0;
	neutralHadronEt = 0;

	// -- trigger object -- //
	_HLT_ntrig = -1;
	_HLT_trigName.clear();
	_HLT_trigPS.clear();

	// -- PU reweight -- //
	PUweight = -1;
	pileUpReweightIn = pileUpReweight = 1.0;
	pileUpReweightPlus = pileUpReweightMinus = 0.0;
	pileUpReweightInMuonPhys = pileUpReweightMuonPhys = 1.0;
	pileUpReweightPlusMuonPhys = pileUpReweightMinusMuonPhys = 0.0;


	CosAngle.clear();
	vtxTrkCkt1Pt.clear();
	vtxTrkCkt2Pt.clear();
	vtxTrkProb.clear();
	vtxTrkNdof.clear();
	vtxTrkChi2.clear();
	vtxTrkDiE1Pt.clear();
	vtxTrkDiE2Pt.clear();
	vtxTrkDiEProb.clear();
	vtxTrkDiENdof.clear();
	vtxTrkDiEChi2.clear();
	vtxTrkEMu1Pt.clear();
	vtxTrkEMu2Pt.clear();
	vtxTrkEMuProb.clear();
	vtxTrkEMuNdof.clear();
	vtxTrkEMuChi2.clear();

	CosAngle_TuneP.clear();
	vtxTrk1Pt_TuneP.clear();
	vtxTrk2Pt_TuneP.clear();
	vtxTrkChi2_TuneP.clear();
	vtxTrkNdof_TuneP.clear();
	vtxTrkProb_TuneP.clear();

	for( int i = 0; i < MPSIZE; i++ )
	{
		// -- Trigger -- //
		_HLT_trigType[i] = -1;
		_HLT_trigFired[i] = -1;
		_HLT_trigPt[i] = _HLT_trigEta[i] = _HLT_trigPhi[i] = -100;

		// -- Jet -- //
		JETbDiscriminant[i] = JETcharge[i] = JETpt[i] = JETeta[i] = JETphi[i] = -100;
		JETflavour[i] = JETntracks[i] = -100;

		Jet_pT[i] = -100; 
		Jet_eta[i] = -100; 
		Jet_phi[i] = -100; 
		Jet_Charge[i] = -100; 
		Jet_Flavor[i] = -100; 

		Jet_bTag[i] = -100; 
		Jet_CHfrac[i] = -100; 
		Jet_NHfrac[i] = -100; 
		Jet_NHEMfrac[i] = -100; 
		Jet_CHEMfrac[i] = -100; 
		Jet_CHmulti[i] = -100; 
		Jet_NHmulti[i] = -100;

		// electron
		Electron_Energy[i] = Electron_et[i] = Electron_pT[i] = Electron_eta[i] = Electron_phi[i] = -100;
		Electron_etCorr[i] = Electron_caloEnergy[i] = -100;
		Electron_Px[i] = Electron_Py[i] = Electron_Pz[i] = -9999;
		Electron_ECalib[i] = Electron_etCalib[i] = Electron_pTCalib[i] = -100;
		Electron_etCorrCalib[i] = Electron_caloEnergyCalib[i] = -100;
		Electron_PxCalib[i] = Electron_PyCalib[i] = Electron_PzCalib[i] = -9999;
		Electron_charge[i] = -100;
		Electron_gsfpT[i] = Electron_gsfEta[i] = Electron_gsfPhi[i] = -100;
		Electron_gsfPx[i] = Electron_gsfPy[i] = Electron_gsfPz[i] = -9999;
		Electron_gsfpTCalib[i] = Electron_gsfPxCalib[i] = Electron_gsfPyCalib[i] = Electron_gsfPzCalib[i] = -9999;
		Electron_gsfCharge[i] = -100;
		Electron_etaSC[i] = -9999;
		Electron_phiSC[i] = -9999;
		Electron_etaWidth[i] = -9999;
		Electron_phiWidth[i] = -9999;
		Electron_dEtaIn[i] = -9999;
		Electron_dPhiIn[i] = -9999;
		Electron_sigmaIEtaIEta[i] = -9999;
		Electron_HoverE[i] = -9999;
		Electron_fbrem[i] = -9999;
		Electron_eOverP[i] = -9999;
		Electron_energyEC[i] = -9999;
		Electron_Pnorm[i] = -9999;
		Electron_InvEminusInvP[i] = -9999;
		Electron_dxyVTX[i] = -9999;
		Electron_dzVTX[i] = -9999;
		Electron_dxy[i] = -9999;
		Electron_dz[i] = -9999;
		Electron_dxyBS[i] = -9999;
		Electron_dzBS[i] = -9999;
		Electron_AEff03[i] = -9999;
		Electron_chIso03[i] = -9999;
		Electron_chIso04[i] = -9999;
		Electron_nhIso03[i] = -9999;
		Electron_nhIso04[i] = -9999;
		Electron_phIso03[i] = -9999;
		Electron_phIso04[i] = -9999;
		Electron_pcIso03[i] = -9999;
		Electron_pcIso04[i] = -9999;
		Electron_relIsoCom03[i] = -9999;
		Electron_relIsoCom04[i] = -9999;
		Electron_relIsoBeta03[i] = -9999;
		Electron_relIsoBeta04[i] = -9999;
		Electron_relIsoRho03[i] = -9999;
		Electron_hasConversion[i] = false;
		Electron_mHits[i] = -9999;

		Electron_crack[i] = -1;
		Electron_ecalDriven[i] = -1;
		Electron_isoEMHADDepth1[i] = -999;
		Electron_25over55[i] = -999;
		Electron_15over55[i] = -999;
		Electron_isoHADDepth2[i] = -999;
		Electron_isoPtTrks[i] = -999;
		Electron_modIsoEMHADDepth1[i] = -999;
		Electron_modIsoPtTrks[i] = -999;
		Electron_modIsoEMHADDepth1Orig[i] = -999;
		Electron_modIsoPtTrksOrig[i] = -999;

		Electron_ambGsf0Pt[i] = Electron_ambGsf0Eta[i] = Electron_ambGsf0Phi[i] = Electron_ambGsf0Charge[i] = -999;
		Electron_ambGsf1Pt[i] = Electron_ambGsf1Eta[i] = Electron_ambGsf1Phi[i] = Electron_ambGsf1Charge[i] = -999;
		Electron_ambGsf2Pt[i] = Electron_ambGsf2Eta[i] = Electron_ambGsf2Phi[i] = Electron_ambGsf2Charge[i] = -999;
		Electron_ambGsf3Pt[i] = Electron_ambGsf3Eta[i] = Electron_ambGsf3Phi[i] = Electron_ambGsf3Charge[i] = -999;

		Electron_r9[i] = -999;
		Electron_EnergySC[i] = -999;
		Electron_preEnergySC[i] = -999;
		Electron_rawEnergySC[i] = -999;
		Electron_etSC[i] = -999;
		Electron_E15[i] = -999;
		Electron_E25[i] = -999;
		Electron_E55[i] = -999;
		Electron_ChIso03FromPU[i] = -999;
		Electron_RelPFIso_dBeta[i] = -999;
		Electron_RelPFIso_Rho[i] = -999;
		Electron_passConvVeto[i] = 0;
		Electron_passVetoID[i] = 0;
		Electron_passLooseID[i] = 0;
		Electron_passMediumID[i] = 0;
		Electron_passTightID[i] = 0;
		Electron_passMVAID_WP80[i] = 0;
		Electron_passMVAID_WP90[i] = 0;
		Electron_passHEEPID[i] = 0;

		// -- PF Isolation -- //
		Muon_PfChargedHadronIsoR05[i] = -1;
		Muon_PfNeutralHadronIsoR05[i] = -1;
		Muon_PfGammaIsoR05[i] = -1;
		Muon_PfChargedHadronIsoR04[i] = -1;
		Muon_PfNeutralHadronIsoR04[i] = -1;
		Muon_PfGammaIsoR04[i] = -1;
		Muon_PFSumPUIsoR04[i] = -1;
		Muon_PfChargedHadronIsoR03[i] = -1;
		Muon_PfNeutralHadronIsoR03[i] = -1;
		Muon_PfGammaIsoR03[i] = -1;
		Muon_PFSumPUIsoR03[i] = -1;

		// muon type
		Muon_muonType[i] = -1;
		isPFmuon[i] = 0;
		isGLBmuon[i] = 0;
		isTRKmuon[i] = 0;
		isSTAmuon[i] = 0;

		// trigger
		Muon_nTrig[i] = Muon_triggerObjectType[i] = Muon_filterName[i] = -1;

		// Muon kinematics
		Muon_phi[i] = Muon_eta[i] = Muon_cktpT[i] = Muon_pT[i] = -100;
		Muon_dB[i] = -999;
		Muon_cktPx[i] = Muon_cktPy[i] = Muon_cktPz[i] = -999;
		Muon_cktpTError[i] = -999;
		Muon_Px[i] = Muon_Py[i] = Muon_Pz[i] = -100;
		Muon_trkiso[i] = Muon_hcaliso[i] = Muon_ecaliso[i] = Muon_chi2dof[i] = -100;
		Muon_trkisoR05[i] = Muon_hcalisoR05[i] = Muon_ecalisoR05[i] = -100;
		Muon_nChambers[i] = Muon_nMatches[i] = Muon_nMatchesRPCLayers[i] = Muon_stationMask[i] = Muon_nSegments[i] =  -1;
		Muon_charge[i] = Muon_nhits[i] = -100;
		Muon_trackerHits[i] = Muon_pixelHits[i] = Muon_muonHits[i] = -1;
		Muon_trackerLayers[i] = -1;
		Muon_trackerHitsGLB[i] = -1;
		Muon_pixelHitsGLB[i] = -1;
		Muon_trackerLayersGLB[i] = -1;

		Muon_qoverp[i] = Muon_theta[i] = Muon_lambda[i] = -100;
		Muon_dxy[i] = Muon_d0[i] = Muon_dsz[i] = Muon_dz[i] = -100;
		Muon_vx[i] = Muon_vy[i] = Muon_vz[i] = -100;
		Muon_dxyBS[i] = Muon_dszBS[i] = Muon_dzBS[i] = -100;
		Muon_dxyVTX[i] = Muon_dszVTX[i] = Muon_dzVTX[i] = -100;
		Muon_dxycktVTX[i] = Muon_dszcktVTX[i] = Muon_dzcktVTX[i] = -100;

		//Various track informations
		//MuonBestTrack
		Muon_Best_pT[i] = -9999;
		Muon_Best_pTError[i] = -9999;
		Muon_Best_Px[i] = -9999;
		Muon_Best_Py[i] = -9999;
		Muon_Best_Pz[i] = -9999;
		Muon_Best_eta[i] = -9999;
		Muon_Best_phi[i] = -9999;
		//Inner Track
		Muon_Inner_pT[i] = -9999;
		Muon_Inner_pTError[i] = -9999;
		Muon_Inner_Px[i] = -9999;
		Muon_Inner_Py[i] = -9999;
		Muon_Inner_Pz[i] = -9999;
		Muon_Inner_eta[i] = -9999;
		Muon_Inner_phi[i] = -9999;
		//Outer Track
		Muon_Outer_pT[i] = -9999;
		Muon_Outer_pTError[i] = -9999;
		Muon_Outer_Px[i] = -9999;
		Muon_Outer_Py[i] = -9999;
		Muon_Outer_Pz[i] = -9999;
		Muon_Outer_eta[i] = -9999;
		Muon_Outer_phi[i] = -9999;
		//Global Track
		Muon_GLB_pT[i] = -9999;
		Muon_GLB_pTError[i] = -9999;
		Muon_GLB_Px[i] = -9999;
		Muon_GLB_Py[i] = -9999;
		Muon_GLB_Pz[i] = -9999;
		Muon_GLB_eta[i] = -9999;
		Muon_GLB_phi[i] = -9999;

		//tuneP MuonBestTrack
		Muon_TuneP_pT[i] = -9999;
		Muon_TuneP_pTError[i] = -9999;
		Muon_TuneP_Px[i] = -9999;
		Muon_TuneP_Py[i] = -9999;
		Muon_TuneP_Pz[i] = -9999;
		Muon_TuneP_eta[i] = -9999;
		Muon_TuneP_phi[i] = -9999;

		// -- LHE -- //
		LHELepton_Px[i] = 0;
		LHELepton_Py[i] = 0;
		LHELepton_Pz[i] = 0;
		LHELepton_E[i] = 0;
		LHELepton_ID[i] = 0;
		LHELepton_status[i] = 0;
		
		// GEN
		GENLepton_phi[i] = GENLepton_eta[i] = GENLepton_pT[i] = GENLepton_mother[i] = -100;
		GENLepton_Px[i] = GENLepton_Py[i] = GENLepton_Pz[i] = -100;
		GENLepton_charge[i] = GENLepton_status[i] = GENLepton_ID[i] = -100;
		GENLepton_isPrompt[i] = 0;
		GENLepton_isPromptFinalState[i] = 0;
		GENLepton_isTauDecayProduct[i] = 0;
		GENLepton_isPromptTauDecayProduct[i] = 0;
		GENLepton_isDirectPromptTauDecayProductFinalState[i] = 0;
		GENLepton_isHardProcess[i] = 0;
		GENLepton_isLastCopy[i] = 0;
		GENLepton_isLastCopyBeforeFSR[i] = 0;
		GENLepton_isPromptDecayed[i] = 0;
		GENLepton_isDecayedLeptonHadron[i] = 0;
		GENLepton_fromHardProcessBeforeFSR[i] = 0;
		GENLepton_fromHardProcessDecayed[i] = 0;
		GENLepton_fromHardProcessFinalState[i] = 0;
		GENLepton_isMostlyLikePythia6Status3[i] = 0;
		GENEvt_weight = 0; //Weights for NLO generated events

		nGenOthers = -1;
		GenOthers_phi[i] = GenOthers_eta[i] = GenOthers_pT[i] = GenOthers_mother[i] = -100;
		GenOthers_Px[i] = GenOthers_Py[i] = GenOthers_Pz[i] = -100;
		GenOthers_charge[i] = GenOthers_status[i] = GenOthers_ID[i] = -100;
		GenOthers_isPrompt[i] = 0;
		GenOthers_isPromptFinalState[i] = 0;
		GenOthers_isTauDecayProduct[i] = 0;
		GenOthers_isPromptTauDecayProduct[i] = 0;
		GenOthers_isDirectPromptTauDecayProductFinalState[i] = 0;
		GenOthers_isHardProcess[i] = 0;
		GenOthers_isLastCopy[i] = 0;
		GenOthers_isLastCopyBeforeFSR[i] = 0;
		GenOthers_isPromptDecayed[i] = 0;
		GenOthers_isDecayedLeptonHadron[i] = 0;
		GenOthers_fromHardProcessBeforeFSR[i] = 0;
		GenOthers_fromHardProcessDecayed[i] = 0;
		GenOthers_fromHardProcessFinalState[i] = 0;
		GenOthers_isMostlyLikePythia6Status3[i] = 0;

		// -- Photon Information -- //
		nPhotons = 0;
		Photon_pT[i] = 0;
		Photon_eta[i] = 0;
		Photon_phi[i] = 0;
		Photon_etaSC[i] = 0;
		Photon_phiSC[i] = 0;
		Photon_HoverE[i] = 0;
		Photon_hasPixelSeed[i] = 0;
		Photon_Full5x5_SigmaIEtaIEta[i] = 0;
		Photon_ChIso[i] = 0;
		Photon_NhIso[i] = 0;
		Photon_PhIso[i] = 0;
		Photon_ChIsoWithEA[i] = 0;
		Photon_NhIsoWithEA[i] = 0;
		Photon_PhIsoWithEA[i] = 0;

		// -- MET -- //
		pfMET_pT = -100;
		pfMET_phi = -100; 
		pfMET_Px = -100; 
		pfMET_Py = -100; 
		pfMET_SumEt = -100;

		pfMET_Type1_pT = -100;
		pfMET_Type1_phi = -100; 
		pfMET_Type1_Px = -100; 
		pfMET_Type1_Py = -100; 
		pfMET_Type1_SumEt = -100;

	} // -- End of "i" iteration -- //

	// cout << "##### Analyze:Initialization #####" << endl;

	nEvt++;
	// -- run number & event number -- //
	runNum = iEvent.id().run();
	evtNum = iEvent.id().event();
	lumiBlock = iEvent.id().luminosityBlock();
	const bool isRD = iEvent.isRealData();

	// edm::Handle<double> weight_;
	// iEvent.getByLabel("PUweight", weight_);

	// if(weight_.isValid())
	// 	PUweight = *weight_;
	// else
	// 	PUweight = 1.0;

	//get the geometry
	edm::ESHandle<GlobalTrackingGeometry> glbTrackingGeometry;
	iSetup.get<GlobalTrackingGeometryRecord>().get(glbTrackingGeometry);

	// -- PileUp Reweighting -- //
	if( isMC )
	{
		edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
		iEvent.getByToken(PileUpInfoToken, PupInfo);
		std::vector<PileupSummaryInfo>::const_iterator PVI;

		int npv = -1;
		// int npvin = -1;

		for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI)
		{
			int BX = PVI->getBunchCrossing();

			if( BX == 0 )
			{
				// npvin = PVI->getPU_NumInteractions(); // in time only
				npv = PVI->getTrueNumInteractions(); // in and out of time
				continue;
			}

		}

		nPileUp = npv;

		//   pileUpReweightIn = LumiWeights_.weight( npvin );
		//   pileUpReweight = LumiWeights_.weight( npv );
		
		//   pileUpReweightPlus  = PShiftUp_.ShiftWeight( npv );
		//   pileUpReweightMinus = PShiftDown_.ShiftWeight( npv );

		//   pileUpReweightInMuonPhys = LumiWeightsMuonPhys_.weight( npvin );
		//   pileUpReweightMuonPhys = LumiWeightsMuonPhys_.weight( npv );
		
		//   pileUpReweightPlusMuonPhys  = PShiftUpMuonPhys_.ShiftWeight( npv );
		//   pileUpReweightMinusMuonPhys = PShiftDownMuonPhys_.ShiftWeight( npv );
	}

	// cout << "##### Analyze:PU-Reweighting #####" << endl;

	// fills
	if( theStoreHLTReportFlag ) hltReport(iEvent);
	if( theStorePriVtxFlag ) fillPrimaryVertex(iEvent);
	if( theStoreJetFlag ) fillJet(iEvent);
	if( theStoreMETFlag ) fillMET(iEvent);
	if( !isRD && theStoreLHEFlag ) fillLHEInfo(iEvent);
	if( !isRD && theStoreGENFlag ) fillGENInfo(iEvent);
	if( !isRD && theStoreGenOthersFlag ) fillGenOthersInfo(iEvent);
	if( theStorePhotonFlag ) fillPhotons(iEvent);
	if( theStoreMuonFlag ) fillMuons(iEvent, iSetup);
	if( theStoreElectronFlag ) fillElectrons(iEvent);
	if( theStoreTTFlag ) fillTT(iEvent);
	DYTree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------ //
void DYntupleMaker::beginJob()
{
	// if( isMC )
	// {
	// 	// Pileup Reweight: 2012, Summer12_S10
	// 	std::vector< float > _PUreweightRun2012 ;
	// 	std::vector< float > _PUreweightRun2012MuonPhys ;
	// 	std::vector< float > _MC2012;

	// 	for( int i = 0; i < 100; ++i)
	// 	{
	// 		_PUreweightRun2012.push_back((float)PileUpRD_[i]);
	// 		_PUreweightRun2012MuonPhys.push_back((float)PileUpRDMuonPhys_[i]);
	// 		_MC2012.push_back((float)PileUpMC_[i]);
	// 	}

	// 	LumiWeights_ = edm::LumiReWeighting(_MC2012, _PUreweightRun2012);
	// 	PShiftDown_ = reweight::PoissonMeanShifter(-0.5);
	// 	PShiftUp_ = reweight::PoissonMeanShifter(0.5);

	// 	LumiWeightsMuonPhys_ = edm::LumiReWeighting(_MC2012, _PUreweightRun2012MuonPhys);
	// 	PShiftDownMuonPhys_ = reweight::PoissonMeanShifter(-0.5);
	// 	PShiftUpMuonPhys_ = reweight::PoissonMeanShifter(0.5);
	// }

	edm::Service<TFileService> fs;
	DYTree = fs->make<TTree>("DYTree","DYTree");

	// -- global event variables -- //
	DYTree->Branch("nTotal",&nEvt,"nTotal/I");
	DYTree->Branch("runNum",&runNum,"runNum/I");
	DYTree->Branch("evtNum",&evtNum,"evtNum/l");
	DYTree->Branch("lumiBlock",&lumiBlock,"lumiBlock/I");
	DYTree->Branch("PUweight",&PUweight,"PUweight/D");
	// DYTree->Branch("sumEt",&sumEt,"sumEt/D");
	// DYTree->Branch("photonEt",&photonEt,"photonEt/D");
	// DYTree->Branch("chargedHadronEt",&chargedHadronEt,"chargedHadronEt/D");
	// DYTree->Branch("neutralHadronEt",&neutralHadronEt,"neutralHadronEt/D");
	DYTree->Branch("nVertices",&nVertices,"nVertices/I");

	DYTree->Branch("Flag_badMuons",&Flag_badMuons,"Flag_badMuons/O");
	DYTree->Branch("Flag_duplicateMuons",&Flag_duplicateMuons,"Flag_duplicateMuons/O");
	DYTree->Branch("Flag_noBadMuons",&Flag_noBadMuons,"Flag_noBadMuons/O");

	if(theStorePriVtxFlag)
	{
		DYTree->Branch("PVtrackSize", &PVtrackSize,"PVtrackSize/I");
		DYTree->Branch("PVchi2", &PVchi2,"PVchi2/D");
		DYTree->Branch("PVndof", &PVndof,"PVndof/D");
		DYTree->Branch("PVnormalizedChi2", &PVnormalizedChi2,"PVnormalizedChi2/D");
		DYTree->Branch("PVx", &PVx,"PVx/D");
		DYTree->Branch("PVy", &PVy,"PVy/D");
		DYTree->Branch("PVz", &PVz,"PVz/D");
	}

	if(theStoreHLTReportFlag)
	{
		DYTree->Branch("HLT_ntrig", &_HLT_ntrig,"HLT_ntrig/I");
		DYTree->Branch("HLT_trigType", &_HLT_trigType,"HLT_trigType[HLT_ntrig]/I");
		DYTree->Branch("HLT_trigFired", &_HLT_trigFired,"HLT_trigFired[HLT_ntrig]/I");
		DYTree->Branch("HLT_trigName", &_HLT_trigName);
		DYTree->Branch("HLT_trigPS", &_HLT_trigPS);
		DYTree->Branch("HLT_trigPt", &_HLT_trigPt,"HLT_trigPt[HLT_ntrig]/D");
		DYTree->Branch("HLT_trigEta", &_HLT_trigEta,"HLT_trigEta[HLT_ntrig]/D");
		DYTree->Branch("HLT_trigPhi", &_HLT_trigPhi,"HLT_trigPhi[HLT_ntrig]/D");
	}

	if(theStoreJetFlag)
	{
		// Jet
		DYTree->Branch("Njets", &Njets,"Njets/I");

		DYTree->Branch("Jet_pT", &Jet_pT,"Jet_pT[Njets]/D");
		DYTree->Branch("Jet_eta", &Jet_eta,"Jet_eta[Njets]/D");
		DYTree->Branch("Jet_phi", &Jet_phi,"Jet_phi[Njets]/D");
		DYTree->Branch("Jet_Charge", &Jet_Charge,"Jet_Charge[Njets]/D");
		DYTree->Branch("Jet_Flavor", &Jet_Flavor,"Jet_Flavor[Njets]/I");

		DYTree->Branch("Jet_bTag", &Jet_bTag,"Jet_bTag[Njets]/D");
		DYTree->Branch("Jet_CHfrac", &Jet_CHfrac,"Jet_CHfrac[Njets]/D");
		DYTree->Branch("Jet_NHfrac", &Jet_NHfrac,"Jet_NHfrac[Njets]/D");
		DYTree->Branch("Jet_NHEMfrac", &Jet_NHEMfrac,"Jet_NHEMfrac[Njets]/D");
		DYTree->Branch("Jet_CHEMfrac", &Jet_CHEMfrac,"Jet_CHEMfrac[Njets]/D");
		DYTree->Branch("Jet_CHmulti", &Jet_CHmulti,"Jet_CHmulti[Njets]/I");
		DYTree->Branch("Jet_NHmulti", &Jet_NHmulti,"Jet_NHmulti[Njets]/I");


		// DYTree->Branch("JETbDiscriminant", &JETbDiscriminant, "JETbDiscriminant[Njets]/D");
		// DYTree->Branch("JETbDiscriminant_alg1", &JETbDiscriminant_alg1, "JETbDiscriminant_alg1[Njets]/D");
		// DYTree->Branch("JETbDiscriminant_alg2", &JETbDiscriminant_alg2, "JETbDiscriminant_alg2[Njets]/D");
		// DYTree->Branch("JETbDiscriminant_alg3", &JETbDiscriminant_alg3, "JETbDiscriminant_alg3[Njets]/D");
		// DYTree->Branch("JETflavour", &JETflavour, "JETflavour[Njets]/I");
		// DYTree->Branch("JETcharge", &JETcharge, "JETcharge[Njets]/D");
		// DYTree->Branch("JETntracks", &JETntracks, "JETntracks[Njets]/I");
		// DYTree->Branch("JETpt", &JETpt, "JETpt[Njets]/D");
		// DYTree->Branch("JETeta", &JETeta, "JETeta[Njets]/D");
		// DYTree->Branch("JETphi", &JETphi, "JETphi[Njets]/D");
		// // b-tagging
		// DYTree->Branch("Nbtagged", &Nbtagged,"Nbtagged/I");
		// DYTree->Branch("NbtaggedCloseMuon", &NbtaggedCloseMuon,"NbtaggedCloseMuon/I");
		// DYTree->Branch("Nbtagged_alg1", &Nbtagged_alg1,"Nbtagged_alg1/I");
		// DYTree->Branch("NbtaggedCloseMuon_alg1", &NbtaggedCloseMuon_alg1,"NbtaggedCloseMuon_alg1/I");
		// DYTree->Branch("Nbtagged_alg2", &Nbtagged_alg2,"Nbtagged_alg2/I");
		// DYTree->Branch("NbtaggedCloseMuon_alg2", &NbtaggedCloseMuon_alg2,"NbtaggedCloseMuon_alg2/I");
		// DYTree->Branch("Nbtagged_alg3", &Nbtagged_alg3,"Nbtagged_alg3/I");
		// DYTree->Branch("NbtaggedCloseMuon_alg3", &NbtaggedCloseMuon_alg3,"NbtaggedCloseMuon_alg3/I");
	}

	// if (theStoreMETFlag)
	// {
	// 	DYTree->Branch("MET_sumEt", &MET_sumEt,"MET_sumEt/D");
	// 	DYTree->Branch("MET_pt", &MET_pt,"MET_pt/D");
	// 	DYTree->Branch("MET_px", &MET_px,"MET_px/D");
	// 	DYTree->Branch("MET_py", &MET_py,"MET_py/D");
	// 	DYTree->Branch("MET_phi", &MET_phi,"MET_phi/D");
	// 	DYTree->Branch("pfMET_sumEt", &pfMET_sumEt,"pfMET_sumEt/D");
	// 	DYTree->Branch("pfMET_pt", &pfMET_pt,"pfMET_pt/D");
	// 	DYTree->Branch("pfMET_px", &pfMET_px,"pfMET_px/D");
	// 	DYTree->Branch("pfMET_py", &pfMET_py,"pfMET_py/D");
	// 	DYTree->Branch("pfMET_phi", &pfMET_phi,"pfMET_phi/D");
	// }

	// Electron
	if( theStoreElectronFlag )
	{
		DYTree->Branch("Nelectrons", &Nelectrons,"Nelectrons/I");
		DYTree->Branch("Electron_Energy", &Electron_Energy, "Electron_Energy[Nelectrons]/D");
		// DYTree->Branch("Electron_et", &Electron_et, "Electron_et[Nelectrons]/D");
		// DYTree->Branch("Electron_etCorr", &Electron_etCorr, "Electron_etCorr[Nelectrons]/D");
		// DYTree->Branch("Electron_caloEnergy", &Electron_caloEnergy, "Electron_caloEnergy[Nelectrons]/D");
		DYTree->Branch("Electron_pT", &Electron_pT, "Electron_pT[Nelectrons]/D");
		// DYTree->Branch("Electron_Px", &Electron_Px, "Electron_Px[Nelectrons]/D");
		// DYTree->Branch("Electron_Py", &Electron_Py, "Electron_Py[Nelectrons]/D");
		// DYTree->Branch("Electron_Pz", &Electron_Pz, "Electron_Pz[Nelectrons]/D");
		// DYTree->Branch("Electron_ECalib", &Electron_ECalib, "Electron_ECalib[Nelectrons]/D");
		// DYTree->Branch("Electron_etCalib", &Electron_etCalib, "Electron_etCalib[Nelectrons]/D");
		// DYTree->Branch("Electron_etCorrCalib", &Electron_etCorrCalib, "Electron_etCorrCalib[Nelectrons]/D");
		// DYTree->Branch("Electron_caloEnergyCalib", &Electron_caloEnergyCalib, "Electron_caloEnergyCalib[Nelectrons]/D");
		// DYTree->Branch("Electron_pTCalib", &Electron_pTCalib, "Electron_pTCalib[Nelectrons]/D");
		// DYTree->Branch("Electron_PxCalib", &Electron_PxCalib, "Electron_PxCalib[Nelectrons]/D");
		// DYTree->Branch("Electron_PyCalib", &Electron_PyCalib, "Electron_PyCalib[Nelectrons]/D");
		// DYTree->Branch("Electron_PzCalib", &Electron_PzCalib, "Electron_PzCalib[Nelectrons]/D");
		DYTree->Branch("Electron_eta", &Electron_eta, "Electron_eta[Nelectrons]/D");
		DYTree->Branch("Electron_phi", &Electron_phi, "Electron_phi[Nelectrons]/D");
		DYTree->Branch("Electron_charge", &Electron_charge, "Electron_charge[Nelectrons]/I");
		DYTree->Branch("Electron_gsfpT", &Electron_gsfpT, "Electron_gsfpT[Nelectrons]/D");
		DYTree->Branch("Electron_gsfPx", &Electron_gsfPx, "Electron_gsfPx[Nelectrons]/D");
		DYTree->Branch("Electron_gsfPy", &Electron_gsfPy, "Electron_gsfPy[Nelectrons]/D");
		DYTree->Branch("Electron_gsfPz", &Electron_gsfPz, "Electron_gsfPz[Nelectrons]/D");
		// DYTree->Branch("Electron_gsfpTCalib", &Electron_gsfpTCalib, "Electron_gsfpTCalib[Nelectrons]/D");
		// DYTree->Branch("Electron_gsfPxCalib", &Electron_gsfPxCalib, "Electron_gsfPxCalib[Nelectrons]/D");
		// DYTree->Branch("Electron_gsfPyCalib", &Electron_gsfPyCalib, "Electron_gsfPyCalib[Nelectrons]/D");
		// DYTree->Branch("Electron_gsfPzCalib", &Electron_gsfPzCalib, "Electron_gsfPzCalib[Nelectrons]/D");
		DYTree->Branch("Electron_gsfEta", &Electron_gsfEta, "Electron_gsfEta[Nelectrons]/D");
		DYTree->Branch("Electron_gsfPhi", &Electron_gsfPhi, "Electron_gsfPhi[Nelectrons]/D");
		DYTree->Branch("Electron_gsfCharge", &Electron_gsfCharge, "Electron_gsfCharge[Nelectrons]/I");
		DYTree->Branch("Electron_etaSC", &Electron_etaSC, "Electron_etaSC[Nelectrons]/D");
		DYTree->Branch("Electron_phiSC", &Electron_phiSC, "Electron_phiSC[Nelectrons]/D");
		DYTree->Branch("Electron_etaWidth", &Electron_etaWidth, "Electron_etaWidth[Nelectrons]/D");
		DYTree->Branch("Electron_phiWidth", &Electron_phiWidth, "Electron_phiWidth[Nelectrons]/D");
		DYTree->Branch("Electron_dEtaIn", &Electron_dEtaIn, "Electron_dEtaIn[Nelectrons]/D");
		DYTree->Branch("Electron_dPhiIn", &Electron_dPhiIn, "Electron_dPhiIn[Nelectrons]/D");
		DYTree->Branch("Electron_sigmaIEtaIEta", &Electron_sigmaIEtaIEta, "Electron_sigmaIEtaIEta[Nelectrons]/D");
		DYTree->Branch("Electron_HoverE", &Electron_HoverE, "Electron_HoverE[Nelectrons]/D");
		DYTree->Branch("Electron_fbrem", &Electron_fbrem, "Electron_fbrem[Nelectrons]/D");
		DYTree->Branch("Electron_eOverP", &Electron_eOverP, "Electron_eOverP[Nelectrons]/D");
		// DYTree->Branch("Electron_energyEC", &Electron_energyEC, "Electron_energyEC[Nelectrons]/D");
		// DYTree->Branch("Electron_Pnorm", &Electron_Pnorm, "Electron_Pnorm[Nelectrons]/D");
		DYTree->Branch("Electron_InvEminusInvP", &Electron_InvEminusInvP, "Electron_InvEminusInvP[Nelectrons]/D");
		DYTree->Branch("Electron_dxyVTX", &Electron_dxyVTX, "Electron_dxyVTX[Nelectrons]/D");
		DYTree->Branch("Electron_dzVTX", &Electron_dzVTX, "Electron_dzVTX[Nelectrons]/D");
		DYTree->Branch("Electron_dxy", &Electron_dxy, "Electron_dxy[Nelectrons]/D");
		DYTree->Branch("Electron_dz", &Electron_dz, "Electron_dz[Nelectrons]/D");
		DYTree->Branch("Electron_dxyBS", &Electron_dxyBS, "Electron_dxyBS[Nelectrons]/D");
		DYTree->Branch("Electron_dzBS", &Electron_dzBS, "Electron_dzBS[Nelectrons]/D");
		// DYTree->Branch("Electron_AEff03", &Electron_AEff03, "Electron_AEff03[Nelectrons]/D");
		DYTree->Branch("Electron_chIso03", &Electron_chIso03, "Electron_chIso03[Nelectrons]/D");
		// DYTree->Branch("Electron_chIso04", &Electron_chIso04, "Electron_chIso04[Nelectrons]/D");
		DYTree->Branch("Electron_nhIso03", &Electron_nhIso03, "Electron_nhIso03[Nelectrons]/D");
		// DYTree->Branch("Electron_nhIso04", &Electron_nhIso04, "Electron_nhIso04[Nelectrons]/D");
		DYTree->Branch("Electron_phIso03", &Electron_phIso03, "Electron_phIso03[Nelectrons]/D");
		// DYTree->Branch("Electron_phIso04", &Electron_phIso04, "Electron_phIso04[Nelectrons]/D");
		// DYTree->Branch("Electron_pcIso03", &Electron_pcIso03, "Electron_pcIso03[Nelectrons]/D");
		// DYTree->Branch("Electron_pcIso04", &Electron_pcIso04, "Electron_pcIso04[Nelectrons]/D");
		// DYTree->Branch("Electron_relIsoCom03", &Electron_relIsoCom03, "Electron_relIsoCom03[Nelectrons]/D");
		// DYTree->Branch("Electron_relIsoCom04", &Electron_relIsoCom04, "Electron_relIsoCom04[Nelectrons]/D");
		// DYTree->Branch("Electron_relIsoBeta03", &Electron_relIsoBeta03, "Electron_relIsoBeta03[Nelectrons]/D");
		// DYTree->Branch("Electron_relIsoBeta04", &Electron_relIsoBeta04, "Electron_relIsoBeta04[Nelectrons]/D");
		// DYTree->Branch("Electron_relIsoRho03", &Electron_relIsoRho03, "Electron_relIsoRho03[Nelectrons]/D");
		DYTree->Branch("Electron_hasConversion", &Electron_hasConversion, "Electron_hasConversion[Nelectrons]/O");
		DYTree->Branch("Electron_mHits", &Electron_mHits, "Electron_mHits[Nelectrons]/I");

		DYTree->Branch("Electron_EnergySC", &Electron_EnergySC, "Electron_EnergySC[Nelectrons]/D");
		DYTree->Branch("Electron_preEnergySC", &Electron_preEnergySC, "Electron_preEnergySC[Nelectrons]/D");
		DYTree->Branch("Electron_rawEnergySC", &Electron_rawEnergySC, "Electron_rawEnergySC[Nelectrons]/D");
		DYTree->Branch("Electron_etSC", &Electron_etSC, "Electron_etSC[Nelectrons]/D");
		DYTree->Branch("Electron_E15", &Electron_E15, "Electron_E15[Nelectrons]/D");
		DYTree->Branch("Electron_E25", &Electron_E25, "Electron_E25[Nelectrons]/D");
		DYTree->Branch("Electron_E55", &Electron_E55, "Electron_E55[Nelectrons]/D");
		DYTree->Branch("Electron_ChIso03FromPU", &Electron_ChIso03FromPU, "Electron_ChIso03FromPU[Nelectrons]/D");
		DYTree->Branch("Electron_RelPFIso_dBeta", &Electron_RelPFIso_dBeta, "Electron_RelPFIso_dBeta[Nelectrons]/D");
		DYTree->Branch("Electron_RelPFIso_Rho", &Electron_RelPFIso_Rho, "Electron_RelPFIso_Rho[Nelectrons]/D");
		DYTree->Branch("Electron_passConvVeto", &Electron_passConvVeto, "Electron_passConvVeto[Nelectrons]/O"); // O (letter; not zero): Boolean 
		DYTree->Branch("Electron_passVetoID", &Electron_passVetoID, "Electron_passVetoID[Nelectrons]/O");
		DYTree->Branch("Electron_passLooseID", &Electron_passLooseID, "Electron_passLooseID[Nelectrons]/O");
		DYTree->Branch("Electron_passMediumID", &Electron_passMediumID, "Electron_passMediumID[Nelectrons]/O");
		DYTree->Branch("Electron_passTightID", &Electron_passTightID, "Electron_passTightID[Nelectrons]/O");
		DYTree->Branch("Electron_passMVAID_WP80", &Electron_passMVAID_WP80, "Electron_passMVAID_WP80[Nelectrons]/O");
		DYTree->Branch("Electron_passMVAID_WP90", &Electron_passMVAID_WP90, "Electron_passMVAID_WP90[Nelectrons]/O");
		DYTree->Branch("Electron_passHEEPID", &Electron_passHEEPID, "Electron_passHEEPID[Nelectrons]/O");
		DYTree->Branch("Electron_r9", &Electron_r9, "Electron_r9[Nelectrons]/D");

		// DYTree->Branch("Electron_crack", &Electron_crack, "Electron_crack[Nelectrons]/I");
		DYTree->Branch("Electron_ecalDriven", &Electron_ecalDriven, "Electron_ecalDriven[Nelectrons]/I");
		// DYTree->Branch("Electron_isoEMHADDepth1", &Electron_isoEMHADDepth1, "Electron_isoEMHADDepth1[Nelectrons]/D");
		// DYTree->Branch("Electron_25over55", &Electron_25over55, "Electron_25over55[Nelectrons]/D");
		// DYTree->Branch("Electron_15over55", &Electron_15over55, "Electron_15over55[Nelectrons]/D");
		// DYTree->Branch("Electron_isoHADDepth2", &Electron_isoHADDepth2, "Electron_isoHADDepth2[Nelectrons]/D");
		// DYTree->Branch("Electron_isoPtTrks", &Electron_isoPtTrks, "Electron_isoPtTrks[Nelectrons]/D");
		// DYTree->Branch("Electron_modIsoEMHADDepth1", &Electron_modIsoEMHADDepth1, "Electron_modIsoEMHADDepth1[Nelectrons]/D");
		// DYTree->Branch("Electron_modIsoPtTrks", &Electron_modIsoPtTrks, "Electron_modIsoPtTrks[Nelectrons]/D");
		// DYTree->Branch("Electron_modIsoEMHADDepth1Orig", &Electron_modIsoEMHADDepth1Orig, "Electron_modIsoEMHADDepth1Orig[Nelectrons]/D");
		// DYTree->Branch("Electron_modIsoPtTrksOrig", &Electron_modIsoPtTrksOrig, "Electron_modIsoPtTrksOrig[Nelectrons]/D");
		DYTree->Branch("Electron_ambGsf0Pt", &Electron_ambGsf0Pt, "Electron_ambGsf0Pt[Nelectrons]/D");
		DYTree->Branch("Electron_ambGsf0Eta", &Electron_ambGsf0Eta, "Electron_ambGsf0Eta[Nelectrons]/D");
		DYTree->Branch("Electron_ambGsf0Phi", &Electron_ambGsf0Phi, "Electron_ambGsf0Phi[Nelectrons]/D");
		DYTree->Branch("Electron_ambGsf0Charge", &Electron_ambGsf0Charge, "Electron_ambGsf0Charge[Nelectrons]/D");
		DYTree->Branch("Electron_ambGsf1Pt", &Electron_ambGsf1Pt, "Electron_ambGsf1Pt[Nelectrons]/D");
		DYTree->Branch("Electron_ambGsf1Eta", &Electron_ambGsf1Eta, "Electron_ambGsf1Eta[Nelectrons]/D");
		DYTree->Branch("Electron_ambGsf1Phi", &Electron_ambGsf1Phi, "Electron_ambGsf1Phi[Nelectrons]/D");
		DYTree->Branch("Electron_ambGsf1Charge", &Electron_ambGsf1Charge, "Electron_ambGsf1Charge[Nelectrons]/D");
		DYTree->Branch("Electron_ambGsf2Pt", &Electron_ambGsf2Pt, "Electron_ambGsf2Pt[Nelectrons]/D");
		DYTree->Branch("Electron_ambGsf2Eta", &Electron_ambGsf2Eta, "Electron_ambGsf2Eta[Nelectrons]/D");
		DYTree->Branch("Electron_ambGsf2Phi", &Electron_ambGsf2Phi, "Electron_ambGsf2Phi[Nelectrons]/D");
		DYTree->Branch("Electron_ambGsf2Charge", &Electron_ambGsf2Charge, "Electron_ambGsf2Charge[Nelectrons]/D");
		DYTree->Branch("Electron_ambGsf3Pt", &Electron_ambGsf3Pt, "Electron_ambGsf3Pt[Nelectrons]/D");
		DYTree->Branch("Electron_ambGsf3Eta", &Electron_ambGsf3Eta, "Electron_ambGsf3Eta[Nelectrons]/D");
		DYTree->Branch("Electron_ambGsf3Phi", &Electron_ambGsf3Phi, "Electron_ambGsf3Phi[Nelectrons]/D");
		DYTree->Branch("Electron_ambGsf3Charge", &Electron_ambGsf3Charge, "Electron_ambGsf3Charge[Nelectrons]/D");
	}

	// -- muon variables -- //
	if( theStoreMuonFlag )
	{
		DYTree->Branch("nMuon",&nMuon,"nMuon/I");
		DYTree->Branch("Nmuons",&Nmuons,"Nmuons/I");
		DYTree->Branch("Muon_muonType", &Muon_muonType,"Muon_muonType[nMuon]/I");
		DYTree->Branch("isPFmuon", &isPFmuon, "isPFmuon[nMuon]/I");
		DYTree->Branch("isGLBmuon", &isGLBmuon, "isGLBmuon[nMuon]/I");
		DYTree->Branch("isTRKmuon", &isTRKmuon, "isTRKmuon[nMuon]/I");
		DYTree->Branch("isSTAmuon", &isSTAmuon, "isSTAmuon[nMuon]/I");
		DYTree->Branch("Muon_nTrig", &Muon_nTrig,"Muon_nTrig[nMuon]/I");
		DYTree->Branch("Muon_triggerObjectType", &Muon_triggerObjectType,"Muon_triggerObjectType[nMuon]/I");
		DYTree->Branch("Muon_filterName", &Muon_filterName,"Muon_filterName[nMuon]/I");
		DYTree->Branch("Muon_phi", &Muon_phi,"Muon_phi[nMuon]/D");
		DYTree->Branch("Muon_dB", &Muon_dB, "Muon_dB[nMuon]/D");
		DYTree->Branch("Muon_eta", &Muon_eta,"Muon_eta[nMuon]/D");
		DYTree->Branch("Muon_pT", &Muon_pT,"Muon_pT[nMuon]/D");
		DYTree->Branch("Muon_cktpT", &Muon_cktpT,"Muon_cktpT[nMuon]/D");
		DYTree->Branch("Muon_cktPx", &Muon_cktPx,"Muon_cktPx[nMuon]/D");
		DYTree->Branch("Muon_cktPy", &Muon_cktPy,"Muon_cktPy[nMuon]/D");
		DYTree->Branch("Muon_cktPz", &Muon_cktPz,"Muon_cktPz[nMuon]/D");
		DYTree->Branch("Muon_cktpTError", &Muon_cktpTError,"Muon_cktpTError[nMuon]/D");
		DYTree->Branch("Muon_Px", &Muon_Px,"Muon_Px[nMuon]/D");
		DYTree->Branch("Muon_Py", &Muon_Py,"Muon_Py[nMuon]/D");
		DYTree->Branch("Muon_Pz", &Muon_Pz,"Muon_Pz[nMuon]/D");
		DYTree->Branch("Muon_trkiso", &Muon_trkiso,"Muon_trkiso[nMuon]/D");
		DYTree->Branch("Muon_hcaliso", &Muon_hcaliso,"Muon_hcaliso[nMuon]/D");
		DYTree->Branch("Muon_ecaliso", &Muon_ecaliso,"Muon_ecaliso[nMuon]/D");
		DYTree->Branch("Muon_trkisoR05", &Muon_trkisoR05,"Muon_trkisoR05[nMuon]/D");
		DYTree->Branch("Muon_hcalisoR05", &Muon_hcalisoR05,"Muon_hcalisoR05[nMuon]/D");
		DYTree->Branch("Muon_ecalisoR05", &Muon_ecalisoR05,"Muon_ecalisoR05[nMuon]/D");

		//Various track informations
		DYTree->Branch("Muon_Best_pT", &Muon_Best_pT, "Muon_Best_pT[nMuon]/D");
		DYTree->Branch("Muon_Best_pTError", &Muon_Best_pTError, "Muon_Best_pTError[nMuon]/D");
		DYTree->Branch("Muon_Best_Px", &Muon_Best_Px, "Muon_Best_Px[nMuon]/D");
		DYTree->Branch("Muon_Best_Py", &Muon_Best_Py, "Muon_Best_Py[nMuon]/D");
		DYTree->Branch("Muon_Best_Pz", &Muon_Best_Pz, "Muon_Best_Pz[nMuon]/D");
		DYTree->Branch("Muon_Best_eta", &Muon_Best_eta, "Muon_Best_eta[nMuon]/D");
		DYTree->Branch("Muon_Best_phi", &Muon_Best_phi, "Muon_Best_phi[nMuon]/D");

		DYTree->Branch("Muon_Inner_pT", &Muon_Inner_pT, "Muon_Inner_pT[nMuon]/D");
		DYTree->Branch("Muon_Inner_pTError", &Muon_Inner_pTError, "Muon_Inner_pTError[nMuon]/D");
		DYTree->Branch("Muon_Inner_Px", &Muon_Inner_Px, "Muon_Inner_Px[nMuon]/D");
		DYTree->Branch("Muon_Inner_Py", &Muon_Inner_Py, "Muon_Inner_Py[nMuon]/D");
		DYTree->Branch("Muon_Inner_Pz", &Muon_Inner_Pz, "Muon_Inner_Pz[nMuon]/D");
		DYTree->Branch("Muon_Inner_eta", &Muon_Inner_eta, "Muon_Inner_eta[nMuon]/D");
		DYTree->Branch("Muon_Inner_phi", &Muon_Inner_phi, "Muon_Inner_phi[nMuon]/D");

		DYTree->Branch("Muon_Outer_pT", &Muon_Outer_pT, "Muon_Outer_pT[nMuon]/D");
		DYTree->Branch("Muon_Outer_pTError", &Muon_Outer_pTError, "Muon_Outer_pTError[nMuon]/D");
		DYTree->Branch("Muon_Outer_Px", &Muon_Outer_Px, "Muon_Outer_Px[nMuon]/D");
		DYTree->Branch("Muon_Outer_Py", &Muon_Outer_Py, "Muon_Outer_Py[nMuon]/D");
		DYTree->Branch("Muon_Outer_Pz", &Muon_Outer_Pz, "Muon_Outer_Pz[nMuon]/D");
		DYTree->Branch("Muon_Outer_eta", &Muon_Outer_eta, "Muon_Outer_eta[nMuon]/D");
		DYTree->Branch("Muon_Outer_phi", &Muon_Outer_phi, "Muon_Outer_phi[nMuon]/D");

		DYTree->Branch("Muon_GLB_pT", &Muon_GLB_pT, "Muon_GLB_pT[nMuon]/D");
		DYTree->Branch("Muon_GLB_pTError", &Muon_GLB_pTError, "Muon_GLB_pTError[nMuon]/D");
		DYTree->Branch("Muon_GLB_Px", &Muon_GLB_Px, "Muon_GLB_Px[nMuon]/D");
		DYTree->Branch("Muon_GLB_Py", &Muon_GLB_Py, "Muon_GLB_Py[nMuon]/D");
		DYTree->Branch("Muon_GLB_Pz", &Muon_GLB_Pz, "Muon_GLB_Pz[nMuon]/D");
		DYTree->Branch("Muon_GLB_eta", &Muon_GLB_eta, "Muon_GLB_eta[nMuon]/D");
		DYTree->Branch("Muon_GLB_phi", &Muon_GLB_phi, "Muon_GLB_phi[nMuon]/D");

		DYTree->Branch("Muon_TuneP_pT", &Muon_TuneP_pT, "Muon_TuneP_pT[nMuon]/D");
		DYTree->Branch("Muon_TuneP_pTError", &Muon_TuneP_pTError, "Muon_TuneP_pTError[nMuon]/D");
		DYTree->Branch("Muon_TuneP_Px", &Muon_TuneP_Px, "Muon_TuneP_Px[nMuon]/D");
		DYTree->Branch("Muon_TuneP_Py", &Muon_TuneP_Py, "Muon_TuneP_Py[nMuon]/D");
		DYTree->Branch("Muon_TuneP_Pz", &Muon_TuneP_Pz, "Muon_TuneP_Pz[nMuon]/D");
		DYTree->Branch("Muon_TuneP_eta", &Muon_TuneP_eta, "Muon_TuneP_eta[nMuon]/D");
		DYTree->Branch("Muon_TuneP_phi", &Muon_TuneP_phi, "Muon_TuneP_phi[nMuon]/D");
		       
		//pf iso
		DYTree->Branch("Muon_PfChargedHadronIsoR05", &Muon_PfChargedHadronIsoR05,"Muon_PfChargedHadronIsoR05[nMuon]/D");
		DYTree->Branch("Muon_PfNeutralHadronIsoR05", &Muon_PfNeutralHadronIsoR05,"Muon_PfNeutralHadronIsoR05[nMuon]/D");
		DYTree->Branch("Muon_PfGammaIsoR05", &Muon_PfGammaIsoR05,"Muon_PfGammaIsoR05[nMuon]/D");
		DYTree->Branch("Muon_PfChargedHadronIsoR04", &Muon_PfChargedHadronIsoR04,"Muon_PfChargedHadronIsoR04[nMuon]/D");
		DYTree->Branch("Muon_PfNeutralHadronIsoR04", &Muon_PfNeutralHadronIsoR04,"Muon_PfNeutralHadronIsoR04[nMuon]/D");
		DYTree->Branch("Muon_PfGammaIsoR04", &Muon_PfGammaIsoR04,"Muon_PfGammaIsoR04[nMuon]/D");
		DYTree->Branch("Muon_PFSumPUIsoR04", &Muon_PFSumPUIsoR04, "Muon_PFSumPUIsoR04[nMuon]/D");
		DYTree->Branch("Muon_PfChargedHadronIsoR03", &Muon_PfChargedHadronIsoR03,"Muon_PfChargedHadronIsoR03[nMuon]/D");
		DYTree->Branch("Muon_PfNeutralHadronIsoR03", &Muon_PfNeutralHadronIsoR03,"Muon_PfNeutralHadronIsoR03[nMuon]/D");
		DYTree->Branch("Muon_PfGammaIsoR03", &Muon_PfGammaIsoR03,"Muon_PfGammaIsoR03[nMuon]/D");
		DYTree->Branch("Muon_PFSumPUIsoR03", &Muon_PFSumPUIsoR03, "Muon_PFSumPUIsoR03[nMuon]/D");
		DYTree->Branch("Muon_charge", &Muon_charge,"Muon_charge[nMuon]/I");
		DYTree->Branch("Muon_nChambers", &Muon_nChambers,"Muon_nChambers[nMuon]/I");
		DYTree->Branch("Muon_nMatches", &Muon_nMatches,"Muon_nMatches[nMuon]/I");
		DYTree->Branch("Muon_nMatchesRPCLayers", &Muon_nMatchesRPCLayers,"Muon_nMatchesRPCLayers[nMuon]/I");
		DYTree->Branch("Muon_stationMask", &Muon_stationMask,"Muon_stationMask[nMuon]/I");
		DYTree->Branch("Muon_nSegments", &Muon_nSegments,"Muon_nSegments[nMuon]/I");
		DYTree->Branch("Muon_chi2dof", &Muon_chi2dof,"Muon_chi2dof[nMuon]/D");
		DYTree->Branch("Muon_nhits", &Muon_nhits,"Muon_nhits[nMuon]/I");
		DYTree->Branch("Muon_trackerHits", &Muon_trackerHits,"Muon_trackerHits[nMuon]/I");
		DYTree->Branch("Muon_trackerLayers", &Muon_trackerLayers,"Muon_trackerLayers[nMuon]/I");
		DYTree->Branch("Muon_pixelHits", &Muon_pixelHits,"Muon_pixelHits[nMuon]/I");
		DYTree->Branch("Muon_trackerHitsGLB", &Muon_trackerHitsGLB,"Muon_trackerHitsGLB[nMuon]/I");
		DYTree->Branch("Muon_trackerLayersGLB", &Muon_trackerLayersGLB,"Muon_trackerLayersGLB[nMuon]/I");
		DYTree->Branch("Muon_pixelHitsGLB", &Muon_pixelHitsGLB,"Muon_pixelHitsGLB[nMuon]/I");
		DYTree->Branch("Muon_muonHits", &Muon_muonHits,"Muon_muonHits[nMuon]/I");
		DYTree->Branch("Muon_qoverp", &Muon_qoverp,"Muon_qoverp[nMuon]/D");
		DYTree->Branch("Muon_theta", &Muon_theta,"Muon_theta[nMuon]/D");
		DYTree->Branch("Muon_lambda", &Muon_lambda,"Muon_lambda[nMuon]/D");
		DYTree->Branch("Muon_dxy", &Muon_dxy,"Muon_dxy[nMuon]/D");
		DYTree->Branch("Muon_d0", &Muon_d0,"Muon_d0[nMuon]/D");
		DYTree->Branch("Muon_dsz", &Muon_dsz,"Muon_dsz[nMuon]/D");
		DYTree->Branch("Muon_dz", &Muon_dz,"Muon_dz[nMuon]/D");
		DYTree->Branch("Muon_dxyBS", &Muon_dxyBS,"Muon_dxyBS[nMuon]/D");
		DYTree->Branch("Muon_dszBS", &Muon_dszBS,"Muon_dszBS[nMuon]/D");
		DYTree->Branch("Muon_dzBS", &Muon_dzBS,"Muon_dzBS[nMuon]/D");
		DYTree->Branch("Muon_dxyVTX", &Muon_dxyVTX,"Muon_dxyVTX[nMuon]/D");
		DYTree->Branch("Muon_dszVTX", &Muon_dszVTX,"Muon_dszVTX[nMuon]/D");
		DYTree->Branch("Muon_dzVTX", &Muon_dzVTX,"Muon_dzVTX[nMuon]/D");
		DYTree->Branch("Muon_dxycktVTX", &Muon_dxycktVTX,"Muon_dxycktVTX[nMuon]/D");
		DYTree->Branch("Muon_dszcktVTX", &Muon_dszcktVTX,"Muon_dszcktVTX[nMuon]/D");
		DYTree->Branch("Muon_dzcktVTX", &Muon_dzcktVTX,"Muon_dzcktVTX[nMuon]/D");
		DYTree->Branch("Muon_vx", &Muon_vx,"Muon_vx[nMuon]/D");
		DYTree->Branch("Muon_vy", &Muon_vy,"Muon_vy[nMuon]/D");
		DYTree->Branch("Muon_vz", &Muon_vz,"Muon_vz[nMuon]/D");
		DYTree->Branch("CosAngle", &CosAngle);
		DYTree->Branch("vtxTrkCkt1Pt", &vtxTrkCkt1Pt);
		DYTree->Branch("vtxTrkCkt2Pt", &vtxTrkCkt2Pt);
		DYTree->Branch("vtxTrkChi2", &vtxTrkChi2);
		DYTree->Branch("vtxTrkProb", &vtxTrkProb);
		DYTree->Branch("vtxTrkNdof", &vtxTrkNdof);
		DYTree->Branch("vtxTrkDiE1Pt", &vtxTrkDiE1Pt);
		DYTree->Branch("vtxTrkDiE2Pt", &vtxTrkDiE2Pt);
		DYTree->Branch("vtxTrkDiEChi2", &vtxTrkDiEChi2);
		DYTree->Branch("vtxTrkDiEProb", &vtxTrkDiEProb);
		DYTree->Branch("vtxTrkDiENdof", &vtxTrkDiENdof);
		DYTree->Branch("vtxTrkEMu1Pt", &vtxTrkEMu1Pt);
		DYTree->Branch("vtxTrkEMu2Pt", &vtxTrkEMu2Pt);
		DYTree->Branch("vtxTrkEMuChi2", &vtxTrkEMuChi2);
		DYTree->Branch("vtxTrkEMuProb", &vtxTrkEMuProb);
		DYTree->Branch("vtxTrkEMuNdof", &vtxTrkEMuNdof);

		DYTree->Branch("CosAngle_TuneP", &CosAngle_TuneP);
		DYTree->Branch("vtxTrk1Pt_TuneP", &vtxTrk1Pt_TuneP);
		DYTree->Branch("vtxTrk2Pt_TuneP", &vtxTrk2Pt_TuneP);
		DYTree->Branch("vtxTrkChi2_TuneP", &vtxTrkChi2_TuneP);
		DYTree->Branch("vtxTrkNdof_TuneP", &vtxTrkNdof_TuneP);
		DYTree->Branch("vtxTrkProb_TuneP", &vtxTrkProb_TuneP);
	}

	// -- LHE info -- //
	if( theStoreLHEFlag )
	{
		DYTree->Branch("nLHEParticle",&nLHEParticle,"nLHEParticle/I");
		DYTree->Branch("LHELepton_Px", &LHELepton_Px,"LHELepton_Px[nLHEParticle]/D");
		DYTree->Branch("LHELepton_Py", &LHELepton_Py,"LHELepton_Py[nLHEParticle]/D");
		DYTree->Branch("LHELepton_Pz", &LHELepton_Pz,"LHELepton_Pz[nLHEParticle]/D");
		DYTree->Branch("LHELepton_E", &LHELepton_E,"LHELepton_E[nLHEParticle]/D");
		DYTree->Branch("LHELepton_ID", &LHELepton_ID,"LHELepton_ID[nLHEParticle]/I");
		DYTree->Branch("LHELepton_status", &LHELepton_status,"LHELepton_status[nLHEParticle]/I");
	}

	// GEN info
	if( theStoreGENFlag )
	{
		DYTree->Branch("GENnPair",&GENnPair,"GENnPair/I");
		DYTree->Branch("GENLepton_phi", &GENLepton_phi,"GENLepton_phi[GENnPair]/D");
		DYTree->Branch("GENLepton_eta", &GENLepton_eta,"GENLepton_eta[GENnPair]/D");
		DYTree->Branch("GENLepton_pT", &GENLepton_pT,"GENLepton_pT[GENnPair]/D");
		DYTree->Branch("GENLepton_Px", &GENLepton_Px,"GENLepton_Px[GENnPair]/D");
		DYTree->Branch("GENLepton_Py", &GENLepton_Py,"GENLepton_Py[GENnPair]/D");
		DYTree->Branch("GENLepton_Pz", &GENLepton_Pz,"GENLepton_Pz[GENnPair]/D");
		DYTree->Branch("GENLepton_mother", &GENLepton_mother,"GENLepton_mother[GENnPair]/D");
		DYTree->Branch("GENLepton_charge", &GENLepton_charge,"GENLepton_charge[GENnPair]/I");
		DYTree->Branch("GENLepton_status", &GENLepton_status,"GENLepton_status[GENnPair]/I");
		DYTree->Branch("GENLepton_ID", &GENLepton_ID,"GENLepton_ID[GENnPair]/I");
		DYTree->Branch("GENLepton_isPrompt", &GENLepton_isPrompt,"GENLepton_isPrompt[GENnPair]/I");
		DYTree->Branch("GENLepton_isPromptFinalState", &GENLepton_isPromptFinalState,"GENLepton_isPromptFinalState[GENnPair]/I");
		DYTree->Branch("GENLepton_isTauDecayProduct", &GENLepton_isTauDecayProduct,"GENLepton_isTauDecayProduct[GENnPair]/I");
		DYTree->Branch("GENLepton_isPromptTauDecayProduct", &GENLepton_isPromptTauDecayProduct,"GENLepton_isPromptTauDecayProduct[GENnPair]/I");
		DYTree->Branch("GENLepton_isDirectPromptTauDecayProductFinalState", &GENLepton_isDirectPromptTauDecayProductFinalState,"GENLepton_isDirectPromptTauDecayProductFinalState[GENnPair]/I");
		DYTree->Branch("GENLepton_isHardProcess",&GENLepton_isHardProcess,"GENLepton_isHardProcess[GENnPair]/I");
		DYTree->Branch("GENLepton_isLastCopy",&GENLepton_isLastCopy,"GENLepton_isLastCopy[GENnPair]/I");
		DYTree->Branch("GENLepton_isLastCopyBeforeFSR",&GENLepton_isLastCopyBeforeFSR,"GENLepton_isLastCopyBeforeFSR[GENnPair]/I");
		DYTree->Branch("GENLepton_isPromptDecayed",&GENLepton_isPromptDecayed,"GENLepton_isPromptDecayed[GENnPair]/I");
		DYTree->Branch("GENLepton_isDecayedLeptonHadron",&GENLepton_isDecayedLeptonHadron,"GENLepton_isDecayedLeptonHadron[GENnPair]/I");
		DYTree->Branch("GENLepton_fromHardProcessBeforeFSR",&GENLepton_fromHardProcessBeforeFSR,"GENLepton_fromHardProcessBeforeFSR[GENnPair]/I");
		DYTree->Branch("GENLepton_fromHardProcessDecayed",&GENLepton_fromHardProcessDecayed,"GENLepton_fromHardProcessDecayed[GENnPair]/I");
		DYTree->Branch("GENLepton_fromHardProcessFinalState",&GENLepton_fromHardProcessFinalState,"GENLepton_fromHardProcessFinalState[GENnPair]/I");
		DYTree->Branch("GENLepton_isMostlyLikePythia6Status3", &GENLepton_isMostlyLikePythia6Status3, "GENLepton_isMostlyLikePythia6Status3[GENnPair]/I");
		DYTree->Branch("GENEvt_weight",&GENEvt_weight,"GENEvt_weight/D");
	}

	if( theStoreGenOthersFlag )
	{
		// -- GEN Others (GenPhoton ...) info -- //
		DYTree->Branch("nGenOthers",&nGenOthers,"nGenOthers/I");
		DYTree->Branch("GenOthers_phi", &GenOthers_phi,"GenOthers_phi[nGenOthers]/D");
		DYTree->Branch("GenOthers_eta", &GenOthers_eta,"GenOthers_eta[nGenOthers]/D");
		DYTree->Branch("GenOthers_pT", &GenOthers_pT,"GenOthers_pT[nGenOthers]/D");
		DYTree->Branch("GenOthers_Px", &GenOthers_Px,"GenOthers_Px[nGenOthers]/D");
		DYTree->Branch("GenOthers_Py", &GenOthers_Py,"GenOthers_Py[nGenOthers]/D");
		DYTree->Branch("GenOthers_Pz", &GenOthers_Pz,"GenOthers_Pz[nGenOthers]/D");
		DYTree->Branch("GenOthers_mother", &GenOthers_mother,"GenOthers_mother[nGenOthers]/D");
		DYTree->Branch("GenOthers_charge", &GenOthers_charge,"GenOthers_charge[nGenOthers]/I");
		DYTree->Branch("GenOthers_status", &GenOthers_status,"GenOthers_status[nGenOthers]/I");
		DYTree->Branch("GenOthers_ID", &GenOthers_ID,"GenOthers_ID[nGenOthers]/I");
		DYTree->Branch("GenOthers_isPrompt", &GenOthers_isPrompt,"GenOthers_isPrompt[nGenOthers]/I");
		DYTree->Branch("GenOthers_isPromptFinalState", &GenOthers_isPromptFinalState,"GenOthers_isPromptFinalState[nGenOthers]/I");
		DYTree->Branch("GenOthers_isTauDecayProduct", &GenOthers_isTauDecayProduct,"GenOthers_isTauDecayProduct[nGenOthers]/I");
		DYTree->Branch("GenOthers_isPromptTauDecayProduct", &GenOthers_isPromptTauDecayProduct,"GenOthers_isPromptTauDecayProduct[nGenOthers]/I");
		DYTree->Branch("GenOthers_isDirectPromptTauDecayProductFinalState", &GenOthers_isDirectPromptTauDecayProductFinalState,"GenOthers_isDirectPromptTauDecayProductFinalState[nGenOthers]/I");
		DYTree->Branch("GenOthers_isHardProcess",&GenOthers_isHardProcess,"GenOthers_isHardProcess[nGenOthers]/I");
		DYTree->Branch("GenOthers_isLastCopy",&GenOthers_isLastCopy,"GenOthers_isLastCopy[nGenOthers]/I");
		DYTree->Branch("GenOthers_isLastCopyBeforeFSR",&GenOthers_isLastCopyBeforeFSR,"GenOthers_isLastCopyBeforeFSR[nGenOthers]/I");
		DYTree->Branch("GenOthers_isPromptDecayed",&GenOthers_isPromptDecayed,"GenOthers_isPromptDecayed[nGenOthers]/I");
		DYTree->Branch("GenOthers_isDecayedLeptonHadron",&GenOthers_isDecayedLeptonHadron,"GenOthers_isDecayedLeptonHadron[nGenOthers]/I");
		DYTree->Branch("GenOthers_fromHardProcessBeforeFSR",&GenOthers_fromHardProcessBeforeFSR,"GenOthers_fromHardProcessBeforeFSR[nGenOthers]/I");
		DYTree->Branch("GenOthers_fromHardProcessDecayed",&GenOthers_fromHardProcessDecayed,"GenOthers_fromHardProcessDecayed[nGenOthers]/I");
		DYTree->Branch("GenOthers_fromHardProcessFinalState",&GenOthers_fromHardProcessFinalState,"GenOthers_fromHardProcessFinalState[nGenOthers]/I");
		DYTree->Branch("GenOthers_isMostlyLikePythia6Status3", &GenOthers_isMostlyLikePythia6Status3, "GenOthers_isMostlyLikePythia6Status3[nGenOthers]/I");
	}

	if( theStorePhotonFlag )
	{
		// -- Photon Information -- //
		DYTree->Branch("nPhotons",&nPhotons,"nPhotons/I");
		DYTree->Branch("Photon_hasPixelSeed",&Photon_hasPixelSeed,"Photon_hasPixelSeed[nPhotons]/I");
		DYTree->Branch("Photon_pT",&Photon_pT,"Photon_pT[nPhotons]/D");
		DYTree->Branch("Photon_eta",&Photon_eta,"Photon_eta[nPhotons]/D");
		DYTree->Branch("Photon_phi",&Photon_phi,"Photon_phi[nPhotons]/D");
		DYTree->Branch("Photon_etaSC",&Photon_etaSC,"Photon_etaSC[nPhotons]/D");
		DYTree->Branch("Photon_phiSC",&Photon_phiSC,"Photon_phiSC[nPhotons]/D");
		DYTree->Branch("Photon_HoverE",&Photon_HoverE,"Photon_HoverE[nPhotons]/D");
		DYTree->Branch("Photon_Full5x5_SigmaIEtaIEta",&Photon_Full5x5_SigmaIEtaIEta,"Photon_Full5x5_SigmaIEtaIEta[nPhotons]/D");
		DYTree->Branch("Photon_ChIso",&Photon_ChIso,"Photon_ChIso[nPhotons]/D");
		DYTree->Branch("Photon_NhIso",&Photon_NhIso,"Photon_NhIso[nPhotons]/D");
		DYTree->Branch("Photon_PhIso",&Photon_PhIso,"Photon_PhIso[nPhotons]/D");
		DYTree->Branch("Photon_ChIsoWithEA",&Photon_ChIsoWithEA,"Photon_ChIsoWithEA[nPhotons]/D");
		DYTree->Branch("Photon_NhIsoWithEA",&Photon_NhIsoWithEA,"Photon_NhIsoWithEA[nPhotons]/D");
		DYTree->Branch("Photon_PhIsoWithEA",&Photon_PhIsoWithEA,"Photon_PhIsoWithEA[nPhotons]/D");
	}


	// Pile-up Reweight
	DYTree->Branch("nPileUp",&nPileUp,"nPileUp/I");
	DYTree->Branch("pileUpReweightIn",&pileUpReweightIn,"pileUpReweightIn/D");
	DYTree->Branch("pileUpReweight",&pileUpReweight,"pileUpReweight/D");
	DYTree->Branch("pileUpReweightPlus",&pileUpReweightPlus,"pileUpReweightPlus/D");
	DYTree->Branch("pileUpReweightMinus",&pileUpReweightMinus,"pileUpReweightMinus/D");
	DYTree->Branch("pileUpReweightInMuonPhys",&pileUpReweightInMuonPhys,"pileUpReweightInMuonPhys/D");
	DYTree->Branch("pileUpReweightMuonPhys",&pileUpReweightMuonPhys,"pileUpReweightMuonPhys/D");
	DYTree->Branch("pileUpReweightPlusMuonPhys",&pileUpReweightPlusMuonPhys,"pileUpReweightPlusMuonPhys/D");
	DYTree->Branch("pileUpReweightMinusMuonPhys",&pileUpReweightMinusMuonPhys,"pileUpReweightMinusMuonPhys/D");

	if( theStoreTTFlag )
	{
		DYTree->Branch("NTT", &NTT,"NTT/I");
		DYTree->Branch("TTrack_dxy", &TTrack_dxy,"TTrack_dxy[NTT]/D");
		DYTree->Branch("TTrack_dxyErr", &TTrack_dxyErr,"TTrack_dxyErr[NTT]/D");
		DYTree->Branch("TTrack_d0", &TTrack_d0,"TTrack_d0[NTT]/D");
		DYTree->Branch("TTrack_d0Err", &TTrack_d0Err,"TTrack_d0Err[NTT]/D");
		DYTree->Branch("TTrack_dsz", &TTrack_dsz,"TTrack_dsz[NTT]/D");
		DYTree->Branch("TTrack_dszErr", &TTrack_dszErr,"TTrack_dszErr[NTT]/D");
		DYTree->Branch("TTrack_dz", &TTrack_dz,"TTrack_dz[NTT]/D");
		DYTree->Branch("TTrack_dzErr", &TTrack_dzErr,"TTrack_dzErr[NTT]/D");
		DYTree->Branch("TTrack_dxyBS", &TTrack_dxyBS,"TTrack_dxyBS[NTT]/D");
		DYTree->Branch("TTrack_dszBS", &TTrack_dszBS,"TTrack_dszBS[NTT]/D");
		DYTree->Branch("TTrack_dzBS", &TTrack_dzBS,"TTrack_dzBS[NTT]/D");
		DYTree->Branch("TTrack_pT", &TTrack_pT,"TTrack_pT[NTT]/D");
		DYTree->Branch("TTrack_Px", &TTrack_Px,"TTrack_Px[NTT]/D");
		DYTree->Branch("TTrack_Py", &TTrack_Py,"TTrack_Py[NTT]/D");
		DYTree->Branch("TTrack_Pz", &TTrack_Pz,"TTrack_Pz[NTT]/D");
		DYTree->Branch("TTrack_eta", &TTrack_eta,"TTrack_eta[NTT]/D");
		DYTree->Branch("TTrack_phi", &TTrack_phi,"TTrack_phi[NTT]/D");
		DYTree->Branch("TTrack_charge", &TTrack_charge,"TTrack_charge[NTT]/D");
	}

	if( theStoreMETFlag )
	{
		// -- MET -- //
		DYTree->Branch("pfMET_pT", &pfMET_pT,"pfMET_pT/D");
		DYTree->Branch("pfMET_phi", &pfMET_phi,"pfMET_phi/D");
		DYTree->Branch("pfMET_Px", &pfMET_Px,"pfMET_Px/D");
		DYTree->Branch("pfMET_Py", &pfMET_Py,"pfMET_Py/D");
		DYTree->Branch("pfMET_SumEt", &pfMET_SumEt,"pfMET_SumEt/D");

		DYTree->Branch("pfMET_Type1_pT", &pfMET_Type1_pT,"pfMET_Type1_pT/D");
		DYTree->Branch("pfMET_Type1_phi", &pfMET_Type1_phi,"pfMET_Type1_phi/D");
		DYTree->Branch("pfMET_Type1_Px", &pfMET_Type1_Px,"pfMET_Type1_Px/D");
		DYTree->Branch("pfMET_Type1_Py", &pfMET_Type1_Py,"pfMET_Type1_Py/D");
		DYTree->Branch("pfMET_Type1_SumEt", &pfMET_Type1_SumEt,"pfMET_Type1_SumEt/D");
	}
}

void DYntupleMaker::beginRun(const Run & iRun, const EventSetup & iSetup)
{
	const int nTrigName = 66;
	string trigs[nTrigName] = 
	{
		// -- single muon triggers -- //
		"HLT_IsoMu18_v*",
		"HLT_IsoMu20_v*",
		"HLT_IsoMu22_v*",
		"HLT_IsoMu22_eta2p1_v*",
		"HLT_IsoMu24_v*",
		"HLT_IsoMu27_v*",
		"HLT_IsoTkMu18_v*",
		"HLT_IsoTkMu20_v*",
		"HLT_IsoTkMu22_v*",
		"HLT_IsoTkMu22_eta2p1_v*",
		"HLT_IsoTkMu24_v*",
		"HLT_IsoTkMu27_v*",

		"HLT_Mu20_v*",
		"HLT_TkMu20_v*",
		"HLT_Mu24_eta2p1_v*",
		"HLT_TkMu24_eta2p1_v*",
		"HLT_Mu27_v*",
		"HLT_TkMu27_v*",
		"HLT_Mu45_eta2p1_v*",
		"HLT_Mu50_v*",
		"HLT_TkMu50_v*",
		"HLT_Mu55_v*",

		// -- double muon triggers -- //
		"HLT_Mu17_Mu8_v*",
		"HLT_Mu17_Mu8_DZ_v*",
		"HLT_Mu20_Mu10_v*",
		"HLT_Mu20_Mu10_DZ_v*",
		"HLT_Mu17_TkMu8_DZ_v*",
		"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*",
		"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
		"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",
		"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",

		"HLT_Mu27_TkMu8_v*",
		"HLT_Mu30_TkMu11_v*",
		"HLT_Mu40_TkMu11_v*",

		// -- for Electrons -- //
		"HLT_Ele22_eta2p1_WP75_Gsf_v*",
		"HLT_Ele22_eta2p1_WPLoose_Gsf_v*",
		"HLT_Ele23_WP75_Gsf_v*",
		"HLT_Ele23_WPLoose_Gsf_v*",
		"HLT_Ele27_WPLoose_Gsf_v*",
		"HLT_Ele27_eta2p1_WPLoose_v*",
		"HLT_Ele27_eta2p1_WPLoose_Gsf_v*"
		"HLT_Ele27_WPTight_Gsf_v*",
		"HLT_Ele27_eta2p1_WPTight_Gsf_v*",
		"HLT_Ele30_WPTight_Gsf_v*",
		"HLT_Ele30_eta2p1_WPLoose_Gsf_v*",
		"HLT_Ele30_eta2p1_WPTight_Gsf_v*",
		"HLT_Ele30WP60_SC4_Mass55_v*",
		"HLT_Ele30WP60_Ele8_Mass55_v*",
		"HLT_Ele32_WPTight_Gsf_v*",
		"HLT_Ele32_eta2p1_WPLoose_Gsf_v*",
		"HLT_Ele32_eta2p1_WPTight_Gsf_v*",
		"HLT_Ele35_WPLoose_Gsf_v*",
		"HLT_Ele45_WPLoose_Gsf_v*",
		"HLT_Ele105_CaloIdVT_GsfTrkIdT_v*",
		"HLT_Ele250_CaloIdVT_GsfTrkIdT_v*",
		"HLT_Ele300_CaloIdVT_GsfTrkIdT_v*",

		"HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*",
		"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
		"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
		"HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v*",
		"HLT_DoubleEle25_CaloIdL_GsfTrkIdVL_v*",
		"HLT_DoubleEle33_CaloIdL_v*",
		"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*",
		"HLT_DoubleEle33_CaloIdL_MW_v*",
		"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v*",
		"HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL_v*",
	};

	MuonHLT.clear();
	MuonHLTPS.clear();

	for( int i = 0; i < nTrigName; i++ )
		MuonHLT.push_back(trigs[i]);

	int listRemoval[nTrigName] = {-1};
	int ntrigName = MuonHLT.size();
	bool changedConfig;
	if (!hltConfig_.init(iRun, iSetup, processName, changedConfig))
	{
		LogError("HLTMuonVal") << "Initialization of HLTConfigProvider failed!!";
		return;
	}
	else
	{
		std::vector<std::string> triggerNames = hltConfig_.triggerNames();

		// -- print all trigger pathes -- // 
		for( size_t i = 0; i < triggerNames.size(); i++)
			cout << "Trigger Path: " << triggerNames[i] << endl;

		/////////////////////////////////////////////
		// -- iteration for each input triggers -- //
		/////////////////////////////////////////////
		for( int itrigName = 0; itrigName < ntrigName; itrigName++ )
		{
			cout << "\t[" << itrigName << "th Input Trigger = " << MuonHLT[itrigName] << "]" << endl;

			listRemoval[itrigName] = 0;

			// check list of triggers
			//cout << "trigger = " << itrigName << " " << MuonHLT[itrigName] << endl;
			// bool isMatched = false;

			// -- find triggers in HLT configuration matched with a input trigger using wild card -- //
			std::vector<std::vector<std::string>::const_iterator> matches = edm::regexMatch(triggerNames, MuonHLT[itrigName]);

			cout << "\t[# matched trigger in HLT configuration for this trigger: " << matches.size() << endl;

			if( !matches.empty() )
			{
				//////////////////////////////////////////////
				// -- iteration for each matched trigger -- //
				//////////////////////////////////////////////
				BOOST_FOREACH(std::vector<std::string>::const_iterator match, matches)
				{
					cout << "\t\t[matched trigger = " << *match << "]" << endl;

					// -- find modules corresponding to a trigger in HLT configuration -- //
					std::vector<std::string> moduleNames = hltConfig_.moduleLabels( *match );

					// -- find prescale value -- //
					int _preScaleValue = hltConfig_.prescaleValue(0, *match);
					MuonHLTPS.push_back(_preScaleValue);

					// cout << "Filter name: " << trigModuleNames[moduleNames.size()-2] << endl;
					// for( size_t j = 0; j < moduleNames.size(); j++)
					// {
					// 	TString name = moduleNames[j];
					// 	cout << "\t  Fliter Name: "<<moduleNames[j] << endl;
					// }

					int nsize = moduleNames.size();

					if( nsize-2 >= 0 )
					{
						//cout << "module names = " << moduleNames[nsize-2] << " " << moduleNames[nsize-3] << endl;
						trigModuleNames.push_back(moduleNames[nsize-2]);
						//cout << "Filter name: " << trigModuleNames[trigModuleNames.size()-1] << endl;

						if( nsize-3 >= 0 )
						{
							trigModuleNames_preFil.push_back(moduleNames[nsize-3]);
						}
						else
						{
							trigModuleNames_preFil.push_back("");
						}
					}

					break; // -- just take into account the case # mathced trigger = 1 -- //

				} // -- end of BOOST_FOREACH(std::vector<std::string>::const_iterator match, matches) -- //

			} // -- end of if( !matches.empty() ) -- //
			else 
				listRemoval[itrigName] = 1;

			cout << endl;			

		} // -- end of for( int itrigName = 0; itrigName < ntrigName; itrigName++ ): trigger iteration -- //

	} // -- end of else of if (!hltConfig_.init(iRun, iSetup, processName, changedConfig)) -- //

	// -- Remove unavailable triggers -- //
	cout << "Check whether there are unavailable triggers ..." << endl;
	int itmp = 0;
	for( vector<string>::iterator iter = MuonHLT.begin(); iter != MuonHLT.end(); )
	{
		if( listRemoval[itmp] > 0 )
		{
			cout << "\t" << *iter << " is not available in this HLT configuration .. remove it" << endl;
			iter = MuonHLT.erase(iter);
		}
		else
			++iter;

		itmp++;
	}
	ntrigName = MuonHLT.size();

	cout << "# triggers after removing un-available triggers: " << nTrigName << " -> " << ntrigName << endl;

	cout << "\n[Prescales]" << endl;
	for( int i = 0; i < ntrigName; i++ )
		cout << "[" << MuonHLT[i] << "]\t\t" << MuonHLTPS[i] << endl;

	// trigger filters
	for( int itrig = 0; itrig < ntrigName; itrig++ )
		cout << "Filter name: " << itrig << " " << MuonHLT[itrig] << " " << trigModuleNames[itrig] << " " << trigModuleNames_preFil[itrig] << endl;

	cout << "##### End of Begin Run #####" << endl;
}

// ------------ method called once each job just after ending the event loop  ------------ //
void DYntupleMaker::endJob()
{
	std::cout <<"++++++++++++++++++++++++++++++++++++++" << std::endl;
	std::cout <<"analyzed " << nEvt << " events: " << std::endl;
	std::cout <<"++++++++++++++++++++++++++++++++++++++" << std::endl;
}

///////////////////////////////////////////////////////
// -- makes hlt report and fills it to the ntuple -- //
///////////////////////////////////////////////////////
void DYntupleMaker::hltReport(const edm::Event &iEvent)
{
	int ntrigName = MuonHLT.size();

	// -- read the whole HLT trigger lists fired in an event -- //
	bool *trigFired = new bool[ntrigName];
	for( int i = 0; i < ntrigName; i++ ) 
		trigFired[i] = false;

	Handle<TriggerResults> trigResult;
	iEvent.getByToken(TriggerToken, trigResult);

	if( !trigResult.failedToGet() )
	{
		int ntrigs = trigResult->size();
		const edm::TriggerNames trigName = iEvent.triggerNames(*trigResult);

		// cout << "trigger names in trigger result (HLT)" << endl;
		// for(int itrig=0; itrig<(int)trigName.size(); itrig++)
		// 	cout << "trigName = " << trigName.triggerName(itrig) << " " << itrig << endl;

		for( int itrigName = 0; itrigName < ntrigName; itrigName++ )
		{
			std::vector<std::vector<std::string>::const_iterator> matches = edm::regexMatch(trigName.triggerNames(), MuonHLT[itrigName]);
			if( !matches.empty() )
			{
				BOOST_FOREACH(std::vector<std::string>::const_iterator match, matches)
				{
					//cout << "trigger match = " << *match << endl;
					if( trigName.triggerIndex(*match) >= (unsigned int)ntrigs ) continue;
					if( trigResult->accept(trigName.triggerIndex(*match)) ) trigFired[itrigName] = true;
					//cout << "trigger fire = " << trigFired[itrigName] << endl;
				}

			}

		} // -- end of for( int itrigName = 0; itrigName < ntrigName; itrigName++ ) -- //
		
	} // -- end of if( !trigResult.failedToGet() ) -- //

	Handle<TriggerResults> trigResultPAT;
	iEvent.getByToken(TriggerTokenPAT, trigResultPAT);

	if( !trigResultPAT.failedToGet() )
	{
		const edm::TriggerNames trigName = iEvent.triggerNames(*trigResultPAT);

		// cout << "trigger names in trigger result (PAT)" << endl;
		// for(int itrig=0; itrig<(int)trigName.size(); itrig++)
		// 	cout << "trigName = " << trigName.triggerName(itrig) << " " << itrig << endl;

		if( trigResultPAT->accept(trigName.triggerIndex("Flag_badMuons")) ) Flag_badMuons = true;
		if( trigResultPAT->accept(trigName.triggerIndex("Flag_duplicateMuons")) ) Flag_duplicateMuons = true;
		if( trigResultPAT->accept(trigName.triggerIndex("Flag_noBadMuons")) ) Flag_noBadMuons = true;

		cout << "Flag_badMuons: " << Flag_badMuons << endl;
		cout << "Flag_duplicateMuons: " << Flag_duplicateMuons << endl;
		cout << "Flag_noBadMuons: " << Flag_noBadMuons << endl;
		cout << endl;
	}

	
	///////////////
	// -- AOD -- //
	///////////////

	// edm::Handle< trigger::TriggerEvent > triggerObject;
	// iEvent.getByLabel(edm::InputTag("hltTriggerSummaryAOD","",processName), triggerObject);

	// const trigger::TriggerObjectCollection & toc(triggerObject->getObjects());
	// int ntrigTot = 0;
	// //cout << "size filter = " << triggerObject->sizeFilters() << endl;

	// // -- Loop for all filters in TriggerEvent -- //
	// for( size_t k = 0; k < triggerObject->sizeFilters(); ++ k )
	// {
	// 	// -- Get the full name of k-th filter -- //
	// 	std::string fullname = triggerObject->filterTag(k).encode();

	// 	std::string filterName;

	// 	// -- Find ":" in the full name -- //
	// 	size_t m = fullname.find_first_of(':');

	// 	// -- if ":" exists in the full name, takes the name before ":" as the filter name -- //
	// 	if( m != std::string::npos )
	// 		filterName = fullname.substr(0, m);
	// 	else
	// 		filterName = fullname;

	// 	if( &toc != 0 )
	// 	{
	// 		// -- Get the collection of "key" of the k-th filter -- //
	// 		const trigger::Keys & it = triggerObject->filterKeys(k);

	// 		// -- Loop for the keys -- //
	// 		for( trigger::Keys::const_iterator ky = it.begin(); ky !=it.end(); ++ky )
	// 		{
	// 			// -- Get the trigger object corresponding to ky-th key -- //
	// 			double hlt_pt = toc[*ky].pt();
	// 			//double hlt_eta = toc[*ky].eta();
	// 			//cout << "hlt kine = " << hlt_pt << " " << hlt_eta << endl;

	// 			// -- Loop for the triggers that a user inserted in this code -- //
	// 			for( int itf = 0; itf < ntrigName; itf++ )
	// 			{
	// 				string names = "";
	// 				//cout << "filterName = " << k << " " << filterName << " " << trigModuleNames[itf] << " " << trigModuleNames_preFil[itf] << " " << _HLT_trigName.size() << endl;

	// 				// -- Store HLT object information only if trigModuleName is equal to this filter name -- //
	// 				if( filterName == trigModuleNames[itf] )
	// 				{
	// 					names = MuonHLT[itf];
	// 					int _ps = MuonHLTPS[itf];
	// 					_HLT_trigType[ntrigTot] = itf;
	// 					_HLT_trigFired[ntrigTot] = trigFired[itf];
	// 					_HLT_trigPt[ntrigTot] = hlt_pt;
	// 					_HLT_trigEta[ntrigTot] = toc[*ky].eta();
	// 					_HLT_trigPhi[ntrigTot] = toc[*ky].phi();
	// 					_HLT_trigName.push_back(names);
	// 					_HLT_trigPS.push_back(_ps);
	// 					ntrigTot++;
	// 				}

	// 			} // -- end of for( int itf = 0; itf < ntrigName; itf++ ) -- //

	// 		} // -- end of for( int itf = 0; itf < ntrigName; itf++ ) -- // 

	// 	} // -- end of if( &toc != 0 ) -- //

	// } // -- end of for( size_t k = 0; k < triggerObject->sizeFilters(); ++ k )-- //

	// _HLT_ntrig = ntrigTot;

	// // for( int i = 0; i < _HLT_ntrig; i++ )
	// // 	cout << "trig = " << i << " " << _HLT_trigType[i] << " " << _HLT_trigPt[i] << " " << _HLT_trigEta[i] << " " << _HLT_trigPhi[i] << " " << _HLT_trigName[i] << " " << _HLT_trigPS[i] <<  endl;


	///////////////////
	// -- MiniAOD -- //
	///////////////////
	// cout << "// -- HLT Report for MINIAOD is used -- //" << endl;


	edm::Handle< std::vector<pat::TriggerObjectStandAlone> > triggerObject;
	iEvent.getByToken(TriggerObjectToken, triggerObject);

	int ntrigTot = 0;

	if( !trigResult.failedToGet() )
	{
		const edm::TriggerNames names = iEvent.triggerNames(*trigResult);

		// cout << "[# of trigger object in this event: " << (*triggerObject).size() << endl;
		for (pat::TriggerObjectStandAlone obj : *triggerObject)
		{
			obj.unpackPathNames(names);

			// cout << "# Filters: " << obj.filterLabels().size() << endl;
			for( size_t i_filter = 0; i_filter < obj.filterLabels().size(); ++i_filter )
			{
				// -- Get the full name of i-th filter -- //
				std::string fullname = obj.filterLabels()[i_filter];

				std::string filterName;

				// -- Find ":" in the full name -- //
				size_t m = fullname.find_first_of(':');

				// -- if ":" exists in the full name, takes the name before ":" as the filter name -- //
				if( m != std::string::npos )
					filterName = fullname.substr(0, m);
				else
					filterName = fullname;

				// cout << "\t[" << i_filter << "th Filter] FullName = " << fullname << ", FilterName = " << filterName << endl;

				// -- Loop for the triggers that a user inserted in this code -- //
				for( int itf = 0; itf < ntrigName; itf++ )
				{
					// cout << "\t\t[" << itf << "th trigger] Name = " << MuonHLT[itf] << ", trigModuleName = " << trigModuleNames[itf] << endl;
					string name = "";

					// -- Store HLT object information only if trigModuleName is equal to this filter name -- //
					if( filterName == trigModuleNames[itf] )
					{
						// cout << "\t\t\t[Matched]: filterName = " << filterName << ", Trigger Name = " << MuonHLT[itf] << endl;
						name = MuonHLT[itf];
						int _ps = MuonHLTPS[itf];
						_HLT_trigType[ntrigTot] = itf;
						_HLT_trigFired[ntrigTot] = trigFired[itf];
						_HLT_trigPt[ntrigTot] = obj.pt();
						_HLT_trigEta[ntrigTot] = obj.eta();
						_HLT_trigPhi[ntrigTot] = obj.phi();
						_HLT_trigName.push_back(name);
						_HLT_trigPS.push_back(_ps);
						ntrigTot++;
					}

					// cout << endl;

				} // -- end of for( int itf = 0; itf < ntrigName; itf++ ) -- //

				// cout << endl;

			} // -- end of filter iteration -- //
			
		} // -- end of trigger object iteration -- //

	} // -- end of !trigResult.failedToGet() -- //

	_HLT_ntrig = ntrigTot;

}

///////////////////////////////////
// -- Get Primary vertex info -- //
///////////////////////////////////
void DYntupleMaker::fillPrimaryVertex(const edm::Event &iEvent)
{
	edm::Handle<reco::VertexCollection> pvHandle;
	iEvent.getByToken(PrimaryVertexToken, pvHandle);
	const reco::VertexCollection vtx = *(pvHandle.product());

	if( vtx.size() > 2 && theDebugLevel > 0) cout << "Reconstructed "<< vtx.size() << " vertices" << endl;
	if (vtx.size() > 0 )
	{
		PVtrackSize = vtx.front().tracksSize();
		PVchi2 = vtx.front().chi2();
		PVndof = vtx.front().ndof();
		PVnormalizedChi2 = vtx.front().normalizedChi2();
		PVx = vtx.front().x();
		PVy = vtx.front().y();
		PVz = vtx.front().z();
		PVprob = TMath::Prob(PVchi2,(int)PVndof);
	}

}

//////////////////////////////
// -- Get Muons info -- //
//////////////////////////////
void DYntupleMaker::fillMuons(const edm::Event &iEvent, const edm::EventSetup& iSetup)
{
	// -- BeamSpot -- //
	edm::Handle<reco::BeamSpot> beamSpotHandle;
	iEvent.getByToken(BeamSpotToken, beamSpotHandle);
	reco::BeamSpot beamSpot = (*beamSpotHandle);

	// -- Vertex -- //
	math::XYZPoint RefVtx;
	RefVtx.SetXYZ(0, 0, 0);

	edm::Handle<reco::VertexCollection> pvHandle;
	iEvent.getByToken(PrimaryVertexToken, pvHandle);
	const reco::VertexCollection &vertices = *pvHandle.product();
	nVertices = pvHandle->size();
	// -- What is the purpose of below line? -- //
	for(reco::VertexCollection::const_iterator it=vertices.begin() ; it!=vertices.end() ; ++it)
	{
		RefVtx = it->position();
		break;
	}
	const reco::Vertex &vtx = pvHandle->front();

	int nLepton = 0;

	// muons
	ESHandle<MagneticField> B;
	iSetup.get<IdealMagneticFieldRecord>().get(B);

	// -- Call PAT muons -- //
	edm::Handle< std::vector<pat::Muon> > muonHandle;
	iEvent.getByToken(MuonToken, muonHandle);
	using reco::MuonCollection;
	MuonCollection::const_iterator imuon;

	int _nMuon = 0;
	Nmuons = muonHandle->size();
	for( unsigned i = 0; i != muonHandle->size(); i++ )
	{
		// cout << "##### Analyze:Start the loop for the muon #####" << endl;
		const pat::Muon imuon = muonHandle->at(i);

		if( imuon.isStandAloneMuon() ) 	isSTAmuon[_nMuon] = 1;
		if( imuon.isGlobalMuon() ) 		isGLBmuon[_nMuon] = 1;	   
		if( imuon.isTrackerMuon() ) 	isTRKmuon[_nMuon] = 1;	
		if( imuon.isPFMuon() ) 			isPFmuon[_nMuon] = 1;

		if( imuon.isStandAloneMuon() )
		{
			if( imuon.isGlobalMuon() )
			{
				if( imuon.isTrackerMuon() )
					Muon_muonType[_nMuon] = 0; // -- STA + GLB + TRK -- //
				else 
					Muon_muonType[_nMuon] = 1; // -- STA + GLB -- //
			}
			else
			{
				if( imuon.isTrackerMuon() ) 
					Muon_muonType[_nMuon] = 2; // -- STA + TM -- //
				else 
					Muon_muonType[_nMuon] = 3; // -- STA -- //
			}
		} // -- End of isStandAloneMuon()
		else
		{
			if( imuon.isTrackerMuon() ) 
				Muon_muonType[_nMuon] = 4; // -- TM -- //
		}

		if( Muon_muonType[_nMuon] == 3 ) continue; // -- Dosen't store STA Muon (Not reconstructed in GLB)-- //

		// -- bits 0-1-2-3 = DT stations 1-2-3-4 -- //
		// -- bits 4-5-6-7 = CSC stations 1-2-3-4 -- //
		int _segments = 0;

		for( int idet = 1; idet < 4; idet++ )
		{
			// -- DT (1), CSC (2), RPC (3) -- //
			for( int istation = 1; istation < 5; istation++ )
			{
				// -- station 1, 2, 3, 4 -- //
				_segments += imuon.numberOfSegments(istation, idet);
			}
		}
		Muon_nSegments[_nMuon] = _segments;

		// cout << "##### Analyze:Muon Type #####" << endl;


		// -- reco track information -- //
		reco::TrackRef trackerTrack = imuon.innerTrack();
		reco::TrackRef muonTrack    = imuon.outerTrack();
		reco::TrackRef glbTrack     = imuon.globalTrack();
		// reco::TrackRef cktTrack 	= (muon::tevOptimized(imuon, 200, 17., 40., 0.25)).first;

		// cout << "##### Analyze:Muon Tracks #####" << endl;


		// -- Global track information -- //
		if( glbTrack.isNonnull() )
		{
			Muon_chi2dof[_nMuon] = glbTrack->normalizedChi2();
			Muon_nhits[_nMuon] = glbTrack->numberOfValidHits();

			Muon_qoverp[_nMuon] = glbTrack->qoverp();
			Muon_theta[_nMuon] = glbTrack->theta();
			Muon_lambda[_nMuon] = glbTrack->lambda();
			Muon_dxy[_nMuon] = glbTrack->dxy();
			Muon_d0[_nMuon] = glbTrack->d0();
			Muon_dsz[_nMuon] = glbTrack->dsz();
			Muon_dz[_nMuon] = glbTrack->dz();
			Muon_dxyBS[_nMuon] = glbTrack->dxy(beamSpot.position());
			Muon_dszBS[_nMuon] = glbTrack->dsz(beamSpot.position());
			Muon_dzBS[_nMuon] = glbTrack->dz(beamSpot.position());

			Muon_vx[_nMuon] = glbTrack->vx();
			Muon_vy[_nMuon] = glbTrack->vy();
			Muon_vz[_nMuon] = glbTrack->vz();

			const reco::HitPattern & glbhit = glbTrack->hitPattern();
			Muon_muonHits[_nMuon] = glbhit.numberOfValidMuonHits();

			Muon_trackerHitsGLB[_nMuon] = glbhit.numberOfValidTrackerHits();
			Muon_pixelHitsGLB[_nMuon] = glbhit.numberOfValidPixelHits();
			Muon_trackerLayersGLB[_nMuon] = glbhit.trackerLayersWithMeasurement();
			
		} // -- end of if( glbTrack.isNonnull() ) -- //
		else
		{
			if( trackerTrack.isNonnull() )
			{
				Muon_chi2dof[_nMuon] = trackerTrack->normalizedChi2();
				Muon_nhits[_nMuon] = trackerTrack->numberOfValidHits();

				Muon_qoverp[_nMuon] = trackerTrack->qoverp();
				Muon_theta[_nMuon] = trackerTrack->theta();
				Muon_lambda[_nMuon] = trackerTrack->lambda();
				Muon_dxy[_nMuon] = trackerTrack->dxy();
				Muon_d0[_nMuon] = trackerTrack->d0();
				Muon_dsz[_nMuon] = trackerTrack->dsz();
				Muon_dz[_nMuon] = trackerTrack->dz();
				Muon_dxyBS[_nMuon] = trackerTrack->dxy(beamSpot.position());
				Muon_dszBS[_nMuon] = trackerTrack->dsz(beamSpot.position());
				Muon_dzBS[_nMuon] = trackerTrack->dz(beamSpot.position());

				Muon_vx[_nMuon] = trackerTrack->vx();
				Muon_vy[_nMuon] = trackerTrack->vy();
				Muon_vz[_nMuon] = trackerTrack->vz();

				if( muonTrack.isNonnull() )
				{
					const reco::HitPattern & muonhit = muonTrack->hitPattern();
					Muon_muonHits[_nMuon] = muonhit.numberOfValidMuonHits();
				}
				else
					Muon_muonHits[_nMuon] = 0;
			}
		} // -- end of else of if( glbTrack.isNonnull() ) -- //

		if( trackerTrack.isNonnull() )
		{
			const reco::HitPattern & inhit = trackerTrack->hitPattern();

			Muon_trackerHits[_nMuon] = inhit.numberOfValidTrackerHits();
			Muon_pixelHits[_nMuon] = inhit.numberOfValidPixelHits();
			Muon_trackerLayers[_nMuon] = inhit.trackerLayersWithMeasurement();
		}

		if( !pvHandle->empty() && !pvHandle->front().isFake() )
		{
			Muon_dxyVTX[_nMuon] = imuon.muonBestTrack()->dxy(vtx.position());
			Muon_dszVTX[_nMuon] = imuon.muonBestTrack()->dsz(vtx.position());
			Muon_dzVTX[_nMuon] = imuon.muonBestTrack()->dz(vtx.position());

			// Muon_dxycktVTX[_nMuon] = cktTrack->dxy(vtx.position());
			// Muon_dszcktVTX[_nMuon] = cktTrack->dsz(vtx.position());
			// Muon_dzcktVTX[_nMuon] = cktTrack->dz(vtx.position());
		}

		// muon1 kinematics
		// Muon_cktpT[_nMuon] = cktTrack->pt();
		// Muon_cktPx[_nMuon] = cktTrack->px();
		// Muon_cktPy[_nMuon] = cktTrack->py();
		// Muon_cktPz[_nMuon] = cktTrack->pz();
		// Muon_cktpTError[_nMuon] = cktTrack->ptError();

		Muon_pT[_nMuon] = imuon.pt();
		Muon_Px[_nMuon] = imuon.px();
		Muon_Py[_nMuon] = imuon.py();
		Muon_Pz[_nMuon] = imuon.pz();
		Muon_eta[_nMuon] = imuon.eta();
		Muon_phi[_nMuon] = imuon.phi();

		Muon_dB[_nMuon] = imuon.dB();

		// -- Various track informations -- //
		// -- MuonBestTrack -- //
		if( imuon.muonBestTrack().isNonnull() )
		{
			Muon_Best_pT[_nMuon] = imuon.muonBestTrack()->pt();
			Muon_Best_pTError[_nMuon] = imuon.muonBestTrack()->ptError();
			Muon_Best_Px[_nMuon] = imuon.muonBestTrack()->px();
			Muon_Best_Py[_nMuon] = imuon.muonBestTrack()->py();
			Muon_Best_Pz[_nMuon] = imuon.muonBestTrack()->pz();
			Muon_Best_eta[_nMuon] = imuon.muonBestTrack()->eta();
			Muon_Best_phi[_nMuon] = imuon.muonBestTrack()->phi();
		}


		// -- Inner Track -- //
		if( imuon.innerTrack().isNonnull() )
		{
			Muon_Inner_pT[_nMuon] = imuon.innerTrack()->pt();
			Muon_Inner_pTError[_nMuon] = imuon.innerTrack()->ptError();
			Muon_Inner_Px[_nMuon] = imuon.innerTrack()->px();
			Muon_Inner_Py[_nMuon] = imuon.innerTrack()->py();
			Muon_Inner_Pz[_nMuon] = imuon.innerTrack()->pz();
			Muon_Inner_eta[_nMuon] = imuon.innerTrack()->eta();
			Muon_Inner_phi[_nMuon] = imuon.innerTrack()->phi();
		}

		// -- Outer Track -- //
		if( imuon.outerTrack().isNonnull() )
		{
			Muon_Outer_pT[_nMuon] = imuon.outerTrack()->pt();
			Muon_Outer_pTError[_nMuon] = imuon.outerTrack()->ptError();
			Muon_Outer_Px[_nMuon] = imuon.outerTrack()->px();
			Muon_Outer_Py[_nMuon] = imuon.outerTrack()->py();
			Muon_Outer_Pz[_nMuon] = imuon.outerTrack()->pz();
			Muon_Outer_eta[_nMuon] = imuon.outerTrack()->eta();
			Muon_Outer_phi[_nMuon] = imuon.outerTrack()->phi();
		}

		// -- Global Track -- //
 		if( imuon.globalTrack().isNonnull() )
		{
			Muon_GLB_pT[_nMuon] = imuon.globalTrack()->pt();
			Muon_GLB_pTError[_nMuon] = imuon.globalTrack()->ptError();
			Muon_GLB_Px[_nMuon] = imuon.globalTrack()->px();
			Muon_GLB_Py[_nMuon] = imuon.globalTrack()->py();
			Muon_GLB_Pz[_nMuon] = imuon.globalTrack()->pz();
			Muon_GLB_eta[_nMuon] = imuon.globalTrack()->eta();
			Muon_GLB_phi[_nMuon] = imuon.globalTrack()->phi();
		}

		// -- tuneP MuonBestTrack -- //
		if( imuon.tunePMuonBestTrack().isNonnull() )
		{
			Muon_TuneP_pT[_nMuon] = imuon.tunePMuonBestTrack()->pt();
			Muon_TuneP_pTError[_nMuon] = imuon.tunePMuonBestTrack()->ptError();
			Muon_TuneP_Px[_nMuon] = imuon.tunePMuonBestTrack()->px();
			Muon_TuneP_Py[_nMuon] = imuon.tunePMuonBestTrack()->py();
			Muon_TuneP_Pz[_nMuon] = imuon.tunePMuonBestTrack()->pz();
			Muon_TuneP_eta[_nMuon] = imuon.tunePMuonBestTrack()->eta();
			Muon_TuneP_phi[_nMuon] = imuon.tunePMuonBestTrack()->phi();
		}

		//-- ISOLATIONS GO HERE -- //
		// -- detector based -- //
		Muon_trkiso[_nMuon] = imuon.isolationR03().sumPt;
		Muon_hcaliso[_nMuon] = imuon.isolationR03().hadEt;
		Muon_ecaliso[_nMuon] = imuon.isolationR03().emEt;
		Muon_trkisoR05[_nMuon] = imuon.isolationR05().sumPt;
		Muon_hcalisoR05[_nMuon] = imuon.isolationR05().hadEt;
		Muon_ecalisoR05[_nMuon] = imuon.isolationR05().emEt; 

		// -- pf isolation -- // 
		Muon_PfChargedHadronIsoR04[_nMuon] = imuon.pfIsolationR04().sumChargedHadronPt;
		Muon_PfNeutralHadronIsoR04[_nMuon] = imuon.pfIsolationR04().sumNeutralHadronEt;
		Muon_PfGammaIsoR04[_nMuon] = imuon.pfIsolationR04().sumPhotonEt;
		Muon_PFSumPUIsoR04[_nMuon] = imuon.pfIsolationR04().sumPUPt;

		Muon_PfChargedHadronIsoR03[_nMuon] = imuon.pfIsolationR03().sumChargedHadronPt;
		Muon_PfNeutralHadronIsoR03[_nMuon] = imuon.pfIsolationR03().sumNeutralHadronEt;
		Muon_PfGammaIsoR03[_nMuon] = imuon.pfIsolationR03().sumPhotonEt;
		Muon_PFSumPUIsoR03[_nMuon] = imuon.pfIsolationR03().sumPUPt;

		// -- Else -- //
		Muon_charge[_nMuon] = imuon.charge();
		Muon_nChambers[_nMuon] = imuon.numberOfChambers(); // -- # of chambers -- //
		Muon_nMatches[_nMuon] = imuon.numberOfMatchedStations(); // -- # of chambers with matched segments -- //
		Muon_nMatchesRPCLayers[_nMuon] = imuon.numberOfMatchedRPCLayers();
		Muon_stationMask[_nMuon] = imuon.stationMask(); // -- bit map of stations with matched segments -- //

		// filter for high pt Tight muon
		// https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideMuonId#New_Version_recommended
		if( theApplyFilter )
		{
			if( theFilterType == 0 )
			{
				if( (Muon_cktpT[_nMuon] > 30.0 || Muon_pT[_nMuon] > 30.0)
					&& fabs(Muon_eta[_nMuon]) < 2.4 
					&& Muon_dxycktVTX[_nMuon] < 0.4 
					&& Muon_dzcktVTX[_nMuon] < 1.0 
					&& Muon_pixelHits[_nMuon] > 0 
					&& Muon_trackerLayers[_nMuon] > 5 ) 
				{ 
					nLepton++;
				}

				//if( Muon_chi2dof[_nMuon] < 0 || Muon_chi2dof[_nMuon] > 10 ) continue;
				//if( Muon_muonHits[_nMuon] <= 0 ) continue;
				//if( Muon_nMatches[_nMuon] <= 1 ) continue;
				//if( Muon_cktpTError[_nMuon]/Muon_cktpT[_nMuon] >= 0.3 ) continue;
				// no isolation cut applied
			}

			if( theFilterType == 1 )
			{
				if( (Muon_cktpT[_nMuon] > 40.0 ||  Muon_pT[_nMuon] > 40.0)
					&& fabs(Muon_eta[_nMuon]) < 2.4 
					&& Muon_dxycktVTX[_nMuon] < 0.4
					&& Muon_dzcktVTX[_nMuon] < 1.0
					)
				{
					nLepton++;
				}
			}
		} // -- end of if( theApplyFilter ) -- //

		_nMuon++;

		// -- dimuon variables -- //
		for( unsigned j = 0; j != muonHandle->size(); j++ )
		{
			if( i <= j ) continue; // -- prevent double-counting -- //

			const pat::Muon imuon2 = muonHandle->at(j);
			int index_type = -1;

			if( imuon2.isStandAloneMuon() )
			{
				if( imuon2.isGlobalMuon() )
				{
					if( imuon2.isTrackerMuon() ) 
						index_type = 0; // -- STA + GLB + TRK -- //
					else
						index_type = 1; // -- STA + GLB -- //
				}
				else
				{
					if( imuon2.isTrackerMuon() ) 
						index_type = 2; // -- STA + TRK -- //
					else 
						index_type = 3; // -- STA -- //
				}
			}
			else
			{
				if( imuon2.isTrackerMuon() )
					index_type = 4; // -- TRK -- //
			}

			if( index_type == 3 ) continue; // -- Don't check when 2nd muon is STA muon -- //

			// -- vertex variables are calculated using InnerTrack information -- //
			reco::TrackRef InnerTrk = imuon.innerTrack();
			reco::TrackRef InnerTrk2 = imuon2.innerTrack();

			if( InnerTrk.isNonnull() && InnerTrk2.isNonnull() )
			{
				reco::TransientTrack muTransient1(InnerTrk, B.product());
				reco::TransientTrack muTransient2(InnerTrk2, B.product());

				vector<reco::TransientTrack> dimuonTracksTrk;
				dimuonTracksTrk.push_back(muTransient1);
				dimuonTracksTrk.push_back(muTransient2);
				KalmanVertexFitter KalmanFitterTrk(true);
				CachingVertex<5> vertexTrk;
				TransientVertex vtxtmpTrk;
				bool isVertexTrk = true;

				try
				{
					vertexTrk = KalmanFitterTrk.vertex(dimuonTracksTrk);
					vtxtmpTrk = KalmanFitterTrk.vertex(dimuonTracksTrk);
				}
				catch( exception & err )
				{
					isVertexTrk = false;
				}

				if( isVertexTrk && vertexTrk.isValid() )
				{
					// inv. mass refit using the dimuon vtx
					InvariantMassFromVertex imfvTrk;
					static const double muon_mass = 0.1056583;
					const CachingVertex<5>& vtxTrk = vertexTrk;
					Measurement1D new_massTrk = imfvTrk.invariantMass(vtxTrk, muon_mass);

					vtxTrkCkt1Pt.push_back(InnerTrk->pt());
					vtxTrkCkt2Pt.push_back(InnerTrk2->pt());
					vtxTrkChi2.push_back(vtxTrk.totalChiSquared());
					vtxTrkNdof.push_back(vtxTrk.degreesOfFreedom());
					vtxTrkProb.push_back( TMath::Prob(vtxTrk.totalChiSquared(),(int)vtxTrk.degreesOfFreedom()) );
				}

				// cosmic variable
				double cosine = acos( -InnerTrk->momentum().Dot( InnerTrk2->momentum() / InnerTrk->p()/InnerTrk2->p()) );
				CosAngle.push_back(cosine);

			} // -- end of if( InnerTrk.isNonnull() && InnerTrk2.isNonnull() ) -- //

			// --vertex variables are calculated using TuneP information -- //
			reco::TrackRef TunePTrk = imuon.tunePMuonBestTrack();
			reco::TrackRef TunePTrk2 = imuon2.tunePMuonBestTrack();

			if( TunePTrk.isNonnull() && TunePTrk2.isNonnull() )
			{
				reco::TransientTrack muTransient1(TunePTrk, B.product());
				reco::TransientTrack muTransient2(TunePTrk2, B.product());

				vector<reco::TransientTrack> dimuonTracksTrk;
				dimuonTracksTrk.push_back(muTransient1);
				dimuonTracksTrk.push_back(muTransient2);
				KalmanVertexFitter KalmanFitterTrk(true);
				CachingVertex<5> vertexTrk;
				TransientVertex vtxtmpTrk;
				bool isVertexTrk = true;

				try
				{
					vertexTrk = KalmanFitterTrk.vertex(dimuonTracksTrk);
					vtxtmpTrk = KalmanFitterTrk.vertex(dimuonTracksTrk);
				}
				catch( exception & err )
				{
					isVertexTrk = false;
				}

				if( isVertexTrk && vertexTrk.isValid() )
				{
					// inv. mass refit using the dimuon vtx
					InvariantMassFromVertex imfvTrk;
					static const double muon_mass = 0.1056583;
					const CachingVertex<5>& vtxTrk = vertexTrk;
					Measurement1D new_massTrk = imfvTrk.invariantMass(vtxTrk, muon_mass);

					vtxTrk1Pt_TuneP.push_back(TunePTrk->pt());
					vtxTrk2Pt_TuneP.push_back(TunePTrk2->pt());
					vtxTrkChi2_TuneP.push_back(vtxTrk.totalChiSquared());
					vtxTrkNdof_TuneP.push_back(vtxTrk.degreesOfFreedom());
					vtxTrkProb_TuneP.push_back( TMath::Prob(vtxTrk.totalChiSquared(),(int)vtxTrk.degreesOfFreedom()) );
				}

				// cosmic variable
				double cosine = acos( -TunePTrk->momentum().Dot( TunePTrk2->momentum() / TunePTrk->p()/TunePTrk2->p()) );
				CosAngle_TuneP.push_back(cosine);

			} // -- end of if( TunePTrk.isNonnull() && InnerTrk2.isNonnull() ) -- //

		} // -- end of for( unsigned j = 0; j != muonHandle->size(); j++ ): iteration for 2nd muon -- //

	} // -- End of imuon iteration -- //

	nMuon = _nMuon;
}

//////////////////////////////
// -- Get Electrons info -- //
//////////////////////////////
void DYntupleMaker::fillElectrons(const edm::Event &iEvent)
{
	// cout << "##### Start of fillElectrons #####" << endl;
	// -- BeamSpot -- //
	edm::Handle<reco::BeamSpot> beamSpotHandle;
	iEvent.getByToken(BeamSpotToken, beamSpotHandle);
	reco::BeamSpot beamSpot = (*beamSpotHandle);

	// -- Primary vertex -- //
	edm::Handle<reco::VertexCollection> pvHandle;
	iEvent.getByToken(PrimaryVertexToken, pvHandle);

	// -- electron -- //
	edm::Handle< edm::View<reco::GsfElectron> > ElecHandle;
	iEvent.getByToken(ElectronToken, ElecHandle);

	// -- Get rho value -- //
	edm::Handle< double > rhoHandle;
	iEvent.getByToken(RhoToken, rhoHandle);

	edm::Handle< std::vector<reco::Conversion> > conversions;
	iEvent.getByToken(ConversionsToken, conversions);

	edm::Handle<edm::ValueMap<bool> > veto_id_decisions;
	iEvent.getByToken(eleVetoIdMapToken, veto_id_decisions);

	edm::Handle<edm::ValueMap<bool> > loose_id_decisions;
	iEvent.getByToken(eleLooseIdMapToken, loose_id_decisions);

	edm::Handle<edm::ValueMap<bool> > medium_id_decisions;
	iEvent.getByToken(eleMediumIdMapToken, medium_id_decisions);

	edm::Handle<edm::ValueMap<bool> > tight_id_decisions;
	iEvent.getByToken(eleTightIdMapToken, tight_id_decisions);

	edm::Handle<edm::ValueMap<bool> > mva_id_wp80_decisions;
	iEvent.getByToken(eleMVAIdWP80MapToken, mva_id_wp80_decisions);

	edm::Handle<edm::ValueMap<bool> > mva_id_wp90_decisions;
	iEvent.getByToken(eleMVAIdWP90MapToken, mva_id_wp90_decisions);

	edm::Handle<edm::ValueMap<bool> > heep_id_decisions;
	iEvent.getByToken(eleHEEPIdMapToken, heep_id_decisions);

	int _nElectron = 0;
	std::vector< double > _ambGsfTrkPt;
	for(int i=0; i< (int)ElecHandle->size(); i++)
	{
		const auto el = ElecHandle->ptrAt(i);

		Electron_pT[_nElectron] = el->pt();
		Electron_eta[_nElectron] = el->eta();
		Electron_phi[_nElectron] = el->phi();
		Electron_Energy[_nElectron] = el->energy();
		Electron_charge[_nElectron] = el->charge();
		Electron_fbrem[_nElectron] = el->fbrem();
		Electron_ecalDriven[_nElectron] = el->ecalDrivenSeed();

		// -- Information from SuperCluster -- //
		Electron_EnergySC[_nElectron] = el->superCluster()->energy();
		Electron_preEnergySC[_nElectron] = el->superCluster()->preshowerEnergy();
		Electron_rawEnergySC[_nElectron] = el->superCluster()->rawEnergy();
		Electron_etaSC[_nElectron] = el->superCluster()->eta();
		Electron_phiSC[_nElectron] = el->superCluster()->phi();
		double R = sqrt(el->superCluster()->x()*el->superCluster()->x() + el->superCluster()->y()*el->superCluster()->y() +el->superCluster()->z()*el->superCluster()->z());
		double Rt = sqrt(el->superCluster()->x()*el->superCluster()->x() + el->superCluster()->y()*el->superCluster()->y());
		Electron_etSC[_nElectron] = el->superCluster()->energy()*(Rt/R);

		// -- Information from ECAL  -- //
		Electron_sigmaIEtaIEta[_nElectron] = el->full5x5_sigmaIetaIeta();
		Electron_E15[_nElectron] = el->e1x5();
		Electron_E25[_nElectron] = el->e2x5Max();
		Electron_E55[_nElectron] = el->e5x5();
		Electron_HoverE[_nElectron] = el->hcalOverEcal();
		Electron_etaWidth[_nElectron] = el->superCluster()->etaWidth();
		Electron_phiWidth[_nElectron] = el->superCluster()->phiWidth();
		Electron_r9[_nElectron] = el->r9();

		// -- Information from ECAL & Track -- //
		Electron_dEtaIn[_nElectron] = el->deltaEtaSuperClusterTrackAtVtx();
		Electron_dPhiIn[_nElectron] = el->deltaPhiSuperClusterTrackAtVtx();

		// -- |1/E-1/p| = |1/E - EoverPinner/E| is computed below. The if protects against ecalEnergy == inf or zero -- //
		if( el->ecalEnergy() == 0 ) 
			Electron_InvEminusInvP[_nElectron] = 1e30;
		else if(  !std::isfinite( el->ecalEnergy() )  ) 
			Electron_InvEminusInvP[_nElectron] = 1e30;
		else 
			Electron_InvEminusInvP[_nElectron] = fabs( 1.0/el->ecalEnergy() - el->eSuperClusterOverP()/el->ecalEnergy() );


		// -- Isolation -- //
		double pfCharged = el->pfIsolationVariables().sumChargedHadronPt;
		double pfNeutral = el->pfIsolationVariables().sumNeutralHadronEt;
		double pfPhoton = el->pfIsolationVariables().sumPhotonEt;
		double pfChargedFromPU = el->pfIsolationVariables().sumPUPt;
		Electron_chIso03[_nElectron] = pfCharged;
		Electron_nhIso03[_nElectron] = pfNeutral;
		Electron_phIso03[_nElectron] = pfPhoton;
		Electron_ChIso03FromPU[_nElectron] = pfChargedFromPU;
		Electron_RelPFIso_dBeta[_nElectron] = (pfCharged + max<float>( 0.0, pfNeutral + pfPhoton - 0.5 * pfChargedFromPU))/(el->pt());
		
		// The effective areas constants file in the local release or default CMSSW, whichever is found
		edm::FileInPath eaConstantsFile("RecoEgamma/ElectronIdentification/data/PHYS14/effAreaElectrons_cone03_pfNeuHadronsAndPhotons.txt");
		EffectiveAreas effectiveAreas(eaConstantsFile.fullPath());
		float abseta = fabs(el->superCluster()->eta());
		float eA = effectiveAreas.getEffectiveArea(abseta);
		Electron_RelPFIso_Rho[_nElectron] = (pfCharged + max<float>( 0.0, pfNeutral + pfPhoton - *rhoHandle * eA))/(el->pt());

		// cout << "##### fillElectrons: Before elecTrk #####" << endl;

		// -- Track - Impact Parameter, Conversion rejection, Converted -- //
		reco::GsfTrackRef elecTrk = el->gsfTrack();

		Electron_mHits[_nElectron] = elecTrk->numberOfLostHits();
		Electron_dxy[_nElectron] = elecTrk->dxy();
		Electron_dz[_nElectron] = elecTrk->dz();
		Electron_dxyBS[_nElectron] = elecTrk->dxy(beamSpot.position());
		Electron_dzBS[_nElectron] = elecTrk->dz(beamSpot.position());

		if( elecTrk.isNonnull() )
		{
			Electron_gsfpT[_nElectron] = elecTrk->pt();
			Electron_gsfPx[_nElectron] = elecTrk->px();
			Electron_gsfPy[_nElectron] = elecTrk->py();
			Electron_gsfPz[_nElectron] = elecTrk->pz();
			Electron_gsfEta[_nElectron] = elecTrk->eta();
			Electron_gsfPhi[_nElectron] = elecTrk->phi();
			Electron_gsfCharge[_nElectron] = elecTrk->charge();
		}

		if( !pvHandle->empty() && !pvHandle->front().isFake() )
		{
			const reco::Vertex &vtx = pvHandle->front();
			Electron_dxyVTX[_nElectron] = elecTrk->dxy(vtx.position());
			Electron_dzVTX[_nElectron] = elecTrk->dz(vtx.position());
		}

		bool passConvVeto = !ConversionTools::hasMatchedConversion(*el, conversions, beamSpot.position());
		Electron_passConvVeto[_nElectron] = passConvVeto;

		// -- for ID variables -- //
		bool isPassVeto  = (*veto_id_decisions)[el];
		bool isPassLoose  = (*loose_id_decisions)[el];
		bool isPassMedium = (*medium_id_decisions)[el];
		bool isPassTight  = (*tight_id_decisions)[el];
		bool isPassMVA_WP80  = (*mva_id_wp80_decisions)[el];
		bool isPassMVA_WP90  = (*mva_id_wp90_decisions)[el];
		bool isPassHEEP  = (*heep_id_decisions)[el];

		// cout << "isPassVeto: " << isPassVeto << ", isPassLoose: " << isPassLoose << ", isPassMedium: " << isPassMedium << ", isPassTight: " << isPassTight << endl;
		Electron_passVetoID[_nElectron] = isPassVeto;
		Electron_passLooseID[_nElectron] = isPassLoose;
		Electron_passMediumID[_nElectron] = isPassMedium;
		Electron_passTightID[_nElectron] = isPassTight;
		Electron_passMVAID_WP80[_nElectron] = isPassMVA_WP80;
		Electron_passMVAID_WP90[_nElectron] = isPassMVA_WP90;
		Electron_passHEEPID[_nElectron] = isPassHEEP;

		// cout << "##### fillElectrons: Start Dielectron Loop #####" << endl;
		// -- Dielectron variables -- //
		for(int j=0; j<(int)ElecHandle->size(); j++)
		{
			if( i <= j ) continue; // -- prevent double-counting -- //

			const auto el2 = ElecHandle->ptrAt(j);
		
			reco::GsfTrackRef elecTrk = el->gsfTrack();
			reco::GsfTrackRef elecTrk2 = el2->gsfTrack();

			if( elecTrk.isNonnull() && elecTrk2.isNonnull() )
			{
				vector<reco::TransientTrack> dielecTracksTrk;
				dielecTracksTrk.push_back(theTTBuilder->build(elecTrk));
				dielecTracksTrk.push_back(theTTBuilder->build(elecTrk2));
				KalmanVertexFitter KalmanFitterTrk(true);
				CachingVertex<5> vertexTrk;
				TransientVertex vtxtmpTrk;
				bool isVertexTrk = true;
				try
				{
					vertexTrk = KalmanFitterTrk.vertex(dielecTracksTrk);
					vtxtmpTrk = KalmanFitterTrk.vertex(dielecTracksTrk);
				}
				catch( exception & err ) 
				{
					isVertexTrk = false;
				}

				if( isVertexTrk && vertexTrk.isValid() )
				{
					// inv. mass refit using the dielec vtx
					InvariantMassFromVertex imfvTrk;
					static const double elec_mass = 0.000511;
					const CachingVertex<5>& vtxTrk = vertexTrk;
					Measurement1D new_massTrk = imfvTrk.invariantMass(vtxTrk, elec_mass);
					
					vtxTrkDiE1Pt.push_back(elecTrk->pt());
					vtxTrkDiE2Pt.push_back(elecTrk2->pt());
					vtxTrkDiEChi2.push_back(vtxTrk.totalChiSquared());
					vtxTrkDiENdof.push_back(vtxTrk.degreesOfFreedom());
					vtxTrkDiEProb.push_back(TMath::Prob(vtxTrk.totalChiSquared(),(int)vtxTrk.degreesOfFreedom()));
				}
			} // -- end of if( elecTrk.isNonnull() && elecTrk2.isNonnull() ) -- // 

		} // -- end of for(int j=0; j<(int)ElecHandle->size(); j++): 2nd electron iteration -- //


		// cout << "##### fillElectrons: Start gsf track associated electron collector #####" << endl;

		// // -- check gsf track associated electron collector -- //
		// for( GsfTrackRefVector::const_iterator igsf = el->ambiguousGsfTracksBegin(); igsf != el->ambiguousGsfTracksEnd(); ++igsf )
		// {
		// 	if( (*igsf)->pt() > 30. )
		// 		_ambGsfTrkPt.push_back((*igsf)->pt());
		// }

		// std::sort(_ambGsfTrkPt.begin(), _ambGsfTrkPt.end(), std::greater<double>());
		// for( GsfTrackRefVector::const_iterator igsf = el->ambiguousGsfTracksBegin(); igsf != el->ambiguousGsfTracksEnd(); ++igsf )
		// {
		// 	if( (*igsf)->pt() > 30. )
		// 	{
		// 		if( fabs(_ambGsfTrkPt[0]-(*igsf)->pt()) < 0.00001 )
		// 		{
		// 			Electron_ambGsf0Pt[_nElectron] = (*igsf)->pt();
		// 			Electron_ambGsf0Eta[_nElectron] = (*igsf)->eta();
		// 			Electron_ambGsf0Phi[_nElectron] = (*igsf)->phi();
		// 			Electron_ambGsf0Charge[_nElectron] = (*igsf)->charge();
		// 		}

		// 		if( fabs(_ambGsfTrkPt[1]-(*igsf)->pt()) < 0.00001 )
		// 		{
		// 			Electron_ambGsf1Pt[_nElectron] = (*igsf)->pt();
		// 			Electron_ambGsf1Eta[_nElectron] = (*igsf)->eta();
		// 			Electron_ambGsf1Phi[_nElectron] = (*igsf)->phi();
		// 			Electron_ambGsf1Charge[_nElectron] = (*igsf)->charge();
		// 		}

		// 		if( fabs(_ambGsfTrkPt[2]-(*igsf)->pt()) < 0.00001 ) 
		// 		{
		// 			Electron_ambGsf2Pt[_nElectron] = (*igsf)->pt();
		// 			Electron_ambGsf2Eta[_nElectron] = (*igsf)->eta();
		// 			Electron_ambGsf2Phi[_nElectron] = (*igsf)->phi();
		// 			Electron_ambGsf2Charge[_nElectron] = (*igsf)->charge();
		// 		}

		// 		if( fabs(_ambGsfTrkPt[3]-(*igsf)->pt()) < 0.00001 ) 
		// 		{
		// 			Electron_ambGsf3Pt[_nElectron] = (*igsf)->pt();
		// 			Electron_ambGsf3Eta[_nElectron] = (*igsf)->eta();
		// 			Electron_ambGsf3Phi[_nElectron] = (*igsf)->phi();
		// 			Electron_ambGsf3Charge[_nElectron] = (*igsf)->charge();
		// 		}

		// 	} // -- end of if( (*igsf)->pt() > 30. ) -- // 

		// } // -- end of for( GsfTrackRefVector::const_iterator igsf = el->ambiguousGsfTracksBegin(); igsf != el->ambiguousGsfTracksEnd(); ++igsf ) -- //

		_nElectron++;

	} // -- end of for(int i=0; i< (int)ElecHandle->size(); i++): 1st electron iteration -- //

	Nelectrons = _nElectron;

	// cout << "##### End of fillElectrons #####" << endl;

}

////////////////////////
// -- Get GEN info -- //
////////////////////////
void DYntupleMaker::fillLHEInfo(const edm::Event &iEvent)
{
	Handle<LHEEventProduct> lheEventProduct;
	iEvent.getByToken(LHEEventProductToken, lheEventProduct);

	const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
	std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;

	Int_t _nLHEParticle = 0;
	for( size_t idxParticle = 0; idxParticle < lheParticles.size(); ++idxParticle )
	{
		Int_t id = lheEvent.IDUP[idxParticle];

		if( fabs(id) == 13 || fabs(id) == 11 || fabs(id) == 15 )
		{
			Double_t Px = lheParticles[idxParticle][0];
			Double_t Py = lheParticles[idxParticle][1];
			Double_t Pz = lheParticles[idxParticle][2];
			Double_t E = lheParticles[idxParticle][3];
			// Double_t M = lheParticles[idxParticle][4];		
			Int_t status = lheEvent.ISTUP[idxParticle];

			LHELepton_ID[_nLHEParticle] = id;
			LHELepton_status[_nLHEParticle] = status;
			LHELepton_Px[_nLHEParticle] = Px;
			LHELepton_Py[_nLHEParticle] = Py;
			LHELepton_Pz[_nLHEParticle] = Pz;
			LHELepton_E[_nLHEParticle] = E;

			_nLHEParticle++;
		}
	}
	nLHEParticle = _nLHEParticle;
}

////////////////////////
// -- Get GEN info -- //
////////////////////////
void DYntupleMaker::fillGENInfo(const edm::Event &iEvent)
{
	edm::Handle < std::vector<reco::GenParticle> > particles;
	iEvent.getByToken(GenParticleToken, particles);

	int _GennPair = 0;
	for( size_t ipar = 0; ipar < particles->size(); ipar++ )
	{
		// const reco::GenParticle &parCand = (*particles)[ipar];
		reco::GenParticle parCand = (*particles)[ipar];
		if( abs(parCand.pdgId()) == 13 || abs(parCand.pdgId()) == 11 || abs(parCand.pdgId()) == 15 ) // -- Electron, Muon and Tau -- //
		{
			GENLepton_ID[_GennPair] = parCand.pdgId(); 
			GENLepton_pT[_GennPair] = parCand.pt(); 
			GENLepton_Px[_GennPair] = parCand.px();
			GENLepton_Py[_GennPair] = parCand.py();
			GENLepton_Pz[_GennPair] = parCand.pz();
			GENLepton_eta[_GennPair] = parCand.eta();
			GENLepton_phi[_GennPair] = parCand.phi();
			GENLepton_charge[_GennPair] = parCand.charge();
			GENLepton_status[_GennPair] = parCand.status();
			GENLepton_mother[_GennPair] = parCand.mother(0)->pdgId();

			//Flags (Ref: https://indico.cern.ch/event/402279/contribution/5/attachments/805964/1104514/mcaod-Jun17-2015.pdf)
			GENLepton_isPrompt[_GennPair] = parCand.statusFlags().isPrompt(); //not from hadron, muon or tau decay
			GENLepton_isPromptFinalState[_GennPair] = parCand.isPromptFinalState(); //isPrompt && final state (status==1)
			GENLepton_isTauDecayProduct[_GennPair] = parCand.statusFlags().isTauDecayProduct(); //is directly or indirectly from a tau decay
			GENLepton_isPromptTauDecayProduct[_GennPair] = parCand.statusFlags().isPromptTauDecayProduct(); //is directly or indirectly from a tau decay, where the tau did not come from a hadron decay
			GENLepton_isDirectPromptTauDecayProductFinalState[_GennPair] = parCand.isDirectPromptTauDecayProductFinalState(); // is the direct decay product from a tau decay (ie no intermediate hadron), where the tau did not come from a hadron decay && final state
			GENLepton_isHardProcess[_GennPair] = parCand.isHardProcess();
			GENLepton_isLastCopy[_GennPair] = parCand.isLastCopy();
			GENLepton_isLastCopyBeforeFSR[_GennPair] = parCand.isLastCopyBeforeFSR();
			GENLepton_isPromptDecayed[_GennPair] = parCand.isPromptDecayed();
			GENLepton_isDecayedLeptonHadron[_GennPair] = parCand.statusFlags().isDecayedLeptonHadron();
			GENLepton_fromHardProcessBeforeFSR[_GennPair] = parCand.fromHardProcessBeforeFSR();
			GENLepton_fromHardProcessDecayed[_GennPair] = parCand.fromHardProcessDecayed();
			GENLepton_fromHardProcessFinalState[_GennPair] = parCand.fromHardProcessFinalState();
			GENLepton_isMostlyLikePythia6Status3[_GennPair] = parCand.isMostlyLikePythia6Status3();

			_GennPair++;
		} 
	} // -- end of for( size_t ipar = 0; ipar < particles->size(); ipar++ ) -- //

	GENnPair = _GennPair;

	edm::Handle<GenEventInfoProduct> genEvtInfo;
	iEvent.getByToken(GenEventInfoToken, genEvtInfo);
	GENEvt_weight = genEvtInfo->weight();

}

///////////////////////////////
// -- Get GEN Others info -- // 
///////////////////////////////
void DYntupleMaker::fillGenOthersInfo(const edm::Event &iEvent)
{
	edm::Handle < std::vector<reco::GenParticle> > particles;
	iEvent.getByToken(GenParticleToken, particles);

	int _nGenOthers = 0;
	for( size_t ipar = 0; ipar < particles->size(); ipar++ )
	{
		reco::GenParticle parCand = (*particles)[ipar];
		if( abs(parCand.pdgId()) == 22 ) // -- Photons -- //
		{
			GenOthers_ID[_nGenOthers] = parCand.pdgId(); 
			GenOthers_pT[_nGenOthers] = parCand.pt(); 
			GenOthers_Px[_nGenOthers] = parCand.px();
			GenOthers_Py[_nGenOthers] = parCand.py();
			GenOthers_Pz[_nGenOthers] = parCand.pz();
			GenOthers_eta[_nGenOthers] = parCand.eta();
			GenOthers_phi[_nGenOthers] = parCand.phi();
			GenOthers_charge[_nGenOthers] = parCand.charge();
			GenOthers_status[_nGenOthers] = parCand.status();
			GenOthers_mother[_nGenOthers] = parCand.mother(0)->pdgId();

			//Flags (Ref: https://indico.cern.ch/event/402279/contribution/5/attachments/805964/1104514/mcaod-Jun17-2015.pdf)
			GenOthers_isPrompt[_nGenOthers] = parCand.statusFlags().isPrompt(); //not from hadron, muon or tau decay
			GenOthers_isPromptFinalState[_nGenOthers] = parCand.isPromptFinalState(); //isPrompt && final state (status==1)
			GenOthers_isTauDecayProduct[_nGenOthers] = parCand.statusFlags().isTauDecayProduct(); //is directly or indirectly from a tau decay
			GenOthers_isPromptTauDecayProduct[_nGenOthers] = parCand.statusFlags().isPromptTauDecayProduct(); //is directly or indirectly from a tau decay, where the tau did not come from a hadron decay
			GenOthers_isDirectPromptTauDecayProductFinalState[_nGenOthers] = parCand.isDirectPromptTauDecayProductFinalState(); // is the direct decay product from a tau decay (ie no intermediate hadron), where the tau did not come from a hadron decay && final state
			GenOthers_isHardProcess[_nGenOthers] = parCand.isHardProcess();
			GenOthers_isLastCopy[_nGenOthers] = parCand.isLastCopy();
			GenOthers_isLastCopyBeforeFSR[_nGenOthers] = parCand.isLastCopyBeforeFSR();
			GenOthers_isPromptDecayed[_nGenOthers] = parCand.isPromptDecayed();
			GenOthers_isDecayedLeptonHadron[_nGenOthers] = parCand.statusFlags().isDecayedLeptonHadron();
			GenOthers_fromHardProcessBeforeFSR[_nGenOthers] = parCand.fromHardProcessBeforeFSR();
			GenOthers_fromHardProcessDecayed[_nGenOthers] = parCand.fromHardProcessDecayed();
			GenOthers_fromHardProcessFinalState[_nGenOthers] = parCand.fromHardProcessFinalState();
			GenOthers_isMostlyLikePythia6Status3[_nGenOthers] = parCand.isMostlyLikePythia6Status3();

			_nGenOthers++;
		}
	} // -- end of for( size_t ipar = 0; ipar < particles->size(); ipar++ ) -- //

	nGenOthers = _nGenOthers;

}

/////////////////////////
// Get Photons info -- // 
/////////////////////////
void DYntupleMaker::fillPhotons(const edm::Event &iEvent)
{
	edm::Handle< edm::View<reco::Photon> > PhotonHandle;
	iEvent.getByToken(PhotonToken, PhotonHandle);

	// Get rho
	edm::Handle< double > rhoH;
	iEvent.getByToken(RhoToken,rhoH);
	float rho_ = *rhoH;

	// Get the full5x5 map
	edm::Handle<edm::ValueMap<float> > full5x5SigmaIEtaIEtaMap;
	iEvent.getByToken(full5x5SigmaIEtaIEtaMapToken, full5x5SigmaIEtaIEtaMap);

	// Get the isolation maps
	edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
	iEvent.getByToken(phoChargedIsolationToken, phoChargedIsolationMap);

	edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
	iEvent.getByToken(phoNeutralHadronIsolationToken, phoNeutralHadronIsolationMap);

	edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
	iEvent.getByToken(phoPhotonIsolationToken, phoPhotonIsolationMap);

	EffectiveAreas effAreaChHadrons_( effAreaChHadronsFile.fullPath() );
	EffectiveAreas effAreaNeuHadrons_( effAreaNeuHadronsFile.fullPath() );
	EffectiveAreas effAreaPhotons_( effAreaPhotonsFile.fullPath() );


	int _nPhotons = 0;
	for(size_t i=0; i< PhotonHandle->size(); ++i)
	{
		const auto pho = PhotonHandle->ptrAt(i);

		Photon_pT[_nPhotons] = pho->pt();
		Photon_eta[_nPhotons] = pho->eta();
		Photon_phi[_nPhotons] = pho->phi();

		Photon_etaSC[_nPhotons] = pho->superCluster()->eta();
		Photon_phiSC[_nPhotons] = pho->superCluster()->phi();

		Photon_HoverE[_nPhotons] = pho->hadTowOverEm();
		Photon_hasPixelSeed[_nPhotons] = (Int_t)pho->hasPixelSeed();
		Photon_Full5x5_SigmaIEtaIEta[_nPhotons] = (*full5x5SigmaIEtaIEtaMap)[ pho ];

		float chIso = (double)(*phoChargedIsolationMap)[pho];
		float nhIso = (double)(*phoNeutralHadronIsolationMap)[pho];
		float phIso = (double)(*phoPhotonIsolationMap)[pho];
		Photon_ChIso[_nPhotons] = chIso;
		Photon_NhIso[_nPhotons] = nhIso;
		Photon_PhIso[_nPhotons] = phIso;

		float abseta = fabs( pho->superCluster()->eta());
		Photon_ChIsoWithEA[_nPhotons] = std::max( (float)0.0, chIso - rho_*effAreaChHadrons_.getEffectiveArea(abseta) );
		Photon_NhIsoWithEA[_nPhotons] = std::max( (float)0.0, nhIso - rho_*effAreaNeuHadrons_.getEffectiveArea(abseta) );
		Photon_PhIsoWithEA[_nPhotons] = std::max( (float)0.0, phIso - rho_*effAreaPhotons_.getEffectiveArea(abseta) );

		_nPhotons++;
	}

	nPhotons = _nPhotons;
}


/////////////////////////
// -- Get METs info -- // 
/////////////////////////
void DYntupleMaker::fillMET(const edm::Event &iEvent)
{
	edm::Handle< std::vector<pat::MET> > metHandle;
	iEvent.getByToken(MetToken,metHandle);

	if( (metHandle->size() > 1) && (theDebugLevel > 0)) cout << "# of METs = " << metHandle->size() << endl;

	pfMET_pT = metHandle->front().uncorPt();
	pfMET_phi = metHandle->front().uncorPhi();
	pfMET_Px = metHandle->front().uncorPx();
	pfMET_Py = metHandle->front().uncorPy();
	pfMET_SumEt = metHandle->front().uncorSumEt();

	// printf("[pfMET] (pT, phi, Px, Py, sumEt) = (%.3lf, %.3lf, %.3lf, %.3lf, %.3lf)\n", pfMET_pT, pfMET_phi, pfMET_Px, pfMET_Py, pfMET_SumEt);

	pfMET_Type1_pT = metHandle->front().pt();
	pfMET_Type1_phi = metHandle->front().phi();
	pfMET_Type1_Px = metHandle->front().px();
	pfMET_Type1_Py = metHandle->front().py();
	pfMET_Type1_SumEt = metHandle->front().sumEt();


	// pat::METCollection::const_iterator iMET = metHandle->begin();
	// MET_sumEt = iMET->sumEt();
	// MET_pt = iMET->pt();
	// MET_px = iMET->px();
	// MET_py = iMET->py();
	// MET_phi = iMET->phi();

	// edm::Handle<reco::PFMETCollection> pfMETcoll;
	// iEvent.getByLabel(pfMetCollection_, pfMETcoll);
	// if( pfMETcoll.isValid() )
	// {
	// 	const PFMETCollection *pfmetcol = pfMETcoll.product();
	// 	const PFMET *pfmet;
	// 	pfmet = &(pfmetcol->front());
	// 	pfMET_sumEt = pfmet->sumEt();
	// 	pfMET_pt = pfmet->pt();
	// 	pfMET_px = pfmet->px();
	// 	pfMET_py = pfmet->py();
	// 	pfMET_phi = pfmet->phi();
	// }
}

/////////////////////////
// -- Get Jets info -- // 
/////////////////////////
void DYntupleMaker::fillJet(const edm::Event &iEvent)
{
	int _njets = 0;
	// int _nbjets_alg1 = 0;
	// int _nbjets_alg2 = 0;
	// int _nbjets_alg3 = 0;
	// int _nbjets_alg1_closemu = 0;
	// int _nbjets_alg2_closemu = 0;
	// int _nbjets_alg3_closemu = 0;

	// edm::Handle<edm::View<pat::Jet> > jetHandle;
	edm::Handle< std::vector<pat::Jet> > jetHandle;
	iEvent.getByToken(JetToken,jetHandle);

	edm::Handle<std::vector<pat::Muon> > muonHandle;
	iEvent.getByToken(MuonToken,muonHandle);

	if( jetHandle->size() > 0 && theDebugLevel > 0) 
		cout << "# of Jets = " << jetHandle->size() << endl;
	
	Njets = jetHandle->size();
	if(Njets == 0) return;

	for (vector<pat::Jet>::const_iterator jets_iter = jetHandle->begin(); jets_iter != jetHandle->end(); ++jets_iter)
	{
		Jet_pT[_njets] = jets_iter->pt();
		Jet_eta[_njets] = jets_iter->eta();
		Jet_phi[_njets] = jets_iter->phi();
		Jet_Charge[_njets] = jets_iter->jetCharge();
		Jet_Flavor[_njets] = jets_iter->partonFlavour();

		Jet_bTag[_njets] = jets_iter->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags");
		Jet_CHfrac[_njets] = jets_iter->chargedHadronEnergyFraction();
		Jet_NHfrac[_njets] = jets_iter->neutralHadronEnergyFraction();
		Jet_NHEMfrac[_njets] = jets_iter->neutralEmEnergyFraction();
		Jet_CHEMfrac[_njets] = jets_iter->chargedEmEnergyFraction();
		Jet_CHmulti[_njets] = jets_iter->chargedMultiplicity();
		Jet_NHmulti[_njets] = jets_iter->neutralMultiplicity();

		_njets++;
	} 


	// for(edm::View<pat::Jet>::const_iterator iJet = jetHandle->begin(); iJet != jetHandle->end(); ++iJet)
	// {
	// 	JETbDiscriminant_alg1[_njets] = iJet->bDiscriminator("trackCountingHighEffBJetTags");

	// 	JETflavour[_njets] = iJet->partonFlavour();
	// 	JETcharge[_njets] = iJet->jetCharge();
	// 	JETntracks[_njets] = iJet->associatedTracks().size();
	// 	JETpt[_njets] = iJet->pt();
	// 	JETeta[_njets] = iJet->eta();
	// 	JETphi[_njets] = iJet->phi();

	// 	bool _isbtagged = false;
	// 	if( JETbDiscriminant_alg1[_njets] > theBDiscriminant_alg1 )
	// 	{
	// 		_nbjets_alg1++;
	// 		_isbtagged = true;
	// 	}

	// 	// -- checking close muon -- //
	// 	if( _isbtagged )
	// 	{
	// 		for(edm::View<pat::Muon>::const_iterator iMuon = muonHandle->begin(); iMuon != muonHandle->end(); ++iMuon)
	// 		{
	// 			double dR_mu_jet = deltaR(iMuon->eta(), iMuon->phi(), iJet->eta(), iJet->phi());
	// 			if( dR_mu_jet < 0.5 )
	// 			{
	// 				_nbjets_alg1_closemu++;
	// 				break;
	// 			}
	// 		}
	// 	}

	// 	// -- track counting high purity algorithm -- //
	// 	JETbDiscriminant_alg2[_njets] = iJet->bDiscriminator("trackCountingHighPurBJetTags");

	// 	_isbtagged = false;
	// 	if( JETbDiscriminant_alg2[_njets] > theBDiscriminant_alg2 )
	// 	{
	// 		_nbjets_alg2++;
	// 		_isbtagged = true;
	// 	}

	// 	// -- checking close muon -- //
	// 	if( _isbtagged )
	// 	{
	// 		for(edm::View<pat::Muon>::const_iterator iMuon = muonHandle->begin(); iMuon != muonHandle->end(); ++iMuon)
	// 		{
	// 			double dR_mu_jet = deltaR(iMuon->eta(), iMuon->phi(), iJet->eta(), iJet->phi());
	// 			if( dR_mu_jet < 0.5 )
	// 			{
	// 				_nbjets_alg2_closemu++;
	// 				break;
	// 			}
	// 		}
	// 	}

	// 	// -- simple secondary vertex algorithm -- //
	// 	JETbDiscriminant_alg3[_njets] = iJet->bDiscriminator("simpleSecondaryVertexBJetTags");

	// 	_isbtagged = false;
	// 	if( JETbDiscriminant_alg3[_njets] > theBDiscriminant_alg3 )
	// 	{
	// 		_nbjets_alg3++;
	// 		_isbtagged = true;
	// 	}

	// 	// -- checking close muon -- //
	// 	if( _isbtagged )
	// 	{
	// 		for(edm::View<pat::Muon>::const_iterator iMuon = muonHandle->begin(); iMuon != muonHandle->end(); ++iMuon)
	// 		{
	// 			double dR_mu_jet = deltaR(iMuon->eta(), iMuon->phi(), iJet->eta(), iJet->phi());
	// 			if( dR_mu_jet < 0.5 )
	// 			{
	// 				_nbjets_alg3_closemu++;
	// 				break;
	// 			}
	// 		}
	// 	}

	// 	_njets++;
	// } // -- end of for(edm::View<pat::Jet>::const_iterator iJet = jetHandle->begin(); iJet != jetHandle->end(); ++iJet): Jet iteration -- //

	// Nbtagged_alg1 = _nbjets_alg1; 
	// NbtaggedCloseMuon_alg1 = _nbjets_alg1_closemu;
	// Nbtagged_alg2 = _nbjets_alg2;
	// NbtaggedCloseMuon_alg2 = _nbjets_alg2_closemu;
	// Nbtagged_alg3 = _nbjets_alg3;
	// NbtaggedCloseMuon_alg3 = _nbjets_alg3_closemu;

	// cout << "# Jets in this event = " << _njets << endl;

	Njets = _njets;
}

//////////////////////////////////////////////////////
// -- Get Tracker track info (all single tracks) -- //
//////////////////////////////////////////////////////
void DYntupleMaker::fillTT(const edm::Event &iEvent)
{
	edm::Handle<edm::View<reco::Track> > trackHandle;
	iEvent.getByToken(TrackToken, trackHandle);

	edm::Handle<reco::BeamSpot> beamSpotHandle;
	iEvent.getByToken(BeamSpotToken, beamSpotHandle);
	reco::BeamSpot beamSpot = (*beamSpotHandle);

	edm::Handle<reco::VertexCollection> _pvHandle;
	iEvent.getByToken(PrimaryVertexToken, _pvHandle);
	const reco::Vertex &vtx = _pvHandle->front();

	// to store gsftracks close to electrons
	edm::Handle< std::vector< reco::GsfTrack > > gsfTracks; 
	iEvent.getByToken(GsfTrackToken, gsfTracks);

	int _nTT = 0;
	for(unsigned igsf = 0; igsf < gsfTracks->size(); igsf++ )
	{
		GsfTrackRef iTT(gsfTracks, igsf);
		//if( iTT->pt() < 1.0 || iTT->pt() > 100000 ) continue;
		bool _isMatch = false;
		for( int i = 0; i < Nelectrons; i++ )
		{
			double dpT = fabs(iTT->pt() - Electron_gsfpT[i]);
			double dR = deltaR(iTT->eta(), iTT->phi(), Electron_gsfEta[i], Electron_gsfPhi[i]);
			//cout << "elec = " << i << " " << Electron_gsfpT[i] << " " << Electron_gsfEta[i] << " " << Electron_gsfPhi[i] << " " << dR << endl;
			if( dR < 0.001 && dpT < 1.0 ) 
				_isMatch = true;
		}
		if( _isMatch ) continue;

		TTrack_dxy[_nTT] = iTT->dxy(vtx.position());
		TTrack_dxyErr[_nTT] = iTT->dxyError();
		TTrack_d0[_nTT] = iTT->d0();
		TTrack_d0Err[_nTT] = iTT->d0Error(); 
		TTrack_dsz[_nTT] = iTT->dsz(vtx.position());
		TTrack_dszErr[_nTT] = iTT->dszError();
		TTrack_dz[_nTT] = iTT->dz(vtx.position());
		TTrack_dzErr[_nTT] = iTT->dzError();
		TTrack_dxyBS[_nTT] = iTT->dxy(beamSpot.position());
		TTrack_dszBS[_nTT] = iTT->dsz(beamSpot.position());
		TTrack_dzBS[_nTT] = iTT->dz(beamSpot.position());
		TTrack_pT[_nTT] = iTT->pt();
		TTrack_Px[_nTT] = iTT->px();
		TTrack_Py[_nTT] = iTT->py();
		TTrack_Pz[_nTT] = iTT->pz();
		TTrack_eta[_nTT] = iTT->eta();
		TTrack_phi[_nTT] = iTT->phi();
		TTrack_charge[_nTT] = iTT->charge();
		_nTT++;

	} // -- end of for(unsigned igsf = 0; igsf < gsfTracks->size(); igsf++ ): GSFTrack iteration -- //

	NTT = _nTT;
}

//define this as a plug-in
DEFINE_FWK_MODULE(DYntupleMaker);
