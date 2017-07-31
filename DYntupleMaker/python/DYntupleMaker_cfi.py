import FWCore.ParameterSet.Config as cms

DYntupleMaker = cms.EDAnalyzer("DYntupleMaker",
	isMC = cms.untracked.bool(True),
	processName = cms.untracked.string("HLT"),
	DebugLevel = cms.untracked.int32(0),

	# -- Object Tags -- #
	Muon = cms.untracked.InputTag("selectedPatMuons"),
	Electron = cms.untracked.InputTag("gedGsfElectrons"),
	Photon = cms.untracked.InputTag("gedPhotons"),
	Jet = cms.untracked.InputTag("selectedPatJets"),
	MET = cms.untracked.InputTag("patMETs"),
	GenParticle = cms.untracked.InputTag("genParticles"),

	# -- electron information -- #
	rho = cms.untracked.InputTag("fixedGridRhoFastjetAll"),
	eleVetoIdMap = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto"),
	eleLooseIdMap = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose"),
	eleMediumIdMap = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium"),
	eleTightIdMap = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight"),
	conversionsInputTag = cms.untracked.InputTag("allConversions"),
	GsfTrack = cms.untracked.InputTag("electronGsfTracks"),

	# -- photon information -- #
	full5x5SigmaIEtaIEtaMap   = cms.untracked.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta"),
	phoChargedIsolation = cms.untracked.InputTag("photonIDValueMapProducer:phoChargedIsolation"),
	phoNeutralHadronIsolation = cms.untracked.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation"),
	phoPhotonIsolation = cms.untracked.InputTag("photonIDValueMapProducer:phoPhotonIsolation"),
	effAreaChHadFile = cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfChargedHadrons_V2.txt"),
	effAreaNeuHadFile= cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfNeutralHadrons_V2.txt"),
	effAreaPhoFile   = cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfPhotons_V2.txt"),

	# -- Jet information -- #
	BDiscriminant_tcheff = cms.untracked.double(0.7),
	BDiscriminant_tchpur = cms.untracked.double(0.7),
	BDiscriminant_ssv = cms.untracked.double(2.05),

	# -- MET information -- #
	pfMET = cms.untracked.InputTag("pfMet"),

	# -- Trigger -- #
	TriggerResults = cms.untracked.InputTag("TriggerResults", "", "HLT"),
	TriggerObject = cms.untracked.InputTag("selectedPatTrigger"),

	# -- Else -- #
	GenEventInfo = cms.untracked.InputTag("generator"),
	BeamSpot = cms.untracked.InputTag("offlineBeamSpot"),
	PrimaryVertex = cms.untracked.InputTag("offlinePrimaryVerticesWithBS"),
	Track = cms.untracked.InputTag("generalTracks"),
	PileUpInfo = cms.untracked.InputTag("addPileupInfo"),

	# -- Store Flags -- #
	StoreElectronFlag = cms.untracked.bool(True),
	StorePhotonFlag = cms.untracked.bool(False),
	StoreJetFlag = cms.untracked.bool(False),
	StoreMETFlag = cms.untracked.bool(False),
	StoreGENFlag = cms.untracked.bool(False),
	StoreGenOthersFlag = cms.untracked.bool(False),
	StorePriVtxFlag = cms.untracked.bool(True),
	StoreTTFlag = cms.untracked.bool(False),
	StoreHLTReportFlag = cms.untracked.bool(True),
	StoreLHEFlag = cms.untracked.bool(True),

	# -- Filters -- #
	ApplyFilter = cms.untracked.bool(False),
	FilterType = cms.untracked.int32(0),

	# -- LHE info -- #
	LHEEventProduct = cms.untracked.InputTag("externalLHEProducer"),
	LHERunInfoProduct = cms.untracked.InputTag("externalLHEProducer"),
)