import FWCore.ParameterSet.Config as cms

isMC = True

process = cms.Process("DYSkim")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

## Options and Output Report
process.options   = cms.untracked.PSet( 
  wantSummary = cms.untracked.bool(True) 
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

## Source
FileName = ""
if isMC == True:
	FileName = "file:/u/user/kplee/scratch/ROOTFiles_Test/80X/ExampleMiniAODv2_ZMuMuPowheg_M120to200_Moriond17.root"
	# FileName = "file:/u/user/dmpai/scratch/ROOTFiles_Test/80X/00312D7A-FEBD-E611-A713-002590DB923E.root"
else:
  FileName = "file:/cms/home/kplee/scratch/ROOTFiles_Test/80X/SingleMuon_Run2016B_v2_Run273450.root"

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring( FileName )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(3000) )

# -- Geometry and Detector Conditions (needed for a few patTuple production steps) -- #
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# -- Global Tags -- #
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
if isMC == True:
  process.GlobalTag.globaltag = cms.string('80X_mcRun2_asymptotic_2016_TrancheIV_v6')
else:
  process.GlobalTag.globaltag = cms.string('80X_dataRun2_Prompt_v8') #prompt-reco global tag


# -- HLT Filters -- #
# import HLTrigger.HLTfilters.hltHighLevel_cfi
# process.dimuonsHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()

# process.dimuonsHLTFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
# process.dimuonsHLTFilter.HLTPaths = ["HLT_Mu*","HLT_DoubleMu*","HLT_IsoMu*"]

process.TFileService = cms.Service("TFileService",
  fileName = cms.string('ntuple_skim.root')
)

# -- FastFilters -- //
process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
   # src = cms.InputTag("offlinePrimaryVertices"),
   src = cms.InputTag("offlineSlimmedPrimaryVertices"), # -- miniAOD -- #
   cut = cms.string("!isFake && ndof > 4 && abs(z) < 24 && position.Rho < 2"), # tracksSize() > 3 for the older cut
   filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)

# process.noscraping = cms.EDFilter("FilterOutScraping",
#    applyfilter = cms.untracked.bool(True),
#    debugOn = cms.untracked.bool(False),
#    numtrack = cms.untracked.uint32(10),
#    thresh = cms.untracked.double(0.25)
# )

# process.FastFilters = cms.Sequence( process.goodOfflinePrimaryVertices + process.noscraping )
process.FastFilters = cms.Sequence( process.goodOfflinePrimaryVertices )

########################
# -- EGM Correction -- #
########################
# -- EGM 80X regression -- #
from EgammaAnalysis.ElectronTools.regressionWeights_cfi import regressionWeights
process = regressionWeights(process)
process.load('EgammaAnalysis.ElectronTools.regressionApplication_cff')

# -- EGM 80X scale and smearing correction -- #
process.load('Configuration.StandardSequences.Services_cff')
process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",
                  calibratedPatElectrons  = cms.PSet( initialSeed = cms.untracked.uint32(8675389),
                                                      engineName = cms.untracked.string('TRandom3'),
                                                      ),
                  calibratedPatPhotons    = cms.PSet( initialSeed = cms.untracked.uint32(8675389),
                                                      engineName = cms.untracked.string('TRandom3'),
                                                      ),
)

process.load('EgammaAnalysis.ElectronTools.calibratedPatElectronsRun2_cfi')
process.load('EgammaAnalysis.ElectronTools.calibratedPatPhotonsRun2_cfi')

process.calibratedPatElectrons.isMC = cms.bool(isMC)
process.calibratedPatPhotons.isMC = cms.bool(isMC)

#########################
# -- for electron ID -- #
#########################
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff']

process.load("RecoEgamma.ElectronIdentification.ElectronIDValueMapProducer_cfi")
process.electronIDValueMapProducer.srcMiniAOD = cms.InputTag('slimmedElectrons')
process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('slimmedElectrons')

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

###################################
# -- (reco) Photon Information -- #
###################################
# -- photon part should be updated! later when it is necessary -- #

#################
# -- DY Tree -- #
#################
from Phys.DYntupleMaker.DYntupleMaker_cfi import *
from Phys.DYntupleMaker.PUreweight2012_cff import *

process.recoTree = DYntupleMaker.clone()
process.recoTree.isMC = isMC

# -- Objects -- #
process.recoTree.Muon = cms.untracked.InputTag("slimmedMuons") # -- miniAOD -- #
process.recoTree.Electron = cms.untracked.InputTag("slimmedElectrons") # -- miniAOD -- #
process.recoTree.CalibElectron = cms.untracked.InputTag("calibratedPatElectrons") # -- miniAOD -- #
process.recoTree.Photon = cms.untracked.InputTag("slimmedPhotons") # -- miniAOD -- #
process.recoTree.Jet = cms.untracked.InputTag("slimmedJets") # -- miniAOD -- #
process.recoTree.MET = cms.untracked.InputTag("slimmedMETs") # -- miniAOD -- #
process.recoTree.GenParticle = cms.untracked.InputTag("prunedGenParticles") # -- miniAOD -- #

# -- for electrons -- #
process.recoTree.rho = cms.untracked.InputTag("fixedGridRhoFastjetAll")
process.recoTree.conversionsInputTag = cms.untracked.InputTag("reducedEgamma:reducedConversions") # -- miniAOD -- #
process.recoTree.eleVetoIdMap = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-veto")
process.recoTree.eleLooseIdMap = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-loose")
process.recoTree.eleMediumIdMap = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-medium")
process.recoTree.eleTightIdMap = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Summer16-80X-V1-tight")
process.recoTree.eleHEEPIdMap = cms.untracked.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70")
process.recoTree.eleMVAIdWP80Map = cms.untracked.InputTag( "egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80" )
process.recoTree.eleMVAIdWP90Map = cms.untracked.InputTag( "egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90" )

# -- for photons -- #
process.recoTree.full5x5SigmaIEtaIEtaMap   = cms.untracked.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta")
process.recoTree.phoChargedIsolation = cms.untracked.InputTag("photonIDValueMapProducer:phoChargedIsolation")
process.recoTree.phoNeutralHadronIsolation = cms.untracked.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation")
process.recoTree.phoPhotonIsolation = cms.untracked.InputTag("photonIDValueMapProducer:phoPhotonIsolation")
process.recoTree.effAreaChHadFile = cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfChargedHadrons_V2.txt")
process.recoTree.effAreaNeuHadFile= cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfNeutralHadrons_V2.txt")
process.recoTree.effAreaPhoFile   = cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfPhotons_V2.txt")

# -- for Track & Vertex -- #
process.recoTree.PrimaryVertex = cms.untracked.InputTag("offlineSlimmedPrimaryVertices") # -- miniAOD -- #

# -- Else -- #
process.recoTree.PileUpInfo = cms.untracked.InputTag("slimmedAddPileupInfo")

# -- Filters -- #
process.recoTree.ApplyFilter = False

# -- Store Flags -- #
process.recoTree.StoreMuonFlag = True
process.recoTree.StoreElectronFlag = True
process.recoTree.StoreCalibElectronFlag = True
# -- photon part should be updated! later when it is necessary -- #
process.recoTree.StorePhotonFlag = False

process.recoTree.StoreLHEFlag = False
process.recoTree.StoreGENFlag = isMC
process.recoTree.StoreGenOthersFlag = isMC
process.recoTree.StoreJetFlag = True
process.recoTree.StoreMETFlag = True

####################
# -- Let it run -- #
####################
process.p = cms.Path(
  process.FastFilters *
  process.regressionApplication *
  process.calibratedPatElectrons *
  #process.calibratedPatPhotons*
  process.egmGsfElectronIDSequence *
  # process.photonIDValueMapProducer *
  #process.egmPhotonIDSequence*
  process.recoTree
)

# process.p.remove(process.makePatPhotons)
# process.p.remove(process.makePatJets)
# process.p.remove(process.makePatTaus)
# process.p.remove(process.makePatMETs)
# process.p.remove(process.patCandidateSummary)

# if isMC == False:
# 	process.p.remove(process.electronMatch)
# 	process.p.remove(process.muonMatch)
