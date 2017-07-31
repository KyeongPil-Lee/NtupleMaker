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
	FileName = "file:/d3/scratch/kplee/ROOTFiles_Test/76X/MiniAOD_DYaMCNLO_M2000to3000.root"
else:
	FileName = "file:/cms/home/kplee/ROOTFiles_Test/ExampleAOD_Run2015Dv4_SingleMuon_Run259891.root"

process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring( FileName )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

# -- Geometry and Detector Conditions (needed for a few patTuple production steps) -- #
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# -- Global Tags -- #
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
if isMC == True:
  process.GlobalTag.globaltag = cms.string('76X_mcRun2_asymptotic_v12')
else:
  process.GlobalTag.globaltag = cms.string('74X_dataRun2_Prompt_v4') #prompt-reco global tag


# -- HLT Filters -- #
import HLTrigger.HLTfilters.hltHighLevel_cfi
process.dimuonsHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()

process.dimuonsHLTFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
process.dimuonsHLTFilter.HLTPaths = ["HLT_Mu*","HLT_DoubleMu*","HLT_IsoMu*"]

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

#########################
# -- for electron ID -- #
#########################
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Spring15_25ns_V1_cff',
                   #'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_PHYS14_PU20bx25_nonTrig_V1_cff',
       'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV60_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

###################################
# -- (reco) Photon Information -- #
###################################
# Several photon variables can not be found inside of a photon object
# and it is easiest to compute them upstream with a dedicated producer,
# such as this standard producer used for photon ID.
#    The producer computes full5x5 cluster shapes and PF isolation variables.
#
# Do not forget to add this producer to the path below!
#
process.load("RecoEgamma/PhotonIdentification/PhotonIDValueMapProducer_cfi")


# -- load the PAT config -- //
# process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
# from PhysicsTools.PatAlgos.tools.coreTools import *
# if isMC==False:
#    removeMCMatching(process, ['All'])

#####################################
# -- FSR weights (for syst.unc.) -- #
#####################################
# Produce event weights to estimate missing QED FSR terms
process.fsrWeight = cms.EDProducer("FSRWeightProducer",
      GenTag = cms.untracked.InputTag("prunedGenParticles"),
)


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
process.recoTree.Photon = cms.untracked.InputTag("slimmedPhotons") # -- miniAOD -- #
process.recoTree.Jet = cms.untracked.InputTag("slimmedJets") # -- miniAOD -- #
process.recoTree.MET = cms.untracked.InputTag("slimmedMETs") # -- miniAOD -- #
process.recoTree.GenParticle = cms.untracked.InputTag("prunedGenParticles") # -- miniAOD -- #
process.recoTree.FSRweight = cms.untracked.InputTag("fsrWeight")

# -- for electrons -- #
process.recoTree.rho = cms.untracked.InputTag("fixedGridRhoFastjetAll")
process.recoTree.conversionsInputTag = cms.untracked.InputTag("reducedEgamma:reducedConversions") # -- miniAOD -- #
process.recoTree.eleVetoIdMap = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-veto")
process.recoTree.eleLooseIdMap = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-loose")
process.recoTree.eleMediumIdMap = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-medium")
process.recoTree.eleTightIdMap = cms.untracked.InputTag("egmGsfElectronIDs:cutBasedElectronID-Spring15-25ns-V1-standalone-tight")

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
process.recoTree.LHEEventProduct = cms.untracked.InputTag("externalLHEProducer")

# -- Filters -- #
process.recoTree.ApplyFilter = False

# -- Store Flags -- #
process.recoTree.StoreElectronFlag = True
process.recoTree.StorePhotonFlag = True
process.recoTree.StoreGENFlag = isMC
process.recoTree.StoreGenOthersFlag = True
process.recoTree.StoreJetFlag = True
process.recoTree.StoreMETFlag = True
process.recoTree.StoreLHEFlag = True

####################
# -- Let it run -- #
####################
process.p = cms.Path(
  process.FastFilters *
  # process.patCandidates *
  process.egmGsfElectronIDSequence *
  process.photonIDValueMapProducer *
  process.fsrWeight *
    # process.patDefaultSequence
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
