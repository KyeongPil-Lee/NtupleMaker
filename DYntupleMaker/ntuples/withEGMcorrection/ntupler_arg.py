#########################################################
#-- usage: cmsRun ntupler_arg.py globalTag=<global tag> inputFile=<inputFile> useSinglePhotonTrigger=<1 or 0> nEvent=<number> isMC=<1 or 0> isSignalMC=<1 or 0> -- #
#########################################################

import FWCore.ParameterSet.Config as cms

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')

options.register('globalTag',
                  "80X_mcRun2_asymptotic_2016_TrancheIV_v6", # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.string,         # string, int, or float
                  "Global tag used for the ntuple production")

# -- GT for DATA (03Feb2017): 80X_dataRun2_2016SeptRepro_v7 (Run2016Hv2,v3: 80X_dataRun2_Prompt_v16)
# -- GT for MC (Moriond17): 80X_mcRun2_asymptotic_2016_TrancheIV_v6

options.register('inputFile',
                  "file:/u/user/kplee/scratch/ROOTFiles_Test/80X/MINIAOD_DYLL_M50toInf_Morind17.root", # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.string,         # string, int, or float
                  "input EDM file location (MINIAOD). Don't forget to add 'file:' for the local file ")

# -- example file for DATA: /u/user/kplee/scratch/ROOTFiles_Test/80X/MINIAOD_SingleMuon_Run2016G_03Feb2017.root
# -- example file for MC: /u/user/kplee/scratch/ROOTFiles_Test/80X/MINIAOD_DYLL_M50toInf_Morind17.root

options.register('useSinglePhotonTrigger',
                  1, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.int,         # string, int, or float
                  "Include single photon triggers if true")

options.register('nEvent',
                  -1, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.int,         # string, int, or float
                  "number of events to run")

options.register('isMC',
                  1, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.int,         # string, int, or float
                  "isMC")

options.register('isSignalMC',
                  1, # default value
                  VarParsing.multiplicity.singleton, # singleton or list
                  VarParsing.varType.int,         # string, int, or float
                  "save LHE information related to PDF systematics. It is always false if isMC = false.")

options.parseArguments()

if not options.isMC: options.isSignalMC = 0

print "input global tag           = ", options.globalTag
print "input file                 = ", options.inputFile
print "use single photon triggers = ", options.useSinglePhotonTrigger
print "number of events to run    = ", options.nEvent
print "isMC                       = ", options.isMC
print "isSignalMC                 = ", options.isSignalMC

process = cms.Process("ntupler")

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

## Options and Output Report
process.options   = cms.untracked.PSet( 
  wantSummary = cms.untracked.bool(True),
  SkipEvent = cms.untracked.vstring('ProductNotFound') # -- a few events have no GenEventInfoProduct or LHEInfoProduct: ignore them
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

## Source
process.source = cms.Source("PoolSource",
	fileNames = cms.untracked.vstring( options.inputFile )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.nEvent) )

# -- Geometry and Detector Conditions (needed for a few patTuple production steps) -- #
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# -- Global Tags -- #
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = cms.string(options.globalTag)

# -- HLT Filters -- #
# import HLTrigger.HLTfilters.hltHighLevel_cfi
# process.dimuonsHLTFilter = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone()

# process.dimuonsHLTFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")
# process.dimuonsHLTFilter.HLTPaths = ["HLT_Mu*","HLT_DoubleMu*","HLT_IsoMu*"]

outputName = ""
if options.isMC: outputName = "ntuple_mc.root"
else:            outputName = "ntuple_data.root"

process.TFileService = cms.Service("TFileService",
  fileName = cms.string(outputName)
)

################################
# -- Level 1 ECAL prefiring -- #
################################
# -- https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe -- #
process.prefiringweight = cms.EDProducer("L1ECALPrefiringWeightProducer",
                                         ThePhotons = cms.InputTag("slimmedPhotons"),
                                         TheJets = cms.InputTag("slimmedJets"),
                                         L1Maps = cms.string("L1PrefiringMaps_new.root"),
                                         DataEra = cms.string("2016BtoH"),
                                         UseJetEMPt = cms.bool(False), #can be set to true to use jet prefiring maps parametrized vs pt(em) instead of pt
                                         PrefiringRateSystematicUncty = cms.double(0.2) #Minimum relative prefiring uncty per object
                                         )

#####################
# -- FastFilters -- #
#####################
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
# -- EGM 80X regression: https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMRegression -- #
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

# -- it makes error if the type is not coverted to boolean by hand
process.calibratedPatElectrons.isMC = cms.bool(bool(options.isMC))
process.calibratedPatPhotons.isMC = cms.bool(bool(options.isMC))

#########################
# -- for electron ID -- #
#########################
from PhysicsTools.SelectorUtils.tools.vid_id_tools import *
# turn on VID producer, indicate data format  to be
# DataFormat.AOD or DataFormat.MiniAOD, as appropriate 
dataFormat = DataFormat.MiniAOD
# -- switchOnVIDElectronIdProducer: load RecoEgamma.ElectronIdentification.egmGsfElectronIDs_cff <- makes egmGsfElectronIDSequence
# -- egmGsfElectronIDSequence = cms.Sequence( electronMVAValueMapProducer * egmGsfElectronIDs * electronRegressionValueMapProducer) -- #
switchOnVIDElectronIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.ElectronIdentification.Identification.cutBasedElectronID_Summer16_80X_V1_cff',
                 'RecoEgamma.ElectronIdentification.Identification.mvaElectronID_Spring16_GeneralPurpose_V1_cff',
                 # 'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff'
                 ]

process.load("RecoEgamma.ElectronIdentification.ElectronIDValueMapProducer_cfi")
# process.electronIDValueMapProducer.srcMiniAOD = cms.InputTag('slimmedElectrons')
# process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('slimmedElectrons')

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)


process.selectedElectrons = cms.EDFilter("PATElectronSelector",
    src = cms.InputTag("calibratedPatElectrons"),
    cut = cms.string("pt>5 && abs(eta)")
)

process.egmGsfElectronIDs.physicsObjectSrc = cms.InputTag('selectedElectrons')
process.electronIDValueMapProducer.srcMiniAOD = cms.InputTag('selectedElectrons')
process.electronRegressionValueMapProducer.srcMiniAOD = cms.InputTag('selectedElectrons')
process.electronMVAValueMapProducer.srcMiniAOD = cms.InputTag('selectedElectrons')

###################################
# -- (reco) Photon Information -- #
###################################

switchOnVIDPhotonIdProducer(process, dataFormat)

# define which IDs we want to produce
my_id_modules = ['RecoEgamma.PhotonIdentification.Identification.cutBasedPhotonID_Spring16_V2p2_cff']

#add them to the VID producer
for idmod in my_id_modules:
    setupAllVIDIdsInModule(process,idmod,setupVIDPhotonSelection)


process.selectedPhotons = cms.EDFilter('PATPhotonSelector',
    src = cms.InputTag('calibratedPatPhotons'),
    cut = cms.string('pt>5 && abs(eta)')
)

process.egmPhotonIDs.physicsObjectSrc = cms.InputTag('selectedPhotons')
process.egmPhotonIsolation.srcToIsolate = cms.InputTag('selectedPhotons')
process.photonIDValueMapProducer.srcMiniAOD = cms.InputTag('selectedPhotons')
process.photonRegressionValueMapProducer.srcMiniAOD = cms.InputTag('selectedPhotons')
process.photonMVAValueMapProducer.srcMiniAOD = cms.InputTag('selectedPhotons')

######################
# MET Phi Correction #
######################
isRealData = True
if options.isMC: isRealData = False
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process, isData=isRealData )  #For MC isData=False

#################
# -- DY Tree -- #
#################
from Phys.DYntupleMaker.DYntupleMaker_cfi import *
from Phys.DYntupleMaker.PUreweight2012_cff import *

process.recoTree = DYntupleMaker.clone()
process.recoTree.isMC = bool(options.isMC)

# -- Objects -- #
process.recoTree.Muon = cms.untracked.InputTag("slimmedMuons") # -- miniAOD -- #
process.recoTree.Electron = cms.untracked.InputTag("selectedElectrons") # -- miniAOD -- #
process.recoTree.UnCorrElectron = cms.untracked.InputTag("slimmedElectrons") # -- miniAOD: before applying energy scale correction -- #
#process.recoTree.Photon = cms.untracked.InputTag("slimmedPhotons") # -- miniAOD -- #
process.recoTree.Photon = cms.untracked.InputTag("selectedPhotons") # -- miniAOD -- #
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
# process.recoTree.eleHEEPIdMap = cms.untracked.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70") # -- HEEP recipe is not working under 80X regression recipe (why?): temporarily disabled -- #
process.recoTree.eleMVAIdWP80Map = cms.untracked.InputTag( "egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp80" )
process.recoTree.eleMVAIdWP90Map = cms.untracked.InputTag( "egmGsfElectronIDs:mvaEleID-Spring16-GeneralPurpose-V1-wp90" )

# -- for photons -- #
#process.recoTree.full5x5SigmaIEtaIEtaMap   = cms.untracked.InputTag("photonIDValueMapProducer:phoFull5x5SigmaIEtaIEta")
process.recoTree.phoChargedIsolation = cms.untracked.InputTag("photonIDValueMapProducer:phoChargedIsolation")
process.recoTree.phoNeutralHadronIsolation = cms.untracked.InputTag("photonIDValueMapProducer:phoNeutralHadronIsolation")
process.recoTree.phoPhotonIsolation = cms.untracked.InputTag("photonIDValueMapProducer:phoPhotonIsolation")
process.recoTree.effAreaChHadFile = cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfChargedHadrons_V2.txt")
process.recoTree.effAreaNeuHadFile= cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfNeutralHadrons_V2.txt")
process.recoTree.effAreaPhoFile   = cms.untracked.FileInPath("RecoEgamma/PhotonIdentification/data/PHYS14/effAreaPhotons_cone03_pfPhotons_V2.txt")
process.recoTree.phoMediumIdMap = cms.untracked.InputTag("egmPhotonIDs:cutBasedPhotonID-Spring16-V2p2-medium")

# -- for Track & Vertex -- #
process.recoTree.PrimaryVertex = cms.untracked.InputTag("offlineSlimmedPrimaryVertices") # -- miniAOD -- #

# -- Else -- #
process.recoTree.PileUpInfo = cms.untracked.InputTag("slimmedAddPileupInfo")

# -- Filters -- #
process.recoTree.ApplyFilter = False

# -- Store Flags -- #
process.recoTree.StoreMuonFlag = True
process.recoTree.StoreElectronFlag = True
#process.recoTree.StorePhotonFlag = False # -- photon part should be updated! later when it is necessary -- #
process.recoTree.StorePhotonFlag = True
process.recoTree.StoreJetFlag = True
process.recoTree.StoreMETFlag = True
process.recoTree.StoreGENFlag = bool(options.isMC)
#process.recoTree.StoreGenOthersFlag = isSignalMC
process.recoTree.StoreGenOthersFlag = bool(options.isMC)
process.recoTree.StoreLHEFlag = bool(options.isSignalMC)

from Phys.DYntupleMaker.HLTList import *
if options.useSinglePhotonTrigger:
  process.recoTree.InputHLTList = cms.untracked.vstring(GetList_HLT_wSinglePhoton2016())
else:
  process.recoTree.InputHLTList = cms.untracked.vstring(GetList_HLT())

####################
# -- Let it run -- #
####################
process.p = cms.Path(
  process.prefiringweight * #Level 1 ECAL prefiring
  process.FastFilters *
  process.regressionApplication *
  process.calibratedPatElectrons *
  process.selectedElectrons *
  process.egmGsfElectronIDSequence *
  process.calibratedPatPhotons *
  process.selectedPhotons *
  process.egmPhotonIDSequence *
  process.fullPatMetSequence *  #This is the phi corrections part
  process.recoTree
)
