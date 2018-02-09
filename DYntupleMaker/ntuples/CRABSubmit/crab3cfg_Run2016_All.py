import os
from shutil import copyfile
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

#version = '_v2p0_'
version = '_v2p2_'

config.General.requestName = ''
config.General.workArea = 'DYntuple'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../withEGMcorrection/DATA_cfg_ReReco.py'

config.Data.inputDataset = ''

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 40
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.outLFNDirBase = '/store/user/%s/%s' % (getUsernameFromSiteDB(), version)
#config.Data.publication = False
config.Data.publication = True

config.Site.storageSite = 'T3_KR_KISTI'


### -- 'MultiCRAB' part -- ###

GoldenJSON = './Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
# MuonPhysJSON = '/u/user/kplee/JSON/Run2016/Cert_271036-277148_13TeV_PromptReco_Collisions16_JSON_MuonPhys.txt'
StartRun = 271036
EndRun = 284044

if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand

    # -- MET phi correction for B to F -- #
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pfMETmultShiftCorrections_B2F_cfi.py')
    dst= os.path.join(os.path.dirname(os.path.abspath(__file__)), '..','..','..','..','JetMETCorrections/Type1MET/python/pfMETmultShiftCorrections_cfi.py')
    copyfile(src,dst)

### -- SingleMuon -- ###

    # -- Run2016B -- #
    config.General.requestName = 'SingleMuon_Run2016B'
    config.Data.inputDataset = '/SingleMuon/Run2016B-03Feb2017_ver2-v2/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016C -- #
    config.General.requestName = 'SingleMuon_Run2016C'
    config.Data.inputDataset = '/SingleMuon/Run2016C-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016D -- #
    config.General.requestName = 'SingleMuon_Run2016D'
    config.Data.inputDataset = '/SingleMuon/Run2016D-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016E -- #
    config.General.requestName = 'SingleMuon_Run2016E'
    config.Data.inputDataset = '/SingleMuon/Run2016E-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016F -- #
    config.General.requestName = 'SingleMuon_Run2016F'
    config.Data.inputDataset = '/SingleMuon/Run2016F-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

### -- DoubleEG -- ###

    # -- Run2016B -- #
    config.General.requestName = 'DoubleEG_Run2016B'
    config.Data.inputDataset = '/DoubleEG/Run2016B-03Feb2017_ver2-v2/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016C -- #
    config.General.requestName = 'DoubleEG_Run2016C'
    config.Data.inputDataset = '/DoubleEG/Run2016C-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016D -- #
    config.General.requestName = 'DoubleEG_Run2016D'
    config.Data.inputDataset = '/DoubleEG/Run2016D-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016E -- #
    config.General.requestName = 'DoubleEG_Run2016E'
    config.Data.inputDataset = '/DoubleEG/Run2016E-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016F -- #
    config.General.requestName = 'DoubleEG_Run2016F'
    config.Data.inputDataset = '/DoubleEG/Run2016F-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

### -- SingleElectron -- ###

    # -- Run2016B -- #
    config.General.requestName = 'SingleElectron_Run2016B'
    config.Data.inputDataset = '/SingleElectron/Run2016B-03Feb2017_ver2-v2/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016C -- #
    config.General.requestName = 'SingleElectron_Run2016C'
    config.Data.inputDataset = '/SingleElectron/Run2016C-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016D -- #
    config.General.requestName = 'SingleElectron_Run2016D'
    config.Data.inputDataset = '/SingleElectron/Run2016D-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016E -- #
    config.General.requestName = 'SingleElectron_Run2016E'
    config.Data.inputDataset = '/SingleElectron/Run2016E-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016F -- #
    config.General.requestName = 'SingleElectron_Run2016F'
    config.Data.inputDataset = '/SingleElectron/Run2016F-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)


    # -- MET phi correction for G to H -- #
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pfMETmultShiftCorrections_GH_cfi.py')
    dst= os.path.join(os.path.dirname(os.path.abspath(__file__)), '..','..','..','..','JetMETCorrections/Type1MET/python/pfMETmultShiftCorrections_cfi.py')
    copyfile(src,dst)

### -- SingleMuon -- ###

    # -- Run2016G -- #
    config.General.requestName = 'SingleMuon_Run2016G'
    config.Data.inputDataset = '/SingleMuon/Run2016G-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016H, v2 -- #
    config.General.requestName = 'SingleMuon_Run2016Hver2'
    config.Data.inputDataset = '/SingleMuon/Run2016H-03Feb2017_ver2-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    config.JobType.psetName = '../withEGMcorrection/DATA_cfg_80X_PromptReco.py'
    crabCommand('submit', config = config)

    # -- Run2016H, v3 -- #
    config.General.requestName = 'SingleMuon_Run2016Hver3'
    config.Data.inputDataset = '/SingleMuon/Run2016H-03Feb2017_ver3-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    config.JobType.psetName = '../withEGMcorrection/DATA_cfg_80X_PromptReco.py'
    crabCommand('submit', config = config)

### -- DoubleEG -- ###

    # -- Run2016G -- #
    config.General.requestName = 'DoubleEG_Run2016G'
    config.Data.inputDataset = '/DoubleEG/Run2016G-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016H, v2 -- #
    config.General.requestName = 'DoubleEG_Run2016Hver2'
    config.Data.inputDataset = '/DoubleEG/Run2016H-03Feb2017_ver2-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    config.JobType.psetName = '../withEGMcorrection/DATA_cfg_80X_PromptReco.py'
    crabCommand('submit', config = config)

    # -- Run2016H, v3 -- #
    config.General.requestName = 'DoubleEG_Run2016Hver3'
    config.Data.inputDataset = '/DoubleEG/Run2016H-03Feb2017_ver3-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    config.JobType.psetName = '../withEGMcorrection/DATA_cfg_80X_PromptReco.py'
    crabCommand('submit', config = config)

### -- SingleElectron -- ###

    # -- Run2016G -- #
    config.General.requestName = 'SingleElectron_Run2016G'
    config.Data.inputDataset = '/SingleElectron/Run2016G-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016H, v2 -- #
    config.General.requestName = 'SingleElectron_Run2016Hver2'
    config.Data.inputDataset = '/SingleElectron/Run2016H-03Feb2017_ver2-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    config.JobType.psetName = '../withEGMcorrection/DATA_cfg_80X_PromptReco.py'
    crabCommand('submit', config = config)

    # -- Run2016H, v3 -- #
    config.General.requestName = 'SingleElectron_Run2016Hver3'
    config.Data.inputDataset = '/SingleElectron/Run2016H-03Feb2017_ver3-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    config.JobType.psetName = '../withEGMcorrection/DATA_cfg_80X_PromptReco.py'
    crabCommand('submit', config = config)



    # config.General.requestName = ''
    # config.Data.inputDataset = ''
    # crabCommand('submit', config = config)
