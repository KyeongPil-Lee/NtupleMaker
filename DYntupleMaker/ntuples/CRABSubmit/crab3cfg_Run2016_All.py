import os
from shutil import copyfile
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = ''
config.General.workArea = 'DYntuple'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../withEGMcorrection/DATA_cfg_ReReco.py'
# config.JobType.maxJobRuntimeMin = 2700 # -- maximum -- #

config.Data.inputDataset = ''

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 40
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False


config.Site.storageSite = 'T3_KR_KISTI'
# config.Site.storageSite = 'T2_KR_KNU'

# config.Data.lumiMask = '/u/user/kplee/JSON/Run2016/Cert_271036-273730_13TeV_PromptReco_Collisions16_JSON.txt'
# config.Data.lumiMask = '/u/user/kplee/JSON/Run2016/Cert_271036-274240_13TeV_PromptReco_Collisions16_JSON.txt'
# config.Data.runRange = '273731-274240'

version = '_v20171022_EGMCorr_'
# 'MultiCRAB' part

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

    # -- Run2016B -- #
    config.General.requestName = 'DYntuple'+version+'SingleMuon_Run2016B_v3_GoldenJSON_%d_to_%d' % (StartRun, EndRun)
    config.Data.inputDataset = '/SingleMuon/Run2016B-03Feb2017_ver2-v2/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016C -- #
    config.General.requestName = 'DYntuple'+version+'SingleMuon_Run2016C_v1_GoldenJSON_%d_to_%d' % (StartRun, EndRun)
    config.Data.inputDataset = '/SingleMuon/Run2016C-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016D -- #
    config.General.requestName = 'DYntuple'+version+'SingleMuon_Run2016D_v1_GoldenJSON_%d_to_%d' % (StartRun, EndRun)
    config.Data.inputDataset = '/SingleMuon/Run2016D-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016E -- #
    config.General.requestName = 'DYntuple'+version+'SingleMuon_Run2016E_v1_GoldenJSON_%d_to_%d' % (StartRun, EndRun)
    config.Data.inputDataset = '/SingleMuon/Run2016E-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016F -- #
    config.General.requestName = 'DYntuple'+version+'SingleMuon_Run2016F_v1_GoldenJSON_%d_to_%d' % (StartRun, EndRun)
    config.Data.inputDataset = '/SingleMuon/Run2016F-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)


    # -- MET phi correction for G to H -- #
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pfMETmultShiftCorrections_GH_cfi.py')
    dst= os.path.join(os.path.dirname(os.path.abspath(__file__)), '..','..','..','..','JetMETCorrections/Type1MET/python/pfMETmultShiftCorrections_cfi.py')
    copyfile(src,dst)

    # -- Run2016G -- #
    config.General.requestName = 'DYntuple'+version+'SingleMuon_Run2016G_v1_GoldenJSON_%d_to_%d' % (StartRun, EndRun)
    config.Data.inputDataset = '/SingleMuon/Run2016G-03Feb2017-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016H, v2 -- #
    config.General.requestName = 'DYntuple'+version+'SingleMuon_Run2016H_v2_GoldenJSON_%d_to_%d' % (StartRun, EndRun)
    config.Data.inputDataset = '/SingleMuon/Run2016H-03Feb2017_ver2-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    config.JobType.psetName = '../withEGMcorrection/DATA_cfg_80X_PromptReco.py'
    crabCommand('submit', config = config)

    # -- Run2016H, v3 -- #
    config.General.requestName = 'DYntuple'+version+'SingleMuon_Run2016H_v3_GoldenJSON_%d_to_%d' % (StartRun, EndRun)
    config.Data.inputDataset = '/SingleMuon/Run2016H-03Feb2017_ver3-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    config.JobType.psetName = '../withEGMcorrection/DATA_cfg_80X_PromptReco.py'
    crabCommand('submit', config = config)

    # config.General.requestName = ''
    # config.Data.inputDataset = ''
    # crabCommand('submit', config = config)