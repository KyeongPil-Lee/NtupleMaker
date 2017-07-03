from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = ''
config.General.workArea = 'DYntuple'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../data_cfg_80X_dataRun2_Prompt_v16.py'
# config.JobType.maxJobRuntimeMin = 2700 # -- maximum -- #

config.Data.inputDataset = ''

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 40
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False


# config.Site.storageSite = 'T3_KR_KISTI'
config.Site.storageSite = 'T2_KR_KNU'

# config.Data.lumiMask = '/u/user/kplee/JSON/Run2016/Cert_271036-273730_13TeV_PromptReco_Collisions16_JSON.txt'
# config.Data.lumiMask = '/u/user/kplee/JSON/Run2016/Cert_271036-274240_13TeV_PromptReco_Collisions16_JSON.txt'
# config.Data.runRange = '273731-274240'

version = '_v20170428_80XReMiniAOD_FixEvtNum_'
# 'MultiCRAB' part

GoldenJSON = '/u/user/kplee/JSON/Run2016/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
# MuonPhysJSON = '/u/user/kplee/JSON/Run2016/Cert_271036-277148_13TeV_PromptReco_Collisions16_JSON_MuonPhys.txt'
StartRun = 271036
EndRun = 284044

if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand

    # -- Run2016H, v2 -- #
    config.General.requestName = 'DYntuple'+version+'SingleMuon_Run2016H_v2_GoldenJSON_%d_to_%d' % (StartRun, EndRun)
    config.Data.inputDataset = '/SingleMuon/Run2016H-03Feb2017_ver2-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016H, v3 -- #
    config.General.requestName = 'DYntuple'+version+'SingleMuon_Run2016H_v3_GoldenJSON_%d_to_%d' % (StartRun, EndRun)
    config.Data.inputDataset = '/SingleMuon/Run2016H-03Feb2017_ver3-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)


    # config.General.requestName = ''
    # config.Data.inputDataset = ''
    # crabCommand('submit', config = config)