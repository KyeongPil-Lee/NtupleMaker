from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'DYntuple_v20160303_SingleMuon_RunD_Rereco_MuonPhys'
config.General.workArea = 'DYntuple'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../data_cfg_76X_dataRun2_v15.py'

config.Data.inputDataset = '/SingleMuon/Run2015D-16Dec2015-v1/MINIAOD'

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 20
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False

config.Site.storageSite = 'T2_KR_KNU'

config.Data.lumiMask = '/cms/home/kplee/JSON/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_MuonPhys.txt'
# config.Data.runRange = '260427-260627'

# config.Data.lumiMask = '/cms/home/kplee/JSON/Cert_246908-259891_13TeV_PromptReco_Collisions15_25ns_JSON_MuonPhys.txt'
# config.Data.runRange = '258751-259891'

