from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = ''
config.General.workArea = 'DYntuple'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../MC_cfg_80X_mcRun2_asymptotic_2016_TrancheIV_v6_SignalMC.py'

config.Data.inputDataset = ''

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
# config.JobType.maxJobRuntimeMin = 2700 # -- 36 hours -- #

config.Site.storageSite = 'T3_KR_KISTI'

version = '_v20170929_'
# 'MultiCRAB' part
if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand

    config.General.requestName = 'DYntuple'+version+'DYLL_M10to50'
    config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'DYLL_M50toInf'
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'DYLL_M100to200'
    config.Data.inputDataset = '/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'DYLL_M200to400'
    config.Data.inputDataset = '/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v2/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'DYLL_M400to500'
    config.Data.inputDataset = '/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'DYLL_M500to700'
    config.Data.inputDataset = '/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'DYLL_M700to800'
    config.Data.inputDataset = '/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'DYLL_M800to1000'
    config.Data.inputDataset = '/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'DYLL_M1000to1500'
    config.Data.inputDataset = '/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'DYLL_M1500to2000'
    config.Data.inputDataset = '/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'DYLL_M2000to3000'
    config.Data.inputDataset = '/DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    # config.General.requestName = ''
    # config.Data.inputDataset = ''
    # crabCommand('submit', config = config)