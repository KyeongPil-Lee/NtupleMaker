import os
from shutil import copyfile
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

#version = '_v2p0_'
version = '_v2p2_'

config.General.requestName = ''
config.General.workArea = 'DYntuple'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../withEGMcorrection/MC_cfg_Others.py'

config.Data.inputDataset = ''

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
#config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.outLFNDirBase = '/store/user/%s/%s' % (getUsernameFromSiteDB(), version)
#config.Data.publication = False
config.Data.publication = True
# config.JobType.maxJobRuntimeMin = 2700 # -- 36 hours -- #

config.Site.storageSite = 'T3_KR_KISTI'


### -- 'MultiCRAB' part -- ###

if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand

    # -- MET phi correction for MC -- #
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pfMETmultShiftCorrections_MC_cfi.py')
    dst= os.path.join(os.path.dirname(os.path.abspath(__file__)), '..','..','..','..','JetMETCorrections/Type1MET/python/pfMETmultShiftCorrections_cfi.py')
    copyfile(src,dst)

    config.General.requestName = 'WToENu_M100toInf'
    config.Data.inputDataset = '/WToENu_M-100_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToENu_M200toInf'
    config.Data.inputDataset = '/WToENu_M-200_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToENu_M500toInf'
    config.Data.inputDataset = '/WToENu_M-500_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToENu_M1000toInf'
    config.Data.inputDataset = '/WToENu_M-1000_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToENu_M2000toInf'
    config.Data.inputDataset = '/WToENu_M-2000_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToENu_M3000toInf'
    config.Data.inputDataset = '/WToENu_M-3000_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToENu_M4000toInf'
    config.Data.inputDataset = '/WToENu_M-4000_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToENu_M5000toInf'
    config.Data.inputDataset = '/WToENu_M-5000_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToENu_M6000toInf'
    config.Data.inputDataset = '/WToENu_M-6000_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToMuNu_M100toInf'
    config.Data.inputDataset = '/WToMuNu_M-100_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToMuNu_M200toInf'
    config.Data.inputDataset = '/WToMuNu_M-200_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToMuNu_M500toInf'
    config.Data.inputDataset = '/WToMuNu_M-500_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToMuNu_M1000toInf'
    config.Data.inputDataset = '/WToMuNu_M-1000_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToMuNu_M2000toInf'
    config.Data.inputDataset = '/WToMuNu_M-2000_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToMuNu_M3000toInf'
    config.Data.inputDataset = '/WToMuNu_M-3000_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToMuNu_M4000toInf'
    config.Data.inputDataset = '/WToMuNu_M-4000_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToMuNu_M5000toInf'
    config.Data.inputDataset = '/WToMuNu_M-5000_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToMuNu_M6000toInf'
    config.Data.inputDataset = '/WToMuNu_M-6000_TuneCUETP8M1_13TeV-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToTauNu_M100toInf'
    config.Data.inputDataset = '/WToTauNu_M-100_TuneCUETP8M1_13TeV-pythia8-tauola/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToTauNu_M200toInf'
    config.Data.inputDataset = '/WToTauNu_M-200_TuneCUETP8M1_13TeV-pythia8-tauola/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToTauNu_M500toInf'
    config.Data.inputDataset = '/WToTauNu_M-500_TuneCUETP8M1_13TeV-pythia8-tauola/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToTauNu_M1000toInf'
    config.Data.inputDataset = '/WToTauNu_M-1000_TuneCUETP8M1_13TeV-pythia8-tauola/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToTauNu_M2000toInf'
    config.Data.inputDataset = '/WToTauNu_M-2000_TuneCUETP8M1_13TeV-pythia8-tauola/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToTauNu_M3000toInf'
    config.Data.inputDataset = '/WToTauNu_M-3000_TuneCUETP8M1_13TeV-pythia8-tauola/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToTauNu_M4000toInf'
    config.Data.inputDataset = '/WToTauNu_M-4000_TuneCUETP8M1_13TeV-pythia8-tauola/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToTauNu_M5000toInf'
    config.Data.inputDataset = '/WToTauNu_M-5000_TuneCUETP8M1_13TeV-pythia8-tauola/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'WToTauNu_M6000toInf'
    config.Data.inputDataset = '/WToTauNu_M-6000_TuneCUETP8M1_13TeV-pythia8-tauola/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)



    # config.General.requestName = ''
    # config.Data.inputDataset = ''
    # crabCommand('submit', config = config)
