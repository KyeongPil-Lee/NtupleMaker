import os
from shutil import copyfile
from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

version = '_v2p0_'

config.General.requestName = ''
config.General.workArea = 'DYntuple'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../withEGMcorrection/MC_cfg_Signal.py'

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

    config.General.requestName = 'DYLL_Pt50to100'
    config.Data.inputDataset = '/DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v3/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_Pt50to100_ext3'
    config.Data.inputDataset = '/DYJetsToLL_Pt-50To100_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext3-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_Pt100to250'
    config.Data.inputDataset = '/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_Pt100to250_ext1'
    config.Data.inputDataset = '/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_Pt100to250_ext2'
    config.Data.inputDataset = '/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_Pt100to250_ext5'
    config.Data.inputDataset = '/DYJetsToLL_Pt-100To250_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext5-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_Pt250to400'
    config.Data.inputDataset = '/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_Pt250to400_ext1'
    config.Data.inputDataset = '/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_Pt250to400_ext2'
    config.Data.inputDataset = '/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_Pt250to400_ext5'
    config.Data.inputDataset = '/DYJetsToLL_Pt-250To400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext5-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_Pt400to650'
    config.Data.inputDataset = '/DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_Pt400to650_ext1'
    config.Data.inputDataset = '/DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_Pt400to650_ext2'
    config.Data.inputDataset = '/DYJetsToLL_Pt-400To650_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_Pt650toInf'
    config.Data.inputDataset = '/DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_Pt650toInf_ext1'
    config.Data.inputDataset = '/DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_Pt650toInf_ext2'
    config.Data.inputDataset = '/DYJetsToLL_Pt-650ToInf_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
    crabCommand('submit', config = config)



    # config.General.requestName = ''
    # config.Data.inputDataset = ''
    # crabCommand('submit', config = config)
