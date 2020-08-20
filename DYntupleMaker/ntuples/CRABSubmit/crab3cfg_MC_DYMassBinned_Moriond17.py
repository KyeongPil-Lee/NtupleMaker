import os
import argparse

parser = argparse.ArgumentParser(description='CRAB3 configuration to submit ntupler jobs. The version should be given by an argument (it will create a subdirectory under SE with version name)')
parser.add_argument('--version', required=True, help="Production version (e.g. v2p8)")
args = parser.parse_args()
print "version: ", args.version

from shutil import copyfile
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = ''
config.General.workArea = 'DYntuple'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../withEGMcorrection/ntupler_arg.py'
# -- arguments for ntupler_arg.py. inputFile & nEvent are not set here (will automiatically be set by CRAB).
config.JobType.pyCfgParams = ['globalTag=80X_mcRun2_asymptotic_2016_TrancheIV_v6', 'useSinglePhotonTrigger=1', 'isMC=1', 'isSignalMC=1']
config.JobType.inputFiles = ["L1PrefiringMaps_new.root"]

config.Data.inputDataset = ''

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
#config.Data.unitsPerJob = 2
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/kplee/%s' % (args.version)
config.Data.publication = False
# config.JobType.maxJobRuntimeMin = 2700 # -- 36 hours -- #

#config.Site.storageSite = 'T3_KR_KISTI'
config.Site.storageSite = 'T2_KR_KNU'


### -- 'MultiCRAB' part -- ###

if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand

    # -- MET phi correction for MC -- #
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pfMETmultShiftCorrections_MC_cfi.py')
    dst= os.path.join(os.path.dirname(os.path.abspath(__file__)), '..','..','..','..','JetMETCorrections/Type1MET/python/pfMETmultShiftCorrections_cfi.py')
    copyfile(src,dst)

    config.General.requestName = 'DYLL_M10to50_v1'
    config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_M10to50_v2'
    config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_M10to50_ext1v1'
    config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_M50toInf'
    config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_M100to200_fixed'
    config.Data.inputDataset = '/DYJetsToLL_M-100to200_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_M200to400_fixed'
    config.Data.inputDataset = '/DYJetsToLL_M-200to400_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_M400to500'
    config.Data.inputDataset = '/DYJetsToLL_M-400to500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_M500to700'
    config.Data.inputDataset = '/DYJetsToLL_M-500to700_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_M700to800'
    config.Data.inputDataset = '/DYJetsToLL_M-700to800_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_M800to1000'
    config.Data.inputDataset = '/DYJetsToLL_M-800to1000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_M1000to1500'
    config.Data.inputDataset = '/DYJetsToLL_M-1000to1500_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_M1500to2000'
    config.Data.inputDataset = '/DYJetsToLL_M-1500to2000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYLL_M2000to3000'
    config.Data.inputDataset = '/DYJetsToLL_M-2000to3000_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    # config.General.requestName = ''
    # config.Data.inputDataset = ''
    # crabCommand('submit', config = config)
