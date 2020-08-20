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
config.JobType.pyCfgParams = ['globalTag=80X_mcRun2_asymptotic_2016_TrancheIV_v6', 'useSinglePhotonTrigger=1', 'isMC=1', 'isSignalMC=0']
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

    config.General.requestName = 'GJets_Pt20to100_Sherpa'
    config.Data.inputDataset = '/GJets_Pt-20To100_13TeV-sherpa/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'GJets_Pt100to200_Sherpa'
    config.Data.inputDataset = '/GJets_Pt-100To200_13TeV-sherpa/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'GJets_Pt200to500_Sherpa'
    config.Data.inputDataset = '/GJets_Pt-200To500_13TeV-sherpa/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'GJets_Pt500to1000_Sherpa'
    config.Data.inputDataset = '/GJets_Pt-500To1000_13TeV-sherpa/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'GJets_Pt1000to2000_Sherpa'
    config.Data.inputDataset = '/GJets_Pt-1000To2000_13TeV-sherpa/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'GJets_Pt2000to5000_Sherpa'
    config.Data.inputDataset = '/GJets_Pt-2000To5000_13TeV-sherpa/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM'
    crabCommand('submit', config = config)




    # config.General.requestName = ''
    # config.Data.inputDataset = ''
    # crabCommand('submit', config = config)
