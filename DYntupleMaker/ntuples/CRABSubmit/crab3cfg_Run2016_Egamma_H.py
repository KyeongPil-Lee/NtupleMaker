# -- Separate configuration for the submission of RunH sample
# ---- The global tag is different with B-F samples
# ---- The change of config.JobType.pyCfgParams in the middle of configuration is not allowed: that's why the configuration is separated
# ---- more info: https://twiki.cern.ch/twiki/bin/view/CMSPublic/CRAB3FAQ#Multiple_submission_fails_with_a

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
config.JobType.pyCfgParams = ['globalTag=80X_dataRun2_Prompt_v16', 'useSinglePhotonTrigger=1', 'isMC=0', 'isSignalMC=0']
config.JobType.inputFiles = ["L1PrefiringMaps_new.root"]

config.Data.inputDataset = ''

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 40
config.Data.outLFNDirBase = '/store/user/kplee/%s' % (args.version)
config.Data.publication = False

#config.Site.storageSite = 'T3_KR_KISTI'
config.Site.storageSite = 'T2_KR_KNU'


### -- 'MultiCRAB' part -- ###

GoldenJSON = './Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
# MuonPhysJSON = '/u/user/kplee/JSON/Run2016/Cert_271036-277148_13TeV_PromptReco_Collisions16_JSON_MuonPhys.txt'
StartRun = 271036
EndRun = 284044

if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand


    # -- MET phi correction for G to H -- #
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'pfMETmultShiftCorrections_GH_cfi.py')
    dst= os.path.join(os.path.dirname(os.path.abspath(__file__)), '..','..','..','..','JetMETCorrections/Type1MET/python/pfMETmultShiftCorrections_cfi.py')
    copyfile(src,dst)

    # -- Run2016H, v2 -- #
    config.General.requestName = 'DoubleEG_Run2016Hver2'
    config.Data.inputDataset = '/DoubleEG/Run2016H-03Feb2017_ver2-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)

    # -- Run2016H, v3 -- #
    config.General.requestName = 'DoubleEG_Run2016Hver3'
    config.Data.inputDataset = '/DoubleEG/Run2016H-03Feb2017_ver3-v1/MINIAOD'
    config.Data.lumiMask = GoldenJSON
    config.Data.runRange = '%d-%d' % (StartRun, EndRun)
    crabCommand('submit', config = config)



    # config.General.requestName = ''
    # config.Data.inputDataset = ''
    # crabCommand('submit', config = config)
