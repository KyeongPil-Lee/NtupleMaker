from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = ''
config.General.workArea = 'DYntuple'

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../MC_cfg_76X_mcRun2_asymptotic_v12.py'

config.Data.inputDataset = ''

config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 2
config.Data.outLFNDirBase = '/store/user/%s/' % (getUsernameFromSiteDB())
config.Data.publication = False

config.Site.storageSite = 'T3_KR_KISTI'
config.JobType.maxJobRuntimeMin = 2160 # -- 36 hours -- #

version = '_v20160525_76X_MINIAODv2_Resubmit_HighMass_'
# 'MultiCRAB' part
if __name__ == '__main__':
    
    from CRABAPI.RawCommand import crabCommand

    # config.General.requestName = 'DYntuple'+version+'ZMuMuPowheg_M50to120_25ns'
    # config.Data.inputDataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_50_120/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    
    # crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'ZMuMuPowheg_M120to200_25ns'
    config.Data.inputDataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_120_200/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'ZMuMuPowheg_M200to400_25ns'
    config.Data.inputDataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_200_400/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'ZMuMuPowheg_M400to800_25ns'
    config.Data.inputDataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_400_800/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'ZMuMuPowheg_M800to1400_25ns'
    config.Data.inputDataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_800_1400/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    # config.General.requestName = 'DYntuple'+version+'ZMuMuPowheg_M1400to2300_25ns'
    # config.Data.inputDataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_1400_2300/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    # crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'ZMuMuPowheg_M2300to3500_25ns'
    config.Data.inputDataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_2300_3500/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'ZMuMuPowheg_M3500to4500_25ns'
    config.Data.inputDataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_3500_4500/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'ZMuMuPowheg_M4500to6000_25ns'
    config.Data.inputDataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_4500_6000/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    config.General.requestName = 'DYntuple'+version+'ZMuMuPowheg_M6000toInf_25ns'
    config.Data.inputDataset = '/ZToMuMu_NNPDF30_13TeV-powheg_M_6000_Inf/RunIIFall15MiniAODv2-PU25nsData2015v1_76X_mcRun2_asymptotic_v12-v1/MINIAODSIM'
    crabCommand('submit', config = config)

    # config.General.requestName = ''
    # config.Data.inputDataset = ''
    # crabCommand('submit', config = config)
