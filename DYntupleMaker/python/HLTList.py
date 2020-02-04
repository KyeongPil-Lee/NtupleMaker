list_HLT = [

# -- single muon triggers -- #
"HLT_IsoMu18_v*",
"HLT_IsoMu20_v*",
"HLT_IsoMu22_v*",
"HLT_IsoMu22_eta2p1_v*",
"HLT_IsoMu24_v*",
"HLT_IsoMu27_v*",
"HLT_IsoTkMu18_v*",
"HLT_IsoTkMu20_v*",
"HLT_IsoTkMu22_v*",
"HLT_IsoTkMu22_eta2p1_v*",
"HLT_IsoTkMu24_v*",
"HLT_IsoTkMu27_v*",

"HLT_Mu20_v*",
"HLT_TkMu20_v*",
"HLT_Mu24_eta2p1_v*",
"HLT_TkMu24_eta2p1_v*",
"HLT_Mu27_v*",
"HLT_TkMu27_v*",
"HLT_Mu45_eta2p1_v*",
"HLT_Mu50_v*",
"HLT_TkMu50_v*",
"HLT_Mu55_v*",

# -- double muon triggers -- #
"HLT_Mu17_Mu8_v*",
"HLT_Mu17_Mu8_DZ_v*",
"HLT_Mu20_Mu10_v*",
"HLT_Mu20_Mu10_DZ_v*",
"HLT_Mu17_TkMu8_DZ_v*",
"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_v*",
"HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v*",
"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_v*",
"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v*",

"HLT_Mu27_TkMu8_v*",
"HLT_Mu30_TkMu11_v*",
"HLT_Mu40_TkMu11_v*",

# -- for Electrons -- #
"HLT_Ele22_eta2p1_WP75_Gsf_v*",
"HLT_Ele22_eta2p1_WPLoose_Gsf_v*",
"HLT_Ele23_WP75_Gsf_v*",
"HLT_Ele23_WPLoose_Gsf_v*",
"HLT_Ele27_WPLoose_Gsf_v*",
"HLT_Ele27_eta2p1_WPLoose_v*",
#"HLT_Ele27_eta2p1_WPLoose_Gsf_v*"
"HLT_Ele27_eta2p1_WPLoose_Gsf_v*",
"HLT_Ele27_WPTight_Gsf_v*",
"HLT_Ele27_eta2p1_WPTight_Gsf_v*",
"HLT_Ele30_WPTight_Gsf_v*",
"HLT_Ele30_eta2p1_WPLoose_Gsf_v*",
"HLT_Ele30_eta2p1_WPTight_Gsf_v*",
"HLT_Ele30WP60_SC4_Mass55_v*",
"HLT_Ele30WP60_Ele8_Mass55_v*",
"HLT_Ele32_WPTight_Gsf_v*",
"HLT_Ele32_eta2p1_WPLoose_Gsf_v*",
"HLT_Ele32_eta2p1_WPTight_Gsf_v*",
"HLT_Ele35_WPLoose_Gsf_v*",
"HLT_Ele45_WPLoose_Gsf_v*",
"HLT_Ele105_CaloIdVT_GsfTrkIdT_v*",
"HLT_Ele250_CaloIdVT_GsfTrkIdT_v*",
"HLT_Ele300_CaloIdVT_GsfTrkIdT_v*",

"HLT_Ele16_Ele12_Ele8_CaloIdL_TrackIdL_v*",
"HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
"HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v*",
"HLT_DoubleEle24_22_eta2p1_WPLoose_Gsf_v*",
"HLT_DoubleEle25_CaloIdL_GsfTrkIdVL_v*",
"HLT_DoubleEle33_CaloIdL_v*",
"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*",
"HLT_DoubleEle33_CaloIdL_MW_v*",
"HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_MW_v*",
"HLT_DoubleEle37_Ele27_CaloIdL_GsfTrkIdVL_v*",
]

# -- for the fake rate measurement in the dielectron channel
# -- method: https://indico.cern.ch/event/622070/contributions/2508978/attachments/1427287/2190433/Preapproval-SMP-17-001-v4.pdf
# -- list of the photon trigger: from the last HLT menu for 2016 data taking
list_HLT_singlePhoton2016 = [
"HLT_Photon22_v*",
"HLT_Photon30_v*",
"HLT_Photon36_v*",
"HLT_Photon50_v*",
"HLT_Photon75_v*",
"HLT_Photon90_v*",
"HLT_Photon120_v*",
"HLT_Photon175_v*",
"HLT_Photon500_v*",
"HLT_Photon600_v*",
]

def GetList_HLT():
    return list_HLT

def GetList_HLT_wSinglePhoton2016():
     return list_HLT + list_HLT_singlePhoton2016

