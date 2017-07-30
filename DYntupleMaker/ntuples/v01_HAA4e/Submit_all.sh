#!bin/bash

cd H_200GeV_A_1GeV
cmsRun MC_cfg_80X_mcRun2_asymptotic_2016_AOD.py >&log&
cd ..

cd H_200GeV_A_50GeV
cmsRun MC_cfg_80X_mcRun2_asymptotic_2016_AOD.py >&log&
cd ..

cd H_800GeV_A_1GeV
cmsRun MC_cfg_80X_mcRun2_asymptotic_2016_AOD.py >&log&
cd ..

cd H_800GeV_A_50GeV
cmsRun MC_cfg_80X_mcRun2_asymptotic_2016_AOD.py >&log&
cd ..

cd H_2000GeV_A_1GeV
cmsRun MC_cfg_80X_mcRun2_asymptotic_2016_AOD.py >&log&
cd ..

cd H_2000GeV_A_50GeV
cmsRun MC_cfg_80X_mcRun2_asymptotic_2016_AOD.py >&log&
cd ..

echo "submission is finished"

