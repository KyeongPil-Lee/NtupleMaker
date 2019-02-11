#!/bin/bash


BaseLocation="/cms/ldap_home/dmpai/CMSSW_8_0_26_patch1/src/Phys/DYntupleMaker/ntuples/CRABSubmit"

## -- Directories which are complete -- ##
Dirs=(
'crab_DYLL_M1500to2000'
'crab_SingleElectron_RunD'
'crab_SingleElectron_RunG'
'crab_WJetsToLNu_HT600to800_ext'
)

echo "============ Move complete directories to [_finished] ============"

for dir in ${Dirs[@]}; do
    echo "work for $dir..."
    mv $BaseLocation/DYntuple/$dir/ $BaseLocation/DYntuple/_finished/
done

echo "======================== [mv] is finished ========================"


