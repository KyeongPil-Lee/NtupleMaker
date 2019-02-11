#!/bin/bash


#CRABdir="DYntuple"
CRABdir=$1

## -- Directories to be resubmitted -- ##
Dirs=(
crab_DYLL_M10to50_v1
crab_DYLL_M10to50_ext1v1
crab_WJetsToLNu_amcatnlo
crab_WJetsToLNu_amcatnlo_ext
crab_WJetsToLNu_amcatnlo_ext2v5
crab_ttbar_M1000toInf
crab_ttbar_M700to1000
)


#for ((iter=0;iter<15;iter++)); do
for ((iter=0;iter<20;iter++)); do

    echo "=================== Start  [iteration : $iter] ==================="
    for dir in ${Dirs[@]}; do
        echo "work for $dir..."
	echo
        crab resubmit -d $CRABdir/$dir/
	echo
    done
    echo "=================== End of [iteration : $iter] ==================="

    echo "Sleep 30m"
    sleep 30m
#    echo "Sleep 1h"
#    sleep 1h

    echo "Sleep is finished"

done

