# NtupleMaker
* Ntuple maker for physics analysis
   * 76X ntuple (for DY analysis with 2015 dataset)


## Environment
	export SCRAM_ARCH=slc6_amd64_gcc493
	export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
	source $VO_CMS_SW_DIR/cmsset_default.sh
	cmsrel CMSSW_7_6_3_patch2
	cd CMSSW_7_6_3_patch2/src
	git clone https://github.com/KyeongPil-Lee/NtupleMaker.git Phys -b 76X