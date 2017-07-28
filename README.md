# NtupleMaker
* Ntuple maker for physics analysis
   * 80X ntuple (for Z'/DY analysis with 2016 dataset)

* IMPORTANT: Please follow the setup in this wiki page to use Egamma ID variables:
	* https://twiki.cern.ch/twiki/bin/view/CMS/HEEPElectronIdentificationRun2#Instructions_to_checkout_HEEPV70 

## Environment
	export SCRAM_ARCH=slc6_amd64_gcc530
	export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
	source $VO_CMS_SW_DIR/cmsset_default.sh
	cmsrel CMSSW_8_0_26_patch1
	cd CMSSW_8_0_26_patch1/src
	cmsenv
	git clone https://github.com/KyeongPil-Lee/NtupleMaker.git Phys -b 80X_AOD
	scram b -j 20 >&log&