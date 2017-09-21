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
	# -- https://twiki.cern.ch/twiki/bin/view/CMS/HEEPElectronIdentificationRun2#Instructions_to_checkout_HEEPV70 -- #
	git cms-merge-topic Sam-Harper:HEEPV70VID_8010_ReducedCheckout  #brings in HEEP V70 into VID
	git cms-merge-topic ikrav:egm_id_80X_v3 #for other E/gamma IDs in VID if you wish to have them
	git cms-merge-topic Sam-Harper:PackedCandNoPuppi #only necessary to run HEEP V70 on AOD (it will crash if this is not present looking for puppi candidates
	mkdir -p ../external/slc6_amd64_gcc530/data/RecoEgamma/ElectronIdentification/ #we need this for the mva weights which runs in VID regardless if you need it or not
	git clone https://github.com/cms-data/RecoEgamma-ElectronIdentification ../external/slc6_amd64_gcc530/data/RecoEgamma/ElectronIdentification/data #we need this for the mva weights which runs in VID regardless if you need it or not

	# -- ntuple maker -- #
	git clone https://github.com/KyeongPil-Lee/NtupleMaker.git Phys -b 80X
	scram b -j 20 >&log&