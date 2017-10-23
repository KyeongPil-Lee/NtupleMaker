# NtupleMaker
* Ntuple maker for physics analysis
   * 80X ntuple (for Z'/DY analysis with 2016 dataset)

## Environment
	export SCRAM_ARCH=slc6_amd64_gcc530
	export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
	source $VO_CMS_SW_DIR/cmsset_default.sh
	cmsrel CMSSW_8_0_26_patch1
	cd CMSSW_8_0_26_patch1/src
	cmsenv

	# -- EGM corrections -- # (https://twiki.cern.ch/twiki/bin/viewauth/CMS/EGMRegression)
	git cms-init
	git cms-merge-topic cms-egamma:EGM_gain_v1
	cd EgammaAnalysis/ElectronTools/data
	git clone -b Moriond17_gainSwitch_unc https://github.com/ECALELFS/ScalesSmearings.git
	cd $CMSSW_BASE/src

	# -- this ntuple maker -- #
	git clone https://github.com/KyeongPil-Lee/NtupleMaker.git Phys -b 80X

	# -- compile -- #
	scram b -j 20 >&log&

	# -- Example for CRAB submission -- #
	source /cvmfs/cms.cern.ch/crab3/crab.sh
	voms-proxy-init --voms cms

	cd Phys/DYntupleMaker/ntuples/CRABSubmit
	python crab3cfg_MC_DYMassBinned_Moriond17.py

	# -- Check CRAB status -- #
	python CRAB_Status.py