# NtupleMaker
* Ntuple maker for physics analysis
   * 80X ntuple (for Z'/DY analysis with 2016 dataset)

## Environment
	export SCRAM_ARCH=slc6_amd64_gcc530
	export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
	source $VO_CMS_SW_DIR/cmsset_default.sh
	cmsrel CMSSW_8_0_32
	cd CMSSW_8_0_32/src
	cmsenv
	git cms-init

	# -- EGM corrections -- # (https://twiki.cern.ch/twiki/bin/view/CMS/EGMSmearer#!ECAL%20scale%20and%20resolution%20corre)
	git cms-merge-topic cms-egamma:EGM_gain_v1
	cd EgammaAnalysis/ElectronTools/data
	# download the txt files with the corrections
	git clone https://github.com/ECALELFS/ScalesSmearings.git
	cd ScalesSmearings
	git checkout Legacy2016_v1 # a tricky part: start
	cp Legacy2016_07Aug2017_pho_unc_*.dat ../ # I can’t understand yet why this dat files are needed
	git checkout Moriond17_23Jan_v2
	mv ../Legacy2016_07Aug2017_pho_unc_*.dat ./ # a tricky part: end
	cd $CMSSW_BASE/src

	# — 80X Photon ID — # (https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Recipe%20for%20regular%20users%20for%208.0)
	git cms-merge-topic ikrav:egm_id_80X_v3_photons

	# -- MET corrections (phi) -- #
	git cms-merge-topic cms-met:METRecipe_8020 -u
	git cms-merge-topic cms-met:METRecipe_80X_part2 -u
	git cms-addpkg JetMETCorrections

	# — L1 pre-firing inefficiency — # (https://twiki.cern.ch/twiki/bin/viewauth/CMS/L1ECALPrefiringWeightRecipe)
	git cms-merge-topic lathomas:L1Prefiring_8_0_32

	# -- this ntuple maker -- #
	git clone https://github.com/KyeongPil-Lee/NtupleMaker.git Phys
	cd Phys
	git checkout v2.6-Muon_ID_Flags
	cd ..

	# -- copy the MET correction recipe (80X) -- #
	cp Phys/DYntupleMaker/python/METCorr/multPhiCorr_*.py JetMETCorrections/Type1MET/python/

	# -- copy the the prefiring map -- #
	cp L1Prefiring/EventWeightProducer/files/L1PrefiringMaps_new.root Phys/DYntupleMaker/ntuples/withEGMcorrection/
	cp L1Prefiring/EventWeightProducer/files/L1PrefiringMaps_new.root Phys/DYntupleMaker/ntuples/CRABSubmit/

	# -- compile -- #
	scram b -j 20 >&log&

	# -- Example for CRAB submission -- #
	source /cvmfs/cms.cern.ch/crab3/crab.sh
	voms-proxy-init --voms cms --valid 168:00

	cd Phys/DYntupleMaker/ntuples/CRABSubmit
	python crab3cfg_MC_DYMassBinned_Moriond17.py

	# -- Check CRAB status -- #
	python CRAB_Status.py
