# NtupleMaker
Ntuple maker for differental DY cross section measurements with full Run-2 data

Twiki: https://twiki.cern.ch/twiki/bin/view/CMS/SNUCMSYooDYntuple



## Recipe (80X)

* This includes a few steps to implement several corrections by hand
  * It would be simplified more in the next version of reprocessed samples (94X or 106X)

### Initial setup & submit a ntupler job to CRAB

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
	git clone https://github.com/KyeongPil-Lee/NtupleMaker.git Phys -b 80X # latest version
	
	# -- copy the MET correction recipe (80X) -- #
	cp Phys/DYntupleMaker/python/METCorr/multPhiCorr_*.py JetMETCorrections/Type1MET/python/
	
	# -- copy the the prefiring map -- #
	cp L1Prefiring/EventWeightProducer/files/L1PrefiringMaps_new.root Phys/DYntupleMaker/ntuples/withEGMcorrection/
	cp L1Prefiring/EventWeightProducer/files/L1PrefiringMaps_new.root Phys/DYntupleMaker/ntuples/CRABSubmit/
	
	# -- compile -- #
	scram b -j 20 >&log& # -- it would take some time
	
	# -- Example for CRAB submission -- #
	source /cvmfs/cms.cern.ch/crab3/crab.sh
	voms-proxy-init --voms cms --valid 168:00
	
	cd Phys/DYntupleMaker/ntuples/CRABSubmit
	python crab3cfg_MC_DYMassBinned_Moriond17.py --version='2p8'



### Python configuration of the ntuple maker

```
/Phys/DYntupleMaker/ntuples/withEGMcorrection/ntupler_arg.py
```

* It requires several arguments
  * ```globalTag```: global tag to use
  * ```inputFile```: input EDM file (MINIAOD). ```file:``` prefix is needed in front of the path of the file if it is local file
  * ```useSinglePhotonTrigger```: inlcude single photon trigger on top of usual lepton triggers
  * ```nEvent```: number of event to run
  * ```isMC```: true if it is MC
  * ```isSignalMC```: true if it is signal MC (DY). Then it will save LHE information, which is used to estimate PDF systematics



* Usage:

  ```
  cmsRun ntupler_arg.py globalTag=<global tag> inputFile=<inputFile> useSinglePhotonTrigger=<1 or 0> nEvent=<number> isMC=<1 or 0> isSignalMC=<1 or 0>
  ```

  * Example (signal MC)

    ```
    cmsRun ntupler_arg.py \
    globalTag=80X_mcRun2_asymptotic_2016_TrancheIV_v6 \
    inputFile=file:/u/user/kplee/scratch/ROOTFiles_Test/80X/MINIAOD_DYLL_M50toInf_Morind17.root \
    useSinglePhotonTrigger=1 \
    nEvent=3000 \
    isMC=1 \
    isSignalMC=1 >&ntupler_arg_signalMC.log&
    ```



### CRAB scripts

Exists under

```
Phys/DYntupleMaker/ntuples/CRABSubmit
```



Common usage

```
python <configuration name> --version=<version>
```

* Example

  ```
  python crab3cfg_MC_DYMassBinned_Moriond17.py --verison='2p8'
  ```



**NOTE**: You may need to update ```config.Site.storageSite``` in the configuration to your site before submission
(Currently it is set with T2_KR_KNU)



* **MC** configurations
  * ```crab3cfg_MC_DYMassBinned_Moriond17.py```: singal DY samples
  * ```crab3cfg_MC_BackgroundMC_Moriond17.py```: usual background MC (ttbar, single top and diboson)
  * ```crab3cfg_MC_WJetsToLNu_Moriond17.py```: W+jets samples (aMC@NLO)
  * ```crab3cfg_MC_GJetsPtBinned_Sherpa_Moriond17.py```: gamma+jet samples
  * ```crab3cfg_MC_QCDPtBinned_Moriond17.py```: QCD mu-enriched and e-enriched samples
* **Data** configurations
  * ```crab3cfg_Run2016_All.py```: SingleMuon and DoubleEG samples
  * ```crab3cfg_Run2016_Egamma.py```: DoubleEG samples only
  * ```crab3cfg_Run2016_SinglePhoton.py```: SinglePhoton samples

