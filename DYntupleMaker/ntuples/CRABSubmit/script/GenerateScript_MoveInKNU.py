#!/usr/bin/env python
import os
import sys
import time
import subprocess

user = 'dmpai'
#xrootdPathBase = '/u/user/dmpai/SE_UserHome/_v2p4_/'
#xrootdPathBase = '/u/user/dmpai/SE_UserHome/_v2p5_/'
# xrootdPathBase = '/u/user/kplee/SE_UserHome/_v2p6_'
xrootdPathBase = '/u/user/kplee/SE_UserHome/_v2p7_'

## -- tamsa2 -- ##
#outputPathBase = '/data9/DATA/DYntuple/v2.3/'
#hostname = '147.47.242.67'

## -- tamsa1 -- ##
# outputPathBase = '/data1/kplee/DYntuple/80X'
# hostname = '147.47.242.42'

## -- prime -- ##
#outputPathBase = '/scratch/DYntuple/v2.3'
#hostname = '147.47.50.219'

## -- knu -- ##
#outputPathBase = '/u/user/dmpai/SE_UserHome/_prime_/DYntuple/v2.4/'
#outputPathBase = '/u/user/dmpai/SE_UserHome/_prime_/DYntuple/v2.5/'
outputPathBase = '/u/user/kplee/SE_UserHome/DYntuple/v2.6'
outputPathBase = '/u/user/kplee/SE_UserHome/DYntuple/v2.7'
hostname = 'cms.knu.ac.kr'


print "=" * 100
print "output directory: " + outputPathBase
CRABDirs = [
'crab_DYLL_M1000to1500',
'crab_DYLL_M100to200_fixed',
'crab_DYLL_M10to50_ext1v1',
'crab_DYLL_M10to50_v2',
'crab_DYLL_M1500to2000',
'crab_DYLL_M2000to3000',
'crab_DYLL_M200to400_fixed',
'crab_DYLL_M400to500',
'crab_DYLL_M500to700',
'crab_DYLL_M50toInf',
'crab_DYLL_M700to800',
'crab_DYLL_M800to1000',
'crab_QCDEMEnriched_Pt170to300',
'crab_QCDEMEnriched_Pt30to50',
'crab_QCDEMEnriched_Pt30to50_ext1',
'crab_QCDEMEnriched_Pt50to80',
'crab_QCDEMEnriched_Pt50to80_ext1',
'crab_QCDEMEnriched_Pt80to120',
'crab_QCDMuEnriched_Pt1000toInf_ext1',
'crab_QCDMuEnriched_Pt120to170',
'crab_QCDMuEnriched_Pt120to170_backup',
'crab_QCDMuEnriched_Pt15to20',
'crab_QCDMuEnriched_Pt170to300',
'crab_QCDMuEnriched_Pt170to300_backup',
'crab_QCDMuEnriched_Pt170to300_ext1',
'crab_QCDMuEnriched_Pt20to30',
'crab_QCDMuEnriched_Pt300to470',
'crab_QCDMuEnriched_Pt300to470_ext1',
'crab_QCDMuEnriched_Pt300to470_ext2',
'crab_QCDMuEnriched_Pt470to600',
'crab_QCDMuEnriched_Pt470to600_ext1',
'crab_QCDMuEnriched_Pt470to600_ext2',
'crab_QCDMuEnriched_Pt50to80',
'crab_QCDMuEnriched_Pt600to800',
'crab_QCDMuEnriched_Pt600to800_backup',
'crab_QCDMuEnriched_Pt600to800_ext1',
'crab_QCDMuEnriched_Pt800to1000',
'crab_QCDMuEnriched_Pt800to1000_ext1',
'crab_QCDMuEnriched_Pt800to1000_ext2',
'crab_QCDMuEnriched_Pt80to120_ext1',
'crab_ST_tW',
'crab_ST_tbarW',
'crab_SinglePhoton_Run2016B',
'crab_SinglePhoton_Run2016C',
'crab_SinglePhoton_Run2016D',
'crab_SinglePhoton_Run2016E',
'crab_SinglePhoton_Run2016F',
'crab_SinglePhoton_Run2016G',
'crab_SinglePhoton_Run2016Hver2',
'crab_WJetsToLNu_amcatnlo_ext',
'crab_WJetsToLNu_amcatnlo_ext2v5',
'crab_WW',
'crab_WZ',
'crab_ZZ',
'crab_ttbar',
'crab_ttbarBackup',
'crab_ttbar_M1000toInf',
]


# FileList = os.listdir(".")
# for File in FileList:
# 	if "crab_" in File and os.path.isdir( File ):
# 		print "Recognized crab directory: " + File
# 		CRABDirs.append( File )

print CRABDirs
print "=" * 100

TIME = time.strftime('%Y%m%d_%H%M%S', time.localtime(time.time()))
f = open( "CRAB_Getoutput_v%s.sh" % (TIME), 'w' )

f.write( "#!bin/bash\n\n" )

for crabDir in CRABDirs:
	outputDir = crabDir.split("crab_")[1]
	print "\toutputDir: " + outputDir

	# -- make the directory in SNU cluster -- #
	cmd_mkdir = "mkdir %s/%s" %( outputPathBase, outputDir)

	result_mkdir = subprocess.Popen(cmd_mkdir, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	(stdout_mkdir, stderr_mkdir) = result_mkdir.communicate()
	print "#" * 100
	print cmd_mkdir+'\n'
	print "[stdout]"
	print stdout_mkdir
	print "[stderr]"
	print stderr_mkdir
	print "#" * 100 + "\n"

	if "File exists" in stderr_mkdir:
		print "directory already exists ... please check"
		sys.exit()

	# -- find corresponding directory in xrootd -- #
	# -- Ref. of os.walk: https://wikidocs.net/39 -- #
	for (path, dir1, files) in os.walk( xrootdPathBase ):
		
		if crabDir in dir1:
			print "find the directory corresponding to %s: \n\t%s/%s" % (crabDir, path, crabDir)

			Set_path = set()
			nFile = 0

			for (path2, dir2, files2) in os.walk( "%s/%s" % (path, crabDir) ):
				
				for filename in files2:
					ext = os.path.splitext(filename)[-1]
					if ext == '.root' and "failed" not in path2:
						nFile += 1
						Set_path.add( path2 )


			f.write( 'echo "[Starting to move %s (# retrieved jobs: %d)]"\n' % (crabDir, nFile) )
			for path_root in Set_path:
				cmd_mv = "mv -v %s/*.root %s/%s" % ( path_root, outputPathBase, outputDir )
				f.write( cmd_mv + "\n" )
				
			f.write( 'echo "[move of %s is finished]"\n\n' % (crabDir) )

			break

print "[move completed crab jobs to _finished directory]"
for crabDir in CRABDirs:
	print "mv "+crabDir+" _finished"

print "\nShell script [CRAB_Getoutput_v%s.sh] is generated" % (TIME)
print "source CRAB_Getoutput_v%s.sh >&log_CRAB_Getoutput_v%s.txt&" % (TIME, TIME)
