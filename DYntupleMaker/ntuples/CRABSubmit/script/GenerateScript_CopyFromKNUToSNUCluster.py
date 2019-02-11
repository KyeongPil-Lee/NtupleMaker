#!/usr/bin/env python
import os
import sys
import time
import subprocess

user = 'dmpai'
#xrootdPathBase = '/xrootd/store/user/dpai/_v2p3_/'
xrootdPathBase = '/xrootd/store/user/dpai/_v2p3_supplement_from_CMSSW_8_0_31_/'

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
outputPathBase = '/u/user/dmpai/SE_UserHome/_prime_/DYntuple/v2.3/'
hostname = 'cms.knu.ac.kr'


print "=" * 100
print "output directory: " + outputPathBase
CRABDirs = [
'crab_DYLL_M50toInf',
'crab_ST_tW',
'crab_ST_tbarW',
'crab_SingleMuon_RunB',
'crab_SingleMuon_RunC',
'crab_SingleMuon_RunD',
'crab_SingleMuon_RunE',
'crab_SingleMuon_RunF',
'crab_SingleMuon_RunG',
'crab_SingleMuon_RunHver2',
'crab_SingleMuon_RunHver3',
'crab_WW',
'crab_WZ',
'crab_ZZ',
'crab_ttbar',
'crab_ttbarBackup',
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
	if hostname is '147.47.50.219':
		cmd_mkdir = "ssh -p 50002 %s@%s 'mkdir %s/%s'" %( user, hostname, outputPathBase, outputDir)
	else:
		cmd_mkdir = "ssh %s@%s 'mkdir %s/%s'" %( user, hostname, outputPathBase, outputDir)

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
		print "directory in SNU-cluster already exists ... please check"
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
				if hostname is '147.47.50.219':
					cmd_scp = "scp -P 50002 -v %s/*.root %s@%s:%s/%s" % ( path_root, user, hostname, outputPathBase, outputDir )
				else:
					cmd_scp = "scp -v %s/*.root %s@%s:%s/%s" % ( path_root, user, hostname, outputPathBase, outputDir )
				f.write( cmd_scp + "\n" )
				
			f.write( 'echo "[copy of %s is finished]"\n\n' % (crabDir) )

			break

print "[move completed crab jobs to _finished directory]"
for crabDir in CRABDirs:
	print "mv "+crabDir+" _finished"

print "\nShell script [CRAB_Getoutput_v%s.sh] is generated" % (TIME)
print "source CRAB_Getoutput_v%s.sh >&log_CRAB_Getoutput_v%s.txt&" % (TIME, TIME)
