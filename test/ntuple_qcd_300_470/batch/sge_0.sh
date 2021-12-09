
#!/bin/sh
#
#(make sure the right shell will be used)
#$ -S /bin/sh
#$ -l site=hh
#$ -l distro=sld6
#
#(the cpu time for this job)
#$ -l h_rt=05:55:00
#
#(the maximum memory usage of this job)
#$ -l h_vmem=4096M
#$ -l cvmfs
#(stderr and stdout are merged together to stdout)
#$ -j y
#$ -m a
#$ -cwd -V
#( -l h_stack=1536M) #try with small stack
#$ -pe local 1 -R y
#$ -P af-cms

export LOGDIR=/afs/cern.ch/work/a/ademoor/CMSSW_11_1_2_patch3/src/DeepNTuples/test/ntuple_qcd_300_470/batch/
export JOB=0
export NTUPLEOUTFILEPATH=/eos/home-a/ademoor/UL17_DV_mini_dataset/111X/Thu_153359_test/ntuple_qcd_300_470/output/ntuple_qcd_300_470_0.root

/afs/cern.ch/work/a/ademoor/CMSSW_11_1_2_patch3/src/DeepNTuples/test/ntuple_qcd_300_470/batchscript.sh /afs/cern.ch/work/a/ademoor/CMSSW_11_1_2_patch3/src/DeepNTuples/test/DeepNtuplizer.py inputScript=QCDPt300to470TuneCP513TeVpythia8RunIISummer19UL17MiniAOD106Xmc2017realisticv6v2MINIAODSIM nJobs=2 job=0 gluonReduction=0.5
            