
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

export LOGDIR=/afs/cern.ch/work/a/ademoor/New_Tagger_Dev/MC_ParT_2024/CMSSW_14_1_0_pre0/src/DeepNTuples/test_2024/ntuple_VBFHToCC_M_125/batch/
export JOB=4
export NTUPLEOUTFILEPATH=/eos/cms/store/group/phys_btag/ParT_2024/Sun_180550_test_2024/ntuple_VBFHToCC_M_125/output/ntuple_VBFHToCC_M_125_4.root

/afs/cern.ch/work/a/ademoor/New_Tagger_Dev/MC_ParT_2024/CMSSW_14_1_0_pre0/src/DeepNTuples/test_2024/ntuple_VBFHToCC_M_125/batchscript.sh /afs/cern.ch/work/a/ademoor/New_Tagger_Dev/MC_ParT_2024/CMSSW_14_1_0_pre0/src/DeepNTuples/test_2024/DeepNtuplizer.py inputScript=VBFHToCCM125TuneCP513p6TeVpowhegpythia8Run3Summer23BPixMiniAODv4130XmcRun32023realisticpostBPixv5v2MINIAODSIM nJobs=5 job=4 
            