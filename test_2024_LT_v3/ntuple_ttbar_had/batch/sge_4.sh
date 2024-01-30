
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

export LOGDIR=/afs/cern.ch/work/a/ademoor/New_Tagger_Dev/Run3_20224/CMSSW_13_0_13/src/DeepNTuples/test_2024_LT_v3/ntuple_ttbar_had/batch/
export JOB=4
export NTUPLEOUTFILEPATH=/eos/cms/store/group/phys_btag/ParticleEdge/s_tagging/Tue_114154_test_2024_LT_v3/ntuple_ttbar_had/output/ntuple_ttbar_had_4.root

/afs/cern.ch/work/a/ademoor/New_Tagger_Dev/Run3_20224/CMSSW_13_0_13/src/DeepNTuples/test_2024_LT_v3/ntuple_ttbar_had/batchscript.sh /afs/cern.ch/work/a/ademoor/New_Tagger_Dev/Run3_20224/CMSSW_13_0_13/src/DeepNTuples/test_2024_LT_v3/DeepNtuplizer.py inputScript=TTto4QTuneCP513p6TeVpowhegpythia8Run3Summer23BPixMiniAODv4130XmcRun32023realisticpostBPixv2v3MINIAODSIM nJobs=5 job=4 
            