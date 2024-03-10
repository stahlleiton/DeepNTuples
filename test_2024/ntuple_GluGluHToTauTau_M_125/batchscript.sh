#!/bin/bash
        
workdir=""
if [ $SGE_CELL ]
then
workdir=`mktemp -d -t DeepNTuplesXXXXXX`        
cd $workdir
workdir=$workdir"/"
fi

exec > "$PWD/stdout.txt" 2>&1
echo "JOBSUB::RUN job running"
trap "echo JOBSUB::FAIL job killed" SIGTERM
export OUTPUT=$workdir"ntuple_GluGluHToTauTau_M_125"
cd /afs/cern.ch/work/a/ademoor/New_Tagger_Dev/MC_ParT_2024/CMSSW_14_1_0_pre0
eval `scramv1 runtime -sh`
export PYTHONPATH=/afs/cern.ch/user/a/ademoor/.deepntuples_scripts_tmp/:$PYTHONPATH
which cmsRun
cd -
cmsRun -n 1 "$@" outputFile=$OUTPUT 
exitstatus=$?
if [ $exitstatus != 0 ]
then
   echo JOBSUB::FAIL job failed with status $exitstatus
else
   pwd
   ls -ltr $OUTPUT*.root
   if [ $SGE_CELL ]
   then
     cp $OUTPUT*.root $NTUPLEOUTFILEPATH
   else
     eoscp $OUTPUT*.root $NTUPLEOUTFILEPATH
   fi
   exitstatus=$?
   rm -f $OUTPUT*.root
   if [ $exitstatus != 0 ]
   then
     echo JOBSUB::FAIL job failed with status $exitstatus
   else
     echo JOBSUB::SUCC job ended sucessfully
   fi
fi
rm -f $OUTPUT*.root
if [ $workdir ]
then
# JOB is only defined for SGE submit
cp $workdir/stdout.txt $LOGDIR/con_out.$JOB.out
rm -rf $workdir
fi
exit $exitstatus
        