#!/bin/bash
# Usage: . test.sh b -->  for batch 
#        . test.sh i -->  for interactive
#        . test.sh cb angle --> for compiling and batch for experiment 1
#        . test.sh ci --> for compiling and interactive
#        . test.sh c -->  for compiling only

if [ $1 == "c" ] || [ $1 == "cb" ] ||  [ $1 == "ci" ] 
then
echo -e  Setting up environment
. muonSetupMac.sh
echo -e  Set up done! 
echo -e  Cleaning last make session
make clean
echo -e  Starting new make session
make 
if [ $1 == "cb" ]
then
echo -e  Batch mode
. run_batch.sh
#cp muonTree01.root muonTree01_angle_$2_w_89.root
#mv muonTree01_angle_$2_w_89.root /Users/ttumuon/hep/g4/g4user/MuonSC8/rp/v5/sim/experimen1_files/
#echo Moved the file!
elif [ $1 == "ci" ]
then
echo -e  Interactive mode:
. interactive_run.sh
fi
fi

if [ $1 == "b" ]
then
echo -e  Batch mode
. runG4CRY_batch.sh
elif [ $1 == "i" ]
then
echo -e  Interactive mode:
. interactive_run.sh
fi
