#####################################
# job script example with GE options
#####################################
#!/bin/bash

#$ -P P_juno
#$ -N muon
#$ -l sps=1
#$ -l ct=18:00:00

#$ -e /sps/hep/juno/pineresr/SUBMIT/muon/MuonTT_eutput_${SGE_TASK_ID}
#$ -o /sps/hep/juno/pineresr/SUBMIT/muon/MuonTT_output_${SGE_TASK_ID}




cd /sps/hep/juno/pineresr/SUBMIT/

echo $JUNOTOP
echo $JUNO_OFFLINE_OFF
unset JUNO_OFFLINE_OFF
source /pbs/throng/juno/SNiPER_OS7.versions/Pre-Release/J18v2r1-Pre1/setup.sh

Muon.exe -n 100000 -mult 0 -seed ${SGE_TASK_ID} -v TT -o muon_${SGE_TASK_ID}.txt -s juno -r Yes -music_dir ${JUNOTOP}/data/Generator/Muon/data


python /pbs/throng/juno/SNiPER_OS7.versions/Pre-Release/J18v2r1-Pre1/offline/Examples/Tutorial/share/tut_detsim.py --no-gdml --start-evtid 0 --no-cd --seed ${SGE_TASK_ID} --no-wp --evtmax 100000 --output MuonTT_${SGE_TASK_ID}.root --user-output MuonTT_${SGE_TASK_ID}_user.root hepevt --exe Muon --file muon_${SGE_TASK_ID}.txt


