#!/bin/bash

for f in `ls TopFCNC_EventSelection_DiMuTrigger_Run2012ABCD_Skim2Mu3Jets_Trees/ | cut -d'_' -f 7- | cut -d'.' -f1`
do
#	echo $f
	cat script.sh | sed "s/XXX/${f}/" > myscript_${f}.sh
	cat TopFCNC_KinTreeMaker.cc | sed "s/XXX/${f}/" > TopFCNC_KinTreeMaker_${f}.cc
#	qsub -q localgrid@cream02 -o script_${f}.stdout -e script_${f}.stderr -l walltime=15:00:00 script_${f}.sh
	qsub -q localgrid@cream02 -o myscript_${f}.stdout -e myscript_${f}.stderr myscript_${f}.sh
done
