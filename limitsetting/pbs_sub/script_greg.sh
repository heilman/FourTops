#!/bin/bash

PROCESS="XXX"

cd /localgrid/ghammad/CMSSW_5_3_11/src/TopBrussels/TopTreeAnalysis/TopFCNC/macros

echo ">> script.sh is checking where it is"
pwd

echo ">> script.sh is checking what the HOME variable is"
echo $HOME

echo ">> script.sh is setting ROOT env variables"
export ROOTSYS=/localgrid/$USER/cmssoft/root_5.30.06/root 
export PATH=$PATH:$ROOTSYS/bin 
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ROOTSYS/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/lib

echo ">> script.sh is compiling TopFCNC_KinTreeMaker.cc "
./CompileAndRun.sh TopFCNC_KinTreeMaker_${PROCESS}.cc

echo ">> script.sh is listing files and directories in the current location"
ls -l

echo ">> script.sh is copying the executable to /scratch"
cp TopFCNC_KinTreeMaker_${PROCESS} /scratch

echo ">> script.sh is copying the needed Tree to /scratch"
cp TopFCNC_EventSelection_DiMuTrigger_Run2012ABCD_Skim2Mu3Jets_Trees/TopFCNC_EventSelection_DiMuTrigger_Run2012ABCD_Skim2Mu3Jets_TTree_${PROCESS}.root /scratch

echo ">> script.sh is copying the needed Btagging infos to /scratch"
cp -r TopFCNC_BtaggingEff_DiMuTrigger_Run2012ABCD_Doc /scratch

echo ">> script.sh is copying the needed resolution files to /scratch"
cp -r ResolutionFiles /scratch

echo ">> script.sh is cd'ing into /scratch partition"
cd /scratch

echo ">> ls of /scratch partition"
ls -l

echo ">> script.sh is executing the TopFCNC_KinTreeMaker executable"
./TopFCNC_KinTreeMaker_${PROCESS}

echo ">> script.sh is copying the output file back to the localgrid sandbox"
cp TopFCNC_Analysis_DiMuTrigger_Run2012ABCD_Skim2Mu3Jets_TTree_${PROCESS}.root /localgrid/ghammad/CMSSW_5_3_11/src/TopBrussels/TopTreeAnalysis/TopFCNC/macros/

echo ">> script.sh is cleaning the /scratch partition"
rm -f TopFCNC_Analysis_DiMuTrigger_Run2012ABCD_Skim2Mu3Jets_TTree_${PROCESS}.root
rm -rf ResolutionFiles
rm -f TopFCNC_KinTreeMaker_${PROCESS}
rm -rf TopFCNC_BtaggingEff_DiMuTrigger_Run2012ABCD_Doc

echo ">> script.sh is cleaning the /localgrid/ghammad/CMSSW_5_3_11/src/TopBrussels/TopTreeAnalysis/TopFCNC/macros dir"
cd /localgrid/ghammad/CMSSW_5_3_11/src/TopBrussels/TopTreeAnalysis/TopFCNC/macros
rm -f TopFCNC_KinTreeMaker_${PROCESS}.cc
rm -f TopFCNC_KinTreeMaker_${PROCESS}
