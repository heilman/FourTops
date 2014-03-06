#Script to run tttt analysis from start to finish



#(A) - Run ttbar systematic samples, these are needed to make the DATA/MC plots.


#1. Scaling uncertainties

    #Scale up
    #edit xml to point to scale up ttbar only:
 #   sed -i 's/add="1"/add="0"/g' config/test_fullsamples.xml


    sed -i 's/name="TTJets_ScaleUp"  title="t\bar{t}+jets"   add="0"/name="TTJets_ScaleUp"  title="t\bar{t}+jets"   add="1"/g' config/test_fullsamples.xml


    #edit macro to turn on scleup bool
    #compile and run macro
 #   makeMacro
 #   ./MACRO


    #Scale down
    #edit xml to point to scale up ttbar only: