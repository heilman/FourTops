#include "TFile.h"
#include "TROOT.h"
#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooRealVar.h"

#include "RooStats/ModelConfig.h"
#include "RooStats/ProfileLikelihoodCalculator.h"
#include "RooStats/LikelihoodInterval.h"
#include "RooStats/LikelihoodIntervalPlot.h"

using namespace RooFit;
using namespace RooStats;



void test(){

  TFile *file = TFile::Open("testmodel_preapp.root");

// get the workspace out of the file
  RooWorkspace* w = (RooWorkspace*) file->Get("combined");
  if(!w){
    cout <<"workspace not found" << endl;
    return;
  }

  // get the modelConfig out of the file
  ModelConfig* mc = (ModelConfig*) w->obj("ModelConfig");

  // get the modelConfig out of the file
  RooAbsData* data = w->data("obsData");

  // make sure ingredients are found
  if(!data || !mc){
    w->Print();
    cout << "data or ModelConfig was not found" <<endl;
    return;
  }


  mc->GetPdf();
  mc->GetPdf()->fitTo(*data,Minimizer("Minuit"), Strategy(2), PrintLevel(3), Save(true) ); 



}
