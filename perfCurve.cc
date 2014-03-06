#include "TH1F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"



void perfCurve(){

  TFile * f = new TFile("FourTop_EventSelection_wMETCut_Mu.root");
  TH1F * hMatched = (TH1F*)f->Get("MultiSamplePlot_Chi2_Matched/Chi2Matched_TTJets");
  TH1F * hUnMatched = (TH1F*)f->Get("MultiSamplePlot_Chi2_UnMatched/Chi2UnMatched_TTJets");

  TGraph * gr = new TGraph();

  int nbins =  hUnMatched->GetNbinsX();

  double scale_matched = 1./hMatched->Integral(); 
  double scale_unmatched = 1./hUnMatched->Integral(); 

  hMatched->Scale(scale_matched);
  hUnMatched->Scale(scale_unmatched);

  for (int a = 0; a < nbins ; a++){

    double eff_matched = (1-hMatched->Integral(0,a));
    double eff_unmatched =(1- hUnMatched->Integral(0,a));  
    gr->SetPoint(a, eff_matched, eff_unmatched);

    }

  gr->Draw("ALP");
}
