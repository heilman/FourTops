#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"

void process(){
  TFile* f_noRW = new TFile("FourTop_EventSelection_wMETCut_Mu_noPTRW.root"); 
  TFile* f = new TFile("FourTop_EventSelection_wMETCut_Mu.root"); 

  TH1F * h_noRW = (TH1F*)f_noRW->Get("MultiSamplePlot_MVA/MVA_TTJets");
  TH1F * h = (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_TTJets");


  double scale = (h_noRW->Integral())/(h->Integral());

  h->Scale(scale);


  h_noRW->SetLineColor(kBlue);
  h->SetLineColor(kRed);
  h_noRW->SetMarkerColor(kBlue);
  h->SetMarkerColor(kRed);

  h_noRW->SetTitle("");

  h_noRW->Draw();
  h->Draw("same");

  TLegend* leg = new TLegend(0.6, 0.73, .97, .95);
  leg->SetFillColor(0);
  leg->AddEntry(h_noRW, "without reweight", "lp");
  leg->AddEntry( h   , "with reweight", "lp");
  leg->Draw();


}
