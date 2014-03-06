#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"



void main(){


  string filename = "FourTopnoRho_EventSelection_El_preApp.root";

  TFile * f = new TFile(filename);



  TH1F * h_Data = new (*TFile)f->Get("MVA_Data");
  TH1F * h_TTJets_ll = new (*TFile)f->Get("MVA_TTJets_ll");
  TH1F * h_TTJets_bb = new (*TFile)f->Get("MVA_TTJets_bb");
  TH1F * h_TTJets_cc = new (*TFile)f->Get("MVA_TTJets_cc");


  TCanvas * c = new TCanvas();

  h_Data->Draw();
  h_TTJets_ll->Draw("same");

  h_TTJets_ll->SetFillColor(kRed);


}
