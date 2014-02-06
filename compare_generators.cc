#include "TCanvas.h"
#include "TH1F.h"
#include "TFile.h"
#include "TLegend.h"

void compare(){

TFile * f = new TFile("FourTop_EventSelection_wMETCut_Mu_compareGen.root");

//TH1F * hMVA_tttt_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_MVA/MVA_NP_overlay_TTTT");


TH1F * f_data = (TH1F*)f->Get("MultiSamplePlot_NbOfSelectedBJets/NbOfSelectedBJets_Data");
TH1F * f_madgraph_pythia = (TH1F*)f->Get("MultiSamplePlot_NbOfSelectedBJets/NbOfSelectedBJets_TTJets");
TH1F * f_powheg_pythia   = (TH1F*)f->Get("MultiSamplePlot_NbOfSelectedBJets/NbOfSelectedBJets_TTJets_powheg_pythia");
TH1F * f_powheg_herwig   = (TH1F*)f->Get("MultiSamplePlot_NbOfSelectedBJets/NbOfSelectedBJets_TTJets_powheg_herwig");
TH1F * f_mcatnlo_herwig = (TH1F*)f->Get("MultiSamplePlot_NbOfSelectedBJets/NbOfSelectedBJets_TTJets_mcatnlo");



 TCanvas * c1 = new TCanvas();

 f_data->Draw("HIST");
 f_madgraph_pythia->Draw("same");
 f_powheg_pythia->Draw("same");
 f_powheg_herwig->Draw("same");
 f_mcatnlo_herwig->Draw("same");

 f_madgraph_pythia->SetLineColor(kRed);
 f_powheg_pythia->SetLineColor(kGreen);
 f_powheg_herwig->SetLineColor(kBlue);
 f_mcatnlo_herwig->SetLineColor(kOrange);



 f_madgraph_pythia->SetMarkerColor(kRed);
 f_powheg_pythia->SetMarkerColor(kGreen);
 f_powheg_herwig->SetMarkerColor(kBlue);
 f_mcatnlo_herwig->SetMarkerColor(kOrange);



 TLegend *templeg = new TLegend(0.6, 0.6, 0.8, 0.8);
 templeg->SetFillColor(0);
 templeg->AddEntry(f_data,"Data","l");
 templeg->AddEntry(f_madgraph_pythia,"MadGraph + Pythia","p");
 templeg->AddEntry(f_powheg_pythia,"Powheg + pythia","p");
 templeg->AddEntry(f_powheg_herwig,"Powheg + Herwig","p");
 templeg->AddEntry(f_mcatnlo_herwig,"MC@NLO + Herwig","p");
 templeg->Draw();




}