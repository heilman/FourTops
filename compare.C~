#include "TFile.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TProfile.h"
#include <cstdlib>
#include <iostream>
#include "TStyle.h"

using namespace std;

void compare(){

  bool debug = true;

  gStyle->SetOptStat(00000);

  if(debug) std::cout <<"here 01"<< std::endl;

  map<string,TProfile*> histoProfile;

  TFile* f_compute_mu = new TFile("FourTop_EventSelection_wMETCut_Mu_forCorrelations.root"); 
  TFile* f_compute_el = new TFile("FourTop_EventSelection_El_forCorrelations.root"); 

  TFile* f_train_mu = new TFile("MasterMVA_Mu.root");
  TFile* f_train_el = new TFile("MasterMVA.root");

  //Comparison 1: mu + jets Vs. el + jets
  TH1F * hMVA_tttt_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_MVA/MVA_NP_overlay_TTTT");
  TH1F * hMVA_tttt_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_MVA/MVA_NP_overlay_TTTT");

  TH1F * hMVA_ttjj_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_MVA/MVA_TTJets");

  //  TH1F * hMVA_ttbb_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_MVA/MVA_TTJets_bb");
  // TH1F * hMVA_ttcc_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_MVA/MVA_TTJets_cc");
  // TH1F * hMVA_ttll_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_MVA/MVA_TTJets_ll");

  TH1F * hMVA_ttbb_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_MVA/MVA_TTJets_bb");
  TH1F * hMVA_ttcc_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_MVA/MVA_TTJets_cc");
  TH1F * hMVA_ttll_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_MVA/MVA_TTJets_ll");


 if(debug) std::cout <<"here 01"<< std::endl;

 // TH1F * hMVA_ttjj_compute_mu = (TH1F*)hMVA_ttbb_compute_mu->Clone();
 // hMVA_ttjj_compute_mu->Reset();
 // hMVA_ttjj_compute_mu->Add(hMVA_ttbb_compute_mu);
 // hMVA_ttjj_compute_mu->Add(hMVA_ttcc_compute_mu);
 // hMVA_ttjj_compute_mu->Add(hMVA_ttll_compute_mu);


  TH1F * hMVA_ttjj_compute_el = (TH1F*)hMVA_ttbb_compute_el->Clone();
  hMVA_ttjj_compute_el->Reset();
  hMVA_ttjj_compute_el->Add(hMVA_ttbb_compute_el);
  hMVA_ttjj_compute_el->Add(hMVA_ttcc_compute_el);
  hMVA_ttjj_compute_el->Add(hMVA_ttll_compute_el);

  hMVA_ttjj_compute_mu->SetMarkerColor(kRed);
  hMVA_ttjj_compute_el->SetMarkerColor(kGreen);

  hMVA_tttt_compute_mu->SetMarkerColor(kBlue);
  hMVA_tttt_compute_el->SetMarkerColor(kBlack);


  hMVA_ttjj_compute_mu->SetTitle("");

  hMVA_ttjj_compute_mu->GetXaxis()->SetTitleOffset(0.85);
  hMVA_ttjj_compute_mu->GetYaxis()->SetTitleOffset(0.85);

  TCanvas * c1 = new TCanvas();
  hMVA_ttjj_compute_mu->DrawNormalized();
  hMVA_ttjj_compute_el->DrawNormalized("same");
  hMVA_tttt_compute_mu->DrawNormalized("same");
  hMVA_tttt_compute_el->DrawNormalized("same");
 
  TLegend *leg = new TLegend(0.6, 0.3, .8, .8);
  leg->SetFillColor(0);
  leg->AddEntry(hMVA_ttjj_compute_mu,"ttjets mu","p");
  leg->AddEntry(hMVA_ttjj_compute_el,"ttjets el","p");
  leg->AddEntry(hMVA_tttt_compute_mu,"tttt mu","p");
  leg->AddEntry(hMVA_tttt_compute_el,"tttt el","p");
  leg->Draw();


  TCanvas * c2 = new TCanvas();


  TH1F * hHTb_tttt_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_HTb_SelectedJets/HTb_SelectedJets_NP_overlay_TTTT");
  TH1F * hHTb_tttt_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTb_SelectedJets/HTb_SelectedJets_NP_overlay_TTTT");

  TH1F * hHTb_ttjj_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTb_SelectedJets/HTb_SelectedJets_TTJets");

  //  TH1F * hHTb_ttbb_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTb_SelectedJets/HTb_SelectedJets_TTJets_bb");
  // TH1F * hHTb_ttcc_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTb_SelectedJets/HTb_SelectedJets_TTJets_cc");
  //TH1F * hHTb_ttll_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTb_SelectedJets/HTb_SelectedJets_TTJets_ll");

  TH1F * hHTb_ttbb_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_HTb_SelectedJets/HTb_SelectedJets_TTJets_bb");
  TH1F * hHTb_ttcc_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_HTb_SelectedJets/HTb_SelectedJets_TTJets_cc");
  TH1F * hHTb_ttll_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_HTb_SelectedJets/HTb_SelectedJets_TTJets_ll");

  //  TH1F * hHTb_ttjj_compute_mu = (TH1F*)hHTb_ttbb_compute_mu->Clone();
  // hHTb_ttjj_compute_mu->Reset();
  //hHTb_ttjj_compute_mu->Add(hHTb_ttbb_compute_mu);
  //hHTb_ttjj_compute_mu->Add(hHTb_ttcc_compute_mu);
  // hHTb_ttjj_compute_mu->Add(hHTb_ttll_compute_mu);


  TH1F * hHTb_ttjj_compute_el = (TH1F*)hHTb_ttbb_compute_el->Clone();
  hHTb_ttjj_compute_el->Reset();
  hHTb_ttjj_compute_el->Add(hHTb_ttbb_compute_el);
  hHTb_ttjj_compute_el->Add(hHTb_ttcc_compute_el);
  hHTb_ttjj_compute_el->Add(hHTb_ttll_compute_el);

  hHTb_ttjj_compute_mu->SetMarkerColor(kRed);
  hHTb_ttjj_compute_el->SetMarkerColor(kGreen);

  hHTb_tttt_compute_mu->SetMarkerColor(kBlue);
  hHTb_tttt_compute_el->SetMarkerColor(kBlack);

  hHTb_ttjj_compute_mu->GetXaxis()->SetTitleOffset(0.85);
  hHTb_ttjj_compute_mu->GetYaxis()->SetTitleOffset(0.85);

  hHTb_ttjj_compute_mu->SetTitle("");
  hHTb_ttjj_compute_mu->DrawNormalized();
  hHTb_ttjj_compute_el->DrawNormalized("same");
  hHTb_tttt_compute_mu->DrawNormalized("same");
  hHTb_tttt_compute_el->DrawNormalized("same");
 
  TLegend *leg2 = new TLegend(0.6, 0.3, .8, .8);
  leg2->SetFillColor(0);
  leg2->AddEntry(hHTb_ttjj_compute_mu,"ttjets mu","p");
  leg2->AddEntry(hHTb_ttjj_compute_el,"ttjets el","p");
  leg2->AddEntry(hHTb_tttt_compute_mu,"tttt mu","p");
  leg2->AddEntry(hHTb_tttt_compute_el,"tttt el","p");
  leg2->Draw();




  TCanvas * c3 = new TCanvas();


  TH1F * hMultiTopness_tttt_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_MultiTopness/MultiTopness_NP_overlay_TTTT");
  TH1F * hMultiTopness_tttt_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_MultiTopness/MultiTopness_NP_overlay_TTTT");

  TH1F * hMultiTopness_ttjj_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_MultiTopness/MultiTopness_TTJets");

  // TH1F * hMultiTopness_ttbb_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_MultiTopness/MultiTopness_TTJets_bb");
  // TH1F * hMultiTopness_ttcc_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_MultiTopness/MultiTopness_TTJets_cc");
  // TH1F * hMultiTopness_ttll_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_MultiTopness/MultiTopness_TTJets_ll");

  TH1F * hMultiTopness_ttbb_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_MultiTopness/MultiTopness_TTJets_bb");
  TH1F * hMultiTopness_ttcc_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_MultiTopness/MultiTopness_TTJets_cc");
  TH1F * hMultiTopness_ttll_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_MultiTopness/MultiTopness_TTJets_ll");

  //  TH1F * hMultiTopness_ttjj_compute_mu = (TH1F*)hMultiTopness_ttbb_compute_mu->Clone();
  //hMultiTopness_ttjj_compute_mu->Reset();
  // hMultiTopness_ttjj_compute_mu->Add(hMultiTopness_ttbb_compute_mu);
  //hMultiTopness_ttjj_compute_mu->Add(hMultiTopness_ttcc_compute_mu);
  //hMultiTopness_ttjj_compute_mu->Add(hMultiTopness_ttll_compute_mu);


  TH1F * hMultiTopness_ttjj_compute_el = (TH1F*)hMultiTopness_ttbb_compute_el->Clone();
  hMultiTopness_ttjj_compute_el->Reset();
  hMultiTopness_ttjj_compute_el->Add(hMultiTopness_ttbb_compute_el);
  hMultiTopness_ttjj_compute_el->Add(hMultiTopness_ttcc_compute_el);
  hMultiTopness_ttjj_compute_el->Add(hMultiTopness_ttll_compute_el);

  hMultiTopness_ttjj_compute_mu->SetMarkerColor(kRed);
  hMultiTopness_ttjj_compute_el->SetMarkerColor(kGreen);
  hMultiTopness_tttt_compute_mu->SetMarkerColor(kBlue);
  hMultiTopness_tttt_compute_el->SetMarkerColor(kBlack);

  hMultiTopness_ttjj_compute_mu->GetXaxis()->SetTitleOffset(0.85);
  hMultiTopness_ttjj_compute_mu->GetYaxis()->SetTitleOffset(0.85);

  hMultiTopness_ttjj_compute_mu->SetTitle("");
  hMultiTopness_ttjj_compute_mu->DrawNormalized();
  hMultiTopness_ttjj_compute_el->DrawNormalized("same");
  hMultiTopness_tttt_compute_mu->DrawNormalized("same");
  hMultiTopness_tttt_compute_el->DrawNormalized("same");
 
  TLegend *leg3 = new TLegend(0.6, 0.3, .8, .8);
  leg3->SetFillColor(0);
  leg3->AddEntry(hMultiTopness_ttjj_compute_mu,"ttjets mu","p");
  leg3->AddEntry(hMultiTopness_ttjj_compute_el,"ttjets el","p");
  leg3->AddEntry(hMultiTopness_tttt_compute_mu,"tttt mu","p");
  leg3->AddEntry(hMultiTopness_tttt_compute_el,"tttt el","p");
  leg3->Draw();




  TCanvas * c4 = new TCanvas();


  TH1F * hHTX_tttt_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_HTX/HTX_NP_overlay_TTTT");
  TH1F * hHTX_tttt_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTX/HTX_NP_overlay_TTTT");


   TH1F * hHTX_ttjj_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTX/HTX_TTJets");

 //  TH1F * hHTX_ttbb_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTX/HTX_TTJets_bb");
 // TH1F * hHTX_ttcc_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTX/HTX_TTJets_cc");
 // TH1F * hHTX_ttll_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTX/HTX_TTJets_ll");

  TH1F * hHTX_ttbb_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_HTX/HTX_TTJets_bb");
  TH1F * hHTX_ttcc_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_HTX/HTX_TTJets_cc");
  TH1F * hHTX_ttll_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_HTX/HTX_TTJets_ll");

  //  TH1F * hHTX_ttjj_compute_mu = (TH1F*)hHTX_ttbb_compute_mu->Clone();
  // hHTX_ttjj_compute_mu->Reset();
  //hHTX_ttjj_compute_mu->Add(hHTX_ttbb_compute_mu);
  //hHTX_ttjj_compute_mu->Add(hHTX_ttcc_compute_mu);
  //hHTX_ttjj_compute_mu->Add(hHTX_ttll_compute_mu);


  TH1F * hHTX_ttjj_compute_el = (TH1F*)hHTX_ttbb_compute_el->Clone();
  hHTX_ttjj_compute_el->Reset();
  hHTX_ttjj_compute_el->Add(hHTX_ttbb_compute_el);
  hHTX_ttjj_compute_el->Add(hHTX_ttcc_compute_el);
  hHTX_ttjj_compute_el->Add(hHTX_ttll_compute_el);

  hHTX_ttjj_compute_mu->SetMarkerColor(kRed);
  hHTX_ttjj_compute_el->SetMarkerColor(kGreen);
  hHTX_tttt_compute_mu->SetMarkerColor(kBlue);
  hHTX_tttt_compute_el->SetMarkerColor(kBlack);

  hHTX_ttjj_compute_mu->GetXaxis()->SetTitleOffset(0.85);
  hHTX_ttjj_compute_mu->GetYaxis()->SetTitleOffset(0.85);

  hHTX_ttjj_compute_mu->SetTitle("");
  hHTX_ttjj_compute_mu->DrawNormalized();
  hHTX_ttjj_compute_el->DrawNormalized("same");
  hHTX_tttt_compute_mu->DrawNormalized("same");
  hHTX_tttt_compute_el->DrawNormalized("same");
 
  TLegend *leg4 = new TLegend(0.6, 0.3, .8, .8);
  leg4->SetFillColor(0);
  leg4->AddEntry(hHTX_ttjj_compute_mu,"ttjets mu","p");
  leg4->AddEntry(hHTX_ttjj_compute_el,"ttjets el","p");
  leg4->AddEntry(hHTX_tttt_compute_mu,"tttt mu","p");
  leg4->AddEntry(hHTX_tttt_compute_el,"tttt el","p");
  leg4->Draw();



  TCanvas * c5 = new TCanvas();


  TH1F * hHTH_tttt_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_HTH/HTH_NP_overlay_TTTT");
  TH1F * hHTH_tttt_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTH/HTH_NP_overlay_TTTT");

  TH1F * hHTH_ttjj_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTH/HTH_TTJets");

  //  TH1F * hHTH_ttbb_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTH/HTH_TTJets_bb");
  //TH1F * hHTH_ttcc_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTH/HTH_TTJets_cc");
  // TH1F * hHTH_ttll_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTH/HTH_TTJets_ll");

  TH1F * hHTH_ttbb_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_HTH/HTH_TTJets_bb");
  TH1F * hHTH_ttcc_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_HTH/HTH_TTJets_cc");
  TH1F * hHTH_ttll_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_HTH/HTH_TTJets_ll");

  //TH1F * hHTH_ttjj_compute_mu = (TH1F*)hHTH_ttbb_compute_mu->Clone();
  // hHTH_ttjj_compute_mu->Reset();
  //hHTH_ttjj_compute_mu->Add(hHTH_ttbb_compute_mu);
  //hHTH_ttjj_compute_mu->Add(hHTH_ttcc_compute_mu);
  // hHTH_ttjj_compute_mu->Add(hHTH_ttll_compute_mu);


  TH1F * hHTH_ttjj_compute_el = (TH1F*)hHTH_ttbb_compute_el->Clone();
  hHTH_ttjj_compute_el->Reset();
  hHTH_ttjj_compute_el->Add(hHTH_ttbb_compute_el);
  hHTH_ttjj_compute_el->Add(hHTH_ttcc_compute_el);
  hHTH_ttjj_compute_el->Add(hHTH_ttll_compute_el);

  hHTH_ttjj_compute_mu->SetMarkerColor(kRed);
  hHTH_ttjj_compute_el->SetMarkerColor(kGreen);
  hHTH_tttt_compute_mu->SetMarkerColor(kBlue);
  hHTH_tttt_compute_el->SetMarkerColor(kBlack);

  hHTH_ttjj_compute_mu->GetXaxis()->SetTitleOffset(0.85);
  hHTH_ttjj_compute_mu->GetYaxis()->SetTitleOffset(0.85);

  hHTH_ttjj_compute_mu->SetTitle("");
  hHTH_ttjj_compute_mu->DrawNormalized();
  hHTH_ttjj_compute_el->DrawNormalized("same");
  hHTH_tttt_compute_mu->DrawNormalized("same");
  hHTH_tttt_compute_el->DrawNormalized("same");
 
  TLegend *leg5 = new TLegend(0.6, 0.3, .8, .8);
  leg5->SetFillColor(0);
  leg5->AddEntry(hHTH_ttjj_compute_mu,"ttjets mu","p");
  leg5->AddEntry(hHTH_ttjj_compute_el,"ttjets el","p");
  leg5->AddEntry(hHTH_tttt_compute_mu,"tttt mu","p");
  leg5->AddEntry(hHTH_tttt_compute_el,"tttt el","p");
  leg5->Draw();




  TCanvas * c6 = new TCanvas();

  TH1F * hHTRat_tttt_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_HTRat/HTRat_NP_overlay_TTTT");
  TH1F * hHTRat_tttt_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTRat/HTRat_NP_overlay_TTTT");

  TH1F * hHTRat_ttjj_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTRat/HTRat_TTJets");

  //  TH1F * hHTRat_ttbb_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTRat/HTRat_TTJets_bb");
  // TH1F * hHTRat_ttcc_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTRat/HTRat_TTJets_cc");
  //TH1F * hHTRat_ttll_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_HTRat/HTRat_TTJets_ll");

  TH1F * hHTRat_ttbb_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_HTRat/HTRat_TTJets_bb");
  TH1F * hHTRat_ttcc_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_HTRat/HTRat_TTJets_cc");
  TH1F * hHTRat_ttll_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_HTRat/HTRat_TTJets_ll");

  //  TH1F * hHTRat_ttjj_compute_mu = (TH1F*)hHTRat_ttbb_compute_mu->Clone();
  // hHTRat_ttjj_compute_mu->Reset();
  // hHTRat_ttjj_compute_mu->Add(hHTRat_ttbb_compute_mu);
  // hHTRat_ttjj_compute_mu->Add(hHTRat_ttcc_compute_mu);
  //hHTRat_ttjj_compute_mu->Add(hHTRat_ttll_compute_mu);


  TH1F * hHTRat_ttjj_compute_el = (TH1F*)hHTRat_ttbb_compute_el->Clone();
  hHTRat_ttjj_compute_el->Reset();
  hHTRat_ttjj_compute_el->Add(hHTRat_ttbb_compute_el);
  hHTRat_ttjj_compute_el->Add(hHTRat_ttcc_compute_el);
  hHTRat_ttjj_compute_el->Add(hHTRat_ttll_compute_el);

  hHTRat_ttjj_compute_mu->SetMarkerColor(kRed);
  hHTRat_ttjj_compute_el->SetMarkerColor(kGreen);
  hHTRat_tttt_compute_mu->SetMarkerColor(kBlue);
  hHTRat_tttt_compute_el->SetMarkerColor(kBlack);

  hHTRat_ttjj_compute_mu->GetXaxis()->SetTitleOffset(0.85);
  hHTRat_ttjj_compute_mu->GetYaxis()->SetTitleOffset(0.85);

  hHTRat_ttjj_compute_mu->SetTitle("");
  hHTRat_ttjj_compute_mu->DrawNormalized();
  hHTRat_ttjj_compute_el->DrawNormalized("same");
  hHTRat_tttt_compute_mu->DrawNormalized("same");
  hHTRat_tttt_compute_el->DrawNormalized("same");
 
  TLegend *leg6 = new TLegend(0.6, 0.3, .8, .8);
  leg6->SetFillColor(0);
  leg6->AddEntry(hHTRat_ttjj_compute_mu,"ttjets mu","p");
  leg6->AddEntry(hHTRat_ttjj_compute_el,"ttjets el","p");
  leg6->AddEntry(hHTRat_tttt_compute_mu,"tttt mu","p");
  leg6->AddEntry(hHTRat_tttt_compute_el,"tttt el","p");
  leg6->Draw();




  TCanvas * c7 = new TCanvas();

  TH1F * hNbOfSelectedJets_tttt_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_NbOfSelectedJets/NbOfSelectedJets_NP_overlay_TTTT");
  TH1F * hNbOfSelectedJets_tttt_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_NbOfSelectedJets/NbOfSelectedJets_NP_overlay_TTTT");

   TH1F * hNbOfSelectedJets_ttjj_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_NbOfSelectedJets/NbOfSelectedJets_TTJets");

  //  TH1F * hNbOfSelectedJets_ttbb_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_NbOfSelectedJets/NbOfSelectedJets_TTJets_bb");
  // TH1F * hNbOfSelectedJets_ttcc_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_NbOfSelectedJets/NbOfSelectedJets_TTJets_cc");
  // TH1F * hNbOfSelectedJets_ttll_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_NbOfSelectedJets/NbOfSelectedJets_TTJets_ll");

  TH1F * hNbOfSelectedJets_ttbb_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_NbOfSelectedJets/NbOfSelectedJets_TTJets_bb");
  TH1F * hNbOfSelectedJets_ttcc_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_NbOfSelectedJets/NbOfSelectedJets_TTJets_cc");
  TH1F * hNbOfSelectedJets_ttll_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_NbOfSelectedJets/NbOfSelectedJets_TTJets_ll");

  // TH1F * hNbOfSelectedJets_ttjj_compute_mu = (TH1F*)hNbOfSelectedJets_ttbb_compute_mu->Clone();
  // hNbOfSelectedJets_ttjj_compute_mu->Reset();
  // hNbOfSelectedJets_ttjj_compute_mu->Add(hNbOfSelectedJets_ttbb_compute_mu);
  // hNbOfSelectedJets_ttjj_compute_mu->Add(hNbOfSelectedJets_ttcc_compute_mu);
  // hNbOfSelectedJets_ttjj_compute_mu->Add(hNbOfSelectedJets_ttll_compute_mu);


  TH1F * hNbOfSelectedJets_ttjj_compute_el = (TH1F*)hNbOfSelectedJets_ttbb_compute_el->Clone();
  hNbOfSelectedJets_ttjj_compute_el->Reset();
  hNbOfSelectedJets_ttjj_compute_el->Add(hNbOfSelectedJets_ttbb_compute_el);
  hNbOfSelectedJets_ttjj_compute_el->Add(hNbOfSelectedJets_ttcc_compute_el);
  hNbOfSelectedJets_ttjj_compute_el->Add(hNbOfSelectedJets_ttll_compute_el);

  hNbOfSelectedJets_ttjj_compute_mu->SetMarkerColor(kRed);
  hNbOfSelectedJets_ttjj_compute_el->SetMarkerColor(kGreen);
  hNbOfSelectedJets_tttt_compute_mu->SetMarkerColor(kBlue);
  hNbOfSelectedJets_tttt_compute_el->SetMarkerColor(kBlack);

  hNbOfSelectedJets_ttjj_compute_mu->GetXaxis()->SetTitleOffset(0.85);
  hNbOfSelectedJets_ttjj_compute_mu->GetYaxis()->SetTitleOffset(0.85);

  hNbOfSelectedJets_ttjj_compute_mu->SetTitle("");
  hNbOfSelectedJets_ttjj_compute_mu->DrawNormalized();
  hNbOfSelectedJets_ttjj_compute_el->DrawNormalized("same");
  hNbOfSelectedJets_tttt_compute_mu->DrawNormalized("same");
  hNbOfSelectedJets_tttt_compute_el->DrawNormalized("same");
 
  TLegend *leg7 = new TLegend(0.6, 0.3, .8, .8);
  leg7->SetFillColor(0);
  leg7->AddEntry(hNbOfSelectedJets_ttjj_compute_mu,"ttjets mu","p");
  leg7->AddEntry(hNbOfSelectedJets_ttjj_compute_el,"ttjets el","p");
  leg7->AddEntry(hNbOfSelectedJets_tttt_compute_mu,"tttt mu","p");
  leg7->AddEntry(hNbOfSelectedJets_tttt_compute_el,"tttt el","p");
  leg7->Draw();





 TCanvas * c8 = new TCanvas();

  TH1F * hNbOfSelectedBJets_tttt_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_NbOfSelectedBJets/NbOfSelectedBJets_NP_overlay_TTTT");
  TH1F * hNbOfSelectedBJets_tttt_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_NbOfSelectedBJets/NbOfSelectedBJets_NP_overlay_TTTT");

  TH1F * hNbOfSelectedBJets_ttjj_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_NbOfSelectedBJets/NbOfSelectedBJets_TTJets");

 //  TH1F * hNbOfSelectedBJets_ttbb_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_NbOfSelectedBJets/NbOfSelectedBJets_TTJets_bb");
 //  TH1F * hNbOfSelectedBJets_ttcc_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_NbOfSelectedBJets/NbOfSelectedBJets_TTJets_cc");
 // TH1F * hNbOfSelectedBJets_ttll_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_NbOfSelectedBJets/NbOfSelectedBJets_TTJets_ll");

  TH1F * hNbOfSelectedBJets_ttbb_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_NbOfSelectedBJets/NbOfSelectedBJets_TTJets_bb");
  TH1F * hNbOfSelectedBJets_ttcc_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_NbOfSelectedBJets/NbOfSelectedBJets_TTJets_cc");
  TH1F * hNbOfSelectedBJets_ttll_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_NbOfSelectedBJets/NbOfSelectedBJets_TTJets_ll");

  //  TH1F * hNbOfSelectedBJets_ttjj_compute_mu = (TH1F*)hNbOfSelectedBJets_ttbb_compute_mu->Clone();
  // hNbOfSelectedBJets_ttjj_compute_mu->Reset();
  // hNbOfSelectedBJets_ttjj_compute_mu->Add(hNbOfSelectedBJets_ttbb_compute_mu);
  // hNbOfSelectedBJets_ttjj_compute_mu->Add(hNbOfSelectedBJets_ttcc_compute_mu);
  //hNbOfSelectedBJets_ttjj_compute_mu->Add(hNbOfSelectedBJets_ttll_compute_mu);


  TH1F * hNbOfSelectedBJets_ttjj_compute_el = (TH1F*)hNbOfSelectedBJets_ttbb_compute_el->Clone();
  hNbOfSelectedBJets_ttjj_compute_el->Reset();
  hNbOfSelectedBJets_ttjj_compute_el->Add(hNbOfSelectedBJets_ttbb_compute_el);
  hNbOfSelectedBJets_ttjj_compute_el->Add(hNbOfSelectedBJets_ttcc_compute_el);
  hNbOfSelectedBJets_ttjj_compute_el->Add(hNbOfSelectedBJets_ttll_compute_el);

  hNbOfSelectedBJets_ttjj_compute_mu->SetMarkerColor(kRed);
  hNbOfSelectedBJets_ttjj_compute_el->SetMarkerColor(kGreen);
  hNbOfSelectedBJets_tttt_compute_mu->SetMarkerColor(kBlue);
  hNbOfSelectedBJets_tttt_compute_el->SetMarkerColor(kBlack);

  hNbOfSelectedBJets_ttjj_compute_mu->GetXaxis()->SetTitleOffset(0.85);
  hNbOfSelectedBJets_ttjj_compute_mu->GetYaxis()->SetTitleOffset(0.85);

  hNbOfSelectedBJets_ttjj_compute_mu->SetTitle("");
  hNbOfSelectedBJets_ttjj_compute_mu->DrawNormalized();
  hNbOfSelectedBJets_ttjj_compute_el->DrawNormalized("same");
  hNbOfSelectedBJets_tttt_compute_mu->DrawNormalized("same");
  hNbOfSelectedBJets_tttt_compute_el->DrawNormalized("same");
 
  TLegend *leg8 = new TLegend(0.6, 0.3, .8, .8);
  leg8->SetFillColor(0);
  leg8->AddEntry(hNbOfSelectedBJets_ttjj_compute_mu,"ttjets mu","p");
  leg8->AddEntry(hNbOfSelectedBJets_ttjj_compute_el,"ttjets el","p");
  leg8->AddEntry(hNbOfSelectedBJets_tttt_compute_mu,"tttt mu","p");
  leg8->AddEntry(hNbOfSelectedBJets_tttt_compute_el,"tttt el","p");
  leg8->Draw();




 TCanvas * c9 = new TCanvas();

  TH1F * hSumJetMassX_tttt_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_SumJetMassX/SumJetMassX_NP_overlay_TTTT");
  TH1F * hSumJetMassX_tttt_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_SumJetMassX/SumJetMassX_NP_overlay_TTTT");


   TH1F * hSumJetMassX_ttjj_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_SumJetMassX/SumJetMassX_TTJets");

  // TH1F * hSumJetMassX_ttbb_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_SumJetMassX/SumJetMassX_TTJets_bb");
  // TH1F * hSumJetMassX_ttcc_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_SumJetMassX/SumJetMassX_TTJets_cc");
  //TH1F * hSumJetMassX_ttll_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_SumJetMassX/SumJetMassX_TTJets_ll");

  TH1F * hSumJetMassX_ttbb_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_SumJetMassX/SumJetMassX_TTJets_bb");
  TH1F * hSumJetMassX_ttcc_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_SumJetMassX/SumJetMassX_TTJets_cc");
  TH1F * hSumJetMassX_ttll_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_SumJetMassX/SumJetMassX_TTJets_ll");

  // TH1F * hSumJetMassX_ttjj_compute_mu = (TH1F*)hSumJetMassX_ttbb_compute_mu->Clone();
  // hSumJetMassX_ttjj_compute_mu->Reset();
  // hSumJetMassX_ttjj_compute_mu->Add(hSumJetMassX_ttbb_compute_mu);
  //hSumJetMassX_ttjj_compute_mu->Add(hSumJetMassX_ttcc_compute_mu);
  // hSumJetMassX_ttjj_compute_mu->Add(hSumJetMassX_ttll_compute_mu);


  TH1F * hSumJetMassX_ttjj_compute_el = (TH1F*)hSumJetMassX_ttbb_compute_el->Clone();
  hSumJetMassX_ttjj_compute_el->Reset();
  hSumJetMassX_ttjj_compute_el->Add(hSumJetMassX_ttbb_compute_el);
  hSumJetMassX_ttjj_compute_el->Add(hSumJetMassX_ttcc_compute_el);
  hSumJetMassX_ttjj_compute_el->Add(hSumJetMassX_ttll_compute_el);

  hSumJetMassX_ttjj_compute_mu->SetMarkerColor(kRed);
  hSumJetMassX_ttjj_compute_el->SetMarkerColor(kGreen);
  hSumJetMassX_tttt_compute_mu->SetMarkerColor(kBlue);
  hSumJetMassX_tttt_compute_el->SetMarkerColor(kBlack);

  hSumJetMassX_ttjj_compute_mu->GetXaxis()->SetTitleOffset(0.85); 
  hSumJetMassX_ttjj_compute_mu->GetYaxis()->SetTitleOffset(0.85);

  hSumJetMassX_ttjj_compute_mu->SetTitle("");
  hSumJetMassX_ttjj_compute_mu->DrawNormalized();
  hSumJetMassX_ttjj_compute_el->DrawNormalized("same");
  hSumJetMassX_tttt_compute_mu->DrawNormalized("same");
  hSumJetMassX_tttt_compute_el->DrawNormalized("same");
 
  TLegend *leg9 = new TLegend(0.6, 0.3, .8, .8);
  leg9->SetFillColor(0);
  leg9->AddEntry(hSumJetMassX_ttjj_compute_mu,"ttjets mu","p");
  leg9->AddEntry(hSumJetMassX_ttjj_compute_el,"ttjets el","p");
  leg9->AddEntry(hSumJetMassX_tttt_compute_mu,"tttt mu","p");
  leg9->AddEntry(hSumJetMassX_tttt_compute_el,"tttt el","p");
  leg9->Draw();






 TCanvas * c10 = new TCanvas();

  TH1F * h5thJetPt_tttt_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_5thJetPt/5thJetPt_NP_overlay_TTTT");
  TH1F * h5thJetPt_tttt_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_5thJetPt/5thJetPt_NP_overlay_TTTT");

  TH1F * h5thJetPt_ttjj_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_5thJetPt/5thJetPt_TTJets");

 //  TH1F * h5thJetPt_ttbb_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_5thJetPt/5thJetPt_TTJets_bb");
 // TH1F * h5thJetPt_ttcc_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_5thJetPt/5thJetPt_TTJets_cc");
 // TH1F * h5thJetPt_ttll_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_5thJetPt/5thJetPt_TTJets_ll");

  TH1F * h5thJetPt_ttbb_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_5thJetPt/5thJetPt_TTJets_bb");
  TH1F * h5thJetPt_ttcc_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_5thJetPt/5thJetPt_TTJets_cc");
  TH1F * h5thJetPt_ttll_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_5thJetPt/5thJetPt_TTJets_ll");

  //  TH1F * h5thJetPt_ttjj_compute_mu = (TH1F*)h5thJetPt_ttbb_compute_mu->Clone();
  // h5thJetPt_ttjj_compute_mu->Reset();
  //h5thJetPt_ttjj_compute_mu->Add(h5thJetPt_ttbb_compute_mu);
  //h5thJetPt_ttjj_compute_mu->Add(h5thJetPt_ttcc_compute_mu);
  //h5thJetPt_ttjj_compute_mu->Add(h5thJetPt_ttll_compute_mu);


  TH1F * h5thJetPt_ttjj_compute_el = (TH1F*)h5thJetPt_ttbb_compute_el->Clone();
  h5thJetPt_ttjj_compute_el->Reset();
  h5thJetPt_ttjj_compute_el->Add(h5thJetPt_ttbb_compute_el);
  h5thJetPt_ttjj_compute_el->Add(h5thJetPt_ttcc_compute_el);
  h5thJetPt_ttjj_compute_el->Add(h5thJetPt_ttll_compute_el);

  h5thJetPt_ttjj_compute_mu->SetMarkerColor(kRed);
  h5thJetPt_ttjj_compute_el->SetMarkerColor(kGreen);
  h5thJetPt_tttt_compute_mu->SetMarkerColor(kBlue);
  h5thJetPt_tttt_compute_el->SetMarkerColor(kBlack);

  h5thJetPt_ttjj_compute_mu->GetXaxis()->SetTitleOffset(0.85); 
  h5thJetPt_ttjj_compute_mu->GetYaxis()->SetTitleOffset(0.85);

  h5thJetPt_ttjj_compute_mu->SetTitle("");
  h5thJetPt_ttjj_compute_mu->DrawNormalized();
  h5thJetPt_ttjj_compute_el->DrawNormalized("same");
  h5thJetPt_tttt_compute_mu->DrawNormalized("same");
  h5thJetPt_tttt_compute_el->DrawNormalized("same");
 
  TLegend *leg10 = new TLegend(0.6, 0.3, .8, .8);
  leg10->SetFillColor(0);
  leg10->AddEntry(h5thJetPt_ttjj_compute_mu,"ttjets mu","p");
  leg10->AddEntry(h5thJetPt_ttjj_compute_el,"ttjets el","p");
  leg10->AddEntry(h5thJetPt_tttt_compute_mu,"tttt mu","p");
  leg10->AddEntry(h5thJetPt_tttt_compute_el,"tttt el","p");
  leg10->Draw();




  TCanvas * c11 = new TCanvas();

  TH1F * h6thJetPt_tttt_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_6thJetPt/6thJetPt_NP_overlay_TTTT");
  TH1F * h6thJetPt_tttt_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_6thJetPt/6thJetPt_NP_overlay_TTTT");


  TH1F * h6thJetPt_ttjj_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_6thJetPt/6thJetPt_TTJets");

  //  TH1F * h6thJetPt_ttbb_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_6thJetPt/6thJetPt_TTJets_bb");
  //TH1F * h6thJetPt_ttcc_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_6thJetPt/6thJetPt_TTJets_cc");
  //TH1F * h6thJetPt_ttll_compute_mu = (TH1F*)f_compute_mu->Get("MultiSamplePlot_6thJetPt/6thJetPt_TTJets_ll");

  TH1F * h6thJetPt_ttbb_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_6thJetPt/6thJetPt_TTJets_bb");
  TH1F * h6thJetPt_ttcc_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_6thJetPt/6thJetPt_TTJets_cc");
  TH1F * h6thJetPt_ttll_compute_el = (TH1F*)f_compute_el->Get("MultiSamplePlot_6thJetPt/6thJetPt_TTJets_ll");

  //  TH1F * h6thJetPt_ttjj_compute_mu = (TH1F*)h6thJetPt_ttbb_compute_mu->Clone();
  // h6thJetPt_ttjj_compute_mu->Reset();
  // h6thJetPt_ttjj_compute_mu->Add(h6thJetPt_ttbb_compute_mu);
  //h6thJetPt_ttjj_compute_mu->Add(h6thJetPt_ttcc_compute_mu);
  // h6thJetPt_ttjj_compute_mu->Add(h6thJetPt_ttll_compute_mu);


  TH1F * h6thJetPt_ttjj_compute_el = (TH1F*)h6thJetPt_ttbb_compute_el->Clone();
  h6thJetPt_ttjj_compute_el->Reset();
  h6thJetPt_ttjj_compute_el->Add(h6thJetPt_ttbb_compute_el);
  h6thJetPt_ttjj_compute_el->Add(h6thJetPt_ttcc_compute_el);
  h6thJetPt_ttjj_compute_el->Add(h6thJetPt_ttll_compute_el);

  h6thJetPt_ttjj_compute_mu->SetMarkerColor(kRed);
  h6thJetPt_ttjj_compute_el->SetMarkerColor(kGreen);
  h6thJetPt_tttt_compute_mu->SetMarkerColor(kBlue);
  h6thJetPt_tttt_compute_el->SetMarkerColor(kBlack);

  h6thJetPt_ttjj_compute_mu->GetXaxis()->SetTitleOffset(0.85); 
  h6thJetPt_ttjj_compute_mu->GetYaxis()->SetTitleOffset(0.85);

  h6thJetPt_ttjj_compute_mu->SetTitle("");
  h6thJetPt_ttjj_compute_mu->DrawNormalized();
  h6thJetPt_ttjj_compute_el->DrawNormalized("same");
  h6thJetPt_tttt_compute_mu->DrawNormalized("same");
  h6thJetPt_tttt_compute_el->DrawNormalized("same");
 
  TLegend *leg11 = new TLegend(0.6, 0.3, .8, .8);
  leg11->SetFillColor(0);
  leg11->AddEntry(h6thJetPt_ttjj_compute_mu,"ttjets mu","p");
  leg11->AddEntry(h6thJetPt_ttjj_compute_el,"ttjets el","p");
  leg11->AddEntry(h6thJetPt_tttt_compute_mu,"tttt mu","p");
  leg11->AddEntry(h6thJetPt_tttt_compute_el,"tttt el","p");
  leg11->Draw();

   std::cout <<"here 0"<< std::endl;



  //Compasion 2: correlations of BDT variables in data and simulation via profile plots
  vector<string> datasets; 
  vector<string> plots; 

  datasets.push_back("Data");
  datasets.push_back("TTJets");
  datasets.push_back("NP_overlay_TTTT");
  
  plots.push_back("HTX_vs_HTH_");
  plots.push_back("HTX_vs_MultiTopness_");
  plots.push_back("HTX_vs_5thJetPt_");
  plots.push_back("HTX_vs_6thJetPt_");
  plots.push_back("HTX_vs_HTb_");
  plots.push_back("HTX_vs_HTRat_");
  plots.push_back("HTX_vs_NbOfSelectedJets_");
  plots.push_back("HTX_vs_NbOfSelectedBJets_");
  plots.push_back("HTX_vs_SumJetMassX_");
  plots.push_back("MultiTopness_vs_HTH_");
  plots.push_back("MultiTopness_vs_5thJetPt_");
  plots.push_back("MultiTopness_vs_6thJetPt_");
  plots.push_back("MultiTopness_vs_HTb_");
  plots.push_back("MultiTopness_vs_HTRat_");
  plots.push_back("MultiTopness_vs_NbOfSelectedJets_");
  plots.push_back("MultiTopness_vs_NbOfSelectedBJets_");
  plots.push_back("MultiTopness_vs_SumJetMassX_");   
  plots.push_back("HTH_vs_5thJetPt_");
  plots.push_back("HTH_vs_6thJetPt_");
  plots.push_back("HTH_vs_HTb_");
  plots.push_back("HTH_vs_HTRat_");
  plots.push_back("HTH_vs_NbOfSelectedJets_");
  plots.push_back("HTH_vs_NbOfSelectedBJets_");
  plots.push_back("HTH_vs_SumJetMassX_");
  plots.push_back("5thJetPt_vs_6thJetPt_");
  plots.push_back("5thJetPt_vs_HTb_");
  plots.push_back("5thJetPt_vs_HTRat_");
  plots.push_back("5thJetPt_vs_NbOfSelectedJets_");
  plots.push_back("5thJetPt_vs_NbOfSelectedBJets_");
  plots.push_back("5thJetPt_vs_SumJetMassX_");
  plots.push_back("6thJetPt_vs_HTb_");
  plots.push_back("6thJetPt_vs_HTRat_");
  plots.push_back("6thJetPt_vs_NbOfSelectedJets_");
  plots.push_back("6thJetPt_vs_NbOfSelectedBJets_");
  plots.push_back("6thJetPt_vs_SumJetMassX_");
  plots.push_back("HTb_vs_HTRat_");
  plots.push_back("HTb_vs_NbOfSelectedJets_");
  plots.push_back("HTb_vs_NbOfSelectedBJets_");
  plots.push_back("HTb_vs_SumJetMassX_");
  plots.push_back("HTRat_vs_NbOfSelectedJets_");
  plots.push_back("HTRat_vs_NbOfSelectedBJets_");
  plots.push_back("HTRat_vs_SumJetMassX_");
  plots.push_back("NbOfSelectedJets_vs_NbOfSelectedBJets_");
  plots.push_back("NbOfSelectedJets_vs_SumJetMassX_");
  plots.push_back("NbOfSelectedBJets_vs_SumJetMassX_");
  
 

 bool first_plot = false;


for (unsigned int p = 0; p < plots.size(); p++){
  TCanvas * ctemp = new TCanvas();

 for (unsigned int d = 0; d < datasets.size(); d++){

 histoProfile[(plots[p]+datasets[d]).c_str()] =  (TProfile*)f_compute_mu->Get(("HistosProfile/"+plots[p]+datasets[d]).c_str());
   cout <<"here 4 "<<  (plots[p]+datasets[d]).c_str()  <<endl;


   if (datasets[d] == "Data")
{
  histoProfile[(plots[p]+datasets[d]).c_str()]->SetMarkerColor(kBlack);
 }else if (datasets[d] == "TTJets"){

 histoProfile[(plots[p]+datasets[d]).c_str()]->SetMarkerColor(kRed);

   }else if (datasets[d] == "NP_overlay_TTTT") {
 histoProfile[(plots[p]+datasets[d]).c_str()]->SetMarkerColor(kGray);

}

   if(!first_plot){
   histoProfile[(plots[p]+datasets[d]).c_str()]->Draw();
   }else{

   histoProfile[(plots[p]+datasets[d]).c_str()]->Draw("same");
}


   if(!first_plot) first_plot = true;

}


 TLegend *templeg = new TLegend(0.8, 0.1, 0.9, 0.4);
  templeg->SetFillColor(0);
  templeg->AddEntry(histoProfile[(plots[p]+"Data").c_str()],"Data","p");
  templeg->AddEntry(histoProfile[(plots[p]+"TTJets").c_str()],"tt + jets","p");
  templeg->AddEntry(histoProfile[(plots[p]+"NP_overlay_TTTT").c_str()],"tttt","p");
  templeg->Draw();



 string canvas_name = "CorrelationPlots/"+plots[p]+"canvas.pdf";

 ctemp->SaveAs(canvas_name.c_str());

 delete ctemp, templeg;

first_plot = false;

 }


  c1->SaveAs("BDT_compare.pdf");
  c2->SaveAs("HTb_compare.pdf");
  c3->SaveAs("MultiTopness_compare.pdf");
  c4->SaveAs("HTX_compare.pdf");
  c5->SaveAs("HTH_compare.pdf");
  c6->SaveAs("HTRat_compare.pdf");
  c7->SaveAs("NbOfSelectedJets_compare.pdf");
  c8->SaveAs("NbOfSelectedBJets_compare.pdf");
  c9->SaveAs("SumJetMassX_compare.pdf");
  c10->SaveAs("5thJetPt_compare.pdf");
  c11->SaveAs("6thJetPt_compare.pdf");



}
