#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>
#include <fstream>

void process(){

 /////////////////////////
 // F I L E S
 /////////////////////////
 TFile* f_Mu = new TFile("FourTop_EventSelection_wMETCut_Mu.root"); 
 TFile* f1 = new TFile("SystematicShapes_Mu.root");


////////////////////////////////////////////////
// P R O C E S S   N O M I N A L   S H A P E S
////////////////////////////////////////////////

 TFile* fnom = new TFile("NominalShapes_Mu.root","RECREATE"); //file for writing


/////////////////////////
//  I N C L U S I V E 
/////////////////////////

  TH1F * hData_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_Data");
  TH1F * hSig_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_NP_overlay_TTTT");

  double sig_scale = 0.92*1.3;//correcting for (very slightly incorrect int. lumi of data histo and correcting to NLO xsection (doesn't make any real difference))
  hSig_Mu->Scale(sig_scale);

  //signal injection test
  TH1F * hdata_100fbtttt = new TH1F();
  TH1F * hdata_150fbtttt = new TH1F();
  TH1F * hdata_200fbtttt = new TH1F();
  hdata_100fbtttt = (TH1F*)hData_Mu->Clone();
  hdata_150fbtttt = (TH1F*)hData_Mu->Clone();
  hdata_200fbtttt = (TH1F*)hData_Mu->Clone();

  hdata_100fbtttt->Add(hSig_Mu, 100/1.3);
  hdata_150fbtttt->Add(hSig_Mu, 150/1.3);
  hdata_200fbtttt->Add(hSig_Mu, 200/1.3);

  hdata_100fbtttt->SetMarkerColor(kBlue);
  hdata_150fbtttt->SetMarkerColor(kGreen);
  hdata_200fbtttt->SetMarkerColor(kRed);
  hSig_Mu->SetLineColor(kOrange);
  hSig_Mu->SetLineWidth(2); 
  hData_Mu->SetLineWidth(2);
  hData_Mu->SetTitle("");
  hData_Mu->SetMinimum(0.001);
  hData_Mu->SetMaximum(100000);

  TCanvas * cx = new TCanvas();
  hData_Mu->Draw("HIST");
  hSig_Mu->Draw("HISTSAME");
  hdata_100fbtttt->Draw("same");
  hdata_150fbtttt->Draw("same");
  hdata_200fbtttt->Draw("same");
  cx->SetLogy();
  cx->Update();

  TLegend* leg = new TLegend(0.6, 0.73, .97, .95);
  leg->SetFillColor(0);
  leg->AddEntry(hData_Mu, "Data", "l");
  leg->AddEntry(hdata_100fbtttt, "Data with signal injection (#sigma = 100 fb)", "p");
  leg->AddEntry(hdata_150fbtttt, "Data with signal injection (#sigma = 150 fb)", "p");
  leg->AddEntry(hdata_200fbtttt, "Data with signal injection (#sigma = 200 fb)", "p");
  leg->AddEntry(hSig_Mu, "Signal", "l");

  leg->Draw();
  cx->SaveAs("SignalInjection_Mu.pdf");

  cout <<"sig injection plot saved..."<<endl;
  TH1F * hTTll_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_TTJets_ll");
  TH1F * hTTcc_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_TTJets_cc");
  TH1F * hTTbb_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_TTJets_bb");
  TH1F * hTT_Mu   = (TH1F*)hTTbb_Mu->Clone();
  hTT_Mu->Reset();
  hTT_Mu->Add(hTTll_Mu);
  hTT_Mu->Add(hTTcc_Mu);
  hTT_Mu->Add(hTTbb_Mu);

  double nExTT_Mu  = hTT_Mu->Integral();
  /*
  hTT_Mu->Write("tt");
  hSig_Mu->Write("tttt");
  hData_Mu->Write("data");
  hdata_100fbtttt->Write("data_w100fbtttt");
  hdata_150fbtttt->Write("data_w150fbtttt");
  hdata_200fbtttt->Write("data_w200fbtttt");

  */
  //  cout << hTTll_Mu->Integral() <<endl;
  //cout << hTTcc_Mu->Integral() <<endl;
  // cout << hTTbb_Mu->Integral() <<endl;
   cout << nExTT_Mu<<endl;

  cout <<"nominal histos 6j..."<<endl;
 /////////////////////////
 // Nominal Histos
 /////////////////////////

  //EW histo
  TH1F * hW4_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_W_4Jets");
  TH1F * hZ4_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_Z_4Jets");
  //  TH1F * hZ4_lm_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_Z_4Jets_lmets");
  TH1F * hSingleTop_t_T_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_SingleTop_t_T");
  TH1F * hSingleTop_t_TBar_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_SingleTop_t_TBar");
  TH1F * hSingleTop_s_T_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_SingleTop_s_T");
  TH1F * hSingleTop_s_TBar_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_SingleTop_s_TBar");
  TH1F * hSingleTop_tW_T_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_SingleTop_tW_T");
  TH1F * hSingleTop_tW_TBar_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_SingleTop_tW_TBar");
  TH1F * hWW_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_WW");
  TH1F * hWZ_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_WZ");
  TH1F * hZZ_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_ZZ");

  TH1F * hEW_Mu = (TH1F*)hW4_Mu->Clone();
  hEW_Mu->Reset();

  hEW_Mu->Add(hW4_Mu);
  hEW_Mu->Add(hZ4_Mu);
  // hEW_Mu->Add(hZ4_lm_Mu);
  hEW_Mu->Add(hSingleTop_t_T_Mu);
  hEW_Mu->Add(hSingleTop_t_TBar_Mu);
  hEW_Mu->Add(hSingleTop_s_T_Mu);
  hEW_Mu->Add(hSingleTop_s_TBar_Mu);
  hEW_Mu->Add(hSingleTop_tW_T_Mu);
  hEW_Mu->Add(hSingleTop_tW_TBar_Mu);
  hEW_Mu->Add(hWW_Mu);
  hEW_Mu->Add(hWZ_Mu);
  hEW_Mu->Add(hZZ_Mu);

  //  hEW_Mu->Write("EW");

  //tt other histo
  TH1F * httW_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_ttW");
  TH1F * httZ_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_ttZ");
  TH1F * httH_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_ttH");
  TH1F * hTTOther_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_TTJets_Other");

  hTTOther_Mu->Add(httW_Mu);
  hTTOther_Mu->Add(httZ_Mu);
  hTTOther_Mu->Add(httH_Mu);

  //  hTTOther_Mu->Write("ttOther");



 /////////////////////////
 // S I X   J E T  B I N
 /////////////////////////

  TH1F * hData_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_Data");
  TH1F * hSig_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_NP_overlay_TTTT");

  sig_scale = 0.92*1.3;//correcting for (very slightly incorrect int. lumi of data histo and correcting to NLO xsection (doesn't make any real difference))
  hSig_Mu_6j->Scale(sig_scale);

  //signal injection test
  TH1F * hdata_100fbtttt_6j = new TH1F();
  TH1F * hdata_150fbtttt_6j = new TH1F();
  TH1F * hdata_200fbtttt_6j = new TH1F();
  hdata_100fbtttt_6j = (TH1F*)hData_Mu_6j->Clone();
  hdata_150fbtttt_6j = (TH1F*)hData_Mu_6j->Clone();
  hdata_200fbtttt_6j = (TH1F*)hData_Mu_6j->Clone();

  hdata_100fbtttt_6j->Add(hSig_Mu_6j, 100/1.3);
  hdata_150fbtttt_6j->Add(hSig_Mu_6j, 150/1.3);
  hdata_200fbtttt_6j->Add(hSig_Mu_6j, 200/1.3);

  hdata_100fbtttt_6j->SetMarkerColor(kBlue);
  hdata_150fbtttt_6j->SetMarkerColor(kGreen);
  hdata_200fbtttt_6j->SetMarkerColor(kRed);
  hSig_Mu_6j->SetLineColor(kOrange);
  hSig_Mu_6j->SetLineWidth(2); 
  hData_Mu_6j->SetLineWidth(2);
  hData_Mu_6j->SetTitle("");
  hData_Mu_6j->SetMinimum(0.001);
  hData_Mu_6j->SetMaximum(100000);

  TCanvas * cx_6j = new TCanvas();
  hData_Mu_6j->Draw("HIST");
  hSig_Mu_6j->Draw("HISTSAME");
  hdata_100fbtttt_6j->Draw("same");
  hdata_150fbtttt_6j->Draw("same");
  hdata_200fbtttt_6j->Draw("same");
  cx_6j->SetLogy();
  cx_6j->Update();

  TLegend* leg_6j = new TLegend(0.6, 0.73, .97, .95);
  leg_6j->SetFillColor(0);
  leg_6j->AddEntry(hData_Mu_6j, "Data", "l");
  leg_6j->AddEntry(hdata_100fbtttt_6j, "Data with signal injection (#sigma = 100 fb)", "p");
  leg_6j->AddEntry(hdata_150fbtttt_6j, "Data with signal injection (#sigma = 150 fb)", "p");
  leg_6j->AddEntry(hdata_200fbtttt_6j, "Data with signal injection (#sigma = 200 fb)", "p");
  leg_6j->AddEntry(hSig_Mu_6j, "Signal", "l");

  leg_6j->Draw();
  cx_6j->SaveAs("SignalInjection_Mu.pdf");

  cout <<"sig injection plot saved..."<<endl;
  TH1F * hTTll_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_TTJets_ll");
  TH1F * hTTcc_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_TTJets_cc");
  TH1F * hTTbb_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_TTJets_bb");
  TH1F * hTT_Mu_6j   = (TH1F*)hTTbb_Mu_6j->Clone();
  hTT_Mu_6j->Reset();
  hTT_Mu_6j->Add(hTTll_Mu_6j);
  hTT_Mu_6j->Add(hTTcc_Mu_6j);
  hTT_Mu_6j->Add(hTTbb_Mu_6j);

  double nExTT_Mu_6j  = hTT_Mu_6j->Integral();
  /*
  hTT_Mu_6j->Write("tt_6j");
  hSig_Mu_6j->Write("tttt_6j");
  hData_Mu_6j->Write("data_6j");
  hdata_100fbtttt_6j->Write("data_w100fbtttt_6j");
  hdata_150fbtttt_6j->Write("data_w150fbtttt_6j");
  hdata_200fbtttt_6j->Write("data_w200fbtttt_6j");


  */
  //  cout << hTTll_Mu->Integral() <<endl;
  //cout << hTTcc_Mu->Integral() <<endl;
  // cout << hTTbb_Mu->Integral() <<endl;
  //  cout << nExTT_Mu<<endl;

  cout <<"nominal histos 6j..."<<endl;
 /////////////////////////
 // Nominal Histos
 /////////////////////////

  //EW histo
  TH1F * hW4_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_W_4Jets");
  TH1F * hZ4_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_Z_4Jets");
  //  TH1F * hZ4_lm_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_Z_4Jets_lm_6Jets");
  TH1F * hSingleTop_t_T_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_SingleTop_t_T");
  TH1F * hSingleTop_t_TBar_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_SingleTop_t_TBar");
  TH1F * hSingleTop_s_T_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_SingleTop_s_T");
  TH1F * hSingleTop_s_TBar_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_SingleTop_s_TBar");
  TH1F * hSingleTop_tW_T_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_SingleTop_tW_T");
  TH1F * hSingleTop_tW_TBar_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_SingleTop_tW_TBar");
  TH1F * hWW_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_WW");
  TH1F * hWZ_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_WZ");
  TH1F * hZZ_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_ZZ");

  TH1F * hEW_Mu_6j = (TH1F*)hW4_Mu_6j->Clone();
  hEW_Mu_6j->Reset();

  hEW_Mu_6j->Add(hW4_Mu_6j);
  hEW_Mu_6j->Add(hZ4_Mu_6j);
  // hEW_Mu_6j->Add(hZ4_lm_Mu_6j);
  hEW_Mu_6j->Add(hSingleTop_t_T_Mu_6j);
  hEW_Mu_6j->Add(hSingleTop_t_TBar_Mu_6j);
  hEW_Mu_6j->Add(hSingleTop_s_T_Mu_6j);
  hEW_Mu_6j->Add(hSingleTop_s_TBar_Mu_6j);
  hEW_Mu_6j->Add(hSingleTop_tW_T_Mu_6j);
  hEW_Mu_6j->Add(hSingleTop_tW_TBar_Mu_6j);
  hEW_Mu_6j->Add(hWW_Mu_6j);
  hEW_Mu_6j->Add(hWZ_Mu_6j);
  hEW_Mu_6j->Add(hZZ_Mu_6j);

  //  hEW_Mu_6j->Write("EW_6j");

  //tt other histo
  TH1F * httW_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_ttW");
  TH1F * httZ_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_ttZ");
  TH1F * httH_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_ttH");
  TH1F * hTTOther_Mu_6j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_TTJets_Other");

  hTTOther_Mu_6j->Add(httW_Mu_6j);
  hTTOther_Mu_6j->Add(httZ_Mu_6j);
  hTTOther_Mu_6j->Add(httH_Mu_6j);

  //  hTTOther_Mu_6j->Write("ttOther_6j");
 
 /////////////////////////
 // S E V E N   J E T  B I N
 /////////////////////////
  cout <<"7j..."<<endl;
  TH1F * hData_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_Data");
  TH1F * hSig_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_NP_overlay_TTTT");

  hSig_Mu_7j->Scale(sig_scale);

  //signal injection test
  TH1F * hdata_100fbtttt_7j = new TH1F();
  TH1F * hdata_150fbtttt_7j = new TH1F();
  TH1F * hdata_200fbtttt_7j = new TH1F();
  hdata_100fbtttt_7j = (TH1F*)hData_Mu_7j->Clone();
  hdata_150fbtttt_7j = (TH1F*)hData_Mu_7j->Clone();
  hdata_200fbtttt_7j = (TH1F*)hData_Mu_7j->Clone();

  hdata_100fbtttt_7j->Add(hSig_Mu_7j, 100/1.3);
  hdata_150fbtttt_7j->Add(hSig_Mu_7j, 150/1.3);
  hdata_200fbtttt_7j->Add(hSig_Mu_7j, 200/1.3);

  hdata_100fbtttt_7j->SetMarkerColor(kBlue);
  hdata_150fbtttt_7j->SetMarkerColor(kGreen);
  hdata_200fbtttt_7j->SetMarkerColor(kRed);
  hSig_Mu_7j->SetLineColor(kOrange);
  hSig_Mu_7j->SetLineWidth(2); 
  hData_Mu_7j->SetLineWidth(2);
  hData_Mu_7j->SetTitle("");
  hData_Mu_7j->SetMinimum(0.001);
  hData_Mu_7j->SetMaximum(100000);

  TCanvas * cx_7j = new TCanvas();
  hData_Mu_7j->Draw("HIST");
  hSig_Mu_7j->Draw("HISTSAME");
  hdata_100fbtttt_7j->Draw("same");
  hdata_150fbtttt_7j->Draw("same");
  hdata_200fbtttt_7j->Draw("same");
  cx_7j->SetLogy();
  cx_7j->Update();

  TLegend* leg_7j = new TLegend(0.6, 0.73, .97, .95);
  leg_7j->SetFillColor(0);
  leg_7j->AddEntry(hData_Mu_7j, "Data", "l");
  leg_7j->AddEntry(hdata_100fbtttt_7j, "Data with signal injection (#sigma = 100 fb)", "p");
  leg_7j->AddEntry(hdata_150fbtttt_7j, "Data with signal injection (#sigma = 150 fb)", "p");
  leg_7j->AddEntry(hdata_200fbtttt_7j, "Data with signal injection (#sigma = 200 fb)", "p");
  leg_7j->AddEntry(hSig_Mu_7j, "Signal", "l");

  leg_7j->Draw();
  cx_7j->SaveAs("SignalInjection_Mu_7j.pdf");

  TH1F * hTTll_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_TTJets_ll");
  TH1F * hTTcc_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_TTJets_cc");
  TH1F * hTTbb_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_TTJets_bb");
  TH1F * hTT_Mu_7j   = (TH1F*)hTTbb_Mu_7j->Clone();
  hTT_Mu_7j->Reset();
  hTT_Mu_7j->Add(hTTll_Mu_7j);
  hTT_Mu_7j->Add(hTTcc_Mu_7j);
  hTT_Mu_7j->Add(hTTbb_Mu_7j);

  double nExTT_Mu_7j  = hTT_Mu_7j->Integral();
  /*
  hTT_Mu_7j->Write("tt_7j");
  hSig_Mu_7j->Write("tttt_7j");
  hData_Mu_7j->Write("data_7j");
  hdata_100fbtttt_7j->Write("data_w100fbtttt_7j");
  hdata_150fbtttt_7j->Write("data_w150fbtttt_7j");
  hdata_200fbtttt_7j->Write("data_w200fbtttt_7j");
  */
  //  cout << hTTll_Mu->Integral() <<endl;
  //cout << hTTcc_Mu->Integral() <<endl;
  // cout << hTTbb_Mu->Integral() <<endl;
  //  cout << nExTT_Mu_7j<<endl;

 cout <<"nominal 6j..."<<endl;
 /////////////////////////
 // Nominal Histos
 /////////////////////////

  //EW histo
  TH1F * hW4_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_W_4Jets");
  TH1F * hZ4_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_Z_4Jets");
  //  TH1F * hZ4_lm_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_Z_4Jets_lm");
  TH1F * hSingleTop_t_T_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_SingleTop_t_T");
  TH1F * hSingleTop_t_TBar_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_SingleTop_t_TBar");
  TH1F * hSingleTop_s_T_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_SingleTop_s_T");
  TH1F * hSingleTop_s_TBar_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_SingleTop_s_TBar");
  TH1F * hSingleTop_tW_T_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_SingleTop_tW_T");
  TH1F * hSingleTop_tW_TBar_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_SingleTop_tW_TBar");
  TH1F * hWW_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_WW");
  TH1F * hWZ_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_WZ");
  TH1F * hZZ_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_ZZ");

  TH1F * hEW_Mu_7j = (TH1F*)hW4_Mu_7j->Clone();
  hEW_Mu_7j->Reset();

  hEW_Mu_7j->Add(hW4_Mu_7j);
  hEW_Mu_7j->Add(hZ4_Mu_7j);
  // hEW_Mu_7j->Add(hZ4_lm_Mu_7j);
  hEW_Mu_7j->Add(hSingleTop_t_T_Mu_7j);
  hEW_Mu_7j->Add(hSingleTop_t_TBar_Mu_7j);
  hEW_Mu_7j->Add(hSingleTop_s_T_Mu_7j);
  hEW_Mu_7j->Add(hSingleTop_s_TBar_Mu_7j);
  hEW_Mu_7j->Add(hSingleTop_tW_T_Mu_7j);
  hEW_Mu_7j->Add(hSingleTop_tW_TBar_Mu_7j);
  hEW_Mu_7j->Add(hWW_Mu_7j);
  hEW_Mu_7j->Add(hWZ_Mu_7j);
  hEW_Mu_7j->Add(hZZ_Mu_7j);

  //  hEW_Mu_7j->Write("EW_7j");

  //tt other histo
  TH1F * httW_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_ttW");
  TH1F * httZ_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_ttZ");
  TH1F * httH_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_ttH");
  TH1F * hTTOther_Mu_7j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_TTJets_Other");

  //  TH1F * hTTOther_Mu = httW_Mu->Clone();
  //hTTOther_Mu->Reset();

  hTTOther_Mu_7j->Add(httW_Mu_7j);
  hTTOther_Mu_7j->Add(httZ_Mu_7j);
  hTTOther_Mu_7j->Add(httH_Mu_7j);

  //  hTTOther_Mu_7j->Write("ttOther_7j");

 /////////////////////////
 // E I G H T   J E T  B I N
 /////////////////////////

 cout <<"8j..."<<endl;


  TH1F * hData_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_Data");
  TH1F * hSig_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_NP_overlay_TTTT");

  hSig_Mu_8j->Scale(sig_scale);

  //signal injection test
  TH1F * hdata_100fbtttt_8j = new TH1F();
  TH1F * hdata_150fbtttt_8j = new TH1F();
  TH1F * hdata_200fbtttt_8j = new TH1F();
  hdata_100fbtttt_8j = (TH1F*)hData_Mu_8j->Clone();
  hdata_150fbtttt_8j = (TH1F*)hData_Mu_8j->Clone();
  hdata_200fbtttt_8j = (TH1F*)hData_Mu_8j->Clone();

  hdata_100fbtttt_8j->Add(hSig_Mu_8j, 100/1.3);
  hdata_150fbtttt_8j->Add(hSig_Mu_8j, 150/1.3);
  hdata_200fbtttt_8j->Add(hSig_Mu_8j, 200/1.3);


  hdata_100fbtttt_8j->SetMarkerColor(kBlue);
  hdata_150fbtttt_8j->SetMarkerColor(kGreen);
  hdata_200fbtttt_8j->SetMarkerColor(kRed);
  hSig_Mu_8j->SetLineColor(kOrange);
  hSig_Mu_8j->SetLineWidth(2); 
  hData_Mu_8j->SetLineWidth(2);
  hData_Mu_8j->SetTitle("");
  hData_Mu_8j->SetMinimum(0.001);
  hData_Mu_8j->SetMaximum(100000);

  TCanvas * cx_8j = new TCanvas();
  hData_Mu_8j->Draw("HIST");
  hSig_Mu_8j->Draw("HISTSAME");
  hdata_100fbtttt_8j->Draw("same");
  hdata_150fbtttt_8j->Draw("same");
  hdata_200fbtttt_8j->Draw("same");
  cx_8j->SetLogy();
  cx_8j->Update();

  TLegend* leg_8j = new TLegend(0.6, 0.73, .97, .95);
  leg_8j->SetFillColor(0);
  leg_8j->AddEntry(hData_Mu_8j, "Data", "l");
  leg_8j->AddEntry(hdata_100fbtttt_8j, "Data with signal injection (#sigma = 100 fb)", "p");
  leg_8j->AddEntry(hdata_150fbtttt_8j, "Data with signal injection (#sigma = 150 fb)", "p");
  leg_8j->AddEntry(hdata_200fbtttt_8j, "Data with signal injection (#sigma = 200 fb)", "p");
  leg_8j->AddEntry(hSig_Mu_8j, "Signal", "l");
  leg_8j->Draw();
  cx_8j->SaveAs("SignalInjection_Mu_8j.pdf");


  TH1F * hTTll_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_TTJets_ll");
  TH1F * hTTcc_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_TTJets_cc");
  TH1F * hTTbb_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_TTJets_bb");
  TH1F * hTT_Mu_8j   = (TH1F*)hTTbb_Mu_8j->Clone();
  hTT_Mu_8j->Reset();
  hTT_Mu_8j->Add(hTTll_Mu_8j);
  hTT_Mu_8j->Add(hTTcc_Mu_8j);
  hTT_Mu_8j->Add(hTTbb_Mu_8j);

  double nExTT_Mu_8j  = hTT_Mu_8j->Integral();
  /*
  hTT_Mu_8j->Write("tt_8j");
  hSig_Mu_8j->Write("tttt_8j");
  hData_Mu_8j->Write("data_8j");
  hdata_100fbtttt_8j->Write("data_w100fbtttt_8j");
  hdata_150fbtttt_8j->Write("data_w150fbtttt_8j");
  hdata_200fbtttt_8j->Write("data_w200fbtttt_8j");
  */
  //  cout << hTTll_Mu->Integral() <<endl;
  //cout << hTTcc_Mu->Integral() <<endl;
  // cout << hTTbb_Mu->Integral() <<endl;
  //  cout << nExTT_Mu_8j<<endl;


 /////////////////////////
 // Nominal Histos
 /////////////////////////
 cout <<"nominal 8j..."<<endl;
  //EW histo
  TH1F * hW4_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_W_4Jets");
  TH1F * hZ4_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_Z_4Jets");
  //  TH1F * hZ4_lm_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_Z_4Jets_lm");
  TH1F * hSingleTop_t_T_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_SingleTop_t_T");
  TH1F * hSingleTop_t_TBar_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_SingleTop_t_TBar");
  TH1F * hSingleTop_s_T_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_SingleTop_s_T");
  TH1F * hSingleTop_s_TBar_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_SingleTop_s_TBar");
  TH1F * hSingleTop_tW_T_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_SingleTop_tW_T");
  TH1F * hSingleTop_tW_TBar_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_SingleTop_tW_TBar");
  TH1F * hWW_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_WW");
  TH1F * hWZ_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_WZ");
  TH1F * hZZ_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_ZZ");

  TH1F * hEW_Mu_8j = (TH1F*)hW4_Mu_8j->Clone();
  hEW_Mu_8j->Reset();

  hEW_Mu_8j->Add(hW4_Mu_8j);
  hEW_Mu_8j->Add(hZ4_Mu_8j);
  // hEW_Mu_8j->Add(hZ4_lm_Mu_8j);
  hEW_Mu_8j->Add(hSingleTop_t_T_Mu_8j);
  hEW_Mu_8j->Add(hSingleTop_t_TBar_Mu_8j);
  hEW_Mu_8j->Add(hSingleTop_s_T_Mu_8j);
  hEW_Mu_8j->Add(hSingleTop_s_TBar_Mu_8j);
  hEW_Mu_8j->Add(hSingleTop_tW_T_Mu_8j);
  hEW_Mu_8j->Add(hSingleTop_tW_TBar_Mu_8j);
  hEW_Mu_8j->Add(hWW_Mu_8j);
  hEW_Mu_8j->Add(hWZ_Mu_8j);
  hEW_Mu_8j->Add(hZZ_Mu_8j);

  //  hEW_Mu_8j->Write("EW_8j");

  //tt other histo
  TH1F * httW_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_ttW");
  TH1F * httZ_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_ttZ");
  TH1F * httH_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_ttH");
  TH1F * hTTOther_Mu_8j = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_TTJets_Other");

  //  TH1F * hTTOther_Mu = httW_Mu->Clone();
  //hTTOther_Mu->Reset();

  hTTOther_Mu_8j->Add(httW_Mu_8j);
  hTTOther_Mu_8j->Add(httZ_Mu_8j);
  hTTOther_Mu_8j->Add(httH_Mu_8j);

  //  hTTOther_Mu_8j->Write("ttOther_8j");
  
  fnom->Close();

  ///////////////////////////////////////////////////////
  //// P R O C E S S   S Y S T E M A T I C   H I S T O S
  ///////////////////////////////////////////////////////

  cout <<"processing sys shapes..."<<endl;

 TFile* f2 = new TFile("SystematicShapes_norm_Mu.root", "RECREATE");

 ///////////////////////////////
 /// I N C L U S I V E 
 /////////////////////////////////

ofstream myfile;
  myfile.open ("chi2_systematics_mu.txt");
  myfile << "Chi2/NDOF of systematics:"<< endl;


  TH1F *  h_Scale_Down = (TH1F*)f1->Get("MVA_Scale_Down");
  TH1F * h_Scale_Up = (TH1F*)f1->Get("MVA_Scale_Up");

 double int_scaledown = h_Scale_Down->Integral();
 double int_scaleup = h_Scale_Up->Integral();

  double s = nExTT_Mu/h_Scale_Up->Integral();
  h_Scale_Up->Scale(s);
  s = nExTT_Mu/h_Scale_Down->Integral();
  h_Scale_Down->Scale(s);
  //  h_Scale_Up->Write("Scale_Up");
  // h_Scale_Down->Write("Scale_Down");

  /*
  TH1F *  h_Scale_Down_1stHalf = (TH1F*)f1->Get("MVA_Scale_Down_1stHalf");
  TH1F * h_Scale_Up_1stHalf = (TH1F*)f1->Get("MVA_Scale_Up_1stHalf");
  s = nExTT_Mu/h_Scale_Up_1stHalf->Integral();
  h_Scale_Up_1stHalf->Scale(s);
  s = nExTT_Mu/h_Scale_Down_1stHalf->Integral();
  h_Scale_Down_1stHalf->Scale(s);
  h_Scale_Up_1stHalf->Write("Scale_Up_1stHalf");
  h_Scale_Down_1stHalf->Write("Scale_Down_1stHalf");

 

  TH1F *  h_Scale_Down_2ndHalf = (TH1F*)f1->Get("MVA_Scale_Down_2ndHalf");
  TH1F * h_Scale_Up_2ndHalf = (TH1F*)f1->Get("MVA_Scale_Up_2ndHalf");
  s = nExTT_Mu/h_Scale_Up_2ndHalf->Integral();
  h_Scale_Up_2ndHalf->Scale(s);
  s = nExTT_Mu/h_Scale_Down_2ndHalf->Integral();
  h_Scale_Down_2ndHalf->Scale(s);
  h_Scale_Up_2ndHalf->Write("Scale_Up_2ndHalf");
  h_Scale_Down_2ndHalf->Write("Scale_Down_2ndHalf");
  */

  TH1F * h_Matching_Down = (TH1F*)f1->Get("MVA_Matching_Down");
  TH1F *h_Matching_Up = (TH1F*)f1->Get("MVA_Matching_Up");

  double int_matchingdown = h_Matching_Down->Integral();
  double int_matchingup = h_Matching_Up->Integral();
  s = nExTT_Mu/h_Matching_Up->Integral();
  h_Matching_Up->Scale(s);
  s = nExTT_Mu/h_Matching_Down->Integral();
  h_Matching_Down->Scale(s);
  //  h_Matching_Up->Write("Matching_Up");
  // h_Matching_Down->Write("Matching_Down");

  TH1F *h_JES_Down = (TH1F*)f1->Get("MVA_JES_Down");
  TH1F *h_JES_Up = (TH1F*)f1->Get("MVA_JES_Up");
 double int_jesdown = h_JES_Down->Integral();
 double int_jesup = h_JES_Up->Integral();
  s = nExTT_Mu/h_JES_Up->Integral();
  h_JES_Up->Scale(s);
  s = nExTT_Mu/h_JES_Down->Integral();
  h_JES_Down->Scale(s);
  //  h_JES_Up->Write("JES_Up");
  // h_JES_Down->Write("JES_Down");


  TH1F * h_ttbb_Down = (TH1F*)f1->Get("MVA_ttbb_Down");
   TH1F *h_ttbb_Up = (TH1F*)f1->Get("MVA_ttbb_Up");

 double int_ttbbdown = h_ttbb_Down->Integral();
 double int_ttbbup = h_ttbb_Up->Integral();

  s = nExTT_Mu/h_ttbb_Up->Integral();
  h_ttbb_Up->Scale(s);
  s = nExTT_Mu/h_ttbb_Down->Integral();
  h_ttbb_Down->Scale(s);
  //  h_ttbb_Up->Write("ttbb_Up");
  // h_ttbb_Down->Write("ttbb_Down");


  TH1F * h_bTag_Down = (TH1F*)f1->Get("MVA_bTag_Down");
   TH1F * h_bTag_Up = (TH1F*)f1->Get("MVA_bTag_Up");
 double int_btagdown = h_bTag_Down->Integral();
 double int_btagup = h_bTag_Up->Integral();
  s = nExTT_Mu/h_bTag_Up->Integral();
  h_bTag_Up->Scale(s);
  s = nExTT_Mu/h_bTag_Down->Integral();
  h_bTag_Down->Scale(s);
  //  h_bTag_Up->Write("bTag_Up");
  // h_bTag_Down->Write("bTag_Down");
 cout <<"processing sys shapes..2."<<endl;




   TH1F * h_misTag_Down = (TH1F*)f1->Get("MVA_misTag_Down");
   TH1F * h_misTag_Up = (TH1F*)f1->Get("MVA_misTag_Up");
 double int_mistagdown = h_misTag_Down->Integral();
 double int_mistagup = h_misTag_Up->Integral();
  s = nExTT_Mu/h_misTag_Up->Integral();
  h_misTag_Up->Scale(s);
  s = nExTT_Mu/h_misTag_Down->Integral();
  h_misTag_Down->Scale(s);
  //  h_misTag_Up->Write("misTag_Up");
  // h_misTag_Down->Write("misTag_Down");


  TH1F * h_leptonSF_Down = (TH1F*)f1->Get("MVA_leptonSF_Down");
   TH1F * h_leptonSF_Up = (TH1F*)f1->Get("MVA_leptonSF_Up");
 double int_leptonsfdown = h_leptonSF_Down->Integral();
 double int_leptonsfup = h_leptonSF_Up->Integral();
  s = nExTT_Mu/h_leptonSF_Up->Integral();
  h_leptonSF_Up->Scale(s);
  s = nExTT_Mu/h_leptonSF_Down->Integral();
  h_leptonSF_Down->Scale(s);
  //  h_leptonSF_Up->Write("leptonSF_Up");
  // h_leptonSF_Down->Write("leptonSF_Down");



  TH1F *h_PU_Down = (TH1F*)f1->Get("MVA_PU_Down");
  TH1F *h_PU_Up = (TH1F*)f1->Get("MVA_PU_Up");
 double int_pudown = h_PU_Down->Integral();
 double int_puup = h_PU_Up->Integral();
  s = nExTT_Mu/h_PU_Up->Integral();
  h_PU_Up->Scale(s);
  s = nExTT_Mu/h_PU_Down->Integral();
  h_PU_Down->Scale(s);
  //  h_PU_Up->Write("PU_Up");
  //  h_PU_Down->Write("PU_Down");


  TH1F *h_JER_Down = (TH1F*)f1->Get("MVA_JER_Down");
  TH1F *h_JER_Up = (TH1F*)f1->Get("MVA_JER_Up");
  double int_jerdown = h_JER_Down->Integral();
  double int_jerup = h_JER_Up->Integral();
 s = nExTT_Mu/h_JER_Up->Integral();
  h_JER_Up->Scale(s);
  s = nExTT_Mu/h_JER_Down->Integral();
  h_JER_Down->Scale(s);
  //  h_JER_Up->Write("JER_Up");
  // h_JER_Down->Write("JER_Down");



 cout <<"processing sys shapes..3."<<endl;
  myfile <<"Scale  & "<< h_Scale_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_Scale_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF")<<" & "<< h_Scale_Down->Chi2Test(hTT_Mu,"WW")<<" & "<< h_Scale_Up->Chi2Test(hTT_Mu,"WW")<<endl;
  myfile <<"Matching &  "<< h_Matching_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_Matching_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF")    << " & "<< h_Matching_Down->Chi2Test(hTT_Mu,"WW") << " & "<< h_Matching_Up->Chi2Test(hTT_Mu,"WW")<<  endl;
 myfile <<"JES  & "<< h_JES_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_JES_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF") << " & "<< h_JES_Down->Chi2Test(hTT_Mu,"WW") << " & "<< h_JES_Up->Chi2Test(hTT_Mu,"WW")<<    endl;
 myfile <<"ttbb  & "<< h_ttbb_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_ttbb_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF") << " & "<< h_ttbb_Down->Chi2Test(hTT_Mu,"WW") << " & "<< h_ttbb_Up->Chi2Test(hTT_Mu,"WW")<<endl;
 myfile <<"btag &  "<< h_bTag_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_bTag_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF")    << " & "<< h_bTag_Down->Chi2Test(hTT_Mu,"WW") << " & "<< h_bTag_Up->Chi2Test(hTT_Mu,"WW")<<endl;
 myfile <<"mistag &  "<< h_misTag_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_misTag_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF") << " & "<< h_misTag_Down->Chi2Test(hTT_Mu,"WW") << " & "<< h_misTag_Up->Chi2Test(hTT_Mu,"WW")   <<endl;
 myfile <<"lepton SF &  "<< h_leptonSF_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_leptonSF_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF") << " & "<< h_leptonSF_Down->Chi2Test(hTT_Mu,"WW") << " & "<< h_leptonSF_Up->Chi2Test(hTT_Mu,"WW")   <<endl;
 myfile <<"PU   &  "<< h_PU_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_PU_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF") << " & "<< h_PU_Down->Chi2Test(hTT_Mu,"WW") << " & "<< h_PU_Up->Chi2Test(hTT_Mu,"WW")   <<endl;
 myfile <<"JER &  "<< h_JER_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_JER_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF") << " & "<< h_JER_Down->Chi2Test(hTT_Mu,"WW") << " & "<< h_JER_Up->Chi2Test(hTT_Mu,"WW")   <<endl;


 //////////////////////////
 //// S I X  J E T  B I N 
 //////////////////////////

  TH1F *  h_Scale_Down_6j = (TH1F*)f1->Get("MVA6Jets_Scale_Down");
  TH1F * h_Scale_Up_6j = (TH1F*)f1->Get("MVA6Jets_Scale_Up");
   s = nExTT_Mu_6j/h_Scale_Up_6j->Integral();
  h_Scale_Up_6j->Scale(s);
  s = nExTT_Mu_6j/h_Scale_Down_6j->Integral();
  h_Scale_Down_6j->Scale(s);
  //  h_Scale_Up_6j->Write("Scale_Up_6j");
  // h_Scale_Down_6j->Write("Scale_Down_6j");


  TH1F * h_Matching_Down_6j = (TH1F*)f1->Get("MVA6Jets_Matching_Down");
  TH1F *h_Matching_Up_6j = (TH1F*)f1->Get("MVA6Jets_Matching_Up");
  s = nExTT_Mu_6j/h_Matching_Up_6j->Integral();
  h_Matching_Up_6j->Scale(s);
  s = nExTT_Mu_6j/h_Matching_Down_6j->Integral();
  h_Matching_Down_6j->Scale(s);
  //  h_Matching_Up_6j->Write("Matching_Up_6j");
  // h_Matching_Down_6j->Write("Matching_Down_6j");


  TH1F *h_JES_Down_6j = (TH1F*)f1->Get("MVA6Jets_JES_Down");
  TH1F *h_JES_Up_6j = (TH1F*)f1->Get("MVA6Jets_JES_Up");
  s = nExTT_Mu_6j/h_JES_Up_6j->Integral();
  h_JES_Up_6j->Scale(s);
  s = nExTT_Mu_6j/h_JES_Down_6j->Integral();
  h_JES_Down_6j->Scale(s);
  //  h_JES_Up_6j->Write("JES_Up_6j");
  // h_JES_Down_6j->Write("JES_Down_6j");


  TH1F * h_ttbb_Down_6j = (TH1F*)f1->Get("MVA6Jets_ttbb_Down");
   TH1F *h_ttbb_Up_6j = (TH1F*)f1->Get("MVA6Jets_ttbb_Up");
  s = nExTT_Mu_6j/h_ttbb_Up_6j->Integral();
  h_ttbb_Up_6j->Scale(s);
  s = nExTT_Mu_6j/h_ttbb_Down_6j->Integral();
  h_ttbb_Down_6j->Scale(s);
  //  h_ttbb_Up_6j->Write("ttbb_Up_6j");
  // h_ttbb_Down_6j->Write("ttbb_Down_6j");


  TH1F * h_bTag_Down_6j = (TH1F*)f1->Get("MVA6Jets_bTag_Down");
   TH1F * h_bTag_Up_6j = (TH1F*)f1->Get("MVA6Jets_bTag_Up");
  s = nExTT_Mu_6j/h_bTag_Up_6j->Integral();
  h_bTag_Up_6j->Scale(s);
  s = nExTT_Mu_6j/h_bTag_Down_6j->Integral();
  h_bTag_Down_6j->Scale(s);
  //  h_bTag_Up_6j->Write("bTag_Up_6j");
  // h_bTag_Down_6j->Write("bTag_Down_6j");


   TH1F * h_misTag_Down_6j = (TH1F*)f1->Get("MVA6Jets_misTag_Down");
   TH1F * h_misTag_Up_6j = (TH1F*)f1->Get("MVA6Jets_misTag_Up");
  s = nExTT_Mu_6j/h_misTag_Up_6j->Integral();
  h_misTag_Up_6j->Scale(s);
  s = nExTT_Mu_6j/h_misTag_Down_6j->Integral();
  h_misTag_Down_6j->Scale(s);
  //  h_misTag_Up_6j->Write("misTag_Up_6j");
  // h_misTag_Down_6j->Write("misTag_Down_6j");


  TH1F * h_leptonSF_Down_6j = (TH1F*)f1->Get("MVA6Jets_leptonSF_Down");
   TH1F * h_leptonSF_Up_6j = (TH1F*)f1->Get("MVA6Jets_leptonSF_Up");
  s = nExTT_Mu_6j/h_leptonSF_Up_6j->Integral();
  h_leptonSF_Up_6j->Scale(s);
  s = nExTT_Mu_6j/h_leptonSF_Down_6j->Integral();
  h_leptonSF_Down_6j->Scale(s);
  //  h_leptonSF_Up_6j->Write("leptonSF_Up_6j");
  // h_leptonSF_Down_6j->Write("leptonSF_Down_6j");


  TH1F *h_PU_Down_6j = (TH1F*)f1->Get("MVA6Jets_PU_Down");
  TH1F *h_PU_Up_6j = (TH1F*)f1->Get("MVA6Jets_PU_Up");
  s = nExTT_Mu_6j/h_PU_Up_6j->Integral();
  h_PU_Up_6j->Scale(s);
  s = nExTT_Mu_6j/h_PU_Down_6j->Integral();
  h_PU_Down_6j->Scale(s);
  //  h_PU_Up_6j->Write("PU_Up_6j");
  // h_PU_Down_6j->Write("PU_Down_6j");



  TH1F *h_JER_Down_6j = (TH1F*)f1->Get("MVA6Jets_JER_Down");
  TH1F *h_JER_Up_6j = (TH1F*)f1->Get("MVA6Jets_JER_Up");
  s = nExTT_Mu_6j/h_JER_Up_6j->Integral();
  h_JER_Up_6j->Scale(s);
  s = nExTT_Mu_6j/h_JER_Down_6j->Integral();
  h_JER_Down_6j->Scale(s);
  //  h_JER_Up_6j->Write("JER_Up_6j");
  // h_JER_Down_6j->Write("JER_Down_6j");


 /////////////////////////
 // S E V E N   J E T  B I N
 /////////////////////////

  TH1F *  h_Scale_Down_7j = (TH1F*)f1->Get("MVA7Jets_Scale_Down");
  TH1F * h_Scale_Up_7j = (TH1F*)f1->Get("MVA7Jets_Scale_Up");
  s = nExTT_Mu_7j/h_Scale_Up_7j->Integral();
  h_Scale_Up_7j->Scale(s);
  s = nExTT_Mu_7j/h_Scale_Down_7j->Integral();
  h_Scale_Down_7j->Scale(s);
  //  h_Scale_Up_7j->Write("Scale_Up_7j");
  // h_Scale_Down_7j->Write("Scale_Down_7j");



  TH1F * h_Matching_Down_7j = (TH1F*)f1->Get("MVA7Jets_Matching_Down");
  TH1F *h_Matching_Up_7j = (TH1F*)f1->Get("MVA7Jets_Matching_Up");
  s = nExTT_Mu_7j/h_Matching_Up_7j->Integral();
  h_Matching_Up_7j->Scale(s);
  s = nExTT_Mu_7j/h_Matching_Down_7j->Integral();
  h_Matching_Down_7j->Scale(s);
  //  h_Matching_Up_7j->Write("Matching_Up_7j");
  // h_Matching_Down_7j->Write("Matching_Down_7j");


  TH1F *h_JES_Down_7j = (TH1F*)f1->Get("MVA7Jets_JES_Down");
  TH1F *h_JES_Up_7j = (TH1F*)f1->Get("MVA7Jets_JES_Up");
  s = nExTT_Mu_7j/h_JES_Up_7j->Integral();
  h_JES_Up_7j->Scale(s);
  s = nExTT_Mu_7j/h_JES_Down_7j->Integral();
  h_JES_Down_7j->Scale(s);
  //  h_JES_Up_7j->Write("JES_Up_7j");
  // h_JES_Down_7j->Write("JES_Down_7j");



  TH1F * h_ttbb_Down_7j = (TH1F*)f1->Get("MVA7Jets_ttbb_Down");
   TH1F *h_ttbb_Up_7j = (TH1F*)f1->Get("MVA7Jets_ttbb_Up");
  s = nExTT_Mu_7j/h_ttbb_Up_7j->Integral();
  h_ttbb_Up_7j->Scale(s);
  s = nExTT_Mu_7j/h_ttbb_Down_7j->Integral();
  h_ttbb_Down_7j->Scale(s);
  //  h_ttbb_Up_7j->Write("ttbb_Up_7j");
  // h_ttbb_Down_7j->Write("ttbb_Down_7j");


  TH1F * h_bTag_Down_7j = (TH1F*)f1->Get("MVA7Jets_bTag_Down");
   TH1F * h_bTag_Up_7j = (TH1F*)f1->Get("MVA7Jets_bTag_Up");
  s = nExTT_Mu_7j/h_bTag_Up_7j->Integral();
  h_bTag_Up_7j->Scale(s);
  s = nExTT_Mu_7j/h_bTag_Down_7j->Integral();
  h_bTag_Down_7j->Scale(s);
  //  h_bTag_Up_7j->Write("bTag_Up_7j");
  // h_bTag_Down_7j->Write("bTag_Down_7j");


   TH1F * h_misTag_Down_7j = (TH1F*)f1->Get("MVA7Jets_misTag_Down");
   TH1F * h_misTag_Up_7j = (TH1F*)f1->Get("MVA7Jets_misTag_Up");
  s = nExTT_Mu_7j/h_misTag_Up_7j->Integral();
  h_misTag_Up_7j->Scale(s);
  s = nExTT_Mu_7j/h_misTag_Down_7j->Integral();
  h_misTag_Down_7j->Scale(s);
  //  h_misTag_Up_7j->Write("misTag_Up_7j");
  // h_misTag_Down_7j->Write("misTag_Down_7j");


  TH1F * h_leptonSF_Down_7j = (TH1F*)f1->Get("MVA7Jets_leptonSF_Down");
   TH1F * h_leptonSF_Up_7j = (TH1F*)f1->Get("MVA7Jets_leptonSF_Up");
  s = nExTT_Mu_7j/h_leptonSF_Up_7j->Integral();
  h_leptonSF_Up_7j->Scale(s);
  s = nExTT_Mu_7j/h_leptonSF_Down_7j->Integral();
  h_leptonSF_Down_7j->Scale(s);
  //  h_leptonSF_Up_7j->Write("leptonSF_Up_7j");
  // h_leptonSF_Down_7j->Write("leptonSF_Down_7j");


  TH1F *h_PU_Down_7j = (TH1F*)f1->Get("MVA7Jets_PU_Down");
  TH1F *h_PU_Up_7j = (TH1F*)f1->Get("MVA7Jets_PU_Up");
  s = nExTT_Mu_7j/h_PU_Up_7j->Integral();
  h_PU_Up_7j->Scale(s);
  s = nExTT_Mu_7j/h_PU_Down_7j->Integral();
  h_PU_Down_7j->Scale(s);
  //  h_PU_Up_7j->Write("PU_Up_7j");
  // h_PU_Down_7j->Write("PU_Down_7j");



  TH1F *h_JER_Down_7j = (TH1F*)f1->Get("MVA7Jets_JER_Down");
  TH1F *h_JER_Up_7j = (TH1F*)f1->Get("MVA7Jets_JER_Up");
  s = nExTT_Mu_7j/h_JER_Up_7j->Integral();
  h_JER_Up_7j->Scale(s);
  s = nExTT_Mu_7j/h_JER_Down_7j->Integral();
  h_JER_Down_7j->Scale(s);
  //  h_JER_Up_7j->Write("JER_Up_7j");
  // h_JER_Down_7j->Write("JER_Down_7j");

  /////////////////////////
 // E I G H T   J E T  B I N
 /////////////////////////

  TH1F *  h_Scale_Down_8j = (TH1F*)f1->Get("MVA8Jets_Scale_Down");
  TH1F * h_Scale_Up_8j = (TH1F*)f1->Get("MVA8Jets_Scale_Up");
  s = nExTT_Mu_8j/h_Scale_Up_8j->Integral();
  h_Scale_Up_8j->Scale(s);
  s = nExTT_Mu_8j/h_Scale_Down_8j->Integral();
  h_Scale_Down_8j->Scale(s);
  //  h_Scale_Up_8j->Write("Scale_Up_8j");
  // h_Scale_Down_8j->Write("Scale_Down_8j");



  TH1F * h_Matching_Down_8j = (TH1F*)f1->Get("MVA8Jets_Matching_Down");
  TH1F *h_Matching_Up_8j = (TH1F*)f1->Get("MVA8Jets_Matching_Up");
  s = nExTT_Mu_8j/h_Matching_Up_8j->Integral();
  h_Matching_Up_8j->Scale(s);
  s = nExTT_Mu_8j/h_Matching_Down_8j->Integral();
  h_Matching_Down_8j->Scale(s);
  //  h_Matching_Up_8j->Write("Matching_Up_8j");
  // h_Matching_Down_8j->Write("Matching_Down_8j");


  TH1F *h_JES_Down_8j = (TH1F*)f1->Get("MVA8Jets_JES_Down");
  TH1F *h_JES_Up_8j = (TH1F*)f1->Get("MVA8Jets_JES_Up");
  s = nExTT_Mu_8j/h_JES_Up_8j->Integral();
  h_JES_Up_8j->Scale(s);
  s = nExTT_Mu_8j/h_JES_Down_8j->Integral();
  h_JES_Down_8j->Scale(s);
  //  h_JES_Up_8j->Write("JES_Up_8j");
  // h_JES_Down_8j->Write("JES_Down_8j");



  TH1F * h_ttbb_Down_8j = (TH1F*)f1->Get("MVA8Jets_ttbb_Down");
   TH1F *h_ttbb_Up_8j = (TH1F*)f1->Get("MVA8Jets_ttbb_Up");
  s = nExTT_Mu_8j/h_ttbb_Up_8j->Integral();
  h_ttbb_Up_8j->Scale(s);
  s = nExTT_Mu_8j/h_ttbb_Down_8j->Integral();
  h_ttbb_Down_8j->Scale(s);
  //  h_ttbb_Up_8j->Write("ttbb_Up_8j");
  // h_ttbb_Down_8j->Write("ttbb_Down_8j");


  TH1F * h_bTag_Down_8j = (TH1F*)f1->Get("MVA8Jets_bTag_Down");
   TH1F * h_bTag_Up_8j = (TH1F*)f1->Get("MVA8Jets_bTag_Up");
  s = nExTT_Mu_8j/h_bTag_Up_8j->Integral();
  h_bTag_Up_8j->Scale(s);
  s = nExTT_Mu_8j/h_bTag_Down_8j->Integral();
  h_bTag_Down_8j->Scale(s);
  //  h_bTag_Up_8j->Write("bTag_Up_8j");
  // h_bTag_Down_8j->Write("bTag_Down_8j");


   TH1F * h_misTag_Down_8j = (TH1F*)f1->Get("MVA8Jets_misTag_Down");
   TH1F * h_misTag_Up_8j = (TH1F*)f1->Get("MVA8Jets_misTag_Up");
  s = nExTT_Mu_8j/h_misTag_Up_8j->Integral();
  h_misTag_Up_8j->Scale(s);
  s = nExTT_Mu_8j/h_misTag_Down_8j->Integral();
  h_misTag_Down_8j->Scale(s);
  //  h_misTag_Up_8j->Write("misTag_Up_8j");
  // h_misTag_Down_8j->Write("misTag_Down_8j");


  TH1F * h_leptonSF_Down_8j = (TH1F*)f1->Get("MVA8Jets_leptonSF_Down");
   TH1F * h_leptonSF_Up_8j = (TH1F*)f1->Get("MVA8Jets_leptonSF_Up");
  s = nExTT_Mu_8j/h_leptonSF_Up_8j->Integral();
  h_leptonSF_Up_8j->Scale(s);
  s = nExTT_Mu_8j/h_leptonSF_Down_8j->Integral();
  h_leptonSF_Down_8j->Scale(s);
  //  h_leptonSF_Up_8j->Write("leptonSF_Up_8j");
  // h_leptonSF_Down_8j->Write("leptonSF_Down_8j");


  TH1F *h_PU_Down_8j = (TH1F*)f1->Get("MVA8Jets_PU_Down");
  TH1F *h_PU_Up_8j = (TH1F*)f1->Get("MVA8Jets_PU_Up");
  s = nExTT_Mu_8j/h_PU_Up_8j->Integral();
  h_PU_Up_8j->Scale(s);
  s = nExTT_Mu_8j/h_PU_Down_8j->Integral();
  h_PU_Down_8j->Scale(s);
  //  h_PU_Up_8j->Write("PU_Up_8j");
  // h_PU_Down_8j->Write("PU_Down_8j");



  TH1F *h_JER_Down_8j = (TH1F*)f1->Get("MVA8Jets_JER_Down");
  TH1F *h_JER_Up_8j = (TH1F*)f1->Get("MVA8Jets_JER_Up");
  s = nExTT_Mu_8j/h_JER_Up_8j->Integral();
  h_JER_Up_8j->Scale(s);
  s = nExTT_Mu_8j/h_JER_Down_8j->Integral();
  h_JER_Down_8j->Scale(s);
  //  h_JER_Up_8j->Write("JER_Up_8j");
  // h_JER_Down_8j->Write("JER_Down_8j");
 


  double k_factor = 1.022727;

cout<<"Comparing effects of shape systematics..test. "<<endl;
//if (h_Scale_Up)cout<<"Comparing effects of shape systematics... "<<h_Scale_Up->Integral()<<endl;
 cout <<"Scale ="<<nExTT_Mu<<" "<< int_scaleup <<" "<<  int_scaledown <<" ==>    "<<100*(int_scaledown - nExTT_Mu )/nExTT_Mu<<"     "<<100*(int_scaleup - nExTT_Mu )/nExTT_Mu   <<  endl;
cout <<"Matching =  "<<nExTT_Mu << " "<< int_matchingup <<" "<< int_matchingdown <<" ==>    "<<100*(int_matchingdown - nExTT_Mu )/nExTT_Mu<<"     "<<100*(int_matchingup - nExTT_Mu )/nExTT_Mu     <<endl;
 cout <<"JES =  "<<nExTT_Mu<< " "<<k_factor*int_jesup <<" "<<k_factor*int_jesdown<<" ==>      "<<100*((k_factor*int_jesdown) - nExTT_Mu )/nExTT_Mu<<"     "<<100*((k_factor*int_jesup) - nExTT_Mu )/nExTT_Mu     <<endl;
 cout <<"JER =  "<<nExTT_Mu<< "   "<< k_factor*int_jerup <<"   "<< k_factor*int_jerdown <<" ==>      "<<100*((k_factor*int_jerdown) - nExTT_Mu )/nExTT_Mu<<"     "<<100*((k_factor*int_jerup) - nExTT_Mu )/nExTT_Mu     <<endl;
 cout <<"ttbb =  "<<  nExTT_Mu         << "   "<<k_factor*int_ttbbup <<"   "<< k_factor*int_ttbbdown<<" ==>      "<<100*((k_factor*int_ttbbdown) - nExTT_Mu )/nExTT_Mu<<"      "<<100*((k_factor*int_ttbbup) - nExTT_Mu )/nExTT_Mu     <<endl;
 cout <<"PU =  "<<nExTT_Mu << " "<<k_factor*int_puup <<"   "<<k_factor*int_pudown<<" ==>  "<<100*((k_factor*int_pudown) - nExTT_Mu )/nExTT_Mu<<"      "<<100*((k_factor*int_puup) - nExTT_Mu )/nExTT_Mu     <<endl;
 cout <<"lep sf = "<<nExTT_Mu << "   "<< k_factor*int_leptonsfup <<"   "<<k_factor*int_leptonsfdown<<" ==>      "<<100*((k_factor*int_leptonsfdown) - nExTT_Mu )/nExTT_Mu<<"     "<<100*((k_factor*int_leptonsfup) - nExTT_Mu )/nExTT_Mu     <<endl;
 cout <<"btag =  "<<nExTT_Mu << "   "<< k_factor*int_btagup <<"   "<< k_factor*int_btagdown <<" ==>      "<<100*((k_factor*int_btagdown) - nExTT_Mu )/nExTT_Mu<<"     "<<100*((k_factor*int_btagup) - nExTT_Mu )/nExTT_Mu     <<endl;
 cout <<"mistag =  "<<  nExTT_Mu<< "   "<< k_factor*int_mistagup <<"   "<< k_factor*int_mistagdown <<" ==>      "<<100*((k_factor*int_mistagdown) - nExTT_Mu )/nExTT_Mu<<"     "<<100*((k_factor*int_mistagup) - nExTT_Mu )/nExTT_Mu     <<endl;








 myfile.close();
  f2->Close();
  f1->Close();



}
