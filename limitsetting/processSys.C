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

 TFile* f_El = new TFile("FourTopnoRho_EventSelection_El.root"); 
 TFile* f1 = new TFile("SystematicShapes_El.root"); 


////////////////////////////////////////////////
// P R O C E S S   N O M I N A L   S H A P E S
////////////////////////////////////////////////

// TFile* fnom = new TFile("NominalShapes_El.root","RECREATE"); 


 //////////////////////
 /// I N C L U S I V E 
 //////////////////////

  TH1F * hData_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_Data");
  TH1F * hSig_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_NP_overlay_TTTT");

  double sig_scale = 0.87*1.3; //correcting for (very slightly incorrect int.lumi of sig histo and correcting to NLO xsection (doesn't make any real difference))


  hSig_El->Scale(sig_scale);

  //signal injection test
  TH1F * hdata_100fbtttt = new TH1F();
  TH1F * hdata_150fbtttt = new TH1F();
  TH1F * hdata_200fbtttt = new TH1F();

  hdata_100fbtttt = (TH1F*)hData_El->Clone();
  hdata_150fbtttt = (TH1F*)hData_El->Clone();
  hdata_200fbtttt = (TH1F*)hData_El->Clone();

  hdata_100fbtttt->Add(hSig_El, 100/1.3);
  hdata_150fbtttt->Add(hSig_El, 150/1.3);
  hdata_200fbtttt->Add(hSig_El, 200/1.3);

  hdata_100fbtttt->SetMarkerColor(kBlue);
  hdata_150fbtttt->SetMarkerColor(kGreen);
  hdata_200fbtttt->SetMarkerColor(kRed);
  hSig_El->SetLineColor(kOrange);
  hSig_El->SetLineWidth(2); 
  hData_El->SetLineWidth(2);
  hData_El->SetTitle("");
  hData_El->SetMinimum(0.001);
  hData_El->SetMaximum(100000);

  TCanvas * cx = new TCanvas();
  hData_El->Draw("HIST");
  hSig_El->Draw("HISTSAME");
  hdata_100fbtttt->Draw("same");
  hdata_150fbtttt->Draw("same");
  hdata_200fbtttt->Draw("same");
  cx->SetLogy();
  cx->Update();

  TLegend* leg = new TLegend(0.6, 0.73, .97, .95);
  leg->SetFillColor(0);
  leg->AddEntry(hData_El, "Data", "l");
  leg->AddEntry(hdata_100fbtttt, "Data with signal injection (#sigma = 100 fb)", "p");
  leg->AddEntry(hdata_150fbtttt, "Data with signal injection (#sigma = 150 fb)", "p");
  leg->AddEntry(hdata_200fbtttt, "Data with signal injection (#sigma = 200 fb)", "p");
  leg->AddEntry(hSig_El, "Signal", "l");
  leg->Draw();
 
  cx->SaveAs("SignalInjection_El.pdf");


  TH1F * hTTll_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_TTJets_ll");
  TH1F * hTTcc_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_TTJets_cc");
  TH1F * hTTbb_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_TTJets_bb");
  TH1F * hTT_El = (TH1F*)hTTbb_El->Clone();
  hTT_El->Reset();
  hTT_El->Add(hTTll_El);
  hTT_El->Add(hTTcc_El);
  hTT_El->Add(hTTbb_El);

  double nExTT_El  = hTT_El->Integral();


  //  hTT_El->Write("tt");
  // hSig_El->Write("tttt");
  // hData_El->Write("data");
  /// hdata_100fbtttt->Write("data_w100fbtttt");
  //hdata_150fbtttt->Write("data_w150fbtttt");
  // hdata_200fbtttt->Write("data_w200fbtttt");

  //  cout << hTTll_El->Integral() <<endl;
  //cout << hTTcc_El->Integral() <<endl;
  //cout << hTTbb_El->Integral() <<endl;
  //  cout << nExTT_El<<endl;

  //EW histo
  TH1F * hW4_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_W_4Jets");
  TH1F * hZ4_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_Z_4Jets");
  TH1F * hZ4_lm_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_Z_4Jets_lm");
  TH1F * hSingleTop_t_T_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_SingleTop_t_T");
  TH1F * hSingleTop_t_TBar_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_SingleTop_t_TBar");
  TH1F * hSingleTop_s_T_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_SingleTop_s_T");
  TH1F * hSingleTop_s_TBar_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_SingleTop_s_TBar");
  TH1F * hSingleTop_tW_T_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_SingleTop_tW_T");
  TH1F * hSingleTop_tW_TBar_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_SingleTop_tW_TBar");
  TH1F * hWW_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_WW");
  TH1F * hWZ_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_WZ");
  TH1F * hZZ_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_ZZ");

  TH1F * hEW_El = (TH1F*)hW4_El->Clone();
  hEW_El->Reset();

  hEW_El->Add(hW4_El);
  hEW_El->Add(hZ4_El);
  hEW_El->Add(hZ4_lm_El);
  hEW_El->Add(hSingleTop_t_T_El);
  hEW_El->Add(hSingleTop_t_TBar_El);
  hEW_El->Add(hSingleTop_s_T_El);
  hEW_El->Add(hSingleTop_s_TBar_El);
  hEW_El->Add(hSingleTop_tW_T_El);
  hEW_El->Add(hSingleTop_tW_TBar_El);
  hEW_El->Add(hWW_El);
  hEW_El->Add(hWZ_El);
  hEW_El->Add(hZZ_El);

  //  hEW_El->Write("EW");

  //tt other histo
  TH1F * httW_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_ttW");
  TH1F * httZ_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_ttZ");
  TH1F * httH_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_ttH");
  TH1F * hTTOther_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_TTJets_Other");

  hTTOther_El->Add(httW_El);
  hTTOther_El->Add(httZ_El);
  hTTOther_El->Add(httH_El);

  //  hTTOther_El->Write("ttOther");





  
 /////////////////////////
 // S I X   J E T  B I N
 /////////////////////////

  TH1F * hData_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_Data");
  TH1F * hSig_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_NP_overlay_TTTT");

 
 sig_scale = 0.87*1.3; //correcting for (very slightly incorrect int. lumi of data histo and correcting to NLO xsection (doesn't make any real difference))
  hSig_El_6j->Scale(sig_scale);

  //signal injection test
  TH1F * hdata_100fbtttt_6j = new TH1F();
  TH1F * hdata_150fbtttt_6j = new TH1F();
  TH1F * hdata_200fbtttt_6j = new TH1F();

  hdata_100fbtttt_6j = (TH1F*)hData_El_6j->Clone();
  hdata_150fbtttt_6j = (TH1F*)hData_El_6j->Clone();
  hdata_200fbtttt_6j = (TH1F*)hData_El_6j->Clone();

  hdata_100fbtttt_6j->Add(hSig_El_6j, 100/1.3);
  hdata_150fbtttt_6j->Add(hSig_El_6j, 150/1.3);
  hdata_200fbtttt_6j->Add(hSig_El_6j, 200/1.3);

  hdata_100fbtttt_6j->SetMarkerColor(kBlue);
  hdata_150fbtttt_6j->SetMarkerColor(kGreen);
  hdata_200fbtttt_6j->SetMarkerColor(kRed);
  hSig_El_6j->SetLineColor(kOrange);
  hSig_El_6j->SetLineWidth(2); 
  hData_El_6j->SetLineWidth(2);
  hData_El_6j->SetTitle("");
  hData_El_6j->SetMinimum(0.001);
  hData_El_6j->SetMaximum(100000);

  TCanvas * cx_6j = new TCanvas();
  hData_El_6j->Draw("HIST");
  hSig_El_6j->Draw("HISTSAME");
  hdata_100fbtttt_6j->Draw("same");
  hdata_150fbtttt_6j->Draw("same");
  hdata_200fbtttt_6j->Draw("same");
  cx_6j->SetLogy();
  cx_6j->Update();

  TLegend* leg_6j = new TLegend(0.6, 0.73, .97, .95);
  leg_6j->SetFillColor(0);
  leg_6j->AddEntry(hData_El_6j, "Data", "l");
  leg_6j->AddEntry(hdata_100fbtttt_6j, "Data with signal injection (#sigma = 100 fb)", "p");
  leg_6j->AddEntry(hdata_150fbtttt_6j, "Data with signal injection (#sigma = 150 fb)", "p");
  leg_6j->AddEntry(hdata_200fbtttt_6j, "Data with signal injection (#sigma = 200 fb)", "p");
  leg_6j->AddEntry(hSig_El_6j, "Signal", "l");
  leg_6j->Draw();
 
  cx_6j->SaveAs("SignalInjection_El_6j.pdf");


  TH1F * hTTll_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_TTJets_ll");
  TH1F * hTTcc_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_TTJets_cc");
  TH1F * hTTbb_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_TTJets_bb");
  TH1F * hTT_El_6j = (TH1F*)hTTbb_El_6j->Clone();
  hTT_El_6j->Reset();
  hTT_El_6j->Add(hTTll_El_6j);
  hTT_El_6j->Add(hTTcc_El_6j);
  hTT_El_6j->Add(hTTbb_El_6j);

  double nExTT_El_6j  = hTT_El_6j->Integral() ;

  //  hTT_El_6j->Write("tt_6j");
  //hSig_El_6j->Write("tttt_6j");
  // hData_El_6j->Write("data_6j");
  //hdata_100fbtttt_6j->Write("data_w100fbtttt_6j");
  //hdata_150fbtttt_6j->Write("data_w150fbtttt_6j");
  //hdata_200fbtttt_6j->Write("data_w200fbtttt_6j");

  //  cout << hTTll_El->Integral() <<endl;
  //cout << hTTcc_El->Integral() <<endl;
  //cout << hTTbb_El->Integral() <<endl;
  //  cout << nExTT_El<<endl;

  //EW histo
  TH1F * hW4_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_W_4Jets");
  TH1F * hZ4_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_Z_4Jets");
  TH1F * hZ4_lm_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_Z_4Jets_lm");
  TH1F * hSingleTop_t_T_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_SingleTop_t_T");
  TH1F * hSingleTop_t_TBar_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_SingleTop_t_TBar");
  TH1F * hSingleTop_s_T_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_SingleTop_s_T");
  TH1F * hSingleTop_s_TBar_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_SingleTop_s_TBar");
  TH1F * hSingleTop_tW_T_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_SingleTop_tW_T");
  TH1F * hSingleTop_tW_TBar_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_SingleTop_tW_TBar");
  TH1F * hWW_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_WW");
  TH1F * hWZ_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_WZ");
  TH1F * hZZ_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_ZZ");

  TH1F * hEW_El_6j = (TH1F*)hW4_El_6j->Clone();
  hEW_El_6j->Reset();

  hEW_El_6j->Add(hW4_El_6j);
  hEW_El_6j->Add(hZ4_El_6j);
  hEW_El_6j->Add(hZ4_lm_El_6j);
  hEW_El_6j->Add(hSingleTop_t_T_El_6j);
  hEW_El_6j->Add(hSingleTop_t_TBar_El_6j);
  hEW_El_6j->Add(hSingleTop_s_T_El_6j);
  hEW_El_6j->Add(hSingleTop_s_TBar_El_6j);
  hEW_El_6j->Add(hSingleTop_tW_T_El_6j);
  hEW_El_6j->Add(hSingleTop_tW_TBar_El_6j);
  hEW_El_6j->Add(hWW_El_6j);
  hEW_El_6j->Add(hWZ_El_6j);
  hEW_El_6j->Add(hZZ_El_6j);

  //  hEW_El_6j->Write("EW_6j");

  //tt other histo
  TH1F * httW_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_ttW");
  TH1F * httZ_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_ttZ");
  TH1F * httH_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_ttH");
  TH1F * hTTOther_El_6j = (TH1F*)f_El->Get("MultiSamplePlot_MVA6Jets/MVA6Jets_TTJets_Other");

  hTTOther_El_6j->Add(httW_El_6j);
  hTTOther_El_6j->Add(httZ_El_6j);
  hTTOther_El_6j->Add(httH_El_6j);

  //  hTTOther_El_6j->Write("ttOther_6j");


 /////////////////////////
 // S E V E N   J E T  B I N
 /////////////////////////

  TH1F * hData_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_Data");
  TH1F * hSig_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_NP_overlay_TTTT");

  //  double sig_scale = 0.87*1.3; //correcting for (very slightly incorrect int. lumi of data histo and correcting to NLO xsection (doesn't make any real difference))
  hSig_El_7j->Scale(sig_scale);
 

  //signal injection test
  TH1F * hdata_100fbtttt_7j = new TH1F();
  TH1F * hdata_150fbtttt_7j = new TH1F();
  TH1F * hdata_200fbtttt_7j = new TH1F();

  hdata_100fbtttt_7j = (TH1F*)hData_El_7j->Clone();
  hdata_150fbtttt_7j = (TH1F*)hData_El_7j->Clone();
  hdata_200fbtttt_7j = (TH1F*)hData_El_7j->Clone();

  hdata_100fbtttt_7j->Add(hSig_El_7j, 100/1.3);
  hdata_150fbtttt_7j->Add(hSig_El_7j, 150/1.3);
  hdata_200fbtttt_7j->Add(hSig_El_7j, 200/1.3);

  hdata_100fbtttt_7j->SetMarkerColor(kBlue);
  hdata_150fbtttt_7j->SetMarkerColor(kGreen);
  hdata_200fbtttt_7j->SetMarkerColor(kRed);
  hSig_El_7j->SetLineColor(kOrange);
  hSig_El_7j->SetLineWidth(2); 
  hData_El_7j->SetLineWidth(2);
  hData_El_7j->SetTitle("");
  hData_El_7j->SetMinimum(0.001);
  hData_El_7j->SetMaximum(100000);

  TCanvas * cx_7j = new TCanvas();
  hData_El_7j->Draw("HIST");
  hSig_El_7j->Draw("HISTSAME");
  hdata_100fbtttt_7j->Draw("same");
  hdata_150fbtttt_7j->Draw("same");
  hdata_200fbtttt_7j->Draw("same");
  cx_7j->SetLogy();
  cx_7j->Update();

  TLegend* leg_7j = new TLegend(0.6, 0.73, .97, .95);
  leg_7j->SetFillColor(0);
  leg_7j->AddEntry(hData_El_7j, "Data", "l");
  leg_7j->AddEntry(hdata_100fbtttt_7j, "Data with signal injection (#sigma = 100 fb)", "p");
  leg_7j->AddEntry(hdata_150fbtttt_7j, "Data with signal injection (#sigma = 150 fb)", "p");
  leg_7j->AddEntry(hdata_200fbtttt_7j, "Data with signal injection (#sigma = 200 fb)", "p");
  leg_7j->AddEntry(hSig_El_7j, "Signal", "l");
  leg_7j->Draw();
 
  cx_7j->SaveAs("SignalInjection_El_7j.pdf");


  TH1F * hTTll_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_TTJets_ll");
  TH1F * hTTcc_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_TTJets_cc");
  TH1F * hTTbb_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_TTJets_bb");
  TH1F * hTT_El_7j = (TH1F*)hTTbb_El_7j->Clone();
  hTT_El_7j->Reset();
  hTT_El_7j->Add(hTTll_El_7j);
  hTT_El_7j->Add(hTTcc_El_7j);
  hTT_El_7j->Add(hTTbb_El_7j);

  double nExTT_El_7j  = hTT_El_7j->Integral() ;

  //  hTT_El_7j->Write("tt_7j");
  // hSig_El_7j->Write("tttt_7j");
  //hData_El_7j->Write("data_7j");
  //hdata_100fbtttt_7j->Write("data_w100fbtttt_7j");
  //hdata_150fbtttt_7j->Write("data_w150fbtttt_7j");
  //hdata_200fbtttt_7j->Write("data_w200fbtttt_7j");

  //  cout << hTTll_El->Integral() <<endl;
  //cout << hTTcc_El->Integral() <<endl;
  //cout << hTTbb_El->Integral() <<endl;
  //  cout << nExTT_El<<endl;

  //EW histo
  TH1F * hW4_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_W_4Jets");
  TH1F * hZ4_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_Z_4Jets");
  TH1F * hZ4_lm_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_Z_4Jets_lm");
  TH1F * hSingleTop_t_T_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_SingleTop_t_T");
  TH1F * hSingleTop_t_TBar_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_SingleTop_t_TBar");
  TH1F * hSingleTop_s_T_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_SingleTop_s_T");
  TH1F * hSingleTop_s_TBar_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_SingleTop_s_TBar");
  TH1F * hSingleTop_tW_T_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_SingleTop_tW_T");
  TH1F * hSingleTop_tW_TBar_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_SingleTop_tW_TBar");
  TH1F * hWW_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_WW");
  TH1F * hWZ_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_WZ");
  TH1F * hZZ_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_ZZ");

  TH1F * hEW_El_7j = (TH1F*)hW4_El_7j->Clone();
  hEW_El_7j->Reset();

  hEW_El_7j->Add(hW4_El_7j);
  hEW_El_7j->Add(hZ4_El_7j);
  hEW_El_7j->Add(hZ4_lm_El_7j);
  hEW_El_7j->Add(hSingleTop_t_T_El_7j);
  hEW_El_7j->Add(hSingleTop_t_TBar_El_7j);
  hEW_El_7j->Add(hSingleTop_s_T_El_7j);
  hEW_El_7j->Add(hSingleTop_s_TBar_El_7j);
  hEW_El_7j->Add(hSingleTop_tW_T_El_7j);
  hEW_El_7j->Add(hSingleTop_tW_TBar_El_7j);
  hEW_El_7j->Add(hWW_El_7j);
  hEW_El_7j->Add(hWZ_El_7j);
  hEW_El_7j->Add(hZZ_El_7j);

  //  hEW_El_7j->Write("EW_7j");

  //tt other histo
  TH1F * httW_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_ttW");
  TH1F * httZ_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_ttZ");
  TH1F * httH_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_ttH");
  TH1F * hTTOther_El_7j = (TH1F*)f_El->Get("MultiSamplePlot_MVA7Jets/MVA7Jets_TTJets_Other");

  hTTOther_El_7j->Add(httW_El_7j);
  hTTOther_El_7j->Add(httZ_El_7j);
  hTTOther_El_7j->Add(httH_El_7j);

  //  hTTOther_El_7j->Write("ttOther_7j");



 /////////////////////////////
 // E I G H T   J E T  B I N
 /////////////////////////////

  TH1F * hData_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_Data");
  TH1F * hSig_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_NP_overlay_TTTT");

  // double sig_scale = 0.87*1.3; //correcting for (very slightly incorrect int. lumi of data histo and correcting to NLO xsection (doesn't make any real difference))
  hSig_El_8j->Scale(sig_scale);
 
  //signal injection test
  TH1F * hdata_100fbtttt_8j = new TH1F();
  TH1F * hdata_150fbtttt_8j = new TH1F();
  TH1F * hdata_200fbtttt_8j = new TH1F();

  hdata_100fbtttt_8j = (TH1F*)hData_El_8j->Clone();
  hdata_150fbtttt_8j = (TH1F*)hData_El_8j->Clone();
  hdata_200fbtttt_8j = (TH1F*)hData_El_8j->Clone();

  hdata_100fbtttt_8j->Add(hSig_El_8j, 100/1.3);
  hdata_150fbtttt_8j->Add(hSig_El_8j, 150/1.3);
  hdata_200fbtttt_8j->Add(hSig_El_8j, 200/1.3);

  hdata_100fbtttt_8j->SetMarkerColor(kBlue);
  hdata_150fbtttt_8j->SetMarkerColor(kGreen);
  hdata_200fbtttt_8j->SetMarkerColor(kRed);
  hSig_El_8j->SetLineColor(kOrange);
  hSig_El_8j->SetLineWidth(2); 
  hData_El_8j->SetLineWidth(2);
  hData_El_8j->SetTitle("");
  hData_El_8j->SetMinimum(0.001);
  hData_El_8j->SetMaximum(100000);

  TCanvas * cx_8j = new TCanvas();
  hData_El_8j->Draw("HIST");
  hSig_El_8j->Draw("HISTSAME");
  hdata_100fbtttt_8j->Draw("same");
  hdata_150fbtttt_8j->Draw("same");
  hdata_200fbtttt_8j->Draw("same");
  cx_8j->SetLogy();
  cx_8j->Update();

  TLegend* leg_8j = new TLegend(0.6, 0.73, .97, .95);
  leg_8j->SetFillColor(0);
  leg_8j->AddEntry(hData_El_8j, "Data", "l");
  leg_8j->AddEntry(hdata_100fbtttt_8j, "Data with signal injection (#sigma = 100 fb)", "p");
  leg_8j->AddEntry(hdata_150fbtttt_8j, "Data with signal injection (#sigma = 150 fb)", "p");
  leg_8j->AddEntry(hdata_200fbtttt_8j, "Data with signal injection (#sigma = 200 fb)", "p");
  leg_8j->AddEntry(hSig_El_8j, "Signal", "l");
  leg_8j->Draw();
 
  cx_8j->SaveAs("SignalInjection_El_8j.pdf");


  TH1F * hTTll_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_TTJets_ll");
  TH1F * hTTcc_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_TTJets_cc");
  TH1F * hTTbb_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_TTJets_bb");
  TH1F * hTT_El_8j = (TH1F*)hTTbb_El_8j->Clone();
  hTT_El_8j->Reset();
  hTT_El_8j->Add(hTTll_El_8j);
  hTT_El_8j->Add(hTTcc_El_8j);
  hTT_El_8j->Add(hTTbb_El_8j);

  double nExTT_El_8j  = hTT_El_8j->Integral() ;

  //  hTT_El_8j->Write("tt_8j");
  // hSig_El_8j->Write("tttt_8j");
  // hData_El_8j->Write("data_8j");
  //hdata_100fbtttt_8j->Write("data_w100fbtttt_8j");
  //hdata_150fbtttt_8j->Write("data_w150fbtttt_8j");
  //hdata_200fbtttt_8j->Write("data_w200fbtttt_8j");

  //  cout << hTTll_El->Integral() <<endl;
  //cout << hTTcc_El->Integral() <<endl;
  //cout << hTTbb_El->Integral() <<endl;
  //  cout << nExTT_El<<endl;

  //EW histo
  TH1F * hW4_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_W_4Jets");
  TH1F * hZ4_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_Z_4Jets");
  TH1F * hZ4_lm_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_Z_4Jets_lm");
  TH1F * hSingleTop_t_T_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_SingleTop_t_T");
  TH1F * hSingleTop_t_TBar_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_SingleTop_t_TBar");
  TH1F * hSingleTop_s_T_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_SingleTop_s_T");
  TH1F * hSingleTop_s_TBar_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_SingleTop_s_TBar");
  TH1F * hSingleTop_tW_T_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_SingleTop_tW_T");
  TH1F * hSingleTop_tW_TBar_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_SingleTop_tW_TBar");
  TH1F * hWW_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_WW");
  TH1F * hWZ_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_WZ");
  TH1F * hZZ_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_ZZ");

  TH1F * hEW_El_8j = (TH1F*)hW4_El_8j->Clone();
  hEW_El_8j->Reset();

  hEW_El_8j->Add(hW4_El_8j);
  hEW_El_8j->Add(hZ4_El_8j);
  hEW_El_8j->Add(hZ4_lm_El_8j);
  hEW_El_8j->Add(hSingleTop_t_T_El_8j);
  hEW_El_8j->Add(hSingleTop_t_TBar_El_8j);
  hEW_El_8j->Add(hSingleTop_s_T_El_8j);
  hEW_El_8j->Add(hSingleTop_s_TBar_El_8j);
  hEW_El_8j->Add(hSingleTop_tW_T_El_8j);
  hEW_El_8j->Add(hSingleTop_tW_TBar_El_8j);
  hEW_El_8j->Add(hWW_El_8j);
  hEW_El_8j->Add(hWZ_El_8j);
  hEW_El_8j->Add(hZZ_El_8j);

  //  hEW_El_8j->Write("EW_8j");

  //tt other histo
  TH1F * httW_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_ttW");
  TH1F * httZ_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_ttZ");
  TH1F * httH_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_ttH");
  TH1F * hTTOther_El_8j = (TH1F*)f_El->Get("MultiSamplePlot_MVA8Jets/MVA8Jets_TTJets_Other");

  hTTOther_El_8j->Add(httW_El_8j);
  hTTOther_El_8j->Add(httZ_El_8j);
  hTTOther_El_8j->Add(httH_El_8j);

  //  hTTOther_El_8j->Write("ttOther_8j");

  

  //  fnom->Close();
  cout<<"processing sytematic shapes..."<<endl;

///////////////////////////////////////////////////////
//// P R O C E S S   S Y S T E M A T I C   H I S T O S
///////////////////////////////////////////////////////

// TFile* f2 = new TFile("SystematicShapes_norm_El.root", "RECREATE"); 

 ///////////////////////
 /// I N C L U S I V E 
 ///////////////////////

 TH1F *h_Scale_Down = (TH1F*)f1->Get("MVA_Scale_Down");
 TH1F *h_Scale_Up = (TH1F*)f1->Get("MVA_Scale_Up");
 double int_scaledown = h_Scale_Down->Integral();
 double int_scaleup = h_Scale_Up->Integral();

ofstream myfile;
  myfile.open ("chi2_systematics_el.txt");
  myfile << "Chi2/NDOF of systematics:"<< endl;
 


 double s = nExTT_El/h_Scale_Up->Integral();
 h_Scale_Up->Scale(s);
 s = nExTT_El/h_Scale_Down->Integral();
 h_Scale_Down->Scale(s);
 //h_Scale_Up->Write("Scale_Up");
 //h_Scale_Down->Write("Scale_Down");


 /*
 TH1F *h_Scale_Down_1stHalf = (TH1F*)f1->Get("MVA_Scale_Down_1stHalf");
 TH1F *h_Scale_Up_1stHalf = (TH1F*)f1->Get("MVA_Scale_Up_1stHalf");
 double int_scaledown_1stHalf = h_Scale_Down_1stHalf->Integral();
 double int_scaleup_1stHalf = h_Scale_Up_1stHalf->Integral();

 s = nExTT_El/h_Scale_Up_1stHalf->Integral();
 h_Scale_Up_1stHalf->Scale(s);
 s = nExTT_El/h_Scale_Down_1stHalf->Integral();
 h_Scale_Down_1stHalf->Scale(s);
 h_Scale_Up_1stHalf->Write("Scale_Up_1stHalf");
 h_Scale_Down_1stHalf->Write("Scale_Down_1stHalf");


 TH1F *h_Scale_Down_2ndHalf = (TH1F*)f1->Get("MVA_Scale_Down_2ndHalf");
 TH1F *h_Scale_Up_2ndHalf = (TH1F*)f1->Get("MVA_Scale_Up_2ndHalf");
 double int_scaledown_2ndHalf = h_Scale_Down_2ndHalf->Integral();
 double int_scaleup_2ndHalf = h_Scale_Up_2ndHalf->Integral();

 s = nExTT_El/h_Scale_Up_2ndHalf->Integral();
 h_Scale_Up_2ndHalf->Scale(s);
 s = nExTT_El/h_Scale_Down_2ndHalf->Integral();
 h_Scale_Down_2ndHalf->Scale(s);
 h_Scale_Up_2ndHalf->Write("Scale_Up_2ndHalf");
 h_Scale_Down_2ndHalf->Write("Scale_Down_2ndHalf");


  cout<<"scale done....."<<endl;
 */



 //  cout<<"scale done....."<<endl;

 TH1F *h_Matching_Down = (TH1F*)f1->Get("MVA_Matching_Down");
 TH1F *h_Matching_Up = (TH1F*)f1->Get("MVA_Matching_Up");
 double int_matchingdown = h_Matching_Down->Integral();
 double int_matchingup = h_Matching_Up->Integral();



 s = nExTT_El/h_Matching_Up->Integral();
 h_Matching_Up->Scale(s);
 s = nExTT_El/h_Matching_Down->Integral();
 h_Matching_Down->Scale(s);
 //h_Matching_Up->Write("Matching_Up");
 //h_Matching_Down->Write("Matching_Down");

 //  myfile<<"matching done....."<<endl;

  TH1F *h_JES_Down = (TH1F*)f1->Get("MVA_JES_Down");
  TH1F *h_JES_Up = (TH1F*)f1->Get("MVA_JES_Up");
  double int_jesdown = h_JES_Down->Integral();
  double int_jesup = h_JES_Up->Integral();

  s = nExTT_El/h_JES_Up->Integral();
  h_JES_Up->Scale(s);
  s = nExTT_El/h_JES_Down->Integral();
  h_JES_Down->Scale(s);
  //h_JES_Up->Write("JES_Up");
  //h_JES_Down->Write("JES_Down");

  //  myfile<<"JES done....."<<endl;
 


  TH1F *h_ttbb_Down = (TH1F*)f1->Get("MVA_ttbb_Down");
  TH1F *h_ttbb_Up = (TH1F*)f1->Get("MVA_ttbb_Up");
  double int_ttbbdown = h_ttbb_Down->Integral();
  double int_ttbbup = h_ttbb_Up->Integral();
  s = nExTT_El/h_ttbb_Up->Integral();
  h_ttbb_Up->Scale(s);
  s = nExTT_El/h_ttbb_Down->Integral();
  h_ttbb_Down->Scale(s);




  //h_ttbb_Up->Write("ttbb_Up");
  //h_ttbb_Down->Write("ttbb_Down");

  //  myfile<<"ttbb done....."<<endl;

  TH1F *h_bTag_Down = (TH1F*)f1->Get("MVA_bTag_Down");
  TH1F *h_bTag_Up = (TH1F*)f1->Get("MVA_bTag_Up");
  double int_btagdown = h_bTag_Down->Integral();
  double int_btagup = h_bTag_Up->Integral();

  s = nExTT_El/h_bTag_Up->Integral();
  h_bTag_Up->Scale(s);
  s = nExTT_El/h_bTag_Down->Integral();
  h_bTag_Down->Scale(s);
  //h_bTag_Up->Write("bTag_Up");
  //h_bTag_Down->Write("bTag_Down");

  //  myfile<<"btag done....."<<endl;



  TH1F *h_misTag_Down = (TH1F*)f1->Get("MVA_misTag_Down");
  TH1F *h_misTag_Up = (TH1F*)f1->Get("MVA_misTag_Up");
 double int_mistagdown = h_misTag_Down->Integral();
  double int_mistagup = h_misTag_Up->Integral();

  s = nExTT_El/h_misTag_Up->Integral();
  h_misTag_Up->Scale(s);
  s = nExTT_El/h_misTag_Down->Integral();
  h_misTag_Down->Scale(s);
  // h_misTag_Up->Write("misTag_Up");
  //h_misTag_Down->Write("misTag_Down");

  //  myfile<<"mistag done....."<<endl;


  TH1F *h_leptonSF_Down = (TH1F*)f1->Get("MVA_LeptonSF_Down");
  TH1F *h_leptonSF_Up = (TH1F*)f1->Get("MVA_LeptonSF_Up");

  double int_leptonsfdown = h_leptonSF_Down->Integral();
  double int_leptonsfup = h_leptonSF_Up->Integral();

  s = nExTT_El/h_leptonSF_Up->Integral();
  h_leptonSF_Up->Scale(s);
  s = nExTT_El/h_leptonSF_Down->Integral();
  h_leptonSF_Down->Scale(s);
  //h_leptonSF_Up->Write("leptonSF_Up");
  //h_leptonSF_Down->Write("leptonSF_Down");

  

  //  myfile<<"lepton sf done....."<<endl;

  TH1F *h_PU_Down = (TH1F*)f1->Get("MVA_PU_Down");
  TH1F *h_PU_Up = (TH1F*)f1->Get("MVA_PU_Up");
  double int_pudown = h_PU_Down->Integral();
  double int_puup = h_PU_Up->Integral();

  s = nExTT_El/h_PU_Up->Integral();
  h_PU_Up->Scale(s);
  s = nExTT_El/h_PU_Down->Integral();
  h_PU_Down->Scale(s);
  //h_PU_Up->Write("PU_Up");
  //h_PU_Down->Write("PU_Down");
 


  //  myfile<<"PU  done....."<<endl;

  TH1F *h_JER_Down = (TH1F*)f1->Get("MVA_JER_Down");
  TH1F *h_JER_Up = (TH1F*)f1->Get("MVA_JER_Up");
  double int_jerdown = h_JER_Down->Integral();
  double int_jerup = h_JER_Up->Integral();
  s = nExTT_El/h_JER_Up->Integral();
  h_JER_Up->Scale(s);
  s = nExTT_El/h_JER_Down->Integral();
  h_JER_Down->Scale(s);
  //h_JER_Up->Write("JER_Up");
  // h_JER_Down->Write("JER_Down");


  //  myfile<<"JER done....."<<endl;
  myfile <<"Scale  & "<< h_Scale_Down->Chi2Test(hTT_El,"WWCHI2/NDF") <<" &  "<< h_Scale_Up->Chi2Test(hTT_El,"WWCHI2/NDF")<<" & "<< h_Scale_Down->Chi2Test(hTT_El,"WW")<<" & "<< h_Scale_Up->Chi2Test(hTT_El,"WW")<<endl;
  myfile <<"Matching &  "<< h_Matching_Down->Chi2Test(hTT_El,"WWCHI2/NDF") <<" &  "<< h_Matching_Up->Chi2Test(hTT_El,"WWCHI2/NDF")    << " & "<< h_Matching_Down->Chi2Test(hTT_El,"WW") << " & "<< h_Matching_Up->Chi2Test(hTT_El,"WW")<<  endl;
 myfile <<"JES  & "<< h_JES_Down->Chi2Test(hTT_El,"WWCHI2/NDF") <<" &  "<< h_JES_Up->Chi2Test(hTT_El,"WWCHI2/NDF") << " & "<< h_JES_Down->Chi2Test(hTT_El,"WW") << " & "<< h_JES_Up->Chi2Test(hTT_El,"WW")<<    endl;
 myfile <<"ttbb  & "<< h_ttbb_Down->Chi2Test(hTT_El,"WWCHI2/NDF") <<" &  "<< h_ttbb_Up->Chi2Test(hTT_El,"WWCHI2/NDF") << " & "<< h_ttbb_Down->Chi2Test(hTT_El,"WW") << " & "<< h_ttbb_Up->Chi2Test(hTT_El,"WW")<<endl;
 myfile <<"btag &  "<< h_bTag_Down->Chi2Test(hTT_El,"WWCHI2/NDF") <<" &  "<< h_bTag_Up->Chi2Test(hTT_El,"WWCHI2/NDF")    << " & "<< h_bTag_Down->Chi2Test(hTT_El,"WW") << " & "<< h_bTag_Up->Chi2Test(hTT_El,"WW")<<endl;
 myfile <<"mistag &  "<< h_misTag_Down->Chi2Test(hTT_El,"WWCHI2/NDF") <<" &  "<< h_misTag_Up->Chi2Test(hTT_El,"WWCHI2/NDF") << " & "<< h_misTag_Down->Chi2Test(hTT_El,"WW") << " & "<< h_misTag_Up->Chi2Test(hTT_El,"WW")   <<endl;
 myfile <<"lepton SF &  "<< h_leptonSF_Down->Chi2Test(hTT_El,"WWCHI2/NDF") <<" &  "<< h_leptonSF_Up->Chi2Test(hTT_El,"WWCHI2/NDF") << " & "<< h_leptonSF_Down->Chi2Test(hTT_El,"WW") << " & "<< h_leptonSF_Up->Chi2Test(hTT_El,"WW")   <<endl;
 myfile <<"PU   &  "<< h_PU_Down->Chi2Test(hTT_El,"WWCHI2/NDF") <<" &  "<< h_PU_Up->Chi2Test(hTT_El,"WWCHI2/NDF") << " & "<< h_PU_Down->Chi2Test(hTT_El,"WW") << " & "<< h_PU_Up->Chi2Test(hTT_El,"WW")   <<endl;
 myfile <<"JER &  "<< h_JER_Down->Chi2Test(hTT_El,"WWCHI2/NDF") <<" &  "<< h_JER_Up->Chi2Test(hTT_El,"WWCHI2/NDF") << " & "<< h_JER_Down->Chi2Test(hTT_El,"WW") << " & "<< h_JER_Up->Chi2Test(hTT_El,"WW")   <<endl;

 cout <<"Scale ="<<nExTT_El<<" "<<h_Scale_Up->Integral() <<" "<<  h_Scale_Down->Integral() <<endl;
 cout <<"Matching ="<<nExTT_El<<" "<<h_Matching_Up->Integral() <<" "<<  h_Matching_Down->Integral() <<endl;
 cout <<"JES ="<<nExTT_El<<" "<<h_JES_Up->Integral() <<" "<<  h_JES_Down->Integral() <<endl;
 cout <<"ttbb ="<<nExTT_El<<" "<<h_ttbb_Up->Integral() <<" "<<  h_ttbb_Down->Integral() <<endl;
 cout <<"btag ="<<nExTT_El<<" "<<h_bTag_Up->Integral() <<" "<<  h_bTag_Down->Integral() <<endl;
 cout <<"mistag ="<<nExTT_El<<" "<<h_misTag_Up->Integral() <<" "<<  h_misTag_Down->Integral() <<endl;
 cout <<"lepton SF ="<<nExTT_El<<" "<<h_leptonSF_Up->Integral() <<" "<<  h_leptonSF_Down->Integral() <<endl;
 cout <<"PU ="<<nExTT_El<<" "<<h_PU_Up->Integral() <<" "<<  h_PU_Down->Integral() <<endl; 
 cout <<"JER ="<<nExTT_El<<" "<<h_JER_Up->Integral() <<" "<<  h_JER_Down->Integral() <<endl; 
 


  /////////////////////////////
  ////  S I X  J E T  B I N 
  /////////////////////////////


 TH1F *h_Scale_Down_6j = (TH1F*)f1->Get("MVA6Jets_Scale_Down");
 TH1F *h_Scale_Up_6j = (TH1F*)f1->Get("MVA6Jets_Scale_Up");
  double int_6j_scaledown = h_Scale_Down_6j->Integral();
  double int_6j_scaleup = h_Scale_Up_6j->Integral();

 s = nExTT_El_6j/h_Scale_Up_6j->Integral();
 h_Scale_Up_6j->Scale(s);
 s = nExTT_El_6j/h_Scale_Down_6j->Integral();
 h_Scale_Down_6j->Scale(s);
 //h_Scale_Up_6j->Write("Scale_Up_6j");
 //h_Scale_Down_6j->Write("Scale_Down_6j");

 //  cout<<"scale done....."<<endl;

 TH1F *h_Matching_Down_6j = (TH1F*)f1->Get("MVA6Jets_Matching_Down");
 TH1F *h_Matching_Up_6j = (TH1F*)f1->Get("MVA6Jets_Matching_Up");
 double int_6j_matchingdown = h_Matching_Down_6j->Integral();
 double int_6j_matchingup = h_Matching_Up_6j->Integral();


 s = nExTT_El_6j/h_Matching_Up_6j->Integral();
 h_Matching_Up_6j->Scale(s);
 s = nExTT_El_6j/h_Matching_Down_6j->Integral();
 h_Matching_Down_6j->Scale(s);
 //h_Matching_Up_6j->Write("Matching_Up_6j");
 //h_Matching_Down_6j->Write("Matching_Down_6j");

 // cout<<"matching done....."<<endl;

  TH1F *h_JES_Down_6j = (TH1F*)f1->Get("MVA6Jets_JES_Down");
  TH1F *h_JES_Up_6j = (TH1F*)f1->Get("MVA6Jets_JES_Up");
  double int_6j_jesdown = h_JES_Down_6j->Integral();
  double int_6j_jesup = h_JES_Up_6j->Integral();

  s = nExTT_El_6j/h_JES_Up_6j->Integral();
  h_JES_Up_6j->Scale(s);
  s = nExTT_El_6j/h_JES_Down_6j->Integral();
  h_JES_Down_6j->Scale(s);
  //h_JES_Up_6j->Write("JES_Up_6j");
  //h_JES_Down_6j->Write("JES_Down_6j");

  //  cout<<"JES done....."<<endl;


  TH1F *h_ttbb_Down_6j = (TH1F*)f1->Get("MVA6Jets_ttbb_Down");
  TH1F *h_ttbb_Up_6j = (TH1F*)f1->Get("MVA6Jets_ttbb_Up");
  double int_6j_ttbbdown = h_ttbb_Down_6j->Integral();
  double int_6j_ttbbup = h_ttbb_Up_6j->Integral();
  s = nExTT_El_6j/h_ttbb_Up_6j->Integral();
  h_ttbb_Up_6j->Scale(s);
  s = nExTT_El_6j/h_ttbb_Down_6j->Integral();
  h_ttbb_Down_6j->Scale(s);
  //h_ttbb_Up_6j->Write("ttbb_Up_6j");
  //h_ttbb_Down_6j->Write("ttbb_Down_6j");

  //  cout<<"ttbb done....."<<endl;

  TH1F *h_bTag_Down_6j = (TH1F*)f1->Get("MVA6Jets_bTag_Down");
  TH1F *h_bTag_Up_6j = (TH1F*)f1->Get("MVA6Jets_bTag_Up");
  double int_6j_btagdown = h_bTag_Down_6j->Integral();
  double int_6j_btagup = h_bTag_Up_6j->Integral();

  s = nExTT_El_6j/h_bTag_Up_6j->Integral();
  h_bTag_Up_6j->Scale(s);
  s = nExTT_El_6j/h_bTag_Down_6j->Integral();
  h_bTag_Down_6j->Scale(s);
  //h_bTag_Up_6j->Write("bTag_Up_6j");
  //h_bTag_Down_6j->Write("bTag_Down_6j");

  //  cout<<"btag done....."<<endl;


  TH1F *h_misTag_Down_6j = (TH1F*)f1->Get("MVA6Jets_misTag_Down");
  TH1F *h_misTag_Up_6j = (TH1F*)f1->Get("MVA6Jets_misTag_Up");
 double int_6j_mistagdown = h_misTag_Down_6j->Integral();
  double int_6j_mistagup = h_misTag_Up_6j->Integral();

  s = nExTT_El_6j/h_misTag_Up_6j->Integral();
  h_misTag_Up_6j->Scale(s);
  s = nExTT_El_6j/h_misTag_Down_6j->Integral();
  h_misTag_Down_6j->Scale(s);
  //h_misTag_Up_6j->Write("misTag_Up_6j");
  //h_misTag_Down_6j->Write("misTag_Down_6j");

  //  cout<<"mistag done....."<<endl;

  TH1F *h_leptonSF_Down_6j = (TH1F*)f1->Get("MVA6Jets_LeptonSF_Down");
  TH1F *h_leptonSF_Up_6j = (TH1F*)f1->Get("MVA6Jets_LeptonSF_Up");

  double int_6j_leptonsfdown = h_leptonSF_Down_6j->Integral();
  double int_6j_leptonsfup = h_leptonSF_Up_6j->Integral();

  s = nExTT_El_6j/h_leptonSF_Up_6j->Integral();
  h_leptonSF_Up_6j->Scale(s);
  s = nExTT_El_6j/h_leptonSF_Down_6j->Integral();
  h_leptonSF_Down_6j->Scale(s);
  //h_leptonSF_Up_6j->Write("leptonSF_Up_6j");
  //h_leptonSF_Down_6j->Write("leptonSF_Down_6j");


  //  cout<<"lepton sf done....."<<endl;

  TH1F *h_PU_Down_6j = (TH1F*)f1->Get("MVA6Jets_PU_Down");
  TH1F *h_PU_Up_6j = (TH1F*)f1->Get("MVA6Jets_PU_Up");
  double int_6j_pudown = h_PU_Down_6j->Integral();
  double int_6j_puup = h_PU_Up_6j->Integral();

  s = nExTT_El_6j/h_PU_Up_6j->Integral();
  h_PU_Up_6j->Scale(s);
  s = nExTT_El_6j/h_PU_Down_6j->Integral();
  h_PU_Down_6j->Scale(s);
  //h_PU_Up_6j->Write("PU_Up_6j");
  //h_PU_Down_6j->Write("PU_Down_6j");


  //  cout<<"PU  done....."<<endl;

  TH1F *h_JER_Down_6j = (TH1F*)f1->Get("MVA6Jets_JER_Down");
  TH1F *h_JER_Up_6j = (TH1F*)f1->Get("MVA6Jets_JER_Up");
  double int_6j_jerdown = h_JER_Down_6j->Integral();
  double int_6j_jerup = h_JER_Up_6j->Integral();
  s = nExTT_El_6j/h_JER_Up_6j->Integral();
  h_JER_Up_6j->Scale(s);
  s = nExTT_El_6j/h_JER_Down_6j->Integral();
  h_JER_Down_6j->Scale(s);
  //h_JER_Up_6j->Write("JER_Up_6j");
  // h_JER_Down_6j->Write("JER_Down_6j");


  //  cout<<"JER done....."<<endl;
  //  cout <<"seven jet bin"<<endl;

 /////////////////////////
 // S E V E N    J E T  B I N
 /////////////////////////

 TH1F *h_Scale_Down_7j = (TH1F*)f1->Get("MVA7Jets_Scale_Down");
 TH1F *h_Scale_Up_7j = (TH1F*)f1->Get("MVA7Jets_Scale_Up");
 s = nExTT_El_7j/h_Scale_Up_7j->Integral();
 h_Scale_Up_7j->Scale(s);
 s = nExTT_El_7j/h_Scale_Down_7j->Integral();
 h_Scale_Down_7j->Scale(s);
 //h_Scale_Up_7j->Write("Scale_Up_7j");
 //h_Scale_Down_7j->Write("Scale_Down_7j");


 TH1F *h_Matching_Down_7j = (TH1F*)f1->Get("MVA7Jets_Matching_Down");
 TH1F *h_Matching_Up_7j = (TH1F*)f1->Get("MVA7Jets_Matching_Up");
 s = nExTT_El_7j/h_Matching_Up_7j->Integral();
 h_Matching_Up_7j->Scale(s);
 s = nExTT_El_7j/h_Matching_Down_7j->Integral();
 h_Matching_Down_7j->Scale(s);
 // h_Matching_Up_7j->Write("Matching_Up_7j");
 //h_Matching_Down_7j->Write("Matching_Down_7j");


  TH1F *h_JES_Down_7j = (TH1F*)f1->Get("MVA7Jets_JES_Down");
  TH1F *h_JES_Up_7j = (TH1F*)f1->Get("MVA7Jets_JES_Up");
  s = nExTT_El_7j/h_JES_Up_7j->Integral();
  h_JES_Up_7j->Scale(s);
  s = nExTT_El_7j/h_JES_Down_7j->Integral();
  h_JES_Down_7j->Scale(s);
  //h_JES_Up_7j->Write("JES_Up_7j");
  // h_JES_Down_7j->Write("JES_Down_7j");


  TH1F *h_ttbb_Down_7j = (TH1F*)f1->Get("MVA7Jets_ttbb_Down");
  TH1F *h_ttbb_Up_7j = (TH1F*)f1->Get("MVA7Jets_ttbb_Up");
  s = nExTT_El_7j/h_ttbb_Up_7j->Integral();
  h_ttbb_Up_7j->Scale(s);
  s = nExTT_El_7j/h_ttbb_Down_7j->Integral();
  h_ttbb_Down_7j->Scale(s);
  //h_ttbb_Up_7j->Write("ttbb_Up_7j");
  // h_ttbb_Down_7j->Write("ttbb_Down_7j");


  TH1F *h_bTag_Down_7j = (TH1F*)f1->Get("MVA7Jets_bTag_Down");
  TH1F *h_bTag_Up_7j = (TH1F*)f1->Get("MVA7Jets_bTag_Up");
  s = nExTT_El_7j/h_bTag_Up_7j->Integral();
  h_bTag_Up_7j->Scale(s);
  s = nExTT_El_7j/h_bTag_Down_7j->Integral();
  h_bTag_Down_7j->Scale(s);
  //  h_bTag_Up_7j->Write("bTag_Up_7j");
  // h_bTag_Down_7j->Write("bTag_Down_7j");


  TH1F *h_misTag_Down_7j = (TH1F*)f1->Get("MVA7Jets_misTag_Down");
  TH1F *h_misTag_Up_7j = (TH1F*)f1->Get("MVA7Jets_misTag_Up");
  s = nExTT_El_7j/h_misTag_Up_7j->Integral();
  h_misTag_Up_7j->Scale(s);
  s = nExTT_El_7j/h_misTag_Down_7j->Integral();
  h_misTag_Down_7j->Scale(s);
  // h_misTag_Up_7j->Write("misTag_Up_7j");
  //h_misTag_Down_7j->Write("misTag_Down_7j");


  TH1F *h_leptonSF_Down_7j = (TH1F*)f1->Get("MVA7Jets_LeptonSF_Down");
  TH1F *h_leptonSF_Up_7j = (TH1F*)f1->Get("MVA7Jets_LeptonSF_Up");
  s = nExTT_El_7j/h_leptonSF_Up_7j->Integral();
  h_leptonSF_Up_7j->Scale(s);
  s = nExTT_El_7j/h_leptonSF_Down_7j->Integral();
  h_leptonSF_Down_7j->Scale(s);
  // h_leptonSF_Up_7j->Write("leptonSF_Up_7j");
  // h_leptonSF_Down_7j->Write("leptonSF_Down_7j");


  TH1F *h_PU_Down_7j = (TH1F*)f1->Get("MVA7Jets_PU_Down");
  TH1F *h_PU_Up_7j = (TH1F*)f1->Get("MVA7Jets_PU_Up");
  s = nExTT_El_7j/h_PU_Up_7j->Integral();
  h_PU_Up_7j->Scale(s);
  s = nExTT_El_7j/h_PU_Down_7j->Integral();
  h_PU_Down_7j->Scale(s);
  // h_PU_Up_7j->Write("PU_Up_7j");
  //h_PU_Down_7j->Write("PU_Down_7j");


  TH1F *h_JER_Down_7j = (TH1F*)f1->Get("MVA7Jets_JER_Down");
  TH1F *h_JER_Up_7j = (TH1F*)f1->Get("MVA7Jets_JER_Up");
  s = nExTT_El_7j/h_JER_Up_7j->Integral();
  h_JER_Up_7j->Scale(s);
  s = nExTT_El_7j/h_JER_Down_7j->Integral();
  h_JER_Down_7j->Scale(s);
  //h_JER_Up_7j->Write("JER_Up_7j");
  // h_JER_Down_7j->Write("JER_Down_7j");

  cout <<"eight jet bin"<<endl;

 /////////////////////////
 // E I G H T  J E T  B I N
 /////////////////////////

 TH1F *h_Scale_Down_8j = (TH1F*)f1->Get("MVA8Jets_Scale_Down");
 TH1F *h_Scale_Up_8j = (TH1F*)f1->Get("MVA8Jets_Scale_Up");
 s = nExTT_El_8j/h_Scale_Up_8j->Integral();
 h_Scale_Up_8j->Scale(s);
 s = nExTT_El_8j/h_Scale_Down_8j->Integral();
 h_Scale_Down_8j->Scale(s);
 // h_Scale_Up_8j->Write("Scale_Up_8j");
 // h_Scale_Down_8j->Write("Scale_Down_8j");


 TH1F *h_Matching_Down_8j = (TH1F*)f1->Get("MVA8Jets_Matching_Down");
 TH1F *h_Matching_Up_8j = (TH1F*)f1->Get("MVA8Jets_Matching_Up");
 s = nExTT_El_8j/h_Matching_Up_8j->Integral();
 h_Matching_Up_8j->Scale(s);
 s = nExTT_El_8j/h_Matching_Down_8j->Integral();
 h_Matching_Down_8j->Scale(s);
 // h_Matching_Up_8j->Write("Matching_Up_8j");
 //h_Matching_Down_8j->Write("Matching_Down_8j");


  TH1F *h_JES_Down_8j = (TH1F*)f1->Get("MVA8Jets_JES_Down");
  TH1F *h_JES_Up_8j = (TH1F*)f1->Get("MVA8Jets_JES_Up");
  s = nExTT_El_8j/h_JES_Up_8j->Integral();
  h_JES_Up_8j->Scale(s);
  s = nExTT_El_8j/h_JES_Down_8j->Integral();
  h_JES_Down_8j->Scale(s);
  // h_JES_Up_8j->Write("JES_Up_8j");
  //h_JES_Down_8j->Write("JES_Down_8j");


  TH1F *h_ttbb_Down_8j = (TH1F*)f1->Get("MVA8Jets_ttbb_Down");
  TH1F *h_ttbb_Up_8j = (TH1F*)f1->Get("MVA8Jets_ttbb_Up");
  s = nExTT_El_8j/h_ttbb_Up_8j->Integral();
  h_ttbb_Up_8j->Scale(s);
  s = nExTT_El_8j/h_ttbb_Down_8j->Integral();
  h_ttbb_Down_8j->Scale(s);
  // h_ttbb_Up_8j->Write("ttbb_Up_8j");
  //h_ttbb_Down_8j->Write("ttbb_Down_8j");


  TH1F *h_bTag_Down_8j = (TH1F*)f1->Get("MVA8Jets_bTag_Down");
  TH1F *h_bTag_Up_8j = (TH1F*)f1->Get("MVA8Jets_bTag_Up");
  s = nExTT_El_8j/h_bTag_Up_8j->Integral();
  h_bTag_Up_8j->Scale(s);
  s = nExTT_El_8j/h_bTag_Down_8j->Integral();
  h_bTag_Down_8j->Scale(s);
  // h_bTag_Up_8j->Write("bTag_Up_8j");
  // h_bTag_Down_8j->Write("bTag_Down_8j");


  TH1F *h_misTag_Down_8j = (TH1F*)f1->Get("MVA8Jets_misTag_Down");
  TH1F *h_misTag_Up_8j = (TH1F*)f1->Get("MVA8Jets_misTag_Up");
  s = nExTT_El_8j/h_misTag_Up_8j->Integral();
  h_misTag_Up_8j->Scale(s);
  s = nExTT_El_8j/h_misTag_Down_8j->Integral();
  h_misTag_Down_8j->Scale(s);
  // h_misTag_Up_8j->Write("misTag_Up_8j");
  // h_misTag_Down_8j->Write("misTag_Down_8j");


  TH1F *h_leptonSF_Down_8j = (TH1F*)f1->Get("MVA8Jets_LeptonSF_Down");
  TH1F *h_leptonSF_Up_8j = (TH1F*)f1->Get("MVA8Jets_LeptonSF_Up");
  s = nExTT_El_8j/h_leptonSF_Up_8j->Integral();
  h_leptonSF_Up_8j->Scale(s);
  s = nExTT_El_8j/h_leptonSF_Down_8j->Integral();
  h_leptonSF_Down_8j->Scale(s);
  // h_leptonSF_Up_8j->Write("leptonSF_Up_8j");
  // h_leptonSF_Down_8j->Write("leptonSF_Down_8j");


  TH1F *h_PU_Down_8j = (TH1F*)f1->Get("MVA8Jets_PU_Down");
  TH1F *h_PU_Up_8j = (TH1F*)f1->Get("MVA8Jets_PU_Up");
  s = nExTT_El_8j/h_PU_Up_8j->Integral();
  h_PU_Up_8j->Scale(s);
  s = nExTT_El_8j/h_PU_Down_8j->Integral();
  h_PU_Down_8j->Scale(s);
  // h_PU_Up_8j->Write("PU_Up_8j");
  // h_PU_Down_8j->Write("PU_Down_8j");


  TH1F *h_JER_Down_8j = (TH1F*)f1->Get("MVA8Jets_JER_Down");
  TH1F *h_JER_Up_8j = (TH1F*)f1->Get("MVA8Jets_JER_Up");
  s = nExTT_El_8j/h_JER_Up_8j->Integral();
  h_JER_Up_8j->Scale(s);
  s = nExTT_El_8j/h_JER_Down_8j->Integral();
  h_JER_Down_8j->Scale(s);
  // h_JER_Up_8j->Write("JER_Up_8j");
  //h_JER_Down_8j->Write("JER_Down_8j");
 

  cout <<"closing files..."<<endl;
 
  // f2->Close();
  f1->Close();

 myfile.close();

    cout <<"end job..."<<endl;


    // cout <<"6J: Scale uncertainty =  "<<  nExTT_El_6j         << "   "<< int_6j_scaleup <<"   "<< int_6j_scaledown<<" ==> " <<(fabs( int_6j_scaledown -  int_6j_scaleup    )/2)/(nExTT_El_6j)  <<endl;
    //  cout <<"7J: Scale uncertainty =  "<<          << "   "<< <<  <<endl;
    //cout <<"8J: Scale uncertainty =  "<<          << "   "<< <<  <<endl;

    //cout <<"6J: Matching uncertainty =  "<<  nExTT_El_6j         << "   "<< int_6j_matchingup <<"   "<< int_6j_matchingdown <<" ==> " <<(fabs( int_6j_matchingdown -  int_6j_matchingup    )/2)/(nExTT_El_6j)  <<endl;
    //  cout <<"7J: Scale uncertainty =  "<<          << "   "<< <<  <<endl;
    //cout <<"8J: Scale uncertainty =  "<<          << "   "<< <<  <<endl;

    //cout <<"6J: JES uncertainty =  "<<  nExTT_El_6j<< "  "<<int_6j_jesup <<" "<<int_6j_jesdown<<" ==> " <<(fabs( int_6j_jesdown -  int_6j_jesup    )/2)/(nExTT_El_6j)  <<endl;


    //cout <<"6J: JER uncertainty =  "<<  nExTT_El_6j<< "   "<< int_6j_jerup <<"   "<< int_6j_jerdown <<" ==> " <<(fabs( int_6j_jerdown -  int_6j_jerup    )/2)/(nExTT_El_6j)  <<endl;

    //cout <<"6J: ttbb uncertainty =  "<<  nExTT_El_6j         << "   "<< int_6j_ttbbup <<"   "<< int_6j_ttbbdown<<" ==> " <<(fabs( int_6j_ttbbdown -  int_6j_ttbbup    )/2)/(nExTT_El_6j)  <<endl;

    //cout <<"6J: PU uncertainty =  "<<  nExTT_El_6j         << "   "<< int_6j_puup <<"   "<< int_6j_pudown<<" ==> " <<(fabs( int_6j_pudown -  int_6j_puup    )/2)/(nExTT_El_6j)  <<endl;

    //cout <<"6J: lepton sf uncertainty =  "<<  nExTT_El_6j         << "   "<< int_6j_leptonsfup <<"   "<< int_6j_leptonsfdown<<" ==> " <<(fabs( int_6j_leptonsfdown -  int_6j_leptonsfup    )/2)/(nExTT_El_6j)  <<endl;

    //cout <<"6J: btag uncertainty =  "<<  nExTT_El_6j         << "   "<< int_6j_btagup <<"   "<< int_6j_btagdown <<" ==> " <<(fabs( int_6j_btagdown -  int_6j_btagup    )/2)/(nExTT_El_6j) <<endl;

    //cout <<"6J: mistag uncertainty =  "<<  nExTT_El_6j         << "   "<< int_6j_mistagup <<"   "<< int_6j_mistagdown <<" ==> " <<(fabs( int_6j_mistagdown -  int_6j_mistagup    )/2)/(nExTT_El_6j)  <<endl;



    double k_factor = 1.022727;

cout<<"Comparing effects of shape systematics..test. "<<endl;
//if (h_Scale_Up)cout<<"Comparing effects of shape systematics... "<<h_Scale_Up->Integral()<<endl;
 cout <<"Scale ="<<nExTT_El<<" "<< int_scaleup <<" "<<  int_scaledown <<" ==>    "<<100*(int_scaledown - nExTT_El )/nExTT_El<<"     "<<100*(int_scaleup - nExTT_El )/nExTT_El   <<  endl;
cout <<"Matching =  "<<nExTT_El << " "<< int_matchingup <<" "<< int_matchingdown <<" ==>    "<<100*(int_matchingdown - nExTT_El )/nExTT_El<<"     "<<100*(int_matchingup - nExTT_El )/nExTT_El     <<endl;
 cout <<"JES =  "<<nExTT_El<< " "<<k_factor*int_jesup <<" "<<k_factor*int_jesdown<<" ==>      "<<100*((k_factor*int_jesdown) - nExTT_El )/nExTT_El<<"     "<<100*((k_factor*int_jesup) - nExTT_El )/nExTT_El     <<endl;
 cout <<"JER =  "<<nExTT_El<< "   "<< k_factor*int_jerup <<"   "<< k_factor*int_jerdown <<" ==>      "<<100*((k_factor*int_jerdown) - nExTT_El )/nExTT_El<<"     "<<100*((k_factor*int_jerup) - nExTT_El )/nExTT_El     <<endl;
 cout <<"ttbb =  "<<  nExTT_El         << "   "<<k_factor*int_ttbbup <<"   "<< k_factor*int_ttbbdown<<" ==>      "<<100*((k_factor*int_ttbbdown) - nExTT_El )/nExTT_El<<"      "<<100*((k_factor*int_ttbbup) - nExTT_El )/nExTT_El     <<endl;
 cout <<"PU =  "<<nExTT_El << " "<<k_factor*int_puup <<"   "<<k_factor*int_pudown<<" ==>  "<<100*((k_factor*int_pudown) - nExTT_El )/nExTT_El<<"      "<<100*((k_factor*int_puup) - nExTT_El )/nExTT_El     <<endl;
 cout <<"lep sf = "<<nExTT_El << "   "<< k_factor*int_leptonsfup <<"   "<<k_factor*int_leptonsfdown<<" ==>      "<<100*((k_factor*int_leptonsfdown) - nExTT_El )/nExTT_El<<"     "<<100*((k_factor*int_leptonsfup) - nExTT_El )/nExTT_El     <<endl;
 cout <<"btag =  "<<nExTT_El << "   "<< k_factor*int_btagup <<"   "<< k_factor*int_btagdown <<" ==>      "<<100*((k_factor*int_btagdown) - nExTT_El )/nExTT_El<<"     "<<100*((k_factor*int_btagup) - nExTT_El )/nExTT_El     <<endl;
 cout <<"mistag =  "<<  nExTT_El<< "   "<< k_factor*int_mistagup <<"   "<< k_factor*int_mistagdown <<" ==>      "<<100*((k_factor*int_mistagdown) - nExTT_El )/nExTT_El<<"     "<<100*((k_factor*int_mistagup) - nExTT_El )/nExTT_El     <<endl;




}
