#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>


void process(){

 /////////////////////////
 // F I L E S
 /////////////////////////

 TFile* f_El = new TFile("FourTopnoRho_EventSelection_El_preApp.root"); 
 TFile* f1 = new TFile("SystematicShapes_El_preApp.root"); 
 // TFile* f2 = new TFile("SystematicShapes_El_preApp_tttt.root"); 


////////////////////////////////////////////////
// P R O C E S S   N O M I N A L   S H A P E S
////////////////////////////////////////////////

 TFile* fnom = new TFile("NominalShapes_El.root","RECREATE"); 


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

  double nExTT_El  = hTT_El->Integral(0,16);

  hTT_El->Scale(0.995);

  hTT_El->Write("tt");
  hSig_El->Write("tttt");
  hData_El->Write("data");
  hdata_100fbtttt->Write("data_w100fbtttt");
  hdata_150fbtttt->Write("data_w150fbtttt");
  hdata_200fbtttt->Write("data_w200fbtttt");

  //  cout << hTTll_El->Integral(0,16) <<endl;
  //cout << hTTcc_El->Integral(0,16) <<endl;
  //cout << hTTbb_El->Integral(0,16) <<endl;
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

  hEW_El->Scale(0.995);

  hEW_El->Write("EW");

  //tt other histo
  TH1F * httW_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_ttW");
  TH1F * httZ_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_ttZ");
  TH1F * httH_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_ttH");
  TH1F * hTTOther_El = (TH1F*)f_El->Get("MultiSamplePlot_MVA/MVA_TTJets_Other");

  hTTOther_El->Add(httW_El);
  hTTOther_El->Add(httZ_El);
  hTTOther_El->Add(httH_El);

  hTTOther_El->Scale(0.995);

  hTTOther_El->Write("ttOther");



  fnom->Close();
  cout<<"processing sytematic shapes..."<<endl;

///////////////////////////////////////////////////////
//// P R O C E S S   S Y S T E M A T I C   H I S T O S
///////////////////////////////////////////////////////

 TFile* f2 = new TFile("SystematicShapes_norm_El.root", "RECREATE"); 

 ///////////////////////
 /// I N C L U S I V E 
 ///////////////////////

 TH1F *h_Scale_Down = (TH1F*)f1->Get("Scale_Down");
 TH1F *h_Scale_Up = (TH1F*)f1->Get("Scale_Up");
 double int_scaledown = h_Scale_Down->Integral(0,16);
 double int_scaleup = h_Scale_Up->Integral(0,16);

 double s = nExTT_El/h_Scale_Up->Integral(0,16);
 h_Scale_Up->Scale(s);
 s = nExTT_El/h_Scale_Down->Integral(0,16);
 h_Scale_Down->Scale(s);
 h_Scale_Up->Write("Scale_Up");
 h_Scale_Down->Write("Scale_Down");

 cout<<"processing sytematic shapes..scale done."<<endl;


 /*
 TH1F *h_Scale_Down_1stHalf = (TH1F*)f1->Get("MVA_Scale_Down_1stHalf");
 TH1F *h_Scale_Up_1stHalf = (TH1F*)f1->Get("MVA_Scale_Up_1stHalf");
 double int_scaledown_1stHalf = h_Scale_Down_1stHalf->Integral(0,16);
 double int_scaleup_1stHalf = h_Scale_Up_1stHalf->Integral(0,16);

 s = nExTT_El/h_Scale_Up_1stHalf->Integral(0,16);
 h_Scale_Up_1stHalf->Scale(s);
 s = nExTT_El/h_Scale_Down_1stHalf->Integral(0,16);
 h_Scale_Down_1stHalf->Scale(s);
 h_Scale_Up_1stHalf->Write("Scale_Up_1stHalf");
 h_Scale_Down_1stHalf->Write("Scale_Down_1stHalf");


 TH1F *h_Scale_Down_2ndHalf = (TH1F*)f1->Get("MVA_Scale_Down_2ndHalf");
 TH1F *h_Scale_Up_2ndHalf = (TH1F*)f1->Get("MVA_Scale_Up_2ndHalf");
 double int_scaledown_2ndHalf = h_Scale_Down_2ndHalf->Integral(0,16);
 double int_scaleup_2ndHalf = h_Scale_Up_2ndHalf->Integral(0,16);

 s = nExTT_El/h_Scale_Up_2ndHalf->Integral(0,16);
 h_Scale_Up_2ndHalf->Scale(s);
 s = nExTT_El/h_Scale_Down_2ndHalf->Integral(0,16);
 h_Scale_Down_2ndHalf->Scale(s);
 h_Scale_Up_2ndHalf->Write("Scale_Up_2ndHalf");
 h_Scale_Down_2ndHalf->Write("Scale_Down_2ndHalf");


  cout<<"scale done....."<<endl;
 */



  cout<<"scale done....."<<endl;

 TH1F *h_Matching_Down = (TH1F*)f1->Get("Matching_Down");
 TH1F *h_Matching_Up = (TH1F*)f1->Get("Matching_Up");
 double int_matchingdown = h_Matching_Down->Integral(0,16);
 double int_matchingup = h_Matching_Up->Integral(0,16);


 s = nExTT_El/h_Matching_Up->Integral(0,16);
 h_Matching_Up->Scale(s);
 s = nExTT_El/h_Matching_Down->Integral(0,16);
 h_Matching_Down->Scale(s);
 h_Matching_Up->Write("Matching_Up");
 h_Matching_Down->Write("Matching_Down");

  cout<<"matching done....."<<endl;

  TH1F *h_JES_Down = (TH1F*)f1->Get("JES_Down");
  TH1F *h_JES_Up = (TH1F*)f1->Get("JES_Up");
  double int_jesdown = h_JES_Down->Integral(0,16);
  double int_jesup = h_JES_Up->Integral(0,16);

  s = nExTT_El/h_JES_Up->Integral(0,16);
  h_JES_Up->Scale(s);
  s = nExTT_El/h_JES_Down->Integral(0,16);
  h_JES_Down->Scale(s);
  h_JES_Up->Write("JES_Up");
  h_JES_Down->Write("JES_Down");

  cout<<"JES done....."<<endl;


  TH1F *h_ttbb_Down = (TH1F*)f1->Get("ttbb_Down");
  TH1F *h_ttbb_Up = (TH1F*)f1->Get("ttbb_Up");
  double int_ttbbdown = h_ttbb_Down->Integral(0,16);
  double int_ttbbup = h_ttbb_Up->Integral(0,16);
  s = nExTT_El/h_ttbb_Up->Integral(0,16);
  h_ttbb_Up->Scale(s);
  s = nExTT_El/h_ttbb_Down->Integral(0,16);
  h_ttbb_Down->Scale(s);
  h_ttbb_Up->Write("ttbb_Up");
  h_ttbb_Down->Write("ttbb_Down");

  cout<<"ttbb done....."<<endl;

  TH1F *h_bTag_Down = (TH1F*)f1->Get("bTag_Down");
  TH1F *h_bTag_Up = (TH1F*)f1->Get("bTag_Up");
  double int_btagdown = h_bTag_Down->Integral(0,16);
  double int_btagup = h_bTag_Up->Integral(0,16);

  s = nExTT_El/h_bTag_Up->Integral(0,16);
  h_bTag_Up->Scale(s);
  s = nExTT_El/h_bTag_Down->Integral(0,16);
  h_bTag_Down->Scale(s);
  h_bTag_Up->Write("bTag_Up");
  h_bTag_Down->Write("bTag_Down");

  cout<<"btag done....."<<endl;


  TH1F *h_misTag_Down = (TH1F*)f1->Get("misTag_Down");
  TH1F *h_misTag_Up = (TH1F*)f1->Get("misTag_Up");
 double int_mistagdown = h_misTag_Down->Integral(0,16);
  double int_mistagup = h_misTag_Up->Integral(0,16);

  s = nExTT_El/h_misTag_Up->Integral(0,16);
  h_misTag_Up->Scale(s);
  s = nExTT_El/h_misTag_Down->Integral(0,16);
  h_misTag_Down->Scale(s);
  h_misTag_Up->Write("misTag_Up");
  h_misTag_Down->Write("misTag_Down");

  cout<<"mistag done....."<<endl;

  TH1F *h_leptonSF_Down = (TH1F*)f1->Get("leptonSF_Down");
  TH1F *h_leptonSF_Up = (TH1F*)f1->Get("leptonSF_Up");

  double int_leptonsfdown = h_leptonSF_Down->Integral(0,16);
  double int_leptonsfup = h_leptonSF_Up->Integral(0,16);

  s = nExTT_El/h_leptonSF_Up->Integral(0,16);
  h_leptonSF_Up->Scale(s);
  s = nExTT_El/h_leptonSF_Down->Integral(0,16);
  h_leptonSF_Down->Scale(s);
  h_leptonSF_Up->Write("leptonSF_Up");
  h_leptonSF_Down->Write("leptonSF_Down");


  cout<<"lepton sf done....."<<endl;

  TH1F *h_PU_Down = (TH1F*)f1->Get("PU_Down");
  TH1F *h_PU_Up = (TH1F*)f1->Get("PU_Up");
  double int_pudown = h_PU_Down->Integral(0,16);
  double int_puup = h_PU_Up->Integral(0,16);

  s = nExTT_El/h_PU_Up->Integral(0,16);
  h_PU_Up->Scale(s);
  s = nExTT_El/h_PU_Down->Integral(0,16);
  h_PU_Down->Scale(s);
  h_PU_Up->Write("PU_Up");
  h_PU_Down->Write("PU_Down");


  cout<<"PU  done....."<<endl;

  TH1F *h_JER_Down = (TH1F*)f1->Get("JER_Down");
  TH1F *h_JER_Up = (TH1F*)f1->Get("JER_Up");
  double int_jerdown = h_JER_Down->Integral(0,16);
  double int_jerup = h_JER_Up->Integral(0,16);
  s = nExTT_El/h_JER_Up->Integral(0,16);
  h_JER_Up->Scale(s);
  s = nExTT_El/h_JER_Down->Integral(0,16);
  h_JER_Down->Scale(s);
  h_JER_Up->Write("JER_Up");
  h_JER_Down->Write("JER_Down");



  cout <<"closing files..."<<endl;
 
  f2->Close();
  f1->Close();

    cout <<"end job..."<<endl;



  double k_factor_1 = 2*1.022727*1.034*1.022525;
  double k_factor_2 = (k_factor_1/2.)*0.9644;
  double k_factor_3 =  1.; 
    // double k_factor = 2.0;


cout<<"Comparing effects of shape systematics..test. "<<endl;
//if (h_Scale_Up)cout<<"Comparing effects of shape systematics... "<<h_Scale_Up->Integral(0,16)<<endl;
 cout <<"Scale ="<<nExTT_El<<" "<< int_scaleup <<" "<<  int_scaledown <<" ==>    "<<100*((k_factor_3*int_scaledown) - nExTT_El )/nExTT_El<<"     "<<100*((k_factor_3*int_scaleup) - nExTT_El )/nExTT_El   <<  endl;
 cout <<"Matching =  "<<nExTT_El << " "<< int_matchingup <<" "<< int_matchingdown <<" ==>    "<<100*((k_factor_3*int_matchingdown) - nExTT_El )/nExTT_El<<"     "<<100*((k_factor_3*int_matchingup) - nExTT_El )/nExTT_El     <<endl;
 cout <<"JES =  "<<nExTT_El<< " "<<k_factor_1*int_jesup <<" "<<k_factor_1*int_jesdown<<" ==>      "<<100*((k_factor_1*int_jesdown) - nExTT_El )/nExTT_El<<"     "<<100*((k_factor_1*int_jesup) - nExTT_El )/nExTT_El     <<endl;
 cout <<"JER =  "<<nExTT_El<< "   "<< k_factor_2*int_jerup <<"   "<< k_factor_2*int_jerdown <<" ==>      "<<100*((k_factor_2*int_jerdown) - nExTT_El )/nExTT_El<<"     "<<100*((k_factor_2*int_jerup) - nExTT_El )/nExTT_El     <<endl;
 cout <<"ttbb =  "<<  nExTT_El         << "   "<<k_factor_1*int_ttbbup <<"   "<< k_factor_1*int_ttbbdown<<" ==>      "<<100*((k_factor_1*int_ttbbdown) - nExTT_El )/nExTT_El<<"      "<<100*((k_factor_1*int_ttbbup) - nExTT_El )/nExTT_El     <<endl;
 cout <<"PU =  "<<nExTT_El << " "<<k_factor_2*int_puup <<"   "<<k_factor_2*int_pudown<<" ==>  "<<100*((k_factor_2*int_pudown) - nExTT_El )/nExTT_El<<"      "<<100*((k_factor_2*int_puup) - nExTT_El )/nExTT_El     <<endl;
 cout <<"lep sf = "<<nExTT_El << "   "<< k_factor_1*int_leptonsfup <<"   "<<k_factor_1*int_leptonsfdown<<" ==>      "<<100*((k_factor_1*int_leptonsfdown) - nExTT_El )/nExTT_El<<"     "<<100*((k_factor_1*int_leptonsfup) - nExTT_El )/nExTT_El     <<endl;
 cout <<"btag =  "<<nExTT_El << "   "<< k_factor_1*int_btagup <<"   "<< k_factor_1*int_btagdown <<" ==>      "<<100*((k_factor_1*int_btagdown) - nExTT_El )/nExTT_El<<"     "<<100*((k_factor_1*int_btagup) - nExTT_El )/nExTT_El     <<endl;
 cout <<"mistag =  "<<  nExTT_El<< "   "<< k_factor_1*int_mistagup <<"   "<< k_factor_1*int_mistagdown <<" ==>      "<<100*((k_factor_1*int_mistagdown) - nExTT_El )/nExTT_El<<"     "<<100*((k_factor_1*int_mistagup) - nExTT_El )/nExTT_El     <<endl;






}
