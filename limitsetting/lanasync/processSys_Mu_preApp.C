#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <iostream>
#include <fstream>
#include "THStack.h"

void process(){

 /////////////////////////
 // F I L E S
 /////////////////////////
 TFile* f_Mu = new TFile("FourTop_EventSelection_wMETCut_Mu_preApp.root"); 
 TFile* f1 = new TFile("SystematicShapes_Mu_preApp.root");


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

  TH1F * ttMu = new TH1F("","", 15,-0.3, 0.4);

  TH1F * hTTll_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_TTJets_ll");
  TH1F * hTTcc_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_TTJets_cc");
  TH1F * hTTbb_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_TTJets_bb");
  TH1F * hTT_Mu   = (TH1F*)hTTbb_Mu->Clone();
  hTT_Mu->Reset();
  hTT_Mu->Add(hTTll_Mu);
  hTT_Mu->Add(hTTcc_Mu);
  hTT_Mu->Add(hTTbb_Mu);


  hTT_Mu->Scale(0.995); //correcting for pixel lumi
 
  cout <<"n TT ll  " << hTTll_Mu->Integral(0,16)<< endl;
  cout <<"n TT cc  " << hTTcc_Mu->Integral(0,16)<< endl;
  cout <<"n TT bb  " << hTTbb_Mu->Integral(0,16)<< endl;
 
 
  double nExTT_Mu  = hTT_Mu->Integral(0,16);

  hTT_Mu->Write("tt");
  hSig_Mu->Write("tttt");
  hData_Mu->Write("data");
  hdata_100fbtttt->Write("data_w100fbtttt");
  hdata_150fbtttt->Write("data_w150fbtttt");
  hdata_200fbtttt->Write("data_w200fbtttt");

  //  cout << hTTll_Mu->Integral(0,16) <<endl;
  //cout << hTTcc_Mu->Integral(0,16) <<endl;
  // cout << hTTbb_Mu->Integral(0,16) <<endl;
  //  cout << nExTT_Mu<<endl;

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


  hEW_Mu->Scale(0.995); //correcting for pixel lumi
  hEW_Mu->Write("EW");



  //tt other histo
  TH1F * httW_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_ttW");
  TH1F * httZ_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_ttZ");
  TH1F * httH_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_ttH");
  TH1F * hTTOther_Mu = (TH1F*)f_Mu->Get("MultiSamplePlot_MVA/MVA_TTJets_Other");

  hTTOther_Mu->Add(httW_Mu);
  hTTOther_Mu->Add(httZ_Mu);
  hTTOther_Mu->Add(httH_Mu);

  hTTOther_Mu->Scale(0.995); //correcting for pixel lumi
  hTTOther_Mu->Write("ttOther");


  fnom->Close();

  ///////////////////////////////////////////////////////
  //// P R O C E S S   S Y S T E M A T I C   H I S T O S
  ///////////////////////////////////////////////////////

 TFile* f2 = new TFile("SystematicShapes_norm_Mu.root", "RECREATE");

 ///////////////////////////////
 /// I N C L U S I V E 
 /////////////////////////////////

ofstream myfile;
  myfile.open ("chi2_systematics_mu_preApp.txt");
  myfile << "Chi2/NDOF of systematics:"<< endl;

  TH1F *  h_Scale_Down = (TH1F*)f1->Get("Scale_Down");
  TH1F * h_Scale_Up = (TH1F*)f1->Get("Scale_Up");

 double int_scaledown = h_Scale_Down->Integral(0,16);
 double int_scaleup = h_Scale_Up->Integral(0,16);


  double s = nExTT_Mu/h_Scale_Up->Integral(0,16);
  h_Scale_Up->Scale(s);
  s = nExTT_Mu/h_Scale_Down->Integral(0,16);
  h_Scale_Down->Scale(s);
  h_Scale_Up->Write("Scale_Up");
  h_Scale_Down->Write("Scale_Down");

  /*
  TH1F *  h_Scale_Down_1stHalf = (TH1F*)f1->Get("Scale_Down_1stHalf");
  TH1F * h_Scale_Up_1stHalf = (TH1F*)f1->Get("Scale_Up_1stHalf");
  s = nExTT_Mu/h_Scale_Up_1stHalf->Integral(0,16);
  h_Scale_Up_1stHalf->Scale(s);
  s = nExTT_Mu/h_Scale_Down_1stHalf->Integral(0,16);
  h_Scale_Down_1stHalf->Scale(s);
  h_Scale_Up_1stHalf->Write("Scale_Up_1stHalf");
  h_Scale_Down_1stHalf->Write("Scale_Down_1stHalf");

 

  TH1F *  h_Scale_Down_2ndHalf = (TH1F*)f1->Get("Scale_Down_2ndHalf");
  TH1F * h_Scale_Up_2ndHalf = (TH1F*)f1->Get("Scale_Up_2ndHalf");
  s = nExTT_Mu/h_Scale_Up_2ndHalf->Integral(0,16);
  h_Scale_Up_2ndHalf->Scale(s);
  s = nExTT_Mu/h_Scale_Down_2ndHalf->Integral(0,16);
  h_Scale_Down_2ndHalf->Scale(s);
  h_Scale_Up_2ndHalf->Write("Scale_Up_2ndHalf");
  h_Scale_Down_2ndHalf->Write("Scale_Down_2ndHalf");
  */

  TH1F * h_Matching_Down = (TH1F*)f1->Get("Matching_Down");
  TH1F *h_Matching_Up = (TH1F*)f1->Get("Matching_Up");

  double int_matchingdown = h_Matching_Down->Integral(0,16);
  double int_matchingup = h_Matching_Up->Integral(0,16);

  s = nExTT_Mu/h_Matching_Up->Integral(0,16);
  h_Matching_Up->Scale(s);
  s = nExTT_Mu/h_Matching_Down->Integral(0,16);
  h_Matching_Down->Scale(s);
  h_Matching_Up->Write("Matching_Up");
  h_Matching_Down->Write("Matching_Down");


  TH1F *h_JES_Down = (TH1F*)f1->Get("JES_Down");
  TH1F *h_JES_Up = (TH1F*)f1->Get("JES_Up");

 double int_jesdown = h_JES_Down->Integral(0,16);
 double int_jesup = h_JES_Up->Integral(0,16);

  s = nExTT_Mu/h_JES_Up->Integral(0,16);
  h_JES_Up->Scale(s);
  s = nExTT_Mu/h_JES_Down->Integral(0,16);
  h_JES_Down->Scale(s);
  h_JES_Up->Write("JES_Up");
  h_JES_Down->Write("JES_Down");


  TH1F * h_ttbb_Down = (TH1F*)f1->Get("ttbb_Down");
   TH1F *h_ttbb_Up = (TH1F*)f1->Get("ttbb_Up");

 double int_ttbbdown = h_ttbb_Down->Integral(0,16);
 double int_ttbbup = h_ttbb_Up->Integral(0,16);

 s = nExTT_Mu/h_ttbb_Up->Integral(0,16);
  h_ttbb_Up->Scale(s);
  s = nExTT_Mu/h_ttbb_Down->Integral(0,16);
  h_ttbb_Down->Scale(s);
  h_ttbb_Up->Write("ttbb_Up");
  h_ttbb_Down->Write("ttbb_Down");


  TH1F * h_bTag_Down = (TH1F*)f1->Get("bTag_Down");
   TH1F * h_bTag_Up = (TH1F*)f1->Get("bTag_Up");

 double int_btagdown = h_bTag_Down->Integral(0,16);
 double int_btagup = h_bTag_Up->Integral(0,16);

  s = nExTT_Mu/h_bTag_Up->Integral(0,16);
  h_bTag_Up->Scale(s);
  s = nExTT_Mu/h_bTag_Down->Integral(0,16);
  h_bTag_Down->Scale(s);
  h_bTag_Up->Write("bTag_Up");
  h_bTag_Down->Write("bTag_Down");


   TH1F * h_misTag_Down = (TH1F*)f1->Get("misTag_Down");
   TH1F * h_misTag_Up = (TH1F*)f1->Get("misTag_Up");

 double int_mistagdown = h_misTag_Down->Integral(0,16);
 double int_mistagup = h_misTag_Up->Integral(0,16);

  s = nExTT_Mu/h_misTag_Up->Integral(0,16);
  h_misTag_Up->Scale(s);
  s = nExTT_Mu/h_misTag_Down->Integral(0,16);
  h_misTag_Down->Scale(s);
  h_misTag_Up->Write("misTag_Up");
  h_misTag_Down->Write("misTag_Down");


  TH1F * h_leptonSF_Down = (TH1F*)f1->Get("leptonSF_Down");
   TH1F * h_leptonSF_Up = (TH1F*)f1->Get("leptonSF_Up");

 double int_leptonsfdown = h_leptonSF_Down->Integral(0,16);
 double int_leptonsfup = h_leptonSF_Up->Integral(0,16);

  s = nExTT_Mu/h_leptonSF_Up->Integral(0,16);
  h_leptonSF_Up->Scale(s);
  s = nExTT_Mu/h_leptonSF_Down->Integral(0,16);
  h_leptonSF_Down->Scale(s);
  h_leptonSF_Up->Write("leptonSF_Up");
  h_leptonSF_Down->Write("leptonSF_Down");


  TH1F *h_PU_Down = (TH1F*)f1->Get("PU_Down");
  TH1F *h_PU_Up = (TH1F*)f1->Get("PU_Up");
 double int_pudown = h_PU_Down->Integral(0,16);
 double int_puup = h_PU_Up->Integral(0,16);

  s = nExTT_Mu/h_PU_Up->Integral(0,16);
  h_PU_Up->Scale(s);
  s = nExTT_Mu/h_PU_Down->Integral(0,16);
  h_PU_Down->Scale(s);
  h_PU_Up->Write("PU_Up");
  h_PU_Down->Write("PU_Down");



  TH1F *h_JER_Down = (TH1F*)f1->Get("JER_Down");
  TH1F *h_JER_Up = (TH1F*)f1->Get("JER_Up");
  double int_jerdown = h_JER_Down->Integral(0,16);
  double int_jerup = h_JER_Up->Integral(0,16);

  s = nExTT_Mu/h_JER_Up->Integral(0,16);
  h_JER_Up->Scale(s);
  s = nExTT_Mu/h_JER_Down->Integral(0,16);
  h_JER_Down->Scale(s);
  h_JER_Up->Write("JER_Up");
  h_JER_Down->Write("JER_Down");

  TCanvas* ctest = new TCanvas();

  // h_Matching_Up->Draw();
  cout<<"HuDmo" << hTT_Mu->Integral()<<endl;

 myfile <<"Systematic  & " <<"#chi^{2} Down  "<< " & #chi^{2} Up  "<<endl;
 myfile <<"Scale  & "<< h_Scale_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_Scale_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF")<<endl;
 myfile <<"Matching &  "<< h_Matching_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_Matching_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF")<<endl;
 myfile <<"JES  & "<< h_JES_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_JES_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF") << endl;
 myfile <<"ttbb  & "<< h_ttbb_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_ttbb_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF") << endl;
 myfile <<"btag &  "<< h_bTag_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_bTag_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF")    << endl;
 myfile <<"mistag &  "<< h_misTag_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_misTag_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF") << endl;
 myfile <<"lepton SF &  "<< h_leptonSF_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_leptonSF_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF") << endl;
 myfile <<"PU   &  "<< h_PU_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_PU_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<  endl;
 myfile <<"JER &  "<< h_JER_Down->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<" &  "<< h_JER_Up->Chi2Test(hTT_Mu,"WWCHI2/NDF") <<  endl;

 //h_Matching_Up->Draw();
 //hTT_Mu->Draw("same");




  //nnlo corr, ttbb corr, int lumi corr, arbitrary corr
 
  double k_factor = 1.023*1.033*1.045*1.01569;

  // double k_factor = 1.0;

cout<<"Comparing effects of shape systematics..test. "<<endl;
//if (h_Scale_Up)cout<<"Comparing effects of shape systematics... "<<h_Scale_Up->Integral(0,16)<<endl;
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

  cout <<"end job..."<<endl;

}
