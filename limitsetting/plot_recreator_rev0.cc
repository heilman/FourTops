#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"
#include "TPad.h"
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <math.h> 

using namespace std;


int main(){

  TCanvas * c = new TCanvas();
  TPad *Canvas_1 = new TPad("Canvas_1", "Canvas_1",0.01,0.14,0.99,0.99);
  Canvas_1->Draw();
  Canvas_1->cd();
  Canvas_1->SetLogy();
  Canvas_1->SetBorderSize(0);
  //  Canvas_1->SetGridx();
  // Canvas_1->SetGridy();
  Canvas_1->SetFrameBorderMode(0);

  string filename = "FourTopnoRho_EventSelection_El_preApp.root";

  TFile * f = new TFile(filename.c_str());
  THStack * ts = new THStack();

  TH1F * h_Data = (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_Data");

  //ttbar
  TH1F * h_TTJets_ll =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_TTJets_ll");
  TH1F * h_TTJets_bb =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_TTJets_bb");
  TH1F * h_TTJets_cc =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_TTJets_cc");
  //  ts->Add(h_TTJets_cc);
  //ts->Add(h_TTJets_bb);
  //ts->Add(h_TTJets_ll);

  h_TTJets_ll->SetFillColor(kRed);
  h_TTJets_cc->SetFillColor(kRed+1);
  h_TTJets_bb->SetFillColor(kRed+2);




  TH1F * h_EW = (TH1F*) h_TTJets_ll->Clone();
  h_EW->Reset();
  h_EW->SetFillColor(kGreen); 

  //EW
  TH1F * h_W4Jets =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_W_4Jets");
  cout <<"w 4 jets "<<  h_W4Jets->Integral()   <<endl;
  h_EW->Add(h_W4Jets);
  TH1F * h_Z4Jets =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_Z_4Jets");
  cout <<"z 4 jets "<<  h_Z4Jets->Integral()   <<endl;
  h_EW->Add(h_Z4Jets);
  TH1F * h_Z4Jets_lm =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_Z_4Jets_lm");
  h_EW->Add(h_Z4Jets_lm);
  cout <<"w 4 jets "<<  h_Z4Jets_lm->Integral()   <<endl;
  TH1F * h_SingleTop_t_T =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_SingleTop_t_T");
  h_EW->Add(h_SingleTop_t_T);
  cout <<"st t T jets "<<  h_SingleTop_t_T->Integral()   <<endl;
  TH1F * h_SingleTop_t_TBar =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_SingleTop_t_TBar");
  cout <<"st t TBar "<<  h_SingleTop_t_TBar->Integral()   <<endl;
  h_EW->Add(h_SingleTop_t_TBar);
  TH1F * h_SingleTop_s_T =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_SingleTop_s_T");
  cout <<"st s T"<<  h_SingleTop_s_T->Integral()   <<endl;
  h_EW->Add(h_SingleTop_s_T);
  TH1F * h_SingleTop_s_TBar =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_SingleTop_s_TBar");
  cout <<"st s TBar"<<  h_SingleTop_s_TBar->Integral()   <<endl;
   h_EW->Add(h_SingleTop_s_TBar);
  TH1F * h_SingleTop_tW_T =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_SingleTop_tW_T");
  cout <<"st tW T"<<  h_SingleTop_tW_T->Integral()   <<endl;
  h_EW->Add(h_SingleTop_tW_T);
  TH1F * h_SingleTop_tW_TBar =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_SingleTop_tW_TBar");
  cout <<"st tW TBar"<<  h_SingleTop_tW_TBar->Integral()   <<endl;
  h_EW->Add(h_SingleTop_tW_TBar);
  TH1F * h_WW =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_WW");
  h_EW->Add(h_WW);
  cout <<"WW "<<  h_WW->Integral()   <<endl;
  TH1F * h_WZ =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_WZ");
  cout <<"WZ "<<  h_WZ->Integral()   <<endl;
  h_EW->Add(h_WZ);
  TH1F * h_ZZ =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_ZZ");
  //h_EW->Add(h_ZZ);
  cout <<"ZZ "<<  h_ZZ->Integral()   <<endl;


  //tt other
  TH1F * h_TTJets_Other = (TH1F*) h_TTJets_ll->Clone();
  h_TTJets_Other->Reset();
  h_TTJets_Other->SetFillColor(kGray);

  TH1F * h_ttZ =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_ttZ");
  h_TTJets_Other->Add(h_ttZ);
  cout <<"ttZ "<<  h_ttZ->Integral()   <<endl;
  TH1F * h_ttW =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_ttW");
  h_TTJets_Other->Add(h_ttW);
  cout <<"ttW "<<  h_ttW->Integral()   <<endl;
  TH1F * h_ttH =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_ttH");
  h_TTJets_Other->Add(h_ttH);
  cout <<"ttH "<<  h_ttH->Integral()   <<endl;
  cout <<"tt other "<<  h_TTJets_Other->Integral()   <<endl;
  TH1F * h_TTJets_other =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_TTJets_Other");
  h_TTJets_Other->Add(h_TTJets_other);
  cout <<"tt other "<<  h_TTJets_Other->Integral()   <<endl;


  //tttt
  TH1F * h_tttt =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_NP_overlay_TTTT");

  h_tttt->Scale(100.);
  h_tttt->SetLineWidth(3);

  ts->Add(h_TTJets_Other);
  ts->Add(h_EW);  
  ts->Add(h_TTJets_bb);
  ts->Add(h_TTJets_cc);
  ts->Add(h_TTJets_ll);

  ts->SetMaximum(100000);
  ts->SetMinimum(1);


  std::ostringstream strs;
  strs<< std::setprecision(5)  << h_TTJets_ll->Integral() ; 
  std::string str = strs.str();
  string ttll_leg = "tt + ll (" + str + " entries)";

  cout <<"tt + ll "<< str    <<endl;
  strs.str(std::string());

  strs<< std::setprecision(4) << h_TTJets_cc->Integral() ; 
  str = strs.str();
  string ttcc_leg = "tt + cc (" + str + " entries)";

  cout <<"tt + cc "<< str    <<endl;
  strs.str(std::string());
  strs<< std::setprecision(4) << h_TTJets_bb->Integral() ; 
  str = strs.str();
  string ttbb_leg = "tt + bb (" + str + " entries)";


  strs.str(std::string());
  strs<< std::setprecision(4) << h_EW->Integral() ; 
  str = strs.str();
  string EW_leg = "EW (" + str + " entries)";

  strs.str(std::string());
  strs<< std::setprecision(4) << h_TTJets_Other->Integral() ; 
  str = strs.str();
  string TTJets_Other_leg = "tt other (" + str + " entries)";


  strs.str(std::string());
  strs<< std::setprecision(4) << h_Data->Integral(0,16) ; 
  str = strs.str();
  string Data_leg = "Data (" + str + " entries)";

  //Add sys error band.

  TFile* f_Sys= new TFile("../ScaleFilesEl_PreApp/Error_MVA.root");

  TH1F * h_Up = (TH1F*)f_Sys->Get("Up");
  TH1F * h_Down = (TH1F*)f_Sys->Get("Down");


  //  h_EW->Draw("HIST");
  // ts->Draw("HIST");
  // h_Data->Draw("same");
  // h_tttt->Draw("HISTsame");

  //  ts->GetHistogram()->GetXaxis()->SetTickLength(0);
  // ts->GetHistogram()->GetXaxis()->SetLabelOffset(999);

   TH1F * h_stack = (TH1F*)ts->GetStack()->Last();
   TGraphAsymmErrors Errors = TGraphAsymmErrors(h_stack);

   //TAxis *xaxis = h->GetXaxis(); etc.
   
   TH1F * h_diff = (TH1F*)h_Data->Clone();
   TGraphAsymmErrors diff_Errors = TGraphAsymmErrors(h_diff);

   for (int i = 1; i<16; i++){

     double binCenter = h_stack->GetXaxis()->GetBinCenter(i);
     double nom = h_stack->GetBinContent(i); 
     double scale_up =  h_Up->GetBinContent(i);
     double scale_down =  h_Down->GetBinContent(i);
     double eyl, eyh =0;
     double data = h_Data->GetBinContent(i);

     double diff =0;
 
     diff_Errors.SetPoint(i-1, binCenter, 0.);

     if (nom!=0 && data !=0) {diff = 100*((data - nom)/(nom));
     h_diff->SetBinContent(i, diff);
     }

     if(scale_up > scale_down){
            eyh = fabs( scale_up - nom );
            eyl = fabs( scale_down - nom );
            }                   
     else{
         eyl = fabs( scale_up - nom );
         eyh = fabs( scale_down - nom );
                      }

     if (nom !=0){
      Errors.SetPointEYlow(i-1, eyl);
      Errors.SetPointEYhigh(i-1, eyh);
     }

      if(data !=0){
      diff_Errors.SetPointEYlow(i-1, 100*(eyl/nom));
      diff_Errors.SetPointEYhigh(i-1, 100*(eyh/nom));
      }

      // Errors.SetPoint(i,binCenter, nom );
      cout <<"i"<< i <<  "  binCenter "<< binCenter << "  nom "<< nom << " eyl "  << eyl << " eyh "  << eyh <<  endl;
    }


     diff_Errors.SetFillColor(1);
     diff_Errors.SetFillStyle(3005);
     Errors.SetFillColor(1);
     Errors.SetFillStyle(3005);
     // Errors.Draw("2");
     Errors.SetMarkerColor(kViolet);

   //  h_EW->Draw("HIST");
 
   ts->Draw("HIST");
   h_Data->Draw("same");
   Errors.Draw("2");
   // h_Data->Draw("same");
   h_tttt->Draw("HISTsame");

   // Errors.Draw("ALP");
   TFile * f_test = new TFile("niceplot.root","RECREATE");
   ts->Write();
   Errors.Write();

   f_test->Write();
   f_test->Close();
   Canvas_1->Draw();

   //   Canvas_1->Modified();
   //Canvas_1->Update();

   //   ts->GetHistogram()->GetXaxis()->SetTickLength(0);
   // ts->GetHistogram()->GetXaxis()->SetLabelOffset(999);

   TLegend * legend = new TLegend(0.6,0.63,0.85,0.88);
   legend->SetFillColor(0);
   legend->AddEntry(h_TTJets_ll,ttll_leg.c_str(), "f");
   legend->AddEntry(h_TTJets_cc,ttcc_leg.c_str()   , "f");
   legend->AddEntry(h_TTJets_bb,ttbb_leg.c_str()  , "f");
   legend->AddEntry(h_EW, EW_leg.c_str()  , "f");
   legend->AddEntry(h_TTJets_Other,TTJets_Other_leg.c_str()  , "f");
   legend->AddEntry(&Errors ,"Scale uncertainty", "f");
   legend->AddEntry(h_tttt,"SM tttt (X 100)", "l"); 
   legend->AddEntry(h_Data,Data_leg.c_str(), "l"); 
  

   legend->Draw();

   string header = "CMS Preliminary, 19.7 fb^{-1} at #sqrt{s} = 8 TeV";

   TLatex text;
   text.SetNDC(true);
   text.SetTextAlign(12);
   text.SetTextFont(42);
   text.SetTextSize(0.05);
   text.DrawLatex(0.11,0.93,header.c_str());


   TLatex text3;
   text3.SetNDC(true);
   text3.SetTextAlign(12);
   text3.SetTextFont(42);
   text3.SetTextSize(0.03);
   text3.DrawLatex(0.15,0.8,"#font[42]{ #splitline{1 iso. e,  #geq 6 Jets, #geq 2 b-tags}{H_{t} #geq 400 GeV/c, E^{miss}_{t} #geq 30 GeV}}");


    c->cd();
    TPad *Canvas_2 = new TPad("Canvas_2", "Canvas_2",0.01,0.04,0.99,0.23);
    Canvas_2->Draw();
    Canvas_2->cd();
    Canvas_2->SetTopMargin(0);
    Canvas_2->SetBottomMargin(0.3);
    Canvas_2->Range(-0.3875,-0.4930139,0.4875,3.98021);
    Canvas_2->SetFillColor(0);
    Canvas_2->SetBorderMode(0);
    Canvas_2->SetBorderSize(0);
    // Canvas_2->SetGridx();
    Canvas_2->SetGridy();
    Canvas_2->SetFrameBorderMode(0);
    h_diff->SetTitle("");

    h_diff->SetMaximum( 170. );
    h_diff->SetMinimum(-170.);
    h_diff->GetXaxis()->SetLabelSize(0.09);
    h_diff->GetYaxis()->SetLabelSize(0.09);
    h_diff->GetXaxis()->SetTitleSize(0.17);
    h_diff->GetXaxis()->SetTitleOffset(0.75);
    h_diff->GetXaxis()->SetTitle("BDT Discriminant");
    h_diff->GetYaxis()->SetTitle("#frac{Data - MC}{MC} (%)");
    h_diff->GetYaxis()->SetTitleSize(0.13);
    h_diff->GetYaxis()->SetTitleOffset(0.25);
    h_diff->SetMarkerSize(1.5);
    h_diff->GetYaxis()->SetNdivisions(6);

    h_diff->Draw("HIST");
    diff_Errors.Draw("2");

    c->SaveAs("plot_recreator.pdf");

}
