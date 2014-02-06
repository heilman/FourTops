#include "TFile.h"
#include "TH1F.h"
#include "THStack.h"
#include "TCanvas.h"
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>

// macro to create shape systeamatic histograms corresponding to
//linear and parabolic deformations of the BDT distributions.

using namespace std;

void test(){


  // string filename = "FourTopnoRho_EventSelection_El_preApp.root";
  //string filename_sys = "SystematicShapes_norm_El.root";

  string filename = "FourTop_EventSelection_wMETCut_Mu_preApp.root";
  string filename_sys = "SystematicShapes_norm_Mu.root";



  TFile * f = new TFile(filename.c_str());
  TFile * f_sys = new TFile(filename_sys.c_str());

  cout <<"here 1 "<<endl;

  TH1F * h_Data =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_Data");
cout <<"here 1.1 "<<endl;
  TH1F * h_TTJets_ll =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_TTJets_ll");
cout <<"here 1.2 "<<endl;
  TH1F * h_TTJets_bb =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_TTJets_bb");
cout <<"here 1.3 "<<endl;
  TH1F * h_TTJets_cc =  (TH1F*)f->Get("MultiSamplePlot_MVA/MVA_TTJets_cc");
cout <<"here 1.4 "<<endl;


  TH1F * h_Scale_Up = (TH1F*)f_sys->Get("Scale_Up");
cout <<"here 1.5 "<<endl;
  TH1F * h_Scale_Down = (TH1F*)f_sys->Get("Scale_Down");

cout <<"here 1.6 "<<endl;

  h_TTJets_ll->Add(h_TTJets_bb);
  h_TTJets_ll->Add(h_TTJets_cc);

  cout <<"here 2.1 "<<endl;

  TCanvas * c = new TCanvas();
  cout <<"here 2.2 "<<endl;

  //  h_Data->Draw();
  cout <<"here 2.3 "<<endl;
  //  h_TTJets_ll->Draw();
  // h_Scale_Up->Draw("same");
  // h_Scale_Down->Draw("same");

  cout <<"here 2.4 "<<endl;

  //    h_TTJets_ll->SetFillColor(kRed);
  h_Scale_Down->SetMarkerColor(kGreen);
  h_Scale_Up->SetMarkerColor(kBlue);

  TH1F * h_linear_up = (TH1F*)h_TTJets_ll->Clone();
  h_linear_up->Reset();

  TH1F * h_linear_down = (TH1F*)h_TTJets_ll->Clone();
  h_linear_down->Reset();

  TH1F * h_parabolic_up = (TH1F*)h_TTJets_ll->Clone();
  h_parabolic_up->Reset();

  TH1F * h_parabolic_down = (TH1F*)h_TTJets_ll->Clone();
  h_parabolic_down->Reset();


  TH1F * h_linear_frac_up = (TH1F*)h_TTJets_ll->Clone();
  h_linear_frac_up->Reset();

  TH1F * h_linear_frac_down = (TH1F*)h_TTJets_ll->Clone();
  h_linear_frac_down->Reset();

  TH1F * h_parabolic_frac_up = (TH1F*)h_TTJets_ll->Clone();
  h_parabolic_frac_up->Reset();

  TH1F * h_parabolic_frac_down = (TH1F*)h_TTJets_ll->Clone();
  h_parabolic_frac_down->Reset();



  //h_linear_down->SetFillColor(kBlue);

  for (int i =1; i< 16; i++){

    double tt = h_TTJets_ll->GetBinContent(i);
    double s_up = h_Scale_Up->GetBinContent(i);
    double s_down = h_Scale_Down->GetBinContent(i);
     if(tt!=0.){
    //    double trend =  100*(fabs(s_up - s_down)/(2*tt));
 
	  //    double trend = (tt*(0.018*i) - (.07*tt))/ tt;
 
    double trend = (1.42*(i-1)) - 10;
    h_linear_frac_up->SetBinContent(i,trend);

    cout<<"bin " << i <<"trend "<< trend  <<" tt "<<tt <<" correction "<< tt + (trend/tt*100)   <<endl;
    h_linear_up->SetBinContent(i,tt + ((tt*trend)/100.));

    //    trend = (tt*(-0.018*i) + (.07*tt))/ tt;
  
   trend =  (-1.42*(i-1)) + 10;

 
   h_linear_down->SetBinContent(i,tt + ((tt*trend)/100.));
   h_linear_frac_down->SetBinContent(i,trend);

   double t = i-7;
   double trend_parabolic = -((t*t)  +(t)  -20);

   h_parabolic_up->SetBinContent(i,tt + ((tt*trend_parabolic)/100.));
   h_parabolic_frac_up->SetBinContent(i,trend_parabolic);
 
   trend_parabolic = (t*t)  +(t)  -20;

   h_parabolic_down->SetBinContent(i,tt + ((tt*trend_parabolic)/100.) );

   h_parabolic_frac_down->SetBinContent(i,trend_parabolic);



     }

    // y = mx + c
    // m = (20-7)/15


}

  double scale_up = h_TTJets_ll->Integral()/h_linear_up->Integral();
  double scale_down = h_TTJets_ll->Integral()/h_linear_down->Integral();

  double scale_up_para = h_TTJets_ll->Integral()/h_parabolic_up->Integral();
  double scale_down_para = h_TTJets_ll->Integral()/h_parabolic_down->Integral();

  h_linear_up->Scale(scale_up);
  h_linear_down->Scale(scale_down);

  h_parabolic_up->Scale(scale_up_para);
  h_parabolic_down->Scale(scale_down_para);

  h_parabolic_up->SetLineColor(kRed);
  h_parabolic_down->SetLineColor(kBlue);

  h_parabolic_frac_up->SetLineColor(kRed);
  h_parabolic_frac_down->SetLineColor(kBlue);

  h_linear_frac_up->SetLineColor(kRed);
  h_linear_frac_down->SetLineColor(kBlue);


  TCanvas * c1 = new TCanvas();

  h_TTJets_ll->Draw();
  h_linear_up->Draw("same");
  h_linear_down->Draw("same");

  h_linear_up->SetLineColor(kRed);
  h_linear_down->SetLineColor(kBlue);

  c1->SaveAs("linear_deform_mu.pdf");

  TCanvas * c2 = new TCanvas();
  h_TTJets_ll->Draw();
  h_parabolic_down->Draw("same");
  h_parabolic_up->Draw("same");

  c2->SaveAs("parabolic_deform_mu.pdf");


  TCanvas * c3 = new TCanvas();

  h_parabolic_frac_down->SetMinimum(-40.0);
  h_parabolic_frac_down->SetMaximum(40.0);
  h_parabolic_frac_down->SetYTitle("% change");
  h_parabolic_frac_down->Draw();
  h_parabolic_frac_up->Draw("same");

  c3->SaveAs("parabolic_frac_mu.pdf");


  TCanvas * c4 = new TCanvas();
  h_linear_frac_down->SetMinimum(-20.0);
  h_linear_frac_down->SetMaximum(20.0);
  h_linear_frac_down->SetYTitle("% change");
  h_linear_frac_down->Draw();
  h_linear_frac_up->Draw("same");

  c4->SaveAs("linear_frac_mu.pdf");



  TFile * f1 = new TFile("extra_shapes_mu.root", "RECREATE");
  h_parabolic_down->Write("parabolic_down");
  h_parabolic_up->Write("parabolic_up");
  h_linear_down->Write("linear_down");
  h_linear_up->Write("linear_up");

  f1->Write();

}
