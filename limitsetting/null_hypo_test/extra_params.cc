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


  string filename = "FourTopnoRho_EventSelection_El_preApp.root";
  string filename_sys = "SystematicShapes_norm_El.root";


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

  h_TTJets_ll->SetFillColor(kRed);
  h_Scale_Down->SetMarkerColor(kGreen);
  h_Scale_Up->SetMarkerColor(kBlue);

  TH1F * h_linear_up = (TH1F*)h_TTJets_ll->Clone();
  h_linear_up->Reset();

  TH1F * h_linear_down = (TH1F*)h_TTJets_ll->Clone();
  h_linear_down->Reset();


h_linear_down->SetFillColor(kBlue);

  for (int i =1; i< 16; i++){

    double tt = h_TTJets_ll->GetBinContent(i);
    double s_up = h_Scale_Up->GetBinContent(i);
    double s_down = h_Scale_Down->GetBinContent(i);
        if(tt!=0.){
    //    double trend =  100*(fabs(s_up - s_down)/(2*tt));
 
    double trend = (tt*(0.018*i) - (.07*tt))/ tt;
    h_linear_up->SetBinContent(i,trend);

    trend = (tt*(-0.018*i) + (.07*tt))/ tt;
     h_linear_down->SetBinContent(i,trend);

      }

    // y = mx + c
    // m = (20-7)/15


}




  h_linear_up->Draw();
  h_linear_down->Draw("same");

  cout <<"here 3 "<<endl;

}
