{

  gStyle->SetOptStat(0000);
  TFile *f1 = new TFile("FourTop_EventSelection_El_bbMinus.root");
  TFile *f2 = new TFile("FourTop_EventSelection_El_bbPlus.root");
    
  h1 = (TH1F*)f1.Get("MultiSamplePlot_NbOfSelectedJets/NbOfSelectedJets_TTJets");
  h2 = (TH1F*)f2.Get("MultiSamplePlot_NbOfSelectedJets/NbOfSelectedJets_TTJets");


  //Colors, text                                                                                                                                           
  //   h1.SetMarkerColor(kRed);
  //  h1.SetLineColor(kRed);

  //  h1.SetFillColor(kRed);
  // h1.SetFillStyle(3005);


    
  h1.SetTitle("");
  h1.SetXTitle("NJets");
  h2.SetLineColor(kRed);
  // h1.SetFillColor(kRed);
  // h1.SetFillStyle(3005);
    
  //  h2.SetLineColor(kBlue);
  // h2.SetFillColor(kBlue);
  //h2.SetFillStyle(3005);
    
    
  TAxis* xax = h1.GetXaxis();
  TAxis* yax = h1.GetYaxis();
  xax->SetTitleSize(0.045);
  xax->SetTitleOffset(1.);

  //6519.8    
double scale1 =   6519.8/h1.Integral();
double scale2 =  6519.8/h2.Integral();

 h1.Scale(scale1);
 h2.Scale(scale2);


  h1.Draw("HIST");
  h2.Draw("HISTsame");

}
