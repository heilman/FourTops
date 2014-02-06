{
  TFile* f_Mu = new TFile("FourTop_EventSelection_Mu.root"); 

  TH1F * hMET_tttt = (TH1F*)f_Mu->Get("MultiSamplePlot_MET/MET_Data");
  TH1F * hMET_ttll = (TH1F*)f_Mu->Get("MultiSamplePlot_MET/MET_NP_overlay_TTTT");




  hMET_tttt->DrawNormalized();
 hMET_ttll->DrawNormalized("same");
 

}
