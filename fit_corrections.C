{

 TFile *f = new TFile("FourTop_EventSelection_Mu.root");
   
 TH1F *h_ttjets = (TH1F*)f->Get("Histos1D/nTagsTTJets");
 TH1F *h_data = (TH1F*)f->Get("Histos1D/nTagsData");

  h_ttjets->Sumw2();
  h_data->Sumw2();

  // h_data->Divide(h_ttjets);

  // h_data->Draw("E1p");

  TH2D *h_TagsHT_ttjets =  (TH2D*)f->Get("Histos2D/nTagsvHTTTJets"); 

  h_Signal =  (TH2D*) h_TagsHT_ttjets->Clone();
  h_Control = (TH2D*) h_TagsHT_ttjets->Clone();

 
  //h_Signal->GetYaxis()->SetRange(500,1000);
  //h_Control->GetYaxis()->SetRange(1,500);
  
  TH1D *h1_sig = h_TagsHT_ttjets->ProjectionX("",0,17);
  h1_sig_clone =  (TH2D*) h1_sig->Clone();


  TH1D *h1_ctrl = h_TagsHT_ttjets->ProjectionX("",17,50);
  h1_ctrl_clone =  (TH2D*) h1_ctrl->Clone();


  //  TH1D* h1_sig =  h_Signal->ProjectionX();
  //TH1D* h1_ctrl =  h_Control->ProjectionX();

   h1_sig_clone->DrawNormalized();
   h1_ctrl_clone->DrawNormalized("same");

   h1_ctrl_clone->SetLineColor(kRed);


}
