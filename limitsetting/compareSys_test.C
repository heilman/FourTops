
{
 TFile* f1_Mu= new TFile("../SystematicShapes_Mu.root");
 TFile* f2_Mu = new TFile("../FourTop_EventSelection_wMETCut_Mu.root");

 TH1F * hSig_Mu = (TH1F*)f2_Mu->Get("MultiSamplePlot_MVA/MVA_NP_overlay_TTTT");
 TH1F * h_leptonSF_Down_tttt = (TH1F*)f1_Mu->Get("MVA_leptonSF_Down_tttt");
 TH1F * h_leptonSF_Up_tttt = (TH1F*)f1_Mu->Get("MVA_leptonSF_Up_tttt");

 cout<<" "<<endl;

 cout<<"tttt "<<endl;
 cout <<"leptonSF  "<< hSig_Mu->Integral(0,16) <<"   "<< h_leptonSF_Down_tttt->Integral(0,16) << "  " <<  h_leptonSF_Up_tttt->Integral(0,16)  <<endl;




}


