{

    gStyle->SetOptStat(0000);
    TFile *f = new TFile("FourTop_EventSelection_Mu.root");
    
    h0  = (TH1F*)f.Get("Histos1D/BDT_Comb_type0");
    h1A = (TH1F*)f.Get("Histos1D/BDT_Comb_type1A");
    h1B = (TH1F*)f.Get("Histos1D/BDT_Comb_type1B");
    h2A = (TH1F*)f.Get("Histos1D/BDT_Comb_type2A");
    h2B = (TH1F*)f.Get("Histos1D/BDT_Comb_type2B");
    h3  = (TH1F*)f.Get("Histos1D/BDT_Comb_type3");
    hAll = (TH1F*)f.Get("Histos1D/BDT_Comb_typeAll");

    h0->SetLineColor(kYellow);
    h1A->SetLineColor(kGreen-5);
    h1B->SetLineColor(kGreen);
    h2A->SetLineColor(kBlue-5);
    h2B->SetLineColor(kBlue);
    h3->SetLineColor(kRed);


    h0->SetLineWidth(3);
    h1A->SetLineWidth(3);
    h1B->SetLineWidth(3);
    h2A->SetLineWidth(3);
    h2B->SetLineWidth(3);
    h3->SetLineWidth(3);



    //    h0->DrawNormalized();
    h1A->DrawNormalized();
    h1B->DrawNormalized("same");
    h2A->DrawNormalized("same");
    h2B->DrawNormalized("same");
    h3->DrawNormalized("same");

    leg = new TLegend(0.1,0.7,0.48,0.9);
    leg->SetHeader("Breakdown of combinations");
    leg->AddEntry(h0,"good combinations","f");
    leg->AddEntry(h1A,"b-jet not in combination","f");
    leg->AddEntry(h1B,"l-jet not in combination","f");
    leg->AddEntry(h2B,"b-jet and l-jet not in combination","f");
    leg->AddEntry(h2A,"both l-jets not in combination","f");
    leg->AddEntry(h3,"all jets not in combination","f");
    leg->AddEntry(hAll,"all combinations","p");
    leg->Draw();


    TCanvas * c2 = new TCanvas();


    h0->SetFillColor(kYellow);
    h1A->SetFillColor(kGreen-5);
    h1B->SetFillColor(kGreen);
    h2A->SetFillColor(kBlue-5);
    h2B->SetFillColor(kBlue);
    h3->SetFillColor(kRed);
    hAll->SetMarkerSize(0.9);


    THStack * Ts = new THStack();

     Ts->Add(h0);
     Ts->Add(h1A);
     Ts->Add(h1B);
     Ts->Add(h2A);
     Ts->Add(h2B);
     Ts->Add(h3);

    Ts->Draw();
    hAll->Draw("samep");
    c1->SetLogy();


    leg->Draw();






}
