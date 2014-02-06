{

/////////////////////////
// F I L E S
//////////////////////////
  TFile* f1 = new TFile("SystematicShapes_norm_El.root");
  TFile* f2 = new TFile("NominalShapes_El.root");

////////////////////////////////////////
// C O M P A R E   S H A P E S
////////////////////////////////////////


      ///////////////////////////////////
      /// E L E C T R O N  C H A N N E L  
      ///////////////////////////////////

            ///////////////////////////
            /// I N C L U S I V E   
            ////////////////////////////

  TCanvas * c = new TCanvas();

  h_Nominal = (TH1F*)f2.Get("tt");
  h_Nominal->SetLineWidth(1);  
  h_Nominal->Draw("HIST");
  h_Nominal->SetTitle("");

  h_Scale_Down = (TH1F*)f1.Get("Scale_Down");
  h_Scale_Up =(TH1F*)f1.Get("Scale_Up");
  h_Scale_Down->SetMarkerColor(kGreen);
  h_Scale_Up->SetMarkerColor(kGreen);
  h_Scale_Down->SetLineColor(kGreen);
  h_Scale_Up->SetLineColor(kGreen);
    
  h_Scale_Down->SetMarkerStyle(23);
  h_Scale_Up->SetMarkerStyle(22);

  h_Scale_Down->Draw("psame");
  h_Scale_Up->Draw("psame");

  h_Matching_Down = (TH1F*)f1.Get("Matching_Down");
  h_Matching_Up =(TH1F*)f1.Get("Matching_Up");
  h_Matching_Down->SetMarkerColor(kViolet);
  h_Matching_Up->SetMarkerColor(kViolet);
  h_Matching_Down->SetLineColor(kViolet);
  h_Matching_Up->SetLineColor(kViolet);
    
  h_Matching_Down->Draw("psame");
  h_Matching_Up->Draw("psame");
  h_Matching_Down->SetMarkerStyle(23);
  h_Matching_Up->SetMarkerStyle(22);

  h_JES_Down = (TH1F*)f1.Get("JES_Down");
  h_JES_Up =(TH1F*)f1.Get("JES_Up");
  h_JES_Down->SetMarkerColor(kYellow);
  h_JES_Up->SetMarkerColor(kYellow);
  h_JES_Down->SetLineColor(kYellow);
  h_JES_Up->SetLineColor(kYellow);
    
  h_JES_Down->Draw("psame");
  h_JES_Up->Draw("psame");
  h_JES_Down->SetMarkerStyle(23);
  h_JES_Up->SetMarkerStyle(22);


  h_ttbb_Down = (TH1F*)f1.Get("ttbb_Down");
  h_ttbb_Up =(TH1F*)f1.Get("ttbb_Up");
  h_ttbb_Down->SetMarkerColor(kRed);
  h_ttbb_Up->SetMarkerColor(kRed);
  h_ttbb_Down->SetLineColor(kRed);
  h_ttbb_Up->SetLineColor(kRed);
  h_ttbb_Down->Draw("psame");
  h_ttbb_Up->Draw("psame");
  h_ttbb_Down->SetMarkerStyle(23);
  h_ttbb_Up->SetMarkerStyle(22);

  h_bTag_Down = (TH1F*)f1.Get("bTag_Down");
  h_bTag_Up =(TH1F*)f1.Get("bTag_Up");
  h_bTag_Down->SetMarkerColor(kBlue);
  h_bTag_Up->SetMarkerColor(kBlue);
  h_bTag_Down->SetLineColor(kBlue);
  h_bTag_Up->SetLineColor(kBlue);
  h_bTag_Down->Draw("psame");
  h_bTag_Up->Draw("psame");
  h_bTag_Down->SetMarkerStyle(23);
  h_bTag_Up->SetMarkerStyle(22);

  h_misTag_Down = (TH1F*)f1.Get("misTag_Down");
  h_misTag_Up =(TH1F*)f1.Get("misTag_Up");
  h_misTag_Down->SetMarkerColor(kCyan);
  h_misTag_Up->SetMarkerColor(kCyan);
  h_misTag_Down->SetMarkerColor(kCyan);
  h_misTag_Up->SetMarkerColor(kCyan);
    
  h_misTag_Down->Draw("psame");
  h_misTag_Up->Draw("psame");
  h_misTag_Down->SetMarkerStyle(23);
  h_misTag_Up->SetMarkerStyle(22);
    
    
    h_leptonSF_Down = (TH1F*)f1.Get("leptonSF_Down");
    h_leptonSF_Up =(TH1F*)f1.Get("leptonSF_Up");
    h_leptonSF_Down->SetMarkerColor(kOrange);
    h_leptonSF_Up->SetMarkerColor(kOrange);
    h_leptonSF_Down->SetMarkerColor(kOrange);
    h_leptonSF_Up->SetMarkerColor(kOrange);
    
    h_leptonSF_Down->Draw("psame");
    h_leptonSF_Up->Draw("psame");
    h_leptonSF_Down->SetMarkerStyle(23);
    h_leptonSF_Up->SetMarkerStyle(22);
    


    h_PU_Down = (TH1F*)f1.Get("PU_Down");
    h_PU_Up =(TH1F*)f1.Get("PU_Up");
    h_PU_Down->SetMarkerColor(kMagenta);
    h_PU_Up->SetMarkerColor(kMagenta);
    h_PU_Down->SetMarkerColor(kMagenta);
    h_PU_Up->SetMarkerColor(kMagenta);
    
    h_PU_Down->Draw("psame");
    h_PU_Up->Draw("psame");
    h_PU_Down->SetMarkerStyle(23);
    h_PU_Up->SetMarkerStyle(22);


    h_JER_Down = (TH1F*)f1.Get("JER_Down");
    h_JER_Up =(TH1F*)f1.Get("JER_Up");
    h_JER_Down->SetMarkerColor(kGray);
    h_JER_Up->SetMarkerColor(kGray);
    h_JER_Down->SetMarkerColor(kGray);
    h_JER_Up->SetMarkerColor(kGray);
    
    h_JER_Down->Draw("psame");
    h_JER_Up->Draw("psame");
    h_JER_Down->SetMarkerStyle(23);
    h_JER_Up->SetMarkerStyle(22);

    TLegend* leg = new TLegend(0.7, 0.25, .95, .95);
    leg->SetFillColor(0);
    leg->AddEntry(h_Scale_Down, "Scaling down", "p");
    leg->AddEntry(h_Scale_Up, "Scaling up", "p");
    leg->AddEntry(h_Matching_Down, "Matching down", "p");
    leg->AddEntry(h_Matching_Up, "Matching up", "p");
    leg->AddEntry(h_JES_Down, "JES down", "p");
    leg->AddEntry(h_JES_Up, "JES up", "p");
    leg->AddEntry(h_bTag_Down, "b-tag down", "p");
    leg->AddEntry(h_bTag_Up, "b-tag up", "p");
    leg->AddEntry(h_misTag_Down, "mis-tag down", "p");
    leg->AddEntry(h_misTag_Up, "mis-tag up", "p");
    leg->AddEntry(h_ttbb_Down, "ttbb down", "p");
    leg->AddEntry(h_ttbb_Up, "ttbb up", "p");
    leg->AddEntry(h_leptonSF_Down, "lepton SF down", "p");
    leg->AddEntry(h_leptonSF_Up, "lepton SF up", "p");
    leg->AddEntry(h_PU_Down, "Pile up down", "p");
    leg->AddEntry(h_PU_Up, "Pile up up", "p");
    leg->AddEntry(h_JER_Down, "JER  down", "p");
    leg->AddEntry(h_JER_Up, "JER up", "p");
    leg->AddEntry(h_Nominal, "Nominal shape", "l");
    leg->Draw();

    TLatex text;
    text.SetNDC(true);
    text.SetTextAlign(12);
    text.SetTextFont(42);
    text.SetTextSize(0.07);
    text.DrawLatex(0.18,0.72,"#font[42]{N_{jets}=6 }");

    cout <<"chi2 test scale : Down -> "<< h_Nominal->Chi2Test(h_Scale_Down, "WWCHI2/NDF") << ",  Up = " << h_Nominal->Chi2Test(h_Scale_Up, "WWCHI2/NDF")  <<endl;
    cout <<"chi2 test matching : Down -> "<< h_Nominal->Chi2Test(h_Matching_Down, "WWCHI2/NDF") << ",  Up = " << h_Nominal->Chi2Test(h_Matching_Up, "WWCHI2/NDF")  <<endl;
    cout <<"chi2 test JES : Down -> "<< h_Nominal->Chi2Test(h_JES_Down, "WWCHI2/NDF") << ",  Up = " << h_Nominal->Chi2Test(h_JES_Up, "WWCHI2/NDF")  <<endl;
    cout <<"chi2 test JER : Down -> "<< h_Nominal->Chi2Test(h_JER_Down, "WWCHI2/NDF") << ",  Up = " << h_Nominal->Chi2Test(h_JER_Up, "WWCHI2/NDF")  <<endl;
    cout <<"chi2 test bTag : Down -> "<< h_Nominal->Chi2Test(h_bTag_Down, "WWCHI2/NDF") << ",  Up = " << h_Nominal->Chi2Test(h_bTag_Up, "WWCHI2/NDF")  <<endl;
    cout <<"chi2 test misTag : Down -> "<< h_Nominal->Chi2Test(h_misTag_Down, "WWCHI2/NDF") << ",  Up = " << h_Nominal->Chi2Test(h_misTag_Up, "WWCHI2/NDF")  <<endl;
    cout <<"chi2 test ttbb : Down -> "<< h_Nominal->Chi2Test(h_ttbb_Down, "WWCHI2/NDF") << ",  Up = " << h_Nominal->Chi2Test(h_ttbb_Up, "WWCHI2/NDF")  <<endl; 
    cout <<"chi2 test leptonSF : Down -> "<< h_Nominal->Chi2Test(h_leptonSF_Down, "WWCHI2/NDF") << ",  Up = " << h_Nominal->Chi2Test(h_leptonSF_Up, "WWCHI2/NDF")  <<endl;
    cout <<"chi2 test PU : Down -> "<< h_Nominal->Chi2Test(h_PU_Down, "WWCHI2/NDF") << ",  Up = " << h_Nominal->Chi2Test(h_PU_Up, "WWCHI2/NDF")  <<endl;


    cout <<"Integrals inclusive "<<endl;
    cout <<"Scale down = "<< h_Scale_Down->Integral()  <<endl;
    cout <<"Scale up = "<< h_Scale_Up->Integral()  <<endl;
    cout <<"Matching down = "<< h_Matching_Down->Integral()  <<endl;
    cout <<"Matching up = "<< h_Matching_Up->Integral()  <<endl;
    cout <<"JES down = "<< h_JES_Down->Integral()  <<endl;
    cout <<"JES up = "<< h_JES_Up->Integral()  <<endl;
    cout <<"bTag down = "<< h_bTag_Down->Integral()  <<endl;
    cout <<"bTag up = "<< h_bTag_Up->Integral()  <<endl;
    cout <<"misTag down = "<< h_misTag_Down->Integral()  <<endl;
    cout <<"misTag up = "<< h_misTag_Up->Integral()  <<endl;
    cout <<"ttbb down = "<< h_ttbb_Down->Integral()  <<endl;
    cout <<"ttbb up = "<< h_ttbb_Up->Integral()  <<endl;
    cout <<"leptonSF down = "<< h_leptonSF_Down->Integral()  <<endl;
    cout <<"leptonSF up = "<< h_leptonSF_Up->Integral()  <<endl;
    cout <<"PU down = "<< h_PU_Down->Integral()  <<endl;
    cout <<"PU up = "<< h_PU_Up->Integral()  <<endl;
    cout <<"JER down = "<< h_JER_Down->Integral()  <<endl;
    cout <<"JER up = "<< h_JER_Up->Integral()  <<endl;
    cout <<"Nominal = "<< h_Nominal->Integral()  <<endl;

    cout <<" "<<endl;


































              ///////////////////////////
              /// S I X  J E T  B I N 
              ////////////////////////////

  TCanvas * c_6j = new TCanvas();

  h_Nominal_6j = (TH1F*)f2.Get("tt_6j");
  h_Nominal_6j->SetLineWidth(1);  
  h_Nominal_6j->Draw("HIST");
  h_Nominal_6j->SetTitle("");

  h_Scale_Down_6j = (TH1F*)f1.Get("Scale_Down_6j");
  h_Scale_Up_6j =(TH1F*)f1.Get("Scale_Up_6j");
  h_Scale_Down_6j->SetMarkerColor(kGreen);
  h_Scale_Up_6j->SetMarkerColor(kGreen);
  h_Scale_Down_6j->SetLineColor(kGreen);
  h_Scale_Up_6j->SetLineColor(kGreen);
    
  h_Scale_Down_6j->SetMarkerStyle(23);
  h_Scale_Up_6j->SetMarkerStyle(22);

  h_Scale_Down_6j->Draw("psame");
  h_Scale_Up_6j->Draw("psame");

  h_Matching_Down_6j = (TH1F*)f1.Get("Matching_Down_6j");
  h_Matching_Up_6j =(TH1F*)f1.Get("Matching_Up_6j");
  h_Matching_Down_6j->SetMarkerColor(kViolet);
  h_Matching_Up_6j->SetMarkerColor(kViolet);
  h_Matching_Down_6j->SetLineColor(kViolet);
  h_Matching_Up_6j->SetLineColor(kViolet);
    
  h_Matching_Down_6j->Draw("psame");
  h_Matching_Up_6j->Draw("psame");
  h_Matching_Down_6j->SetMarkerStyle(23);
  h_Matching_Up_6j->SetMarkerStyle(22);

  h_JES_Down_6j = (TH1F*)f1.Get("JES_Down_6j");
  h_JES_Up_6j =(TH1F*)f1.Get("JES_Up_6j");
  h_JES_Down_6j->SetMarkerColor(kYellow);
  h_JES_Up_6j->SetMarkerColor(kYellow);
  h_JES_Down_6j->SetLineColor(kYellow);
  h_JES_Up_6j->SetLineColor(kYellow);
    
  h_JES_Down_6j->Draw("psame");
  h_JES_Up_6j->Draw("psame");
  h_JES_Down_6j->SetMarkerStyle(23);
  h_JES_Up_6j->SetMarkerStyle(22);


  h_ttbb_Down_6j = (TH1F*)f1.Get("ttbb_Down_6j");
  h_ttbb_Up_6j =(TH1F*)f1.Get("ttbb_Up_6j");
  h_ttbb_Down_6j->SetMarkerColor(kRed);
  h_ttbb_Up_6j->SetMarkerColor(kRed);
  h_ttbb_Down_6j->SetLineColor(kRed);
  h_ttbb_Up_6j->SetLineColor(kRed);
  h_ttbb_Down_6j->Draw("psame");
  h_ttbb_Up_6j->Draw("psame");
  h_ttbb_Down_6j->SetMarkerStyle(23);
  h_ttbb_Up_6j->SetMarkerStyle(22);

  h_bTag_Down_6j = (TH1F*)f1.Get("bTag_Down_6j");
  h_bTag_Up_6j =(TH1F*)f1.Get("bTag_Up_6j");
  h_bTag_Down_6j->SetMarkerColor(kBlue);
  h_bTag_Up_6j->SetMarkerColor(kBlue);
  h_bTag_Down_6j->SetLineColor(kBlue);
  h_bTag_Up_6j->SetLineColor(kBlue);
  h_bTag_Down_6j->Draw("psame");
  h_bTag_Up_6j->Draw("psame");
  h_bTag_Down_6j->SetMarkerStyle(23);
  h_bTag_Up_6j->SetMarkerStyle(22);

  h_misTag_Down_6j = (TH1F*)f1.Get("misTag_Down_6j");
  h_misTag_Up_6j =(TH1F*)f1.Get("misTag_Up_6j");
  h_misTag_Down_6j->SetMarkerColor(kCyan);
  h_misTag_Up_6j->SetMarkerColor(kCyan);
  h_misTag_Down_6j->SetMarkerColor(kCyan);
  h_misTag_Up_6j->SetMarkerColor(kCyan);
    
  h_misTag_Down_6j->Draw("psame");
  h_misTag_Up_6j->Draw("psame");
  h_misTag_Down_6j->SetMarkerStyle(23);
  h_misTag_Up_6j->SetMarkerStyle(22);
    
    
    h_leptonSF_Down_6j = (TH1F*)f1.Get("leptonSF_Down_6j");
    h_leptonSF_Up_6j =(TH1F*)f1.Get("leptonSF_Up_6j");
    h_leptonSF_Down_6j->SetMarkerColor(kOrange);
    h_leptonSF_Up_6j->SetMarkerColor(kOrange);
    h_leptonSF_Down_6j->SetMarkerColor(kOrange);
    h_leptonSF_Up_6j->SetMarkerColor(kOrange);
    
    h_leptonSF_Down_6j->Draw("psame");
    h_leptonSF_Up_6j->Draw("psame");
    h_leptonSF_Down_6j->SetMarkerStyle(23);
    h_leptonSF_Up_6j->SetMarkerStyle(22);
    


    h_PU_Down_6j = (TH1F*)f1.Get("PU_Down_6j");
    h_PU_Up_6j =(TH1F*)f1.Get("PU_Up_6j");
    h_PU_Down_6j->SetMarkerColor(kMagenta);
    h_PU_Up_6j->SetMarkerColor(kMagenta);
    h_PU_Down_6j->SetMarkerColor(kMagenta);
    h_PU_Up_6j->SetMarkerColor(kMagenta);
    
    h_PU_Down_6j->Draw("psame");
    h_PU_Up_6j->Draw("psame");
    h_PU_Down_6j->SetMarkerStyle(23);
    h_PU_Up_6j->SetMarkerStyle(22);


    h_JER_Down_6j = (TH1F*)f1.Get("JER_Down_6j");
    h_JER_Up_6j =(TH1F*)f1.Get("JER_Up_6j");
    h_JER_Down_6j->SetMarkerColor(kGray);
    h_JER_Up_6j->SetMarkerColor(kGray);
    h_JER_Down_6j->SetMarkerColor(kGray);
    h_JER_Up_6j->SetMarkerColor(kGray);
    
    h_JER_Down_6j->Draw("psame");
    h_JER_Up_6j->Draw("psame");
    h_JER_Down_6j->SetMarkerStyle(23);
    h_JER_Up_6j->SetMarkerStyle(22);

    TLegend* leg_6j = new TLegend(0.7, 0.25, .95, .95);
    leg_6j->SetFillColor(0);
    leg_6j->AddEntry(h_Scale_Down_6j, "Scaling down", "p");
    leg_6j->AddEntry(h_Scale_Up_6j, "Scaling up", "p");
    leg_6j->AddEntry(h_Matching_Down_6j, "Matching down", "p");
    leg_6j->AddEntry(h_Matching_Up_6j, "Matching up", "p");
    leg_6j->AddEntry(h_JES_Down_6j, "JES down", "p");
    leg_6j->AddEntry(h_JES_Up_6j, "JES up", "p");
    leg_6j->AddEntry(h_bTag_Down_6j, "b-tag down", "p");
    leg_6j->AddEntry(h_bTag_Up_6j, "b-tag up", "p");
    leg_6j->AddEntry(h_misTag_Down_6j, "mis-tag down", "p");
    leg_6j->AddEntry(h_misTag_Up_6j, "mis-tag up", "p");
    leg_6j->AddEntry(h_ttbb_Down_6j, "ttbb down", "p");
    leg_6j->AddEntry(h_ttbb_Up_6j, "ttbb up", "p");
    leg_6j->AddEntry(h_leptonSF_Down_6j, "lepton SF down", "p");
    leg_6j->AddEntry(h_leptonSF_Up_6j, "lepton SF up", "p");
    leg_6j->AddEntry(h_PU_Down_6j, "Pile up down", "p");
    leg_6j->AddEntry(h_PU_Up_6j, "Pile up up", "p");
    leg_6j->AddEntry(h_JER_Down_6j, "JER  down", "p");
    leg_6j->AddEntry(h_JER_Up_6j, "JER up", "p");
    leg_6j->AddEntry(h_Nominal_6j, "Nominal shape", "l");
    leg_6j->Draw();

    TLatex text_6j;
    text_6j.SetNDC(true);
    text_6j.SetTextAlign(12);
    text_6j.SetTextFont(42);
    text_6j.SetTextSize(0.07);
    text_6j.DrawLatex(0.18,0.72,"#font[42]{N_{jets}=6 }");


    cout <<"Integrals 6j "<<endl;
    cout <<"Scale down = "<< h_Scale_Down_6j->Integral()  <<endl;
    cout <<"Scale up = "<< h_Scale_Up_6j->Integral()  <<endl;
    cout <<"Matching down = "<< h_Matching_Down_6j->Integral()  <<endl;
    cout <<"Matching up = "<< h_Matching_Up_6j->Integral()  <<endl;
    cout <<"JES down = "<< h_JES_Down_6j->Integral()  <<endl;
    cout <<"JES up = "<< h_JES_Up_6j->Integral()  <<endl;
    cout <<"bTag down = "<< h_bTag_Down_6j->Integral()  <<endl;
    cout <<"bTag up = "<< h_bTag_Up_6j->Integral()  <<endl;
    cout <<"misTag down = "<< h_misTag_Down_6j->Integral()  <<endl;
    cout <<"misTag up = "<< h_misTag_Up_6j->Integral()  <<endl;
    cout <<"ttbb down = "<< h_ttbb_Down_6j->Integral()  <<endl;
    cout <<"ttbb up = "<< h_ttbb_Up_6j->Integral()  <<endl;
    cout <<"leptonSF down = "<< h_leptonSF_Down_6j->Integral()  <<endl;
    cout <<"leptonSF up = "<< h_leptonSF_Up_6j->Integral()  <<endl;
    cout <<"PU down = "<< h_PU_Down_6j->Integral()  <<endl;
    cout <<"PU up = "<< h_PU_Up_6j->Integral()  <<endl;
    cout <<"JER down = "<< h_JER_Down_6j->Integral()  <<endl;
    cout <<"JER up = "<< h_JER_Up_6j->Integral()  <<endl;
    cout <<"Nominal = "<< h_Nominal_6j->Integral()  <<endl;

    cout <<" "<<endl;



    //////////////////////////////
    ///// S E V E N  J E T  B I N
    //////////////////////////////

  TCanvas * c_7j = new TCanvas();

  h_Nominal_7j = (TH1F*)f2.Get("tt_7j");
  h_Nominal_7j->SetLineWidth(1);  
  h_Nominal_7j->Draw("HIST");
  h_Nominal_7j->SetTitle("");

  h_Scale_Down_7j = (TH1F*)f1.Get("Scale_Down_7j");
  h_Scale_Up_7j =(TH1F*)f1.Get("Scale_Up_7j");
  h_Scale_Down_7j->SetMarkerColor(kGreen);
  h_Scale_Up_7j->SetMarkerColor(kGreen);
  h_Scale_Down_7j->SetLineColor(kGreen);
  h_Scale_Up_7j->SetLineColor(kGreen);
    
  h_Scale_Down_7j->SetMarkerStyle(23);
  h_Scale_Up_7j->SetMarkerStyle(22);

  h_Scale_Down_7j->Draw("psame");
  h_Scale_Up_7j->Draw("psame");

  h_Matching_Down_7j = (TH1F*)f1.Get("Matching_Down_7j");
  h_Matching_Up_7j =(TH1F*)f1.Get("Matching_Up_7j");
  h_Matching_Down_7j->SetMarkerColor(kViolet);
  h_Matching_Up_7j->SetMarkerColor(kViolet);
  h_Matching_Down_7j->SetLineColor(kViolet);
  h_Matching_Up_7j->SetLineColor(kViolet);
    
  h_Matching_Down_7j->Draw("psame");
  h_Matching_Up_7j->Draw("psame");
  h_Matching_Down_7j->SetMarkerStyle(23);
  h_Matching_Up_7j->SetMarkerStyle(22);

  h_JES_Down_7j = (TH1F*)f1.Get("JES_Down_7j");
  h_JES_Up_7j =(TH1F*)f1.Get("JES_Up_7j");
  h_JES_Down_7j->SetMarkerColor(kYellow);
  h_JES_Up_7j->SetMarkerColor(kYellow);
  h_JES_Down_7j->SetLineColor(kYellow);
  h_JES_Up_7j->SetLineColor(kYellow);
    
  h_JES_Down_7j->Draw("psame");
  h_JES_Up_7j->Draw("psame");
  h_JES_Down_7j->SetMarkerStyle(23);
  h_JES_Up_7j->SetMarkerStyle(22);


  h_ttbb_Down_7j = (TH1F*)f1.Get("ttbb_Down_7j");
  h_ttbb_Up_7j =(TH1F*)f1.Get("ttbb_Up_7j");
  h_ttbb_Down_7j->SetMarkerColor(kRed);
  h_ttbb_Up_7j->SetMarkerColor(kRed);
  h_ttbb_Down_7j->SetLineColor(kRed);
  h_ttbb_Up_7j->SetLineColor(kRed);
  h_ttbb_Down_7j->Draw("psame");
  h_ttbb_Up_7j->Draw("psame");
  h_ttbb_Down_7j->SetMarkerStyle(23);
  h_ttbb_Up_7j->SetMarkerStyle(22);

  h_bTag_Down_7j = (TH1F*)f1.Get("bTag_Down_7j");
  h_bTag_Up_7j =(TH1F*)f1.Get("bTag_Up_7j");
  h_bTag_Down_7j->SetMarkerColor(kBlue);
  h_bTag_Up_7j->SetMarkerColor(kBlue);
  h_bTag_Down_7j->SetLineColor(kBlue);
  h_bTag_Up_7j->SetLineColor(kBlue);
  h_bTag_Down_7j->Draw("psame");
  h_bTag_Up_7j->Draw("psame");
  h_bTag_Down_7j->SetMarkerStyle(23);
  h_bTag_Up_7j->SetMarkerStyle(22);

  h_misTag_Down_7j = (TH1F*)f1.Get("misTag_Down_7j");
  h_misTag_Up_7j =(TH1F*)f1.Get("misTag_Up_7j");
  h_misTag_Down_7j->SetMarkerColor(kCyan);
  h_misTag_Up_7j->SetMarkerColor(kCyan);
  h_misTag_Down_7j->SetMarkerColor(kCyan);
  h_misTag_Up_7j->SetMarkerColor(kCyan);
    
  h_misTag_Down_7j->Draw("psame");
  h_misTag_Up_7j->Draw("psame");
  h_misTag_Down_7j->SetMarkerStyle(23);
  h_misTag_Up_7j->SetMarkerStyle(22);
    
    
    h_leptonSF_Down_7j = (TH1F*)f1.Get("leptonSF_Down_7j");
    h_leptonSF_Up_7j =(TH1F*)f1.Get("leptonSF_Up_7j");
    h_leptonSF_Down_7j->SetMarkerColor(kOrange);
    h_leptonSF_Up_7j->SetMarkerColor(kOrange);
    h_leptonSF_Down_7j->SetMarkerColor(kOrange);
    h_leptonSF_Up_7j->SetMarkerColor(kOrange);
    
    h_leptonSF_Down_7j->Draw("psame");
    h_leptonSF_Up_7j->Draw("psame");
    h_leptonSF_Down_7j->SetMarkerStyle(23);
    h_leptonSF_Up_7j->SetMarkerStyle(22);
    


    h_PU_Down_7j = (TH1F*)f1.Get("PU_Down_7j");
    h_PU_Up_7j =(TH1F*)f1.Get("PU_Up_7j");
    h_PU_Down_7j->SetMarkerColor(kMagenta);
    h_PU_Up_7j->SetMarkerColor(kMagenta);
    h_PU_Down_7j->SetMarkerColor(kMagenta);
    h_PU_Up_7j->SetMarkerColor(kMagenta);
    
    h_PU_Down_7j->Draw("psame");
    h_PU_Up_7j->Draw("psame");
    h_PU_Down_7j->SetMarkerStyle(23);
    h_PU_Up_7j->SetMarkerStyle(22);


    h_JER_Down_7j = (TH1F*)f1.Get("JER_Down_7j");
    h_JER_Up_7j =(TH1F*)f1.Get("JER_Up_7j");
    h_JER_Down_7j->SetMarkerColor(kGray);
    h_JER_Up_7j->SetMarkerColor(kGray);
    h_JER_Down_7j->SetMarkerColor(kGray);
    h_JER_Up_7j->SetMarkerColor(kGray);
    
    h_JER_Down_7j->Draw("psame");
    h_JER_Up_7j->Draw("psame");
    h_JER_Down_7j->SetMarkerStyle(23);
    h_JER_Up_7j->SetMarkerStyle(22);

    TLegend* leg_7j = new TLegend(0.7, 0.25, 0.95, 0.95);
    leg_7j->SetFillColor(0);
    leg_7j->AddEntry(h_Scale_Down_7j, "Scaling down", "p");
    leg_7j->AddEntry(h_Scale_Up_7j, "Scaling up", "p");
    leg_7j->AddEntry(h_Matching_Down_7j, "Matching down", "p");
    leg_7j->AddEntry(h_Matching_Up_7j, "Matching up", "p");
    leg_7j->AddEntry(h_JES_Down_7j, "JES down", "p");
    leg_7j->AddEntry(h_JES_Up_7j, "JES up", "p");
    leg_7j->AddEntry(h_bTag_Down_7j, "b-tag down", "p");
    leg_7j->AddEntry(h_bTag_Up_7j, "b-tag up", "p");
    leg_7j->AddEntry(h_misTag_Down_7j, "mis-tag down", "p");
    leg_7j->AddEntry(h_misTag_Up_7j, "mis-tag up", "p");
    leg_7j->AddEntry(h_ttbb_Down_7j, "ttbb down", "p");
    leg_7j->AddEntry(h_ttbb_Up_7j, "ttbb up", "p");
    leg_7j->AddEntry(h_leptonSF_Down_7j, "lepton SF down", "p");
    leg_7j->AddEntry(h_leptonSF_Up_7j, "lepton SF up", "p");
    leg_7j->AddEntry(h_PU_Down_7j, "Pile up down", "p");
    leg_7j->AddEntry(h_PU_Up_7j, "Pile up up", "p");
    leg_7j->AddEntry(h_JER_Down_7j, "JER  down", "p");
    leg_7j->AddEntry(h_JER_Up_7j, "JER up", "p");
    leg_7j->AddEntry(h_Nominal_7j, "Nominal shape", "l");
    leg_7j->Draw();



    TLatex text_7j;
    text_7j.SetNDC(true);
    text_7j.SetTextAlign(12);
    text_7j.SetTextFont(42);
    text_7j.SetTextSize(0.07);
    text_7j.DrawLatex(0.18,0.72,"#font[42]{N_{jets}=7 }");



  //////////////////////////////
  ///// E I G H T    J E T  B I N
  //////////////////////////////

  TCanvas * c_8j = new TCanvas();

  h_Nominal_8j = (TH1F*)f2.Get("tt_8j");
  h_Nominal_8j->SetLineWidth(1);  
  h_Nominal_8j->Draw("HIST");
  h_Nominal_8j->SetTitle("");

  h_Scale_Down_8j = (TH1F*)f1.Get("Scale_Down_8j");
  h_Scale_Up_8j =(TH1F*)f1.Get("Scale_Up_8j");
  h_Scale_Down_8j->SetMarkerColor(kGreen);
  h_Scale_Up_8j->SetMarkerColor(kGreen);
  h_Scale_Down_8j->SetLineColor(kGreen);
  h_Scale_Up_8j->SetLineColor(kGreen);
    
  h_Scale_Down_8j->SetMarkerStyle(23);
  h_Scale_Up_8j->SetMarkerStyle(22);

  h_Scale_Down_8j->Draw("psame");
  h_Scale_Up_8j->Draw("psame");

  h_Matching_Down_8j = (TH1F*)f1.Get("Matching_Down_8j");
  h_Matching_Up_8j =(TH1F*)f1.Get("Matching_Up_8j");
  h_Matching_Down_8j->SetMarkerColor(kViolet);
  h_Matching_Up_8j->SetMarkerColor(kViolet);
  h_Matching_Down_8j->SetLineColor(kViolet);
  h_Matching_Up_8j->SetLineColor(kViolet);
    
  h_Matching_Down_8j->Draw("psame");
  h_Matching_Up_8j->Draw("psame");
  h_Matching_Down_8j->SetMarkerStyle(23);
  h_Matching_Up_8j->SetMarkerStyle(22);

  h_JES_Down_8j = (TH1F*)f1.Get("JES_Down_8j");
  h_JES_Up_8j =(TH1F*)f1.Get("JES_Up_8j");
  h_JES_Down_8j->SetMarkerColor(kYellow);
  h_JES_Up_8j->SetMarkerColor(kYellow);
  h_JES_Down_8j->SetLineColor(kYellow);
  h_JES_Up_8j->SetLineColor(kYellow);
    
  h_JES_Down_8j->Draw("psame");
  h_JES_Up_8j->Draw("psame");
  h_JES_Down_8j->SetMarkerStyle(23);
  h_JES_Up_8j->SetMarkerStyle(22);


  h_ttbb_Down_8j = (TH1F*)f1.Get("ttbb_Down_8j");
  h_ttbb_Up_8j =(TH1F*)f1.Get("ttbb_Up_8j");
  h_ttbb_Down_8j->SetMarkerColor(kRed);
  h_ttbb_Up_8j->SetMarkerColor(kRed);
  h_ttbb_Down_8j->SetLineColor(kRed);
  h_ttbb_Up_8j->SetLineColor(kRed);
  h_ttbb_Down_8j->Draw("psame");
  h_ttbb_Up_8j->Draw("psame");
  h_ttbb_Down_8j->SetMarkerStyle(23);
  h_ttbb_Up_8j->SetMarkerStyle(22);

  h_bTag_Down_8j = (TH1F*)f1.Get("bTag_Down_8j");
  h_bTag_Up_8j =(TH1F*)f1.Get("bTag_Up_8j");
  h_bTag_Down_8j->SetMarkerColor(kBlue);
  h_bTag_Up_8j->SetMarkerColor(kBlue);
  h_bTag_Down_8j->SetLineColor(kBlue);
  h_bTag_Up_8j->SetLineColor(kBlue);
  h_bTag_Down_8j->Draw("psame");
  h_bTag_Up_8j->Draw("psame");
  h_bTag_Down_8j->SetMarkerStyle(23);
  h_bTag_Up_8j->SetMarkerStyle(22);

  h_misTag_Down_8j = (TH1F*)f1.Get("misTag_Down_8j");
  h_misTag_Up_8j =(TH1F*)f1.Get("misTag_Up_8j");
  h_misTag_Down_8j->SetMarkerColor(kCyan);
  h_misTag_Up_8j->SetMarkerColor(kCyan);
  h_misTag_Down_8j->SetMarkerColor(kCyan);
  h_misTag_Up_8j->SetMarkerColor(kCyan);
    
  h_misTag_Down_8j->Draw("psame");
  h_misTag_Up_8j->Draw("psame");
  h_misTag_Down_8j->SetMarkerStyle(23);
  h_misTag_Up_8j->SetMarkerStyle(22);
    
    
    h_leptonSF_Down_8j = (TH1F*)f1.Get("leptonSF_Down_8j");
    h_leptonSF_Up_8j =(TH1F*)f1.Get("leptonSF_Up_8j");
    h_leptonSF_Down_8j->SetMarkerColor(kOrange);
    h_leptonSF_Up_8j->SetMarkerColor(kOrange);
    h_leptonSF_Down_8j->SetMarkerColor(kOrange);
    h_leptonSF_Up_8j->SetMarkerColor(kOrange);
    
    h_leptonSF_Down_8j->Draw("psame");
    h_leptonSF_Up_8j->Draw("psame");
    h_leptonSF_Down_8j->SetMarkerStyle(23);
    h_leptonSF_Up_8j->SetMarkerStyle(22);
    


    h_PU_Down_8j = (TH1F*)f1.Get("PU_Down_8j");
    h_PU_Up_8j =(TH1F*)f1.Get("PU_Up_8j");
    h_PU_Down_8j->SetMarkerColor(kMagenta);
    h_PU_Up_8j->SetMarkerColor(kMagenta);
    h_PU_Down_8j->SetMarkerColor(kMagenta);
    h_PU_Up_8j->SetMarkerColor(kMagenta);
    
    h_PU_Down_8j->Draw("psame");
    h_PU_Up_8j->Draw("psame");
    h_PU_Down_8j->SetMarkerStyle(23);
    h_PU_Up_8j->SetMarkerStyle(22);


    h_JER_Down_8j = (TH1F*)f1.Get("JER_Down_8j");
    h_JER_Up_8j =(TH1F*)f1.Get("JER_Up_8j");
    h_JER_Down_8j->SetMarkerColor(kGray);
    h_JER_Up_8j->SetMarkerColor(kGray);
    h_JER_Down_8j->SetMarkerColor(kGray);
    h_JER_Up_8j->SetMarkerColor(kGray);
    
    h_JER_Down_8j->Draw("psame");
    h_JER_Up_8j->Draw("psame");
    h_JER_Down_8j->SetMarkerStyle(23);
    h_JER_Up_8j->SetMarkerStyle(22);

    TLegend* leg_8j = new TLegend(0.7, 0.25, 0.95, 0.95);
    leg_8j->SetFillColor(0);
    leg_8j->AddEntry(h_Scale_Down_8j, "Scaling down", "p");
    leg_8j->AddEntry(h_Scale_Up_8j, "Scaling up", "p");
    leg_8j->AddEntry(h_Matching_Down_8j, "Matching down", "p");
    leg_8j->AddEntry(h_Matching_Up_8j, "Matching up", "p");
    leg_8j->AddEntry(h_JES_Down_8j, "JES down", "p");
    leg_8j->AddEntry(h_JES_Up_8j, "JES up", "p");
    leg_8j->AddEntry(h_bTag_Down_8j, "b-tag down", "p");
    leg_8j->AddEntry(h_bTag_Up_8j, "b-tag up", "p");
    leg_8j->AddEntry(h_misTag_Down_8j, "mis-tag down", "p");
    leg_8j->AddEntry(h_misTag_Up_8j, "mis-tag up", "p");
    leg_8j->AddEntry(h_ttbb_Down_8j, "ttbb down", "p");
    leg_8j->AddEntry(h_ttbb_Up_8j, "ttbb up", "p");
    leg_8j->AddEntry(h_leptonSF_Down_8j, "lepton SF down", "p");
    leg_8j->AddEntry(h_leptonSF_Up_8j, "lepton SF up", "p");
    leg_8j->AddEntry(h_PU_Down_8j, "Pile up down", "p");
    leg_8j->AddEntry(h_PU_Up_8j, "Pile up up", "p");
    leg_8j->AddEntry(h_JER_Down_8j, "JER  down", "p");
    leg_8j->AddEntry(h_JER_Up_8j, "JER up", "p");
    leg_8j->AddEntry(h_Nominal_8j, "Nominal shape", "l");
    leg_8j->Draw();

    TLatex text_8j;
    text_8j.SetNDC(true);
    text_8j.SetTextAlign(12);
    text_8j.SetTextFont(42);
    text_8j.SetTextSize(0.07);
    text_8j.DrawLatex(0.18,0.72,"#font[42]{N_{jets}>=8 }");




}
