using namespace RooFit ;

void ex3b()
{
  TFile* f = TFile::Open("model.root") ;
  RooWorkspace* w = f->Get("w") ;
  w->exportToCint() ;

  RooPlot* frame = w::x.frame() ;
  w::data.plotOn(frame) ;
  w::model.plotOn(frame) ;
  w::model.plotOn(frame,Components("bkg"),LineStyle(kDashed)) ;
  frame->Draw() ;
}
