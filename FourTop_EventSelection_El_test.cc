//////////////////////////////////////////////////////////////////////////////////
////         Analysis code for search for Four Top Production.                  ///
////////////////////////////////////////////////////////////////////////////////////

#define _USE_MATH_DEFINES
#include "TStyle.h"
#include "TPaveText.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "TRandom3.h"
#include "TRandom.h"
#include <iostream>
#include <map>
#include <cstdlib> 

//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
#include "TopTreeAnalysisBase/Selection/interface/FourTopSelectionTable.h"

#include "TopTreeAnalysisBase/Content/interface/AnalysisEnvironment.h"
#include "TopTreeAnalysisBase/Content/interface/Dataset.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"
#include "TopTreeAnalysisBase/Tools/interface/PlottingTools.h"
#include "TopTreeAnalysisBase/Tools/interface/MultiSamplePlot.h"
#include "TopTreeAnalysisBase/Tools/interface/TTreeLoader.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"

#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MakeBinning.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MEzCalculator.h"

#include "TopTreeAnalysisBase/Reconstruction/interface/TTreeObservables.h"

//This header file is taken directly from the BTV wiki. It contains
// to correctly apply an event level Btag SF. It is not yet on CVS
// as I hope to merge the functionality into BTagWeigtTools.h
#include "TopTreeAnalysisBase/Tools/interface/BTagSFUtil.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"

#include "TopTreeAnalysisBase/Tools/interface/JetCombiner.h"
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"

#include "TopTreeAnalysis/macros/Style.C"

using namespace std;
using namespace TopTree;
using namespace reweight;

bool split_ttbar = true;

pair<float, vector<unsigned int> > MVAvals1;
pair<float, vector<unsigned int> > MVAvals2;
pair<float, vector<unsigned int> > MVAvals2ndPass;
pair<float, vector<unsigned int> > MVAvals3rdPass;

int nMVASuccesses=0;
int nMatchedEvents=0;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

/// MultiPadPlot
map<string,MultiSamplePlot*> MultiPadPlot;

struct HighestTCHEBtag{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const{
    	return j1->btag_trackCountingHighEffBJetTags() > j2->btag_trackCountingHighEffBJetTags();
    }
};
struct HighestCVSBtag{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const{
    	return j1->btag_combinedSecondaryVertexBJetTags() > j2->btag_combinedSecondaryVertexBJetTags();
    }
};

bool match;

//To cout the Px, Py, Pz, E and Pt of objects
void coutObjectsFourVector(vector < TRootMuon* > init_muons, vector < TRootElectron* > init_electrons, vector < TRootJet* > init_jets, vector < TRootMET* > mets, string Comment);

int Factorial(int N);


int main (int argc, char *argv[])
{

  ofstream eventlist;
  eventlist.open ("interesting_events_mu.txt");


    TRandom3* rand = new TRandom3();

    //    BTagWeightTools * bTool = new BTagWeightTools("SFb-pt_payload_Moriond13.txt", "CSVM") ;
    
 BTagWeightTools * bTool = new BTagWeightTools("SFb-pt_NOttbar_payload_EPS13.txt", "CSVM") ;


  int doJESShift = 0; // 0: off 1: minus 2: plus
  cout << "doJESShift: " << doJESShift << endl;

  int doJERShift = 0; // 0: off (except nominal scalefactor for jer) 1: minus 2: plus
  cout << "doJERShift: " << doJERShift << endl;

  int dobTagEffShift = 0; //0: off (except nominal scalefactor for btag eff) 1: minus 2: plus
  cout << "dobTagEffShift: " << dobTagEffShift << endl;

  int domisTagEffShift = 0; //0: off (except nominal scalefactor for mistag eff) 1: minus 2: plus
  cout << "domisTagEffShift: " << domisTagEffShift << endl;

  int doPUShift = 0; //0: off (except nominal PU reweighting) 1: minus 2: plus
  cout << "doPUShift: " << doPUShift << endl;
    
  int doScaleShift = 0;
  cout << "doScaleShift: " << doScaleShift << endl;

 int doMatchingShift = 0;  //0: off (except nominal PU reweighting) 1: minus 2: plus
  cout << "doMatchingShift: " << doMatchingShift << endl;

 int dottbbShift = 0;  //0: off (except nominal PU reweighting) 1: minus 2: plus
  cout << "dottbbShift: " << dottbbShift << endl;

 int doLeptonSFShift = 0;
 cout << "doLeptonSFShift: " << doLeptonSFShift << endl;

 string leptonsyst = "Nominal";


	if (doLeptonSFShift==1){

	  leptonsyst = "Minus";
	}else if(doLeptonSFShift==2){
	  leptonsyst = "Plus";
}


  string btagger = "CSVM";
// b-tag scalefactor => TCHEL: data/MC scalefactor = 0.95 +- 0.10,    TCHEM: data/MC scalefactor = 0.94 +- 0.09
// mistag scalefactor => TCHEL: data/MC scalefactor = 1.11 +- 0.12,    TCHEM: data/MC scalefactor = 1.21 +- 0.17
  float scalefactorbtageff, mistagfactor;
  if(btagger == "TCHPM"  || btagger == "TCHET"  ||  btagger == "SSV" ){
    cout<<"This tagger ("<< btagger <<")is not commisioned in 2012, please use CSV, TCHP or JetProb"<<endl;
    exit(1);
}
  else if(btagger == "TCHEM") //redundant for now, but will use as skeleton for CSVM
  {
  	  if(dobTagEffShift == 0)
		scalefactorbtageff = 0.94;
	  if(dobTagEffShift == -1)
		scalefactorbtageff = 0.85;
	  if(dobTagEffShift == 1)
		scalefactorbtageff = 1.03;
		
	  if(domisTagEffShift == 0)
		mistagfactor = 1.21;
	  if(domisTagEffShift == 1)
		mistagfactor = 1.04;
	  if(domisTagEffShift == 2)
		mistagfactor = 1.38;
  }
  float workingpointvalue = 9999; //working points updated to 2012 BTV-POG recommendations.
 
  if(btagger == "TCHPM"  || btagger == "TCHET"  ||  btagger == "SSV" ){
    cout<<"This tagger ("<< btagger <<")is not commisioned in 2012, please use CSV, TCHP or JetProb"<<endl;
    exit(1); 
  }
  else if(btagger == "TCHPL")
     workingpointvalue = 1.470;
  else if(btagger == "TCHPT")
    workingpointvalue = 3.42;
  else if(btagger == "CSVL")
     workingpointvalue = .244;	
  else if(btagger == "CSVM")
    workingpointvalue = .679;
  else if(btagger == "CSVT")
    workingpointvalue = .898;

  clock_t start = clock();

 cout << "*************************************************************" << endl;
  cout << " Beginning of the program for the FourTop search ! "           << endl;
  cout << "*************************************************************" << endl;

  //SetStyle if needed
 // setTDRStyle();
  setGregStyle();
  //setMyStyle();

  string postfix = "noRho_EventSelection"; // to relabel the names of the output file

  if (doJESShift == 1)
    postfix= postfix+"_JESMinus";
  if (doJESShift == 2)
    postfix= postfix+"_JESPlus";
  if (doJERShift == 1)
    postfix= postfix+"_JERMinus";
  if (doJERShift == 2)
    postfix= postfix+"_JERPlus";
  if (dobTagEffShift == -1)
    postfix= postfix+"_bTagMinus";
  if (dobTagEffShift == 1)
    postfix= postfix+"_bTagPlus";
  if(domisTagEffShift == -1)
    postfix= postfix+"_misTagMinus";
  if(domisTagEffShift == 1)
    postfix= postfix+"_misTagPlus";
 if(doPUShift == 1)
    postfix= postfix+"_PUMinus";
  if(doPUShift == 2)
    postfix= postfix+"_PUPlus";

  ///////////////////////////////////////
  // Configuration
  ///////////////////////////////////////

  string channelpostfix = "";
  string xmlFileName = "";

  bool Electron = true; // use Electron channel?
  bool Muon = false; // use Muon channel?
  if(Electron && Muon){
	cout << "  --> Using both Muon and Electron channel? Choose only one ( since different samples/skims are required)!" << endl;
	exit(1);
  }

  if(Muon){
	cout << " --> Using the Muon channel..." << endl;
	channelpostfix = "_Mu";
	//xmlFileName = "../config/myTopFCNCconfig_Muon.xml";
      xmlFileName = "config/test_fullsamples.xml";

  }
  else if(Electron){
	cout << " --> Using the Electron channel..." << endl;
	channelpostfix = "_El";
//	xmlFileName = "../config/myTopFCNCconfig_Electron.xml";
      xmlFileName = "config/test_fullsamples_el.xml";

  }

//    xmlFileName = "config/test_fullsamples.xml";

   // xmlFileName = "config/test_quicktest.xml";
     // xmlFileName = "config/test_dimuon.xml";
     //xmlFileName = "config/test_2.xml";
      // xmlFileName = "config/refsel.xml";


  bool TrainMVA = false; // If false, the previously trained MVA will be used to calculate stuff
  bool trainEventMVA = false; // If false, the previously trained MVA will be used to calculate stuff
  bool computeEventMVA = true;

 if (doScaleShift==1){
  xmlFileName = "config/test_ScaleDown_el.xml";
}
  else if (doScaleShift==2){ 
 xmlFileName = "config/test_ScaleUp_el.xml";
  }
  else  if (doMatchingShift==1){
  xmlFileName = "config/test_MatchingDown_el.xml";
  }
  else  if (doMatchingShift==2){
  xmlFileName = "config/test_MatchingUp_el.xml";
  }

 else  if (doJESShift!=0){
   xmlFileName = "config/test_mconly_el.xml";
  }

 else  if (doLeptonSFShift!=0){
   xmlFileName = "config/test_mconly_el.xml";
 }


 else  if (doJERShift!=0){
   xmlFileName = "config/test_mconly_el.xml";
  }

else  if (dobTagEffShift!=0){
  xmlFileName = "config/test_mconly_el.xml";
  }

else  if (domisTagEffShift!=0){
  xmlFileName = "config/test_mconly_el.xml";
  }

else  if (doPUShift!=0){
  xmlFileName = "config/test_mconly_el.xml";
  }

else  if (dottbbShift!=0){
  xmlFileName = "config/test_mconly_el.xml";
  }

else  if (trainEventMVA!=0){
  xmlFileName = "config/train_fullsamples_el_rereco.xml";
  }

  else{
 xmlFileName = "config/test_fullsamples_el_rereco.xml";
  }


  const char *xmlfile = xmlFileName.c_str();
  cout << "used config file: " << xmlfile << endl;    

/////////////////////////////
/// AnalysisEnvironment
/////////////////////////////

  AnalysisEnvironment anaEnv;
  cout<<" - Loading environment ..."<<endl;
  AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
  int verbose = 2;//anaEnv.Verbose;

////////////////////////////////
//  Load datasets
////////////////////////////////
    
  TTreeLoader treeLoader;
  vector < Dataset* > datasets;
  cout << " - Loading datasets ..." << endl;
  treeLoader.LoadDatasets (datasets, xmlfile);
  float Luminosity = 19721.734; //pb^-1??
  vector<string> MVAvars;

  //A few bools to steer the MassReco and Event MVAs
  string MVAmethod = "BDT"; // MVAmethod to be used to get the good jet combi calculation (not for training! this is chosen in the jetcombiner class)
  //  bool TrainMVA = false; // If false, the previously trained MVA will be used to calculate stuff
  //bool trainEventMVA = false; // If false, the previously trained MVA will be used to calculate stuff
  //bool computeEventMVA = true;
    cout <<"Instantiating jet combiner..."<<endl;

  JetCombiner* jetCombiner = new JetCombiner(TrainMVA, Luminosity, datasets, MVAmethod, false);
    cout <<"Instantiated jet combiner..."<<endl;
    double bestTopMass1 =0.;
    double bestTopMass2 = 0.;
    double bestTopMass2ndPass = 0.;
    double bestTopPt =0.;

    //uncomment these for training
    //      MVAComputer* Eventcomputer_ =0;
    //MVATrainer* Eventtrainer_ = new MVATrainer("BDT","MasterMVA", "MasterMVA.root");

    //comment this for training
      MVATrainer* Eventtrainer_ = 0;
        
    if(trainEventMVA){
        cout<<"instantiating trainer..."<<endl;
        
    Eventtrainer_->bookWeight("Weight");
        
   // Eventtrainer_->bookInputVar("H");
   // Eventtrainer_->bookInputVar("HT");
    Eventtrainer_->bookInputVar("HTX");
    Eventtrainer_->bookInputVar("HTH");
    Eventtrainer_->bookInputVar("HTRat");
    Eventtrainer_->bookInputVar("HTb");
    // Eventtrainer_->bookInputVar("HTXHX");
  //  Eventtrainer_->bookInputVar("EventMass");
    // Eventtrainer_->bookInputVar("EventMassX");
  //  Eventtrainer_->bookInputVar("SumJetMass");
    Eventtrainer_->bookInputVar("SumJetMassX");
 //   Eventtrainer_->bookInputVar("PTBalTopEventX");
 //   Eventtrainer_->bookInputVar("PTBalTopSumJetX");
    Eventtrainer_->bookInputVar("MultiTopness");
    //  Eventtrainer_->bookInputVar("nJets");
    Eventtrainer_->bookInputVar("nTags");
 //    Eventtrainer_->bookInputVar("bTagRatio13");
 //  Eventtrainer_->bookInputVar("bTagRatio14");
   // Eventtrainer_->bookInputVar("Jet1Pt");
   // Eventtrainer_->bookInputVar("Jet2Pt");
   // Eventtrainer_->bookInputVar("Jet3Pt");
   // Eventtrainer_->bookInputVar("Jet4Pt");
    Eventtrainer_->bookInputVar("Jet5Pt");
    Eventtrainer_->bookInputVar("Jet6Pt");


    }else if (computeEventMVA){
    
  //  MVAvars.push_back("H");
  //  MVAvars.push_back("HT");
      MVAvars.push_back("HTX");
      MVAvars.push_back("HTH");
      MVAvars.push_back("HTRat");
      MVAvars.push_back("HTb");
      //      MVAvars.push_back("HTXHX");
      // MVAvars.push_back("EventMass");
      // MVAvars.push_back("EventMassX");
      // MVAvars.push_back("SumJetMass");
      MVAvars.push_back("SumJetMassX");
      // MVAvars.push_back("PTBalTopEventX");
      // MVAvars.push_back("PTBalTopSumJetX");
      MVAvars.push_back("MultiTopness");
      //MVAvars.push_back("nJets");
      MVAvars.push_back("nTags");
      // MVAvars.push_back("bTagRatio13");
      //MVAvars.push_back("bTagRatio14");
      //  MVAvars.push_back("Jet1Pt");
      //  MVAvars.push_back("Jet2Pt");
      //  MVAvars.push_back("Jet3Pt");
      //  MVAvars.push_back("Jet4Pt");
        MVAvars.push_back("Jet5Pt");
        MVAvars.push_back("Jet6Pt");

    cout << " Initialized Eventcomputer_" << endl;

    // MVAComputer* Eventcomputer_ = new MVAComputer("BDT","MasterMVA.root","MasterMVA",MVAvars, "test_El");

    }
    
  MVAComputer* Eventcomputer_ = new MVAComputer("BDT","MasterMVA.root","MasterMVA",MVAvars, "test_El");

    
  for (unsigned int d = 0; d < datasets.size (); d++)
  {
     cout <<"found sample with equivalent lumi "<<  datasets[d]->EquivalentLumi() <<endl;
     string dataSetName = datasets[d]->Name();
   if(dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0)
	{
		  Luminosity = datasets[d]->EquivalentLumi();
		  cout <<"found DATA sample with equivalent lumi "<<  datasets[d]->EquivalentLumi() <<endl;
		  break;
	 }
  }

  cout << "Rescaling to an integrated luminosity of "<< Luminosity <<" pb^-1" << endl;
    int ndatasets = datasets.size() - 1 ;

    // frac is the factor by which the TTJets sample is divided, it is minimally 1.4
    double frac =2.;
    double currentLumi;
    double newlumi;
    for (unsigned int d = 0; d < datasets.size (); d++)
    {
    if( datasets[d]->Name()=="TTJets" ||datasets[d]->Name()=="TTJets_AllHad" || datasets[d]->Name()=="TTJets_Other"  ){
     currentLumi = datasets[d]->EquivalentLumi();
    cout <<"Old lumi =   "<< currentLumi  <<endl;
     newlumi = currentLumi/frac;
    datasets[d]->SetEquivalentLuminosity(newlumi);
    }
    }
    
    // for splitting the ttbar sample, it is essential to have the ttjets sample as the last
    //dataset loaded
    if (split_ttbar){
        cout << " - splitting TTBar dataset ..." << ndatasets   << endl;
        vector<string> ttbar_filenames = datasets[ndatasets]->Filenames();
        cout <<"ttbar filenames =  "<< ttbar_filenames[0] <<endl;
        
        Dataset* ttbar_ll = new Dataset("TTJets_ll","tt + ll" , true, 633, 2, 2, 1, 213.4,ttbar_filenames );
        Dataset* ttbar_cc = new Dataset("TTJets_cc","tt + cc" , true, 633, 2, 2, 1, 6.9, ttbar_filenames );
        Dataset* ttbar_bb = new Dataset("TTJets_bb","tt + bb" , true, 633, 2, 2, 1, 4.8, ttbar_filenames );
        
	///heavy flav re-weight -> scaling ttbb up and ttjj down so that ttbb/ttjj matches CMS measurement.
	double ll_rw = 0.976;
	double bb_rw = 3.;

         ttbar_ll->SetEquivalentLuminosity(newlumi/ll_rw);
	 ttbar_cc->SetEquivalentLuminosity(newlumi/ll_rw);
	 ttbar_bb->SetEquivalentLuminosity(newlumi/bb_rw);
        
	 //   ttbar_ll->SetEquivalentLuminosity(newlumi);
	 // ttbar_cc->SetEquivalentLuminosity(newlumi);
	 // ttbar_bb->SetEquivalentLuminosity(newlumi);


        ttbar_ll->SetColor(kRed);
        ttbar_cc->SetColor(kRed-3);
        ttbar_bb->SetColor(kRed+2);
        
        
        datasets.pop_back();
        datasets.push_back(ttbar_bb);
        datasets.push_back(ttbar_cc);
        datasets.push_back(ttbar_ll);     
    }     
    
  //Output ROOT file
  string rootFileName ("FourTop"+postfix+channelpostfix+".root");
  TFile *fout = new TFile (rootFileName.c_str(), "RECREATE");

  //vector of objects
  cout << " - Variable declaration ..." << endl;
  vector < TRootVertex* >   vertex;
  vector < TRootMuon* >     init_muons;
  vector < TRootElectron* > init_electrons;
  vector < TRootJet* >      init_jets;
  vector < TRootMET* >      mets;

  //Global variable
  TRootEvent* event = 0;
    
    ////////////////////////////////////////////////////////////////////
    ////////////////// MultiSample plots  //////////////////////////////
    ////////////////////////////////////////////////////////////////////
    
    MSPlot["RhoCorrection"]              = new MultiSamplePlot(datasets, "RhoCorrection", 100, 0, 100, "#rho");
    MSPlot["NbOfVertices"]               = new MultiSamplePlot(datasets, "NbOfVertices", 40, 0, 40, "Nb. of vertices");
    //Muons
    MSPlot["NbOfIsolatedMuons"]          = new MultiSamplePlot(datasets, "NbOfIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
    MSPlot["NbOfIsolatedElectrons"]      = new MultiSamplePlot(datasets, "NbOfIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
    MSPlot["NbOfExtraIsolatedMuons"]     = new MultiSamplePlot(datasets, "NbOfExtraIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
    MSPlot["NbOfExtraIsolatedElectrons"] = new MultiSamplePlot(datasets, "NbOfExtraIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
  
    MSPlot["ElectronRelIsolation"] = new MultiSamplePlot(datasets, "ElectronRelIsolation", 50, 0, .25, "RelIso");
    MSPlot["ElectronPt"]              = new MultiSamplePlot(datasets, "ElectronPt", 50, 0, 300, "PT_{e}");
    MSPlot["ElectronEta"]              = new MultiSamplePlot(datasets, "ElectronEta", 25, -2.4, 2.4, "#eta_{e}");
    MSPlot["ElectronPhi"]              = new MultiSamplePlot(datasets, "ElectronPhi", 50, -4, 4, "#phi_{e}");
    MSPlot["Electrond0"]              = new MultiSamplePlot(datasets, "Electrond0", 50, -0.05, 0.05, "d0_{#mu}");
    MSPlot["ElectronMissingHits"]              = new MultiSamplePlot(datasets, "ElectronMissingHits", 11, -0.05,    10.5, "Missing Hits");
    MSPlot["ElectronMVATrigID"]              = new MultiSamplePlot(datasets, "ElectronMVATrigID", 70, 0.7,    1, "MVA Trig ID");

    //    MSPlot["ElectrondZPVz"]              = new MultiSamplePlot(datasets, "ElectrondZPVz", 50, 0, .5, "dZPVZ_{#mu}");

    MSPlot["ElectrondRJets"]              = new MultiSamplePlot(datasets, "ElectrondRJets", 50, 0, 10, "dRJets_{#mu}");   
    //    MSPlot["ElectronDistVzPVz"]              = new MultiSamplePlot(datasets, "ElectronDistVzPVz", 50, 0 ,.3, "DistVzPVz_{#mu}");
    //    MSPlot["ElectronDz"]              = new MultiSamplePlot(datasets, "ElectronDz", 25, -.6 ,.6, "Dz_{#mu}");
   
    MSPlot["NbOfLooseElectron"]     = new MultiSamplePlot(datasets, "NbOfLooseElectron", 10, 0, 10, "Nb. of loose Electrons");
    
    //B-tagging discriminators
    MSPlot["BdiscBJetCand_CSV"]  = new MultiSamplePlot(datasets, "BdiscBJetCand_CSV", 75, 0, 1, "CSV b-disc.");
      
    MSPlot["bTagRatio13"]            = new MultiSamplePlot(datasets, "bTagRatio13", 35, 0, 20, "bTagRatio13");
    MSPlot["bTagRatio14"]            = new MultiSamplePlot(datasets, "bTagRatio14", 35, 0, 20, "bTagRatio14");

    //Jets
    MSPlot["NbOfSelectedJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedJets", 15, 0, 15, "Nb. of jets");
    MSPlot["NbOfSelectedLightJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedLightJets", 10, 0, 10, "Nb. of jets");
    MSPlot["NbOfSelectedBJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedBJets", 8, 0, 8, "Nb. of tags");
    MSPlot["JetEta"]                  = new MultiSamplePlot(datasets, "JetEta", 30,-3, 3, "Jet #eta");
    MSPlot["JetPhi"]                  = new MultiSamplePlot(datasets, "JetPhi", 50, -4,4 , "Jet #phi");
    
    MSPlot["NbOfBadTrijets"]                  = new MultiSamplePlot(datasets, "NbOfBadTriJets", 150, 0, 150, "Nb. of Bad Combs");
    
    MSPlot["TriJetMass_Matched"] = new MultiSamplePlot(datasets, "TriJetMassMatched", 100, 0, 1000, "m_{bjj}");
    MSPlot["TriJetMass_UnMatched"] = new MultiSamplePlot(datasets, "TriJetMassUnMatched", 100, 0, 1000, "m_{bjj}");

    MSPlot["MultiTopness"] = new MultiSamplePlot(datasets, "MultiTopness", 35, -1., 0.3, "MultiTopness");
    
    MSPlot["MVA1TriJetMass"] = new MultiSamplePlot(datasets, "MVA1TriJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA1DiJetMass"] = new MultiSamplePlot(datasets, "MVA1DiJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA1PtRat"] = new MultiSamplePlot(datasets, "MVA1PtRat", 25, 0, 2, "P_{t}^{Rat}");
    MSPlot["MVA1BTag"] = new MultiSamplePlot(datasets, "MVA1BTag", 35, 0, 1, "BTag");
    MSPlot["MVA1AnThBh"] = new MultiSamplePlot(datasets, "MVA1AnThBh", 35, 0, 3.14, "AnThBh");
    MSPlot["MVA1AnThWh"] = new MultiSamplePlot(datasets, "MVA1AnThWh", 35, 0, 3.14, "AnThWh");

    MSPlot["MVA1TriJet"] = new MultiSamplePlot(datasets, "MVA1TriJet", 35, -1., 0, "BDT Discriminator");
    MSPlot["MVA2TriJetMass"] = new MultiSamplePlot(datasets, "MVA2TriJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA2ndPassTriJetMass"] = new MultiSamplePlot(datasets, "MVA2ndPassTriJetMass", 30, 0, 1000, "m_{bjj}");
    MSPlot["MVA2ndPassTriJet"] = new MultiSamplePlot(datasets, "MVA2ndPassTriJet", 35, -1., 0, "BDT Discriminator");
    MSPlot["MVA3rdPassTriJet"] = new MultiSamplePlot(datasets, "MVA3rdPassTriJet", 35, -1., 0, "BDT Discriminator");  
    MSPlot["MVA1TriJetMassMatched"] = new MultiSamplePlot(datasets, "MVA1TriJetMassMatched", 75, 0, 500, "m_{bjj}");
    MSPlot["BDisc_Asym"] = new MultiSamplePlot(datasets, "BDisc_Asym", 100, -5, 5, "BDiscAsym");
    MSPlot["HTComb"] = new MultiSamplePlot(datasets, "HTComb", 100, 0, 800, "HTComb");
    MSPlot["WMassComb"] = new MultiSamplePlot(datasets, "WMassComb", 75, 0, 650, "WMassComb");

    MSPlot["WMt"] = new MultiSamplePlot(datasets, "WMt", 50, 0, 250, "W Transverse Mass");
    MSPlot["LepWMass"] = new MultiSamplePlot(datasets, "LepWMass", 50, 0, 200, "MuMET");
    MSPlot["LepWMass_g"] = new MultiSamplePlot(datasets, "LepWMass_g", 50, 0, 200, "MuMET");
    MSPlot["SelectedJetPt"] = new MultiSamplePlot(datasets, "JetPt", 50, 0, 300, "PT_{jet}");
  
    MSPlot["1stJetPt"] = new MultiSamplePlot(datasets, "1stJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["2ndJetPt"] = new MultiSamplePlot(datasets, "2ndJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["3rdJetPt"] = new MultiSamplePlot(datasets, "3rdJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["4thJetPt"] = new MultiSamplePlot(datasets, "4thJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["5thJetPt"] = new MultiSamplePlot(datasets, "5thJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["6thJetPt"] = new MultiSamplePlot(datasets, "6thJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["7thJetPt"] = new MultiSamplePlot(datasets, "7thJetPt", 60, 0, 400, "PT_{jet}");

    MSPlot["SelectedJetPt_light"] = new MultiSamplePlot(datasets, "JetPt_light", 50, 0, 1000, "PT_{lightjet}");
    MSPlot["SelectedJetPt_b"] = new MultiSamplePlot(datasets, "JetPt_b", 50, 0, 1000, "PT_{bjet}");
    MSPlot["HT_SelectedJets"] = new MultiSamplePlot(datasets, "HT_SelectedJets", 50, 0, 1500, "HT");
    MSPlot["HTb_SelectedJets"] = new MultiSamplePlot(datasets, "HTb_SelectedJets", 50, 0, 1500, "HTb");
    
    MSPlot["H"] = new MultiSamplePlot(datasets, "H", 50, 0, 3000, "H");
    MSPlot["HX"] = new MultiSamplePlot(datasets, "HX", 50, 0, 1500, "HX");
    
    MSPlot["HTH"] = new MultiSamplePlot(datasets, "HTH", 50, 0, 1, "HT/H");
    MSPlot["HTXHX"] = new MultiSamplePlot(datasets, "HTXHX", 50, 0, 1, "HTX/HX");
    /*
    MSPlot["DeltaPhi_LepTop"] = new MultiSamplePlot(datasets, "DeltaPhi_LepTop", 25, -3.14, 3.14, "DeltaPhi_LepTop");
    MSPlot["HT_Asym"] = new MultiSamplePlot(datasets, "HT_Asym", 25, -1., 1., "");
    
    
    MSPlot["MHT_SelectedJets"] = new MultiSamplePlot(datasets, "MHT_SelectedJets", 75, 0, 1000, "MHT");
    MSPlot["MHTSig_SelectedJets"] = new MultiSamplePlot(datasets, "MHTSig_SelectedJets", 75, 0, 30, "MHTSig");*/
    MSPlot["MET"] = new MultiSamplePlot(datasets, "MET", 75, 0, 700, "MET");
    MSPlot["MET_MHT"]= new MultiSamplePlot(datasets, "MET_MHT", 75, 0, 200, "MET_MHT");
    MSPlot["STLep"] = new MultiSamplePlot(datasets, "STLep", 40, 0, 2500, "STLep");
    MSPlot["STJet"] = new MultiSamplePlot(datasets, "STJet", 40, 0, 2500, "STJet");
    
    MSPlot["EventMass"] = new MultiSamplePlot(datasets, "EventMass", 40, 0, 3000, "EventMass");
    MSPlot["EventMassX"] = new MultiSamplePlot(datasets, "EventMassX", 30, 0, 3000, "EventMassX");
    MSPlot["SumJetMass"] = new MultiSamplePlot(datasets, "SumJetMass", 40, 0, 3000, "SumJetMass");
    MSPlot["SumJetMassX"] = new MultiSamplePlot(datasets, "SumJetMassX", 30, 0, 3000, "SumJetMassX");
    MSPlot["HTX"] = new MultiSamplePlot(datasets, "HTX", 20, 0, 1000, "HTX");
    MSPlot["MVA"] = new MultiSamplePlot(datasets, "MVA", 17, -0.5, 0.4, "BDT Discriminator");
  
    MSPlot["HTRat"]                  = new MultiSamplePlot(datasets, "HTRat", 50, 0, 20, "HTRat");
    MSPlot["H_SelectedJets"]                  = new MultiSamplePlot(datasets, "H_SelectedJets", 50, 0, 2000, "H");
  


   //Plots after loose cut on BDT discriminator
    MSPlot["NbOfSelectedJets_BDTCut"]                  = new MultiSamplePlot(datasets, "NbOfSelectedJets_BDTCut", 15, 0, 15, "Nb. of jets");
    MSPlot["NbOfSelectedBJets_BDTCut"]                  = new MultiSamplePlot(datasets, "NbOfSelectedBJets_BDTCut", 8, 0, 8, "Nb. of tags");
    MSPlot["HTX_BDTCut"] = new MultiSamplePlot(datasets, "HTX_BDTCut", 20, 0, 1000, "HTX");
    MSPlot["HTH_BDTCut"] = new MultiSamplePlot(datasets, "HTH_BDTCut", 50, 0, 1, "HT/H");
    MSPlot["MultiTopness_BDTCut"] = new MultiSamplePlot(datasets, "MultiTopness_BDTCut", 35, -1., 0.5, "MultiTopness");
    MSPlot["HTb_SelectedJets_BDTCut"] = new MultiSamplePlot(datasets, "HTb_SelectedJets_BDTCut", 50, 0, 1500, "HTb");

 histo1D["btag_weight"] = new TH1F("btag_weight","btag_weight;btag_weight;#events",150,-10,10);  



 //Declare arrays of MSPlots
    Int_t minNJets=6, maxNJets=8, minNBJets=0, maxNBJets=3;
    Int_t q =0;

    for (Int_t p = minNJets; p<= maxNJets; p++){

    string NJets_str = static_cast<ostringstream*>( &(ostringstream() << p) )->str();
    string MVA_Name = "MVA"+NJets_str+"Jets" ;        

    MSPlot[MVA_Name.c_str()] = new MultiSamplePlot(datasets, MVA_Name.c_str(), 17, -0.5, 0.4, "BDT Discriminator");

}

    for (Int_t p = minNBJets; p<= maxNBJets; p++){

    string NBJets_str = static_cast<ostringstream*>( &(ostringstream() << p) )->str();
    string MVA_Name = "MVA"+NBJets_str+"Tags";        

    MSPlot[MVA_Name.c_str()] = new MultiSamplePlot(datasets, MVA_Name.c_str(), 17, -0.5, 0.4, "BDT Discriminator");

 }



  //Plots
  string pathPNG_MVA = "MVAPlots_"+postfix+channelpostfix;
  string pathPNG = "FourTop"+postfix+channelpostfix;
  pathPNG += "_MSPlots/"; 	
  //pathPNG = pathPNG +"/";
  mkdir(pathPNG.c_str(),0777);
    
    cout <<"Making directory :"<< pathPNG  <<endl;
    
 

  /////////////////////////////
  // Selection table: e + jets
  /////////////////////////////

  vector<string> CutsSelecTableEl;
  // CutsSelecTableEl.push_back(string("initial"));
  //  CutsSelecTableEl.push_back(string("PU reweighting"));
  CutsSelecTableEl.push_back(string("Trigger and PV"));
  CutsSelecTableEl.push_back(string("1 iso. e"));
  CutsSelecTableEl.push_back(string("Loose e veto"));
  CutsSelecTableEl.push_back(string("Loose mu veto"));
  CutsSelecTableEl.push_back(string("$>$= 6 Jets"));
  CutsSelecTableEl.push_back(string("$>$= 2 b-tags"));
  CutsSelecTableEl.push_back(string("HT $>$= 400 GeV"));
  CutsSelecTableEl.push_back(string("E^{Miss}_{T} > 30 GeV"));

  FourTopSelectionTable selecTableEl(CutsSelecTableEl, datasets);
  selecTableEl.SetLuminosity(Luminosity);
  selecTableEl.SetPrecision(1);
  

////////////////////////////
//  Pile up reweighting
////////////////////////////
//NEW METHOD (TRUE INTERACTIONS)

    LumiReWeighting LumiWeights;

    if(doPUShift==1){

 LumiWeights = LumiReWeighting("PUReweighting/pileup_MC_S10.root", "PUReweighting/PURewtoptree_id_3427_3429_3446_3444__PileupHistogram_Systematic_Down_5perc.root", "pileup", "pileup");


}
    else if (doPUShift==2){

 LumiWeights = LumiReWeighting("PUReweighting/pileup_MC_S10.root", "PUReweighting/PURewtoptree_id_3427_3429_3446_3444__PileupHistogram_Systematic_Up_5perc.root", "pileup", "pileup");

    }else if(doPUShift==0){

    LumiWeights = LumiReWeighting("PUReweighting/pileup_MC_S10.root", "PUReweighting/PURewtoptree_id_2217_2227_2262_2423_2435_2417_PileupHistogram.root", "pileup", "pileup");
    }


    //reweight::PoissonMeanShifter PShiftDown_ = reweight::PoissonMeanShifter(-0.6);
    //reweight::PoissonMeanShifter PShiftUp_ = reweight::PoissonMeanShifter(0.6);

   cout << " Initialized LumiReWeighting stuff" << endl;

    
    // initialize lepton SF
    LeptonTools* leptonTools = new LeptonTools(false);
    // leptonTools->readMuonSF("LeptonSF/Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.root", "LeptonSF/MuonEfficiencies_Run_2012A_2012B_53X.root", "LeptonSF/MuonEfficiencies_Run_2012C_53X.root", "LeptonSF/TriggerMuonEfficiencies_Run_2012D_53X.root");
    leptonTools->readElectronSF();
    
    
   vector<TRootMCParticle*> mcParticles_flav;
  TRootGenEvent* genEvt_flav = 0;
  
/////////////////////////////////
//////// Loop on datasets
/////////////////////////////////
    
    TGraph *scaleGraph, *BTagGraph_B,*BTagGraph_C, *BTagGraph_L_eta008, *BTagGraph_L_eta0816, *BTagGraph_L_eta1624; ;

    
    TFile * btageffFile;
    btageffFile = new TFile("FourTop_MCStudy_El.root","READ");
    

  cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
    
  for (unsigned int d = 0; d < datasets.size(); d++) //d < datasets.size()
  {
    if (verbose > 1){
      cout << "   Dataset " << d << " name : " << datasets[d]->Name () << " / title : " << datasets[d]->Title () << endl;
      cout << " - Cross section = " << datasets[d]->Xsection() << endl;
      cout << " - IntLumi = " << datasets[d]->EquivalentLumi() << "  NormFactor = " << datasets[d]->NormFactor() << endl;
      cout << " - Nb of events : " << datasets[d]->NofEvtsToRunOver() << endl;
    }
    //open files and load
    cout<<"Load Dataset"<<endl;
    treeLoader.LoadDataset (datasets[d], anaEnv);
		
    string previousFilename = "";
    int iFile = -1;
    
    string dataSetName = datasets[d]->Name();	
    
//////////////////////////////////////////////////
/// Initialize JEC factors ///////////////////////
//////////////////////////////////////////////////
		
    vector<JetCorrectorParameters> vCorrParam;
    JetCorrectorParameters *L1JetCorPar ,*L2JetCorPar, *L3JetCorPar, *L2L3ResJetCorPar;

    if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) // Data!
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L1FastJet_AK5PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L2Relative_AK5PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L3Absolute_AK5PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
      JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L2L3Residual_AK5PFchs.txt");
      vCorrParam.push_back(*L2L3ResJetCorPar);
    }
    else
    {
      JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_L1FastJet_AK5PFchs.txt");
      vCorrParam.push_back(*L1JetCorPar);
      JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_L2Relative_AK5PFchs.txt");
      vCorrParam.push_back(*L2JetCorPar);
      JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_L3Absolute_AK5PFchs.txt");
      vCorrParam.push_back(*L3JetCorPar);
    }
    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("../../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_Uncertainty_AK5PFchs.txt");
//    JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Fall12_V7_DATA_UncertaintySources_AK5PFchs.txt", "SubTotalMC")));
//    JetCorrectionUncertainty *jecUncTotal = new JetCorrectionUncertainty(*(new JetCorrectorParameters("JECFiles/Fall12_V7_DATA_UncertaintySources_AK5PFchs.txt", "Total")));
    
    JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);		


//////////////////////////////////////////////////
//      Loop on events
/////////////////////////////////////////////////

    int itrigger = -1, previousRun = -1;
   
    int start = 0;
    unsigned int end = datasets[d]->NofEvtsToRunOver();

    cout <<"Number of events = "<<  end  <<endl;
      
    bool debug = false;
    int event_start;
    if (dataSetName == "Data")TrainMVA=false, trainEventMVA=false;
      
     

    if (verbose > 1) cout << " - Loop over events " << endl;
      
      int nBBBar, nCCBar, nLLBar;
      nBBBar=  nCCBar = nLLBar = 0;
      
      double MHT, MHTSig, STJet, EventMass, EventMassX , SumJetMass, SumJetMassX,H,HX , HT, HTX,HTH,HTXHX, sumpx_X, sumpy_X, sumpz_X, sume_X, sumpx, sumpy, sumpz, sume, jetpt,PTBalTopEventX,PTBalTopSumJetX , PTBalTopMuMet;
        
      double currentfrac =0.;
      double end_d = end;
      
      double fakeSig_xs = 100.; //xsection of fake signal in fb.
      
      // to control what fraction of the ttjets we run over(remember to alter the int. lumi of the sample)
      if(dataSetName=="TTJets" || dataSetName=="TTJets_ll" || dataSetName=="TTJets_cc" || dataSetName=="TTJets_bb" ) {

	if(trainEventMVA){

	event_start = 200000;
	 end_d =      1000000.;
	}else if (computeEventMVA){

	event_start = 1000000;
        end_d = 1000000. + end_d/frac;

}
 
      }

      else if (datasets[d]->Name()=="TTJets_AllHad" || datasets[d]->Name()=="TTJets_Other"  ) {
          event_start = 0;
          end_d =  end_d/frac;
          cout <<"N ttjets to run over =   "<< end_d/frac <<endl;
      }
else if(dataSetName=="NP_overlay_TTTT"){

	if(trainEventMVA){

	event_start = 0;
	 end_d = 20000;
	}else if (computeEventMVA){

	event_start = 20000.;
	 end_d = end;
}

      }else if(dataSetName=="Data"){
	event_start = 0;
          end_d = end;
      }
      
      else{
          event_start = 0.;
          end_d = end;
      }
      cout <<"Will run over "<<  end_d<< " events..."<<endl;
      cout <<"Starting event =  "<< event_start  << "Ending event =  "<< end_d   << endl; 
    for (unsigned int ievt = event_start; ievt < end_d ; ievt++)
    {
        
MHT = 0.,MHTSig = 0., STJet = 0., EventMass =0., EventMassX =0., SumJetMass = 0., SumJetMassX=0.  ,H = 0., HX =0., HT = 0., HTX = 0.,HTH=0.,HTXHX=0., sumpx_X = 0., sumpy_X= 0., sumpz_X =0., sume_X= 0. , sumpx =0., sumpy=0., sumpz=0., sume=0., jetpt =0., PTBalTopEventX = 0., PTBalTopSumJetX =0.;
        
        double ievt_d = ievt;
        currentfrac = ievt_d/end_d;
        if (debug)cout <<"event loop 1"<<end <<"  "<<start<<endl;
        
	if(ievt%10000 == 0)
		std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";

        if (debug)cout <<"event loop 2"<<endl;

	// scale factor for the event
	float scaleFactor = 1.;
        
	//load event
	event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);

//cout <<" "<<endl;
	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"First cout after loading event:");

////////////////////////////////////////////////////////////////////////////
/// Splitting TTBar sample into tt +ll, tt+ cc and tt + bb /////////////////
////////////////////////////////////////////////////////////////////////////
        
        if (debug)cout <<"event loop 4"<<endl;
        
    //load mcparticles to check jet flavour for ttjets events
	//  vector<TRootMCParticle*> mcParticles_flav;
    Int_t ttbar_flav = -1;

        double nExB,nExC,nExL;
        nExB = nExC = nExL = 0.;
        
  
    genEvt_flav = treeLoader.LoadGenEvent(ievt,false);
    treeLoader.LoadMCEvent(ievt, genEvt_flav, 0, mcParticles_flav,false); 
        
   
 /////////////////////////////
 /////  Generating shape sytematic by reweighting tt+bb events
 /////////////////////////////
 //	int hf_rw = 2;

	double rw = 1.;

	//	if (dottbbShift != 0){

  if(dataSetName == "TTJets" )
    {

  //load mcparticles to check jet flavour for ttjets events

        double nExB,nExC,nExL;
        nExB = nExC = nExL = 0.;
        
	//  TRootGenEvent* genEvt_flav = 0;
    genEvt_flav = treeLoader.LoadGenEvent(ievt,false);
    treeLoader.LoadMCEvent(ievt, genEvt_flav, 0, mcParticles_flav,false);


       for(unsigned int p=0; p<mcParticles_flav.size(); p++) {
          
            if(mcParticles_flav[p]->status()==3 && abs(mcParticles_flav[p]->type())==5 && abs(mcParticles_flav[p]->motherType())!=6) {
	      // ttbar_flav=2;
                nExB++;  
            }
            
            else if (mcParticles_flav[p]->status()==3 && abs(mcParticles_flav[p]->type())==4 && abs(mcParticles_flav[p]->motherType())!=6
                && abs(mcParticles_flav[p]->motherType())!=5 && abs(mcParticles_flav[p]->motherType())!=24
                ){
	      //ttbar_flav=1;
                 nExC++; 
            }
            
            else if (mcParticles_flav[p]->status()==3 && abs(mcParticles_flav[p]->type())<4 && abs(mcParticles_flav[p]->motherType())!=6){
	      // ttbar_flav=1;
                nExL++; 
            }

        }
            
     //   cout <<"TTBar flav composition : " << nExL  <<"  light, " << nExC <<"  C, " << nExB<< "  B" <<  endl;            
       
            if (nExB >= 2.){
            ttbar_flav =2; 
            nBBBar++ ; //  bbbar
            }
            else if ( nExC >=2.) {
            ttbar_flav =1; 
            nCCBar++ ; //ccbar
            }
            else{
            ttbar_flav =0.; 
                nLLBar++;  //llbar   
            }


            if (ttbar_flav ==0 && (dataSetName == "TTJets_cc"  || dataSetName == "TTJets_bb"))  continue;
            if (ttbar_flav ==1 && (dataSetName == "TTJets_ll"  || dataSetName == "TTJets_bb" ))  continue;
            if (ttbar_flav ==2 && (dataSetName == "TTJets_ll"  || dataSetName == "TTJets_cc" ))  continue;


	    if (ttbar_flav ==2){

	      //try scaling ttbb up by 300%
	      
	    if (dottbbShift == 1 ) {

	      rw = 1.5;

	    }else if (dottbbShift == 2){

               rw = 0.5;

	    }else{
	      rw = 1.;
}

    }

}


  //    cout <<"ll  "<< nExL  <<"  cc  "  <<   nExC  << "  bb  "  <<  nExB  <<  " ttbar flav "<< ttbar_flav  <<endl;


 if(dataSetName != "Data"){
	scaleFactor = scaleFactor*rw;
 }
 

//////////////////
//Top pt reweighting
//////////////////
	double SF_t = 1;
	double SF_tbar = 1;
	double t_pt;
	double tbar_pt;
	double a = 0.159;
	double b = -0.00141;
	double t_pt_rw =1.;

	bool applyTopPtRW = false;
  if(applyTopPtRW){
 for ( int i=0; i<mcParticles_flav.size(); i++){
                 if( mcParticles_flav[i]->status() != 3) continue;
                 
                 if(  mcParticles_flav[i]->type()== 6){
		   t_pt = mcParticles_flav[i]->Pt();
		   SF_t =exp(a + (b*t_pt ));
                 }else if (  mcParticles_flav[i]->type()== -6){
                   tbar_pt = mcParticles_flav[i]->Pt();
		   SF_tbar =exp(a + (b*tbar_pt ));
}


 }

 t_pt_rw  =sqrt(SF_t*SF_tbar);


 // cout <<" SF PT "<< t_pt_rw   <<endl;

 if(dataSetName != "Data"){
 scaleFactor = scaleFactor*t_pt_rw;
 }

  }

//////////////////
//Loading Gen jets
//////////////////
        if (debug)cout <<"event loop 5"<<endl;

        
	vector<TRootGenJet*> genjets;
	if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
	{
	  // loading GenJets as I need them for JER
	  		genjets = treeLoader.LoadGenJet(ievt);
	}

	// check which file in the dataset it is to have the HLTInfo right
	string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
	if(previousFilename != currentFilename)
	{
		previousFilename = currentFilename;
        	iFile++;
		cout<<"File changed!!! => iFile = "<<iFile<<endl;
	}
        
///////////////////////////////////////////
//  Trigger
///////////////////////////////////////////
	bool trigged = false;
    std::string filterName = "";
	int currentRun = event->runId();
	if(previousRun != currentRun)
	{
	 // cout <<"What run? "<< currentRun<<endl;
      		previousRun = currentRun;
	
			if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
			{
                if (debug)cout <<"event loop 6a"<<endl;

		// cout << " RUN " << event->runId() << endl;
			
                if( event->runId() <= 190738 )
                    itrigger = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v8"), currentRun, iFile);
                else if( event->runId() >= 190782 && event->runId() <= 191411 )
                    itrigger = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v9"), currentRun, iFile);
                else if( event->runId() >= 191695 && event->runId() <= 196531)
                    itrigger = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v10"), currentRun, iFile);
                else if( event->runId() >= 198049 && event->runId() <= 208686)
                    itrigger = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v11"), currentRun, iFile);
		// else if( event->runId() > 208686)
		//     itrigger = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v11"), currentRun, iFile);
       	    		

	  if(itrigger == 9999)
				{
    		  cout << "NO VALID TRIGGER FOUND FOR THIS EVENT (DATA) IN RUN " << event->runId() << endl;
		 //   exit(1);
 	 			}
			}
	   		else 
			  {
				 itrigger = treeLoader.iTrigger (string ("HLT_Ele27_WP80_v10"), currentRun, iFile);
    
  				if(itrigger == 9999)
				{
    			  		cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;
						//exit(1);
				}
			}
	
	} //end previousRun != currentRun

      
////////////////////////////////////////////////////////////////////////////////////
//  JES Corrections: The nominal corrections are already applied at PAT level     //
//    so these tools should only be used for studies of the effect of systematics //
////////////////////////////////////////////////////////////////////////////////////

    // Apply Jet Corrections on-the-fly
      if( dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
      {
        jetTools->unCorrectMETTypeOne(init_jets, mets[0], true);
        jetTools->correctJets(init_jets, event->kt6PFJets_rho(), true);
        jetTools->correctMETTypeOne(init_jets, mets[0], true);
      }
      else
      {
        jetTools->unCorrectMETTypeOne(init_jets, mets[0], false);
        jetTools->correctJets(init_jets, event->kt6PFJets_rho(), false);
        jetTools->correctMETTypeOne(init_jets, mets[0], false);
      }

  
        ///////////////////////
        // JER smearing
        //////////////////////
        
        if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
        {
            //JER
            doJERShift == 0;
            if(doJERShift == 1)
	      jetTools->correctJetJER(init_jets, genjets, mets[0], "minus",false);
            else if(doJERShift == 2)
	      jetTools->correctJetJER(init_jets, genjets, mets[0], "plus",false);
            else
	      jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal", false);
            
            //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JER correction:");
            
            // JES systematic!
            if (doJESShift == 1)
                jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
            else if (doJESShift == 2)
                jetTools->correctJetJESUnc(init_jets, mets[0], "plus");
            
            //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction:");
            
        }
        
    
    
////////////////////////////////////////
//  Beam scraping and PU reweighting
////////////////////////////////////////
        
	// scale factor for the event
	//	float scaleFactor = 1.;

	if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
	{
        if (debug)cout <<"event loop 6bb"<<endl;

		// Apply the scraping veto. (Is it still needed?)
        	bool isBeamBG = true;
        	if(event->nTracks() > 10)
        	{
                if (debug)cout <<"event loop 6bcxx"<<endl;
                
			if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
			isBeamBG = false;
		}
      		if(isBeamBG) continue;
        if (debug)cout <<"event loop 6bcc"<<endl;

        
	}
	else{
	double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
	double lumiWeightOLD=lumiWeight;
        if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0){
	  lumiWeight=1;
        }
	scaleFactor = scaleFactor*lumiWeight;
        if (debug)cout <<"event loop 6c"<<endl;

	}

	//histo1D["lumiWeights"]->Fill(scaleFactor);
			
///////////////////////////////////////////////////////////
//   Event selection
///////////////////////////////////////////////////////////
	       
	// Apply trigger selection
	trigged = treeLoader.EventTrigged (itrigger);
	if (debug)cout<<"triggered? Y/N?  "<< trigged  <<endl;

        if (debug)cout <<"event loop 6bbb"<<endl;

	//Applying trigger selection again with 2012 Electron triggers.
	if(!trigged)		   continue;

	// Declare selection instance    
	Selection selection(init_jets, init_muons, init_electrons, mets, event->kt6PFJets_rho());

	// Define object selection cuts
	selection.setJetCuts(30.,2.4,0.01,1.,0.98,0,0);//Pt, Eta, EMF, n90Hits, fHPD, dRJetElectron, DRJets

	//    selection.setElectronCuts();//	selection.setElectronCuts(30.,2.5,0.1,0.02,0.,999.,0,1); //Pt,Eta,RelIso,d0,MVAId, DRJets, MaxMissingHits
        
        selection.setElectronCuts(30.,2.5,0.1,0.02,0.9,999.,0);
        
	selection.setLooseElectronCuts(20,2.5,0.15,0.);

	//selection.setMuonCuts(26.,2.1,.12,0,0.2,0,1,0.5,5,1 ); //Pt,Eta,RelIso,NValidMuHits,d0, dRJets, NMatchedStations,DistVzPVz,NTrackerLayersWithMeas
    selection.setMuonCuts(26.0,2.1,.12,0.2,0.3,1,0.5,5,0 );
    selection.setLooseMuonCuts(10,2.5,0.2);
	  
	//Select objects 
	//	vector<TRootElectron*> selectedElectrons_NoIso = selection.GetSelectedElectrons(20,2.4,999.,vertex[0]);
	//vector<TRootElectron*> selectedElectrons        = selection.GetSelectedElectrons(vertex[0]);
   
    vector<TRootElectron*> selectedElectrons        = selection.GetSelectedElectrons();
	vector<TRootElectron*> selectedExtraElectrons;
	vector<TRootMuon*>     selectedMuons_NoIso      = selection.GetSelectedMuons(26,2.4,999.); 
	vector<TRootMuon*>     selectedMuons            = selection.GetSelectedMuons(vertex[0]);
	vector<TRootMuon*>     selectedExtraMuons;
	vector<TRootJet*>      selectedJets             = selection.GetSelectedJets(true); // ApplyJetId
    vector<TRootJet*>      selectedJets2ndPass;
    vector<TRootJet*>      selectedJets3rdPass;
    vector<TRootJet*>   MVASelJets1;

    vector<TRootElectron*> selectedElectrons_NoIso  = selection.GetSelectedElectrons(30,2.4,999.);
        
    vector<TRootJet*>      selectedSoftJets         = selection.GetSelectedJets(20.,2.5, selectedMuons, 0., true); // ApplyJetId
	vector<TRootMuon*>     selectedLooseMuons       = selection.GetSelectedLooseMuons(10., 2.5, 0.2);
    vector<TRootElectron*> selectedLooseElectrons   = selection.GetSelectedLooseElectrons(20.,2.5,0.15); // VBTF ID
    vector<TRootJet*>      selectedBJets; // B-Jets
    vector<TRootJet*>      selectedLightJets; // light-Jets
    vector<TRootJet*>       selectedCSVOrderedJets     = selection.GetSelectedJets(true); //CSV ordered Jet collection added by JH

	//order jets wrt to Pt, then set bool corresponding to RefSel cuts.
    sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order muons wrt Pt.                                                                    
    sort(selectedCSVOrderedJets.begin(), selectedCSVOrderedJets.end(), HighestCVSBtag()); //order Jets wrt CSVtag

    int JetCut =0;
    int nMu = selectedMuons.size();
    int nEl = selectedElectrons.size();
  
  selecTableEl.Fill(d,0,1.);
	// Apply primary vertex selection
	bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
        if(!isGoodPV) continue;

	double temp_HT = 0.;
	double HT_leading = 0.;
	double HT_lagging = 0.;
	double HTRat = 0;
	for (Int_t seljet0 =0; seljet0 < selectedJets.size(); seljet0++ ){
	  temp_HT += selectedJets[seljet0]->Pt();
	  if (seljet0 < 4){
	    HT_leading += selectedJets[seljet0]->Pt();
	  }else{

	    HT_lagging += selectedJets[seljet0]->Pt();
	  }
        }

        if (debug)cout <<"event loop 6"<<endl;

        bool isTagged =false;
        int seljet;
        
        
        //boolean which controls the application of the Btag scalefactor
        int applyBTagSF = 1; // 1 = EventWeigthing using SF and MC efficiencies.
 if (debug)cout <<"getting btag eff graphs..."<<endl;

   btageffFile->GetObject("TGraph/BTagEfficiency_e_b_TTJets", BTagGraph_B);
    btageffFile->GetObject("TGraph/BTagEfficiency_e_c_TTJets", BTagGraph_C);
    btageffFile->GetObject("TGraph/BTagEfficiency_e_l_TTJets_eta0-08", BTagGraph_L_eta008);
    btageffFile->GetObject("TGraph/BTagEfficiency_e_l_TTJets_eta08-16", BTagGraph_L_eta0816);
    btageffFile->GetObject("TGraph/BTagEfficiency_e_l_TTJets_eta16-24", BTagGraph_L_eta1624);

          if (applyBTagSF ==1){

 if (debug)cout <<"applying event weighting for Btagging...."<<endl;

         for ( seljet =0; seljet < selectedJets.size(); seljet++ ){
                if (selectedJets[seljet]->btag_combinedSecondaryVertexBJetTags() > workingpointvalue   ){
                    
                    selectedBJets.push_back(selectedJets[seljet]);
                }
                
                else{
                    selectedLightJets.push_back(selectedJets[seljet]);
                }

            }
	 ///start datacheck here...

	 if (dataSetName != "Data"){

   double p_tags_tagged_mc = 1.;
   double p_tags_untagged_mc = 1.;
   double p_tags_tagged_data = 1.;
   double p_tags_untagged_data = 1.;
   double p_mc = 1., p_data = 1.;
   int jet_flavor;
   float eff=1   ;
   float scaled_eff=1 ;
   float a_eff = 1;
   float sf_a_eff = 1;
   double LightJeteff;
   double JetPt, JetEta;
   double SF_tag =1.;
   double event_weight = 1.;
   for ( seljet =0; seljet < selectedJets.size(); seljet++ ){
     //for each jet, get the flavour, PT, efficiency and SF
      jet_flavor = selectedJets[seljet]->partonFlavour();
      JetPt = selectedJets[seljet]->Pt() ;
      JetEta = selectedJets[seljet]->Eta() ;
      if (JetPt > 800.) JetPt = 800;
      if (JetEta > 2.4) {JetEta = 2.4;
      }else if (JetEta < -2.4) {JetEta = -2.4;
      }

      if(fabs(jet_flavor) == 5 || fabs(jet_flavor) == 4  ){
      SF_tag =  bTool->getWeight(selectedJets[seljet]->Pt(),selectedJets[seljet]->Eta(),jet_flavor,dobTagEffShift );
		 
      //cout <<" b jet... "<<endl;
}
else{
  //  cout <<" light jet... "<<endl;
  SF_tag =  bTool->getWeight(selectedJets[seljet]->Pt(),selectedJets[seljet]->Eta(),jet_flavor,domisTagEffShift);
}


      //      cout <<"SFl mean "<<   bTool->getWeight(selectedJets[seljet]->Pt(),selectedJets[seljet]->Eta(),1,0)  <<endl;
      // cout <<"SFl min "<<   bTool->getWeight(selectedJets[seljet]->Pt(),selectedJets[seljet]->Eta(),1,-1)  <<endl;
      //cout <<"SFl max "<<   bTool->getWeight(selectedJets[seljet]->Pt(),selectedJets[seljet]->Eta(),1,1)  <<endl;
      //cout <<" "<< endl;


      if (debug) cout <<"got SFs"<<endl;

      if (fabs(jet_flavor) == 5  ) eff =  BTagGraph_B->Eval(JetPt,0,"");
      else  if( fabs(jet_flavor) == 4  ) eff =  BTagGraph_C->Eval(JetPt,0,"");
      else if(fabs(JetEta) < 0.8) eff =  BTagGraph_L_eta008->Eval(JetPt,0,"");
      else if(fabs(JetEta) >= 0.8  &&  fabs(JetEta) < 1.6 ) eff =  BTagGraph_L_eta0816->Eval(JetPt,0,"");
      else if(fabs(JetEta) >= 1.6 &&  fabs(JetEta) <= 2.4)  eff =  BTagGraph_L_eta1624->Eval(JetPt,0,"");

      a_eff = 1 - eff;
      scaled_eff = SF_tag*eff;
      sf_a_eff = 1 - scaled_eff;

      //  cout <<"eff  "<< eff <<endl;
      // cout <<"scaled_eff  "<< scaled_eff <<endl;
      //cout <<"sf_a_eff  "<< sf_a_eff <<endl;
      //cout <<"1 - eff  "<< a_eff <<endl;


      if(eff==1.)    cout <<"jet flavor "<<jet_flavor <<" jet pt "<< JetPt  <<"jet Eta = "<< JetEta  <<"  eff = " << eff  << " SF tag = "<<  SF_tag  << endl;
     if (selectedJets[seljet]->btag_combinedSecondaryVertexBJetTags() > workingpointvalue   ){
       p_tags_tagged_mc =  p_tags_tagged_mc*eff;
       p_tags_tagged_data = p_tags_tagged_data*(SF_tag*eff) ;

}
     else{

       p_tags_untagged_mc = (p_tags_untagged_mc)*(a_eff);
       p_tags_untagged_data = (p_tags_untagged_data)*(sf_a_eff) ;
}

}

   p_mc = (p_tags_tagged_mc)*(p_tags_untagged_mc);
   p_data = (p_tags_tagged_data)*( p_tags_untagged_data);


   if (p_mc != 0.)event_weight = (p_data)/(p_mc);
  
 histo1D["btag_weight"]->Fill(event_weight);

 if (debug)cout <<" filled weight histo"<<endl;

 if(event_weight < 0.){

   cout <<"NEGATIVE EVENT WEIGHT FROM BTAGGING!!!!! MUST BE FIXED!"<<endl;
       cout <<"p_tags_tagged_mc =   "<<  p_tags_tagged_mc <<endl;
       cout <<"p_tags_untagged_mc =   "<<  p_tags_untagged_mc <<endl;
      cout <<"p_tags_tagged_data =   "<<  p_tags_tagged_data <<endl;
      cout <<"p_tags_untagged_data =   "<<  p_tags_untagged_data <<endl;
      cout <<"p_mc = "<< p_mc <<endl;
      cout <<"p_data = "<< p_data <<endl;
      cout <<"event weight = = "<< event_weight <<endl;
      cout <<" "<<endl;
      cout <<" "<<endl;
      cout <<" "<<endl;
      cout <<" "<<endl;
      cout <<" "<<endl;
 }

 if (debug)cout <<"Btag Event weight applied"<<endl;
 

 if (dataSetName != "Data" && event_weight > 0. && event_weight < 10.) scaleFactor = scaleFactor*event_weight; 

}
	}
        else{ //not applying SF
            if (debug)cout <<"NOT APPLYING SF"<<endl;
            
            
            for ( seljet =0; seljet < selectedJets.size(); seljet++ ){
                if (selectedJets[seljet]->btag_combinedSecondaryVertexBJetTags() > workingpointvalue   ){
                    
                    selectedBJets.push_back(selectedJets[seljet]);
                }
                
                else{
                    selectedLightJets.push_back(selectedJets[seljet]);
                }
                
            }
            
        }
        

	HTRat = HT_leading/HT_lagging;

   double HTb = 0.;
        for ( Int_t sljt =0; sljt < selectedBJets.size(); sljt++ ){
            HTb += selectedBJets[sljt]->Pt();
        } 
      


	if (dataSetName != "Data"&&  selectedElectrons.size() ==1 ) {
   scaleFactor = scaleFactor*leptonTools->getElectronSF(selectedElectrons[0]->Eta(), selectedElectrons[0]->Pt(),leptonsyst );
 }


	if(isGoodPV && trigged){
	selecTableEl.Fill(d,0,scaleFactor);
	if (nEl==1){  selecTableEl.Fill(d,1,scaleFactor);
	  if(selectedLooseElectrons.size() == 1){selecTableEl.Fill(d,2,scaleFactor);
	    if(selectedLooseMuons.size() == 0){selecTableEl.Fill(d,3, scaleFactor);
	      if (selectedJets.size() >= 6) {  selecTableEl.Fill(d,4,scaleFactor) ;
		if (selectedBJets.size() >= 2){ selecTableEl.Fill(d,5,scaleFactor) ;
		  if (temp_HT >= 400.){ selecTableEl.Fill(d,6,scaleFactor) ;
		  
		  if (mets[0]->Et() >= 30.){ selecTableEl.Fill(d,7,scaleFactor) ;
		  }                  
}
                }
	      }
	    }
	  }
        }
      }
	/////////////////////////////////
	/////Applying baseline selection
	/////////////////////////////////

	  //Apply the lepton, btag and HT selections
        if  (  !( selectedJets.size() >= 6 && selectedBJets.size() >= 2 && selectedElectrons.size() == 1 && selectedLooseElectrons.size() ==1 && selectedLooseMuons.size() ==0 && temp_HT >= 400. && mets[0]->Et()> 30. )) continue;

	//	cout<<"event passed...."<<endl;

  
        vector<TLorentzVector*> selectedMuonTLV_JC;
        selectedMuonTLV_JC.push_back(selectedElectrons[0]);
        

        vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV, mcMuonsTLV;
        vector<TRootMCParticle*> mcParticlesMatching_;
        bool muPlusFromTop = false, muMinusFromTop = false;
        bool elPlusFromTop = false, elMinusFromTop = false;
        pair<unsigned int, unsigned int> leptonicBJet_, hadronicBJet_, hadronicWJet1_, hadronicWJet2_; //First index is the JET number, second one is the parton
        leptonicBJet_ = hadronicBJet_ = hadronicWJet1_ = hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
        
////////////////////////////////////////////////////////////////////////////////////
//// Getting Gen Event
////////////////////////////////////////////////////////////////////////////////////
        vector<TRootMCParticle*> mcParticles;

         if(dataSetName != "data" && dataSetName != "Data" && dataSetName != "Data"){
        vector<TRootMCParticle*> mcParticles;
        vector<TRootMCParticle*> mcTops;
             
        mcParticlesMatching_.clear();
        mcParticlesTLV.clear();
        selectedJetsTLV.clear();
                     
        int leptonPDG, muonPDG = 13, electronPDG = 11;
        leptonPDG = electronPDG;
             
        TRootGenEvent* genEvt = 0;
        genEvt = treeLoader.LoadGenEvent(ievt,false);
        sort(selectedJets.begin(),selectedJets.end(),HighestPt()); 
        treeLoader.LoadMCEvent(ievt, genEvt, 0, mcParticles,false);  
    
        if (debug) cout <<"size   "<< mcParticles.size()<<endl;

         
             //Pick out MCParticles of interest
             for (Int_t i=0; i<mcParticles.size(); i++){
                 if( mcParticles[i]->status() != 3) continue;
                 
                 if( (abs(mcParticles[i]->type())< 6||mcParticles[i]->type()==21)&&mcParticles[i]->status() ==3){
                     mcParticlesTLV.push_back(*mcParticles[i]);
                     mcParticlesMatching_.push_back(mcParticles[i]);
                 }
                 
                 // check if there is a mu +/- or a el +/-
                 if( mcParticles[i]->type() == leptonPDG && mcParticles[i]->motherType() == -24 && mcParticles[i]->grannyType() == -6 )
                 {
                     //if(muMinusFromTop) cerr<<"muMinusFromTop was already true"<<endl;
                     if(leptonPDG==muonPDG) muMinusFromTop = true;
                     else if(leptonPDG==electronPDG) elMinusFromTop = true;
                 }
                 if( mcParticles[i]->type() == -leptonPDG && mcParticles[i]->motherType() == 24 && mcParticles[i]->grannyType() == 6 )
                 {
                     //if(muPlusFromTop) cerr<<"muPlusFromTop was already true"<<endl;
                     if(leptonPDG==muonPDG) muPlusFromTop = true;
                     else if(leptonPDG==electronPDG) elPlusFromTop = true;
                 }
                 
             }
         
         }
        
        MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);

        if (debug)cout <<"event loop 8"<<endl;
        
        //Electrons
        MSPlot["NbOfIsolatedElectrons"]->Fill(selectedElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["ElectronPt"]->Fill(selectedElectrons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["ElectronEta"]->Fill(selectedElectrons[0]->Eta(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["ElectronPhi"]->Fill(selectedElectrons[0]->Phi(), datasets[d], true, Luminosity*scaleFactor);
                                    


    double EffectiveArea = 0.;
    double rho_ = event->kt6PFJets_rho();    
    double isocorr = 0;


    // HCP 2012 updated for electron conesize = 0.3, taken from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/UserCode/EGamma/EGammaAnalysisTools/interface/ElectronEffectiveArea.h?revision=1.4&view=markup

    if (fabs(selectedElectrons[0]->superClusterEta()) >= 0.0   && fabs(selectedElectrons[0]->superClusterEta()) < 1.0   ) EffectiveArea = 0.130;
    if (fabs(selectedElectrons[0]->superClusterEta()) >= 1.0   && fabs(selectedElectrons[0]->superClusterEta()) < 1.479 ) EffectiveArea = 0.137;
    if (fabs(selectedElectrons[0]->superClusterEta()) >= 1.479 && fabs(selectedElectrons[0]->superClusterEta()) < 2.0   ) EffectiveArea = 0.067;
    if (fabs(selectedElectrons[0]->superClusterEta()) >= 2.0   && fabs(selectedElectrons[0]->superClusterEta()) < 2.2   ) EffectiveArea = 0.089;
    if (fabs(selectedElectrons[0]->superClusterEta()) >= 2.2   && fabs(selectedElectrons[0]->superClusterEta()) < 2.3   ) EffectiveArea = 0.107;
    if (fabs(selectedElectrons[0]->superClusterEta()) >= 2.3   && fabs(selectedElectrons[0]->superClusterEta()) < 2.4   ) EffectiveArea = 0.110;
    if (fabs(selectedElectrons[0]->superClusterEta()) >= 2.4) EffectiveArea = 0.138;


    isocorr = rho_*EffectiveArea;


 float reliso = (selectedElectrons[0]->chargedHadronIso() + max(0.0 , selectedElectrons[0]->neutralHadronIso() + selectedElectrons[0]->photonIso() - isocorr) )/ selectedElectrons[0]->Pt();


        MSPlot["ElectronRelIsolation"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["ElectronMissingHits"]->Fill(selectedElectrons[0]->missingHits(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["ElectronMVATrigID"]->Fill(selectedElectrons[0]->mvaTrigId(), datasets[d], true, Luminosity*scaleFactor);


/////////////
////Jets
/////////////
	  MSPlot["NbOfSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["NbOfSelectedLightJets"]->Fill(selectedLightJets.size(), datasets[d], true, Luminosity*scaleFactor);
      MSPlot["NbOfSelectedBJets"]->Fill(selectedBJets.size(), datasets[d], true, Luminosity*scaleFactor);
	  //  MSPlot["RhoCorrection"]->Fill(event->kt6PFJetsPF2PAT_rho(), datasets[d], true, Luminosity*scaleFactor);
	  if (debug) cout <<"per jet plots.."<<endl;


       MSPlot["1stJetPt"]->Fill(selectedJets[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
       MSPlot["2ndJetPt"]->Fill(selectedJets[1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
       MSPlot["3rdJetPt"]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
       MSPlot["4thJetPt"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
       MSPlot["5thJetPt"]->Fill(selectedJets[4]->Pt(), datasets[d], true, Luminosity*scaleFactor);

      if (selectedJets.size()>=6) MSPlot["6thJetPt"]->Fill(selectedJets[5]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      if (selectedJets.size()>=7) MSPlot["7thJetPt"]->Fill(selectedJets[6]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	
    if (debug) cout <<"got muons and mets"<<endl;
        
/////////////////////////////////
/// Find indices of jets from Tops
////////////////////////////////
            
        for(unsigned int i=0; i<selectedJets.size(); i++) selectedJetsTLV.push_back(*selectedJets[i]);
        JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);
        
        vector< pair<unsigned int, unsigned int> > JetPartonPair;
        for(unsigned int i=0; i<mcParticlesTLV.size(); i++) //loop through mc particles and find matched jets
        {
            int matchedJetNumber = matching.getMatchForParton(i, 0);
            if(matchedJetNumber != -1)
                JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );// Jet index, MC Particle index
        }
        
       if (debug) cout <<"n sel jets  "<<selectedJets.size()  << "   n mc particles tlv : "<< mcParticlesTLV.size() << " jet parton pari size :   "<< JetPartonPair.size()<<"  "<< muPlusFromTop<<muMinusFromTop<<endl;
        
        for(unsigned int i=0; i<JetPartonPair.size(); i++)//looping through matched jet-parton pairs
        {
            unsigned int j = JetPartonPair[i].second;	  //get index of matched mc particle
        
            if( fabs(mcParticlesMatching_[j]->type()) < 5 )
            {
                if( ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == -24 && mcParticlesMatching_[j]->grannyType() == -6 )
                   || ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == 24 && mcParticlesMatching_[j]->grannyType() == 6 ) )
                {
                    if(hadronicWJet1_.first == 9999)
                    {
                        hadronicWJet1_ = JetPartonPair[i];
                       // MCPermutation[0] = JetPartonPair[i].first;
                    }
                    else if(hadronicWJet2_.first == 9999)
                    {
                        hadronicWJet2_ = JetPartonPair[i];
                        //MCPermutation[1] = JetPartonPair[i].first;
                    }
                    //else cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
                }
            }
            else if( fabs(mcParticlesMatching_[j]->type()) == 5 )
            {
               
                if(  ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == -6) || ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == 6 ) )
                {
                    hadronicBJet_ = JetPartonPair[i];
                    //MCPermutation[2] = JetPartonPair[i].first;
                }
                else if((muPlusFromTop && mcParticlesMatching_[j]->motherType() == 6) ||  ( muMinusFromTop &&mcParticlesMatching_[j]->motherType() == -6) )
                {
                    leptonicBJet_ = JetPartonPair[i];
                    //MCPermutation[3] = JetPartonPair[i].first;
                }
            }
        }
 
  //  cout <<"  "<<endl;
   if (debug) cout <<"Indices of matched jets are :  "<< hadronicBJet_.first<<"  "<< hadronicWJet1_.first  <<" " << hadronicWJet2_.first <<endl;
        
/////////////////////////////////
/// TMVA for mass reconstruction
////////////////////////////////

jetCombiner->ProcessEvent_SingleHadTop(datasets[d], mcParticles_flav, selectedJets, selectedMuonTLV_JC[0], genEvt_flav, scaleFactor);
if (debug) cout <<"Processing event with jetcombiner :  "<< endl;

 double MultiTopness;
      
if(!TrainMVA){
        MVAvals1 = jetCombiner->getMVAValue(MVAmethod, 1); // 1 means the highest MVA value
    
    MSPlot["MVA1TriJet"]->Fill(MVAvals1.first   ,  datasets[d], true, Luminosity*scaleFactor );

    selectedJets2ndPass.clear();
    selectedJets3rdPass.clear();
    MVASelJets1.clear();    


    //make vector of jets excluding thise selected by 1st pass of mass reco
    for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ ){
      if (seljet1 == MVAvals1.second[0] || seljet1 == MVAvals1.second[1] || seljet1 == MVAvals1.second[2]){ 
     MVASelJets1.push_back(selectedJets[seljet1]);
     continue;
      }
        selectedJets2ndPass.push_back(selectedJets[seljet1]);
         }
    
            jetCombiner->ProcessEvent_SingleHadTop(datasets[d], mcParticles_flav, selectedJets2ndPass, selectedMuonTLV_JC[0], genEvt_flav, scaleFactor);
            MVAvals2ndPass = jetCombiner->getMVAValue(MVAmethod, 1);


 MultiTopness = MVAvals2ndPass.first;

   MSPlot["MultiTopness"]->Fill(MultiTopness,  datasets[d], true, Luminosity*scaleFactor );

//check data-mc agreement of kin. reco. variables.
    float mindeltaR =100.;
    float mindeltaR_temp =100.;
    int wj1;
    int wj2;
    int bj1;
    
    //define the jets from W as the jet pair with smallest deltaR
    for (int m=0; m<MVASelJets1.size(); m++) {
        for (int n=0; n<MVASelJets1.size(); n++) {
            if(n==m) continue;
            TLorentzVector lj1  = *MVASelJets1[m];
            TLorentzVector lj2  = *MVASelJets1[n];
            mindeltaR_temp  = lj1.DeltaR(lj2);
            if (mindeltaR_temp < mindeltaR){
                mindeltaR = mindeltaR_temp;
                wj1 = m;
                wj2 = n;
            }
        }
    }
    // find the index of the jet not chosen as a W-jet
    for (unsigned int p=0; p<MVASelJets1.size(); p++) {
        if(p!=wj1 && p!=wj2) bj1 = p;
    }
    
    //now that putative b and W jets are chosen, calculate the six kin. variables.
    TLorentzVector Wh = *MVASelJets1[wj1]+*MVASelJets1[wj2];
    TLorentzVector Bh = *MVASelJets1[bj1];
    TLorentzVector Th = Wh+Bh;
    
    double TriJetMass = Th.M();

    double DiJetMass = Wh.M();
    //DeltaR
    float AngleThWh = fabs(Th.DeltaPhi(Wh));
    float AngleThBh = fabs(Th.DeltaPhi(Bh));

    float btag = MVASelJets1[bj1]->btag_combinedSecondaryVertexBJetTags();
    
    double PtRat = (( *MVASelJets1[0] + *MVASelJets1[1] + *MVASelJets1[2] ).Pt())/( MVASelJets1[0]->Pt() + MVASelJets1[1]->Pt() + MVASelJets1[2]->Pt() );
    
      MSPlot["MVA1TriJetMass"]->Fill(TriJetMass,  datasets[d], true, Luminosity*scaleFactor );
      MSPlot["MVA1DiJetMass"]->Fill(DiJetMass,  datasets[d], true, Luminosity*scaleFactor );
      MSPlot["MVA1BTag"]->Fill(btag,  datasets[d], true, Luminosity*scaleFactor );
      MSPlot["MVA1PtRat"]->Fill(PtRat,  datasets[d], true, Luminosity*scaleFactor );
      MSPlot["MVA1AnThWh"]->Fill(AngleThWh,  datasets[d], true, Luminosity*scaleFactor );
      MSPlot["MVA1AnThBh"]->Fill(AngleThBh,  datasets[d], true, Luminosity*scaleFactor );




    //make vector of jets excluding thise selected by 2nd pass of mass reco
    for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ ){
        if (seljet1 == MVAvals1.second[0] || seljet1 == MVAvals1.second[1] || seljet1 == MVAvals1.second[2] || seljet1 == MVAvals2ndPass.second[0] || seljet1 == MVAvals2ndPass.second[1] || seljet1 == MVAvals2ndPass.second[2] ) continue;
        selectedJets3rdPass.push_back(selectedJets[seljet1]);
    }
    

    if(selectedJets3rdPass.size()>= 3){
    jetCombiner->ProcessEvent_SingleHadTop(datasets[d], mcParticles_flav, selectedJets3rdPass, selectedMuonTLV_JC[0], genEvt_flav, scaleFactor);
    MVAvals3rdPass = jetCombiner->getMVAValue(MVAmethod, 1);

    MSPlot["MVA3rdPassTriJet"]->Fill(MVAvals3rdPass.first,  datasets[d], true, Luminosity*scaleFactor );
    }
    
    
    MSPlot["MVA2ndPassTriJet"]->Fill(MVAvals2ndPass.first,  datasets[d], true, Luminosity*scaleFactor );

    
           // cout <<"sel jets 1st pass = "<< selectedJets.size() << "sel jets 2nd pass =  " << selectedJets2ndPass.size() << endl;
            
         bestTopMass1 =( *selectedJets[MVAvals1.second[0]] + *selectedJets[MVAvals1.second[1]] + *selectedJets[MVAvals1.second[2]]).M();
         bestTopMass2ndPass =( *selectedJets[MVAvals2ndPass.second[0]] + *selectedJets[MVAvals2ndPass.second[1]] + *selectedJets[MVAvals2ndPass.second[2]]).M();

         bestTopPt =( *selectedJets[MVAvals1.second[0]] + *selectedJets[MVAvals1.second[1]] + *selectedJets[MVAvals1.second[2]]).Pt();

       // cout <<"Indices of best MVA jets are :  "<< MVAvals1.second[0] <<"  "<< MVAvals1.second[1]  <<" " << MVAvals1.second[2]<<endl;
            
         //   cout <<"MVA Mass 1 = "<< bestTopMass1 << " MVA Mass 2 = "<< bestTopMass2ndPass << endl;
    
   // cout <<"   "<<endl;
    
        MSPlot["MVA1TriJetMass"]->Fill(bestTopMass1,  datasets[d], true, Luminosity*scaleFactor );
        MSPlot["MVA2ndPassTriJetMass"]->Fill(bestTopMass2ndPass,  datasets[d], true, Luminosity*scaleFactor );

       if (debug)  cout <<"MVA Mass 1 = "<< bestTopMass1 << " MVA Mass 2 = "<< bestTopMass2 << endl;

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /////  Calculating how well the MVA jet selection is doing: Fraction of ttbar events            ////
        ////    where the jets selected by the TMVA massReco match the true jets from the hadronic top) ////
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        
        if(   ( hadronicBJet_.first == MVAvals1.second[0] || hadronicBJet_.first == MVAvals1.second[1] || hadronicBJet_.first == MVAvals1.second[2]   )  && ( hadronicWJet1_.first == MVAvals1.second[0] || hadronicWJet1_.first == MVAvals1.second[1] || hadronicWJet1_.first == MVAvals1.second[2]   )    && ( hadronicWJet2_.first == MVAvals1.second[0] || hadronicWJet2_.first == MVAvals1.second[1] || hadronicWJet2_.first == MVAvals1.second[2]   )      ){
            nMVASuccesses++;
        }
           
        }
        
        
////////////////////////////////////////////////////////////////////////
//          Plotting jet and event-level variables for rejecting ttbar + X
////////////////////////////////////////////////////////////////////////
        
        if (debug) cout <<"Caculating event level variables...  "<< endl;

        MEzCalculator NuPz;
        NuPz.SetMET(*mets[0]);
        NuPz.SetMuon(*selectedElectrons[0]);
        if (debug) cout <<"Created NuPz object "<< endl;
        
        for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ ){
            MSPlot["BdiscBJetCand_CSV"]->Fill(selectedJets[seljet1]->btag_combinedSecondaryVertexBJetTags(),datasets[d], true, Luminosity*scaleFactor);
            MSPlot["SelectedJetPt"]->Fill(selectedJets[seljet1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
            MSPlot["JetEta"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
            MSPlot["JetPhi"]->Fill(selectedJets[seljet1]->Phi() , datasets[d], true, Luminosity*scaleFactor);
            
            //Event-level variables
            jetpt = selectedJets[seljet1]->Pt();
            HT = HT + jetpt;
            H = H + selectedJets[seljet1]->P();
            sumpx = sumpx + selectedJets[seljet1]->Px();
            sumpy = sumpy + selectedJets[seljet1]->Py();
            sumpz = sumpz + selectedJets[seljet1]->Pz();
            sume = sume + selectedJets[seljet1]->E();  
            
            if(!TrainMVA){
            if(seljet1 == MVAvals1.second[0]  || seljet1 == MVAvals1.second[1] || seljet1 == MVAvals1.second[2]  ) continue;
            HTX = HTX + jetpt;
            HX = HX + selectedJets[seljet1]->P();
            sumpx_X = sumpx_X + selectedJets[seljet1]->Px();
            sumpy_X = sumpy_X + selectedJets[seljet1]->Py();
            sumpz_X = sumpz_X + selectedJets[seljet1]->Pz();
            sume_X = sume_X + selectedJets[seljet1]->E();
            }
        }
        
        sort(selectedJets2ndPass.begin(),selectedJets2ndPass.end(),HighestCVSBtag()); //order muons wrt dsicriminator

        if (debug) cout <<"Creating sumjets "<< endl;

        TRootJet sumjet (TLorentzVector (sumpx, sumpy, sumpz,sume )); //Object representing all the jets summed
        TRootJet sumjet_X (TLorentzVector (sumpx_X, sumpy_X, sumpz_X,sume_X )); //Object representing all the jets summed minus the hadronic system

        //now extend sumjet object to form event and reduced event
        //1. add muon
        sumpx_X = sumpx_X + selectedElectrons[0]->Px();
        sumpy_X = sumpy_X + selectedElectrons[0]->Py();
        sumpz_X = sumpz_X + selectedElectrons[0]->Pz();
        sume_X  = sume_X  + selectedElectrons[0]->E();
        
        //2. add MET
        sumpx_X = sumpx_X + mets[0]->Px();
        sumpy_X = sumpy_X + mets[0]->Py();
        sumpz_X = sumpz_X + NuPz.Calculate();
        sume_X = sume_X + mets[0]->E();
        
        
        //now extend sumjet object to form event
        //1. add muon
        sumpx = sumpx + selectedElectrons[0]->Px();
        sumpy = sumpy + selectedElectrons[0]->Py();
        sumpz = sumpz + selectedElectrons[0]->Pz();
        sume  = sume  + selectedElectrons[0]->E();
        
        //2. add MET
        sumpx = sumpx + mets[0]->Px();
        sumpy = sumpy + mets[0]->Py();
        sumpz = sumpz + NuPz.Calculate();
        sume = sume + mets[0]->E();
        double HT_top =1.;
        TLorentzVector bestMVATop;
        double ang_leptop = 0.;
        
        if(!TrainMVA){ bestMVATop =( *selectedJets[MVAvals1.second[0]] + *selectedJets[MVAvals1.second[1]] + *selectedJets[MVAvals1.second[2]]);
            HT_top = selectedJets[MVAvals1.second[0]]->Pt() + selectedJets[MVAvals1.second[1]]->Pt() + selectedJets[MVAvals1.second[2]]->Pt();
            ang_leptop = bestMVATop.DeltaPhi(*selectedElectrons[0]);
            
        }
        TRootJet Event (TLorentzVector (sumpx, sumpy, sumpz,sume )); //Object representing the whole event

        TRootJet Event_X (TLorentzVector (sumpx_X, sumpy_X, sumpz_X,sume_X )); //Object representing the event minus the hadronic top

        EventMass = Event.M();
        EventMassX = Event_X.M();
        SumJetMass = sumjet.M();
        SumJetMassX = sumjet_X.M();
        
        PTBalTopEventX = bestMVATop.Pt() - Event_X.Pt()  ;
        PTBalTopSumJetX = bestMVATop.Pt() - sumjet_X.Pt()  ;
        
        HTH = HT/H;
        HTXHX = HTX/HX;
        
        double HT_Asym = (HT_top - HT )/(HT_top + HT);
        
	//  MSPlot["HT_Asym"]->Fill(HT_Asym,datasets[d], true, Luminosity*scaleFactor);
	// MSPlot["DeltaPhi_LepTop"]->Fill(ang_leptop,datasets[d], true, Luminosity*scaleFactor);
        
         TLorentzVector mumet;
        
        //Form Lep W
         double   mumpx =   selectedElectrons[0]->Px() + mets[0]->Px();
         double   mumpy =   selectedElectrons[0]->Py() + mets[0]->Py();
         double   mumpz =   selectedElectrons[0]->Pz() + NuPz.Calculate();
         double   mumpe =   selectedElectrons[0]->E()  + mets[0]->E();
         mumet.SetPx(mumpx);
         mumet.SetPy(mumpy);
         mumet.SetPz(mumpz);
         mumet.SetE(mumpe);
        
        TLorentzVector muMETHadTop;
        TLorentzVector muMETB;
        TLorentzVector muMETBHadTop;
        
        if(!TrainMVA){
        muMETHadTop =( *selectedJets[MVAvals1.second[0]] + *selectedJets[MVAvals1.second[1]] + *selectedJets[MVAvals1.second[2]]  + mumet);
        muMETB =( *selectedJets2ndPass[0]  + mumet);
        muMETBHadTop =( *selectedJets[MVAvals1.second[0]] + *selectedJets[MVAvals1.second[1]] + *selectedJets[MVAvals1.second[2]]  + mumet + *selectedJets2ndPass[0] );
        }
        
        double deltaMTopMuMetB = muMETB.M() - bestMVATop.M() ;

        double deltaMTopMuMet = mumet.M() - bestMVATop.M() ;

        //cout <<"  mumetB mass  "<<  muMETB.M()   << "  MVA Tops Mass "<< bestMVATop.M()   <<  "  difference : "<< deltaMTopMuMetB <<endl;
        
        //cout<<" H  = "<< H  << " HT  = "<< HT <<" EventMass =   "<< EventMass<< " Event mass X  "<< EventMassX <<  "SumJet mass " <<  SumJetMass<< "SumJet mass X " << SumJetMassX  <<endl;
        
        if (debug) cout <<"arrays  "<<endl;
        
        MHT = sumjet.Pt();
        MHTSig = MHT/sqrt(HT);
        STJet = HT + mets[0]->Et() + selectedElectrons[0]->Pt();
        EventMass = sumjet.M();
    
        sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag()); //order muons wrt dsicriminator
        double bTagRatio13 =0.;
        double bTagRatio14 =0.;

        if(selectedJets[2]->btag_combinedSecondaryVertexBJetTags() != 0 && selectedJets[3]->btag_combinedSecondaryVertexBJetTags() != 0){
        bTagRatio13 = selectedJets[0]->btag_combinedSecondaryVertexBJetTags()/selectedJets[2]->btag_combinedSecondaryVertexBJetTags();
        bTagRatio14 = selectedJets[0]->btag_combinedSecondaryVertexBJetTags()/selectedJets[3]->btag_combinedSecondaryVertexBJetTags();
        }
            
        MSPlot["bTagRatio13"]->Fill(bTagRatio13,datasets[d], true, Luminosity*scaleFactor);
        MSPlot["bTagRatio14"]->Fill(bTagRatio14,datasets[d], true, Luminosity*scaleFactor);

        sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order muons wrt dsicriminator
        
        if (debug) cout <<"jets sorted...  "<<endl;

        
        if(trainEventMVA){
        
        if(dataSetName=="NP_overlay_TTTT"){
        Eventtrainer_->FillWeight("S","Weight", scaleFactor);

       // Eventtrainer_->Fill("S","H", H);
       // Eventtrainer_->Fill("S","HT", HT);
        Eventtrainer_->Fill("S","HTX", HTX);
        Eventtrainer_->Fill("S","HTH", HTH);  
        Eventtrainer_->Fill("S","HTb", HTb);
        Eventtrainer_->Fill("S","HTRat", HTRat);
	// Eventtrainer_->Fill("S","HTXHX", HTXHX);
       // Eventtrainer_->Fill("S","EventMass", EventMass);
       // Eventtrainer_->Fill("S","EventMassX", EventMassX);
       // Eventtrainer_->Fill("S","SumJetMass", SumJetMass);
        Eventtrainer_->Fill("S","SumJetMassX", SumJetMassX);
       // Eventtrainer_->Fill("S","PTBalTopEventX", PTBalTopEventX);
       // Eventtrainer_->Fill("S","PTBalTopSumJetX", PTBalTopSumJetX);
        Eventtrainer_->Fill("S","MultiTopness", MultiTopness);
            
	//  Eventtrainer_->Fill("S","bTagRatio13", bTagRatio13);
	// Eventtrainer_->Fill("S","bTagRatio14", bTagRatio14);
	//  Eventtrainer_->Fill("S","nJets",selectedJets.size()  );
       Eventtrainer_->Fill("S","nTags",selectedBJets.size()  );
       // Eventtrainer_->Fill("S","Jet1Pt",selectedJets[0]->Pt());
       // Eventtrainer_->Fill("S","Jet2Pt",selectedJets[1]->Pt());
       // Eventtrainer_->Fill("S","Jet3Pt",selectedJets[2]->Pt());
       // Eventtrainer_->Fill("S","Jet4Pt",selectedJets[3]->Pt());
        Eventtrainer_->Fill("S","Jet5Pt",selectedJets[4]->Pt());
        Eventtrainer_->Fill("S","Jet6Pt",selectedJets[5]->Pt());
        }
        
        if(dataSetName=="TTJets"){
         //   cout <<"jets sorted.  ttjets..  "<<endl;

              Eventtrainer_->FillWeight("B","Weight", scaleFactor);

           // Eventtrainer_->Fill("B","H", H);
           // cout <<"jets sorted.  ttjets.. cc0 "<<endl;

           // Eventtrainer_->Fill("B","HT", HT);
            Eventtrainer_->Fill("B","HTX", HTX);
            Eventtrainer_->Fill("B","HTH", HTH);
            Eventtrainer_->Fill("B","HTb", HTb);
            Eventtrainer_->Fill("B","HTRat", HTRat);
	    //            Eventtrainer_->Fill("B","HTXHX", HTXHX);
        //    Eventtrainer_->Fill("B","EventMass", EventMass);
	    //   Eventtrainer_->Fill("B","EventMassX", EventMassX);
          //  Eventtrainer_->Fill("B","SumJetMass", SumJetMass);
            Eventtrainer_->Fill("B","SumJetMassX", SumJetMassX);
         //   Eventtrainer_->Fill("B","PTBalTopEventX", PTBalTopEventX);
         //   Eventtrainer_->Fill("B","PTBalTopSumJetX", PTBalTopSumJetX);
	    // Eventtrainer_->Fill("B","bestTopMass2ndPass", bestTopMass2ndPass);
            
            Eventtrainer_->Fill("B","MultiTopness", MultiTopness);
	    // Eventtrainer_->Fill("B","bTagRatio13", bTagRatio13);
	    // Eventtrainer_->Fill("B","bTagRatio14", bTagRatio14);
	    //  Eventtrainer_->Fill("B","nJets", selectedJets.size() );
       Eventtrainer_->Fill("B","nTags", selectedBJets.size() ); 
        //    Eventtrainer_->Fill("B","Jet1Pt", selectedJets[0]->Pt() );
        //    Eventtrainer_->Fill("B","Jet2Pt", selectedJets[1]->Pt() );
        //    Eventtrainer_->Fill("B","Jet3Pt", selectedJets[2]->Pt() );
        //    Eventtrainer_->Fill("B","Jet4Pt", selectedJets[3]->Pt() );
            Eventtrainer_->Fill("B","Jet5Pt", selectedJets[4]->Pt() );
            Eventtrainer_->Fill("B","Jet6Pt", selectedJets[5]->Pt() );
        
            
        }
            
        }
        
        else if (computeEventMVA){
            if (debug) cout <<"filling computer...."<<H <<endl;

            if (Eventcomputer_ == 0) cout <<"null computer...."<<H <<endl;

            
       // Eventcomputer_->FillVar("H", H);
            
        if (debug) cout <<"filling computer...."<<endl;

	//        Eventcomputer_->FillVar("HT", HT);
        Eventcomputer_->FillVar("HTX", HTX);
        Eventcomputer_->FillVar("HTH", HTH);
        Eventcomputer_->FillVar("HTb", HTb);
        Eventcomputer_->FillVar("HTRat", HTRat);
        //Eventcomputer_->FillVar("EventMass", EventMass);
	//        Eventcomputer_->FillVar("EventMassX", EventMassX);
        //Eventcomputer_->FillVar("SumJetMass", SumJetMass);
        Eventcomputer_->FillVar("SumJetMassX", SumJetMassX);
        //Eventcomputer_->FillVar("PTBalTopEventX", PTBalTopEventX);
        //Eventcomputer_->FillVar("PTBalTopSumJetX", PTBalTopSumJetX);
        Eventcomputer_->FillVar("MultiTopness", MultiTopness);
        Eventcomputer_->FillVar("nTags", selectedBJets.size());
	//        Eventcomputer_->FillVar("bTagRatio14", bTagRatio14);
	//  Eventcomputer_->FillVar("nJets", selectedJets.size() );
            
            if (debug) cout <<"pre jets file///."<<endl;

      //  Eventcomputer_->FillVar("Jet1Pt", selectedJets[0]->Pt() );
      //  Eventcomputer_->FillVar("Jet2Pt", selectedJets[1]->Pt() );
      //  Eventcomputer_->FillVar("Jet3Pt", selectedJets[2]->Pt() );
      //  Eventcomputer_->FillVar("Jet4Pt", selectedJets[3]->Pt() );
        Eventcomputer_->FillVar("Jet5Pt", selectedJets[4]->Pt() );
        Eventcomputer_->FillVar("Jet6Pt", selectedJets[5]->Pt() );
            if (debug) cout <<"filled inside loop...."<<endl;

            
            }
        if (debug) cout <<"computer filled...."<<endl;
        
        double BDTscore;
        
        if(computeEventMVA){  std::map<std::string,Float_t> MVAVals = Eventcomputer_->GetMVAValues();
        
        for (std::map<std::string,Float_t>::const_iterator it = MVAVals.begin(); it != MVAVals.end(); ++it){
        
          //  cout <<"MVA Method : "<< it->first    <<" Score : "<< it->second <<endl;
            BDTscore = it->second;

        }
        
    }
       // cout <<"  "<<endl;
                
        MSPlot["MVA"]->Fill(BDTscore, datasets[d], true, Luminosity*scaleFactor);
        
     
       
	if((BDTscore > 0.2) && (dataSetName=="Data")){
	  eventlist <<currentRun  << " " << event->lumiBlockId() <<" " <<event->eventId() << " "  << BDTscore <<endl;        

	}


	if((BDTscore > 0.05)){

	  MSPlot["NbOfSelectedJets_BDTCut"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);                
    MSPlot["NbOfSelectedBJets_BDTCut"]->Fill(selectedBJets.size(), datasets[d], true, Luminosity*scaleFactor);             
    MSPlot["HTX_BDTCut"]->Fill(HTX, datasets[d], true, Luminosity*scaleFactor); 
    MSPlot["HTH_BDTCut"]->Fill(HTH, datasets[d], true, Luminosity*scaleFactor); 
    MSPlot["MultiTopness_BDTCut"]->Fill(MultiTopness, datasets[d], true, Luminosity*scaleFactor); 
    MSPlot["HTb_SelectedJets_BDTCut"]->Fill(HTb, datasets[d], true, Luminosity*scaleFactor); 


}



       if (debug) cout <<"filling event level histos  "<<endl;

        MSPlot["EventMass"]->Fill(EventMass, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["EventMassX"]->Fill(EventMassX, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["SumJetMass"]->Fill(SumJetMass, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["SumJetMassX"]->Fill(SumJetMassX, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["HTX"]->Fill(HTX, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["H"]->Fill(H, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["HX"]->Fill(HX, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["HTXHX"]->Fill(HTXHX, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["HTH"]->Fill(HTH, datasets[d], true, Luminosity*scaleFactor);
        
        MSPlot["HT_SelectedJets"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["HTb_SelectedJets"]->Fill(HTb, datasets[d], true, Luminosity*scaleFactor);

        MSPlot["HTRat"]->Fill(HTRat, datasets[d], true, Luminosity*scaleFactor);

        MSPlot["MET"]->Fill(mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor);
	 
        if (debug) cout <<"filling event level histos  "<<endl;  
        
    
//////////////////////////////////////////////////////////////////////////////////
/////          Filling NJet Vs NBJet arrays for                             //////
/////                                                                       //////
///// 1. HT, 2. STLep, 3.STJet, 4.HT(w/o ttbar), 5.InvMass (w/o) ttbar      //////
//////////////////////////////////////////////////////////////////////////////////

      for (Int_t c = minNJets; c<= maxNJets; c++){
          string NJets_str = static_cast<ostringstream*>( &(ostringstream() << c) )->str();

          string MVA_Name = "MVA"+NJets_str+"Jets" ;

          
          
        if(c<=7){
            if(selectedJets.size() == c  ) {
        MSPlot[MVA_Name.c_str()]->Fill(BDTscore,  datasets[d], true, Luminosity*scaleFactor );

 }
                        }
        else if ( c > 7) {
	  if( selectedJets.size() >= c  ) {

                MSPlot[MVA_Name.c_str()]->Fill(BDTscore,  datasets[d], true, Luminosity*scaleFactor );

        }
}
}



  for (Int_t c = minNBJets; c <= maxNBJets; c++){

       string NBJets_str = static_cast<ostringstream*>( &(ostringstream() << c) )->str();
       string MVA_Name = "MVA"+NBJets_str+"Tags" ;
          
        if(c<=2){
            if(selectedBJets.size() == c  ) {
        MSPlot[MVA_Name.c_str()]->Fill(BDTscore,  datasets[d], true, Luminosity*scaleFactor );

 }
                        }
        else if ( c > 2) {
	  if( selectedBJets.size() >= c  ) {
	
                MSPlot[MVA_Name.c_str()]->Fill(BDTscore,  datasets[d], true, Luminosity*scaleFactor );

        }
}

}


    }//loop on events
    
    if (debug)cout <<"N BBar = = " << nBBBar <<"  N CCBar = = " << nCCBar <<"  N LLBar = =  " << nLLBar << endl;
      

    if(jetTools) delete jetTools;
    if(L1JetCorPar) delete L1JetCorPar;
    if(L2JetCorPar) delete L2JetCorPar;
    if(L3JetCorPar) delete L3JetCorPar;
    if(L2L3ResJetCorPar) delete L2L3ResJetCorPar;



      //important: free memory
      treeLoader.UnLoadDataset();
      
      double nMVASuccessesd = nMVASuccesses;
      double nMatchedEventsd = nMatchedEvents;
      
      cout <<"Efficiency of MVA jet selection = = "<<  nMVASuccessesd/nMatchedEventsd   <<endl;
      
  } //loop on datasets
   

eventlist.close();

	if (leptonTools)	delete leptonTools;

    if(trainEventMVA) Eventtrainer_->TrainMVA("Block","",0,0,"",0,0,"test_El");
    
  //Once everything is filled ...
  cout << " We ran over all the data ;-)" << endl;
  
  ///////////////////
  // Writing
  //////////////////
  cout << " - Writing outputs to the files ..." << endl;

  //Selection tables
  if(Muon){
	//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
      
  //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
  }
    else if(Electron){
	//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
	selecTableEl.TableCalculator(  false, true, true, true, true);
   //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
    selecTableEl.Write(  "FourTop"+postfix+"Table_El.tex",  false,true,true,true,false,false,true);

  }
    
    fout->cd();
   TFile *foutmva = new TFile ("foutMVA.root","RECREATE");
   //foutmva->cd();
    cout <<" after cd .."<<endl;
    
    //Output ROOT file
    string mvarootFileName ("MVA"+postfix+channelpostfix+".root");
  //  TFile * foutmva;
  //  TFile* foutmva = new TFile(mvarootFileName.c_str(), "RECREATE");
    //foutmva->cd();
    
   // TH1F * h1 = new TH1F();
   // h1->Write();
    //foutmva->Write();
  //  foutmva->Close();
    
    string pathPNGJetCombi = pathPNG + "JetCombination/";
    mkdir(pathPNGJetCombi.c_str(),0777);
    
   // if (foutmva) cout <<"fout "<<endl;
    //cout <<"pointer to file -->"<< foutmva << endl;
    //foutmva->Print();
    
    if(TrainMVA)jetCombiner->Write(foutmva, true, pathPNGJetCombi.c_str());
   
//   jetCombiner->Write(fout, true, pathPNGJetCombi.c_str());
    
//   TFile* sysFile = new TFile("SystematicShapes.root","RECREATE");

  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        string name = it->first;
        
        MultiSamplePlot *temp = it->second;
        TH1F *tempHisto_data;
        TH1F *tempHisto_TTTT;


      //Option to write ROOT files containing histograms with systematics varied +/-
      //TFile * tempErrorFile;
        if(doScaleShift == 1){
	  cout <<"    "<<endl;
	  cout <<"    "<<endl;
	  cout <<"    "<<endl;
	  cout <<"    "<<endl;
	  cout <<"#########################################"<<endl;
          cout <<"Scaling sys down!!!!!!!!!"<<endl;
	  cout <<"#########################################"<<endl;
	  cout <<"    "<<endl;
	  cout <<"    "<<endl;
	  cout <<"    "<<endl;
	  cout <<"    "<<endl;

            string filename = "ScaleFilesEl/Error_"+name+".root";
            TFile* tempErrorFile = new TFile(filename.c_str(),"RECREATE");
            TH1F * tempHisto = temp->getTH1F("TTJets_ScaleDown");

            tempHisto->Write("Down");
            tempErrorFile->Write();
            tempErrorFile->Close();
	    if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	     string hist_name =  name + "_Scale_Down";
             tempHisto->Write(hist_name.c_str());
	    sysFile->Close();
	    }
	    delete tempErrorFile;
        }else if(doScaleShift == 2){
            cout <<"Scaling sys up"<<endl;
            string filename = "ScaleFilesEl/Error_"+name+".root";
            TFile* tempErrorFile = new TFile(filename.c_str(),"UPDATE");
            TH1F * tempHisto = temp->getTH1F("TTJets_ScaleUp");
            tempHisto->Write("Up");
            tempErrorFile->Write();
            tempErrorFile->Close();
	     if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	     string hist_name =  name + "_Scale_Up";
             tempHisto->Write(hist_name.c_str());
	    sysFile->Close();
	    }

	    delete tempErrorFile;
        }
        else if(doScaleShift == 0){
            cout <<"Scaling sys nominal"<<endl;
        }
        

	     //Option to write ROOT files containing histograms with systematics varied +/-
      //TFile * tempErrorFile;
        if(doMatchingShift == 1){
	  cout <<"    "<<endl;
	  cout <<"    "<<endl;
	  cout <<"    "<<endl;
	  cout <<"    "<<endl;
	  cout <<"#########################################"<<endl;
          cout <<"Matching sys down!!!!!!!!!"<<endl;
	  cout <<"#########################################"<<endl;
	  cout <<"    "<<endl;
	  cout <<"    "<<endl;
	  cout <<"    "<<endl;
	  cout <<"    "<<endl;

            string filename = "MatchingFilesEl/Error_"+name+".root";
            TFile* tempErrorFile = new TFile(filename.c_str(),"RECREATE");
            TH1F * tempHisto = temp->getTH1F("TTJets_MatchingDown");
            tempHisto->Write("Down");
            tempErrorFile->Write();
            tempErrorFile->Close();
	    if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	     string hist_name =  name + "_Matching_Down";
             tempHisto->Write(hist_name.c_str());
	    sysFile->Close();
	    }
	    delete tempErrorFile;
        }else if(doMatchingShift == 2){
            cout <<"Matching sys up"<<endl;
            string filename = "MatchingFilesEl/Error_"+name+".root";
            TFile* tempErrorFile = new TFile(filename.c_str(),"UPDATE");
            TH1F * tempHisto = temp->getTH1F("TTJets_MatchingUp");
            tempHisto->Write("Up");
            tempErrorFile->Write();
            tempErrorFile->Close();
	    if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	     string hist_name =  name + "_Matching_Up";
             tempHisto->Write(hist_name.c_str());
	    sysFile->Close();
	    }

	    delete tempErrorFile;
        }
	else if (doMatchingShift == 0){

 cout <<"Matching sys nominal"<<endl;
}

	if (dottbbShift == 1){
	  if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TH1F * tempHisto = temp->getTH1F("TTJets");
             TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	     string hist_name =  name + "_ttbb_Down";
             tempHisto->Write(hist_name.c_str());
	    sysFile->Close();
 }
}
	else if(dottbbShift == 2){
	  if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TH1F * tempHisto = temp->getTH1F("TTJets");
             TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	     string hist_name =  name + "_ttbb_Up";
             tempHisto->Write(hist_name.c_str());
	    sysFile->Close();
 }
}
	else if (dottbbShift == 0) {
         cout <<"ttbb sys nominal"<<endl;
}



        if (doLeptonSFShift == 1){
	  if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
	    TH1F * tempHisto = temp->getTH1F("TTJets");
	    TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	    string hist_name =  name + "_LeptonSF_Down";
	    tempHisto->Write(hist_name.c_str());
	   sysFile->Close();
	  }
	}
        else if(doLeptonSFShift == 2){
	  if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
	    TH1F * tempHisto = temp->getTH1F("TTJets");
	    TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	    string hist_name =  name + "_LeptonSF_Up";
	    tempHisto->Write(hist_name.c_str());
	   sysFile->Close();
	  }
	}
        else if (doLeptonSFShift == 0) {
	  cout <<"leptonSF sys nominal"<<endl;
	}


	if (doJERShift == 1){
	  if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
            TH1F * tempHisto = temp->getTH1F("TTJets");
            TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	    string hist_name =  name + "_JER_Down";
	    tempHisto->Write(hist_name.c_str());
	   sysFile->Close();
          }
	}
        else if(doJERShift == 2){
	  if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
            TH1F * tempHisto = temp->getTH1F("TTJets");
            TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	    string hist_name =  name + "_JER_Up";
	    tempHisto->Write(hist_name.c_str());
	   sysFile->Close();
          }
        }
        else if (doJERShift == 0) {
          cout <<"JER sys nominal"<<endl;
        }


	if (doPUShift == 1){
	  if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
            TH1F * tempHisto = temp->getTH1F("TTJets");
            TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	    string hist_name =  name + "_PU_Down";
	    tempHisto->Write(hist_name.c_str());
	   sysFile->Close();
          }
	}
        else if(doPUShift == 2){
	  if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
            TH1F * tempHisto = temp->getTH1F("TTJets");
            TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	    string hist_name =  name + "_PU_Up";
	    tempHisto->Write(hist_name.c_str());
	   sysFile->Close();
          }
        }
        else if (doPUShift == 0) {
          cout <<"PU sys nominal"<<endl;
        }



	if (domisTagEffShift == -1){
	  if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TH1F * tempHisto = temp->getTH1F("TTJets");
             TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	     string hist_name =  name + "_misTag_Down";
             tempHisto->Write(hist_name.c_str());
	    sysFile->Close();
 }
}
	else if(domisTagEffShift == 1){
	  if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TH1F * tempHisto = temp->getTH1F("TTJets");
             TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	     string hist_name =  name + "_misTag_Up";
             tempHisto->Write(hist_name.c_str());
	    sysFile->Close(); 
}
}
	else if (domisTagEffShift == 0) {
         cout <<"misTag sys nominal"<<endl;
}




	if (dobTagEffShift == -1){
	  if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TH1F * tempHisto = temp->getTH1F("TTJets");
             TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	     string hist_name =  name + "_bTag_Down";
             tempHisto->Write(hist_name.c_str());
	    sysFile->Close();
 }
}
	else if(dobTagEffShift == 1){
	  if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TH1F * tempHisto = temp->getTH1F("TTJets");
             TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	     string hist_name =  name + "_bTag_Up";
             tempHisto->Write(hist_name.c_str());
	    sysFile->Close();
 }
}
	else if (dobTagEffShift == 0) {
         cout <<"bTag sys nominal"<<endl;
}

           if (doJESShift == 1){
            string filename = "ErrorBandsEl/Error_"+name+".root";
            TFile* tempErrorFile = new TFile(filename.c_str(),"RECREATE");
            TH1F * tempHisto = temp->getTH1F("TTJets");
            tempHisto->Write("Minus");
            tempErrorFile->Write();
            tempErrorFile->Close();

	    if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	     string hist_name =  name + "_JES_Down";
             tempHisto->Write(hist_name.c_str());
	    sysFile->Close();
	    }

	    delete tempErrorFile;
        }
        else if  (doJESShift ==2){
            
            string filename = "ErrorBandsEl/Error_"+name+".root";
            TFile* tempErrorFile = new TFile(filename.c_str(),"UPDATE");
            TH1F * tempHisto = temp->getTH1F("TTJets");
            tempHisto->Write("Plus");
            tempErrorFile->Write();
            tempErrorFile->Close();

	    if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TFile* sysFile = new TFile("SystematicShapes_El.root","UPDATE");
	     string hist_name =  name + "_JES_Up";
             tempHisto->Write(hist_name.c_str());
	    sysFile->Close();
	    }
	    delete tempErrorFile;

            cout <<"JES sys down"<<endl;
            
        }
        else if  (doJESShift ==0){
            cout <<"JES sys off "<<endl;
      
   // temp->addText("CMS preliminary");
	    temp->Draw_wSysUnc(false,"ScaleFilesEl" , name, true, true, false, false, false,100.,true, false, false, true); // merge TT/QCD/W/Z/ST/
	//temp->Draw(false, name, true, true, true, true, true,100.,false); // merge TT/QCD/W/Z/ST/

            
            //Draw(bool addRandomPseudoData = false, string label = string("CMSPlot"), bool mergeTT = false, bool mergeQCD = false, bool mergeW = false, bool mergeZ = false, bool mergeST = false, int scaleNPSignal = 1, bool addRatio = false, bool mergeVV = false, bool mergeTTV = false);
      cout <<" Writing plot:  name = "<< name<<"  path = "<< pathPNG  <<endl;
      temp->Write(fout, name, true, pathPNG, "pdf");
    }
  }

  TDirectory* th1dir = fout->mkdir("Histos1D");
  th1dir->cd();
  for(map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {

	TH1F *temp = it->second;
	temp->Write();
	//TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
	//tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }
  TDirectory* th2dir = fout->mkdir("Histos2D");
  th2dir->cd();
   for(map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
  {

	TH2F *temp = it->second;
	temp->Write();
	//TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
	//tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }
    
  delete fout;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
}

int Factorial(int N = 1)
{
    int fact = 1;
    for( int i=1; i<=N; i++ )
        fact = fact * i;  // OR fact *= i;
    return fact;
}


//To cout the Px, Py, Pz, E and Pt of objects
void coutObjectsFourVector(vector < TRootMuon* > init_muons, vector < TRootElectron* > init_electrons, vector < TRootJet* > init_jets, vector < TRootMET* > mets, string Comment)
{
     cout<<Comment<<endl;
     
     for(unsigned int k=0; k<init_muons.size(); k++)
     {
	   cout<<" init_muons["<<k<<"] -> Px() = "<<init_muons[k]->Px()<<endl;
	   cout<<"              -> Py() = "<<init_muons[k]->Py()<<endl;
	   cout<<"              -> Pz() = "<<init_muons[k]->Pz()<<endl;
	   cout<<"                -> Pt() = "<<init_muons[k]->Pt()<<endl;
	   cout<<"              -> E() = "<<init_muons[k]->E()<<endl;   
     }
     for(unsigned int k=0; k<init_electrons.size(); k++)
     {
	   cout<<" init_electrons["<<k<<"] -> Px() = "<<init_electrons[k]->Px()<<endl;
	   cout<<"              -> Py() = "<<init_electrons[k]->Py()<<endl;
	   cout<<"              -> Pz() = "<<init_electrons[k]->Pz()<<endl;
	   cout<<"                -> Pt() = "<<init_electrons[k]->Pt()<<endl;
	   cout<<"              -> E() = "<<init_electrons[k]->E()<<endl;   
     }         
     for(unsigned int k=0; k<init_jets.size(); k++) //init_jets.size()
     {
	   cout<<" init_jets["<<k<<"] -> Px() = "<<init_jets[k]->Px()<<endl;
	   cout<<"              -> Py() = "<<init_jets[k]->Py()<<endl;
	   cout<<"              -> Pz() = "<<init_jets[k]->Pz()<<endl;
	   cout<<"                -> Pt() = "<<init_jets[k]->Pt()<<endl;
	   cout<<"              -> E() = "<<init_jets[k]->E()<<endl;	   
     }
     for(unsigned int k=0; k<mets.size(); k++)
     {
           cout<<" mets["<<k<<"] -> Px() = "<<mets[k]->Px()<<endl;
           cout<<"         ->  Py() = "<<mets[k]->Py()<<endl;
	   cout<<"         ->  Pz() = "<<mets[k]->Pz()<<endl;
	   cout<<"              -> Pt() = "<<mets[k]->Pt()<<endl;
	   cout<<"         ->  E() = "<<mets[k]->E()<<endl;
	   cout<<"              -> Et() = "<<mets[k]->Et()<<endl;
     }
};
