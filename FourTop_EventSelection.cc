//////////////////////////////////////////////////////////////////////////////
////         Analysis code for search for Four Top Production.                  ////
////////////////////////////////////////////////////////////////////////////////////

// ttbar @ NLO : 
//all-had ->234 * .46 = 107.64
//semi-lep ->234 *.45 = 105.3
//di-lep-> 234* .09 = 21.06

//ttbar @ NNLO: 
//all-had -> 245.8 * .46 = 113.068
//semi-lep-> 245.8 * .45 = 110.61
//di-lep ->  245.8 * .09 = 22.122

#define _USE_MATH_DEFINES
#include "TStyle.h"
#include "TPaveText.h"
#include "TTree.h"
#include "TNtuple.h"
#include <ctime>

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "TRandom3.h"
#include "TRandom.h"
#include "TProfile.h"
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
#include "TopTreeAnalysisBase/Tools/interface/AnalysisEnvironmentLoader.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectorParameters.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/JetCorrectionUncertainty.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MakeBinning.h"
#include "TopTreeAnalysisBase/MCInformation/interface/LumiReWeighting.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MEzCalculator.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"

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

bool split_ttbar = false;

pair<float, vector<unsigned int> > MVAvals1;
pair<float, vector<unsigned int> > MVAvals2;
pair<float, vector<unsigned int> > MVAvals2ndPass;
pair<float, vector<unsigned int> > MVAvals3rdPass;

int nMVASuccesses=0;
int nMatchedEvents=0;

/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;
map<string,TProfile*> histoProfile;

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

  int passed = 0;
  int ndefs =0;

  //  TRandom3* rand = new TRandom3();

 // BTagWeightTools * bTool = new BTagWeightTools("SFb-pt_payload_Moriond13.txt", "CSVM") ;
   
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

 int dottbbShift = 2;  //0: off (except nominal PU reweighting) 1: minus 2: plus
  cout << "dottbbShift: " << dottbbShift << endl;

 int doLeptonSFShift = 0;
 cout << "doLeptonSFShift: " << doLeptonSFShift << endl;

 string leptonsyst = "Nominal";

	if (doLeptonSFShift==1){

	  leptonsyst = "Minus";
	}else if(doLeptonSFShift==2){
	  leptonsyst = "Plus";
}


//bool to control whether to look at SIGNAL (HT > 400 GeV) or CONTROL (HT < 400 GeV) regions
   bool sig_region = true;

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
	  if(dobTagEffShift == 1)
		scalefactorbtageff = 0.85;
	  if(dobTagEffShift == 2)
		scalefactorbtageff = 1.03;
		
	  if(domisTagEffShift == 0)
		mistagfactor = 1.21;
	  if(domisTagEffShift == 1)
		mistagfactor = 1.04;
	  if(domisTagEffShift == 2)
		mistagfactor = 1.38;
  }

  float workingpointvalue = 0.679; //working points updated to 2012 BTV-POG recommendations.
 
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
  //  setGregStyle();
  //setMyStyle();

  string postfix = "_EventSelection_synclana_noNJets"; // to relabel the names of the output file

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
 if(dottbbShift == 1)
    postfix= postfix+"_ttbbMinus";
  if(dottbbShift == 2)
    postfix= postfix+"_ttbbPlus";
 if(doLeptonSFShift == 1)
    postfix= postfix+"_leptonSFMinus";
  if(doLeptonSFShift == 2)
    postfix= postfix+"_leptonSFPlus";



  ///////////////////////////////////////
  // Configuration
  ///////////////////////////////////////

  string channelpostfix = "";
  string xmlFileName = "";

  bool Muon = true; // use Muon channel?
  bool Electron = false;

  if(Muon){
	cout << " --> Using the Muon channel..." << endl;
	channelpostfix = "_Mu";
	//xmlFileName = "../config/myTopFCNCconfig_Muon.xml";
	//    xmlFileName = "config/test_fullsamples.xml";

  }


  bool TrainMVA = false; // If false, the previously trained MVA will be used to calculate stuff                
  bool trainEventMVA = false; // If false, the previously trained MVA will be used to calculate stuff             
  bool computeEventMVA = true; 


 if (doScaleShift==1){
  xmlFileName = "config/test_ScaleDown.xml";
}
  else if (doScaleShift==2){ 
 xmlFileName = "config/test_ScaleUp.xml";
  }
  else  if (doMatchingShift==1){
  xmlFileName = "config/test_MatchingDown.xml";
  }
  else  if (doMatchingShift==2){
  xmlFileName = "config/test_MatchingUp.xml";
  }

 else  if (doJESShift!=0){
   xmlFileName = "config/test_mconly.xml";
  }

 else  if (doLeptonSFShift!=0){
   xmlFileName = "config/test_mconly.xml";
 }

 else  if (doJERShift!=0){
   xmlFileName = "config/test_mconly.xml";
  }

else  if (dobTagEffShift!=0){
  xmlFileName = "config/test_mconly.xml";
  }

else  if (domisTagEffShift!=0){
  xmlFileName = "config/test_mconly.xml";
  }

else  if (doPUShift!=0){
  xmlFileName = "config/test_mconly.xml";
  }

else  if (dottbbShift!=0){
  xmlFileName = "config/test_mconly.xml";
  }
 else if(trainEventMVA==1){
   xmlFileName = "config/train_fullsamples_rereco.xml";
  }
  else{
         xmlFileName = "config/test_fullsamples_rereco.xml";
    // xmlFileName = "config/test_mconly.xml";
  }

 //    xmlFileName = "config/test_mconly.xml";

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
  float Luminosity = 18848.367; //pb^-1??
  vector<string> MVAvars;

  //A few bools to steer the MassReco and Event MVAs
  string MVAmethod = "BDT"; // MVAmethod to be used to get the good jet combi calculation (not for training! this is chosen in the jetcombiner class)
  // bool TrainMVA = false; // If false, the previously trained MVA will be used to calculate stuff
  //bool trainEventMVA = true; // If false, the previously trained MVA will be used to calculate stuff
  //bool computeEventMVA = false;
    cout <<"Instantiating jet combiner..."<<endl;
    
    
   JetCombiner* jetCombiner = new JetCombiner(TrainMVA, Luminosity, datasets, MVAmethod, false);
    cout <<"Instantiated jet combiner..."<<endl;
    double bestTopMass1 =0.;
    double bestTopMass2 = 0.;
    double bestTopMass2ndPass = 0.;
    double bestTopPt =0.;

    //uncomment these two lines for training
    // MVAComputer* Eventcomputer_ =0;
    //MVATrainer* Eventtrainer_ = new MVATrainer("BDT","MasterMVA_Mu_25thFeb", "MasterMVA_Mu_25thFeb.root");

    cout <<"Instantiated jet combiner..."<<endl;
   
    //comment these lines for training 
    MVATrainer* Eventtrainer_ = 0;
    
    
    if(trainEventMVA){
        cout<<"instantiating trainer..."<<endl;
        
    Eventtrainer_->bookWeight("Weight");    
    Eventtrainer_->bookInputVar("HTX");
    Eventtrainer_->bookInputVar("HTH");
    Eventtrainer_->bookInputVar("HTRat");
    Eventtrainer_->bookInputVar("HTb");
    Eventtrainer_->bookInputVar("SumJetMassX");
    Eventtrainer_->bookInputVar("MultiTopness");
    Eventtrainer_->bookInputVar("nTags"); 
    Eventtrainer_->bookInputVar("nJets");
    Eventtrainer_->bookInputVar("Jet5Pt");
    Eventtrainer_->bookInputVar("Jet6Pt");
    }else if (computeEventMVA){
    
    MVAvars.push_back("HTX");
    MVAvars.push_back("HTH");
    MVAvars.push_back("HTRat");
    MVAvars.push_back("HTb");
    MVAvars.push_back("SumJetMassX");
    MVAvars.push_back("MultiTopness");
    MVAvars.push_back("nTags");
    MVAvars.push_back("nJets");
    MVAvars.push_back("Jet5Pt");
    MVAvars.push_back("Jet6Pt");
    cout << " Initialized Eventcomputer_" << endl;

   // MVAComputer* Eventcomputer_ = new MVAComputer("BDT","MasterMVA.root","MasterMVA",MVAvars, "test");

    }
      //comment these lines for training 
  MVAComputer* Eventcomputer_ = new MVAComputer("BDT","MasterMVA_Mu_25thFeb.root","MasterMVA_Mu_25thFeb",MVAvars, "_25thFeb");

    
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
    double frac = 2.;
    double currentLumi;
    double newlumi;
    for (unsigned int d = 0; d < datasets.size (); d++)
    {
    if( datasets[d]->Name()=="TTJets"  ||datasets[d]->Name()=="TTJets_AllHad" || datasets[d]->Name()=="TTJets_Other"){
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


	///heavy flav re-weight
	double ll_rw = 0.976;
	double bb_rw = 3.;
        
	 ttbar_ll->SetEquivalentLuminosity(newlumi/ll_rw);
         ttbar_cc->SetEquivalentLuminosity(newlumi/ll_rw);
	 ttbar_bb->SetEquivalentLuminosity(newlumi/bb_rw);
        

	 // ttbar_ll->SetEquivalentLuminosity(newlumi);                                                
         //ttbar_cc->SetEquivalentLuminosity(newlumi);                                                      
         //ttbar_bb->SetEquivalentLuminosity(newlumi);  


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
    MSPlot["MuonRelIsolation"] = new MultiSamplePlot(datasets, "MuonRelIsolation", 50, 0, .25, "RelIso");
    MSPlot["MuonRelIsolation_PreTrig"] = new MultiSamplePlot(datasets, "MuonRelIsolation_PreTrig", 15, 0, .25, "RelIso");
    MSPlot["MuonRelIsolation_PostTrig"] = new MultiSamplePlot(datasets, "MuonRelIsolation_PreTrig", 15, 0, .25, "RelIso");

    MSPlot["MuonPt"]              = new MultiSamplePlot(datasets, "MuonPt", 35, 0, 300, "PT_{#mu}");
    MSPlot["MuonEta"]              = new MultiSamplePlot(datasets, "MuonEta", 25, -2.4, 2.4, "#eta_{#mu}");
    MSPlot["MuonPhi"]              = new MultiSamplePlot(datasets, "MuonPhi", 50, -4, 4, "#phi_{#mu}");
    MSPlot["MuonNValidHits"]              = new MultiSamplePlot(datasets, "MuonNValidHits", 30, 0, 30, "NValidHits_{#mu}");
    MSPlot["Muond0"]              = new MultiSamplePlot(datasets, "Muond0", 50, -0.05, 0.05, "d0_{#mu}");
    MSPlot["MuondZPVz"]              = new MultiSamplePlot(datasets, "MuondZPVz", 50, 0, .5, "dZPVZ_{#mu}");

    MSPlot["MuondRJets"]              = new MultiSamplePlot(datasets, "MuondRJets", 50, 0, 10, "dRJets_{#mu}");
    MSPlot["MuonNMatchedStations"]              = new MultiSamplePlot(datasets, "MuonNMatchedStations", 10, 0, 10, "NMatchedStations_{#mu}");
    MSPlot["MuonDistVzPVz"]              = new MultiSamplePlot(datasets, "MuonDistVzPVz", 50, 0 ,.3, "DistVzPVz_{#mu}");
    MSPlot["MuonDz"]              = new MultiSamplePlot(datasets, "MuonDz", 25, -.6 ,.6, "Dz_{#mu}");

    MSPlot["MuonTrackerLayersWithMeasurement"]    = new MultiSamplePlot(datasets, "MuonTrackerLayersWithMeasurement", 25, 0, 25, "nLayers");
    MSPlot["DiMuon_InvMass"]     = new MultiSamplePlot(datasets, "DiMuon_InvMass", 60, 0, 120, "DiMuon_InvMass");
    MSPlot["NbOfLooseMuon"]     = new MultiSamplePlot(datasets, "NbOfLooseMuon", 10, 0, 10, "Nb. of loose muons");
    
    //B-tagging discriminators
    MSPlot["BdiscBJetCand_CSV"]  = new MultiSamplePlot(datasets, "BdiscBJetCand_CSV", 75, 0, 1, "CSV b-disc.");

    //Jets
    MSPlot["NbOfSelectedJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedJets", 15, 0, 15, "Nb. of jets");
    MSPlot["NbOfSelectedLightJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedLightJets", 10, 0, 10, "Nb. of jets");
    MSPlot["NbOfSelectedBJets"]                  = new MultiSamplePlot(datasets, "NbOfSelectedBJets", 8, 0, 8, "Nb. of tags");
    MSPlot["JetEta"]                  = new MultiSamplePlot(datasets, "JetEta", 30,-3, 3, "Jet #eta");
    MSPlot["JetPhi"]                  = new MultiSamplePlot(datasets, "JetPhi", 50, -4,4 , "Jet #phi");
    
    MSPlot["NbOfBadTrijets"]                  = new MultiSamplePlot(datasets, "NbOfBadTriJets", 150, 0, 150, "Nb. of Bad Combs");
    
    MSPlot["MVA1TriJetMass"] = new MultiSamplePlot(datasets, "MVA1TriJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA1DiJetMass"] = new MultiSamplePlot(datasets, "MVA1DiJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA1PtRat"] = new MultiSamplePlot(datasets, "MVA1PtRat", 25, 0, 2, "P_{t}^{Rat}");
    MSPlot["MVA1BTag"] = new MultiSamplePlot(datasets, "MVA1BTag", 35, 0, 1, "BTag");
    MSPlot["MVA1AnThBh"] = new MultiSamplePlot(datasets, "MVA1AnThBh", 35, 0, 3.14, "AnThBh");
    MSPlot["MVA1AnThWh"] = new MultiSamplePlot(datasets, "MVA1AnThWh", 35, 0, 3.14, "AnThWh");

    MSPlot["TriJetMass_Matched"] = new MultiSamplePlot(datasets, "TriJetMassMatched", 100, 0, 1000, "m_{bjj}");
    MSPlot["TriJetMass_UnMatched"] = new MultiSamplePlot(datasets, "TriJetMassUnMatched", 100, 0, 1000, "m_{bjj}");
    MSPlot["Chi2_Matched"] = new MultiSamplePlot(datasets, "Chi2Matched", 100, 0, 100, "#chi^{2}");
    MSPlot["Chi2_UnMatched"] = new MultiSamplePlot(datasets, "Chi2UnMatched", 100, 0, 100, "#chi^{2}");
 
    MSPlot["MVA1TriJet"] = new MultiSamplePlot(datasets, "MVA1TriJet", 35, -1., 0.5, "BDT Discriminator");

    //    MSPlot["MVA1TriJet_1A"] = new MultiSamplePlot(datasets, "MVA1TriJet_1A", 35, -1., 0.5, "BDT Discriminator");
    //MSPlot["MVA1TriJet_1B"] = new MultiSamplePlot(datasets, "MVA1TriJet_1B", 35, -1., 0.5, "BDT Discriminator");
    //MSPlot["MVA1TriJet_2A"] = new MultiSamplePlot(datasets, "MVA1TriJet_2A", 35, -1., 0.5, "BDT Discriminator");
    //MSPlot["MVA1TriJet_2B"] = new MultiSamplePlot(datasets, "MVA1TriJet_2B", 35, -1., 0.5, "BDT Discriminator");
    //MSPlot["MVA1TriJet_3"] = new MultiSamplePlot(datasets, "MVA1TriJet_3", 35, -1., 0.5, "BDT Discriminator");

    MSPlot["MVA2TriJetMass"] = new MultiSamplePlot(datasets, "MVA2TriJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA2ndPassTriJetMass"] = new MultiSamplePlot(datasets, "MVA2ndPassTriJetMass", 30, 0, 1000, "m_{bjj}");
    MSPlot["MultiTopness"] = new MultiSamplePlot(datasets, "MultiTopness", 35, -1., 0.5, "MultiTopness");

    //    MSPlot["MVA3rdPassTriJet"] = new MultiSamplePlot(datasets, "MVA3rdPassTriJet", 35, -1., 0, "BDT Discriminator");
    MSPlot["MVA1TriJetMassMatched"] = new MultiSamplePlot(datasets, "MVA1TriJetMassMatched", 75, 0, 500, "m_{bjj}");
    
    //MSPlot["BDisc_Asym"] = new MultiSamplePlot(datasets, "BDisc_Asym", 100, -5, 5, "BDiscAsym");

    MSPlot["SelectedJetPt"] = new MultiSamplePlot(datasets, "JetPt", 50, 0, 300, "PT_{jet}");
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
    MSPlot["MET"] = new MultiSamplePlot(datasets, "MET", 75, 0, 700, "MET");
   
    MSPlot["SumJetMass"] = new MultiSamplePlot(datasets, "SumJetMass", 40, 0, 3000, "SumJetMass");
    MSPlot["SumJetMassX"] = new MultiSamplePlot(datasets, "SumJetMassX", 30, 0, 3000, "SumJetMassX");
    MSPlot["HTX"] = new MultiSamplePlot(datasets, "HTX", 20, 0, 1000, "HTX");
    MSPlot["MVA"] = new MultiSamplePlot(datasets, "MVA", 15, -0.3, 0.4, "BDT Discriminator");
    MSPlot["HTRat"]                  = new MultiSamplePlot(datasets, "HTRat", 50, 0, 20, "HTRat");

    //    //Plots after loose cut on BDT discriminator
      // MSPlot["NbOfSelectedJets_BDTCut"]                  = new MultiSamplePlot(datasets, "NbOfSelectedJets_BDTCut", 15, 0, 15, "Nb. of jets");
    //MSPlot["NbOfSelectedBJets_BDTCut"]                  = new MultiSamplePlot(datasets, "NbOfSelectedBJets_BDTCut", 8, 0, 8, "Nb. of tags");
    //MSPlot["HTX_BDTCut"] = new MultiSamplePlot(datasets, "HTX_BDTCut", 20, 0, 1000, "HTX");
    //MSPlot["HTH_BDTCut"] = new MultiSamplePlot(datasets, "HTH_BDTCut", 50, 0, 1, "HT/H");
    //MSPlot["MultiTopness_BDTCut"] = new MultiSamplePlot(datasets, "MultiTopness_BDTCut", 35, -1., 0.5, "MultiTopness");
    //MSPlot["HTb_SelectedJets_BDTCut"] = new MultiSamplePlot(datasets, "HTb_SelectedJets_BDTCut", 50, 0, 1500, "HTb");

    //Declare arrays of MSPlots
    Int_t minNJets=6, maxNJets=8, minNBJets=2, maxNBJets=3;
    Int_t q =0;

    //  for (Int_t q = minNBJets; q <= maxNBJets; q++){
    for (Int_t p = minNJets; p<= maxNJets; p++){

    string NJets_str = static_cast<ostringstream*>( &(ostringstream() << p) )->str();
    string NBJets_str = static_cast<ostringstream*>( &(ostringstream() << q) )->str();

    string HT_Name = "HT_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string HTX_Name = "HTX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string H_Name = "H_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string HX_Name = "HX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string HTH_Name = "HTH_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string HTXHX_Name = "HTXHX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;

    string EventMass_Name = "EventMass_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string EventMassX_Name = "EventMassX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string SumJetMass_Name = "SumJetMass_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string SumJetMassX_Name = "SumJetMassX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string PTBalTopEventX_Name = "PTBalTopEventX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string PTBalTopSumJetX_Name = "PTBalTopSumJetX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
        
    string FirstTriJetMass_1BJet_g_Name = "FirstTriJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string FirstDiJetMass_1BJet_g_Name = "FirstDiJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string SecondDiJetMass_1BJet_g_Name = "SecondDiJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;

    string SecondTriJetMass_1BJet_g_Name = "SecondTriJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string SecondTriJetMass_1BJet_g_chi2_Name = "SecondTriJetMass_1BJet_g_chi2_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string MuMetBMasses_g_Name = "MuMetBMasses_g_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string MuMetMasses_g_Name = "MuMetMasses_g_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string MET_Name = "MET"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string MVA1TriJetMass_Name = "MVA1TriJetMass"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string MVA2TriJetMass_Name = "MVA2TriJetMass"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
    string MVA_Name = "MVA"+NJets_str+"Jets" ;        

    MSPlot[MVA_Name.c_str()] = new MultiSamplePlot(datasets, MVA_Name.c_str(), 17, -0.5, 0.4, "BDT Discriminator");

    //    MSPlot[MET_Name.c_str() ] = new MultiSamplePlot(datasets, MET_Name.c_str() , 50, 0, 700, "MET");
    //  MSPlot[MuMetMasses_g_Name.c_str() ] = new MultiSamplePlot(datasets, MuMetMasses_g_Name.c_str() , 50, 0, 700, "M_{muMET}");
    // MSPlot[MuMetBMasses_g_Name.c_str() ] = new MultiSamplePlot(datasets, MuMetBMasses_g_Name.c_str() , 50, 0, 700, "M_{muMETb}");
    // MSPlot[HT_Name.c_str() ] = new MultiSamplePlot(datasets, HT_Name.c_str() , 20, 0, 1700, "HT");
    // MSPlot[HTX_Name.c_str() ] = new MultiSamplePlot(datasets, HTX_Name.c_str() , 20, 0, 1200, "HTX");
    // MSPlot[H_Name.c_str() ] = new MultiSamplePlot(datasets, H_Name.c_str() , 20, 0, 1700, "H");
    // MSPlot[HX_Name.c_str() ] = new MultiSamplePlot(datasets, HX_Name.c_str() , 20, 0, 1700, "HX");
    //MSPlot[HTH_Name.c_str() ] = new MultiSamplePlot(datasets, HTH_Name.c_str() , 20, 0, 1, "HT/H");
    //MSPlot[HTXHX_Name.c_str() ] = new MultiSamplePlot(datasets, HTXHX_Name.c_str() , 20, 0, 1, "HTX/HX");
        
    // MSPlot[EventMass_Name.c_str() ] = new MultiSamplePlot(datasets, EventMass_Name.c_str() , 15, 0, 2500, "EventMass");
    //MSPlot[EventMassX_Name.c_str() ] = new MultiSamplePlot(datasets, EventMassX_Name.c_str() , 15, 0, 1700, "EventMassX");
    
    //MSPlot[SumJetMass_Name.c_str() ] = new MultiSamplePlot(datasets, SumJetMass_Name.c_str() , 15, 0, 2500, "SumJetMass");
    //MSPlot[SumJetMassX_Name.c_str() ] = new MultiSamplePlot(datasets, SumJetMassX_Name.c_str() , 15, 0, 1700, "SumJetMassX");
    // MSPlot[PTBalTopEventX_Name.c_str() ] = new MultiSamplePlot(datasets, PTBalTopEventX_Name.c_str() , 15, 0, 500, "PTBalTopEvent_Name");
    //MSPlot[PTBalTopSumJetX_Name.c_str() ] = new MultiSamplePlot(datasets, PTBalTopSumJetX_Name.c_str() , 15, 0, 500, "PTBalTopSumJetX_Name");
        
  
    //    MSPlot[MVA1TriJetMass_Name.c_str() ] = new MultiSamplePlot(datasets, MVA1TriJetMass_Name.c_str() , 40, 20 ,500 , "M_{bjj}");
    //MSPlot[MVA2TriJetMass_Name.c_str() ] = new MultiSamplePlot(datasets, MVA2TriJetMass_Name.c_str() , 40, 20 ,500 , "M_{bjj}");
}
    //}

    for (Int_t p = minNBJets; p<= maxNBJets; p++){

    string NBJets_str = static_cast<ostringstream*>( &(ostringstream() << p) )->str();

    string MVA_Name = "MVA"+NBJets_str+"Tags" ;    
    MSPlot[MVA_Name.c_str()] = new MultiSamplePlot(datasets, MVA_Name.c_str(), 17, -0.5, 0.4, "BDT Discriminator");

 }


  ///////////////////
  // 1D histograms
  ///////////////////
    
     histo1D["btag_weight"] = new TH1F("btag_weight","btag_weight;btag_weight;#events",150,-10,10);
     histo1D["leptonScales"] = new TH1F("leptonScales","leptonScaleFactor;leptonScaleFactor;#events",100,0,4);
     histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);
     histo1D["scaleFactor"] = new TH1F("scaleFactor","scaleFactor;scaleFactor;#events",100,0,4);

     /*


    histo1D["RelIso_PreTrig"] = new TH1F("RelIso_PreTrig","RelIso_PreTrig;RelIso_PreTrig;#events",40,0,0.75);
    histo1D["RelIso_PostTrig"] = new TH1F("RelIso_PostTrig","RelIso_PostTrig;RelIso_PostTrig;#events",40,0,0.75);


 histo1D["4thJetPt"] = new TH1F("4thJetPt","4thJetPt;4thJetPt;#events",25,0,500);

    histo1D["Pt_PreTrig"] = new TH1F("Pt_PreTrig","Pt_PreTrig;Pt_PreTrig;#events",25,0,500);
    histo1D["Pt_PostTrig"] = new TH1F("Pt_PostTrig","Pt_PostTrig;Pt_PostTrig;#events",25,0,500);
    histo1D["d0_PreTrig"] = new TH1F("d0_PreTrig","d0_PreTrig;d0_PreTrig;#events",20,0,0.02);
    histo1D["d0_PostTrig"] = new TH1F("d0_PostTrig","d0_PostTrig;d0_PostTrig;#events",20,0,0.02);
    histo1D["Eta_PreTrig"] = new TH1F("Eta_PreTrig","Eta_PreTrig;Eta_PreTrig;#events",20,-2.5,2.5);
    histo1D["Eta_PostTrig"] = new TH1F("Eta_PostTrig","Eta_PostTrig;Eta_PostTrig;#events",20,-2.5,2.5);
    histo1D["leptonScales"] = new TH1F("leptonScales","leptonScaleFactor;leptonScaleFactor;#events",100,0,4);
    histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);

    histo1D["BDT_Comb_typeAll"] = new TH1F("BDT_Comb_typeAll","BDT_Comb_typeAll;BDT_Comb_typeAll;#combinations",35,-1.,1);
    histo1D["BDT_Comb_type0"] = new TH1F("BDT_Comb_type0","BDT_Comb_type0;BDT_Comb_type0;#combinations",35,-1.,1);
    histo1D["BDT_Comb_type1A"] = new TH1F("BDT_Comb_type1A","BDT_Comb_type1A;BDT_Comb_type1A;#combinations",35,-1.,1);
    histo1D["BDT_Comb_type1B"] = new TH1F("BDT_Comb_type1B","BDT_Comb_type1B;BDT_Comb_type1B;#combinations",35,-1.,1);
    histo1D["BDT_Comb_type2A"] = new TH1F("BDT_Comb_type2A","BDT_Comb_type2A;BDT_Comb_type2A;#combinations",35,-1.,1);
    histo1D["BDT_Comb_type2B"] = new TH1F("BDT_Comb_type2B","BDT_Comb_type2B;BDT_Comb_type2B;#combinations",35,-1.,1);
    histo1D["BDT_Comb_type3"] = new TH1F("BDT_Comb_type3","BDT_Comb_type3;BDT_Comb_type3;#combinations",35,-1.,1);
    
 histo1D["BDT_Comb_0def"] = new TH1F("BDT_Comb_0def","BDT_Comb_0def;BDT_Comb_0def;#combinations",35,-1.,1);
 histo1D["BDT_Comb_1def"] = new TH1F("BDT_Comb_1def","BDT_Comb_1def;BDT_Comb_1def;#combinations",35,-1.,1);
 histo1D["BDT_Comb_2def"] = new TH1F("BDT_Comb_2def","BDT_Comb_2def;BDT_Comb_2def;#combinations",35,-1.,1);
 histo1D["BDT_Comb_3def"] = new TH1F("BDT_Comb_3def","BDT_Comb_3def;BDT_Comb_3def;#combinations",35,-1.,1);

    
  for (unsigned int d = 0; d < datasets.size(); d++){
      
  histoProfile[("HTX_vs_MultiTopness_"+datasets[d]->Name()).c_str()] = new TProfile(("HTX_vs_MultiTopness_"+datasets[d]->Name()).c_str(),"HTX:MultiTopness",30,-1,1, 0,1000);
  histoProfile[("HTX_vs_HTH_"+datasets[d]->Name()).c_str()] = new TProfile(("HTX_vs_HTH_"+datasets[d]->Name()).c_str(),"HTX:HTH",30,0,1, 0,1000);
  histoProfile[("HTX_vs_5thJetPt_"+datasets[d]->Name()).c_str()] = new TProfile(("HTX_vs_5thJetPt_"+datasets[d]->Name()).c_str(),"HTX:5thJetPt",50,0,300, 0,1000);
  histoProfile[("HTX_vs_6thJetPt_"+datasets[d]->Name()).c_str()] = new TProfile(("HTX_vs_6thJetPt_"+datasets[d]->Name()).c_str(),"HTX:6thJetPt",50,0,300, 0,1000);
  histoProfile[("HTX_vs_HTb_"+datasets[d]->Name()).c_str()] = new TProfile(("HTX_vs_HTb_"+datasets[d]->Name()).c_str(),"HTX:HTb",70,0,1400, 0,1000);
  histoProfile[("HTX_vs_HTRat_"+datasets[d]->Name()).c_str()] = new TProfile(("HTX_vs_HTRat_"+datasets[d]->Name()).c_str(),"HTX:HTRat",50,0,20, 0,1000);
  histoProfile[("HTX_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str()] = new TProfile(("HTX_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str(),"HTX:NbOfSelectedJets",15,0,15, 0,1000);
  histoProfile[("HTX_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str()] = new TProfile(("HTX_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str(),"HTX:NbOfSelectedBJets",7,0,7, 0,1000);
  histoProfile[("HTX_vs_SumJetMassX_"+datasets[d]->Name()).c_str()] = new TProfile(("HTX_vs_SumJetMassX_"+datasets[d]->Name()).c_str(),"HTX:SumJetMassX",70,0,1500, 0,1000);
  histoProfile[("MultiTopness_vs_HTH_"+datasets[d]->Name()).c_str()] = new TProfile(("MultiTopness_vs_HTH_"+datasets[d]->Name()).c_str(),"MultiTopness:HTH",30,0,1, -1,1);
  histoProfile[("MultiTopness_vs_5thJetPt_"+datasets[d]->Name()).c_str()] = new TProfile(("MultiTopness_vs_5thJetPt_"+datasets[d]->Name()).c_str(),"MultiTopness:5thJetPt",50,0, 300, -1,1);
  histoProfile[("MultiTopness_vs_6thJetPt_"+datasets[d]->Name()).c_str()] = new TProfile(("MultiTopness_vs_6thJetPt_"+datasets[d]->Name()).c_str(),"MultiTopness:6thJetPt",50,0, 300, -1,1);
  histoProfile[("MultiTopness_vs_HTb_"+datasets[d]->Name()).c_str()] = new TProfile(("MultiTopness_vs_HTb_"+datasets[d]->Name()).c_str(),"MultiTopness:HTb",70,0, 1400, -1,1);
  histoProfile[("MultiTopness_vs_HTRat_"+datasets[d]->Name()).c_str()] = new TProfile(("MultiTopness_vs_HTRat_"+datasets[d]->Name()).c_str(),"MultiTopness:HTRat",50,0, 20, -1,1);
  histoProfile[("MultiTopness_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str()] = new TProfile(("MultiTopness_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str(),"MultiTopness:NbOfSelectedJets",15,0, 15, -1,1);
  histoProfile[("MultiTopness_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str()] = new TProfile(("MultiTopness_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str(),"MultiTopness:NbOfSelectedBJets",7,0, 7, -1,1);
  histoProfile[("MultiTopness_vs_SumJetMassX_"+datasets[d]->Name()).c_str()] = new TProfile(("MultiTopness_vs_SumJetMassX_"+datasets[d]->Name()).c_str(),"MultiTopness:SumJetMassX",70,0, 1500, -1,1);
   
  histoProfile[("HTH_vs_5thJetPt_"+datasets[d]->Name()).c_str()] = new TProfile(("HTH_vs_5thJetPt_"+datasets[d]->Name()).c_str(),"HTH:5thJetPt",50,0,300,0,1);
  histoProfile[("HTH_vs_6thJetPt_"+datasets[d]->Name()).c_str()] = new TProfile(("HTH_vs_6thJetPt_"+datasets[d]->Name()).c_str(),"HTH:6thJetPt",50,0,300,0,1);
  histoProfile[("HTH_vs_HTb_"+datasets[d]->Name()).c_str()] = new TProfile(("HTH_vs_HTb_"+datasets[d]->Name()).c_str(),"HTH:HTb",70,0,1400,0,1);
  histoProfile[("HTH_vs_HTRat_"+datasets[d]->Name()).c_str()] = new TProfile(("HTH_vs_HTRat_"+datasets[d]->Name()).c_str(),"HTH:HTRat",50,0,20,0,1);
  histoProfile[("HTH_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str()] = new TProfile(("HTH_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str(),"HTH:NbOfSelectedJets",15,0,15, 0,1);
  histoProfile[("HTH_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str()] = new TProfile(("HTH_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str(),"HTH:NbOfSelectedBJets",7,0,7, 0,1);
  histoProfile[("HTH_vs_SumJetMassX_"+datasets[d]->Name()).c_str()] = new TProfile(("HTH_vs_SumJetMassX_"+datasets[d]->Name()).c_str(),"HTH:SumJetMassX",70,0,1500, 0,1);
  histoProfile[("5thJetPt_vs_6thJetPt_"+datasets[d]->Name()).c_str()] = new TProfile(("5thJetPt_vs_6thJetPt_"+datasets[d]->Name()).c_str(),"5thJetPt:6thJetPt",50,0,300,0,300);
  histoProfile[("5thJetPt_vs_HTb_"+datasets[d]->Name()).c_str()] = new TProfile(("5thJetPt_vs_HTb_"+datasets[d]->Name()).c_str(),"5thJetPt:HTb",70,0,1400,0,300);
  histoProfile[("5thJetPt_vs_HTRat_"+datasets[d]->Name()).c_str()] = new TProfile(("5thJetPt_vs_HTRat_"+datasets[d]->Name()).c_str(),"5thJetPt:HTRat",50,0,20,0,300);
  histoProfile[("5thJetPt_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str()] = new TProfile(("5thJetPt_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str(),"5thJetPt:NbOfSelectedJets",15,0,15,0,300);
  histoProfile[("5thJetPt_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str()] = new TProfile(("5thJetPt_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str(),"5thJetPt:NbOfSelectedBJets",7,0,7,0,300);
  histoProfile[("5thJetPt_vs_SumJetMassX_"+datasets[d]->Name()).c_str()] = new TProfile(("5thJetPt_vs_SumJetMassX_"+datasets[d]->Name()).c_str(),"5thJetPt:SumJetMassX",70,0,1500,0,300);

  histoProfile[("6thJetPt_vs_HTb_"+datasets[d]->Name()).c_str()] = new TProfile(("6thJetPt_vs_HTb_"+datasets[d]->Name()).c_str(),"6thJetPt:HTb",70,0,1400,0,300);
  histoProfile[("6thJetPt_vs_HTRat_"+datasets[d]->Name()).c_str()] = new TProfile(("6thJetPt_vs_HTRat_"+datasets[d]->Name()).c_str(),"6thJetPt:HTRat",50,0,20,0,300);
  histoProfile[("6thJetPt_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str()] = new TProfile(("6thJetPt_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str(),"6thJetPt:NbOfSelectedJets",15,0,15,0,300);
  histoProfile[("6thJetPt_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str()] = new TProfile(("6thJetPt_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str(),"6thJetPt:NbOfSelectedBJets",7,0,7,0,300);
  histoProfile[("6thJetPt_vs_SumJetMassX_"+datasets[d]->Name()).c_str()] = new TProfile(("6thJetPt_vs_SumJetMassX_"+datasets[d]->Name()).c_str(),"6thJetPt:SumJetMassX",70,0,1500,0,300);

  histoProfile[("HTb_vs_HTRat_"+datasets[d]->Name()).c_str()] = new TProfile(("HTb_vs_HTRat_"+datasets[d]->Name()).c_str(),"HTb:HTRat",50,0,20,0,1500);
  histoProfile[("HTb_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str()] = new TProfile(("HTb_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str(),"HTb:NbOfSelectedJets",15,0,15,0,1500);
  histoProfile[("HTb_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str()] = new TProfile(("HTb_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str(),"HTb:NbOfSelectedJets",7,0,7,0,1500);
  histoProfile[("HTb_vs_SumJetMassX_"+datasets[d]->Name()).c_str()] = new TProfile(("HTb_vs_SumJetMassX_"+datasets[d]->Name()).c_str(),"HTb:SumJetMassX",70,0,1500,0,1500);
  histoProfile[("HTRat_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str()] = new TProfile(("HTRat_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str(),"HTRat:NbOfSelectedJets",15,0,15,0,20);
  histoProfile[("HTRat_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str()] = new TProfile(("HTRat_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str(),"HTRat:NbOfSelectedBJets",7,0,7,0,20);
  histoProfile[("HTRat_vs_SumJetMassX_"+datasets[d]->Name()).c_str()] = new TProfile(("HTRat_vs_SumJetMassX_"+datasets[d]->Name()).c_str(),"HTRat:SumJetMassX",70,0,1500,0,20);
  histoProfile[("NbOfSelectedJets_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str()] = new TProfile(("NbOfSelectedJets_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str(),"NbOfSelectedJets:NbOfSelectedBJets",7,0,7,0,15);
  histoProfile[("NbOfSelectedJets_vs_SumJetMassX_"+datasets[d]->Name()).c_str()] = new TProfile(("NbOfSelectedJets_vs_SumJetMassX_"+datasets[d]->Name()).c_str(),"NbOfSelectedJets:NbOfSelectedBJets",70,0,1500,0,15);
  histoProfile[("NbOfSelectedBJets_vs_SumJetMassX_"+datasets[d]->Name()).c_str()] = new TProfile(("NbOfSelectedBJets_vs_SumJetMassX_"+datasets[d]->Name()).c_str(),"NbOfSelectedBJets:SumJetMassX",70,0,1500,0,7);


histo2D[("RelIso_vs_MET_"+datasets[d]->Name()).c_str()] = new TH2F(("RelIso_vs_MET_"+datasets[d]->Name()).c_str(),"RelIso:MET",100,0,1000, 100, 0,1);
histo2D[("5thJetPt_vs_6thJetPt_"+datasets[d]->Name()).c_str()] = new TH2F(("5thJetPt_vs_6thJetPt_"+datasets[d]->Name()).c_str(),"5thJetPt:6thJetPt",100,0,500, 100, 0,500);

  }
*/    
  
  //Plots
  string pathPNG_MVA = "MVAPlots_"+postfix+channelpostfix;
  string pathPNG = "FourTop"+postfix+channelpostfix;
  pathPNG += "_MSPlots/"; 	
  //pathPNG = pathPNG +"/";
  mkdir(pathPNG.c_str(),0777);
    
    cout <<"Making directory :"<< pathPNG  <<endl;
    
  /////////////////////
  // Selection table: µ + jets
  /////////////////////////////////
  vector<string> CutsSelecTableMu;
  //  CutsSelecTableMu.push_back(string("initial"));
  //  CutsSelecTableMu.push_back(string("PU reweighting"));
    CutsSelecTableMu.push_back(string("Event cleaning and Trigger"));
    CutsSelecTableMu.push_back(string("Exactly 1 isolated muon"));
    CutsSelecTableMu.push_back(string("Loose muon veto"));
    CutsSelecTableMu.push_back(string("Loose electron veto"));
    CutsSelecTableMu.push_back(string("$>$= 6 Jets (30 GeV)"));
    CutsSelecTableMu.push_back(string("$>$= 2 CSVM Tags"));
    CutsSelecTableMu.push_back(string("H_{T} $>$=  400 GeV")); 
    CutsSelecTableMu.push_back(string("E^{Miss}_{T} $>$=  30 GeV"));

    
//  CutsSelecTableMu.push_back(string("Event cleaning and Trigger"));
//  CutsSelecTableMu.push_back(string("Exactly 1 isolated muon"));
//  CutsSelecTableMu.push_back(string("Loose muon veto"));
//  CutsSelecTableMu.push_back(string("Electron veto"));
//  CutsSelecTableMu.push_back(string("$\\geq$ 1 jet"));
//  CutsSelecTableMu.push_back(string("$\\geq$ 2 jet"));
//  CutsSelecTableMu.push_back(string("$\\geq$ 3 jet"));
//  CutsSelecTableMu.push_back(string("$\\geq$ 4 jet"));
//  CutsSelecTableMu.push_back(string("$b\\texttt{-}disc \\geq 0.679$ (CSVM)"));

  FourTopSelectionTable selecTableMu(CutsSelecTableMu, datasets);
  selecTableMu.SetLuminosity(Luminosity);
  selecTableMu.SetPrecision(1);

 FourTopSelectionTable selecTableEl(CutsSelecTableMu, datasets);

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

    }else if(doPUShift==0) {

    LumiWeights = LumiReWeighting("PUReweighting/pileup_MC_S10.root", "PUReweighting/PURewtoptree_id_2217_2227_2262_2423_2435_2417_PileupHistogram.root", "pileup", "pileup");
    }


    //    LumiWeights = LumiReWeighting("PUReweighting/pileup_MC_S10.root", "PUReweighting/PURewtoptree_id_2217_2227_2262_2423_2435_2417_PileupHistogram.root", "pileup", "pileup");
    
    //    if (doPUShift==1)  reweight::PoissonMeanShifter PShiftDown_ = reweight::PoissonMeanShifter(-0.6);
    // else if (doPUShift==2) reweight::PoissonMeanShifter PShiftUp_ = reweight::PoissonMeanShifter(0.6);

   cout << " Initialized LumiReWeighting stuff" << endl;  
    
    // initialize lepton SF
    LeptonTools* leptonTools = new LeptonTools(false);
    leptonTools->readMuonSF("LeptonSF/Muon_ID_iso_Efficiencies_Run_2012ABCD_53X.root", "LeptonSF/MuonEfficiencies_Run_2012A_2012B_53X.root", "LeptonSF/MuonEfficiencies_Run_2012C_53X.root", "LeptonSF/TriggerMuonEfficiencies_Run_2012D_53X.root");
   // leptonTools->readElectronSF();
    
    
/////////////////////////////////
//////// Loop on datasets
/////////////////////////////////

  cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;
    
    TGraph *scaleGraph, *BTagGraph_B,*BTagGraph_C, *BTagGraph_L_eta008, *BTagGraph_L_eta0816, *BTagGraph_L_eta1624;
    
     TFile * btageffFile;
     btageffFile = new TFile("FourTop_MCStudy_Mu.root","READ");
    
    
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
    
     time_t t = time(0);   // get time now
     struct tm * now = localtime( & t );

     int year = now->tm_year + 1900;
     int month =  now->tm_mon + 1;
     int day = now->tm_mday;
     int hour = now->tm_hour;
     int min = now->tm_min; 
     int sec = now->tm_sec;

     string year_str;
     string month_str;
     string day_str;
     string hour_str;
     string min_str;
     string sec_str;

     ostringstream convert;   // stream used for the conversion
     convert << year;      // insert the textual representation of 'Number' in the characters in the stream
     year_str = convert.str();
     convert.str("");
     convert.clear();
     convert << month;      // insert the textual representation of 'Number' in the characters in the stream
     month_str = convert.str();
     convert.str("");
     convert.clear();
     convert << day;      // insert the textual representation of 'Number' in the characters in the stream
     day_str = convert.str();
     convert.str("");
     convert.clear();
     convert << hour;      // insert the textual representation of 'Number' in the characters in the stream
     hour_str = convert.str();
     convert.str("");
     convert.clear();
     convert << min;      // insert the textual representation of 'Number' in the characters in the stream
     min_str = convert.str();
     convert.str("");
     convert.clear();
     convert << day;      // insert the textual representation of 'Number' in the characters in the stream
     sec_str = convert.str();
     convert.str("");
     convert.clear();


     string date_str = day_str + "_" + month_str + "_" + year_str;

     cout <<"DATE STRING   "<<date_str << endl;

     string dataSetName = datasets[d]->Name();	
     string date_dir = "Craneens_Mu/Craneens" + date_str  +"/";
     mkdir(date_dir.c_str(),0777);


     //     string Ntupname = "Craneens/Craneen_" + dataSetName +postfix + "_" + date_str+  ".root";   

     string Ntupname = "Craneens_Mu/Craneens"+ date_str  +"/Craneen_" + dataSetName +postfix + ".root";
     string Ntuptitle = "Craneen_" + dataSetName;

     TFile * tupfile = new TFile(Ntupname.c_str(),"RECREATE");

     TNtuple * tup = new TNtuple(Ntuptitle.c_str(),Ntuptitle.c_str(),"MVA:Njets:Ntags:MultiTopness:HTb:HTH:HTRat:HTX:SumJetMassX:Jet5thPt:Jet6thPt:ScaleFactor:NormFactor:Luminosity");

    
//////////////////////////////////////////////////
/// Initialize JEC factors ///////////////////////
//////////////////////////////////////////////////
		
    vector<JetCorrectorParameters> vCorrParam;

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
      
      
      /*
      if (TrainMVA){
          cout <<"Training Mass Reconstruction MVA"<<endl;
        event_start = 0;
        end = 200000;
      }
      else{
      
       if (trainEventMVA){
          if(  (dataSetName == "NP_overlay_TTTT" ) ){
          event_start = 0;
          end = 9000;
      //  end = datasets[d]->NofEvtsToRunOver();

          }
          else{
              event_start = 200000;
              end = 250000;
            //end = datasets[d]->NofEvtsToRunOver();
          }

      }
      
      else if (computeEventMVA){
      
           if (dataSetName == "Data"){
              event_start = 0;
              
              // end = 100;
              end = datasets[d]->NofEvtsToRunOver();
          }
           else if(  (dataSetName == "NP_overlay_TTTT" ) ){
               event_start = 9000;
               end = datasets[d]->NofEvtsToRunOver();
           }
           else if(  (dataSetName == "TTJets" ) ){
               event_start = 250000;
               end = datasets[d]->NofEvtsToRunOver();
           }
          
      }
      else{
          
          if (dataSetName == "Data"){
              event_start = 0;
              end = datasets[d]->NofEvtsToRunOver();
          }
          
         else if(  (dataSetName == "NP_overlay_TTTT" ) ){
              event_start = 0;
              end = datasets[d]->NofEvtsToRunOver();
             
          }
          else{
              event_start = 200000;
              end = 250000;
          }
      
      }
      }
       
       */
      
      
    if (verbose > 1) cout << " - Loop over events " << endl;
      
      int nBBBar, nCCBar, nLLBar;
      nBBBar=  nCCBar = nLLBar = 0;
      
      double MHT, MHTSig, STJet, EventMass, EventMassX , SumJetMass, SumJetMassX,H,HX , HT, HTX,HTH,HTXHX, sumpx_X, sumpy_X, sumpz_X, sume_X, sumpx, sumpy, sumpz, sume, jetpt,PTBalTopEventX,PTBalTopSumJetX , PTBalTopMuMet;
      
      double currentfrac =0.;
      double end_d = end;
      
      double fakeSig_xs = 100.; //xsection of fake signal in fb.
          
      // to control what fraction of the ttjets we run over(remember to alter the int. lumi of the sample)
      if(dataSetName=="TTJets" || dataSetName=="TTJets_ll" || dataSetName=="TTJets_cc" || dataSetName=="TTJets_bb") {

	if(trainEventMVA){

	event_start = 200000;
        end_d = 250000;
	}else if (computeEventMVA){

	event_start = 250000;
	 end_d =  end_d = 250000. + end_d/frac;

}
          cout <<"N ttjets to run over =   "<< end_d/frac <<endl;
 
      }
  else if (datasets[d]->Name()=="TTJets_AllHad" || datasets[d]->Name()=="TTJets_Other" ) {
             event_start = 0;
         end_d =  end_d/frac;

    //     event_start = 44723;
    // end_d =  end_d;

          cout <<"N ttjets to run over =   "<< end_d/frac <<endl;
      }
else if(dataSetName=="NP_overlay_TTTT"){
	if(trainEventMVA){

	event_start = 0;
	 end_d = 9000;
	}else if (computeEventMVA){

	event_start = 9000;
	 end_d = end;

}

      }else if(dataSetName=="Data"){
          //end_d =  103947. + 250.;
         // end_d = 104197.;
          event_start = 0;
          end_d = end;
      }
      
      else{
          event_start = 0;
          end_d = end;
      }
      cout <<"Will run over "<<  end_d<< " events..."<<endl;

      cout <<"Starting event = = = = "<< event_start  << endl;



      //define object containers
   
        vector<TRootElectron*> selectedElectrons;    
        vector<TRootJet*>      selectedJets;         
        vector<TRootJet*>      selectedJets2ndPass;
        vector<TRootJet*>   MVASelJets1;
	vector<TRootMuon*>     selectedMuons; 
        vector<TRootMuon*>     selectedLooseMuons; 
        vector<TRootElectron*> selectedLooseElectrons;
        vector<TRootJet*>      selectedBJets;
        vector<TRootJet*>      selectedLightJets; 


        vector<TRootJet*>      selectedJets3rdPass;
      	vector<TRootElectron*> selectedExtraElectrons;
       	vector<TRootMuon*>     selectedMuons_NoIso; 
       	vector<TRootMuon*>     selectedExtraMuons;


        selectedElectrons.reserve(10);

        //selectedMuons_NoIso;
	// selectedJets.reserve(20);       
	selectedMuons.reserve(10);      
	selectedLooseMuons.reserve(10);
        selectedLooseElectrons.reserve(10);






    for (unsigned int ievt = event_start; ievt < end_d; ievt++)
    {
        



MHT = 0.,MHTSig = 0., STJet = 0., EventMass =0., EventMassX =0., SumJetMass = 0., SumJetMassX=0.  ,H = 0., HX =0., HT = 0., HTX = 0.,HTH=0.,HTXHX=0., sumpx_X = 0., sumpy_X= 0., sumpz_X =0., sume_X= 0. , sumpx =0., sumpy=0., sumpz=0., sume=0., jetpt =0., PTBalTopEventX = 0., PTBalTopSumJetX =0.;
        
        double ievt_d = ievt;
        currentfrac = ievt_d/end_d;
        if (debug)cout <<"event loop 1"<<endl;
        
	if(ievt%1000 == 0)
		std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";

        if (debug)cout <<"event loop 2"<<endl;
        
     
  //   selecTableEl.Fill(d,0,1.);


	// scale factor for the event
	float scaleFactor = 1.;

	//	init_muons.clear();
	//	init_electrons.clear();
	//init_jets.clear();
	//	mets.clear();
	//	if (vertex) delete vertex;
	//if (ievt) delete 


	//	cout <<"test object "<<  ievt<< vertex << init_muons << init_electrons << init_jets << mets<< endl;

	//load event
	event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);

        if (debug)cout <<"event loop 2.5"<<endl;

	//cout <<" "<<endl;
	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"First cout after loading event:");


        if (debug)cout <<"event loop 3"<<endl;
        
        float rho = event->kt6PFJets_rho();
        string graphName;
        

////////////////////////////////////////////////////////////////////////////
/// Splitting TTBar sample into tt +ll, tt+ cc and tt + bb /////////////////
////////////////////////////////////////////////////////////////////////////
        
        if (debug)cout <<"event loop 4"<<endl;
        
    //load mcparticles to check jet flavour for ttjets events
    vector<TRootMCParticle*> mcParticles_flav;
    Int_t ttbar_flav = -1;

        double nExB,nExC,nExL;
        nExB = nExC = nExL = 0.;
        
    TRootGenEvent* genEvt_flav = 0;
    genEvt_flav = treeLoader.LoadGenEvent(ievt,false);
    treeLoader.LoadMCEvent(ievt, genEvt_flav, 0, mcParticles_flav,false); 
        
        if(  (dataSetName == "TTJets_ll" || dataSetName == "TTJets_cc" || dataSetName == "TTJets_bb" ) )
        {
        for(unsigned int p=0; p<mcParticles_flav.size(); p++) {
          
            if(mcParticles_flav[p]->status()==3 && abs(mcParticles_flav[p]->type())==5 && abs(mcParticles_flav[p]->motherType())!=6) {
               // ttbar_flav=2;
                nExB++;  
            }
            
            else if (mcParticles_flav[p]->status()==3 && abs(mcParticles_flav[p]->type())==4 && abs(mcParticles_flav[p]->motherType())!=6
                && abs(mcParticles_flav[p]->motherType())!=5 && abs(mcParticles_flav[p]->motherType())!=24
                ){
           // ttbar_flav=1;
                 nExC++; 
            }
            
            else if (mcParticles_flav[p]->status()==3 && abs(mcParticles_flav[p]->type())<4 && abs(mcParticles_flav[p]->motherType())!=6){
                // ttbar_flav=1;
                nExL++; 
            }

        }
            
     //   cout <<"TTBar flav composition : " << nExL  <<"  light, " << nExC <<"  C, " << nExB<< "  B" <<  endl;
            
     //   if (ttbar_flav != 1 && ttbar_flav != 2 ) ttbar_flav = 0;
       
            if (nExB >= 2.){
            ttbar_flav =2; 
            nBBBar++ ; //  bbbar
            }
            else if ( nExC >=2.) {
            ttbar_flav =1; 
            nCCBar++ ; //ccbar
            }
            else{
            ttbar_flav =0; 
                nLLBar++;  //llbar   
            }
            
            if (ttbar_flav ==0 && (dataSetName == "TTJets_cc"  || dataSetName == "TTJets_bb"))  continue;
            if (ttbar_flav ==1 && (dataSetName == "TTJets_ll"  || dataSetName == "TTJets_bb" ))  continue;
            if (ttbar_flav ==2 && (dataSetName == "TTJets_ll"  || dataSetName == "TTJets_cc" ))  continue;
        
        }

        if (debug)cout <<"event loop 5"<<endl;

	/////////////////////////////
	/////  Reweighting tt+bb events
	/////////////////////////////
	//	int hf_rw = 2;
	double rw = 1.;
	// Int_t ttbar_flav = -1;

	//	if (dottbbShift != 0){

  if(  (dataSetName == "TTJets" || dataSetName == "TTJets_Other" || dataSetName == "TTJets_AllHad" ) )
    {

  //load mcparticles to check jet flavour for ttjets events
    vector<TRootMCParticle*> mcParticles_flav;


        double nExB,nExC,nExL;
        nExB = nExC = nExL = 0.;
        
    TRootGenEvent* genEvt_flav = 0;
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
           // ttbar_flav=1;
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

	    if (ttbar_flav ==2){
	    if (dottbbShift == 1 ) {

	      rw = 0.5;

	    }else if (dottbbShift == 2){

               rw = 1.5;

	    }else{
	      rw = 1.;
}

    }


}

        if (debug)cout <<"event loop 6"<<endl;
 if(dataSetName != "Data"){
 	scaleFactor = scaleFactor*rw;
 }

	//	cout <<"ttbar flav "<<  ttbar_flav <<"  rw " << rw <<endl;


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

 if(dataSetName != "Data"){

 scaleFactor = scaleFactor*t_pt_rw;
 }

 // cout <<" SF PT "<< t_pt_rw   <<endl;

  }

        if (debug)cout <<"event loop 4"<<endl;

//////////////////
//Loading Gen jets
//////////////////
        
	vector<TRootGenJet*> genjets;
	if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
	{
	  // loading GenJets as I need them for JER
	  		genjets = treeLoader.LoadGenJet(ievt);
	}


	//	cout <<"NGenJets"<<  genjets.size() <<endl;

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
		if(Muon)
		{
                // semi-muon
                if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
                {
                    if( event->runId() >= 190456 && event->runId() <= 190738 ){
                        itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v11"), currentRun, iFile);
                        filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10";
                    }
                    else if( event->runId() >= 190782 && event->runId() <= 193621){
                        itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v12"), currentRun, iFile);
                        filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10";
                    }
                    else if(event->runId() >= 193834  && event->runId() <= 196531 ){
                        itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
                        filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";
                    }
                    else if( event->runId() >= 198022  && event->runId() <= 199608){
                        itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v14"), currentRun, iFile);
                        filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";
                    }
                    else if( event->runId() >= 199698 && event->runId() <= 209151){
                        itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v15"), currentRun, iFile);
                        filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";
                    }
                    else{
                        cout << "Unknown run for SemiMu HLTpath selection: " << event->runId() << endl;
                        filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";

                    }
                    if( itrigger == 9999 )
                    {
                        cout << "itriggerSemiMu == 9999 for SemiMu HLTpath selection: " << event->runId() << endl;
                        //exit(-1);
                    }
                }
    
            else 
	   		{
                  if(dataSetName == "TTJets" || dataSetName == "TTJets_AllHad" ||  dataSetName == "TTJetsOther"  ||   dataSetName == "TTJets_DiLep"    ||dataSetName == "TTJets_ll" || dataSetName == "TTJets_cc"  || dataSetName == "TTJets_bb" ||  dataSetName == "TTJets_powheg_pythia" ||  dataSetName == "TTJets_powheg_herwig" || dataSetName == "TTJets_mcatnlo" || dataSetName == "TTJets_ScaleUp" ||  dataSetName == "TTJets_ScaleDown" ||  dataSetName == "TTJets_MatchingUp"  ||  dataSetName == "TTJets_MatchingDown" ) itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
                
                else  if(dataSetName == "NP_overlay_TTTT") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"),currentRun, iFile);
                
		else  if(dataSetName == "NP_overlay_HG_M500_gt05") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"),currentRun, iFile);

                else  if(dataSetName == "T1TTTT") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"),currentRun, iFile);

                else if (dataSetName == "DY_4Jets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);

                else if (dataSetName == "Z_3Jets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);

                else if (dataSetName == "Z_4Jets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);

                else if (dataSetName == "Z_4Jets_lowmass") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);

      		else if (dataSetName == "W_1Jets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
                
                else if (dataSetName == "W_2Jets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
                
                else if (dataSetName == "W_3Jets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
                
                else if (dataSetName == "W_4Jets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);

                else if (dataSetName == "WW") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
                else if (dataSetName == "WZ") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
                else if (dataSetName == "ZZ") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
                else if (dataSetName == "ttW") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
                else if (dataSetName == "ttZ") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);


                else if (dataSetName == "ZJets") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);

			    else if (dataSetName == "SingleTop_t_T") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
                else if (dataSetName == "SingleTop_t_TBar") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
			    else if (dataSetName == "SingleTop_s_T") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
                else if (dataSetName == "SingleTop_s_TBar") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
       			else if (dataSetName == "SingleTop_tW_T") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
       			else if (dataSetName == "SingleTop_tW_TBar") itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);

				else if (dataSetName == "MultiJet") {
                                                                   
	                 	  itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile); 
				
							  }
				else {

                                        itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
                               }
				       
  				if(itrigger == 9999)
				{
    			  		cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;
					//	exit(1);
				}
			}

        if (debug)cout <<"event loop 7"<<endl;
		} //end if Muon
		else if(Electron)
		{
	
		} //end if Electron
	} //end previousRun != currentRun

////////////////////////////////////////////////////////////////////////////////////
//  JES Corrections: The nominal corrections are already applied at PAT level     //
//    so these tools should only be used for studies of the effect of systematics //
////////////////////////////////////////////////////////////////////////////////////

	// Apply Jet Corrections on-the-fly
	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction on the fly:");
	//	if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
	//		jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),true); //last boolean: isData (needed for L2L3Residual...)
	//	else
	//		jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),false); //last boolean: isData (needed for L2L3Residual...)
	//	coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JES correction on the fly:");

        
	//coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction on the fly:");
  

  // Apply Jet Corrections on-the-fly
      if( dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
      {
        jetTools->unCorrectMETTypeOne(init_jets, mets[0], true);
        jetTools->correctJets(init_jets, event->kt6PFJets_rho(), true);
        jetTools->correctMETTypeOne(init_jets, mets[0], true);
      }
      else
      {
        if (debug)cout <<"event loop 8"<<endl;
        jetTools->unCorrectMETTypeOne(init_jets, mets[0], false);
        jetTools->correctJets(init_jets, event->kt6PFJets_rho(), false);
        jetTools->correctMETTypeOne(init_jets, mets[0], false);
      }


      // coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JES correction on the fly:");



        ///////////////////////
        // JER smearing
        //////////////////////
        
        if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
        {
            //JER
            doJERShift == 0;
            if(doJERShift == 1)
                jetTools->correctJetJER(init_jets, genjets, mets[0], "minus");
            else if(doJERShift == 2)
                jetTools->correctJetJER(init_jets, genjets, mets[0], "plus");
            else
                jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal");
            
	    //     coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JER correction:");
            
            // JES sysematic!
            if (doJESShift == 1)
                jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
            else if (doJESShift == 2)
                jetTools->correctJetJESUnc(init_jets, mets[0], "plus");
            
	    //            coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction:");
            
        }
        
        
    
////////////////////////////////////////
//  Beam scraping and PU reweighting
////////////////////////////////////////
        
	// scale factor for the event
	//float scaleFactor = 1.;

	if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA")
	{
		// Apply the scraping veto. (Is it still needed?)
        	bool isBeamBG = true;
        	if(event->nTracks() > 10)
        	{
			if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
			isBeamBG = false;
		}
      		if(isBeamBG) continue;
	}
	else{
	double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
	double lumiWeightOLD=lumiWeight;
        if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0){
	  lumiWeight=1;
        }
	scaleFactor = scaleFactor*lumiWeight;

	}

		histo1D["lumiWeights"]->Fill(scaleFactor);	
			
///////////////////////////////////////////////////////////
//   Event selection
///////////////////////////////////////////////////////////
	       
	// Apply trigger selection
	trigged = treeLoader.EventTrigged (itrigger);
	if (debug)cout<<"triggered? Y/N?  "<< trigged  <<endl;


	//Applying trigger selection again with 2012 Muon+Jets triggers.
	if(!trigged)		   continue;

	// Declare selection instance    
	Selection selection(init_jets, init_muons, init_electrons, mets);

	// Define object selection cuts
	selection.setJetCuts(30.,2.5,0.01,1.,0.98,0,0);//Pt, Eta, EMF, n90Hits, fHPD, dRJetElectron, DRJets

        selection.setElectronCuts();//	selection.setElectronCuts(30.,2.5,0.1,0.02,0.,999.,0,1); //Pt,Eta,RelIso,d0,MVAId,DistVzPVz, DRJets, MaxMissingHits
	selection.setLooseElectronCuts(20,2.5,0.15,0.);

	//selection.setMuonCuts(26.,2.1,.12,0,0.2,0,1,0.5,5,1 ); //Pt,Eta,RelIso,NValidMuHits,d0, dRJets, NMatchedStations,DistVzPVz,NTrackerLayersWithMeas
	//	selection.setMuonCuts(26.0,2.1,0.12,0.2,999.,1,0.5,5,0 );  // Pt,  Eta, RelIso, d0,  DRJets, NMatchedStations, Dz,  NTrackerLayersWithMeas, NValidPixelHits
  
	selection.setMuonCuts();

      selection.setLooseMuonCuts(10,2.5,0.2);
	  
	//Select objects 
      /*
	//	vector<TRootElectron*> selectedElectrons_NoIso = selection.GetSelectedElectrons(20,2.4,999.,vertex[0]);
	//vector<TRootElectron*> selectedElectrons        = selection.GetSelectedElectrons(vertex[0]);
   
        vector<TRootElectron*> selectedElectrons        = selection.GetSelectedElectrons();
	//	vector<TRootElectron*> selectedExtraElectrons;
	//	vector<TRootMuon*>     selectedMuons_NoIso      = selection.GetSelectedMuons(26,2.4,999.); 
	//	vector<TRootMuon*>     selectedExtraMuons;
	vector<TRootJet*>      selectedJets             = selection.GetSelectedJets(true); // ApplyJetId
        vector<TRootJet*>      selectedJets2ndPass;
        vector<TRootJet*>   MVASelJets1;
	// vector<TRootJet*>      selectedJets3rdPass;
	// vector<TRootElectron*> selectedElectrons_NoIso  = selection.GetSelectedElectrons(30,2.4,999.);
	vector<TRootMuon*>     selectedMuons = selection.GetSelectedMuons();
        //vector<TRootJet*>      selectedSoftJets         = selection.GetSelectedJets(20.,2.5, selectedMuons, 0., true); // ApplyJetId
        vector<TRootMuon*>     selectedLooseMuons       = selection.GetSelectedLooseMuons(10., 2.5, 0.2);
        vector<TRootElectron*> selectedLooseElectrons   = selection.GetSelectedLooseElectrons(20.,2.5,0.15); // VBTF ID
        vector<TRootJet*>      selectedBJets; // B-Jets
        vector<TRootJet*>      selectedLightJets; // light-Jets
	// vector<TRootJet*>       selectedCSVOrderedJets     = selection.GetSelectedJets(true); //CSV ordered Jet collection added by JH

      
        selectedElectrons.reserve(10);
        //selectedMuons_NoIso;
        selectedJets.reserve(20);       
	selectedMuons.reserve(10);      
	selectedLooseMuons.reserve(10);
        selectedLooseElectrons.reserve(10);
      */
        selectedElectrons        = selection.GetSelectedElectrons();
	//	selectedMuons_NoIso      = selection.GetSelectedMuons(26,2.4,999.); 
	selectedJets             = selection.GetSelectedJets(true); // ApplyJetId
        selectedMuons            = selection.GetSelectedMuons();
        selectedLooseMuons       = selection.GetSelectedLooseMuons(10., 2.5, 0.2);
        selectedLooseElectrons   = selection.GetSelectedLooseElectrons(20.,2.5,0.15); // VBTF ID

	//	cout <<"SEL JETS SIZE = A   "<<selectedJets.size() << endl;
         

	//order jets wrt to Pt, then set bool corresponding to RefSel cuts.
	// sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order muons wrt Pt.                                                                    
	//sort(selectedCSVOrderedJets.begin(), selectedCSVOrderedJets.end(), HighestCVSBtag()); //order Jets wrt CSVtag

    int JetCut =0;
    int nMu = selectedMuons.size();
    int nEl = selectedElectrons.size();


    //////////////////////////////////////////////////
    ///// Calculate all Scalefactors first..
    //////////////////////////////////////////////////
        
        bool isTagged =false;
        int seljet;
        
        //boolean which controls the application of the Btag scalefactor
    int applyBTagSF = 1; // 1 = EventWeigthing using SF and MC efficiencies.
 
    //    btageffFile->GetObject("TGraph/BTagEfficiency_e_b_TTJets", BTagGraph_B);
    //btageffFile->GetObject("TGraph/BTagEfficiency_e_c_TTJets", BTagGraph_C);
    //btageffFile->GetObject("TGraph/BTagEfficiency_e_l_TTJets_eta0-08", BTagGraph_L_eta008);
    //btageffFile->GetObject("TGraph/BTagEfficiency_e_l_TTJets_eta08-16", BTagGraph_L_eta0816);
    //btageffFile->GetObject("TGraph/BTagEfficiency_e_l_TTJets_eta16-24", BTagGraph_L_eta1624);
       

    TGraph * BTagGraph_B = (TGraph*)btageffFile->Get("TGraph/BTagEfficiency_e_b_TTJets");
    TGraph * BTagGraph_C = (TGraph*)btageffFile->Get("TGraph/BTagEfficiency_e_c_TTJets");
    TGraph * BTagGraph_L_eta008 = (TGraph*)btageffFile->Get("TGraph/BTagEfficiency_e_l_TTJets_eta0-08");
    TGraph * BTagGraph_L_eta0816 = (TGraph*)btageffFile->Get("TGraph/BTagEfficiency_e_l_TTJets_eta08-16");
    TGraph * BTagGraph_L_eta1624 = (TGraph*)btageffFile->Get("TGraph/BTagEfficiency_e_l_TTJets_eta16-24");

    //	cout <<"SEL JETS SIZE = B   "<<selectedJets.size() << endl;

  if (applyBTagSF ==1){

    //cout <<"SEL JETS SIZE = C   "<<selectedJets.size() << endl;

    
         for ( seljet =0; seljet < selectedJets.size(); seljet++ ){
	   //  cout <<"seljet  "<< seljet  << endl;
                if (selectedJets[seljet]->btag_combinedSecondaryVertexBJetTags() > workingpointvalue   ){

		  //cout <<"tagged  "  << endl;
                    selectedBJets.push_back(selectedJets[seljet]);
                }
                
                else{
		  //	   cout <<"untagged  " << endl;
                    selectedLightJets.push_back(selectedJets[seljet]);

                }
		//         cout <<"running tag total  = "<< selectedBJets.size()  <<endl;

            }

	 //	 cout <<"Ntags = "<< selectedBJets.size()  <<endl;



	 //cout <<"SEL JETS SIZE = D   "<<selectedJets.size() << endl;

	 //	 cout <<"  "<<endl;

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
  //  cout <<" "<<endl;		 
  ////cout <<"jet SF nom "<< bTool->getWeight(selectedJets[seljet]->Pt(),selectedJets[seljet]->Eta(),jet_flavor,0 )    <<endl;
  //cout <<"jet SF minus "<< bTool->getWeight(selectedJets[seljet]->Pt(),selectedJets[seljet]->Eta(),jet_flavor,-1 )    <<endl;
  //cout <<"jet SF plus "<< bTool->getWeight(selectedJets[seljet]->Pt(),selectedJets[seljet]->Eta(),jet_flavor,1 )    <<endl;
}
else{
  //  cout <<" light jet... "<<endl;
  SF_tag =  bTool->getWeight(selectedJets[seljet]->Pt(),selectedJets[seljet]->Eta(),jet_flavor,domisTagEffShift);
}

      if (fabs(jet_flavor) == 5  ) eff =  BTagGraph_B->Eval(JetPt,0,"");
      else  if( fabs(jet_flavor) == 4  ) eff =  BTagGraph_C->Eval(JetPt,0,"");
      else if(fabs(JetEta) < 0.8) eff =  BTagGraph_L_eta008->Eval(JetPt,0,"");
      else if(fabs(JetEta) >= 0.8  &&  fabs(JetEta) < 1.6 ) eff =  BTagGraph_L_eta0816->Eval(JetPt,0,"");
      else if(fabs(JetEta) >= 1.6 &&  fabs(JetEta) <= 2.4)  eff =  BTagGraph_L_eta1624->Eval(JetPt,0,"");

      a_eff = 1 - eff;
      scaled_eff = SF_tag*eff;
      sf_a_eff = 1 - scaled_eff;

      //      cout <<"eff  "<< eff <<endl;
      // cout <<"scaled_eff  "<< scaled_eff <<endl;
      // cout <<"sf_a_eff  "<< sf_a_eff <<endl;
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

   event_weight = (p_data)/(p_mc);

// cout <<"event_weight  "<< event_weight  <<endl;




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
      cout <<" "<<endl;
      cout <<" "<<endl;
      cout <<" "<<endl;
 }
 
 // cout<<"applying event weight... "<<event_weight <<endl; 

 if (dataSetName != "Data" && event_weight > 0. && event_weight < 10.) scaleFactor = scaleFactor*event_weight; 
 histo1D["btag_weight"]->Fill(event_weight);
}
        else{ //not applying SF
            if (debug)cout <<"NOT APPLYING SF"<<endl;
            
            
            for ( seljet =0; seljet < selectedJets.size(); seljet++ ){
                if (selectedJets[seljet]->btag_combinedSecondaryVertexBJetTags() > workingpointvalue   ){
                    
		  cout<<"pushing back 2"<<endl;
                    selectedBJets.push_back(selectedJets[seljet]);
                }
                
                else{
                    selectedLightJets.push_back(selectedJets[seljet]);
                }

            }
        
        }
        
 if (dataSetName != "Data" && selectedMuons.size() ==1 ) {
	   scaleFactor = scaleFactor*leptonTools->getMuonSF(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), leptonsyst);
	   	   histo1D["leptonScales"]->Fill(leptonTools->getMuonSF(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), leptonsyst));
	   // cout << "lepton SF " << leptonTools->getMuonSF(selectedMuons[0]->Eta(), selectedMuons[0]->Pt(), leptonsyst)  <<endl; 

        }


    int nbadcombs = Factorial(selectedJets.size())    /  (  (Factorial(3)) * (Factorial(  selectedJets.size() - 3     ) )       ) ;
         
            
        int njets = selectedJets.size();
        int ntags = selectedBJets.size();
	selectedBJets.clear();

	//       cout <<"***********"<<endl;
	// cout<<"ntags = = " << ntags << endl;
	//cout <<"***********"<<endl;
	//cout <<"   "<<endl;

        double temp_HT = 0.;
        double HTRat = 0.;
        double HT_leading = 0.;
        double HT_lagging = 0.;
        for (Int_t seljet0 =0; seljet0 < selectedJets.size(); seljet0++ ){
            temp_HT += selectedJets[seljet0]->Pt();
	    if (seljet0 < 4){
HT_leading += selectedJets[seljet0]->Pt();
	    }else{

HT_lagging += selectedJets[seljet0]->Pt();
}
        }

	HTRat = HT_leading/HT_lagging;

         double HTb = 0.;
        for (Int_t seljet1 =0; seljet1 < ntags; seljet1++ ){
            HTb += selectedBJets[seljet1]->Pt();
        }

        ////////////////////
       // Sync'ing cutflow
       /////////////////////

   if (debug)	cout <<" applying baseline event selection..."<<endl;
      	// Apply primary vertex selection
	bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
        if(isGoodPV && trigged){

 selecTableMu.Fill(d,0,scaleFactor);

 if (nMu==1){
        selecTableMu.Fill(d,1,scaleFactor);
        if(selectedLooseMuons.size() ==1 ){
               selecTableMu.Fill(d,2,scaleFactor);
                     if(selectedLooseElectrons.size() ==0 ){
                          selecTableMu.Fill(d,3,scaleFactor);
                          if (selectedJets.size() >= 6){
                             selecTableMu.Fill(d,4,scaleFactor);
                                if (ntags >= 2){
                                     selecTableMu.Fill(d,5,scaleFactor);
                                     if (temp_HT >= 400.){
                                          selecTableMu.Fill(d,6,scaleFactor);
                                          if (mets[0]->Et() >= 30.){
                                          selecTableMu.Fill(d,7,scaleFactor);
                                                        }

                                                        }
                                                        }
                                                        }
                                                        }
                                                        }
                                                        }
                                                        }

	///////////////////////////
	/////Applying baseline selection
	/////
	///////////////////////////
	if (!trigged) continue;
	if (!isGoodPV) continue;

       if (!(selectedJets.size() >= 6)) continue;

	if(Muon){

    if (debug) cout <<"Number of Muons, Jets, BJets, JetCut  ===>  "<< selectedMuons.size() <<"  "  << selectedJets.size()   <<"  " <<  ntags   <<"  "<<JetCut  <<endl;
        int nTightLeptons, nVetoLeptonsSF, nVetoLeptonsOF;
        
   if (debug)	cout <<" applying baseline event selection..."<<endl;
      	// Apply primary vertex selection
//Apply the lepton, btag and HT selections
   if  (  !( nMu == 1 && selectedLooseMuons.size() ==1 && ntags >= 2    &&selectedLooseElectrons.size()   ==0 && temp_HT > 400. && mets[0]->Et() > 30.    )) continue;
 passed++;

        
 if (debug) cout <<"Event passed..22."<<endl;

          
        vector<TLorentzVector*> selectedMuonTLV_JC;
        selectedMuonTLV_JC.push_back(selectedMuons[0]);        

         
         vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV, mcMuonsTLV;
        vector<TRootMCParticle*> mcParticlesMatching_;
        bool muPlusFromTop = false, muMinusFromTop = false;
        bool elPlusFromTop = false, elMinusFromTop = false;

	//comment these when not using gen event
	//        pair<unsigned int, unsigned int> leptonicBJet_, hadronicBJet_, hadronicWJet1_, hadronicWJet2_; //First index is the JET number, second one is the parton
        //leptonicBJet_ = hadronicBJet_ = hadronicWJet1_ = hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
        
        
             if (debug) cout <<"trigger objects..1"<<endl;

        
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
        mcParticles.clear();                     
        mcTops.clear();

        int leptonPDG, muonPDG = 13, electronPDG = 11;
        leptonPDG = muonPDG;

	//////
	///// Get Gen Event (can be switched off)
	//////
	/*
        TRootGenEvent* genEvt = 0;
        genEvt = treeLoader.LoadGenEvent(ievt,false);
        sort(selectedJets.begin(),selectedJets.end(),HighestPt()); 
        treeLoader.LoadMCEvent(ievt, genEvt, 0, mcParticles,false);  
    
        if (debug) cout <<"size   "<< mcParticles.size()<<endl;

                  if (debug) cout <<"getting gen event."<<endl;

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
         
	*/
	//End of get gen event
         }
        

         if (debug) cout <<"filling histos..."<<endl;
	 //	 cout <<"Scale factor ==> " <<    scaleFactor <<"   luminosity  "<< Luminosity   <<endl;

 
//////////////////////////////////////
////Filling histograms / plotting
//////////////////////////////////////

//	 cout <<" here...."<<  Luminosity <<endl;
//	 cout <<" here...."<<  scaleFactor <<endl;
//	 cout <<" "<< endl;

        histo1D["scaleFactor"]->Fill(scaleFactor);

        MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);
        


/////////////////
////Muons
/////////////////
        
	    MSPlot["NbOfIsolatedMuons"]->Fill(selectedMuons.size(), datasets[d], true, Luminosity*scaleFactor);

        if (debug) cout <<"filling muid   "<< selectedMuons[0]->Pt()<<endl;
 
       float reliso = (selectedMuons[0]->chargedHadronIso() + max( 0.0, selectedMuons[0]->neutralHadronIso() + selectedMuons[0]->photonIso() - 0.5*selectedMuons[0]->puChargedHadronIso() ) ) / selectedMuons[0]->Pt();

        double muzPVz = fabs(selectedMuons[0]->vz() - vertex[0]->Z());
        MSPlot["MuonDz"]->Fill( selectedMuons[0]->dz() , datasets[d], true, Luminosity*scaleFactor );
        MSPlot["MuonPt"]->Fill(selectedMuons[0]->Pt(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["MuonEta"]->Fill(selectedMuons[0]->Eta(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["MuonPhi"]->Fill(selectedMuons[0]->Phi(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["MuonNValidHits"]->Fill(selectedMuons[0]->nofValidHits(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["Muond0"]->Fill(selectedMuons[0]->d0(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["MuonDistVzPVz"]->Fill(muzPVz, datasets[d], true, Luminosity*scaleFactor );
        MSPlot["MuonNMatchedStations"]->Fill(selectedMuons[0]->nofMatchedStations(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["MuonTrackerLayersWithMeasurement"]->Fill(selectedMuons[0]->nofTrackerLayersWithMeasurement(), datasets[d], true, Luminosity*scaleFactor);
        MSPlot["MuonRelIsolation"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
        
      if (debug) cout <<"filled all muID plots  .."<<endl;
        
      /*
       
        if(dataSetName != "data" && dataSetName != "DATA" && dataSetName != "Data") {
            for(unsigned int i=0; i<mcParticles.size(); i++) {
                if( abs(mcParticles[i]->type()) == 13 && mcParticles[i]->status() == 3) { //Matrix Element Level Muon
                    mcMuonsTLV.push_back(*mcParticles[i]);
                }
            }
            vector<TLorentzVector> selectedMuonTLV;
            selectedMuonTLV.push_back(*selectedMuons[0]);
            JetPartonMatching muonMatching = JetPartonMatching(mcMuonsTLV, selectedMuonTLV, 2, true, true, 0.3);
            
            for(unsigned int i=0; i<mcMuonsTLV.size(); i++) {
                int matchedMuonNumber = muonMatching.getMatchForParton(i, 0);
                if(matchedMuonNumber != -1) MSPlot["LeptonTruth"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
                else MSPlot["LeptonTruth"]->Fill(1, datasets[d], true, Luminosity*scaleFactor);
            }
        }


      */        

/////////////
////Jets
/////////////
	  MSPlot["NbOfSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
	  MSPlot["NbOfSelectedLightJets"]->Fill(selectedLightJets.size(), datasets[d], true, Luminosity*scaleFactor);
          MSPlot["NbOfSelectedBJets"]->Fill(ntags, datasets[d], true, Luminosity*scaleFactor);
	  //  MSPlot["RhoCorrection"]->Fill(event->kt6PFJetsPF2PAT_rho(), datasets[d], true, Luminosity*scaleFactor);
	  if (debug) cout <<"per jet plots.."<<endl;

	//plots to to inspire staggered Jet Pt selection
	  //	  histo1D["4thJetPt"]->Fill(selectedJets[3]->Pt());
      if (selectedJets.size()>=4) MSPlot["4thJetPt"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      if (selectedJets.size()>=5) MSPlot["5thJetPt"]->Fill(selectedJets[4]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      if (selectedJets.size()>=6) MSPlot["6thJetPt"]->Fill(selectedJets[5]->Pt(), datasets[d], true, Luminosity*scaleFactor);
      if (selectedJets.size()>=7) MSPlot["7thJetPt"]->Fill(selectedJets[6]->Pt(), datasets[d], true, Luminosity*scaleFactor);
	
    if (debug) cout <<"got muons and mets"<<endl;
        
    //    histo2D[("5thJetPt_vs_6thJetPt_"+datasets[d]->Name()).c_str()]->Fill(selectedJets[5]->Pt(), selectedJets[4]->Pt());


/////////////////////////////////
/// Find indices of jets from Tops (optional)
////////////////////////////////
/*
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
*/
 
  //  cout <<"  "<<endl;
    //   if (debug) cout <<"Indices of matched jets are :  "<< hadronicBJet_.first<<"  "<< hadronicWJet1_.first  <<" " << hadronicWJet2_.first <<endl;
        
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


	    selectedMuonTLV_JC.clear();


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
      //   for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ ){
      //  if (seljet1 == MVAvals1.second[0] || seljet1 == MVAvals1.second[1] || seljet1 == MVAvals1.second[2] || seljet1 == MVAvals2ndPass.second[0] || seljet1 == MVAvals2ndPass.second[1] || seljet1 == MVAvals2ndPass.second[2] ) continue;
      //  selectedJets3rdPass.push_back(selectedJets[seljet1]);
      // }
    

    //    if(selectedJets3rdPass.size()>= 3){
    // jetCombiner->ProcessEvent_SingleHadTop(datasets[d], mcParticles_flav, selectedJets3rdPass, selectedMuonTLV_JC[0], genEvt_flav, scaleFactor);
    // MVAvals3rdPass = jetCombiner->getMVAValue(MVAmethod, 1);
    //
    // MSPlot["MVA3rdPassTriJet"]->Fill(MVAvals3rdPass.first,  datasets[d], true, Luminosity*scaleFactor );
    // }

    MultiTopness = MVAvals2ndPass.first;
    
    MSPlot["MultiTopness"]->Fill(MultiTopness,  datasets[d], true, Luminosity*scaleFactor );


          // cout <<"sel jets 1st pass = "<< selectedJets.size() << "sel jets 2nd pass =  " << selectedJets2ndPass.size() << endl;
            
    //         bestTopMass1 =( *selectedJets[MVAvals1.second[0]] + *selectedJets[MVAvals1.second[1]] + *selectedJets[MVAvals1.second[2]]).M();
    //     bestTopMass2ndPass =( *selectedJets[MVAvals2ndPass.second[0]] + *selectedJets[MVAvals2ndPass.second[1]] + *selectedJets[MVAvals2ndPass.second[2]]).M();

    //   bestTopPt =( *selectedJets[MVAvals1.second[0]] + *selectedJets[MVAvals1.second[1]] + *selectedJets[MVAvals1.second[2]]).Pt();

       // cout <<"Indices of best MVA jets are :  "<< MVAvals1.second[0] <<"  "<< MVAvals1.second[1]  <<" " << MVAvals1.second[2]<<endl;
            
         //   cout <<"MVA Mass 1 = "<< bestTopMass1 << " MVA Mass 2 = "<< bestTopMass2ndPass << endl;
    
   // cout <<"   "<<endl;
    
    //   MSPlot["MVA1TriJetMass"]->Fill(bestTopMass1,  datasets[d], true, Luminosity*scaleFactor );
    //  MSPlot["MVA2ndPassTriJetMass"]->Fill(bestTopMass2ndPass,  datasets[d], true, Luminosity*scaleFactor );

       if (debug)  cout <<"MVA Mass 1 = "<< bestTopMass1 << " MVA Mass 2 = "<< bestTopMass2 << endl;

        ////////////////////////////////////////////////////////////////////////////////////////////////////
        /////  Calculating how well the MVA jet selection is doing: Fraction of ttbar events            ////
        ////    where the jets selected by the TMVA massReco match the true jets from the hadronic top) ////
        ////////////////////////////////////////////////////////////////////////////////////////////////////
        
//dividing bad combinations into sub-categories to understand structure in BDT discriminant
		      // Class 1: exactly 1 incorrect jet
		      //          1A: the b-jet is not selected
		      //          1B: one of the jets from the W is not selected.
		      // Class 2: exactly 2 incorrect jets
		      //          2A: both jets from the W are incorrect
		      //          2B: 1 jet from the b and one jet from the W are incorrect
		      // Class3:  All jets are incorrect.
	
	   // now make appropriate combinations of random bad jets with good jets to model
	   //each sub-category
       //histo1D["BDT_Comb_typeAll"]->Fill(MVAvals1.first);



//histo1D["BDT_Comb_1def"]

	//START OF OPTIONAL ANLAYIS OF COMBINATORICS

/*
  ndefs =0;

 for (int def = 0;  def < 3; def++){
   if(selectedJets[MVAvals1.second[def]]->btag_combinedSecondaryVertexBJetTags()  <  0) ndefs++;
}


  if(   ( hadronicBJet_.first != MVAvals1.second[0] || hadronicBJet_.first != MVAvals1.second[1] || hadronicBJet_.first != MVAvals1.second[2]   )  || ( hadronicWJet1_.first != MVAvals1.second[0] || hadronicWJet1_.first != MVAvals1.second[1] || hadronicWJet1_.first != MVAvals1.second[2]   )  || ( hadronicWJet2_.first != MVAvals1.second[0] || hadronicWJet2_.first != MVAvals1.second[1] || hadronicWJet2_.first != MVAvals1.second[2]   )      ){
           
    if(ndefs ==0) histo1D["BDT_Comb_0def"]->Fill(MVAvals1.first);
   else if(ndefs ==1) histo1D["BDT_Comb_1def"]->Fill(MVAvals1.first);
   else if(ndefs ==2) histo1D["BDT_Comb_2def"]->Fill(MVAvals1.first);
   else if(ndefs ==3) histo1D["BDT_Comb_3def"]->Fill(MVAvals1.first);

    }





        if(   ( hadronicBJet_.first == MVAvals1.second[0] || hadronicBJet_.first == MVAvals1.second[1] || hadronicBJet_.first == MVAvals1.second[2]   )  && ( hadronicWJet1_.first == MVAvals1.second[0] || hadronicWJet1_.first == MVAvals1.second[1] || hadronicWJet1_.first == MVAvals1.second[2]   )    && ( hadronicWJet2_.first == MVAvals1.second[0] || hadronicWJet2_.first == MVAvals1.second[1] || hadronicWJet2_.first == MVAvals1.second[2]   )      ){
            nMVASuccesses++;
histo1D["BDT_Comb_type0"]->Fill(MVAvals1.first);
        }

else if (   ( hadronicBJet_.first != MVAvals1.second[0] || hadronicBJet_.first != MVAvals1.second[1] || hadronicBJet_.first != MVAvals1.second[2]   ) && ( hadronicWJet1_.first == MVAvals1.second[0] || hadronicWJet1_.first == MVAvals1.second[1] || hadronicWJet1_.first == MVAvals1.second[2]   )    && ( hadronicWJet2_.first == MVAvals1.second[0] || hadronicWJet2_.first == MVAvals1.second[1] || hadronicWJet2_.first == MVAvals1.second[2]   )   ){

	  //This is a type-1A combination 
MSPlot["MVA1TriJet_1A"]->Fill(MVAvals1.first   ,  datasets[d], true, Luminosity*scaleFactor );
histo1D["BDT_Comb_type1A"]->Fill(MVAvals1.first);


	}else if (  ( hadronicBJet_.first == MVAvals1.second[0] || hadronicBJet_.first == MVAvals1.second[1] || hadronicBJet_.first == MVAvals1.second[2]   )  && ( hadronicWJet1_.first != MVAvals1.second[0] || hadronicWJet1_.first != MVAvals1.second[1] || hadronicWJet1_.first != MVAvals1.second[2]   )    && ( hadronicWJet2_.first == MVAvals1.second[0] || hadronicWJet2_.first == MVAvals1.second[1] || hadronicWJet2_.first == MVAvals1.second[2]   )          ){

 //This is a type-1B combination 

MSPlot["MVA1TriJet_1B"]->Fill(MVAvals1.first   ,  datasets[d], true, Luminosity*scaleFactor );
histo1D["BDT_Comb_type1B"]->Fill(MVAvals1.first);

}
	else if ( ( hadronicBJet_.first == MVAvals1.second[0] || hadronicBJet_.first == MVAvals1.second[1] || hadronicBJet_.first == MVAvals1.second[2]   )  && ( hadronicWJet1_.first != MVAvals1.second[0] || hadronicWJet1_.first != MVAvals1.second[1] || hadronicWJet1_.first != MVAvals1.second[2]   )    && ( hadronicWJet2_.first != MVAvals1.second[0] || hadronicWJet2_.first != MVAvals1.second[1] || hadronicWJet2_.first != MVAvals1.second[2]   )     ){

//This is a type-2A combination 
MSPlot["MVA1TriJet_2A"]->Fill(MVAvals1.first   ,  datasets[d], true, Luminosity*scaleFactor );
histo1D["BDT_Comb_type2A"]->Fill(MVAvals1.first);
	}

else if (   ( hadronicBJet_.first != MVAvals1.second[0] || hadronicBJet_.first != MVAvals1.second[1] || hadronicBJet_.first != MVAvals1.second[2]   )  && ( hadronicWJet1_.first != MVAvals1.second[0] || hadronicWJet1_.first != MVAvals1.second[1] || hadronicWJet1_.first != MVAvals1.second[2]   )    && ( hadronicWJet2_.first == MVAvals1.second[0] || hadronicWJet2_.first == MVAvals1.second[1] || hadronicWJet2_.first == MVAvals1.second[2]   )         ){

//This is a type-2B combination 
MSPlot["MVA1TriJet_2B"]->Fill(MVAvals1.first   ,  datasets[d], true, Luminosity*scaleFactor );
histo1D["BDT_Comb_type2B"]->Fill(MVAvals1.first);
	}else if ( ( hadronicBJet_.first != MVAvals1.second[0] || hadronicBJet_.first != MVAvals1.second[1] || hadronicBJet_.first != MVAvals1.second[2]   )  && ( hadronicWJet1_.first != MVAvals1.second[0] || hadronicWJet1_.first != MVAvals1.second[1] || hadronicWJet1_.first != MVAvals1.second[2]   )    && ( hadronicWJet2_.first != MVAvals1.second[0] || hadronicWJet2_.first != MVAvals1.second[1] || hadronicWJet2_.first != MVAvals1.second[2]   )     ){

//This is a type-3 combination 
MSPlot["MVA1TriJet_3"]->Fill(MVAvals1.first   ,  datasets[d], true, Luminosity*scaleFactor );
histo1D["BDT_Comb_type3"]->Fill(MVAvals1.first);
}
	//now look at the 2nd highest scoring combinatons, these should be more enhanced in wrong combinations

else if (   ( hadronicBJet_.first != MVAvals2ndPass.second[0] || hadronicBJet_.first != MVAvals2ndPass.second[1] || hadronicBJet_.first != MVAvals2ndPass.second[2]   )  && ( hadronicWJet1_.first == MVAvals2ndPass.second[0] || hadronicWJet1_.first == MVAvals2ndPass.second[1] || hadronicWJet1_.first == MVAvals2ndPass.second[2]   )    && ( hadronicWJet2_.first == MVAvals2ndPass.second[0] || hadronicWJet2_.first == MVAvals2ndPass.second[1] || hadronicWJet2_.first == MVAvals2ndPass.second[2]   )   ){

	  //This is a type-1A combination 
MSPlot["MVA1TriJet_1A"]->Fill(MVAvals2ndPass.first   ,  datasets[d], true, Luminosity*scaleFactor );
histo1D["BDT_Comb_type1A"]->Fill(MVAvals2ndPass.first);


	}else if (  ( hadronicBJet_.first == MVAvals2ndPass.second[0] || hadronicBJet_.first == MVAvals2ndPass.second[1] || hadronicBJet_.first == MVAvals2ndPass.second[2]   )  && ( hadronicWJet1_.first != MVAvals2ndPass.second[0] || hadronicWJet1_.first != MVAvals2ndPass.second[1] || hadronicWJet1_.first != MVAvals2ndPass.second[2]   )    && ( hadronicWJet2_.first == MVAvals2ndPass.second[0] || hadronicWJet2_.first == MVAvals2ndPass.second[1] || hadronicWJet2_.first == MVAvals2ndPass.second[2]   )          ){

 //This is a type-1B combination 

MSPlot["MVA1TriJet_1B"]->Fill(MVAvals2ndPass.first   ,  datasets[d], true, Luminosity*scaleFactor );

}
	else if ( ( hadronicBJet_.first == MVAvals2ndPass.second[0] || hadronicBJet_.first == MVAvals2ndPass.second[1] || hadronicBJet_.first == MVAvals2ndPass.second[2]   )  && ( hadronicWJet1_.first != MVAvals2ndPass.second[0] || hadronicWJet1_.first != MVAvals2ndPass.second[1] || hadronicWJet1_.first != MVAvals2ndPass.second[2]   )    && ( hadronicWJet2_.first != MVAvals2ndPass.second[0] || hadronicWJet2_.first != MVAvals2ndPass.second[1] || hadronicWJet2_.first != MVAvals2ndPass.second[2]   )     ){

//This is a type-2A combination 
MSPlot["MVA1TriJet_2A"]->Fill(MVAvals2ndPass.first   ,  datasets[d], true, Luminosity*scaleFactor );
histo1D["BDT_Comb_type2A"]->Fill(MVAvals2ndPass.first);

	}

else if (   ( hadronicBJet_.first != MVAvals2ndPass.second[0] || hadronicBJet_.first != MVAvals2ndPass.second[1] || hadronicBJet_.first != MVAvals2ndPass.second[2]   )  && ( hadronicWJet1_.first != MVAvals2ndPass.second[0] || hadronicWJet1_.first != MVAvals2ndPass.second[1] || hadronicWJet1_.first != MVAvals2ndPass.second[2]   )    && ( hadronicWJet2_.first == MVAvals2ndPass.second[0] || hadronicWJet2_.first == MVAvals2ndPass.second[1] || hadronicWJet2_.first == MVAvals2ndPass.second[2]   )         ){

//This is a type-2B combination 
MSPlot["MVA1TriJet_2B"]->Fill(MVAvals2ndPass.first   ,  datasets[d], true, Luminosity*scaleFactor );
histo1D["BDT_Comb_type2B"]->Fill(MVAvals2ndPass.first);
	}else if ( ( hadronicBJet_.first != MVAvals2ndPass.second[0] || hadronicBJet_.first != MVAvals2ndPass.second[1] || hadronicBJet_.first != MVAvals2ndPass.second[2]   )  && ( hadronicWJet1_.first != MVAvals2ndPass.second[0] || hadronicWJet1_.first != MVAvals2ndPass.second[1] || hadronicWJet1_.first != MVAvals2ndPass.second[2]   )    && ( hadronicWJet2_.first != MVAvals2ndPass.second[0] || hadronicWJet2_.first != MVAvals2ndPass.second[1] || hadronicWJet2_.first != MVAvals2ndPass.second[2]   )     ){

//This is a type-3 combination 
MSPlot["MVA1TriJet_3"]->Fill(MVAvals2ndPass.first   ,  datasets[d], true, Luminosity*scaleFactor );
histo1D["BDT_Comb_type3"]->Fill(MVAvals2ndPass.first);

}

*/
	//END OF OPTIONAL ANLAYIS OF COMBINATORICS

	   /*
	   
double unmatchedMass_C1A =( *selectedJets[badjet1] + *selectedJets[hadronicWJet1_.first] + *selectedJets[hadronicWJet2_.first]).M();
double unmatchedMass_C1B =( *selectedJets[hadronicBJet_.first] + *selectedJets[badjet2] + *selectedJets[hadronicWJet2_.first]).M();

double unmatchedMass_C2A =( *selectedJets[hadronicBJet_.first] + *selectedJets[badjet2] + *selectedJets[badjet2]).M();
double unmatchedMass_C2B =( *selectedJets[badjet1] + *selectedJets[badjet2] + *selectedJets[hadronicWJet2_.first]).M();


	   */

/*
            if (debug)  cout <<"checked mva success..."<< endl;
            double matchedMass, unmatchedMass;
            
            // Pick some random jets to investigate wrong combinations...
            int badjet1 =-1 , badjet2 =-1 , badjet3=-1;
            
	    vector<int> badjets;
	    badjets.clear();

           while((badjet1==-1)||(badjet2==-1)||(badjet3==-1)) {
                Int_t rand_jet_index  = rand->Integer(selectedJets.size());
                if ((rand_jet_index != hadronicBJet_.first) ||(rand_jet_index != hadronicWJet1_.first ) || (rand_jet_index != hadronicWJet2_.first)){

                    if (badjet1== -1) {
                        badjet1 = rand_jet_index;
			badjets.push_back(badjet1);
                        continue;
                    }
                    if ((badjet2== -1) && (rand_jet_index!= badjet1)) {
                        badjet2 = rand_jet_index;
			badjets.push_back(badjet2);
                        continue;
    
                    }
                    if ((badjet3== -1)  && (rand_jet_index!= badjet1) &&( rand_jet_index!= badjet2)) {
                        badjet3 = rand_jet_index;
			badjets.push_back(badjet3);
                        continue;
   
                    }
                }
            } 

	   double delR = 999.9;

	   double dijetmass_unmatched;
	   for(int bj1 = 0; bj1 < badjets.size(); bj1++ ){
                   for(int bj2 = 0; bj2 < badjets.size(); bj2++ ){

		     if (bj1 == bj2) continue;
		     TLorentzVector bj1_tlv = *selectedJets[badjets[bj1]];
		     TLorentzVector bj2_tlv = *selectedJets[badjets[bj2]];

                     double  delR_temp =   bj1_tlv.DeltaR(bj2_tlv);
		     if (delR_temp < delR) {
		       delR = delR_temp;
                       dijetmass_unmatched =  ( *selectedJets[badjets[bj1]] + *selectedJets[badjets[bj2]]).M();
     
		     }
                     }
                   }


            unmatchedMass =( *selectedJets[badjet1] + *selectedJets[badjet2] + *selectedJets[badjet3]).M();

	   

	  double topmass = 176.0;
	  double Wmass = 80.9;
          double trijetmass_unmatched  =( *selectedJets[badjet1] + *selectedJets[badjet2] + *selectedJets[badjet3]).M();
	  // double dijetmass_unmatched  =( *selectedJets[badjet1] + *selectedJets[badjet2]).M();
          double chi2_unmatched = (pow( dijetmass_unmatched - Wmass,2)/pow(10.38,2)) + (pow(trijetmass_unmatched - topmass,2)/pow(16.95,2));

          MSPlot["Chi2_UnMatched"]->Fill(chi2_unmatched,  datasets[d], true, Luminosity*scaleFactor );

          MSPlot["TriJetMass_UnMatched"]->Fill(unmatchedMass,  datasets[d], true, Luminosity*scaleFactor );

           // cout <<"Indices of matched jets are :  "<< hadronicBJet_.first<<"  "<< hadronicWJet1_.first  <<" " << hadronicWJet2_.first <<endl;
           // cout <<"Indices of random  jets are :  "<<badjet1 <<"  "<< badjet2  <<" " << badjet3  <<endl;
            //cout<<" "<<endl;
        
        if(   ( hadronicBJet_.first != 9999 )  && ( hadronicWJet1_.first != 9999   )    && ( hadronicWJet2_.first != 9999  )      ){
            nMatchedEvents++;
            if (debug) cout <<"matched..." <<endl;
            
            matchedMass =( *selectedJets[hadronicBJet_.first] + *selectedJets[hadronicWJet1_.first] + *selectedJets[hadronicWJet2_.first]).M();
          double trijetmass_matched  = ( *selectedJets[hadronicBJet_.first] + *selectedJets[hadronicWJet1_.first] + *selectedJets[hadronicWJet2_.first]).M();
          double dijetmass_matched  =  ( *selectedJets[hadronicWJet1_.first] + *selectedJets[hadronicWJet2_.first]).M();

          double chi2_matched = (pow( dijetmass_matched - Wmass,2)/pow(10.38,2)) + (pow(trijetmass_matched - topmass,2)/pow(16.95,2));

          MSPlot["Chi2_Matched"]->Fill(chi2_matched,  datasets[d], true, Luminosity*scaleFactor );

            MSPlot["MVA1TriJetMassMatched"]->Fill(bestTopMass1,  datasets[d], true, Luminosity*scaleFactor );
            MSPlot["TriJetMass_Matched"]->Fill(matchedMass,  datasets[d], true, Luminosity*scaleFactor );
        }
*/
        }
        
//end of optional checks of good/bad combinations
	       
////////////////////////////////////////////////////////////////////////
//          Plotting jet and event-level variables for rejecting ttbar + X
////////////////////////////////////////////////////////////////////////
        
        if (debug) cout <<"Caculating event level variables...  "<< endl;

	///        MEzCalculator NuPz;
        //NuPz.SetMET(*mets[0]);
	// NuPz.SetMuon(*selectedMuons[0]);
        //if (debug) cout <<"Created NuPz object "<< endl;
        
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
        


	//        sort(selectedJets2ndPass.begin(),selectedJets2ndPass.end(),HighestCVSBtag()); //order muons wrt dsicriminator

        if (debug) cout <<"Creating sumjets "<< endl;

	//        TRootJet sumjet (TLorentzVector (sumpx, sumpy, sumpz,sume )); //Object representing all the jets summed
        TRootJet sumjet_X (TLorentzVector (sumpx_X, sumpy_X, sumpz_X,sume_X )); //Object representing all the jets summed minus the hadronic system

        //now extend sumjet object to form event and reduced event
        //1. add muon
	//  sumpx_X = sumpx_X + selectedMuons[0]->Px();
        //sumpy_X = sumpy_X + selectedMuons[0]->Py();
        //sumpz_X = sumpz_X + selectedMuons[0]->Pz();
        //sume_X  = sume_X  + selectedMuons[0]->E();
        
        //2. add MET
        //sumpx_X = sumpx_X + mets[0]->Px();
	// sumpy_X = sumpy_X + mets[0]->Py();
        //sumpz_X = sumpz_X + NuPz.Calculate();
        //sume_X = sume_X + mets[0]->E();
        
	//   
        //now extend sumjet object to form event
        //1. add muon
	// sumpx = sumpx + selectedMuons[0]->Px();
        //sumpy = sumpy + selectedMuons[0]->Py();
        //sumpz = sumpz + selectedMuons[0]->Pz();
        //sume  = sume  + selectedMuons[0]->E();
        
        //2. add MET
	// sumpx = sumpx + mets[0]->Px();
        //sumpy = sumpy + mets[0]->Py();
        //sumpz = sumpz + NuPz.Calculate();
        //sume = sume + mets[0]->E();
        //double HT_top =1.;
        //TLorentzVector bestMVATop;
        //double ang_leptop = 0.;
        
	// if(!TrainMVA){ bestMVATop =( *selectedJets[MVAvals1.second[0]] + *selectedJets[MVAvals1.second[1]] + *selectedJets[MVAvals1.second[2]]);
        //    HT_top = selectedJets[MVAvals1.second[0]]->Pt() + selectedJets[MVAvals1.second[1]]->Pt() + selectedJets[MVAvals1.second[2]]->Pt();
        //    ang_leptop = bestMVATop.DeltaPhi(*selectedMuons[0]);  
	//        }

	//        TRootJet Event (TLorentzVector (sumpx, sumpy, sumpz,sume )); //Object representing the whole event
	//
	// TRootJet Event_X (TLorentzVector (sumpx_X, sumpy_X, sumpz_X,sume_X )); //Object representing the event minus the hadronic top

	//        EventMass = Event.M();
	// EventMassX = Event_X.M();
	// SumJetMass = sumjet.M();
        SumJetMassX = sumjet_X.M();
        
	// PTBalTopEventX = bestMVATop.Pt() - Event_X.Pt()  ;
	// PTBalTopSumJetX = bestMVATop.Pt() - sumjet_X.Pt()  ;
        
        HTH = HT/H;
        HTXHX = HTX/HX;
        
        //double HT_Asym = (HT_top - HT )/(HT_top + HT);
        
	/*
        //Form Lep W

         TLorentzVector mumet;      
         double   mumpx =   selectedMuons[0]->Px() + mets[0]->Px();
         double   mumpy =   selectedMuons[0]->Py() + mets[0]->Py();
         double   mumpz =   selectedMuons[0]->Pz() + NuPz.Calculate();
         double   mumpe =   selectedMuons[0]->E()  + mets[0]->E();
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
	*/


        //cout <<"  mumetB mass  "<<  muMETB.M()   << "  MVA Tops Mass "<< bestMVATop.M()   <<  "  difference : "<< deltaMTopMuMetB <<endl;
        
        //cout<<" H  = "<< H  << " HT  = "<< HT <<" EventMass =   "<< EventMass<< " Event mass X  "<< EventMassX <<  "SumJet mass " <<  SumJetMass<< "SumJet mass X " << SumJetMassX  <<endl;
        
        if (debug) cout <<"arrays  "<<endl;
        
	//        MHT = sumjet.Pt();
	// MHTSig = MHT/sqrt(HT);
        //STJet = HT + mets[0]->Et() + selectedMuons[0]->Pt();
        //EventMass = sumjet.M();
    
	//  sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag()); //order muons wrt dsicriminator
	//        sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order muons wrt dsicriminator
        
        if (debug) cout <<"jets sorted...  "<<endl;

	/*
	histoProfile[("HTX_vs_MultiTopness_"+datasets[d]->Name()).c_str()]->Fill(MultiTopness, HTX);
        histoProfile[("HTX_vs_HTH_"+datasets[d]->Name()).c_str()]->Fill(HTH,HTX);
        histoProfile[("HTX_vs_5thJetPt_"+datasets[d]->Name()).c_str()]->Fill(selectedJets[4]->Pt(),HTX);
        histoProfile[("HTX_vs_6thJetPt_"+datasets[d]->Name()).c_str()]->Fill(selectedJets[5]->Pt(),HTX);
        histoProfile[("HTX_vs_HTb_"+datasets[d]->Name()).c_str()]->Fill(HTb,HTX);
        histoProfile[("HTX_vs_HTRat_"+datasets[d]->Name()).c_str()]->Fill(HTRat,HTX);
        histoProfile[("HTX_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str()]->Fill(selectedJets.size(),HTX);        
        histoProfile[("HTX_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str()]->Fill(ntags,HTX);        
        histoProfile[("HTX_vs_SumJetMassX_"+datasets[d]->Name()).c_str()]->Fill(SumJetMassX,HTX);        

        histoProfile[("MultiTopness_vs_HTH_"+datasets[d]->Name()).c_str()]->Fill(HTH,MultiTopness);
        histoProfile[("MultiTopness_vs_5thJetPt_"+datasets[d]->Name()).c_str()]->Fill(selectedJets[4]->Pt(),MultiTopness);
        histoProfile[("MultiTopness_vs_6thJetPt_"+datasets[d]->Name()).c_str()]->Fill(selectedJets[5]->Pt(),MultiTopness);
        histoProfile[("MultiTopness_vs_HTb_"+datasets[d]->Name()).c_str()]->Fill(HTb,MultiTopness);
        histoProfile[("MultiTopness_vs_HTRat_"+datasets[d]->Name()).c_str()]->Fill(HTRat,MultiTopness);
        histoProfile[("MultiTopness_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str()]->Fill(selectedJets.size(),MultiTopness);
        histoProfile[("MultiTopness_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str()]->Fill(ntags,MultiTopness);
        histoProfile[("MultiTopness_vs_SumJetMassX_"+datasets[d]->Name()).c_str()]->Fill(SumJetMassX,MultiTopness);

        histoProfile[("HTH_vs_5thJetPt_"+datasets[d]->Name()).c_str()]->Fill(selectedJets[4]->Pt(),HTH);
        histoProfile[("HTH_vs_6thJetPt_"+datasets[d]->Name()).c_str()]->Fill(selectedJets[5]->Pt(),HTH);
        histoProfile[("HTH_vs_HTb_"+datasets[d]->Name()).c_str()]->Fill(HTb,HTH);
        histoProfile[("HTH_vs_HTRat_"+datasets[d]->Name()).c_str()]->Fill(HTRat,HTH);
        histoProfile[("HTH_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str()]->Fill(selectedJets.size(),HTH);
        histoProfile[("HTH_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str()]->Fill(ntags,HTH);
        histoProfile[("HTH_vs_SumJetMassX_"+datasets[d]->Name()).c_str()]->Fill(SumJetMassX ,HTH);


        histoProfile[("5thJetPt_vs_6thJetPt_"+datasets[d]->Name()).c_str()]->Fill(selectedJets[5]->Pt(),selectedJets[4]->Pt());
        histoProfile[("5thJetPt_vs_HTb_"+datasets[d]->Name()).c_str()]->Fill(HTb,selectedJets[4]->Pt());
        histoProfile[("5thJetPt_vs_HTRat_"+datasets[d]->Name()).c_str()]->Fill(HTRat,selectedJets[4]->Pt());
        histoProfile[("5thJetPt_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str()]->Fill(selectedJets.size() ,selectedJets[4]->Pt());
        histoProfile[("5thJetPt_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str()]->Fill(ntags ,selectedJets[4]->Pt());
        histoProfile[("5thJetPt_vs_SumJetMassX_"+datasets[d]->Name()).c_str()]->Fill(SumJetMassX ,selectedJets[4]->Pt());


        histoProfile[("6thJetPt_vs_HTb_"+datasets[d]->Name()).c_str()]->Fill(HTb,selectedJets[5]->Pt());
        histoProfile[("6thJetPt_vs_HTRat_"+datasets[d]->Name()).c_str()]->Fill(HTRat,selectedJets[5]->Pt());
        histoProfile[("6thJetPt_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str()]->Fill(selectedJets.size(),selectedJets[5]->Pt());
        histoProfile[("6thJetPt_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str()]->Fill(ntags,selectedJets[5]->Pt());
        histoProfile[("6thJetPt_vs_SumJetMassX_"+datasets[d]->Name()).c_str()]->Fill(SumJetMassX ,selectedJets[5]->Pt());


        histoProfile[("HTb_vs_HTRat_"+datasets[d]->Name()).c_str()]->Fill(HTRat,HTb);
        histoProfile[("HTb_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str()]->Fill(selectedJets.size(),HTb);
        histoProfile[("HTb_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str()]->Fill(ntags,HTb);
        histoProfile[("HTb_vs_SumJetMassX_"+datasets[d]->Name()).c_str()]->Fill(SumJetMassX,HTb);


        histoProfile[("HTRat_vs_NbOfSelectedJets_"+datasets[d]->Name()).c_str()]->Fill(selectedJets.size(),HTRat);
        histoProfile[("HTRat_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str()]->Fill(ntags,HTRat);
        histoProfile[("HTRat_vs_SumJetMassX_"+datasets[d]->Name()).c_str()]->Fill(SumJetMassX,HTRat);

        histoProfile[("NbOfSelectedJets_vs_NbOfSelectedBJets_"+datasets[d]->Name()).c_str()]->Fill(ntags,selectedJets.size());
        histoProfile[("NbOfSelectedJets_vs_SumJetMassX_"+datasets[d]->Name()).c_str()]->Fill(SumJetMassX,selectedJets.size());

        histoProfile[("NbOfSelectedBJets_vs_SumJetMassX_"+datasets[d]->Name()).c_str()]->Fill(SumJetMassX,ntags);


	*/


	float temp_weight = 1.0;
        if(trainEventMVA){
        
        if(dataSetName=="NP_overlay_TTTT"){

        Eventtrainer_->FillWeight("S","Weight",scaleFactor);
            
	//  	cout <<"scalefactor  = = "<< scaleFactor << endl;
        Eventtrainer_->Fill("S","HTX", HTX);
        Eventtrainer_->Fill("S","HTH", HTH);
        Eventtrainer_->Fill("S","HTb", HTb);
        Eventtrainer_->Fill("S","HTRat", HTRat);
        Eventtrainer_->Fill("S","SumJetMassX", SumJetMassX);
        Eventtrainer_->Fill("S","MultiTopness", MultiTopness);
        Eventtrainer_->Fill("S","nTags",ntags  );
        Eventtrainer_->Fill("S","nJets",selectedJets.size()  );
        Eventtrainer_->Fill("S","Jet5Pt",selectedJets[4]->Pt());
        Eventtrainer_->Fill("S","Jet6Pt",selectedJets[5]->Pt());
        }
        
        if(dataSetName=="TTJets"){

	   Eventtrainer_->FillWeight("B","Weight", scaleFactor);

            Eventtrainer_->Fill("B","HTX", HTX);
            Eventtrainer_->Fill("B","HTH", HTH);
            Eventtrainer_->Fill("B","HTb", HTb);
            Eventtrainer_->Fill("B","HTRat", HTRat);
            Eventtrainer_->Fill("B","SumJetMassX", SumJetMassX);
            Eventtrainer_->Fill("B","MultiTopness", MultiTopness);
            Eventtrainer_->Fill("B","nTags", ntags );
	    Eventtrainer_->Fill("B","nJets", selectedJets.size() );
            Eventtrainer_->Fill("B","Jet5Pt", selectedJets[4]->Pt() );
            Eventtrainer_->Fill("B","Jet6Pt", selectedJets[5]->Pt() );        
            
        }
            
        }
        
        else if (computeEventMVA){
            if (debug) cout <<"filling computer...."<<H <<endl;

            if (Eventcomputer_ == 0) cout <<"null computer...."<<H <<endl;

           
            
        if (debug) cout <<"filling computer...."<<endl;

        Eventcomputer_->FillVar("HTX", HTX);
        Eventcomputer_->FillVar("HTH", HTH);
        Eventcomputer_->FillVar("HTb", HTb);
        Eventcomputer_->FillVar("HTRat",HTRat);
        Eventcomputer_->FillVar("SumJetMassX", SumJetMassX);
        Eventcomputer_->FillVar("MultiTopness", MultiTopness);
        Eventcomputer_->FillVar("nTags",ntags );     
        Eventcomputer_->FillVar("nJets", selectedJets.size() );
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
       

	//testing super-lite TNtuples (to be compared with ttree)
      tup->Fill(BDTscore,njets,ntags,MultiTopness,HTb,HTH,HTRat,HTX,SumJetMassX,selectedJets[4]->Pt(),selectedJets[5]->Pt(),scaleFactor,datasets[d]->NormFactor(),Luminosity);
      //   cout <<BDTscore <<" " <<njets <<" "  <<ntags <<" " <<MultiTopness  <<" " << HTb <<" " << HTRat <<" "    << HTX <<  SumJetMassX  << "  "     <<selectedJets[4]->Pt()   <<  " " << selectedJets[5]->Pt()  <<endl;

	if((BDTscore > 0.2) && (dataSetName=="Data")){
	  eventlist <<currentRun  << " " << event->lumiBlockId() <<" " <<event->eventId() << " "  << BDTscore <<" N Jets = " << selectedJets.size()<<" N Tags = " << ntags<< " MultiTopness => "  <<MultiTopness<< " HTX = >  "<<HTX<< " HTRat = > "<< HTRat << " HTb =>  "<< HTb  << endl;        

	}

	//	if((BDTscore > 0.05)){
	  //    MSPlot["NbOfSelectedJets_BDTCut"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);                
	  //MSPlot["NbOfSelectedBJets_BDTCut"]->Fill(ntags, datasets[d], true, Luminosity*scaleFactor);             
	  //MSPlot["HTX_BDTCut"]->Fill(HTX, datasets[d], true, Luminosity*scaleFactor); 
	  //MSPlot["HTH_BDTCut"]->Fill(HTH, datasets[d], true, Luminosity*scaleFactor); 
	  //MSPlot["MultiTopness_BDTCut"]->Fill(MultiTopness, datasets[d], true, Luminosity*scaleFactor); 
	  //MSPlot["HTb_SelectedJets_BDTCut"]->Fill(HTb, datasets[d], true, Luminosity*scaleFactor); 
	  //}

        if (debug) cout <<"filling event level histos  "<<endl;

        
	// MSPlot["EventMass"]->Fill(EventMass, datasets[d], true, Luminosity*scaleFactor);
	// MSPlot["EventMassX"]->Fill(EventMassX, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["SumJetMass"]->Fill(SumJetMass, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["SumJetMassX"]->Fill(SumJetMassX, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["HTX"]->Fill(HTX, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["H"]->Fill(H, datasets[d], true, Luminosity*scaleFactor);
	// MSPlot["HX"]->Fill(HX, datasets[d], true, Luminosity*scaleFactor);
	//        MSPlot["HTXHX"]->Fill(HTXHX, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["HTH"]->Fill(HTH, datasets[d], true, Luminosity*scaleFactor);       
        MSPlot["HT_SelectedJets"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["HTb_SelectedJets"]->Fill(HTb, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["HTRat"]->Fill(HTRat, datasets[d], true, Luminosity*scaleFactor);
        MSPlot["MET"]->Fill(mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor);
        
        if (debug) cout <<"filled event level histos  "<<endl;
        

//////////////////////////////////////////////////////////////////////////////////
/////          Filling NJet Vs NBJet arrays for                             //////
/////                                                                       //////
///// 1. HT, 2. STLep, 3.STJet, 4.HT(w/o ttbar), 5.InvMass (w/o) ttbar      //////
//////////////////////////////////////////////////////////////////////////////////
	Int_t b =0;
//  for (Int_t b = minNJets; b <= maxNJets; b++){
      for (Int_t c = minNJets; c<= maxNJets; c++){
          string NJets_str = static_cast<ostringstream*>( &(ostringstream() << c) )->str();
	  // string NBJets_str = static_cast<ostringstream*>( &(ostringstream() << c) )->str();
	  // string HT_Name = "HT_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string HTX_Name = "HTX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string H_Name = "H_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string HX_Name = "HX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string HTH_Name = "HTH_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string HTXHX_Name = "HTXHX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
	  //          string EventMass_Name = "EventMass_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
	  // string EventMassX_Name = "EventMassX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string FirstTriJetMass_1BJet_g_Name = "FirstTriJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string FirstDiJetMass_1BJet_g_Name = "FirstDiJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string SecondDiJetMass_1BJet_g_Name = "SecondDiJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;

	  // string SumJetMass_Name = "SumJetMass_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string SumJetMassX_Name = "SumJetMassX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string PTBalTopEventX_Name = "PTBalTopEventX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string PTBalTopSumJetX_Name = "PTBalTopSumJetX_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;

	  // string MuMetBMasses_g_Name = "MuMetBMasses_g_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string MuMetMasses_g_Name = "MuMetMasses_g_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
	  // string MET_Name = "MET"+NJets_str+"Jets_"+NBJets_str+"Tags" ;

          //string SecondTriJetMass_1BJet_g_Name = "SecondTriJetMass_1BJet_g"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string SecondTriJetMass_1BJet_g_chi2_Name = "SecondTriJetMass_1BJet_g_chi2_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string MVA1TriJetMass_Name = "MVA1TriJetMass"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
          //string MVA2TriJetMass_Name = "MVA2TriJetMass"+NJets_str+"Jets_"+NBJets_str+"Tags" ;

          string MVA_Name = "MVA"+NJets_str+"Jets" ;

          
          
        if(c<=7){
            if(selectedJets.size() == c  ) {

	      //	      cout <<"filling "<< c <<endl;
	      //    MSPlot[HT_Name.c_str() ]->Fill(HT,datasets[d], true, Luminosity*scaleFactor);
              //  MSPlot[HTX_Name.c_str() ]->Fill(HTX,datasets[d], true, Luminosity*scaleFactor);
              //  MSPlot[H_Name.c_str() ]->Fill(H,datasets[d], true, Luminosity*scaleFactor);
              //  MSPlot[HX_Name.c_str() ]->Fill(HX,datasets[d], true, Luminosity*scaleFactor);
              //  MSPlot[HTH_Name.c_str() ]->Fill(HTH,datasets[d], true, Luminosity*scaleFactor);
              //  MSPlot[HTXHX_Name.c_str() ]->Fill(HTXHX,datasets[d], true, Luminosity*scaleFactor);
              //  MSPlot[EventMass_Name.c_str() ]->Fill(EventMass,datasets[d], true, Luminosity*scaleFactor);
              //  MSPlot[EventMassX_Name.c_str() ]->Fill(EventMassX,datasets[d], true, Luminosity*scaleFactor);
              //  MSPlot[SumJetMass_Name.c_str() ]->Fill(SumJetMass,datasets[d], true, Luminosity*scaleFactor);
              //  MSPlot[SumJetMassX_Name.c_str() ]->Fill(SumJetMassX,datasets[d], true, Luminosity*scaleFactor);
		//    MSPlot[PTBalTopEventX_Name.c_str() ]->Fill(PTBalTopEventX,datasets[d], true, Luminosity*scaleFactor);
		// MSPlot[PTBalTopSumJetX_Name.c_str() ]->Fill(PTBalTopSumJetX,datasets[d], true, Luminosity*scaleFactor);
            
	      // MSPlot[MET_Name.c_str()]->Fill(mets[0]->Et() ,  datasets[d], true, Luminosity*scaleFactor );       
	      // MSPlot[MVA1TriJetMass_Name.c_str()]->Fill(bestTopMass1,  datasets[d], true, Luminosity*scaleFactor);
	      //MSPlot[MVA2TriJetMass_Name.c_str()]->Fill(bestTopMass2,  datasets[d], true, Luminosity*scaleFactor);
        MSPlot[MVA_Name.c_str()]->Fill(BDTscore,  datasets[d], true, Luminosity*scaleFactor );


 }
                        }
        else if ( c > 7) {
	  if( selectedJets.size() >= c  ) {
	      // MSPlot[HT_Name.c_str() ]->Fill(HT,datasets[d], true, Luminosity*scaleFactor);
	      // MSPlot[HTX_Name.c_str() ]->Fill(HTX,datasets[d], true, Luminosity*scaleFactor);
              //  MSPlot[H_Name.c_str() ]->Fill(H,datasets[d], true, Luminosity*scaleFactor);
              //  MSPlot[HX_Name.c_str() ]->Fill(HX,datasets[d], true, Luminosity*scaleFactor);
              //  MSPlot[HTH_Name.c_str() ]->Fill(HTH,datasets[d], true, Luminosity*scaleFactor);
              //  MSPlot[HTXHX_Name.c_str() ]->Fill(HTXHX,datasets[d], true, Luminosity*scaleFactor);
                
	      //MSPlot[EventMass_Name.c_str() ]->Fill(EventMass,datasets[d], true, Luminosity*scaleFactor);
	      // MSPlot[EventMassX_Name.c_str() ]->Fill(EventMassX,datasets[d], true, Luminosity*scaleFactor);

              //  MSPlot[MET_Name.c_str()]->Fill(mets[0]->Et() ,  datasets[d], true, Luminosity*scaleFactor );
                
              //  MSPlot[SumJetMass_Name.c_str() ]->Fill(SumJetMass,datasets[d], true, Luminosity*scaleFactor);
	      // MSPlot[SumJetMassX_Name.c_str() ]->Fill(SumJetMassX,datasets[d], true, Luminosity*scaleFactor);
		// MSPlot[PTBalTopEventX_Name.c_str() ]->Fill(PTBalTopEventX,datasets[d], true, Luminosity*scaleFactor);
                //MSPlot[PTBalTopSumJetX_Name.c_str() ]->Fill(PTBalTopSumJetX,datasets[d], true, Luminosity*scaleFactor);
                
      
                //MSPlot[MVA1TriJetMass_Name.c_str()]->Fill(bestTopMass1,  datasets[d], true, Luminosity*scaleFactor);
                //MSPlot[MVA2TriJetMass_Name.c_str()]->Fill(bestTopMass2,  datasets[d], true, Luminosity*scaleFactor);

                MSPlot[MVA_Name.c_str()]->Fill(BDTscore,  datasets[d], true, Luminosity*scaleFactor );

		//          }
        }
}
}
 


  for (Int_t b = minNBJets; b <= maxNBJets; b++){

       string NBJets_str = static_cast<ostringstream*>( &(ostringstream() << b) )->str();
       string MVA_Name = "MVA"+NBJets_str+"Tags" ;
          
        if(b<=2){
            if(ntags == b  ) {
        MSPlot[MVA_Name.c_str()]->Fill(BDTscore,  datasets[d], true, Luminosity*scaleFactor );

 }
                        }
        else if ( b > 2) {
	  if( ntags >= b  ) {
	
                MSPlot[MVA_Name.c_str()]->Fill(BDTscore,  datasets[d], true, Luminosity*scaleFactor );

        }
}

}

  
	}


      
	//	cout <<"clearing b jets vector..."<< ntags  << endl;


        if (debug) cout <<"freeing memory  "<<endl;





	/*
      selectedElectrons.clear();     
      selectedExtraElectrons.clear();
      selectedMuons_NoIso.clear(); 
      selectedExtraMuons.clear();
      selectedJets.clear();  
      selectedJets2ndPass.clear();
      MVASelJets1.clear();
      //   selectedJets3rdPass.clear();
      // selectedElectrons_NoIso.clear();
      selectedMuons.clear();  
      //selectedSoftJets.clear();
      selectedLooseMuons.clear();
      selectedLooseElectrons.clear();
      selectedBJets.clear();
      selectedLightJets.clear();
*/
      //selectedCSVOrderedJets.clear();
      // if (genEvt) delete genEvt;
	// if(genEvt_flav) delete genEvt_flav;

 
	// if (BTagGraph_B) delete BTagGraph_B; 
	// if (BTagGraph_C) delete BTagGraph_C; 

	//if (BTagGraph_L_eta008)  delete BTagGraph_L_eta008;
	// if (BTagGraph_L_eta0816) delete BTagGraph_L_eta0816;
	//if (BTagGraph_L_eta1624) delete BTagGraph_L_eta1624;

      //      if( btageffFile) delete  btageffFile;

        if (debug) cout <<" memory  freed "<<endl;



      




    }//loop on events
    


     tup->Write();
     tupfile->Close();



    cout <<"n events passed  =  "<<passed <<endl;

    if (debug)cout <<"N BBar = = " << nBBBar <<"  N CCBar = = " << nCCBar <<"  N LLBar = =  " << nLLBar << endl;
      
    
    // if(jetTools) delete jetTools;
    
    //    if(L1JetCorPar) delete L1JetCorPar;
    // if(L2JetCorPar) delete L2JetCorPar;
    //if(L3JetCorPar) delete L3JetCorPar;
    // if(jecUnc) delete jecUnc;    

      //important: free memory
      treeLoader.UnLoadDataset();
      
      double nMVASuccessesd = nMVASuccesses;
      double nMatchedEventsd = nMatchedEvents;
      
      cout <<"Efficiency of MVA jet selection = = "<<  nMVASuccessesd/nMatchedEventsd   <<endl;
      
  } //loop on datasets
    
eventlist.close();

    if(trainEventMVA) Eventtrainer_->TrainMVA("Block","",0,0,"",0,0,"_25thFeb");
    
  //Once everything is filled ...
  cout << " We ran over all the data ;-)" << endl;
  
  ///////////////////
  // Writing
  //////////////////
  cout << " - Writing outputs to the files ..." << endl;

  //Selection tables
  if(Muon){
	//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
	selecTableMu.TableCalculator(  true, true, true, true, true);
      
  //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawsyNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
  	selecTableMu.Write(  "FourTop"+postfix+"Table_Mu.tex",    false,true,true,true,false,false,false);
  }
    else if(Electron){
	//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
	selecTableEl.TableCalculator(  false, true, true, true, true);
   //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
    selecTableEl.Write(  "FourTop"+postfix+"Table_El.tex",  true,true,true,true,false,false,true);

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
    
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
        string name = it->first;
        
        MultiSamplePlot *temp = it->second;
        TH1F *tempHisto_data;
        TH1F *tempHisto_TTTT;

      //Option to write ROOT files containing histograms with systematics varied +/-
      //TFile * tempErrorFile;
        if(doScaleShift == 1){
            cout <<"Scaling sys down"<<endl;

            string filename = "ScaleFilesMu/Error_"+name+".root";
            TFile* tempErrorFile = new TFile(filename.c_str(),"RECREATE");
            TH1F * tempHisto = temp->getTH1F("TTJets_ScaleDown");
            tempHisto->Write("Down");
            tempErrorFile->Write();
            tempErrorFile->Close();

	    //  TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
	    //	    TFile* sysFile; 


            if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
	    TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
             string hist_name =  name + "_Scale_Down";
	     //  tempHisto->Write(hist_name.c_str());
	     sysFile->Close();

	    }
            delete tempErrorFile, tempHisto;
        }else if(doScaleShift == 2){
            cout <<"Scaling sys up"<<endl;
            string filename = "ScaleFilesMu/Error_"+name+".root";
            TFile* tempErrorFile = new TFile(filename.c_str(),"UPDATE");
            TH1F * tempHisto = temp->getTH1F("TTJets_ScaleUp");
            tempHisto->Write("Up");
            tempErrorFile->Write();
            tempErrorFile->Close();
     if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
	      string hist_name =  name + "_Scale_Up";
	      // tempHisto->Write(hist_name.c_str());
             sysFile->Close();
	    }
            delete tempErrorFile, tempHisto;
        }
        else if(doScaleShift == 0){
	  //cout <<"Scaling sys nominal"<<endl;
        }
        

      if(doMatchingShift == 1){
            cout <<"Matching sys down"<<endl;

            string filename = "MatchingFiles/Error_"+name+".root";
            TFile* tempErrorFile = new TFile(filename.c_str(),"UPDATE");
            TH1F * tempHisto = temp->getTH1F("TTJets_MatchingDown");
            tempHisto->Write("Down");
            tempErrorFile->Write();
            tempErrorFile->Close();

   if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){

             TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
             string hist_name =  name + "_Matching_Down";
             tempHisto->Write(hist_name.c_str());
             sysFile->Close();
	    }
	    delete tempErrorFile, tempHisto;
        }else if(doMatchingShift == 2){
            cout <<"Matching sys up"<<endl;
            string filename = "MatchingFiles/Error_"+name+".root";
            TFile* tempErrorFile = new TFile(filename.c_str(),"UPDATE");
            TH1F * tempHisto = temp->getTH1F("TTJets_MatchingUp");
            tempHisto->Write("Up");
            tempErrorFile->Write();
            tempErrorFile->Close();
      if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
             string hist_name =  name + "_Matching_Up";             
             tempHisto->Write(hist_name.c_str());
             sysFile->Close();
	    }

            delete tempErrorFile, tempHisto;

        }
        else if(doMatchingShift == 0){
	  //cout <<"Matching sys nominal"<<endl;
        }


	if (dottbbShift == 1){
               if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TH1F * tempHisto = temp->getTH1F("TTJets");
             TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
             string hist_name =  name + "_ttbb_Down"; 
             tempHisto->Write(hist_name.c_str());
             sysFile->Close();
 }
}
	else if(dottbbShift == 2){
     if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TH1F * tempHisto = temp->getTH1F("TTJets");
             TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
             string hist_name =  name + "_ttbb_Up"; 
             tempHisto->Write(hist_name.c_str());
             sysFile->Close();
 }
}
	else if (dottbbShift == 0) {
         cout <<"ttbb sys nominal"<<endl;
}


	if (doPUShift == 1){
               if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TH1F * tempHisto = temp->getTH1F("TTJets");
             TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
             string hist_name =  name + "_PU_Down"; 
             tempHisto->Write(hist_name.c_str());
             sysFile->Close();
 }
}
	else if(doPUShift == 2){
     if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TH1F * tempHisto = temp->getTH1F("TTJets");
             TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
             string hist_name =  name + "_PU_Up"; 
             tempHisto->Write(hist_name.c_str());
             sysFile->Close();
 }
}
	else if (doPUShift == 0) {
         cout <<"PU sys nominal"<<endl;
}


    if (doLeptonSFShift == 1){


	     if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
	    TH1F * tempHisto = temp->getTH1F("TTJets");
            TH1F * tempHisto_tttt = temp->getTH1F("NP_overlay_TTTT");
	    TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
            string hist_name =  name + "_leptonSF_Down"; 
            string hist_name_tttt =  name + "_leptonSF_Down_tttt";
	    tempHisto->Write(hist_name.c_str());
      cout <<"writing histos...."<<endl;
            tempHisto_tttt->Write(hist_name_tttt.c_str());
            sysFile->Close();
	  }
	}
        else if(doLeptonSFShift == 2){
 cout <<"writing histos...."<<endl;
	      if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
	    TH1F * tempHisto = temp->getTH1F("TTJets");
            TH1F * tempHisto_tttt = temp->getTH1F("NP_overlay_TTTT");
	    TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
            string hist_name =  name + "_leptonSF_Up"; 
            string hist_name_tttt =  name + "_leptonSF_Up_tttt";
	    tempHisto->Write(hist_name.c_str());
            tempHisto_tttt->Write(hist_name_tttt.c_str());
            sysFile->Close();
	  }
	}
        else if (doLeptonSFShift == 0) {
	  cout <<"leptonSF sys nominal"<<endl;
	}


	if (doJERShift == 1){
              if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
            TH1F * tempHisto = temp->getTH1F("TTJets");
            TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
            string hist_name =  name + "_JER_Down"; 
            tempHisto->Write(hist_name.c_str());
             sysFile->Close();
          }
	}
        else if(doJERShift == 2){
              if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
            TH1F * tempHisto = temp->getTH1F("TTJets");
            TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
            string hist_name =  name + "_JER_Up"; 
            tempHisto->Write(hist_name.c_str());
             sysFile->Close();
          }
        }
        else if (doJERShift == 0) {
          cout <<"JER sys nominal"<<endl;
        }

	if (domisTagEffShift == -1){
               if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TH1F * tempHisto = temp->getTH1F("TTJets");
             TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
             string hist_name =  name + "_misTag_Down"; 
             tempHisto->Write(hist_name.c_str());
             sysFile->Close();
 }
}
	else if(domisTagEffShift == 1){
     if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TH1F * tempHisto = temp->getTH1F("TTJets");
             TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
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
             TH1F * tempHisto_tttt = temp->getTH1F("NP_overlay_TTTT");
             TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
             string hist_name =  name + "_bTag_Down";
             string hist_name_tttt =  name + "_leptonSF_Down_tttt"; 
             tempHisto->Write(hist_name.c_str());
             tempHisto_tttt->Write(hist_name_tttt.c_str());
             sysFile->Close();
 }
}
	else if(dobTagEffShift == 1){
     if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TH1F * tempHisto = temp->getTH1F("TTJets");
             TH1F * tempHisto_tttt = temp->getTH1F("NP_overlay_TTTT");
             TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
             string hist_name =  name + "_bTag_Up"; 
             string hist_name_tttt =  name + "_leptonSF_Down_tttt"; 
             tempHisto->Write(hist_name.c_str());
             tempHisto_tttt->Write(hist_name_tttt.c_str());
             sysFile->Close();
 }
}
	else if (dobTagEffShift == 0) {
         cout <<"bTag sys nominal"<<endl;
}

             
        if (doJESShift == 1){
            string filename = "JESFiles/Error_"+name+".root";
            TFile* tempErrorFile = new TFile(filename.c_str(),"UPDATE");
            TH1F * tempHisto = temp->getTH1F("TTJets");
            tempHisto->Write("Minus");
            tempErrorFile->Write();
            tempErrorFile->Close();
       if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
             TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
             string hist_name =  name + "_JES_Down";            
             tempHisto->Write(hist_name.c_str());
             sysFile->Close();
	    }
	    delete tempErrorFile, tempHisto;
        }
        else if  (doJESShift ==2){
            
            string filename = "JESFiles/Error_"+name+".root";
            TFile* tempErrorFile = new TFile(filename.c_str(),"UPDATE");
            TH1F * tempHisto = temp->getTH1F("TTJets");
            tempHisto->Write("Plus");
            tempErrorFile->Write();
            tempErrorFile->Close();
               if(name=="MVA" || name=="MVA6Jets" || name=="MVA7Jets" || name=="MVA8Jets" ){
            TFile* sysFile = new TFile("SystematicShapes_Mu.root","UPDATE");
            string hist_name =  name + "_JES_Up";            
            tempHisto->Write(hist_name.c_str());
             sysFile->Close();
	    }
	    delete tempErrorFile, tempHisto;

	    //            cout <<"JES sys down"<<endl;
            
        }
        else if  (doJESShift ==0){
	  // cout <<"JES sys off "<<endl;
      
   // temp->addText("CMS preliminary");
	  //	temp->Draw(false, name, true, true, true, true, true,100.,true); // merge TT/QCD/W/Z/ST/
	  temp->Draw_wSysUnc(false,"ScaleFilesMu", name, true, true, false, false, false,100,true, false, false, true); // merge TT/QCD/W/Z/ST/

            
            //Draw(bool addRandomPseudoData = false, string label = string("CMSPlot"), bool mergeTT = false, bool mergeQCD = false, bool mergeW = false, bool mergeZ = false, bool mergeST = false, int scaleNPSignal = 1, bool addRatio = false, bool mergeVV = false, bool mergeTTV = false);
	  //   cout <<" Writing plot:  name = "<< name<<"  path = "<< pathPNG  <<endl;
      temp->Write(fout, name, true, pathPNG, "pdf");
    }
  }

  TDirectory* th1dir = fout->mkdir("Histos1D");
  th1dir->cd();
 
  TH1F *ttJets_forCorr_Njets =   MSPlot["NbOfSelectedJets"]->getTH1F("TTJets");
  TH1F *data_forCorr_Njets =   MSPlot["NbOfSelectedJets"]->getTH1F("Data");
  TH1F *ttJets_forCorr_Ntags =   MSPlot["NbOfSelectedBJets"]->getTH1F("TTJets");
  TH1F *data_forCorr_Ntags =   MSPlot["NbOfSelectedBJets"]->getTH1F("Data");
    
  ttJets_forCorr_Njets->Write();
  ttJets_forCorr_Ntags->Write();
 
  cout <<" integral ttjets njets = =  "<< ttJets_forCorr_Njets->Integral() <<endl;
  cout <<" integral ttjets ntags = =  "<< ttJets_forCorr_Ntags->Integral() <<endl;


 for(map<std::string,TH1F*>::const_iterator it = histo1D.begin(); it != histo1D.end(); it++)
  {

	TH1F *temp = it->second;
	temp->Write();
	//TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
	//tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }
    
    double corr_factor=1.;
    double ttjets_exp = 1;
    double data_obs = 1;


   // for(int jet_bin = 1; jet_bin < 16; jet_bin++){
    //    for(int tag_bin = 1; tag_bin < 9; tag_bin++){
    //        corr_factor=1.
   //         ttjets_exp = ttJets_forCorr_Ntags->GetBinContent(tag_bin);
   //         data_obs = data_forCorr_Ntags->GetBinContent(tag_bin);
   //         if(ttjets_exp !=0) corr_factor = (data_obs)/(ttjets_exp);
   //         cout <<"corr factor == " << corr_factor << endl;
   //     }
  //  }
    
    
    
  TDirectory* th2dir = fout->mkdir("Histos2D");
  th2dir->cd();
  
    
   for(map<std::string,TH2F*>::const_iterator it = histo2D.begin(); it != histo2D.end(); it++)
  {

	TH2F *temp = it->second;
	temp->Write();
	//TCanvas* tempCanvas = TCanvasCreator(temp, it->first);
	//tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );
  }
    
   TDirectory* tProfiledir = fout->mkdir("HistosProfile");
   tProfiledir->cd();


   for(map<std::string,TProfile*>::const_iterator it = histoProfile.begin(); it != histoProfile.end(); it++)
     {

       TProfile *temp = it->second;
       temp->Write();
       TCanvas* tempCanvas = TCanvasCreator(temp, it->first);                                                  
       tempCanvas->SaveAs( (pathPNG+it->first+".png").c_str() );                                               
     }

   if (leptonTools) delete leptonTools;
   if(Eventtrainer_) delete Eventtrainer_;
   if(jetCombiner) delete jetCombiner;
   if (bTool) delete bTool;
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
