#include "TStyle.h"
#include "TPaveText.h"

#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include "TRandom3.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"


//user code
#include "TopTreeProducer/interface/TRootRun.h"
#include "TopTreeProducer/interface/TRootEvent.h"
#include "TopTreeAnalysisBase/Selection/interface/SelectionTable.h"
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
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
//#include "TopTreeAnalysisBase/InclFourthGenSearch/interface/InclFourthGenSearchTools.h"
#include "TopTreeAnalysisBase/MCInformation/interface/ResolutionFit.h"

#include "TopTreeAnalysisBase/Tools/interface/JetCombiner.h"
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"


#include "TLorentzVector.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/TTreeObservables.h"
//#include "interface/FourTopTree.h"
#include "../macros/Style.C"

using namespace std;
using namespace TopTree;
using namespace reweight;

float workingpointvalue = 0.679;

Int_t n_true_bjets = 0, n_tagged_true_bjets = 0, n_true_cjets = 0, n_tagged_true_cjets = 0, n_true_ljets = 0, n_tagged_true_ljets = 0;

string TreeFileName = "FourTopTree.root";

struct HighestCSVBtag{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const{
    	return j1->btag_combinedSecondaryVertexBJetTags() > j2->btag_combinedSecondaryVertexBJetTags();
    }
};


/// Normal Plots (TH1F* and TH2F*)
map<string,TH1F*> histo1D;
map<string,TH2F*> histo2D;

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

/// MultiPadPlot
map<string,MultiSamplePlot*> MultiPadPlot;

//TGraphAsymmErrors
map<string,TGraphAsymmErrors*> tgraph;

bool match;

int main (int argc, char *argv[])
{

  clock_t start = clock();

cout << "*************************************************************" << endl;
cout << " Beginning of the program for the FourTop search ! "           << endl;
cout << "*************************************************************" << endl;

    //SetStyle if needed
setTDRStyle();

string postfix = "_MCStudy"; // to relabel the names of the output file

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////// Configuration ///////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
string channelpostfix = "";
string xmlFileName = "";

bool Electron = false; // use Electron channel?
bool Muon = true; // use Muon channel?

if(Electron && Muon){
    cout << "  --> Using both Muon and Electron channel? Choose only one ( since different samples/skims are required)!" << endl;
	exit(1);
  }


  if(Muon){
	cout << " --> Using the Muon channel..." << endl;
	channelpostfix = "_Mu";
	xmlFileName = "config/test_mcstudy_mu.xml";
  }
  else if(Electron){
	cout << " --> Using the Electron channel..." << endl;
	channelpostfix = "_El";
	xmlFileName = "config/test_mcstudy_el.xml";
  }


const char *xmlfile = xmlFileName.c_str();
cout << "used config file: " << xmlfile << endl;    


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////// AnalysisEnvironment
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

AnalysisEnvironment anaEnv;
cout<<" - Loading environment ..."<<endl;
AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
int verbose = 2;//anaEnv.Verbose;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////// Load Datasets
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TTreeLoader treeLoader;
vector < Dataset* > datasets;
cout << " - Loading datasets ..." << endl;
treeLoader.LoadDatasets (datasets, xmlfile);
float Luminosity = 5343.64; //pb^-1??

bool Tprime = false; // If false, regular variables are used in MVA    

string MVAmethod = "Likelihood"; // MVAmethod to be used to get the good jet combi calculation (not for training! this is chosen in the jetcombiner class)

bool TrainMVA = true; // If false, the previously trained MVA will be used to calculate stuff
//JetCombiner* jetCombiner = new JetCombiner(TrainMVA, Luminosity, datasets, MVAmethod, Tprime);

    

int ndatasets = datasets.size() - 1 ;
cout << " - splitting TTBar dataset ..." << ndatasets   << endl;
vector<string> ttbar_filenames = datasets[ndatasets]->Filenames();
cout <<"ttbar filenames =  "<< ttbar_filenames[0] <<endl;

Dataset* ttbar_right = new Dataset("TTJets_right","right combinations" , true, 413, 2, 2, 1, 213.4,ttbar_filenames );
Dataset* ttbar_wrong = new Dataset("TTJets_wrong","wrong combinations" , true, 633, 2, 2, 1, 213.4, ttbar_filenames );    
    
    
MSPlot["TriJetMass"] = new MultiSamplePlot(datasets, "TriJetMass", 15, 0, 1000, "m_{bjj}");
//MSPlot["TriJetMass_wrong"] = new MultiSamplePlot(datasets, "TriJetMass_wrong", 15, 0, 1000, "m_{bjj}");

    
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

MSPlot["DeltaR_ClosestTop"] = new MultiSamplePlot(datasets, "DeltaR_ClosestTop ", 40, -10,10 , "Delta R Closest_top");
    
////////////////////////////////////////////////////////////////////
////////////////// 1D histograms  //////////////////////////////////
////////////////////////////////////////////////////////////////////

    
    
    
    tgraph["e_b_TTJets"] = new TGraphAsymmErrors();
    tgraph["e_b_TTTT"] = new TGraphAsymmErrors();

    tgraph["e_c_TTJets"] = new TGraphAsymmErrors();
    tgraph["e_c_TTTT"] = new TGraphAsymmErrors();
    
    // need light efficiency in the following eta bins: 0-0.8, 0.8-1.6, 1.6-2.4 
    
    tgraph["e_l_TTJets_eta0-08"] = new TGraphAsymmErrors();
    tgraph["e_l_TTTT_eta0-08"] = new TGraphAsymmErrors();
    
    tgraph["e_l_TTJets_eta08-16"] = new TGraphAsymmErrors();
    tgraph["e_l_TTTT_eta08-16"] = new TGraphAsymmErrors();
    
    tgraph["e_l_TTJets_eta16-24"] = new TGraphAsymmErrors();
    tgraph["e_l_TTTT_eta16-24"] = new TGraphAsymmErrors();
    
    //float ptmin[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600};
    //  float ptmax[] = {30, 40, 50, 60, 70, 80,100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
    int nbins = 16;
    float xbins[] = {20, 30, 40, 50, 60, 70, 80, 100, 120, 160, 210, 260, 320, 400, 500, 600, 800};
    
    
    
    
histo1D["hadronicPartonWMass"]= new TH1F("hadronicPartonWMass","hadronicPartonWMass;hadronicPartonWMass;#combinations",80,0,200);
histo1D["hadronicRecoWMass"]= new TH1F("hadronicRecoWMass","hadronicRecoWMass;hadronicRecoWMass;#combinations",80,30,310);    
histo1D["leptonicPartonWMass"]= new TH1F("leptonicPartonWMass","leptonicPartonWMass;#combinations",80,0,200);
histo1D["leptonicRecoWMass"]= new TH1F("leptonicRecoWMass","leptonicRecoWMass;#combinations",80,0,250);   

histo1D["hadronicPartonTopMass"]= new TH1F("hadronicPartonTopMass","hadronicPartonTopMass;hadronicPartonTopMass;#combinations",80,30,310);
histo1D["hadronicRecoTopMass"]= new TH1F("hadronicRecoTopMass","hadronicRecoTopMass;hadronicRecoTopMass;#combinations",80,30,310);        
histo1D["hadronicRecoTopMass_right"]= new TH1F("hadronicRecoTopMass_right","hadronicRecoTopMass_right;hadronicRecoTopMass_right;#combinations",80,30,550);        
histo1D["hadronicRecoTopMass_wrong"]= new TH1F("hadronicRecoTopMass_wrong","hadronicRecoTopMass_wrong;hadronicRecoTopMass_wrong;#combinations",80,30,550);        
histo1D["hadronicRecoTopMass_all"]= new TH1F("hadronicRecoTopMass_all","hadronicRecoTopMass_all;hadronicRecoTopMass_all;#combinations",80,30,550);        
    
histo1D["leptonicPartonTopMass"]= new TH1F("leptonicPartonTopMass","leptonicPartonTopMass;leptonicPartonTopMass;#combinations",80,100,250);
histo1D["leptonicRecoTopMass"]= new TH1F("leptonicRecoTopMass","leptonicRecoTopMass;leptonicRecoTopMass;#combinations",80,0,350); 

histo1D["lumiWeights"] = new TH1F("lumiWeights","lumiWeights;lumiWeight;#events",100,0,4);

histo1D["BJetPt_TTJets"]= new TH1F("BJetPt_TTJets","BJetPt;BJetPt;#combinations",nbins,xbins);
histo1D["BJetPt_tagged_TTJets"]= new TH1F("BJetPt_tagged_TTJets","BJetPt_tagged;BJetPt_tagged;#combinations",nbins,xbins);

histo1D["BJetPt_TTTT"]= new TH1F("BJetPt_TTTT","BJetPt;BJetPt;#combinations",nbins,xbins);
histo1D["BJetPt_tagged_TTTT"]= new TH1F("BJetPt_tagged_TTTT","BJetPt_tagged;BJetPt_tagged;#combinations",nbins,xbins);
    
    
histo1D["CJetPt_TTJets"]= new TH1F("CJetPt_TTJets","CJetPt;CJetPt;#combinations",nbins,xbins);
histo1D["CJetPt_tagged_TTJets"]= new TH1F("CJetPt_tagged_TTJets","CJetPt_tagged;CJetPt_tagged;#combinations",nbins,xbins);
    
histo1D["CJetPt_TTTT"]= new TH1F("CJetPt_TTTT","CJetPt;CJetPt;#combinations",nbins,xbins);
histo1D["CJetPt_tagged_TTTT"]= new TH1F("CJetPt_tagged_TTTT","CJetPt_tagged;CJetPt_tagged;#combinations",nbins,xbins);
    

histo1D["LJetPt_TTJets_eta0-08"]= new TH1F("LJetPt_TTJets_eta0-08","LJetPt;LJetPt;#combinations",nbins,xbins);
histo1D["LJetPt_tagged_TTJets_eta0-08"]= new TH1F("LJetPt_tagged_TTJets_eta0-08","LJetPt_tagged;LJetPt_tagged;#combinations",nbins,xbins);

histo1D["LJetPt_TTJets_eta08-16"]= new TH1F("LJetPt_TTJets_eta08-16","LJetPt;LJetPt;#combinations",nbins,xbins);
histo1D["LJetPt_tagged_TTJets_eta08-16"]= new TH1F("LJetPt_tagged_TTJets_eta08-16","LJetPt_tagged;LJetPt_tagged;#combinations",nbins,xbins);
    
histo1D["LJetPt_TTJets_eta16-24"]= new TH1F("LJetPt_TTJets_eta16-24","LJetPt;LJetPt;#combinations",nbins,xbins);
histo1D["LJetPt_tagged_TTJets_eta16-24"]= new TH1F("LJetPt_tagged_TTJets_eta16-24","LJetPt_tagged;LJetPt_tagged;#combinations",nbins,xbins);

    
histo1D["LJetPt_TTTT_eta0-08"]= new TH1F("LJetPt_TTTT_eta0-08","LJetPt;LJetPt;#combinations",nbins,xbins);
histo1D["LJetPt_tagged_TTTT_eta0-08"]= new TH1F("LJetPt_tagged_TTTT_eta0-08","LJetPt_tagged;LJetPt_tagged;#combinations",nbins,xbins);
    
histo1D["LJetPt_TTTT_eta08-16"]= new TH1F("LJetPt_TTTT_eta08-16","LJetPt;LJetPt;#combinations",nbins,xbins);
histo1D["LJetPt_tagged_TTTT_eta08-16"]= new TH1F("LJetPt_tagged_TTTT_eta08-16","LJetPt_tagged;LJetPt_tagged;#combinations",nbins,xbins);
    
histo1D["LJetPt_TTTT_eta16-24"]= new TH1F("LJetPt_TTTT_eta16-24","LJetPt;LJetPt;#combinations",nbins,xbins);
histo1D["LJetPt_tagged_TTTT_eta16-24"]= new TH1F("LJetPt_tagged_TTTT_eta16-24","LJetPt_tagged;LJetPt_tagged;#combinations",nbins,xbins);
    
    
for (unsigned int d = 0; d < datasets.size(); d++){
    histo2D[("RelIso_vs_MET_"+datasets[d]->Name()).c_str()] = new TH2F(("RelIso_vs_MET_"+datasets[d]->Name()).c_str(),"RelIso:MET",100,0,1000, 100, 0,1);
    histo2D[("ThirdTopMass_vs_HT_"+datasets[d]->Name()).c_str()] = new TH2F(("ThirdTopMass_vs_HT_"+datasets[d]->Name()).c_str(),"ThirdTopMass_vs_HT",30,0,1000, 30, 0,1500);
    histo2D[("ThirdTopMass_vs_SecondTopMass"+datasets[d]->Name()).c_str()] = new TH2F(("ThirdTopMass_vs_SecondTopMass"+datasets[d]->Name()).c_str(),"ThirdTopMass_vs_SecondTopMass",30,0,1000, 30, 0,1000);
      
    histo2D[("MassChi2_vs_HT"+datasets[d]->Name()).c_str()] = new TH2F(("MassChi2_vs_HT"+datasets[d]->Name()).c_str(),"MassChi2_vs_HT",15,0,400, 15, 0,1400);
    histo2D[("EventMassX_vs_HT"+datasets[d]->Name()).c_str()] = new TH2F(("EventMassX_vs_HT"+datasets[d]->Name()).c_str(),"EventMassX_vs_HT",20,0,1500, 15, 0,1200);
    histo2D[("EventMassX_vs_HTX"+datasets[d]->Name()).c_str()] = new TH2F(("EventMassX_vs_HTX"+datasets[d]->Name()).c_str(),"EventMassX_vs_HTX",20,0,1500, 15, 0,1200);

  }
 
//  cout << " - Declared histograms ..." <<  endl;
	
////////////////////////////////////////////////////////////////////
////////////////// Plots  //////////////////////////////////////////
////////////////////////////////////////////////////////////////////

string pathPNG = "FourTop"+postfix+channelpostfix;
pathPNG += "_MSPlots/"; 	
mkdir(pathPNG.c_str(),0777);

////////////////////////////////////////////////////////////////////
///////////////////// Tree for parton level studies    ////////////////////////
////////////////////////////////////////////////////////////////////    
    Float_t SumPtg, RPerg, RPerllg, RPerblg, VecSumPtg,BDiscRatg, WMassg , TopMassg ;
    Float_t SumPtb, RPerb, RPerllb, RPerblb, VecSumPtb,BDiscRatb, WMassb, TopMassb;
    TFile* treeFile;
    treeFile = new TFile(TreeFileName.c_str(),"RECREATE");
   // TTree* myFourTopTree;
    TTree* combinationsTree_good;
    TTree* combinationsTree_bad;
    combinationsTree_good = new TTree("CombinationsTree_good","Tree containing the good combinations");
    combinationsTree_bad = new TTree("CombinationsTree_bad","Tree containing the bad combinations");

    TBranch* good =0 ;

    combinationsTree_good->Branch("SumPt",&SumPtg, "SumPtg/F");   
    combinationsTree_good->Branch("RPer",&RPerg, "RPerg/F"); 
    combinationsTree_good->Branch("RPerll",&RPerllg, "RPerllg/F");   
    combinationsTree_good->Branch("RPerbl",&RPerblg, "RPerblg/F");   

    combinationsTree_good->Branch("VecSumPt",&VecSumPtg, "VecSumPtg/F");   
    combinationsTree_good->Branch("BDiscRat",&BDiscRatg, "BDiscRatg/F");   
    combinationsTree_good->Branch("WMass",&WMassg, "WMassg/F"); 
    combinationsTree_good->Branch("TopMass",&TopMassg, "TopMassg/F");   

    combinationsTree_bad->Branch("SumPt",&SumPtb, "SumPtb/F");   
    combinationsTree_bad->Branch("RPer",&RPerb, "RPerb/F");   
    combinationsTree_bad->Branch("RPerll",&RPerllb, "RPerllb/F");   
    combinationsTree_bad->Branch("RPerbl",&RPerblb, "RPerblb/F"); 
    combinationsTree_bad->Branch("VecSumPt",&VecSumPtb, "VecSumPtb/F");   
    combinationsTree_bad->Branch("BDiscRat",&BDiscRatb, "BDiscRatb/F");   
    combinationsTree_bad->Branch("WMass",&WMassb, "WMassb/F"); 
    combinationsTree_bad->Branch("TopMass",&TopMassb, "TopMassb/F");   


    
   // FourTopTree* myBranch_selectedEvents = 0;		
    bool datadriven = false;
    //if(!datadriven){
     //   treeFile = new TFile(TreeFileName.c_str(),"RECREATE");
    //    myFourTopTree = new TTree("FourTopTree","Tree containing the four top information");
      //  myFourTopTree->Branch("FourTopBranch_selectedEvents","FourTopTree",&myBranch_selectedEvents);
    //}
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////// Loop on dataset ////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
      
      
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ///////////////////////////////////////////////// Initialize JEC factors /////////////////////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      
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
      
      
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////                      Loop on events                                                    ///////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    double scaleFactor = 1.;
    int itrigger = -1, previousRun = -1;
    int start = 0;
      
    unsigned int end = datasets[d]->NofEvtsToRunOver();
    cout <<"Number of events = "<<  end  <<endl;

    bool debug = false;

    if (verbose > 1) cout << " - Loop over events " << endl;
      
      double end_d =end;
      
      if (dataSetName=="TTJets"){
     //    end_d = end;
          end_d = end_d/1.;
      }
      
    for (unsigned int ievt = start; ievt < 100000; ievt++)
    {
	if(ievt%1000 == 0)
		std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";

	//load event
	event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);

	vector<TRootGenJet*> genjets;
	if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
	{
	  // loading GenJets as I need them for JER
	  		genjets = treeLoader.LoadGenJet(ievt);
	}

     if (debug)   cout << " got gen jets " << endl;

        
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////// Jet energy scale corrections     ///////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        // Apply Jet Corrections on-the-fly
        //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction on the fly:");
        //if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
           // jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),true); //last boolean: isData (needed for L2L3Residual...)
       // else
            //jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),false); //last boolean: isData (needed for L2L3Residual...)
        //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JES correction on the fly:");  
        
        
// Declare selection instance    
Selection selection(init_jets, init_muons, init_electrons, mets);     
        
// Define object selection cuts
selection.setJetCuts(40.,2.5,0.01,1.,0.98,0,0);//Pt, Eta, EMF, n90Hits, fHPD, dRJetElectron, DRJets
//selection.setElectronCuts(30.,2.5,0.1,0.02,0.,999.,0,1); //Pt,Eta,RelIso,d0,MVAId,DistVzPVz, DRJets
        selection.setLooseElectronCuts(20,2.5,0.2,0.);

        selection.setLooseElectronCuts(20,2.5,0.2,0.);
//selection.setMuonCuts(26.,2.1,.12,0,0.2,0,1,0.5,5,1 );
        
selection.setMuonCuts(26,2.1,0.12,0.2,0.3,1,0.5,5,0);
 //Pt,Eta,RelIso,NValidMuHits,d0, dRJets, NMatchedStations,DistVzPVz, 
                //float, float, float, int, float, float, int, float, int, int)
        
selection.setLooseMuonCuts(10,2.5,0.2);

//Select objects 
//	vector<TRootElectron*> selectedElectrons_NoIso = selection.GetSelectedElectrons(20,2.4,999.,vertex[0]);
vector<TRootElectron*> selectedElectrons       = selection.GetSelectedElectrons();
vector<TRootMuon*>     selectedMuons       = selection.GetSelectedMuons();
vector<TRootJet*>      selectedJets        = selection.GetSelectedJets(true); // ApplyJetId
vector<TRootMuon*>     selectedLooseMuons     = selection.GetSelectedLooseMuons();
vector<TRootElectron*> selectedLooseElectrons = selection.GetSelectedElectrons(); // VBTF ID
vector<TRootJet*>      selectedBJets; // B-Jets
vector<TRootJet*>      selectedLightJets; // light-Jets      
//vector<TRootJet*>      selectedJetsTLV;
//order jets wrt to Pt, then set bool corresponding to RefSel cuts.             
sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order muons wrt Pt. 
        
       // cout<< "sel jets size :" << selectedJets.size()  <<endl;

        
        selectedBJets.clear();
        
        for(unsigned int i=0; i<selectedJets.size(); i++){
            selectedBJets.push_back(selectedJets[i]);   
        }

        
	if(Muon){
	  if (debug) cout<< "Event passed..."<<endl;
   
        //applying baseline
        
	  //   if (!(selectedMuons.size() == 1   && selectedJets.size() ==6  && selectedBJets.size() >=2   )) continue;
        
     if ( !(selectedElectrons.size() == 1  && selectedLooseMuons.size() ==0  && selectedLooseElectrons.size() ==1  && mets[0]->Et()> 30.   ) ) continue;
        
////////////////////////////////////////////////////////////////////////////////////
//// Getting Gen Event
////////////////////////////////////////////////////////////////////////////////////
if(dataSetName != "data" && dataSetName != "Data" && dataSetName != "Data"){
             
bool Wbosonpartonsmatched = false; // True if the Wboson ttbar semi-mu or semi-el partons are matched to their 2 jets (not necessarily the 4 highest pt jets)
float WMassmatched_had = -9999;
float TopMassmatched_had = -9999;
float WMassmatched_lep = -9999;
float TopMassmatched_lep = -9999;
    
    
    
pair<unsigned int, unsigned int> leptonicBJet_, hadronicBJet_, hadronicWJet1_, hadronicWJet2_; //First index is the JET number, second one is the parton
leptonicBJet_ = hadronicBJet_ = hadronicWJet1_ = hadronicWJet2_ = pair<unsigned int,unsigned int>(9999,9999);
vector<TRootMCParticle*> mcParticles;
vector<TRootMCParticle*> mcTops;
vector<TRootMCParticle*>  mcBs;
vector<TRootMCParticle*> mcLights;
vector<TRootMCParticle*> mcLeps;
vector<TRootMCParticle*> mcNus;
vector<TRootMCParticle*> mcQGs;
vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV;
vector<TRootMCParticle*> mcParticlesMatching_; 
vector<TRootJet*> selectedJets_;
//vector<TRootJet*> selectedJets;
    
    
mcParticlesMatching_.clear();
selectedJets_.clear();
mcParticlesTLV.clear(); 
selectedJetsTLV.clear();


vector<TLorentzVector*> mcTops_vec;
vector<TLorentzVector*>  mcBs_vec;
vector<TLorentzVector*> mcLights_vec;
vector<TLorentzVector*> mcLeps_vec;
vector<TLorentzVector*> mcNus_vec;
             
int MCPermutation[4]; for (unsigned int i=0;i<4;i++) MCPermutation[i] = -1;


//load Gen event
TRootGenEvent* genEvt = 0;
genEvt = treeLoader.LoadGenEvent(ievt,false);
treeLoader.LoadMCEvent(ievt, genEvt, 0, mcParticles,true);  
        
if (debug) cout <<"size   "<< mcParticles.size()<<endl;
        
//myBranch_selectedEvents = new FourTopTree();
//    cout <<" "<< endl;
            
TLorentzVector met(1.,1.,1.,1.);
TLorentzVector* mcTop1_vec = new TLorentzVector();
TLorentzVector* mcTop2_vec = new TLorentzVector();
TLorentzVector* mcTop3_vec = new TLorentzVector();
TLorentzVector* mcTop4_vec = new TLorentzVector();
             
TLorentzVector* mcB1_vec = new TLorentzVector();
TLorentzVector* mcB2_vec = new TLorentzVector();
TLorentzVector* mcB3_vec = new TLorentzVector();
TLorentzVector* mcB4_vec = new TLorentzVector();     
             
TLorentzVector* mcLight1_vec = new TLorentzVector();
TLorentzVector* mcLight2_vec = new TLorentzVector();
TLorentzVector* mcLight3_vec = new TLorentzVector();
TLorentzVector* mcLight4_vec = new TLorentzVector();    
TLorentzVector* mcLight5_vec = new TLorentzVector();
TLorentzVector* mcLight6_vec = new TLorentzVector();         

TLorentzVector* mcLep_vec = new TLorentzVector();
TLorentzVector* mcNu_vec = new TLorentzVector();  
    
    bool muPlusFromTop = false, muMinusFromTop = false;
    bool elPlusFromTop = false, elMinusFromTop = false;
    int leptonPDG, muonPDG = 13, electronPDG = 11;
    leptonPDG = muonPDG; 
             
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
    
        Int_t dau1Id = mcParticles[i]->dauOneId();            
        Int_t dau2Id = mcParticles[i]->dauTwoId();            
        Int_t dau3Id = mcParticles[i]->dauThreeId();          
                              
        if(abs(mcParticles[i]->type()) == 13 && abs(mcParticles[i]->motherType())==24) mcLeps.push_back(mcParticles[i]);
        if(abs(mcParticles[i]->type()) == 14 && abs(mcParticles[i]->motherType())==24) mcNus.push_back(mcParticles[i]);
        if(abs(mcParticles[i]->type())< 5 && abs(mcParticles[i]->motherType())==24 ) mcLights.push_back(mcParticles[i]);
        if(abs(mcParticles[i]->type())==5 && abs(mcParticles[i]->motherType())==6  ) mcBs.push_back(mcParticles[i]);
        if(abs(mcParticles[i]->type())==6 ) mcTops.push_back(mcParticles[i]);
        }
             
             if (debug) cout<<" filled mc vectors...  "<<endl;
     //   if (mcLeps.size() != 1) continue;
      
             if (debug) cout<<" jet matching...  "<<endl;             
////////////////////
/// Jet Matching ///    
////////////////////
     //cout<<"  jet-parts size...  "<< selectedJets.size() <<endl;
    
    for(unsigned int i=0; i<selectedJets.size(); i++) selectedJetsTLV.push_back(*selectedJets[i]);
    if (debug) cout<<" jet matching... 1 "<<endl;

    
    JetPartonMatching matching = JetPartonMatching(mcParticlesTLV, selectedJetsTLV, 2, true, true, 0.3);
            
        vector< pair<unsigned int, unsigned int> > JetPartonPair; // First one is jet number, second one is mcParticle number
            

             for(unsigned int i=0; i<mcParticlesTLV.size(); i++) //loop through mc particles and find matched jets
             {
                 int matchedJetNumber = matching.getMatchForParton(i, 0);
                 if(matchedJetNumber != -1)
                 JetPartonPair.push_back( pair<unsigned int, unsigned int> (matchedJetNumber, i) );
             }
    //cout <<"  "<<endl;
    //cout<<" N jets.....  "<< selectedJets.size()  <<endl;             
    //cout<<" N partons.....  "<< mcParticlesTLV.size() <<endl;             
    //cout<<" N jetpartonpairs.....  "<< JetPartonPair.size() <<endl; 
    
    
             
             for(unsigned int i=0; i<JetPartonPair.size(); i++)//looping through matched jet-parton pairs
             {
                 unsigned int j = JetPartonPair[i].second;	  //get index of matched mc particle
                 if( fabs(mcParticlesMatching_[j]->type()) < 5 )
                 {
                     //cout << "found true light-jet " << JetPartonPair[i].first <<endl;
                     if( ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == -24 && mcParticlesMatching_[j]->grannyType() == -6 )
		  				|| ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == 24 && mcParticlesMatching_[j]->grannyType() == 6 ) )
                     {
                         if(hadronicWJet1_.first == 9999)
                         {
                             hadronicWJet1_ = JetPartonPair[i];
                             MCPermutation[0] = JetPartonPair[i].first;
                         }
                         else if(hadronicWJet2_.first == 9999)
                         {
                             hadronicWJet2_ = JetPartonPair[i];
                             MCPermutation[1] = JetPartonPair[i].first;
                         } 
                         else cerr<<"Found a third jet coming from a W boson which comes from a top quark..."<<endl;
                     }
                     
                     if (debug) cout<<" jet matching...2  "<<endl;

                    // cout << "found true light-jet 2"  <<endl;

                 }
                 
                 if( fabs(mcParticlesMatching_[j]->type()) == 5 )
                 {
//                cout << "found true b-jet " << JetPartonPair[i].first <<endl;

                     n_true_bjets++;
                     if(dataSetName=="TTJets"){ histo1D["BJetPt_TTJets"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                        // cout <<"filling plots...ttjets"<<endl;
                     }
                     
                     if (debug) cout<<" jet matching...3  "<<endl;

                     
                     if(dataSetName=="TTTT"){ histo1D["BJetPt_TTTT"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());

//                        cout << "found true b-jet " << JetPartonPair[i].first <<endl;
   
                     }
                     if (selectedJets[JetPartonPair[i].first]->btag_combinedSecondaryVertexBJetTags() > workingpointvalue){
                         n_tagged_true_bjets++;
                         if(dataSetName=="TTJets") {  histo1D["BJetPt_tagged_TTJets"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                           //  cout <<"filling plots...ttjets 2"<<endl;
   
                         }
                         if(dataSetName=="TTTT")   histo1D["BJetPt_tagged_TTTT"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                        
                     }
                     
                     
                     if (debug) cout<<" jet matching...4  "<<endl;

                     if(  ( muPlusFromTop && mcParticlesMatching_[j]->motherType() == -6) || ( muMinusFromTop && mcParticlesMatching_[j]->motherType() == 6 ) )
                     {
                         hadronicBJet_ = JetPartonPair[i];
                         MCPermutation[2] = JetPartonPair[i].first;
                     }
                     else if((muPlusFromTop && mcParticlesMatching_[j]->motherType() == 6) ||  ( muMinusFromTop &&mcParticlesMatching_[j]->motherType() == -6) ) 
                     {
                         leptonicBJet_ = JetPartonPair[i];
                         MCPermutation[3] = JetPartonPair[i].first;
                     }
                 }  else if ( fabs(mcParticlesMatching_[j]->type()) == 4 )     {
                 
                     if (debug) cout<<" jet matching...c1  "<<endl;

                     n_true_cjets++;
                     if(dataSetName=="TTJets"){ histo1D["CJetPt_TTJets"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                         // cout <<"filling plots...ttjets"<<endl;
                     }
                     
                     if(dataSetName=="TTTT"){ histo1D["CJetPt_TTTT"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                         
                         //                        cout << "found true b-jet " << JetPartonPair[i].first <<endl;
                         
                     }
                     if (selectedJets[JetPartonPair[i].first]->btag_combinedSecondaryVertexBJetTags() > workingpointvalue){
                         n_tagged_true_bjets++;
                         if(dataSetName=="TTJets") {  histo1D["CJetPt_tagged_TTJets"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                             //  cout <<"filling plots...ttjets 2"<<endl;
                             
                         }
                         if(dataSetName=="TTTT")   histo1D["CJetPt_tagged_TTTT"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                         
                     }
                 
                 }
                 
                 else{
                     if (debug) cout<<" jet matching...L1  "<<endl;

                 
                     n_true_ljets++;
                     if(dataSetName=="TTJets"){
                        // cout <<"filled plots...ttjets 1"<<endl;
                         
    if (fabs(selectedJets[JetPartonPair[i].first]->Eta()) < 0.8) histo1D["LJetPt_TTJets_eta0-08"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                         //cout <<"filled plots...ttjets 2"<<endl;

    if (fabs(selectedJets[JetPartonPair[i].first]->Eta()) > 0.8  && fabs(selectedJets[JetPartonPair[i].first]->Eta()) < 1.6 )  histo1D["LJetPt_TTJets_eta08-16"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                         //cout <<"filled plots...ttjets 3"<<endl;

    if (fabs(selectedJets[JetPartonPair[i].first]->Eta()) > 1.6  && fabs(selectedJets[JetPartonPair[i].first]->Eta()) < 2.4 )histo1D["LJetPt_TTJets_eta16-24"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                        // cout <<"filled plots...ttjets 4"<<endl;

                
                     }
                     
                     if(dataSetName=="TTTT"){
                        // cout <<"filling plots...tttt"<<endl;
                
                         
                         if (fabs(selectedJets[JetPartonPair[i].first]->Eta()) < 0.8) histo1D["LJetPt_TTTT_eta0-08"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                         if (fabs(selectedJets[JetPartonPair[i].first]->Eta()) > 0.8  && fabs(selectedJets[JetPartonPair[i].first]->Eta()) < 1.6 )  histo1D["LJetPt_TTTT_eta08-16"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                         if (fabs(selectedJets[JetPartonPair[i].first]->Eta()) > 1.6  && fabs(selectedJets[JetPartonPair[i].first]->Eta()) < 2.4 )histo1D["LJetPt_TTTT_eta16-24"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                         
                         
                         //                        cout << "found true b-jet " << JetPartonPair[i].first <<endl;
                         
                     }
                     if (selectedJets[JetPartonPair[i].first]->btag_combinedSecondaryVertexBJetTags() > workingpointvalue){
                         n_tagged_true_bjets++;
                         if(dataSetName=="TTJets") {
                           //  cout <<"filling tagged plots...ttjets 1"<<endl;

                             
                        if (fabs(selectedJets[JetPartonPair[i].first]->Eta()) < 0.8) histo1D["LJetPt_tagged_TTJets_eta0-08"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                             //cout <<"filling tagged plots...ttjets 2"<<endl;
                             
                        if (fabs(selectedJets[JetPartonPair[i].first]->Eta()) > 0.8  && fabs(selectedJets[JetPartonPair[i].first]->Eta()) < 1.6 )
                            histo1D["LJetPt_tagged_TTJets_eta08-16"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                             
                             
                            // cout <<"filling tagged plots...ttjets 3"<<endl;

                             
                        if (fabs(selectedJets[JetPartonPair[i].first]->Eta()) > 1.6  && fabs(selectedJets[JetPartonPair[i].first]->Eta()) < 2.4 ) histo1D["LJetPt_tagged_TTJets_eta16-24"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                             
                             //cout <<"filling tagged plots...ttjets 4"<<endl;

                             
                             
                             //  cout <<"filling plots...ttjets 2"<<endl;
                             
                         }
                         if(dataSetName=="TTTT"){
                             
                             
                             if (fabs(selectedJets[JetPartonPair[i].first]->Eta()) < 0.8) histo1D["LJetPt_tagged_TTTT_eta0-08"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                             if (fabs(selectedJets[JetPartonPair[i].first]->Eta()) > 0.8  && fabs(selectedJets[JetPartonPair[i].first]->Eta()) < 1.6 )  histo1D["LJetPt_tagged_TTTT_eta08-16"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                             if (fabs(selectedJets[JetPartonPair[i].first]->Eta()) > 1.6  && fabs(selectedJets[JetPartonPair[i].first]->Eta()) < 2.4 )histo1D["LJetPt_tagged_TTTT_eta16-24"]->Fill(selectedJets[JetPartonPair[i].first]->Pt());
                             
                             
                         }
                     }
                 
                 
                 
                 }
                 
                 
             }
    
   // cout <<"filling plots...ttjets 3"<<endl;

    
             if (debug) cout<<" filling histos.....  "<<endl;
             

             if(hadronicWJet1_.first != 9999 && hadronicWJet2_.first != 9999 && hadronicWJet1_.first != hadronicWJet2_.first && hadronicBJet_.first != 9999 && leptonicBJet_.first != 9999  && hadronicBJet_.first != leptonicBJet_.first  && mcLeps.size() !=0 )
             {
                 
                 if (debug) cout<<" matching.....  "<<endl;

                 
                 histo1D["hadronicPartonWMass"]->Fill((*mcParticlesMatching_[hadronicWJet1_.second]+*mcParticlesMatching_[hadronicWJet2_.second]).M(),scaleFactor);
                 histo1D["hadronicPartonTopMass"]->Fill((*mcParticlesMatching_[hadronicWJet1_.second]+*mcParticlesMatching_[hadronicWJet2_.second]+*mcParticlesMatching_[hadronicBJet_.second]).M(),scaleFactor);

            //     histo1D["leptonicPartonWMass"]->Fill((*mcNus[0]+*mcLeps[0]).M(),scaleFactor);
             //    histo1D["leptonicPartonTopMass"]->Fill((*mcNus[0]+*mcLeps[0] +  *mcParticlesMatching_[leptonicBJet_.second]).M() ,scaleFactor);
                 if (debug) cout<<" matched.....  "<<endl;

                 Wbosonpartonsmatched = true;
             }					
             
             if(Wbosonpartonsmatched)
             {
if (debug)     cout <<"partons matched to W"<<endl;
                 
    TopMassmatched_had = (*selectedJets[hadronicWJet1_.first]+*selectedJets[hadronicWJet2_.first]+*selectedJets[hadronicBJet_.first]).M();
    WMassmatched_had = (*selectedJets[hadronicWJet1_.first]+*selectedJets[hadronicWJet2_.first]).M();
                
      if (debug)   cout <<"filling 1 "<< selectedMuons.size()  <<endl;
                 if(selectedMuons.size() !=0){
        TopMassmatched_lep = (*selectedMuons[0]+*selectedJets[leptonicBJet_.first]+*mets[0]).M();    
        WMassmatched_lep = (*selectedMuons[0]+*mets[0]).M();    
                 }
        histo1D["leptonicRecoWMass"]->Fill(WMassmatched_lep,scaleFactor);
        histo1D["leptonicRecoTopMass"]->Fill(TopMassmatched_lep,scaleFactor);
                 
        histo1D["hadronicRecoWMass"]->Fill(WMassmatched_had,scaleFactor);
        histo1D["hadronicRecoTopMass"]->Fill(TopMassmatched_had,scaleFactor);
                 
      //  histo1D["hadronicRecoTopMass_right"]->Fill(TopMassmatched_had,scaleFactor);

                if (debug) cout <<"filling 2 "<< endl;
                 
             }
             else{
                 histo1D["hadronicRecoTopMass_wrong"]->Fill(TopMassmatched_had,scaleFactor);
                             
             }
  if (debug) cout <<"filling plots...ttjets neste"<<endl;

    
        for(unsigned int i=0; i<selectedJets.size(); i++){
            for(unsigned int j=0; j<selectedJets.size(); j++){
                for(unsigned int k=0; k<selectedJets.size(); k++){
                    if (k==j || k==i || i==j) continue;
                    histo1D["hadronicRecoTopMass_all"]->Fill((*selectedJets[i]+*selectedJets[j]+*selectedJets[k]).M()  ,scaleFactor);
            
               if (debug)     cout <<"filling plots...ttjets 1A"<<endl;
                    
                    if ( !((i == hadronicWJet1_.first || i == hadronicWJet2_.first ) && (j == hadronicWJet1_.first || j == hadronicWJet2_.first ) && ( k == hadronicBJet_.first ) ))
                    {
                        
                 if (debug)       cout <<"filling plots...ttjets 1B"<<endl;
                        //JetPermutation.clear();
                        //JetPermutation.push_back(selectedJets[i]);
                        //JetPermutation.push_back(selectedJets[j]);
                        //JetPermutation.push_back(selectedJets[k]);
                        //sort(JetPermutation.begin(), JetPermutation.end(),HighestCSVBtag());
                        BDiscRatb = (selectedJets[k]->btag_combinedSecondaryVertexBJetTags()/ (selectedJets[i]->btag_combinedSecondaryVertexBJetTags() + selectedJets[j]->btag_combinedSecondaryVertexBJetTags()));
                   if (debug)    cout <<"filling plots...ttjets 11C"<<endl;

                        WMassb = ( *selectedJets[i] + *selectedJets[j]).M();
                     if (debug)   cout <<"filling plots...ttjets 11D"<<endl;

                        TopMassb = ( *selectedJets[i] + *selectedJets[j] + *selectedJets[k]).M();
                       if (debug) cout <<"filling plots...ttjets 11E"<<endl;

                        SumPtb = selectedJets[i]->Pt() + selectedJets[j]->Pt() + selectedJets[k]->Pt();
                     if (debug)   cout <<"filling plots...ttjets 11F"<<endl;

            
                        RPerb =  sqrt(   pow(selectedJets[i]->Eta() - selectedJets[j]->Eta(),2) + pow(selectedJets[i]->Phi() - selectedJets[j]->Phi(),2 )   )  +  sqrt(   pow(selectedJets[j]->Eta() - selectedJets[k]->Eta(),2) + pow(selectedJets[j]->Phi() - selectedJets[k]->Phi(),2 )   )  + sqrt(   pow(selectedJets[k]->Eta() - selectedJets[i]->Eta(),2) + pow(selectedJets[k]->Phi() - selectedJets[i]->Phi(),2 )) ;
                        
                if (debug) cout <<"filling plots...ttjets 11G"<<endl;

                        RPerllb =  sqrt(   pow(selectedJets[i]->Eta() - selectedJets[j]->Eta(),2) + pow(selectedJets[i]->Phi() - selectedJets[j]->Phi(),2 )   ) ;
                    if (debug)    cout <<"filling plots...ttjets 11H"<<endl;

                        RPerblb = sqrt(   pow(selectedJets[k]->Eta() - selectedJets[i]->Eta(),2) + pow(selectedJets[k]->Phi() - selectedJets[i]->Phi(),2 ))  ;
                    if (debug)    cout <<"filling plots...ttjets 11I"<<endl;
                        
                        
    if (SumPtb !=0.)  VecSumPtb = ((*selectedJets[i] + *selectedJets[j] + *selectedJets[k]).Pt() / SumPtb   );
                //        combinationsTree_bad->Fill();
                histo1D["hadronicRecoTopMass_wrong"]->Fill((*selectedJets[i]+*selectedJets[j]+*selectedJets[k]).M()  ,scaleFactor);
                       if (debug) cout <<"filling plots...ttjets 11J"<<endl;

                        MSPlot["TriJetMass"]->Fill((*selectedJets[i]+*selectedJets[j]+*selectedJets[k]).M(),  ttbar_wrong , true, Luminosity*scaleFactor );
                    if (debug)   cout <<"filling plots...ttjets 11K"<<endl;
                    }
                    else{
                         histo1D["hadronicRecoTopMass_right"]->Fill((*selectedJets[i]+*selectedJets[j]+*selectedJets[k]).M()  ,scaleFactor);
                        //JetPermutation.clear();
                        //JetPermutation.push_back(selectedJets[i]);
                        //JetPermutation.push_back(selectedJets[j]);
                        //JetPermutation.push_back(selectedJets[k]);
                        //sort(JetPermutation.begin(), JetPermutation.end(),HighestCSVBtag());
                    
                        BDiscRatg = (selectedJets[k]->btag_combinedSecondaryVertexBJetTags()/ (selectedJets[i]->btag_combinedSecondaryVertexBJetTags() + selectedJets[j]->btag_combinedSecondaryVertexBJetTags()));
                        WMassg = ( *selectedJets[i] + *selectedJets[j]).M();
                        TopMassg = ( *selectedJets[i] + *selectedJets[j] + *selectedJets[k]).M();
                        SumPtg = selectedJets[i]->Pt() + selectedJets[j]->Pt() + selectedJets[k]->Pt();
                        VecSumPtg = ((*selectedJets[i] + *selectedJets[j] + *selectedJets[k]).Pt() / SumPtg   );
                        
                        RPerg =  sqrt(   pow(selectedJets[i]->Eta() - selectedJets[j]->Eta(),2) + pow(selectedJets[i]->Phi() - selectedJets[j]->Phi(),2 )   )  +  sqrt(   pow(selectedJets[j]->Eta() - selectedJets[k]->Eta(),2) + pow(selectedJets[j]->Phi() - selectedJets[k]->Phi(),2 )   )  + sqrt(   pow(selectedJets[k]->Eta() - selectedJets[i]->Eta(),2) + pow(selectedJets[k]->Phi() - selectedJets[i]->Phi(),2 )) ;
                        
                        RPerllg =  sqrt(   pow(selectedJets[i]->Eta() - selectedJets[j]->Eta(),2) + pow(selectedJets[i]->Phi() - selectedJets[j]->Phi(),2 )   ) ;
                        
                        RPerblg = sqrt(   pow(selectedJets[k]->Eta() - selectedJets[i]->Eta(),2) + pow(selectedJets[k]->Phi() - selectedJets[i]->Phi(),2 ))  ;
                        
                       if (debug) cout <<"filling plots...ttjets 111"<<endl;

                        MSPlot["TriJetMass"]->Fill((*selectedJets[i]+*selectedJets[j]+*selectedJets[k]).M(),  ttbar_right , true, Luminosity*scaleFactor );
                       
                    }
                }

                
            }

    }
    if (debug) cout <<"filling plots...ttjets rper"<<endl;
    if (RPerg > 10000.) RPerg =-1.;
   // combinationsTree_good->Fill();
             //MC Tree stuff
             //cout <<mcTops.size()  << "  "<<mcBs.size()  << "  "<<mcLights.size() << "  "<<mcLeps.size()<<"  " << mcNus.size() << endl;    
             
            if(mcTops.size() > 0)mcTop1_vec->SetPxPyPzE(mcTops[0]->Px(),mcTops[0]->Py(),mcTops[0]->Pz(),mcTops[0]->E());
            if(mcBs.size() > 0) mcB1_vec->SetPxPyPzE(mcBs[0]->Px(),mcBs[0]->Py(),mcBs[0]->Pz(),mcBs[0]->E());
             
            if(mcTops.size() > 1) mcTop2_vec->SetPxPyPzE(mcTops[1]->Px(),mcTops[1]->Py(),mcTops[1]->Pz(),mcTops[1]->E());
            if(mcBs.size() > 1) mcB2_vec->SetPxPyPzE(mcBs[1]->Px(),mcBs[1]->Py(),mcBs[1]->Pz(),mcBs[1]->E());
    
             
            if(mcTops.size() > 2)  mcTop3_vec->SetPxPyPzE(mcTops[2]->Px(),mcTops[2]->Py(),mcTops[2]->Pz(),mcTops[2]->E());   
            if(mcBs.size() > 2)    mcB3_vec->SetPxPyPzE(mcBs[2]->Px(),mcBs[2]->Py(),mcBs[2]->Pz(),mcBs[2]->E());   
             
            if(mcTops.size() > 3) mcTop4_vec->SetPxPyPzE(mcTops[3]->Px(),mcTops[3]->Py(),mcTops[3]->Pz(),mcTops[3]->E());   
            if(mcBs.size() > 3)  mcB4_vec->SetPxPyPzE(mcBs[3]->Px(),mcBs[3]->Py(),mcBs[3]->Pz(),mcBs[3]->E());   
           
             
            if (mcLights.size() >0)mcLight1_vec->SetPxPyPzE(mcLights[0]->Px(),mcLights[0]->Py(),mcLights[0]->Pz(),mcLights[0]->E());
            if (mcLights.size() >1)mcLight2_vec->SetPxPyPzE(mcLights[1]->Px(),mcLights[1]->Py(),mcLights[1]->Pz(),mcLights[1]->E());
             
            if (mcLights.size() >2) mcLight3_vec->SetPxPyPzE(mcLights[2]->Px(),mcLights[2]->Py(),mcLights[2]->Pz(),mcLights[2]->E());
            if (mcLights.size() >3) mcLight4_vec->SetPxPyPzE(mcLights[3]->Px(),mcLights[3]->Py(),mcLights[3]->Pz(),mcLights[3]->E());
             
            if (mcLights.size() >4)mcLight5_vec->SetPxPyPzE(mcLights[4]->Px(),mcLights[4]->Py(),mcLights[4]->Pz(),mcLights[4]->E());
            if (mcLights.size() >5)mcLight6_vec->SetPxPyPzE(mcLights[5]->Px(),mcLights[5]->Py(),mcLights[5]->Pz(),mcLights[5]->E());
             
            if (mcLeps.size() >0)mcLep_vec->SetPxPyPzE(mcLeps[0]->Px(),mcLeps[0]->Py(),mcLeps[0]->Pz(),mcLeps[0]->E());   
            if (mcNus.size() >0)mcNu_vec->SetPxPyPzE(mcNus[0]->Px(),mcNus[0]->Py(),mcNus[0]->Pz(),mcNus[0]->E());   

                
         }
       if (debug) cout <<"filling plots...ttjets out"<<endl;

        
        
        if (debug) cout <<"filling tree  "<<endl;        
     //  if (datasets.size() == 1)  myFourTopTree->Fill(); 
        
        if (debug) cout <<"filling tree 2 "<<endl;
    //   delete myBranch_selectedEvents;
        
	}
	else if(Electron){
	  //	MSPlot["1stLeadingElectronRelIsolation"]->Fill(selectedElectrons_NoIso[0]->relativePfIso(), datasets[d], true, Luminosity*scaleFactor);

	}
        
    }//loop on events
    
      
     if (debug) cout <<"setting tree"<<endl;
   if (datasets.size() == 1)   // myFourTopTree->Write();
    
      combinationsTree_bad->Write();
      combinationsTree_good->Write();
      treeFile->Write();
      treeFile->Close();
    //  delete treeFile;
      
    
    
      //important: free memory
      treeLoader.UnLoadDataset();
      
      
  } //loop on datasets
    
    double n_true_bjets_d = n_true_bjets;
    double n_tagged_true_bjets_d = n_tagged_true_bjets;
    
    double Eff = n_tagged_true_bjets_d/n_true_bjets_d;
    
    double sigEff = (1/n_true_bjets_d)*sqrt(n_tagged_true_bjets_d*(1-Eff));
    
    cout <<" Eff = "  << Eff << " +/- "<< sigEff << endl;
    
    cout <<"N true bjets =   "<< n_true_bjets  <<    "   N tagged true b-jets =  "<<  n_tagged_true_bjets<< endl;

  //Once everything is filled ...
  cout << " We ran over all the data ;-)" << endl;
  
  ///////////////////
  // Writing
  //////////////////
  cout << " - Writing outputs to the files ..." << endl;

  //Selection tables
  if(Muon){ 
	//(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)

      cout <<"  "<<endl;

  //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
  //	selecTableMu.Write(  "FourTop"+postfix+"Table_Mu.tex",    true,true,true,true,false,false,true);
  }
    else if(Electron){
	
  }

    fout->cd();
    
    cout <<"1"<<endl;
    if(histo1D["BJetPt_tagged_TTJets"] && histo1D["BJetPt_TTJets"] && histo1D["BJetPt_tagged_TTTT"] && histo1D["BJetPt_TTTT"]  ){
        cout <<"2"<<endl;
        tgraph["e_b_TTJets"]->BayesDivide(histo1D["BJetPt_tagged_TTJets"],histo1D["BJetPt_TTJets"] );
        tgraph["e_b_TTTT"]->BayesDivide(histo1D["BJetPt_tagged_TTTT"],histo1D["BJetPt_TTTT"] );
        
        tgraph["e_c_TTJets"]->BayesDivide(histo1D["CJetPt_tagged_TTJets"],histo1D["CJetPt_TTJets"] );
        tgraph["e_c_TTTT"]->BayesDivide(histo1D["CJetPt_tagged_TTTT"],histo1D["CJetPt_TTTT"] );
        
       // tgraph["e_l_TTJets"]->BayesDivide(histo1D["LJetPt_tagged_TTJets"],histo1D["LJetPt_TTJets"] );
       // tgraph["e_l_TTTT"]->BayesDivide(histo1D["LJetPt_tagged_TTTT"],histo1D["LJetPt_TTTT"] );
        
        tgraph["e_l_TTJets_eta0-08"]->BayesDivide(histo1D["LJetPt_tagged_TTJets_eta0-08"],histo1D["LJetPt_TTJets_eta0-08"] );
        tgraph["e_l_TTTT_eta0-08"]->BayesDivide(histo1D["LJetPt_tagged_TTTT_eta0-08"],histo1D["LJetPt_TTTT_eta0-08"] );
        
        tgraph["e_l_TTJets_eta08-16"]->BayesDivide(histo1D["LJetPt_tagged_TTJets_eta08-16"],histo1D["LJetPt_TTJets_eta08-16"] );
        tgraph["e_l_TTTT_eta08-16"]->BayesDivide(histo1D["LJetPt_tagged_TTTT_eta08-16"],histo1D["LJetPt_TTTT_eta08-16"] );
    
        tgraph["e_l_TTJets_eta16-24"]->BayesDivide(histo1D["LJetPt_tagged_TTJets_eta16-24"],histo1D["LJetPt_TTJets_eta16-24"] );
        tgraph["e_l_TTTT_eta16-24"]->BayesDivide(histo1D["LJetPt_tagged_TTTT_eta16-24"],histo1D["LJetPt_TTTT_eta16-24"] );
    
    }
    
    
    cout <<"3"<<endl;
    
     TDirectory* tEffdir = fout->mkdir("TGraph");
     tEffdir->cd();
    cout <<"4"<<endl;
    
    for(map<string,TGraphAsymmErrors*>::const_iterator ig = tgraph.begin(); ig != tgraph.end(); ig++)
    {
        
        TGraphAsymmErrors* test = ig->second ;
        
        string grname = ig->first;
        
        string grtitle = "BTagEfficiency_" + grname;
        
        test->SetTitle(grtitle.c_str());
        
        test->Write(grtitle.c_str());
    }

  
  for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++)
    {
	MultiSamplePlot *temp = it->second;
	TH1F *tempHisto_data;
	TH1F *tempHisto_TTTT;
	//	temp->addText("CMS preliminary");
	string name = it->first;
//	temp->Draw(false, name, true, true, true, true, true,1,false); // merge TT/QCD/W/Z/ST/
	//Draw(bool addRandomPseudoData = false, string label = string("CMSPlot"), bool mergeTT = false, bool mergeQCD = false, bool mergeW = false, bool mergeZ = false, bool mergeST = false, int scaleNPSignal = 1, bool addRatio = false, bool mergeVV = false, bool mergeTTV = false);
      
        temp->Write(fout, name, true, pathPNG, "pdf");
      
  }

  cout <<"1D  "<< histo1D.size()  <<"2D   "  <<  histo2D.size() <<endl;
 
   // histo1D["BJetPt"]->Write();
    
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

  //delete  
  delete fout;

  cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;

  cout << "********************************************" << endl;
  cout << "           End of the program !!            " << endl;
  cout << "********************************************" << endl;

  return 0;
}

