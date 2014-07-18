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
//#include "TopTreeAnalysisBase/Selection/interface/FourTopSelectionTable.h"

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
#include "TopTreeAnalysisBase/MCInformation/interface/JetPartonMatching.h"
#include "TopTreeAnalysisBase/Reconstruction/interface/MEzCalculator.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"

#include "TopTreeAnalysisBase/Reconstruction/interface/TTreeObservables.h"

//This header file is taken directly from the BTV wiki. It contains
// to correctly apply an event level Btag SF. It is not yet on CVS
// as I hope to merge the functionality into BTagWeigtTools.h
//#include "TopTreeAnalysisBase/Tools/interface/BTagSFUtil.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"


#include "TopTreeAnalysisBase/Tools/interface/JetCombiner.h"
#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"

using namespace std;
using namespace TopTree;
using namespace reweight;

bool split_ttbar = false;
bool debug = false;


pair<float, vector<unsigned int> > MVAvals1;
pair<float, vector<unsigned int> > MVAvals2;
pair<float, vector<unsigned int> > MVAvals2ndPass;
pair<float, vector<unsigned int> > MVAvals3rdPass;

int nMVASuccesses=0;
int nMatchedEvents=0;

/// Normal Plots (TH1F* and TH2F*)

/// MultiSamplePlot
map<string,MultiSamplePlot*> MSPlot;

/// MultiPadPlot
map<string,MultiSamplePlot*> MultiPadPlot;

struct HighestTCHEBtag {
    bool operator()( TRootJet* j1, TRootJet* j2 ) const {
        return j1->btag_trackCountingHighEffBJetTags() > j2->btag_trackCountingHighEffBJetTags();
        }
    };
struct HighestCVSBtag {
    bool operator()( TRootJet* j1, TRootJet* j2 ) const {
        return j1->btag_combinedSecondaryVertexBJetTags() > j2->btag_combinedSecondaryVertexBJetTags();
        }
    };

bool match;

//To cout the Px, Py, Pz, E and Pt of objects
int Factorial(int N);

int main (int argc, char *argv[]) {

    ofstream eventlist;
    eventlist.open ("interesting_events_mu.txt");

    int passed = 0;
    int ndefs =0;

    string btagger = "CSVM";
    float scalefactorbtageff, mistagfactor;
    float workingpointvalue = 0.679; //working points updated to 2012 BTV-POG recommendations.

    if(btagger == "CSVL")
        workingpointvalue = .244;
    else if(btagger == "CSVM")
        workingpointvalue = .679;
    else if(btagger == "CSVT")
        workingpointvalue = .898;

    clock_t start = clock();

    cout << "*************************************************************" << endl;
    cout << " Beginning of the program for the FourTop search ! "           << endl;
    cout << "*************************************************************" << endl;


    string postfix = "_EventSelection_Dilepton_Unskimmed"; // to relabel the names of the output file

    ///////////////////////////////////////
    // Configuration
    ///////////////////////////////////////

    string channelpostfix = "";
    string xmlFileName = "";

    //Setting Lepton Channels (Setting both flags true will select Muon-Electron Channel)
    bool Muon = true;
    bool Electron = false;

    if(Muon && !Electron) {
        cout << " --> Using the di-Muon channel..." << endl;
        channelpostfix = "_MuMu";
        xmlFileName = "config/test_fullsamples.xml";
        }
    else if(Muon && Electron) {
        cout << " --> Using the Muon-Electron channel..." << endl;
        channelpostfix = "_MuEl";
        xmlFileName = "config/test_fullsamples.xml";
        }
    else if(!Muon && Electron) {
        cout << " --> Using the di-Electron channel..." << endl;
        channelpostfix = "_ElEl";
        xmlFileName = "config/test_fullsamples.xml";
        }
    else {
        cerr<<"Correct Di-lepton Channel not selected."<<endl;
        exit(1);
        }

    bool TrainMVA = false; // If false, the previously trained MVA will be used to calculate stuff
    bool trainEventMVA = false; // If false, the previously trained MVA will be used to calculate stuff
    bool computeEventMVA = false;

    const char *xmlfile = xmlFileName.c_str();
    cout << "used config file: " << xmlfile << endl;

    /////////////////////////////
    //  Set up AnalysisEnvironment
    /////////////////////////////

    AnalysisEnvironment anaEnv;
    cout<<" - Loading environment ..."<<endl;
    AnalysisEnvironmentLoader anaLoad(anaEnv,xmlfile);
    int verbose = 2;//anaEnv.Verbose;

    LumiReWeighting LumiWeights;
    LumiWeights = LumiReWeighting("PUReweighting/pileup_MC_S10.root", "PUReweighting/PURewtoptree_id_2217_2227_2262_2423_2435_2417_PileupHistogram.root", "pileup", "pileup");


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

    cout <<"Instantiating jet combiner..."<<endl;

    JetCombiner* jetCombiner = new JetCombiner(TrainMVA, Luminosity, datasets, MVAmethod, false);
    cout <<"Instantiated jet combiner..."<<endl;
    double bestTopMass1 =0.;
    double bestTopMass2 = 0.;
    double bestTopMass2ndPass = 0.;
    double bestTopPt =0.;

    //uncomment these two lines for training
    //  MVAComputer* Eventcomputer_ =0;
    // MVATrainer* Eventtrainer_ = new MVATrainer("BDT","MasterMVA_El_25thFeb", "MasterMVA_El_25thFeb.root");

    //comment this for training
    MVATrainer* Eventtrainer_ = 0;
    string dataSetName;
    int doJERShift = -1;
    int doJESShift = 0;

    /////////////////////////////////
    //  Loop over Datasets
    /////////////////////////////////

    for (unsigned int d = 0; d < datasets.size (); d++) {
        cout <<"found sample with equivalent lumi "<<  datasets[d]->EquivalentLumi() <<endl;
        dataSetName = datasets[d]->Name();
        if(dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0) {
            Luminosity = datasets[d]->EquivalentLumi();
            cout <<"found DATA sample with equivalent lumi "<<  datasets[d]->EquivalentLumi() <<endl;
            break;
            }
        }

    cout << "Rescaling to an integrated luminosity of "<< Luminosity <<" pb^-1" << endl;
    int ndatasets = datasets.size() - 1 ;

    double frac = 2.;  // frac is the factor by which the TTJets sample is divided, it is minimally 1.4
    double currentLumi;
    double newlumi;
    for (unsigned int d = 0; d < datasets.size (); d++) {
        if( datasets[d]->Name()=="TTJets"  ||datasets[d]->Name()=="TTJets_AllHad" || datasets[d]->Name()=="TTJets_Other") {
            currentLumi = datasets[d]->EquivalentLumi();
            cout <<"Old lumi =   "<< currentLumi  <<endl;
            newlumi = currentLumi/frac;
            datasets[d]->SetEquivalentLuminosity(newlumi);
            }
        }


    // for splitting the ttbar sample, it is essential to have the ttjets sample as the last
    //dataset loaded
    if (split_ttbar) {
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

    MSPlot["RhoCorrection"]                                 = new MultiSamplePlot(datasets, "RhoCorrection", 100, 0, 100, "#rho");
    MSPlot["NbOfVertices"]                                  = new MultiSamplePlot(datasets, "NbOfVertices", 40, 0, 40, "Nb. of vertices");
    //Muons
    MSPlot["NbOfMuonsPreSel"]                               = new MultiSamplePlot(datasets, "NbOfMuonsPreSel", 5, 0, 5, "Nb. of  muons");
    MSPlot["NbOfMuonsPostSel"]                              = new MultiSamplePlot(datasets, "NbOfMuonsPostSel", 5, 0, 5, "Nb. of  muons");
    MSPlot["MuonDRMin"]                                     = new MultiSamplePlot(datasets, "MuonDRMin", 40, 0, 4, "dR_{min}");
    MSPlot["IsolatedMuonDRMin"]                             = new MultiSamplePlot(datasets, "IsolatedMuonDRMin", 40, 0, 4, "dR_{min}");
    MSPlot["MuonRelIsolation"]                              = new MultiSamplePlot(datasets, "MuonRelIsolation", 50, 0, .25, "RelIso");
    MSPlot["MuonPt"]                                        = new MultiSamplePlot(datasets, "MuonPt", 35, 0, 300, "PT_{#mu}");
    MSPlot["MuonPtRel"]                                     = new MultiSamplePlot(datasets, "MuonPtRel", 35, 0, 300, "PT_{#mu}");
    MSPlot["MuonEta"]                                       = new MultiSamplePlot(datasets, "MuonEta", 25, -2.4, 2.4, "#eta_{#mu}");
    MSPlot["MuonPhi"]                                       = new MultiSamplePlot(datasets, "MuonPhi", 24, -M_PI, M_PI, "#phi_{#mu}");
    MSPlot["MuonNValidHits"]                                = new MultiSamplePlot(datasets, "MuonNValidHits", 30, 0, 30, "NValidHits_{#mu}");
    MSPlot["Muond0"]                                        = new MultiSamplePlot(datasets, "Muond0", 50, -0.05, 0.05, "d0_{#mu}");
    MSPlot["MuondZPVz"]                                     = new MultiSamplePlot(datasets, "MuondZPVz", 50, 0, .5, "dZPVZ_{#mu}");
    MSPlot["IsolatedMuMudR"]                                = new MultiSamplePlot(datasets, "IsolatedMuMudR", 25, 0, 5, "dR_{#mu#mu}");
    MSPlot["BadPairMuondR"]                                 = new MultiSamplePlot(datasets, "BadPairMuondR", 20, 0, 0.4, "dR_{#mu#mu}");
    MSPlot["MuonNMatchedStations"]                          = new MultiSamplePlot(datasets, "MuonNMatchedStations", 10, 0, 10, "NMatchedStations_{#mu}");
    MSPlot["MuonDistVzPVz"]                                 = new MultiSamplePlot(datasets, "MuonDistVzPVz", 50, 0 ,.3, "DistVzPVz_{#mu}");
    MSPlot["MuonDz"]                                        = new MultiSamplePlot(datasets, "MuonDz", 12, -.6 ,.6, "Dz_{#mu}");
    MSPlot["MuonTrackerLayersWithMeasurement"]              = new MultiSamplePlot(datasets, "MuonTrackerLayersWithMeasurement", 25, 0, 25, "nLayers");
    MSPlot["DiMuon_InvMass"]                                = new MultiSamplePlot(datasets, "DiMuon_InvMass", 60, 0, 120, "DiMuon_InvMass");
    MSPlot["BadDiMuon_InvMass"]                             = new MultiSamplePlot(datasets, "BadDiMuon_InvMass", 20, 0, 10, "DiMuon_InvMass");
    MSPlot["BadPairMuonDEta"]                               = new MultiSamplePlot(datasets, "BadPairMuonDEta", 20, 0, 0.2, "abs(#dEta_{#mu})");
    MSPlot["BadPairMuonDPhi"]                               = new MultiSamplePlot(datasets, "BadPairMuonDPhi", 20, 0, 0.2, "abs(#dPhi_{#mu})");
    MSPlot["NbOfLooseMuons"]                                = new MultiSamplePlot(datasets, "NbOfLooseMuons", 10, 0, 10, "Nb. of loose muons");
    MSPlot["MuonPairCharge"]                                = new MultiSamplePlot(datasets, "MuonPairCharge", 5, -2.5, 2.5, "Total charge");

    //Electrons
    MSPlot["NbOfIsolatedElectrons"]                         = new MultiSamplePlot(datasets, "NbOfIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
    MSPlot["NbOfExtraIsolatedElectrons"]                    = new MultiSamplePlot(datasets, "NbOfExtraIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
    MSPlot["ElectronRelIsolation"]                          = new MultiSamplePlot(datasets, "ElectronRelIsolation", 50, 0, .25, "RelIso");
    MSPlot["ElectronPt"]                                    = new MultiSamplePlot(datasets, "ElectronPt", 35, 0, 300, "PT_{#mu}");
    MSPlot["ElectronEta"]                                   = new MultiSamplePlot(datasets, "ElectronEta", 25, -2.4, 2.4, "#eta_{#mu}");
    MSPlot["ElectronPhi"]                                   = new MultiSamplePlot(datasets, "ElectronPhi", 50, -4, 4, "#phi_{#mu}");
    MSPlot["Electrond0"]                                    = new MultiSamplePlot(datasets, "Electrond0", 50, -0.05, 0.05, "d0_{#mu}");
    MSPlot["ElectrondZPVz"]                                 = new MultiSamplePlot(datasets, "ElectrondZPVz", 50, 0, .5, "dZPVZ_{#mu}");
    MSPlot["ElectrondRJets"]                                = new MultiSamplePlot(datasets, "ElectrondRJets", 50, 0, 10, "dRJets_{#mu}");
    MSPlot["ElectronDistVzPVz"]                             = new MultiSamplePlot(datasets, "ElectronDistVzPVz", 50, 0 ,.3, "DistVzPVz_{#mu}");
    MSPlot["ElectronDz"]                                    = new MultiSamplePlot(datasets, "ElectronDz", 25, -.6 ,.6, "dZ_{#mu}");
    MSPlot["DiElectron_InvMass"]                            = new MultiSamplePlot(datasets, "DiElectron_InvMass", 60, 0, 120, "DiElectron_InvMass");
    MSPlot["ElectronPairCharge"]                            = new MultiSamplePlot(datasets, "ElectronPairCharge", 5, -2.5, 2.5, "Total charge");

    //Plots Specific to MuEl channel
    MSPlot["MuElPairCharge"]                                = new MultiSamplePlot(datasets, "MuElPairCharge", 5, -2.5, 2.5, "Total charge");

    //B-tagging discriminators
    MSPlot["BdiscBJetCand_CSV"]                             = new MultiSamplePlot(datasets, "BdiscBJetCand_CSV", 75, 0, 1, "CSV b-disc.");
    //Jets
    MSPlot["NbOfSelectedJets"]                              = new MultiSamplePlot(datasets, "NbOfSelectedJets", 15, 0, 15, "Nb. of jets");
    MSPlot["NbOfSelectedLightJets"]                         = new MultiSamplePlot(datasets, "NbOfSelectedLightJets", 10, 0, 10, "Nb. of jets");
    MSPlot["NbOfSelectedMBJets"]                            = new MultiSamplePlot(datasets, "NbOfSelectedMBJets", 8, 0, 8, "Nb. of CSVM tags");
    MSPlot["JetEta"]                                        = new MultiSamplePlot(datasets, "JetEta", 30,-3, 3, "Jet #eta");
    MSPlot["JetPhi"]                                        = new MultiSamplePlot(datasets, "JetPhi", 50, -M_PI,M_PI , "Jet #phi");
    MSPlot["SelectedJetPt"]                                 = new MultiSamplePlot(datasets, "SelectedJetPt", 50, 0, 300, "PT_{jet}");
    MSPlot["3rdJetPt"]                                      = new MultiSamplePlot(datasets, "3rdJetPt", 40, 0, 400, "PT_{jet}");
    MSPlot["4thJetPt"]                                      = new MultiSamplePlot(datasets, "4thJetPt", 40, 0, 400, "PT_{jet}");
    MSPlot["5thJetPt"]                                      = new MultiSamplePlot(datasets, "5thJetPt", 40, 0, 400, "PT_{jet}");
    MSPlot["6thJetPt"]                                      = new MultiSamplePlot(datasets, "6thJetPt", 40, 0, 400, "PT_{jet}");
    MSPlot["7thJetPt"]                                      = new MultiSamplePlot(datasets, "7thJetPt", 40, 0, 400, "PT_{jet}");
    MSPlot["SelectedJetPt_light"]                           = new MultiSamplePlot(datasets, "SelectedJetPt_light", 50, 0, 1000, "PT_{lightjet}");
    MSPlot["HT_SelectedJets"]                               = new MultiSamplePlot(datasets, "HT_SelectedJets", 30, 0, 1500, "HT");
    MSPlot["HTb_SelectedJets"]                              = new MultiSamplePlot(datasets, "HTb_SelectedJets", 30, 0, 1500, "HTb");
    MSPlot["H"]                                             = new MultiSamplePlot(datasets, "H", 30, 0, 3000, "H");
    MSPlot["HTH"]                                           = new MultiSamplePlot(datasets, "HTH", 25, 0, 1, "HT/H");
    MSPlot["HTRat"]                                         = new MultiSamplePlot(datasets, "HTRat", 50, 0, 20, "HTRat");
    MSPlot["HTExcess2M"]                                    = new MultiSamplePlot(datasets, "HTExcess2M", 50, 0, 1000, "HT");
    MSPlot["HExcess2M"]                                     = new MultiSamplePlot(datasets, "HExcess2M", 50, 0, 3000, "H");
    MSPlot["HTExcess1M2L"]                                  = new MultiSamplePlot(datasets, "HTExcess1M2L", 50, 0, 1000, "HT");
    MSPlot["HExcess1M2L"]                                   = new MultiSamplePlot(datasets, "HExcess1M2L", 50, 0, 3000, "H");

    //MET
    MSPlot["MET"]                                           = new MultiSamplePlot(datasets, "MET", 75, 0, 700, "MET");

    //Truth Plots
    MSPlot["LeptonTruth"]                                   = new MultiSamplePlot(datasets, "LeptonTruth", 6, -1, 5, "LeptonTruth");
    MSPlot["SameSignMuonType"]                              = new MultiSamplePlot(datasets, "SameSignMuonType", 26, -1, 25, "Muon Type Code");
    MSPlot["SameSignMuonParentType"]                        = new MultiSamplePlot(datasets, "SameSignMuonParentType", 26, -1, 25, "Parent Type Code");
    MSPlot["SameSignMuonGrannyType"]                        = new MultiSamplePlot(datasets, "SameSignMuonGrannyType", 26, -1, 25, "Granny Type Code");
    MSPlot["SameSignTotalCharge"]                           = new MultiSamplePlot(datasets, "SameSignTotalCharge", 7, -3.5, 3.5, "Total Charge in Muons");
    MSPlot["SameSignNbOfMuons"]                             = new MultiSamplePlot(datasets, "SameSignNbOfMuons", 6, 0, 6, "Number of Muons");
    MSPlot["OpSignMuonType"]                                = new MultiSamplePlot(datasets, "OpSignMuonType", 26, -1, 25, "Muon Type Code");
    MSPlot["OpSignMuonParentType"]                          = new MultiSamplePlot(datasets, "OpSignMuonParentType", 26, -1, 25, "Parent Type Code");
    MSPlot["OpSignMuonGrannyType"]                          = new MultiSamplePlot(datasets, "OpSignMuonGrannyType", 26, -1, 25, "Granny Type Code");

    //MVA Top Roconstruction Plots
    MSPlot["MVA1TriJet"]                                    = new MultiSamplePlot(datasets, "MVA1TriJet", 30, -1.0, 0.2, "MVA1TriJet");
    MSPlot["MVA1TriJetMass"]                                = new MultiSamplePlot(datasets, "MVA1TriJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA1DiJetMass"]                                 = new MultiSamplePlot(datasets, "MVA1DiJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA1PtRat"]                                     = new MultiSamplePlot(datasets, "MVA1PtRat", 25, 0, 2, "P_{t}^{Rat}");
    MSPlot["MVA1BTag"]                                      = new MultiSamplePlot(datasets, "MVA1BTag", 35, 0, 1, "BTag");
    MSPlot["MVA1AnThBh"]                                    = new MultiSamplePlot(datasets, "MVA1AnThBh", 35, 0, 3.14, "AnThBh");
    MSPlot["MVA1AnThWh"]                                    = new MultiSamplePlot(datasets, "MVA1AnThWh", 35, 0, 3.14, "AnThWh");

    //Declare arrays of MSPlots
    Int_t minNJets=6, maxNJets=8, minNBJets=2, maxNBJets=3;
    Int_t q =0;
    for (Int_t p = minNJets; p<= maxNJets; p++) {
        for (Int_t q = minNBJets; q<= maxNBJets; q++) {
            string NJets_str = static_cast<ostringstream*>( &(ostringstream() << p) )->str();
            string NBJets_str = static_cast<ostringstream*>( &(ostringstream() << q) )->str();
            string HT_Name = "HT_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
            string H_Name = "H_"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
            string MET_Name = "MET"+NJets_str+"Jets_"+NBJets_str+"Tags" ;
            }
        }


    ///////////////////
    // 1D histograms
    ///////////////////

    //Plots
    string pathPNG = "FourTop"+postfix+channelpostfix;
    pathPNG += "_MSPlots/";
    //pathPNG = pathPNG +"/";
    mkdir(pathPNG.c_str(),0777);

    cout <<"Making directory :"<< pathPNG  <<endl;

    /////////////////////////////////
    // Selection table: Dilepton + jets
    /////////////////////////////////
    vector<string> CutsselecTable;
    if(Muon && !Electron) {
        CutsselecTable.push_back(string("initial"));
        CutsselecTable.push_back(string("Event cleaning and Trigger"));
        CutsselecTable.push_back(string("At least 2 Loose Isolated Muons (relIso04$\\leq 0.20$)"));
        CutsselecTable.push_back(string("Highest 2 pT Muon Mass $\\geq 20 GeV$"));
//        CutsselecTable.push_back(string("At least 1 CSVM Jet"));
//        CutsselecTable.push_back(string("At least 1 extra CSVL Jet"));
//        CutsselecTable.push_back(string("At least 2 extra CSVL Jet"));
        CutsselecTable.push_back(string("At least 1  Jet"));
        CutsselecTable.push_back(string("At least 2 Jets"));
        CutsselecTable.push_back(string("At least 3 Jets"));
        CutsselecTable.push_back(string("Same Sign Leading Leptons"));
        CutsselecTable.push_back(string("Same Sign $+$ Extra Lepton"));
        }
    if(Muon && Electron) {
        CutsselecTable.push_back(string("initial"));
        CutsselecTable.push_back(string("Event cleaning and Trigger"));
        CutsselecTable.push_back(string("At least 1 Loose Isolated Muon (relIso04$\\leq 0.20$)"));
        CutsselecTable.push_back(string("At least 1 Loose Electron"));
        CutsselecTable.push_back(string("At least 1 CSVM Jet"));
        CutsselecTable.push_back(string("At least 1 extra CSVL Jet"));
        CutsselecTable.push_back(string("At least 2 extra CSVL Jet"));
        CutsselecTable.push_back(string("Same Sign Leading Leptons"));
        }
    if(!Muon && Electron) {
        CutsselecTable.push_back(string("initial"));
        CutsselecTable.push_back(string("Event cleaning and Trigger"));
        CutsselecTable.push_back(string("Exactly 0 Loose Muon"));
        CutsselecTable.push_back(string("Exactly 2 Loose Electron"));
        CutsselecTable.push_back(string("At least 1 CSVM Jet"));
        CutsselecTable.push_back(string("At least 1 extra CSVL Jet"));
        CutsselecTable.push_back(string("At least 2 extra CSVL Jet"));
        CutsselecTable.push_back(string("Same Sign Leading Leptons"));
        }

    SelectionTable selecTable(CutsselecTable, datasets);
    selecTable.SetLuminosity(Luminosity);
    selecTable.SetPrecision(1);

    /////////////////////////////////
    // Loop on datasets
    /////////////////////////////////

    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;

    for (unsigned int d = 0;
            d < datasets.size();
            d++) {
        cout<<"Load Dataset"<<endl;
        treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset
        string previousFilename = "";
        int iFile = -1;
        dataSetName = datasets[d]->Name();

        //////////////////////////////////////////////////
        // Initialize JEC factors ///////////////////////
        //////////////////////////////////////////////////

        vector<JetCorrectorParameters> vCorrParam;

        if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) { // Data!
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L1FastJet_AK5PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L2Relative_AK5PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L3Absolute_AK5PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L2L3Residual_AK5PFchs.txt");
            vCorrParam.push_back(*L2L3ResJetCorPar);
            }
        else {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_L1FastJet_AK5PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_L2Relative_AK5PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_L3Absolute_AK5PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            }
        JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_Uncertainty_AK5PFchs.txt");
        JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);


        //////////////////////////////////////////////////
        // Loop on events
        /////////////////////////////////////////////////

        int itrigger = -1, previousRun = -1;

        int start = 0;
        unsigned int end = datasets[d]->NofEvtsToRunOver();

        cout <<"Number of events = "<<  end  <<endl;

        int event_start = 0;
        if (verbose > 1) cout << " - Loop over events " << endl;

        double MHT, MHTSig, STJet, EventMass, EventMassX , SumJetMass, SumJetMassX,H,HX , HT, HTX,HTH,HTXHX, sumpx_X, sumpy_X, sumpz_X, sume_X, sumpx, sumpy, sumpz, sume, jetpt,PTBalTopEventX,PTBalTopSumJetX , PTBalTopMuMet;

        double currentfrac =0.;
        double end_d = end;

        cout <<"Will run over "<<  end_d<< " events..."<<endl;
        cout <<"Starting event = = = = "<< event_start  << endl;

        //define object containers
        vector<TRootElectron*> selectedElectrons;
        vector<TRootJet*>      selectedJets;
        vector<TRootMuon*>     selectedMuons;
        vector<TRootMuon*>     selectedLooseMuons;
        vector<TRootElectron*> selectedLooseElectrons;
        vector<TRootElectron*> selectedExtraElectrons;
        vector<TRootMuon*>     selectedMuons_NoIso;
        vector<TRootMuon*>     selectedExtraMuons;
        selectedElectrons.reserve(10);
        selectedMuons.reserve(10);
        selectedLooseMuons.reserve(10);
        selectedLooseElectrons.reserve(10);

        //////////////////////////////////////
        // Begin Event Loop
        //////////////////////////////////////

        for (unsigned int ievt = event_start; ievt < end_d; ievt++) {
            MHT = 0.,MHTSig = 0., STJet = 0., EventMass =0., EventMassX =0., SumJetMass = 0., SumJetMassX=0.  ,H = 0., HX =0., HT = 0., HTX = 0.,HTH=0.,HTXHX=0., sumpx_X = 0., sumpy_X= 0., sumpz_X =0., sume_X= 0. , sumpx =0., sumpy=0., sumpz=0., sume=0., jetpt =0., PTBalTopEventX = 0., PTBalTopSumJetX =0.;

            double ievt_d = ievt;
            currentfrac = ievt_d/end_d;
            if (debug)cout <<"event loop 1"<<endl;

            if(ievt%1000 == 0)
                std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r"<<endl;

            float scaleFactor = 1.;  // scale factor for the event
            event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets);  //load event

            float rho = event->kt6PFJets_rho();
            string graphName;

            /////////////////////////
            // Loading Gen jets
            /////////////////////////

            vector<TRootGenJet*> genjets;
            if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) ) {
                // loading GenJets as I need them for JER
                genjets = treeLoader.LoadGenJet(ievt);
                }
            // check which file in the dataset it is to have the HLTInfo right
            string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
            if(previousFilename != currentFilename) {
                previousFilename = currentFilename;
                iFile++;
                cout<<"File changed!!! => iFile = "<<iFile<<endl;
                }

            ///////////////////////////////////////////
            // Trigger
            ///////////////////////////////////////////
            bool trigged = false;
            std::string filterName = "";
            int currentRun = event->runId();
            if(previousRun != currentRun) {
                // cout <<"What run? "<< currentRun<<endl;
                previousRun = currentRun;
                if(Muon && !Electron) { // di-muon
                    if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) {
                        if( event->runId() >= 190456 && event->runId() <= 190738 ) {
                            itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v11"), currentRun, iFile);
                            filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10";
                            }
                        else if( event->runId() >= 190782 && event->runId() <= 193621) {
                            itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v12"), currentRun, iFile);
                            filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10";
                            }
                        else if(event->runId() >= 193834  && event->runId() <= 196531 ) {
                            itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
                            filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";
                            }
                        else if( event->runId() >= 198022  && event->runId() <= 199608) {
                            itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v14"), currentRun, iFile);
                            filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";
                            }
                        else if( event->runId() >= 199698 && event->runId() <= 209151) {
                            itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v15"), currentRun, iFile);
                            filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";
                            }
                        else {
                            cout << "Unknown run for SemiMu HLTpath selection: " << event->runId() << endl;
                            filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";

                            }
                        if( itrigger == 9999 ) {
                            cout << "itriggerSemiMu == 9999 for SemiMu HLTpath selection: " << event->runId() << endl;
                            //exit(-1);
                            }
                        }
                    else if(Muon && Electron) { // Muon-Electron
                        }
                    else if(!Muon && Electron) { // di-Electron
                        }
                    else {
                        if(dataSetName == "TTJets" ) itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);

                        else {

                            itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
                            }

                        if(itrigger == 9999) {
                            cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;
                            //exit(1);
                            }
                        }
                    } //end Muon Triggering
                } //end previousRun != currentRun

            ////////////////////////////////////////////////////////////////////////////////////
            // JES Corrections: The nominal corrections are already applied at PAT level     //
            // so these tools should only be used for studies of the effect of systematics //
            ////////////////////////////////////////////////////////////////////////////////////

            // Apply Jet Corrections on-the-fly
            //coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction on the fly:");
            //	if( dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" )
            //		jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),true); //last boolean: isData (needed for L2L3Residual...)
            //	else
            //		jetTools->correctJets(init_jets,event->kt6PFJetsPF2PAT_rho(),false); //last boolean: isData (needed for L2L3Residual...)
            //	coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JES correction on the fly:");

            // Apply Jet Corrections on-the-fly
//            if( dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) {
//                jetTools->unCorrectMETTypeOne(init_jets, mets[0], true);
//                jetTools->correctJets(init_jets, event->kt6PFJets_rho(), true);
//                jetTools->correctMETTypeOne(init_jets, mets[0], true);
//                }
//            else {
//                if (debug)cout <<"event loop 8"<<endl;
//                jetTools->unCorrectMETTypeOne(init_jets, mets[0], false);
//                jetTools->correctJets(init_jets, event->kt6PFJets_rho(), false);
//                jetTools->correctMETTypeOne(init_jets, mets[0], false);
//                }

            ///////////////////////
            // JER smearing
            //////////////////////

            if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) ) {
                //JER
                if (debug) cout << "Doing JER Corrections: ";
                doJERShift = -1;
                if(doJERShift == 1)
                    jetTools->correctJetJER(init_jets, genjets, mets[0], "minus");
                else if(doJERShift == 2)
                    jetTools->correctJetJER(init_jets, genjets, mets[0], "plus");
                else if(doJERShift == 0) {
                    if (debug) cout << "Nominal" << endl;
                    jetTools->correctJetJER(init_jets, genjets, mets[0], "nominal");
                    }
                else if (debug) cout << "JER Disabled." << endl;

                //     coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"After JER correction:");

                // JES sysematic!
                if (doJESShift == 1)
                    jetTools->correctJetJESUnc(init_jets, mets[0], "minus");
                else if (doJESShift == 2)
                    jetTools->correctJetJESUnc(init_jets, mets[0], "plus");

                //            coutObjectsFourVector(init_muons,init_electrons,init_jets,mets,"Before JES correction:");

                }
            if (debug) {
                cout<<"Debugging Jet Momentum"<<endl;
                for(int curJet = 0; curJet<init_jets.size(); curJet++) {
                    cout << "Px: "<<init_jets[curJet]->Px()<<" Py: "<<init_jets[curJet]->Py()<<" Pz: "<<init_jets[curJet]->Pz()<<" PT: "<<init_jets[curJet]->Pt()<<endl;
                    cout<<"Eta: "<<init_jets[curJet]->Eta()<<endl;
                    }
                cout<<"Debugging Gen Jet Momentum"<<endl;
                for(int curJet = 0; curJet<genjets.size(); curJet++) {
                    cout << "Px: "<<genjets[curJet]->Px()<<" Py: "<<genjets[curJet]->Py()<<" Pz: "<<genjets[curJet]->Pz()<<" E: "<<genjets[curJet]->E()<<" PT: "<<genjets[curJet]->Pt()<<endl;
                    cout<<"Eta: "<<genjets[curJet]->Eta()<<endl;
                    }
                cin.get();
                }



            ////////////////////////////////////////
            // Beam scraping and PU reweighting
            ////////////////////////////////////////

            // scale factor for the event
            //float scaleFactor = 1.;

            if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") {
                // Apply the scraping veto. (Is it still needed?)
                bool isBeamBG = true;
                if(event->nTracks() > 10) {
                    if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
                        isBeamBG = false;
                    }
                if(isBeamBG) continue;
                }
            else {
                double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
                double lumiWeightOLD=lumiWeight;
                if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) {
                    lumiWeight=1;
                    }
                scaleFactor = scaleFactor*lumiWeight;

                }

            ///////////////////////////////////////////////////////////
            // Event selection
            ///////////////////////////////////////////////////////////

            // Apply trigger selection
//            trigged = treeLoader.EventTrigged (itrigger);
            trigged = true;  // Disabling the HLT requirement
            if (debug)cout<<"triggered? Y/N?  "<< trigged  <<endl;
            if(!trigged)		   continue;  //If an HLT condition is not present, skip this event in the loop.
            // Declare selection instance
            Selection selection(init_jets, init_muons, init_electrons, mets);

            // Define object selection cuts
            selection.setJetCuts(30.,2.5,0.01,1.,0.98,0,0);//Pt, Eta, EMF, n90Hits, fHPD, dRJetElectron, DRJets
            selection.setDiElectronCuts(20,2.5,.15,.04,.5,1,0,0);//	selection.setElectronCuts(30.,2.5,0.1,0.02,0.,999.,0,1); //Pt,Eta,RelIso,d0,MVAId,DistVzPVz, DRJets, MaxMissingHits
            selection.setLooseElectronCuts(10,2.5,0.15,0);
            selection.setDiMuonCuts();
            selection.setLooseMuonCuts();
            if (debug)cout<<"Getting Tight Electrons"<<endl;
            selectedElectrons        = selection.GetSelectedDiElectrons();
            if (debug)cout<<"Getting Jets"<<endl;
            selectedJets             = selection.GetSelectedJets(true); // ApplyJetId
            if (debug)cout<<"Getting Tight Muons"<<endl;
            selectedMuons            = selection.GetSelectedDiMuons(20,2.4,0.20);
            if (debug)cout<<"Getting Loose Muons"<<endl;
            selectedLooseMuons       = selection.GetSelectedLooseMuons();
            if (debug)cout<<"Getting Loose Electrons"<<endl;
            selectedLooseElectrons   = selection.GetSelectedLooseElectrons(10,2.5,0.15); // VBTF ID

            vector<TRootJet*>      selectedMBJets;
            vector<TRootJet*>      selectedLBJets;
            vector<TRootJet*>      selectedLightJets;
            vector<TRootJet*>       MVASelJets1;

            //order jets wrt to Pt, then set bool corresponding to RefSel cuts.
            //sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order muons wrt Pt.
            //sort(selectedCSVOrderedJets.begin(), selectedCSVOrderedJets.end(), HighestCVSBtag()); //order Jets wrt CSVtag

            ////////////////////////////////////////////////////////////
            // Calculating dRMin and pTRel for each muon in the event //
            // This is done here for muon isolation selection         //
            ////////////////////////////////////////////////////////////

            double dR = 0, dPhi = 0, dEta = 0, lepSq = 0, projSq = 0, jetSq = 0, diMuMass = 0;
            double dRTemp = 0, pTRelTemp = 0;
            float pairCharge = -9999, totalCharge = 0;
            vector<double>  MuondRMin;
            vector<double>  MuonPTRel;
            int jetPosMin = 0, nIsoMu = 0;

            vector<TLorentzVector> selectedMuonsTLV_JC, selectedElectronsTLV_JC;
            vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV, mcMuonsTLV;
            vector<TRootMCParticle*> mcParticlesMatching_;
            vector<int> mcMuonIndex;
            JetPartonMatching muonMatching;

            for (Int_t selmu =0; selmu < selectedMuons.size(); selmu++ ) {
                MuondRMin.push_back(99999);
                selectedMuonsTLV_JC.push_back(*selectedMuons[selmu]);
                for (Int_t seljet =0; seljet<selectedJets.size(); seljet++) {
                    dRTemp = abs(selectedMuons[selmu]->DeltaR(*selectedJets[seljet]));
                    if (dRTemp <= MuondRMin[selmu]) {
                        MuondRMin[selmu] = dRTemp;
                        jetPosMin = seljet;
                        }
                    }
                lepSq = pow(selectedMuons[selmu]->Px(),2) + pow(selectedMuons[selmu]->Py(),2) + pow(selectedMuons[selmu]->Pz(),2); //square of lepton momentum
                projSq = pow(selectedMuons[selmu]->Px() * selectedJets[jetPosMin]->Px() + selectedMuons[selmu]->Py() * selectedJets[jetPosMin]->Py() + selectedMuons[selmu]->Pz() * selectedJets[jetPosMin]->Pz(),2); //square of lepton momentum projected onto jet axis
                jetSq = pow(selectedJets[jetPosMin]->Px(),2) + pow(selectedJets[jetPosMin]->Py(),2) + pow(selectedJets[jetPosMin]->Pz(),2); //square of jet momentum
                pTRelTemp = sqrt( lepSq - (projSq/jetSq) );
                MuonPTRel.push_back(pTRelTemp);
                }
            if(selectedMuons.size() >=2) {
                TLorentzVector diMuon = selectedMuonsTLV_JC[0] + selectedMuonsTLV_JC[1];
                diMuMass = diMuon.M();
                }

            int JetCut =0;
            int nMu, nEl;
            if(Muon && !Electron) {
                nMu = selectedMuons.size(); //Number of Muons in Event
                nEl = selectedElectrons.size(); //Number of Electrons in Event
                }
            if(Muon && Electron) {
                nMu = selectedMuons.size(); //Number of Muons in Event
                nEl = selectedElectrons.size(); //Number of Electrons in Event
                }
            if(!Muon && Electron) {
                nMu = selectedLooseMuons.size(); //Number of Muons in Event
                nEl = selectedElectrons.size(); //Number of Electrons in Event
                }

            bool isTagged =false;
            int seljet;

            ////////////////////////////////////////////////////
            // Looping over Jets to populate Bjet Collections //
            ////////////////////////////////////////////////////
            for ( seljet =0; seljet < selectedJets.size(); seljet++ ) {
                if (selectedJets[seljet]->btag_combinedSecondaryVertexBJetTags() > 0.244   ) {
                    selectedLBJets.push_back(selectedJets[seljet]);
                    if (selectedJets[seljet]->btag_combinedSecondaryVertexBJetTags() > workingpointvalue) {
                        selectedMBJets.push_back(selectedJets[seljet]);
                        }
                    }
                else {
                    selectedLightJets.push_back(selectedJets[seljet]);
                    }
                }
            int nJets = selectedJets.size(); //Number of Jets in Event
            int nMtags = selectedMBJets.size(); //Number of CSVM tags in Event
            int nLtags = selectedLBJets.size(); //Number of CSVL tags in Event (includes jets that pass CSVM)
            ///////////////////////////////////////////////////////////////////////////////////
            // Summing HT and calculating leading, lagging, and ratio for Selected and BJets //
            ///////////////////////////////////////////////////////////////////////////////////
            double temp_HT = 0.;
            double HTRat = 0.;
            double HT_leading = 0.;
            double HT_lagging = 0.;
            for (Int_t seljet0 =0; seljet0 < selectedJets.size(); seljet0++ ) {
                temp_HT += selectedJets[seljet0]->Pt();
                if (seljet0 < 4) { //Defines the leading Jets and the first 4
                    HT_leading += selectedJets[seljet0]->Pt();
                    }
                else {
                    HT_lagging += selectedJets[seljet0]->Pt();
                    }
                }
            HTRat = HT_leading/HT_lagging;
            double HTb = 0.;
            for (Int_t seljet1 =0; seljet1 < nMtags; seljet1++ ) {
                HTb += selectedMBJets[seljet1]->Pt();
                }

            //////////////////////
            // Sync'ing cutflow //
            //////////////////////

            if (debug)	cout <<" applying baseline event selection for cut table..."<<endl;
            // Apply primary vertex selection
            bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
            if (debug)	cout <<"PrimaryVertexBit: " << isGoodPV << " TriggerBit: " << trigged <<endl;
            if (debug) cin.get();
            selecTable.Fill(d,0,1);
            if(Muon && !Electron) { //Di-Muon Selection Table
                if(isGoodPV && trigged) {
                    selecTable.Fill(d,1,1);
                    if (nMu>=2 ) {
                        selecTable.Fill(d,2,1);
                        if(diMuMass >= 20) {
                            selecTable.Fill(d,3,1);
//                            if(nMtags>=1) {
//                                selecTable.Fill(d,4,1);
//                                if(nLtags>=2) {
//                                    selecTable.Fill(d,5,1);
//                                    if(nLtags>=3) {
//                                        selecTable.Fill(d,6,1);
                            if(nJets>=1) {
                                selecTable.Fill(d,4,1);
                                if(nJets>=2) {
                                    selecTable.Fill(d,5,1);
                                    if(nJets>=3) {
                                        selecTable.Fill(d,6,1);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            if(Muon && Electron) { //Muon-Electron Selection Table
                if(isGoodPV && trigged) {
                    selecTable.Fill(d,1,1);
                    if (nMu>=1) {
                        selecTable.Fill(d,2,1);
                        if(nEl>=1) {
                            selecTable.Fill(d,3,1);
                            if(nMtags>=1) {
                                selecTable.Fill(d,4,1);
                                if(nLtags>=2) {
                                    selecTable.Fill(d,5,1);
                                    if(nLtags>=3) {
                                        selecTable.Fill(d,6,1);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            if(!Muon && Electron) { //Di-Electron Selection Table
                if(isGoodPV && trigged) {
                    selecTable.Fill(d,1,1);
                    if (nMu==0) {
                        selecTable.Fill(d,2,1);
                        if(nEl==2) {
                            selecTable.Fill(d,3,1);
                            if(nMtags>=1) {
                                selecTable.Fill(d,4,1);
                                if(nLtags>=2) {
                                    selecTable.Fill(d,5,1);
                                    if(nLtags>=3) {
                                        selecTable.Fill(d,6,1);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            /////////////////////////////////
            // Applying baseline selection //
            /////////////////////////////////
//            if (!trigged) continue;  // Redunant check that an HLT was triggered
            if (!isGoodPV) continue; // Check that there is a good Primary Vertex
////            if (!(selectedJets.size() >= 6)) continue; //Selection of a minimum of 6 Jets in Event
//
//            if (debug) cout <<"Number of Muons, Electrons, Jets, BJets, JetCut, MuonChannel, ElectronChannel ===>  "<< nMu <<"  "  <<nEl<<" "<< selectedJets.size()   <<"  " <<  nMtags   <<"  "<<JetCut  <<"  "<<Muon<<" "<<Electron<<endl;
//            if (debug) cin.get();
//            int nTightLeptons, nVetoLeptonsSF, nVetoLeptonsOF;
//            if (debug)	cout <<" applying baseline event selection..."<<endl;
//            //Apply the lepton, btag and HT selections
//            if (Muon && !Electron) {
//                MSPlot["NbOfMuonsPreSel"]->Fill(nMu, datasets[d], true, Luminosity*scaleFactor);
//                if  (  !( nMu >= 2 && diMuMass >= 20) ) continue; // Di-Muon Channel Selection with dirty muon iso.
//                }
//            else if (Muon && Electron) {
//                if  (  !( nMu >= 1 && nEl >= 1 )) continue; // Muon-Electron Channel Selection
//                }
//            else if (!Muon && Electron) {
//                if  (  !( nEl == 2 )) continue; // Di-Electron Channel Selection
//                }
//            else {
//                cerr<<"Correct Di-lepton Channel not selected."<<endl;
//                exit(1);
//                }
            sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());

            if(nMtags>=2) {
                for (Int_t seljet1 =2; seljet1 < selectedJets.size(); seljet1++ ) {

                    jetpt = selectedJets[seljet1]->Pt();
                    HT = HT + jetpt;
                    H = H + selectedJets[seljet1]->P();
                    }
                MSPlot["HTExcess2M"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["HExcess2M"]->Fill(H, datasets[d], true, Luminosity*scaleFactor);
                }


            //if (!(nMtags>=1 && nLtags>=3)) continue; //Jet Tag Event Selection Requirements
//            if (!(nJets>=3)) continue; //Jet Tag Event Selection Requirements
//            if(debug) {
//                cout<<"Selection Passed."<<endl;
//                cin.get();
//                }
            passed++;


            ///////////////////////
            // Getting Gen Event //
            ///////////////////////

            TRootGenEvent* genEvt = 0;

            if(dataSetName != "data" && dataSetName != "Data" && dataSetName != "Data") {
                vector<TRootMCParticle*> mcParticles;
                vector<TRootMCParticle*> mcTops;
                mcParticlesMatching_.clear();
                mcParticlesTLV.clear();
                selectedJetsTLV.clear();
                mcParticles.clear();
                mcTops.clear();

                int leptonPDG, muonPDG = 13, electronPDG = 11;
                leptonPDG = muonPDG;

                genEvt = treeLoader.LoadGenEvent(ievt,false);
                treeLoader.LoadMCEvent(ievt, genEvt, 0, mcParticlesMatching_,false);
                if (debug) cout <<"size   "<< mcParticlesMatching_.size()<<endl;
                }

            ///////////////////////////
            // Event Level Variables //
            ///////////////////////////

//            //////////////////////////////////////
//            // MVA Hadronic Top Reconstructions //
//            //////////////////////////////////////
//
//            jetCombiner->ProcessEvent_SingleHadTop(datasets[d], mcParticlesMatching_, selectedJets, &selectedMuonsTLV_JC[0], genEvt, scaleFactor);
//
//            if(!TrainMVA) {
//                MVAvals1 = jetCombiner->getMVAValue(MVAmethod, 1); // 1 means the highest MVA value
//                MSPlot["MVA1TriJet"]->Fill(MVAvals1.first, datasets[d], true, Luminosity*scaleFactor );
//
//                for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ ) {
//                    if (seljet1 == MVAvals1.second[0] || seljet1 == MVAvals1.second[1] || seljet1 == MVAvals1.second[2]) {
//                        MVASelJets1.push_back(selectedJets[seljet1]);
//                        }
//
//                    }
//
//                //check data-mc agreement of kin. reco. variables.
//                float mindeltaR =100.;
//                float mindeltaR_temp =100.;
//                int wj1;
//                int wj2;
//                int bj1;
//
//                //define the jets from W as the jet pair with smallest deltaR
//                for (int m=0; m<MVASelJets1.size(); m++) {
//                    for (int n=0; n<MVASelJets1.size(); n++) {
//                        if(n==m) continue;
//                        TLorentzVector lj1  = *MVASelJets1[m];
//                        TLorentzVector lj2  = *MVASelJets1[n];
//                        mindeltaR_temp  = lj1.DeltaR(lj2);
//                        if (mindeltaR_temp < mindeltaR) {
//                            mindeltaR = mindeltaR_temp;
//                            wj1 = m;
//                            wj2 = n;
//                            }
//                        }
//                    }
//                // find the index of the jet not chosen as a W-jet
//                for (unsigned int p=0; p<MVASelJets1.size(); p++) {
//                    if(p!=wj1 && p!=wj2) bj1 = p;
//                    }
//
//                if (debug) cout <<"Processing event with jetcombiner : 3 "<< endl;
//
//                //now that putative b and W jets are chosen, calculate the six kin. variables.
//                TLorentzVector Wh = *MVASelJets1[wj1]+*MVASelJets1[wj2];
//                TLorentzVector Bh = *MVASelJets1[bj1];
//                TLorentzVector Th = Wh+Bh;
//
//                double TriJetMass = Th.M();
//
//                double DiJetMass = Wh.M();
//                //DeltaR
//                float AngleThWh = fabs(Th.DeltaPhi(Wh));
//                float AngleThBh = fabs(Th.DeltaPhi(Bh));
//
//                float btag = MVASelJets1[bj1]->btag_combinedSecondaryVertexBJetTags();
//
//                double PtRat = (( *MVASelJets1[0] + *MVASelJets1[1] + *MVASelJets1[2] ).Pt())/( MVASelJets1[0]->Pt() + MVASelJets1[1]->Pt() + MVASelJets1[2]->Pt() );
//                if (debug) cout <<"Processing event with jetcombiner : 4 "<< endl;
//
//                MSPlot["MVA1TriJetMass"]->Fill(TriJetMass,  datasets[d], true, Luminosity*scaleFactor );
//                MSPlot["MVA1DiJetMass"]->Fill(DiJetMass,  datasets[d], true, Luminosity*scaleFactor );
//                MSPlot["MVA1BTag"]->Fill(btag,  datasets[d], true, Luminosity*scaleFactor );
//                MSPlot["MVA1PtRat"]->Fill(PtRat,  datasets[d], true, Luminosity*scaleFactor );
//                MSPlot["MVA1AnThWh"]->Fill(AngleThWh,  datasets[d], true, Luminosity*scaleFactor );
//                MSPlot["MVA1AnThBh"]->Fill(AngleThBh,  datasets[d], true, Luminosity*scaleFactor );
//
//
//                if (debug) cout <<"Processing event with jetcombiner : 8 "<< endl;
//
//
//            }

                ///////////////////////////////////
                // Filling histograms / plotting //
                ///////////////////////////////////

                MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);




                //////////////////////
                // Muon Based Plots //
                //////////////////////
                MSPlot["NbOfMuonsPostSel"]->Fill(nMu, datasets[d], true, Luminosity*scaleFactor);
                if(nMu>=2 && Muon && !Electron) {
                    pairCharge = selectedMuons[0]->charge() + selectedMuons[1]->charge();
                    MSPlot["MuonPairCharge"]->Fill(pairCharge, datasets[d], true, Luminosity*scaleFactor);
                    }

                for (Int_t selmu =0; selmu < selectedMuons.size(); selmu++ ) {
                    float reliso = selectedMuons[selmu]->relPfIso(4, 0.5);
                    double muzPVz = fabs(selectedMuons[selmu]->vz() - vertex[0]->Z());
                    totalCharge += selectedMuons[selmu]->charge();
                    if(debug) cout << "Reliso: "<< reliso << endl;
                    if(debug) cout << " chISO: " << selectedMuons[selmu]->chargedHadronIso() << endl;
                    MSPlot["MuonDz"]->Fill(selectedMuons[selmu]->dz() , datasets[d], true, Luminosity*scaleFactor );
                    MSPlot["MuonPt"]->Fill(selectedMuons[selmu]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["MuonPtRel"]->Fill(MuonPTRel[selmu], datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["MuonEta"]->Fill(selectedMuons[selmu]->Eta(), datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["MuonPhi"]->Fill(selectedMuons[selmu]->Phi(), datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["MuonNValidHits"]->Fill(selectedMuons[selmu]->nofValidHits(), datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["Muond0"]->Fill(selectedMuons[selmu]->d0(), datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["MuonDistVzPVz"]->Fill(muzPVz, datasets[d], true, Luminosity*scaleFactor );
                    MSPlot["MuonNMatchedStations"]->Fill(selectedMuons[selmu]->nofMatchedStations(), datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["MuonTrackerLayersWithMeasurement"]->Fill(selectedMuons[selmu]->nofTrackerLayersWithMeasurement(), datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["MuonDRMin"]->Fill(MuondRMin[selmu], datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["MuonRelIsolation"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
                    }

                if(Muon && !Electron) {
                    //Truth matching in Di-Muon channel
                    for(unsigned int i=0; i<mcParticlesMatching_.size(); i++) {
                        if( abs(mcParticlesMatching_[i]->type()) == 13 && mcParticlesMatching_[i]->status() == 1) { //Final State Muon
                            mcMuonsTLV.push_back(*mcParticlesMatching_[i]);
                            mcMuonIndex.push_back(i);
                            }
                        }
                    muonMatching = JetPartonMatching(selectedMuonsTLV_JC, mcMuonsTLV, 2, true, true, 0.3);

                    for(unsigned int i=0; i<selectedMuonsTLV_JC.size(); i++) {
                        int matchedMuonNumber = muonMatching.getMatchForParton(i, 0); //Gives the index of mcMuonsTLV where the match is
                        if(matchedMuonNumber < 0) MSPlot["LeptonTruth"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
                        else MSPlot["LeptonTruth"]->Fill(1, datasets[d], true, Luminosity*scaleFactor);
                        if(abs(pairCharge!=0) && i<2 ) {
                            if(matchedMuonNumber >=0) {
                                MSPlot["SameSignMuonType"]->Fill(abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->type()), datasets[d], true, Luminosity*scaleFactor);
                                if(abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->type()) == 13) {
                                    int motherType = abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->motherType());
                                    if(motherType<25) MSPlot["SameSignMuonParentType"]->Fill(motherType, datasets[d], true, Luminosity*scaleFactor);
                                    else MSPlot["SameSignMuonParentType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                                    int grannyType = abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->grannyType());
                                    if(grannyType<25) MSPlot["SameSignMuonGrannyType"]->Fill(grannyType, datasets[d], true, Luminosity*scaleFactor);
                                    else MSPlot["SameSignMuonGrannyType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                                    }
                                }
                            else {
                                MSPlot["SameSignMuonType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                                }
                            }
                        if(abs(pairCharge==0) && i<2 ) {
                            if(matchedMuonNumber >=0) {
                                MSPlot["OpSignMuonType"]->Fill(abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->type()), datasets[d], true, Luminosity*scaleFactor);
                                if(abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->type()) == 13) {
                                    int motherType = abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->motherType());
                                    MSPlot["OpSignMuonParentType"]->Fill(motherType, datasets[d], true, Luminosity*scaleFactor);
                                    int grannyType = abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->grannyType());
                                    MSPlot["OpSignMuonGrannyType"]->Fill(grannyType, datasets[d], true, Luminosity*scaleFactor);
                                    }
                                else {
                                    MSPlot["OpSignMuonType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                                    MSPlot["OpSignMuonParentType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                                    MSPlot["OpSignMuonGrannyType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                                    }
                                }
                            else {
                                MSPlot["OpSignMuonType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                                MSPlot["OpSignMuonParentType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                                MSPlot["OpSignMuonGrannyType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                                }
                            }
                        }
                    }


                if(pairCharge!=0 && Muon && !Electron) {
                    MSPlot["SameSignTotalCharge"]->Fill(totalCharge, datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["SameSignNbOfMuons"]->Fill(nMu, datasets[d], true, Luminosity*scaleFactor);
                    selecTable.Fill(d,7,1);
                    if(nMu>2) selecTable.Fill(d,8,1);
                    }
//            if(diMuMass < 10.0) {
//                MSPlot["BadDiMuon_InvMass"]->Fill(diMuMass, datasets[d], true, Luminosity*scaleFactor);
//                }
                MSPlot["DiMuon_InvMass"]->Fill(diMuMass, datasets[d], true, Luminosity*scaleFactor);


                //////////////////////////
                // Electron Based Plots //
                //////////////////////////

                MSPlot["NbOfIsolatedElectrons"]->Fill(nEl, datasets[d], true, Luminosity*scaleFactor);
                for (Int_t selel =0; selel < selectedElectrons.size(); selel++ ) {
                    float reliso = (selectedElectrons[selel]->chargedHadronIso() + max( 0.0, selectedElectrons[selel]->neutralHadronIso() + selectedElectrons[selel]->photonIso() - 0.5*selectedElectrons[selel]->puChargedHadronIso() ) ) / selectedElectrons[selel]->Pt();
                    double elzPVz = fabs(selectedElectrons[selel]->vz() - vertex[0]->Z());
                    selectedElectronsTLV_JC.push_back(*selectedElectrons[selel]);
                    MSPlot["ElectronDz"]->Fill(selectedElectrons[selel]->dz() , datasets[d], true, Luminosity*scaleFactor );
                    MSPlot["ElectronPt"]->Fill(selectedElectrons[selel]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["ElectronEta"]->Fill(selectedElectrons[selel]->Eta(), datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["ElectronPhi"]->Fill(selectedElectrons[selel]->Phi(), datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["Electrond0"]->Fill(selectedElectrons[selel]->d0(), datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["ElectronDistVzPVz"]->Fill(elzPVz, datasets[d], true, Luminosity*scaleFactor );
                    MSPlot["ElectronRelIsolation"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
                    }
                if(nEl>=2 && !Muon && Electron) {
                    float pairCharge = selectedElectrons[0]->charge() + selectedElectrons[1]->charge();
                    MSPlot["MuonPairCharge"]->Fill(pairCharge, datasets[d], true, Luminosity*scaleFactor);
                    TLorentzVector diElectron = selectedElectronsTLV_JC[0] + selectedElectronsTLV_JC[1];
                    MSPlot["DiElectron_InvMass"]->Fill(diElectron.M(), datasets[d], true, Luminosity*scaleFactor);
                    }


                ////////////////////////////
                // Plots for MuEl Channel //
                ////////////////////////////

                if(nEl>=1 && nMu>=1 && Muon && Electron) {
                    float pairCharge = selectedMuons[0]->charge() + selectedElectrons[0]->charge();
                    MSPlot["MuElPairCharge"]->Fill(pairCharge, datasets[d], true, Luminosity*scaleFactor);
                    }

                //////////////////////
                // Jets Based Plots //
                //////////////////////

                MSPlot["NbOfSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["NbOfSelectedLightJets"]->Fill(selectedLightJets.size(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["NbOfSelectedMBJets"]->Fill(nMtags, datasets[d], true, Luminosity*scaleFactor);
                if (selectedJets.size()>=3) MSPlot["3rdJetPt"]->Fill(selectedJets[2]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                if (selectedJets.size()>=4) MSPlot["4thJetPt"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                if (selectedJets.size()>=5) MSPlot["5thJetPt"]->Fill(selectedJets[4]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                if (selectedJets.size()>=6) MSPlot["6thJetPt"]->Fill(selectedJets[5]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                if (selectedJets.size()>=7) MSPlot["7thJetPt"]->Fill(selectedJets[6]->Pt(), datasets[d], true, Luminosity*scaleFactor);


                ////////////////////////////////////
                // Plotting event-level variables //
                ////////////////////////////////////

                HT = 0;
                H = 0;
                double HT1M2L=0, H1M2L=0;

                for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ ) {
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
                    if(seljet1>=3) {
                        HT1M2L += jetpt;
                        H1M2L += selectedJets[seljet1]->P();
                        }
                    }
                MSPlot["HTExcess1M2L"]->Fill(HT1M2L, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["HExcess1M2L"]->Fill(H1M2L, datasets[d], true, Luminosity*scaleFactor);
                HTH = HT/H;
                MSPlot["H"]->Fill(H, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["HTH"]->Fill(HTH, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["HT_SelectedJets"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["HTRat"]->Fill(HTRat, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MET"]->Fill(mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor);
                } //End Loop on Events
            cout <<"n events passed  =  "<<passed <<endl;
            //important: free memory
            treeLoader.UnLoadDataset();
            } //End Loop on Datasets

        eventlist.close();

        /////////////
        // Writing //
        /////////////

        cout << " - Writing outputs to the files ..." << endl;

        //////////////////////
        // Selection tables //
        //////////////////////

        //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
        selecTable.TableCalculator(  true, true, true, true, true);

        //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawsyNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
        selecTable.Write(  "FourTop"+postfix+"_Table"+channelpostfix+".tex",    false,true,true,true,false,false,false);

        fout->cd();
        TFile *foutmva = new TFile ("foutMVA.root","RECREATE");
        cout <<" after cd .."<<endl;

        string pathPNGJetCombi = pathPNG + "JetCombination/";
        mkdir(pathPNGJetCombi.c_str(),0777);
        if(TrainMVA)jetCombiner->Write(foutmva, true, pathPNGJetCombi.c_str());

//Output ROOT file
        for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin();
                it != MSPlot.end();
                it++) {
            string name = it->first;
            MultiSamplePlot *temp = it->second;
            temp->Write(fout, name, true, pathPNG, "pdf");
            }
        delete fout;

        cout << "It took us " << ((double)clock() - start) / CLOCKS_PER_SEC << " to run the program" << endl;
        cout << "********************************************" << endl;
        cout << "           End of the program !!            " << endl;
        cout << "********************************************" << endl;

        return 0;
        }

    int Factorial(int N = 1) {
        int fact = 1;
        for( int i=1; i<=N; i++ )
            fact = fact * i;  // OR fact *= i;
        return fact;
        }
























