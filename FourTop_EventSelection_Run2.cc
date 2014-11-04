//////////////////////////////////////////////////////////////////////////////
////         Analysis code for search for Four Top Production.                  ////
////////////////////////////////////////////////////////////////////////////////////

// ttbar @ NLO 13 TeV:
//all-had ->679 * .46 = 312.34
//semi-lep ->679 *.45 = 305.55
//di-lep-> 679* .09 = 61.11

//ttbar @ NNLO 8 TeV:
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

struct HighestTCHEBtag
{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const
    {
        return j1->btag_trackCountingHighEffBJetTags() > j2->btag_trackCountingHighEffBJetTags();
    }
};
struct HighestCVSBtag
{
    bool operator()( TRootJet* j1, TRootJet* j2 ) const
    {
        return j1->btag_combinedSecondaryVertexBJetTags() > j2->btag_combinedInclusiveSecondaryVertexV2BJetTags();
    }
};

bool match;

//To cout the Px, Py, Pz, E and Pt of objects
int Factorial(int N);

int main (int argc, char *argv[])
{

    ofstream eventlist;
    eventlist.open ("interesting_events_mu.txt");

    int passed = 0;
    int ndefs =0;

    string btagger = "CSVM";
    float scalefactorbtageff, mistagfactor;
    float workingpointvalue = 0.679; //working points updated to 2012 BTV-POG recommendations.
    bool bx25 = false;

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


    string postfix = "_Run2_TopTree_Study"; // to relabel the names of the output file

    ///////////////////////////////////////
    // Configuration
    ///////////////////////////////////////

    string channelpostfix = "";
    string xmlFileName = "";

    //Setting Lepton Channels (Setting both flags true will select Muon-Electron Channel when dilepton is also true)
    bool dilepton = true;
    bool Muon = true;
    bool Electron = true;

    if(Muon && !Electron && dilepton)
    {
        cout << " --> Using the di-Muon channel..." << endl;
        channelpostfix = "_MuMu";
        xmlFileName = "config/Run2_Samples.xml";
    }
    else if(Muon && Electron && dilepton)
    {
        cout << " --> Using the Muon-Electron channel..." << endl;
        channelpostfix = "_MuEl";
        xmlFileName = "config/Run2_Samples.xml";
    }
    else if(!Muon && Electron && dilepton)
    {
        cout << " --> Using the di-Electron channel..." << endl;
        channelpostfix = "_ElEl";
        xmlFileName = "config/Run2_Samples.xml";
    }
    else if(Muon && !Electron && !dilepton)
    {
        cout << " --> Using the Single Muon channel..." << endl;
        channelpostfix = "_Mu";
        xmlFileName = "config/Run2_Samples.xml";
//        xmlFileName = "config/test_mconly.xml";
    }
    else if(!Muon && Electron && !dilepton)
    {
        cout << " --> Using the Single Electron channel..." << endl;
        channelpostfix = "_El";
        xmlFileName = "config/Run2_Samples.xml";
//        xmlFileName = "config/test_mconly.xml";
    }
    else
    {
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
    float Luminosity = 5000.0; //pb^-1??
    vector<string> MVAvars;

    //A few bools to steer the MassReco and Event MVAs
    string MVAmethod = "BDT"; // MVAmethod to be used to get the good jet combi calculation (not for training! this is chosen in the jetcombiner class)

    cout <<"Instantiating jet combiner..."<<endl;

//    JetCombiner* jetCombiner = new JetCombiner(TrainMVA, Luminosity, datasets, MVAmethod, false);
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

    for (unsigned int d = 0; d < datasets.size (); d++)
    {
        cout <<"found sample with equivalent lumi "<<  datasets[d]->EquivalentLumi() <<endl;
        dataSetName = datasets[d]->Name();
        if(dataSetName.find("Data")<=0 || dataSetName.find("data")<=0 || dataSetName.find("DATA")<=0)
        {
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
//    for (unsigned int d = 0; d < datasets.size (); d++)
//    {
//        if( datasets[d]->Name()=="TTJets"  ||datasets[d]->Name()=="TTJets_AllHad" || datasets[d]->Name()=="TTJets_Other")
//        {
//            currentLumi = datasets[d]->EquivalentLumi();
//            cout <<"Old lumi =   "<< currentLumi  <<endl;
//            newlumi = currentLumi/frac;
//            datasets[d]->SetEquivalentLuminosity(newlumi);
//        }
//    }


    // for splitting the ttbar sample, it is essential to have the ttjets sample as the last
    //dataset loaded

//    cout << " - splitting TTBar dataset ..." << ndatasets   << endl;
//    vector<string> ttbar_filenames = datasets[ndatasets]->Filenames();
//    cout <<"ttbar filenames =  "<< ttbar_filenames[0] <<endl;
//
//    Dataset* ttbar_loose = new Dataset("TTJets_LooseSel","tt loose selection" , true, 633, 2, 2, 1, 213.4,ttbar_filenames );
//
//    ttbar_loose->SetEquivalentLuminosity(newlumi);
//
//    ttbar_loose->SetColor(kRed);

    //datasets.pop_back();
//    datasets.push_back(ttbar_loose);

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
    MSPlot["NbOfVertices"]                                  = new MultiSamplePlot(datasets, "NbOfVertices", 60, 0, 60, "Nb. of vertices");
    MSPlot["NbOfVerticesPreSel"]                            = new MultiSamplePlot(datasets, "NbOfVerticesPreSel", 60, 0, 60, "Nb. of vertices");
    //Muons
    MSPlot["NbOfGenMuons"]                                  = new MultiSamplePlot(datasets, "NbOfGenMuons", 5, 0, 5, "Nb. of gen muons");
    MSPlot["NbOfMuonsPreSel"]                               = new MultiSamplePlot(datasets, "NbOfMuonsPreSel", 5, 0, 5, "Nb. of muons");
    MSPlot["NbOfMuonsPostSel"]                              = new MultiSamplePlot(datasets, "NbOfMuonsPostSel", 5, 0, 5, "Nb. of muons");
    MSPlot["MuonDRMin"]                                     = new MultiSamplePlot(datasets, "MuonDRMin", 40, 0, 4, "dR_{min}");
    MSPlot["IsolatedMuonDRMin"]                             = new MultiSamplePlot(datasets, "IsolatedMuonDRMin", 40, 0, 4, "dR_{min}");
    MSPlot["MuonRelIsolation"]                              = new MultiSamplePlot(datasets, "MuonRelIsolation", 10, 0, .25, "RelIso");
    MSPlot["MuonPt"]                                        = new MultiSamplePlot(datasets, "MuonPt", 30, 0, 300, "PT_{#mu}");
    MSPlot["MuonPtRel"]                                     = new MultiSamplePlot(datasets, "MuonPtRel", 30, 0, 300, "PT_{#mu}");
    MSPlot["MuonEta"]                                       = new MultiSamplePlot(datasets, "MuonEta", 24, -2.4, 2.4, "#eta_{#mu}");
    MSPlot["MuonPhi"]                                       = new MultiSamplePlot(datasets, "MuonPhi", 35, -3.5, 3.5, "#phi_{#mu}");
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
    MSPlot["MuonChi2"]                                      = new MultiSamplePlot(datasets, "MuonChi2", 11, 0, 11, "#chi^{2}/dof");
    MSPlot["MuonNPixelHits"]                                = new MultiSamplePlot(datasets, "MuonNPixelHits", 6, 0, 6, "NPixelHits_{#mu}");
    MSPlot["MuonRelIsolationPreSel"]                        = new MultiSamplePlot(datasets, "MuonRelIsolationPreSel", 10, 0, .25, "RelIso");
    MSPlot["MuonPtPreSel"]                                  = new MultiSamplePlot(datasets, "MuonPtPreSel", 30, 0, 300, "PT_{#mu}");
    MSPlot["MuonEtaPreSel"]                                 = new MultiSamplePlot(datasets, "MuonEtaPreSel", 24, -2.4, 2.4, "#eta_{#mu}");
    MSPlot["MuonRelNHI4PreSel"]                             = new MultiSamplePlot(datasets, "MuonRelNHI4PreSel", 50, 0,1, "Relative NHad p_{T}");
    MSPlot["MuonRelCHI4PreSel"]                             = new MultiSamplePlot(datasets, "MuonRelCHI4PreSel", 50, 0,1, "Relative CHad p_{T}");
    MSPlot["MuonRelPhI4PreSel"]                             = new MultiSamplePlot(datasets, "MuonRelPhI4PreSel", 50, 0,1, "Relative Photon p_{T}");
    MSPlot["MuonRelPUCHI4PreSel"]                           = new MultiSamplePlot(datasets, "MuonRelPUCHI4PreSel", 50, 0,1, "Relative CHad p_{T} from PU");
    MSPlot["MuonRelTotN4PreSel"]                            = new MultiSamplePlot(datasets, "MuonRelTotN4PreSel", 50, 0,1, "Relative Neu p_{T}");
    MSPlot["MuonRelIsoPreSel"]                              = new MultiSamplePlot(datasets, "MuonRelIsoPreSel", 50, 0, 1, "RelIso");
    MSPlot["MuonRelNminusPUCH"]                             = new MultiSamplePlot(datasets, "MuonRelNminusPUCH", 50, -1, 1, "Relative N-PUCH");
    MSPlot["MuonRelNminusBPUCH"]                            = new MultiSamplePlot(datasets, "MuonRelNminusBPUCH", 50, -1, 1, "Relative N-#beta * PUCH");
    MSPlot["MuonBeta"]                                      = new MultiSamplePlot(datasets, "MuonBeta", 50, -1, 1, "#Delta#beta");

    //Different dBeta RelIso Plots
    MSPlot["MuonRelIso0"]                                   = new MultiSamplePlot(datasets, "MuonRelIso0", 50, 0, 1, "RelIso(#Delta#beta = 0.1)");
    MSPlot["MuonRelIso1"]                                   = new MultiSamplePlot(datasets, "MuonRelIso1", 50, 0, 1, "RelIso(#Delta#beta = 0.1)");
    MSPlot["MuonRelIso2"]                                   = new MultiSamplePlot(datasets, "MuonRelIso2", 50, 0, 1, "RelIso(#Delta#beta = 0.2)");
    MSPlot["MuonRelIso3"]                                   = new MultiSamplePlot(datasets, "MuonRelIso3", 50, 0, 1, "RelIso(#Delta#beta = 0.3)");
    MSPlot["MuonRelIso4"]                                   = new MultiSamplePlot(datasets, "MuonRelIso4", 50, 0, 1, "RelIso(#Delta#beta = 0.4)");
    MSPlot["MuonRelIso5"]                                   = new MultiSamplePlot(datasets, "MuonRelIso5", 50, 0, 1, "RelIso(#Delta#beta = 0.5)");
    MSPlot["MuonRelIso6"]                                   = new MultiSamplePlot(datasets, "MuonRelIso6", 50, 0, 1, "RelIso(#Delta#beta = 0.6)");
    MSPlot["MuonRelIso7"]                                   = new MultiSamplePlot(datasets, "MuonRelIso7", 50, 0, 1, "RelIso(#Delta#beta = 0.7)");
    MSPlot["MuonRelIso8"]                                   = new MultiSamplePlot(datasets, "MuonRelIso8", 50, 0, 1, "RelIso(#Delta#beta = 0.8)");
    MSPlot["MuonRelIso9"]                                   = new MultiSamplePlot(datasets, "MuonRelIso9", 50, 0, 1, "RelIso(#Delta#beta = 0.9)");
    MSPlot["MuonRelIso10"]                                  = new MultiSamplePlot(datasets, "MuonRelIso10", 50, 0, 1, "RelIso(#Delta#beta = 1.0)");



    //Electrons
    MSPlot["NbOfIsolatedElectrons"]                         = new MultiSamplePlot(datasets, "NbOfIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
    MSPlot["NbOfExtraIsolatedElectrons"]                    = new MultiSamplePlot(datasets, "NbOfExtraIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
    MSPlot["ElectronRelIsolation"]                          = new MultiSamplePlot(datasets, "ElectronRelIsolation", 10, 0, .25, "RelIso");
    MSPlot["ElectronPt"]                                    = new MultiSamplePlot(datasets, "ElectronPt", 30, 0, 300, "PT_{e}");
    MSPlot["ElectronEta"]                                   = new MultiSamplePlot(datasets, "ElectronEta", 24, -2.4, 2.4, "#eta_{e}");
    MSPlot["ElectronPhi"]                                   = new MultiSamplePlot(datasets, "ElectronPhi", 35, -3.5, 3.5, "#phi_{e}");
    MSPlot["Electrond0"]                                    = new MultiSamplePlot(datasets, "Electrond0", 50, -0.05, 0.05, "d0_{e}");
    MSPlot["ElectrondZPVz"]                                 = new MultiSamplePlot(datasets, "ElectrondZPVz", 50, 0, .5, "dZPVZ_{e}");
    MSPlot["ElectrondRJets"]                                = new MultiSamplePlot(datasets, "ElectrondRJets", 50, 0, 10, "dRJets_{e}");
    MSPlot["ElectronDistVzPVz"]                             = new MultiSamplePlot(datasets, "ElectronDistVzPVz", 50, 0 ,.3, "DistVzPVz_{e}");
    MSPlot["Electrondz"]                                    = new MultiSamplePlot(datasets, "Electrondz", 25, -.6 ,.6, "dZ_{e}");
    MSPlot["DiElectron_InvMass"]                            = new MultiSamplePlot(datasets, "DiElectron_InvMass", 60, 0, 120, "DiElectron_InvMass");
    MSPlot["ElectronPairCharge"]                            = new MultiSamplePlot(datasets, "ElectronPairCharge", 5, -2.5, 2.5, "Total charge");
    MSPlot["ElectronTrackChi2"]                             = new MultiSamplePlot(datasets, "ElectronTrackChi2", 20, 0, 20, "#chi^{2}/dof");
    MSPlot["ElectronGSFChi2"]                               = new MultiSamplePlot(datasets, "ElectronGSFChi2", 20, 0, 20, "#chi^{2}/dof");
    MSPlot["ElectronPtPreSel"]                              = new MultiSamplePlot(datasets, "ElectronPtPreSel", 30, 0, 300, "PT_{e}");
    MSPlot["ElectronEtaPreSel"]                             = new MultiSamplePlot(datasets, "ElectronEtaPreSel", 24, -2.4, 2.4, "#eta_{e}");
    MSPlot["ElectronRelIsolationPreSel"]                    = new MultiSamplePlot(datasets, "ElectronRelIsolationPreSel", 10, 0, .25, "RelIso");



    //Plots Specific to MuEl channel
    MSPlot["MuElPairCharge"]                                = new MultiSamplePlot(datasets, "MuElPairCharge", 5, -2.5, 2.5, "Total charge");

    //B-tagging discriminators
    MSPlot["BdiscBJetCand_CSV"]                             = new MultiSamplePlot(datasets, "BdiscBJetCand_CSV", 20, 0, 1, "CSV b-disc.");
    MSPlot["BdiscBJetCand_CSV_LowPU"]                       = new MultiSamplePlot(datasets, "BdiscBJetCand_CSV_LowPU", 20, 0, 1, "CSV b-disc.");
    MSPlot["BdiscBJetCand_CSV_HighPU"]                      = new MultiSamplePlot(datasets, "BdiscBJetCand_CSV_HighPU", 20, 0, 1, "CSV b-disc.");
    //Jets
    MSPlot["NbOfSelectedJets"]                              = new MultiSamplePlot(datasets, "NbOfSelectedJets", 20, 0, 20, "Nb. of jets");
    MSPlot["NbOfSelectedLightJets"]                         = new MultiSamplePlot(datasets, "NbOfSelectedLightJets", 10, 0, 10, "Nb. of jets");
    MSPlot["NbOfSelectedMBJets"]                            = new MultiSamplePlot(datasets, "NbOfSelectedMBJets", 8, 0, 8, "Nb. of CSVM tags");
    MSPlot["JetEta"]                                        = new MultiSamplePlot(datasets, "JetEta", 40,-4, 4, "Jet #eta");
    MSPlot["JetPhi"]                                        = new MultiSamplePlot(datasets, "JetPhi", 35, -3.5, 3.5, "Jet #phi");
    MSPlot["JetCEF"]                                        = new MultiSamplePlot(datasets, "JetCEF", 50, 0, 1, "Jet CHEF");
    MSPlot["JetNEF"]                                        = new MultiSamplePlot(datasets, "JetNEF", 50, 0, 1, "Jet NEEF");
    MSPlot["SelectedJetPt"]                                 = new MultiSamplePlot(datasets, "SelectedJetPt", 30, 0, 300, "PT_{jet}");
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
    MSPlot["HTRat"]                                         = new MultiSamplePlot(datasets, "HTRat", 40, 0, 20, "HTRat");
    MSPlot["HTExcess2M"]                                    = new MultiSamplePlot(datasets, "HTExcess2M", 50, 0, 1000, "HT");
    MSPlot["HExcess2M"]                                     = new MultiSamplePlot(datasets, "HExcess2M", 50, 0, 3000, "H");
    MSPlot["HTExcess1M2L"]                                  = new MultiSamplePlot(datasets, "HTExcess1M2L", 50, 0, 1000, "HT");
    MSPlot["HExcess1M2L"]                                   = new MultiSamplePlot(datasets, "HExcess1M2L", 50, 0, 3000, "H");
    MSPlot["SelectedJetPtPreSel"]                           = new MultiSamplePlot(datasets, "SelectedJetPtPreSel", 30, 0, 300, "PT_{jet}");
    MSPlot["JetEtaPreSel"]                                  = new MultiSamplePlot(datasets, "JetEtaPreSel", 40,-4, 4, "Jet #eta");
    MSPlot["JetEnergyRatio"]                                = new MultiSamplePlot(datasets, "JetEnergyRatio", 50,0, 10, "Charge/Neutral");
    MSPlot["JetMultiplicityRatio"]                          = new MultiSamplePlot(datasets, "JetMultiplicityRatio", 50,0, 10, "Charge/Neutral");



    //MET
    MSPlot["MET"]                                           = new MultiSamplePlot(datasets, "MET", 70, 0, 700, "MET");

    //Truth Plots
    MSPlot["LeptonTruth"]                                   = new MultiSamplePlot(datasets, "LeptonTruth", 6, -1, 5, "LeptonTruth");
    MSPlot["CSVLEfficiency"]                                = new MultiSamplePlot(datasets, "CSVLEfficiency", 2, 0, 2, "True B-Jet");
    MSPlot["CSVMEfficiency"]                                = new MultiSamplePlot(datasets, "CSVMEfficiency", 2, 0, 2, "True B-Jet");
    MSPlot["CSVTEfficiency"]                                = new MultiSamplePlot(datasets, "CSVTEfficiency", 2, 0, 2, "True B-Jet");
    MSPlot["TrueBJetPt"]                                    = new MultiSamplePlot(datasets, "TrueBJetPt", 30, 0, 300, "PT_{jet}");
    MSPlot["TrueBJetEta"]                                   = new MultiSamplePlot(datasets, "TrueBJetEta", 40, -4, 4, "#eta_{jet}");
    MSPlot["TrueBJetCSV"]                                   = new MultiSamplePlot(datasets, "TrueBJetCSV", 20, 0, 1, "CSV_{jet}");
    MSPlot["TrueCJetPt"]                                    = new MultiSamplePlot(datasets, "TrueCJetPt", 30, 0, 300, "PT_{jet}");
    MSPlot["TrueCJetEta"]                                   = new MultiSamplePlot(datasets, "TrueCJetEta", 40, -4, 4, "#eta_{jet}");
    MSPlot["TrueCJetCSV"]                                   = new MultiSamplePlot(datasets, "TrueCJetCSV", 20, 0, 1, "CSV_{jet}");
    MSPlot["TrueLJetPt"]                                    = new MultiSamplePlot(datasets, "TrueLJetPt", 30, 0, 300, "PT_{jet}");
    MSPlot["TrueLJetEta"]                                   = new MultiSamplePlot(datasets, "TrueLJetEta", 40, -4, 4, "#eta_{jet}");
    MSPlot["TrueLJetCSV"]                                   = new MultiSamplePlot(datasets, "TrueLJetCSV", 20, 0, 1, "CSV_{jet}");
    MSPlot["TrueBJetPtLowPU"]                               = new MultiSamplePlot(datasets, "TrueBJetPtLowPU", 30, 0, 300, "PT_{jet}");
    MSPlot["TrueBJetEtaLowPU"]                              = new MultiSamplePlot(datasets, "TrueBJetEtaLowPU", 40, -4, 4, "#eta_{jet}");
    MSPlot["TrueBJetCSVLowPU"]                              = new MultiSamplePlot(datasets, "TrueBJetCSVLowPU", 20, 0, 1, "CSV_{jet}");
    MSPlot["TrueBJetPtHighPU"]                              = new MultiSamplePlot(datasets, "TrueBJetPtHighPU", 30, 0, 300, "PT_{jet}");
    MSPlot["TrueBJetEtaHighPU"]                             = new MultiSamplePlot(datasets, "TrueBJetEtaHighPU", 40, -4, 4, "#eta_{jet}");
    MSPlot["TrueBJetCSVHighPU"]                             = new MultiSamplePlot(datasets, "TrueBJetCSVHighPU", 20, 0, 1, "CSV_{jet}");
    MSPlot["CSVMTrueBJetPt"]                                = new MultiSamplePlot(datasets, "CSVMTrueBJetPt", 30, 0, 300, "PT_{jet}");
    MSPlot["CSVMTrueCJetPt"]                                = new MultiSamplePlot(datasets, "CSVMTrueCJetPt", 30, 0, 300, "PT_{jet}");
    MSPlot["CSVMTrueLJetPt"]                                = new MultiSamplePlot(datasets, "CSVMTrueLJetPt", 30, 0, 300, "PT_{jet}");
    MSPlot["CSVMTrueBJetEta"]                               = new MultiSamplePlot(datasets, "CSVMTrueBJetEta", 40, -4, 4, "#eta_{jet}");
    MSPlot["CSVMTrueCJetEta"]                               = new MultiSamplePlot(datasets, "CSVMTrueCJetEta", 40, -4, 4, "#eta_{jet}");
    MSPlot["CSVMTrueLJetEta"]                               = new MultiSamplePlot(datasets, "CSVMTrueLJetEta", 40, -4, 4, "#eta_{jet}");
    MSPlot["MuonType"]                                      = new MultiSamplePlot(datasets, "MuonType", 26, -1, 25, "Muon Type Code");
    MSPlot["MuonParentType"]                                = new MultiSamplePlot(datasets, "MuonParentType", 26, -1, 25, "Parent Type Code");
    MSPlot["MuonGrannyType"]                                = new MultiSamplePlot(datasets, "MuonGrannyType", 26, -1, 25, "Granny Type Code");
    MSPlot["SameSignMuonType"]                              = new MultiSamplePlot(datasets, "SameSignMuonType", 26, -1, 25, "Muon Type Code");
    MSPlot["SameSignMuonParentType"]                        = new MultiSamplePlot(datasets, "SameSignMuonParentType", 26, -1, 25, "Parent Type Code");
    MSPlot["SameSignMuonGrannyType"]                        = new MultiSamplePlot(datasets, "SameSignMuonGrannyType", 26, -1, 25, "Granny Type Code");
    MSPlot["SameSignTotalCharge"]                           = new MultiSamplePlot(datasets, "SameSignTotalCharge", 7, -3.5, 3.5, "Total Charge in Muons");
    MSPlot["SameSignNbOfMuons"]                             = new MultiSamplePlot(datasets, "SameSignNbOfMuons", 6, 0, 6, "Number of Muons");
    MSPlot["OpSignMuonType"]                                = new MultiSamplePlot(datasets, "OpSignMuonType", 26, -1, 25, "Muon Type Code");
    MSPlot["OpSignMuonParentType"]                          = new MultiSamplePlot(datasets, "OpSignMuonParentType", 26, -1, 25, "Parent Type Code");
    MSPlot["OpSignMuonGrannyType"]                          = new MultiSamplePlot(datasets, "OpSignMuonGrannyType", 26, -1, 25, "Granny Type Code");

    //N-1 Plots
    MSPlot["ElectronLoosePt"]                               = new MultiSamplePlot(datasets, "ElectronLoosePt", 30, 0, 300, "PT_{e}");
    MSPlot["ElectronLooseEta"]                              = new MultiSamplePlot(datasets, "ElectronLooseEta", 50, -5, 5, "#eta_{e}");
    MSPlot["ElectronLooseIso"]                              = new MultiSamplePlot(datasets, "ElectronLooseIso", 40, 0, 1, "PFRelIso_{e}");

    MSPlot["MuonLoosePt"]                                   = new MultiSamplePlot(datasets, "MuonLoosePt", 30, 0, 300, "PT_{#mu}");
    MSPlot["MuonLooseEta"]                                  = new MultiSamplePlot(datasets, "MuonLooseEta", 50, -5, 5, "#eta_{#mu}");
    MSPlot["MuonLooseIso"]                                  = new MultiSamplePlot(datasets, "MuonLooseIso", 40, 0, 1, "PFRelIso_{#mu}");

    MSPlot["JetLoosePt"]                                    = new MultiSamplePlot(datasets, "JetLoosePt", 30, 0, 300, "PT_{e}");
    MSPlot["JetLooseEta"]                                   = new MultiSamplePlot(datasets, "JetLooseEta", 50, -5, 5, "#eta_{e}");

    //MVA Top Roconstruction Plots
    MSPlot["MVA1TriJet"]                                    = new MultiSamplePlot(datasets, "MVA1TriJet", 30, -1.0, 0.2, "MVA1TriJet");
    MSPlot["MVA1TriJetMass"]                                = new MultiSamplePlot(datasets, "MVA1TriJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA1DiJetMass"]                                 = new MultiSamplePlot(datasets, "MVA1DiJetMass", 75, 0, 500, "m_{bjj}");
    MSPlot["MVA1PtRat"]                                     = new MultiSamplePlot(datasets, "MVA1PtRat", 25, 0, 2, "P_{t}^{Rat}");
    MSPlot["MVA1BTag"]                                      = new MultiSamplePlot(datasets, "MVA1BTag", 35, 0, 1, "BTag");
    MSPlot["MVA1AnThBh"]                                    = new MultiSamplePlot(datasets, "MVA1AnThBh", 35, 0, 3.14, "AnThBh");
    MSPlot["MVA1AnThWh"]                                    = new MultiSamplePlot(datasets, "MVA1AnThWh", 35, 0, 3.14, "AnThWh");


    ///////////////////
    // 1D histograms
    ///////////////////

    //Plots
    string pathPNG = "FourTop"+postfix+channelpostfix;
    pathPNG += "_MSPlots/";
    //pathPNG = pathPNG +"/";
    mkdir(pathPNG.c_str(),0777);

    cout <<"Making directory :"<< pathPNG  <<endl;
    vector<string> CutsselecTable;
    if(dilepton)
    {
        /////////////////////////////////
        // Selection table: Dilepton + jets
        /////////////////////////////////
        if(Muon && !Electron)
        {
            CutsselecTable.push_back(string("initial"));
            CutsselecTable.push_back(string("Event cleaning and Trigger"));
            CutsselecTable.push_back(string("At least 2 op sign Muons"));
            CutsselecTable.push_back(string("Muons Isolated (PFRelIso04$\\leq 0.2$)"));
            CutsselecTable.push_back(string("Highest 2 pT Muon Mass $\\geq 20 GeV$"));
            CutsselecTable.push_back(string("Z Veto ($m_{Z}\\pm 15 GeV$)"));
            CutsselecTable.push_back(string("At least 2 Jets"));
            CutsselecTable.push_back(string("MET $\\geq 40 GeV$"));
            CutsselecTable.push_back(string("At least 1 CSVL Jet"));
            CutsselecTable.push_back(string("At least 2 CSVL Jets"));
        }
        if(Muon && Electron)
        {
            CutsselecTable.push_back(string("initial"));
            CutsselecTable.push_back(string("Event cleaning and Trigger"));
            CutsselecTable.push_back(string("At least 1 Loose Isolated Muon"));
            CutsselecTable.push_back(string("At least 1 Loose Electron"));
            CutsselecTable.push_back(string("At least 4 Jets"));
            CutsselecTable.push_back(string("At least 1 CSVM Jet"));
            CutsselecTable.push_back(string("At least 2 CSVM Jets"));
//            CutsselecTable.push_back(string("Same Sign Leading Leptons"));
        }
        if(!Muon && Electron)
        {
            CutsselecTable.push_back(string("initial"));
            CutsselecTable.push_back(string("Event cleaning and Trigger"));
            CutsselecTable.push_back(string("Exactly 1 Isolated Electron (PFRelIso03$\\leq 0.15$)"));
//        CutsselecTable.push_back(string("Exactly 2 Loose Electron"));
            CutsselecTable.push_back(string("At least 2 CSVM Jets"));
            CutsselecTable.push_back(string("At least 2 Jets"));
            CutsselecTable.push_back(string("At least 3 Jets"));
            CutsselecTable.push_back(string("At least 4 Jets"));
        }
    }
    else
    {
        /////////////////////////////////
        // Selection table: Lepton + jets
        /////////////////////////////////
        if(Muon && !Electron)
        {
            CutsselecTable.push_back(string("initial"));
            CutsselecTable.push_back(string("Event cleaning and Trigger"));
            CutsselecTable.push_back(string("Exactly 1 Isolated Muon (relIso04$\\leq 0.12$)"));
            CutsselecTable.push_back(string("Loose Muon Veto"));
            CutsselecTable.push_back(string("Electron Veto"));
            CutsselecTable.push_back(string("At least 1 Jet"));
            CutsselecTable.push_back(string("At least 2 Jets"));
            CutsselecTable.push_back(string("At least 3 Jets"));
            CutsselecTable.push_back(string("At least 4 Jets"));
            CutsselecTable.push_back(string("At least 1 CSVM Jet"));
        }
        if(!Muon && Electron)
        {
            CutsselecTable.push_back(string("initial"));
            CutsselecTable.push_back(string("Event cleaning and Trigger"));
            CutsselecTable.push_back(string("Exactly 1 Isolated Electron"));
            CutsselecTable.push_back(string("Loose Electron Veto"));
            CutsselecTable.push_back(string("Loose Muon Veto"));
            CutsselecTable.push_back(string("At least 1 Jet"));
            CutsselecTable.push_back(string("At least 2 Jets"));
            CutsselecTable.push_back(string("At least 3 Jets"));
            CutsselecTable.push_back(string("At least 4 Jets"));
            CutsselecTable.push_back(string("At least 1 CSVM Jet"));
        }
    }

    SelectionTable selecTable(CutsselecTable, datasets);
    selecTable.SetLuminosity(Luminosity);
    selecTable.SetPrecision(1);

    /////////////////////////////////
    // Loop on datasets
    /////////////////////////////////

    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;

    for (unsigned int d = 0; d < datasets.size(); d++)
    {
        cout<<"Load Dataset"<<endl;
        treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset
        string previousFilename = "";
        int iFile = -1;
        dataSetName = datasets[d]->Name();
        if(datasets[d]->Name()=="TTJetsPU20bx25") bx25 = true;
        if(datasets[d]->Name()=="TTJetsPU40bx50") bx25 = false;
        if(bx25) cout << "Dataset with 25ns Bunch Spacing!" <<endl;
        else cout << "Dataset with 50ns Bunch Spacing!" <<endl;

        //////////////////////////////////////////////////
        // Initialize JEC factors ///////////////////////
        //////////////////////////////////////////////////

        vector<JetCorrectorParameters> vCorrParam;

        if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )   // Data!
        {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L1FastJet_AK5PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L2Relative_AK5PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L3Absolute_AK5PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L2L3Residual_AK5PFchs.txt");
            vCorrParam.push_back(*L2L3ResJetCorPar);
        }
        else
        {
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
        vector<TRootPFJet*>    selectedJets, selectedLoosePtJets, selectedLooseEtaJets;
        vector<TRootMuon*>     selectedMuons;
        vector<TRootMuon*>     selectedLooseIsoMuons, selectedLoosePtMuons, selectedLooseEtaMuons;
        vector<TRootElectron*> selectedLooseIsoElectrons, selectedLoosePtElectrons, selectedLooseEtaElectrons;
        vector<TRootElectron*> selectedExtraElectrons;
        vector<TRootMuon*>     selectedMuons_NoIso;
        vector<TRootMuon*>     selectedExtraMuons;
        selectedElectrons.reserve(10);
        selectedMuons.reserve(10);

        //////////////////////////////////////
        // Begin Event Loop
        //////////////////////////////////////

        for (unsigned int ievt = event_start; ievt < end_d; ievt++)
        {
            MHT = 0.,MHTSig = 0., STJet = 0., EventMass =0., EventMassX =0., SumJetMass = 0., SumJetMassX=0.  ,H = 0., HX =0., HT = 0., HTX = 0.,HTH=0.,HTXHX=0., sumpx_X = 0., sumpy_X= 0., sumpz_X =0., sume_X= 0. , sumpx =0., sumpy=0., sumpz=0., sume=0., jetpt =0., PTBalTopEventX = 0., PTBalTopSumJetX =0.;

            double ievt_d = ievt;
            currentfrac = ievt_d/end_d;
            if (debug)cout <<"event loop 1"<<endl;

            if(ievt%1000 == 0)
            {
                std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r"<<endl;
            }

            float scaleFactor = 1.;  // scale factor for the event
            event = treeLoader.LoadEvent (ievt, vertex, init_muons, init_electrons, init_jets, mets, debug);  //load event
            if (debug)cout <<"Number of Electrons Loaded: " << init_electrons.size() <<endl;

            float rho = event->fixedGridRhoFastjetAll();
            string graphName;

            /////////////////////////
            // Loading Gen jets
            /////////////////////////

            vector<TRootGenJet*> genjets;
//            if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) ) {
//                // loading GenJets as I need them for JER
//                genjets = treeLoader.LoadGenJet(ievt);
//                }
            // check which file in the dataset it is to have the HLTInfo right
            string currentFilename = datasets[d]->eventTree()->GetFile()->GetName();
            if(previousFilename != currentFilename)
            {
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
            if(previousRun != currentRun)
            {
                // cout <<"What run? "<< currentRun<<endl;
                previousRun = currentRun;
                if(Muon && !Electron)   // di-muon
                {
                    if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 )
                    {
                        if( event->runId() >= 190456 && event->runId() <= 190738 )
                        {
                            itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v11"), currentRun, iFile);
                            filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10";
                        }
                        else if( event->runId() >= 190782 && event->runId() <= 193621)
                        {
                            itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v12"), currentRun, iFile);
                            filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10";
                        }
                        else if(event->runId() >= 193834  && event->runId() <= 196531 )
                        {
                            itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
                            filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";
                        }
                        else if( event->runId() >= 198022  && event->runId() <= 199608)
                        {
                            itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v14"), currentRun, iFile);
                            filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";
                        }
                        else if( event->runId() >= 199698 && event->runId() <= 209151)
                        {
                            itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v15"), currentRun, iFile);
                            filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";
                        }
                        else
                        {
                            cout << "Unknown run for SemiMu HLTpath selection: " << event->runId() << endl;
                            filterName = "hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15";

                        }
                        if( itrigger == 9999 )
                        {
                            cout << "itriggerSemiMu == 9999 for SemiMu HLTpath selection: " << event->runId() << endl;
                            //exit(-1);
                        }
                    }
                    else if(Muon && Electron)   // Muon-Electron
                    {
                    }
                    else if(!Muon && Electron)   // di-Electron
                    {
                    }
                    else
                    {
                        if(dataSetName == "TTJets" ) itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);

                        else
                        {

                            itrigger = treeLoader.iTrigger (string ("HLT_IsoMu24_eta2p1_v13"), currentRun, iFile);
                        }

                        if(itrigger == 9999)
                        {
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

            if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) )
            {
                //JER
                if (debug) cout << "Doing JER Corrections: ";
                doJERShift = -1;
                if(doJERShift == 1)
                    jetTools->correctJetJER(init_jets, genjets, mets[0], "minus");
                else if(doJERShift == 2)
                    jetTools->correctJetJER(init_jets, genjets, mets[0], "plus");
                else if(doJERShift == 0)
                {
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
            if (debug)
            {
                cout<<"Debugging Jet Momentum"<<endl;
                for(int curJet = 0; curJet<init_jets.size(); curJet++)
                {
                    cout << "Px: "<<init_jets[curJet]->Px()<<" Py: "<<init_jets[curJet]->Py()<<" Pz: "<<init_jets[curJet]->Pz()<<" PT: "<<init_jets[curJet]->Pt()<<endl;
                    cout<<"Eta: "<<init_jets[curJet]->Eta()<<endl;
                }
                cout<<"Debugging Gen Jet Momentum"<<endl;
                for(int curJet = 0; curJet<genjets.size(); curJet++)
                {
                    cout << "Px: "<<genjets[curJet]->Px()<<" Py: "<<genjets[curJet]->Py()<<" Pz: "<<genjets[curJet]->Pz()<<" E: "<<genjets[curJet]->E()<<" PT: "<<genjets[curJet]->Pt()<<endl;
                    cout<<"Eta: "<<genjets[curJet]->Eta()<<endl;
                }
                cin.get();
            }



            ////////////////////////////////////////
            // Beam scraping and PU reweighting
            // Turned off on Aug 11, 2014 by JH for CSA14 studies
            ////////////////////////////////////////

            // scale factor for the event
            //float scaleFactor = 1.;

//            if(dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA") {
//                // Apply the scraping veto. (Is it still needed?)
//                bool isBeamBG = true;
//                if(event->nTracks() > 10) {
//                    if( ( (float) event->nHighPurityTracks() ) / ( (float) event->nTracks() ) > 0.25 )
//                        isBeamBG = false;
//                    }
//                if(isBeamBG) continue;
//                }
//            else {
//                double lumiWeight = LumiWeights.ITweight( (int)event->nTruePU() );
//                double lumiWeightOLD=lumiWeight;
//                if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0) {
//                    lumiWeight=1;
//                    }
//                scaleFactor = scaleFactor*lumiWeight;
//
//                }

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
            if(debug)
            {
                for(int selel = 0; selel < init_electrons.size(); selel++)
                {
                    cout << "PreSel Electron Pt: " << init_electrons[selel]->Pt() << endl;
                    cout << "PreSel Electron Eta: " << init_electrons[selel]->Eta() << endl;
                    cout << "PreSel Electron d0: " << init_electrons[selel]->d0() << endl;
                    cout << "PreSel Electron reliso: " << init_electrons[selel]->relPfIso(3, 0.5) << endl;
                    cout << "PreSel Electron MVA: " << init_electrons[selel]->mvaTrigId() << endl;
                    cout << "PreSel Electron dz: " << init_electrons[selel]->dz() << endl;
                    cout << "PreSel Electron MissHits: " << init_electrons[selel]->missingHits() << endl;
                    cout << "PreSel Electron Passed Conversion: " << init_electrons[selel]->passConversion() << endl;
                }
            }
            if(debug)
            {
                for(int selmu = 0; selmu < init_muons.size(); selmu++)
                {
                    cout << "PreSel Muon Pt: " << init_muons[selmu]->Pt() << endl;
                    cout << "PreSel Muon Eta: " << init_muons[selmu]->Eta() << endl;
                    cout << "PreSel Muon d0: " << init_muons[selmu]->d0() << endl;
                    cout << "PreSel Muon reliso: " << init_muons[selmu]->relPfIso(3, 0.5) << endl;
                    cout << "PreSel Muon Trasker Layers: " << init_muons[selmu]->nofTrackerLayersWithMeasurement() << endl;
                    cout << "PreSel Muon dz: " << init_muons[selmu]->dz() << endl;
                    cout << "PreSel Muon MatchedStations: " << init_muons[selmu]->nofMatchedStations() << endl;
                    cout << "PreSel Muon Pixel hits: " << init_muons[selmu]->nofValidPixelHits() << endl;
                }
            }
            // Define object selection cuts
            if(Muon && !Electron)   //Muon Channel Cuts
            {
                if(dilepton)   //Dilepton Cuts
                {
                    selection.setJetCuts(30.,2.5,0.01,1.,0.98,0,0); //Pt, Eta, EMF, n90Hits, fHPD, dRJetElectron, DRJets
                    selection.setElectronCuts(20,2.5,.15,0.04,.5,1,0); //Pt,Eta,RelIso,d0,MVAId,DistVzPVz, DRJets, MaxMissingHits
                    selection.setMuonCuts(20, 2.4, 0.20, 1, 0.3, 1, 1, 0, 0);

                    if (debug)cout<<"Getting Jets"<<endl;
                    selectedJets                                        = selection.GetSelectedPFJets(true); // ApplyJetId
                    selectedLoosePtJets                                 = selection.GetSelectedPFJets(0, 2.5, true);
                    selectedLooseEtaJets                                = selection.GetSelectedPFJets(30, 5, true);
                    if (debug)cout<<"Getting Tight Muons"<<endl;
                    selectedMuons                                       = selection.GetSelectedMuons();
                    selectedLoosePtMuons                                = selection.GetSelectedMuons(0, 2.4, 0.2);
                    selectedLooseEtaMuons                               = selection.GetSelectedMuons(20, 5, 0.2);
                    selectedLooseIsoMuons                               = selection.GetSelectedMuons(20, 2.4, 1);
                    if (debug)cout<<"Getting Loose Electrons"<<endl;
                    selectedElectrons                                   = selection.GetSelectedElectrons(20,2.5,0.15); // VBTF ID
                    selectedLoosePtElectrons                            = selection.GetSelectedElectrons(0, 2.5, 0.15);
                    selectedLooseEtaElectrons                           = selection.GetSelectedElectrons(20, 5, 0.15);
                    selectedLooseIsoElectrons                           = selection.GetSelectedElectrons(20, 2.5, 1);
                }
                else    //Semi-Leptonic Cuts
                {
                    selection.setJetCuts(30.,2.5,0.01,1.,0.98,0,0); //Pt, Eta, EMF, n90Hits, fHPD, dRJetElectron, DRJets
                    selection.setElectronCuts(20,2.5,.15,0.04,.5,1,0); //Pt,Eta,RelIso,d0,MVAId,DistVzPVz, DRJets, MaxMissingHits
                    selection.setMuonCuts(26, 2.1, 0.12, 0.02, 0.3, 1, 0.5, 5, 0);

                    if (debug)cout<<"Getting Jets"<<endl;
                    selectedJets                                        = selection.GetSelectedPFJets(true); // ApplyJetId
                    selectedLoosePtJets                                 = selection.GetSelectedPFJets(0, 2.5, true);
                    selectedLooseEtaJets                                = selection.GetSelectedPFJets(30, 5, true);
                    if (debug)cout<<"Getting Tight Muons"<<endl;
                    selectedMuons                                       = selection.GetSelectedMuons();
                    selectedLoosePtMuons                                = selection.GetSelectedMuons(0, 2.1, 0.12);
                    selectedLooseEtaMuons                               = selection.GetSelectedMuons(26, 5, 0.12);
                    selectedLooseIsoMuons                               = selection.GetSelectedMuons(26, 2.1, 1);
//                selection.setMuonCuts(10, 2.5, 0.20, 1, 0.3, -10000, 1, 0, 0);
                    selectedExtraMuons                                  =selection.GetSelectedMuons(10, 2.5, 0.20);
                    if (debug)cout<<"Getting Loose Electrons"<<endl;
                    selectedElectrons                                   = selection.GetSelectedElectrons(20,2.5,0.15, false, bx25, false); // VBTF ID
                    selectedLoosePtElectrons                            = selection.GetSelectedElectrons(0, 2.5, 0.15, false, bx25, false);
                    selectedLooseEtaElectrons                           = selection.GetSelectedElectrons(20, 5, 0.15, false, bx25, false);
                    selectedLooseIsoElectrons                           = selection.GetSelectedElectrons(20, 2.5, 1), false, bx25, false;
                }
            }

            else if(!Muon && Electron)   //Electron Channel Cuts
            {
                if(dilepton)   //Dilepton Cuts
                {
                    selection.setJetCuts(30.,2.5,0.01,1.,0.98,0,0); //Pt, Eta, EMF, n90Hits, fHPD, dRJetElectron, DRJets
                    selection.setElectronCuts(20,2.5,.15,0.04,.5,1,0); //Pt,Eta,RelIso,d0,MVAId,DistVzPVz, DRJets, MaxMissingHits
                    selection.setMuonCuts(20, 2.4, 0.20, 1, 0.3, 1, 1, 0, 0);

                    if (debug)cout<<"Getting Jets"<<endl;
                    selectedJets                                        = selection.GetSelectedPFJets(true); // ApplyJetId
                    selectedLoosePtJets                                 = selection.GetSelectedPFJets(0, 2.5, true);
                    selectedLooseEtaJets                                = selection.GetSelectedPFJets(30, 5, true);
                    if (debug)cout<<"Getting Tight Muons"<<endl;
                    selectedMuons                                       = selection.GetSelectedMuons();
                    selectedLoosePtMuons                                = selection.GetSelectedMuons(0, 2.4, 0.2);
                    selectedLooseEtaMuons                               = selection.GetSelectedMuons(20, 5, 0.2);
                    selectedLooseIsoMuons                               = selection.GetSelectedMuons(20, 2.4, 1);
                    if (debug)cout<<"Getting Loose Electrons"<<endl;
                    selectedElectrons                                   = selection.GetSelectedElectrons(20,2.5,0.15); // VBTF ID
                    selectedLoosePtElectrons                            = selection.GetSelectedElectrons(0, 2.5, 0.15);
                    selectedLooseEtaElectrons                           = selection.GetSelectedElectrons(20, 5, 0.15);
                    selectedLooseIsoElectrons                           = selection.GetSelectedElectrons(20, 2.5, 1);
                }
                else    //Semi-Leptonic Cuts
                {
                    selection.setJetCuts(30.,2.5,0.01,1.,0.98,0,0); //Pt, Eta, EMF, n90Hits, fHPD, dRJetElectron, DRJets
                    selection.setElectronCuts(20,2.5,.15,0.04,.5,1,0); //Pt,Eta,RelIso,d0,MVAId,DistVzPVz, DRJets, MaxMissingHits
                    selection.setMuonCuts(26, 2.1, 0.12, 0.02, 0.3, 1, 0.5, 5, 0);

                    if (debug)cout<<"Getting Jets"<<endl;
                    selectedJets                                        = selection.GetSelectedPFJets(true); // ApplyJetId
                    selectedLoosePtJets                                 = selection.GetSelectedPFJets(0, 2.5, true);
                    selectedLooseEtaJets                                = selection.GetSelectedPFJets(30, 5, true);
                    if (debug)cout<<"Getting Tight Muons"<<endl;
                    selectedMuons                                       = selection.GetSelectedMuons(20, 2.5, 0.2);
                    selectedLoosePtMuons                                = selection.GetSelectedMuons(0, 2.5, 0.2);
                    selectedLooseEtaMuons                               = selection.GetSelectedMuons(20, 5, 0.2);
                    selectedLooseIsoMuons                               = selection.GetSelectedMuons(20, 2.5, 1);
//                selection.setMuonCuts(10, 2.5, 0.20, 1, 0.3, -10000, 1, 0, 0);
                    selectedExtraMuons                                  =selection.GetSelectedMuons(10, 2.5, 0.20);
                    if (debug)cout<<"Getting Loose Electrons"<<endl;
                    selectedElectrons                                   = selection.GetSelectedElectrons(30,2.5,0.15, false, bx25, false); // Cut Based ID for CSA14
                    selectedExtraElectrons                              = selection.GetSelectedElectrons(20,2.5,0.15, false, bx25, false);
                    selectedLoosePtElectrons                            = selection.GetSelectedElectrons(0, 2.5, 0.15, false, bx25, false);
                    selectedLooseEtaElectrons                           = selection.GetSelectedElectrons(30, 5, 0.15, false, bx25, false);
                    selectedLooseIsoElectrons                           = selection.GetSelectedElectrons(30, 2.5, 1), false, bx25, false;
                }
            }
            else if (Muon && Electron && dilepton)
            {
                selection.setJetCuts(30.,2.5,0.01,1.,0.98,0,0); //Pt, Eta, EMF, n90Hits, fHPD, dRJetElectron, DRJets
                selection.setElectronCuts(20,2.5,.15,0.04,.5,1,0); //Pt,Eta,RelIso,d0,MVAId,DistVzPVz, DRJets, MaxMissingHits
                selection.setMuonCuts(20, 2.4, 0.20, 1, 0.3, 1, 1, 0, 0);

                if (debug)cout<<"Getting Jets"<<endl;
                selectedJets                                        = selection.GetSelectedPFJets(true); // ApplyJetId
                selectedLoosePtJets                                 = selection.GetSelectedPFJets(0, 2.5, true);
                selectedLooseEtaJets                                = selection.GetSelectedPFJets(30, 5, true);
                if (debug)cout<<"Getting Tight Muons"<<endl;
                selectedMuons                                       = selection.GetSelectedMuons();
                selectedLoosePtMuons                                = selection.GetSelectedMuons(0, 2.4, 0.2);
                selectedLooseEtaMuons                               = selection.GetSelectedMuons(20, 5, 0.2);
                selectedLooseIsoMuons                               = selection.GetSelectedMuons(20, 2.4, 1);
                if (debug)cout<<"Getting Loose Electrons"<<endl;
                selectedElectrons                                   = selection.GetSelectedElectrons(20,2.5,0.15); // VBTF ID
                selectedLoosePtElectrons                            = selection.GetSelectedElectrons(0, 2.5, 0.15);
                selectedLooseEtaElectrons                           = selection.GetSelectedElectrons(20, 5, 0.15);
                selectedLooseIsoElectrons                           = selection.GetSelectedElectrons(20, 2.5, 1);
            }






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

            double dR = 0, dPhi = 0, dEta = 0, lepSq = 0, projSq = 0, jetSq = 0, diMuMass = 0, looseIsoDiMuMass = 0, deltaBeta = 0;
            double dRTemp = 0, pTRelTemp = 0, sumpTTemp = 0, sumpT=0, NE = 0, CE = 0, PUCE = 0;
            float pairCharge = -9999, totalCharge = 0;
            vector<double>  MuondRMin;
            vector<double>  MuonPTRel;
            int jetPosMin = 0, nIsoMu = 0, muPos1 = 0, muPos2 = 0;
            bool pairFlag = false, looseIsoPairFlag = false;

            vector<TLorentzVector> selectedMuonsTLV_JC, selectedElectronsTLV_JC, selectedLooseIsoMuonsTLV;
            vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV, mcMuonsTLV;
            vector<TRootMCParticle*> mcParticlesMatching_;
            vector<int> mcMuonIndex;
            JetPartonMatching muonMatching;

            for (Int_t selmu =0; selmu < selectedLooseIsoMuons.size(); selmu++ )
            {
                NE = selectedLooseIsoMuons[selmu]->neutralHadronIso(4)+selectedLooseIsoMuons[selmu]->photonIso(4);
                CE = selectedLooseIsoMuons[selmu]->chargedHadronIso(4);
                PUCE = selectedLooseIsoMuons[selmu]->puChargedHadronIso(4);
                if(pow(NE,2) >= (4*PUCE*CE))
                {
                    double discrim = pow(NE,2) - (4*PUCE*CE);
                    if(sqrt(discrim) >= NE)
                    {
                        deltaBeta = (sqrt(discrim) - NE)/(2*PUCE);
                    }
                    else
                    {
                        deltaBeta = -1;
                    }
                }
                else
                {
                    deltaBeta = -1;
                }
                MSPlot["MuonBeta"]->Fill(deltaBeta, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelNHI4PreSel"]->Fill((selectedLooseIsoMuons[selmu]->neutralHadronIso(4)/selectedLooseIsoMuons[selmu]->Pt()), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelCHI4PreSel"]->Fill((CE/selectedLooseIsoMuons[selmu]->Pt()), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelPhI4PreSel"]->Fill((selectedLooseIsoMuons[selmu]->photonIso(4)/selectedLooseIsoMuons[selmu]->Pt()), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelPUCHI4PreSel"]->Fill((PUCE/selectedLooseIsoMuons[selmu]->Pt()), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelIsoPreSel"]->Fill(selectedLooseIsoMuons[selmu]->relPfIso(4, 0.5), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelTotN4PreSel"]->Fill((NE/selectedLooseIsoMuons[selmu]->Pt()), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelIso0"]->Fill(selectedLooseIsoMuons[selmu]->relPfIso(4, 0.0), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelIso1"]->Fill(selectedLooseIsoMuons[selmu]->relPfIso(4, 0.1), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelIso2"]->Fill(selectedLooseIsoMuons[selmu]->relPfIso(4, 0.2), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelIso3"]->Fill(selectedLooseIsoMuons[selmu]->relPfIso(4, 0.3), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelIso4"]->Fill(selectedLooseIsoMuons[selmu]->relPfIso(4, 0.4), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelIso5"]->Fill(selectedLooseIsoMuons[selmu]->relPfIso(4, 0.5), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelIso6"]->Fill(selectedLooseIsoMuons[selmu]->relPfIso(4, 0.6), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelIso7"]->Fill(selectedLooseIsoMuons[selmu]->relPfIso(4, 0.7), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelIso8"]->Fill(selectedLooseIsoMuons[selmu]->relPfIso(4, 0.8), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelIso9"]->Fill(selectedLooseIsoMuons[selmu]->relPfIso(4, 0.9), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelIso10"]->Fill(selectedLooseIsoMuons[selmu]->relPfIso(4, 1.0), datasets[d], true, Luminosity*scaleFactor);

                MSPlot["MuonRelNminusPUCH"]->Fill(((NE-PUCE)/selectedLooseIsoMuons[selmu]->Pt()), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelNminusBPUCH"]->Fill(((NE-0.5*PUCE)/selectedLooseIsoMuons[selmu]->Pt()), datasets[d], true, Luminosity*scaleFactor);
            }

            for (Int_t selmu =0; selmu < selectedMuons.size(); selmu++ )
            {
                MuondRMin.push_back(99999);
                selectedMuonsTLV_JC.push_back(*selectedMuons[selmu]);
                for (Int_t seljet =0; seljet<selectedJets.size(); seljet++)
                {
                    dRTemp = abs(selectedMuons[selmu]->DeltaR(*selectedJets[seljet]));
                    if (dRTemp <= MuondRMin[selmu])
                    {
                        MuondRMin[selmu] = dRTemp;
                        jetPosMin = seljet;
                    }
                }
                lepSq = pow(selectedMuons[selmu]->Px(),2) + pow(selectedMuons[selmu]->Py(),2) + pow(selectedMuons[selmu]->Pz(),2); //square of lepton momentum
                projSq = pow(selectedMuons[selmu]->Px() * selectedJets[jetPosMin]->Px() + selectedMuons[selmu]->Py() * selectedJets[jetPosMin]->Py() + selectedMuons[selmu]->Pz() * selectedJets[jetPosMin]->Pz(),2); //square of lepton momentum projected onto jet axis
                jetSq = pow(selectedJets[jetPosMin]->Px(),2) + pow(selectedJets[jetPosMin]->Py(),2) + pow(selectedJets[jetPosMin]->Pz(),2); //square of jet momentum
                pTRelTemp = sqrt( lepSq - (projSq/jetSq) );
                MuonPTRel.push_back(pTRelTemp);
                for (Int_t selmu1 =0; selmu1 < selectedMuons.size(); selmu1++ )
                {
                    if(selectedMuons[selmu]->charge() != selectedMuons[selmu1]->charge())
                    {
                        sumpTTemp = selectedMuons[selmu]->Pt() + selectedMuons[selmu1]->Pt();
                        pairFlag = true;
                        if(sumpTTemp > sumpT)
                        {
                            sumpT = sumpTTemp;
                            muPos1 = selmu;
                            muPos2 = selmu1;
                        }
                    }
                }
            }
            for (Int_t selmu =0; selmu < selectedLooseIsoMuons.size(); selmu++ )
            {
                selectedLooseIsoMuonsTLV.push_back(*selectedLooseIsoMuons[selmu]);
                for (Int_t selmu1 =0; selmu1 < selectedLooseIsoMuons.size(); selmu1++ )
                {
                    if(selectedLooseIsoMuons[selmu]->charge() != selectedLooseIsoMuons[selmu1]->charge())
                    {
                        looseIsoPairFlag = true;
                    }
                }
            }
            if(selectedMuons.size() >=2 && pairFlag)
            {
                TLorentzVector diMuon = selectedMuonsTLV_JC[muPos1] + selectedMuonsTLV_JC[muPos2];
                diMuMass = diMuon.M();
            }


            int JetCut =0;
            int nMu, nEl, nLooseIsoMu;
            if(Muon && !Electron)
            {
                nMu = selectedMuons.size(); //Number of Muons in Event
                nEl = selectedElectrons.size(); //Number of Electrons in Event
                nLooseIsoMu = selectedLooseIsoMuons.size();
            }
            if(Muon && Electron)
            {
                nMu = selectedMuons.size(); //Number of Muons in Event
                nEl = selectedElectrons.size(); //Number of Electrons in Event
            }
            if(!Muon && Electron)
            {
                nMu = selectedExtraMuons.size(); //Number of Muons in Event
                nEl = selectedElectrons.size(); //Number of Electrons in Event
            }

            bool isTagged =false;
            int seljet;

            ////////////////////////////////////////////////////
            // Looping over Jets to populate Bjet Collections //
            ////////////////////////////////////////////////////
            for ( seljet =0; seljet < selectedJets.size(); seljet++ )
            {
                if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > 0.244   )
                {
                    selectedLBJets.push_back(selectedJets[seljet]);
                    if (selectedJets[seljet]->btag_combinedInclusiveSecondaryVertexV2BJetTags() > workingpointvalue)
                    {
                        selectedMBJets.push_back(selectedJets[seljet]);
                    }
                }
                else
                {
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
            double sumCH = 0, sumNeu = 0, JERatio = 0., JMR = 0.;
            for (Int_t seljet0 =0; seljet0 < selectedJets.size(); seljet0++ )
            {
                sumCH = selectedJets[seljet0]->chargedHadronEnergyFraction() + selectedJets[seljet0]->chargedEmEnergyFraction();
                sumNeu = selectedJets[seljet0]->neutralHadronEnergyFraction() + selectedJets[seljet0]->neutralEmEnergyFraction();
                JERatio = sumCH/sumNeu;
                JMR = selectedJets[seljet0]->chargedMultiplicity()/selectedJets[seljet0]->neutralMultiplicity();
                MSPlot["JetCEF"]->Fill(sumCH, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["JetNEF"]->Fill(sumNeu, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["JetEnergyRatio"]->Fill(JERatio, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["JetMultiplicityRatio"]->Fill(JMR, datasets[d], true, Luminosity*scaleFactor);
                temp_HT += selectedJets[seljet0]->Pt();
                if (seljet0 < 4)   //Defines the leading Jets and the first 4
                {
                    HT_leading += selectedJets[seljet0]->Pt();
                }
                else
                {
                    HT_lagging += selectedJets[seljet0]->Pt();
                }
            }
            HTRat = HT_leading/HT_lagging;
            double HTb = 0.;
            for (Int_t seljet1 =0; seljet1 < nMtags; seljet1++ )
            {
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
            if(Muon && !Electron && dilepton)   //Di-Muon Selection Table
            {
                if(isGoodPV && trigged)
                {
                    selecTable.Fill(d,1,1);
                    if (nLooseIsoMu>=2 && looseIsoPairFlag)
                    {
                        selecTable.Fill(d,2,1);
                        if (nMu>=2 && pairFlag)
                        {
                            selecTable.Fill(d,3,1);
                            if(diMuMass>=20)
                            {
                                selecTable.Fill(d,4,1);
                                if(diMuMass>=106 || diMuMass<=76)
                                {
                                    selecTable.Fill(d,5,1);
                                    if(nJets>=2)
                                    {
                                        selecTable.Fill(d,6,1);
                                        if(mets[0]->Et()>=40)
                                        {
                                            selecTable.Fill(d,7,1);
                                            if(nLtags>=1)
                                            {
                                                selecTable.Fill(d,8,1);
                                                if(nLtags>=2)
                                                {
                                                    selecTable.Fill(d,9,1);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if(Muon && Electron && dilepton)   //Muon-Electron Selection Table
            {
                if(isGoodPV && trigged)
                {
                    selecTable.Fill(d,1,1);
                    if (nMu>=1)
                    {
                        selecTable.Fill(d,2,1);
                        if(nEl>=1)
                        {
                            selecTable.Fill(d,3,1);
                            if(nJets>=4)
                            {
                                selecTable.Fill(d,4,1);
                                if(nMtags>=1)
                                {
                                    selecTable.Fill(d,5,1);
                                    if(nMtags>=2)
                                    {
                                        selecTable.Fill(d,6,1);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if(!Muon && Electron && dilepton)   //Di-Electron Selection Table
            {
                if(isGoodPV && trigged)
                {
                    selecTable.Fill(d,1,1);
                    if (nEl==1)
                    {
                        selecTable.Fill(d,2,1);
                        if(nMtags==2)
                        {
                            selecTable.Fill(d,3,1);
                            if(nJets>=2)
                            {
                                selecTable.Fill(d,4,1);
                                if(nJets>=3)
                                {
                                    selecTable.Fill(d,5,1);
                                    if(nJets>=4)
                                    {
                                        selecTable.Fill(d,6,1);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if(Muon && !Electron && !dilepton)   //Single Muon Selection Table
            {
                if(isGoodPV && trigged)
                {
                    selecTable.Fill(d,1,1);
                    if (nMu==1 )
                    {
                        selecTable.Fill(d,2,1);
                        if (selectedExtraMuons.size() == 1)
                        {
                            selecTable.Fill(d,3,1);
                            if(nEl == 0)
                            {
                                selecTable.Fill(d,4,1);
                                if(nJets>=1)
                                {
                                    selecTable.Fill(d,5,1);
                                    if(nJets>=2)
                                    {
                                        selecTable.Fill(d,6,1);
                                        if(nJets>=3)
                                        {
                                            selecTable.Fill(d,7,1);
                                            if(nJets>=4)
                                            {
                                                selecTable.Fill(d,8,1);
                                                if(nMtags>=1)
                                                {
                                                    selecTable.Fill(d,9,1);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if(!Muon && Electron && !dilepton)   //Single Electron Selection Table
            {
                if(isGoodPV && trigged)
                {
                    selecTable.Fill(d,1,1);
                    if (nEl==1 )
                    {
                        selecTable.Fill(d,2,1);
                        if (selectedExtraElectrons.size() == 1)
                        {
                            selecTable.Fill(d,3,1);
                            if(nMu == 0)
                            {
                                selecTable.Fill(d,4,1);
                                if(nJets>=1)
                                {
                                    selecTable.Fill(d,5,1);
                                    if(nJets>=2)
                                    {
                                        selecTable.Fill(d,6,1);
                                        if(nJets>=3)
                                        {
                                            selecTable.Fill(d,7,1);
                                            if(nJets>=4)
                                            {
                                                selecTable.Fill(d,8,1);
                                                if(nMtags>=1)
                                                {
                                                    selecTable.Fill(d,9,1);
                                                }
                                            }
                                        }
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

            //Filling Histogram of the number of vertices before Event Selection
            MSPlot["NbOfVerticesPreSel"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);

//            if (!trigged) continue;  // Redunant check that an HLT was triggered
            if (!isGoodPV) continue; // Check that there is a good Primary Vertex
////            if (!(selectedJets.size() >= 6)) continue; //Selection of a minimum of 6 Jets in Event
//
            if (debug) cout <<"Number of Muons, Electrons, Jets, BJets, JetCut, MuonChannel, ElectronChannel ===>  "<< nMu <<"  "  <<nEl<<" "<< selectedJets.size()   <<"  " <<  nMtags   <<"  "<<JetCut  <<"  "<<Muon<<" "<<Electron<<endl;
            if (debug) cin.get();

            ///////////////////////
            // Filling N-1 Plots //
            ///////////////////////

            for (Int_t selmu =0; selmu < selectedLoosePtMuons.size(); selmu++ )
            {
                MSPlot["MuonLoosePt"]->Fill(selectedLoosePtMuons[selmu]->Pt(), datasets[d], true, Luminosity*scaleFactor);
            }
            for (Int_t selmu =0; selmu < selectedLooseEtaMuons.size(); selmu++ )
            {
                MSPlot["MuonLooseEta"]->Fill(selectedLooseEtaMuons[selmu]->Eta(), datasets[d], true, Luminosity*scaleFactor);
            }
            for (Int_t selmu =0; selmu < selectedLooseIsoMuons.size(); selmu++ )
            {
                MSPlot["MuonLooseIso"]->Fill(selectedLooseIsoMuons[selmu]->relPfIso(4, 0.5), datasets[d], true, Luminosity*scaleFactor);
            }

            for (Int_t selel =0; selel < selectedLoosePtElectrons.size(); selel++ )
            {
                MSPlot["ElectronLoosePt"]->Fill(selectedLoosePtElectrons[selel]->Pt(), datasets[d], true, Luminosity*scaleFactor);
            }
            for (Int_t selel =0; selel < selectedLooseEtaElectrons.size(); selel++ )
            {
                MSPlot["ElectronLooseEta"]->Fill(selectedLooseEtaElectrons[selel]->Eta(), datasets[d], true, Luminosity*scaleFactor);
            }
            for (Int_t selel =0; selel < selectedLooseIsoElectrons.size(); selel++ )
            {
                MSPlot["ElectronLooseIso"]->Fill(selectedLooseIsoElectrons[selel]->relPfIso(3), datasets[d], true, Luminosity*scaleFactor);
            }

            for (Int_t seljet =0; seljet < selectedLoosePtJets.size(); seljet++ )
            {
                MSPlot["JetLoosePt"]->Fill(selectedLoosePtJets[seljet]->Pt(), datasets[d], true, Luminosity*scaleFactor);
            }
            for (Int_t seljet =0; seljet < selectedLooseEtaJets.size(); seljet++ )
            {
                MSPlot["JetLooseEta"]->Fill(selectedLooseEtaJets[seljet]->Eta(), datasets[d], true, Luminosity*scaleFactor);
            }

            //Preselection Plots
            for (Int_t selmu =0; selmu < selectedMuons.size(); selmu++ )
            {
                MSPlot["MuonPtPreSel"]->Fill(selectedMuons[selmu]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelIsolationPreSel"]->Fill(selectedMuons[selmu]->relPfIso(4, 0.5), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonEtaPreSel"]->Fill(selectedMuons[selmu]->Eta(), datasets[d], true, Luminosity*scaleFactor);
            }
            for (Int_t selel =0; selel < selectedElectrons.size(); selel++ )
            {
                MSPlot["ElectronPtPreSel"]->Fill(selectedElectrons[selel]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["ElectronRelIsolationPreSel"]->Fill(selectedElectrons[selel]->relPfIso(4, 0.5), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["ElectronEtaPreSel"]->Fill(selectedElectrons[selel]->Eta(), datasets[d], true, Luminosity*scaleFactor);
            }
            for (Int_t seljet =0; seljet < selectedJets.size(); seljet++ )
            {
                MSPlot["SelectedJetPtPreSel"]->Fill(selectedJets[seljet]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["JetEtaPreSel"]->Fill(selectedJets[seljet]->Eta(), datasets[d], true, Luminosity*scaleFactor);
            }

            int nTightLeptons, nVetoLeptonsSF, nVetoLeptonsOF;
            if (debug)	cout <<" applying baseline event selection..."<<endl;
            //Apply the lepton, btag and HT selections
            if (Muon && !Electron && dilepton)
            {
                MSPlot["NbOfMuonsPreSel"]->Fill(nMu, datasets[d], true, Luminosity*scaleFactor);
                if  (  !( nMu >= 2) ) continue; // Di-Muon Channel Selection
            }
            else if (Muon && Electron && dilepton)
            {
                if  (  !( nMu >= 1 && nEl >= 1 )) continue; // Muon-Electron Channel Selection
            }
            else if (!Muon && Electron && dilepton)
            {
                if  (  !( nEl == 1 )) continue; // Di-Electron Channel Selection
            }
            else if (Muon && !Electron && !dilepton)
            {
                if ( !( nMu ==1 && selectedExtraMuons.size() == 1 && nEl == 0) ) continue;  //Single Muon Selection
            }
            else if (!Muon && Electron && !dilepton)
            {
                if ( !( nEl ==1 && selectedExtraElectrons.size() == 1 && nMu == 0) ) continue;  //Single Electron Selection
            }
            else
            {
                cerr<<"Correct Channel not selected."<<endl;
                exit(1);
            }
            sort(selectedJets.begin(),selectedJets.end(),HighestCVSBtag());

            if(nMtags>=2)
            {
                for (Int_t seljet1 =2; seljet1 < selectedJets.size(); seljet1++ )
                {

                    jetpt = selectedJets[seljet1]->Pt();
                    HT = HT + jetpt;
                    H = H + selectedJets[seljet1]->P();
                }
                MSPlot["HTExcess2M"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["HExcess2M"]->Fill(H, datasets[d], true, Luminosity*scaleFactor);
            }

            if (dilepton && Muon && !Electron)
            {
                if (!(nLtags>=2 && nJets>=2 && diMuMass >= 20 && (diMuMass>=106 || diMuMass<=76) && mets[0]->Et()>=40 )) continue; //Jet Tag Event Selection Requirements for dilepton
            }
            else if (dilepton && Muon && Electron)
            {
                if (!(nJets>=4 && nMtags >=2 )) continue; //Jet Tag Event Selection Requirements for Mu-El dilepton channel
            }
            else
            {
                if (!(nMtags>=1 && nJets>=4 )) continue; //Jet Tag Event Selection Requirements for single lepton
            }
            if(debug)
            {
                cout<<"Selection Passed."<<endl;
                cin.get();
            }
            passed++;


            ///////////////////////
            // Getting Gen Event //
            ///////////////////////

            TRootGenEvent* genEvt = 0;

            if(dataSetName != "data" && dataSetName != "Data" && dataSetName != "Data")
            {
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
//                float btag = MVASelJets1[bj1]->btag_combinedInclusiveSecondaryVertexV2BJetTags();
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
            if(nMu>=2 && Muon && !Electron)
            {
                pairCharge = selectedMuons[0]->charge() + selectedMuons[1]->charge();
                MSPlot["MuonPairCharge"]->Fill(pairCharge, datasets[d], true, Luminosity*scaleFactor);
            }

            for (Int_t selmu =0; selmu < selectedMuons.size(); selmu++ )
            {
                float reliso = selectedMuons[selmu]->relPfIso(4, 0.5);
                double muzPVz = fabs(selectedMuons[selmu]->vz() - vertex[0]->Z());
                totalCharge += selectedMuons[selmu]->charge();
                if(debug) cout << "Reliso: "<< reliso << endl;
                if(debug) cout << " chISO: " << selectedMuons[selmu]->chargedHadronIso() << endl;
                MSPlot["MuonDz"]->Fill(selectedMuons[selmu]->dz() , datasets[d], true, Luminosity*scaleFactor );
                MSPlot["MuonChi2"]->Fill(selectedMuons[selmu]->chi2() , datasets[d], true, Luminosity*scaleFactor );
                MSPlot["MuonPt"]->Fill(selectedMuons[selmu]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonPtRel"]->Fill(MuonPTRel[selmu], datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonEta"]->Fill(selectedMuons[selmu]->Eta(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonPhi"]->Fill(selectedMuons[selmu]->Phi(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonNValidHits"]->Fill(selectedMuons[selmu]->nofValidHits(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["Muond0"]->Fill(selectedMuons[selmu]->d0(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonDistVzPVz"]->Fill(muzPVz, datasets[d], true, Luminosity*scaleFactor );
                MSPlot["MuonNMatchedStations"]->Fill(selectedMuons[selmu]->nofMatchedStations(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonNPixelHits"]->Fill(selectedMuons[selmu]->nofValidPixelHits(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonTrackerLayersWithMeasurement"]->Fill(selectedMuons[selmu]->nofTrackerLayersWithMeasurement(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonDRMin"]->Fill(MuondRMin[selmu], datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MuonRelIsolation"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
//                if(reliso == 0) cout << "Zero Muon pfRelIso.  Components: " << "nHadIso: " << selectedMuons[selmu]->neutralHadronIso(4) << " phoIso: " << selectedMuons[selmu]->photonIso(4) << " PuChHadIso: " << selectedMuons[selmu]->puChargedHadronIso(4) << " chHadIso: " << selectedMuons[selmu]->chargedHadronIso(4) << endl;
//                if(selectedMuons[selmu]->chargedHadronIso(4) == 0) cout << "Zero Muon CHIso. nHadIso: " << selectedMuons[selmu]->neutralHadronIso(4) << " phoIso: " << selectedMuons[selmu]->photonIso(4) << " PuChHadIso: " << selectedMuons[selmu]->puChargedHadronIso(4) << " pfRelIso: " << reliso << endl;

            }

            if(Muon && !Electron)
            {
                //Truth matching in Di-Muon channel
                for(unsigned int i=0; i<mcParticlesMatching_.size(); i++)
                {
                    if( abs(mcParticlesMatching_[i]->type()) == 13 && mcParticlesMatching_[i]->status() == 1)   //Final State Muon
                    {
                        mcMuonsTLV.push_back(*mcParticlesMatching_[i]);
                        mcMuonIndex.push_back(i);
                    }
                }
                muonMatching = JetPartonMatching(selectedMuonsTLV_JC, mcMuonsTLV, 2, true, true, 0.3);

                MSPlot["NbOfGenMuons"]->Fill(mcMuonsTLV.size(), datasets[d], true, Luminosity*scaleFactor);

                for(unsigned int i=0; i<selectedMuonsTLV_JC.size(); i++)
                {
                    int matchedMuonNumber = muonMatching.getMatchForParton(i, 0); //Gives the index of mcMuonsTLV where the match is
                    if(matchedMuonNumber < 0) MSPlot["LeptonTruth"]->Fill(2, datasets[d], true, Luminosity*scaleFactor);
                    else MSPlot["LeptonTruth"]->Fill(1, datasets[d], true, Luminosity*scaleFactor);
                    if(abs(pairCharge!=0) && i<2 )
                    {
                        if(matchedMuonNumber >=0)
                        {
                            MSPlot["SameSignMuonType"]->Fill(abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->type()), datasets[d], true, Luminosity*scaleFactor);
                            if(abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->type()) == 13)
                            {
                                int motherType = abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->motherType());
                                if(motherType<25) MSPlot["SameSignMuonParentType"]->Fill(motherType, datasets[d], true, Luminosity*scaleFactor);
                                else MSPlot["SameSignMuonParentType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                                int grannyType = abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->grannyType());
                                if(grannyType<25) MSPlot["SameSignMuonGrannyType"]->Fill(grannyType, datasets[d], true, Luminosity*scaleFactor);
                                else MSPlot["SameSignMuonGrannyType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                            }
                        }
                        else
                        {
                            MSPlot["SameSignMuonType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                        }
                    }
                    if(abs(pairCharge==0) && i<2 )
                    {
                        if(matchedMuonNumber >=0)
                        {
                            MSPlot["OpSignMuonType"]->Fill(abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->type()), datasets[d], true, Luminosity*scaleFactor);
                            if(abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->type()) == 13)
                            {
                                int motherType = abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->motherType());
                                MSPlot["OpSignMuonParentType"]->Fill(motherType, datasets[d], true, Luminosity*scaleFactor);
                                int grannyType = abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->grannyType());
                                MSPlot["OpSignMuonGrannyType"]->Fill(grannyType, datasets[d], true, Luminosity*scaleFactor);
                            }
                            else
                            {
                                MSPlot["OpSignMuonType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                                MSPlot["OpSignMuonParentType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                                MSPlot["OpSignMuonGrannyType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                            }
                        }
                        else
                        {
                            MSPlot["OpSignMuonType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                            MSPlot["OpSignMuonParentType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                            MSPlot["OpSignMuonGrannyType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                        }
                    }

                    if(matchedMuonNumber >=0)
                    {
                        MSPlot["MuonType"]->Fill(abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->type()), datasets[d], true, Luminosity*scaleFactor);
                        int motherType = abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->motherType());
                        if(motherType<25) MSPlot["MuonParentType"]->Fill(motherType, datasets[d], true, Luminosity*scaleFactor);
                        else MSPlot["MuonParentType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                        int grannyType = abs(mcParticlesMatching_[mcMuonIndex[matchedMuonNumber]]->grannyType());
                        if(grannyType<25) MSPlot["MuonGrannyType"]->Fill(grannyType, datasets[d], true, Luminosity*scaleFactor);
                        else MSPlot["MuonGrannyType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                    }
                    else
                    {
                        MSPlot["MuonType"]->Fill(-1, datasets[d], true, Luminosity*scaleFactor);
                    }
                }
            }


//            if(pairCharge!=0 && Muon && !Electron) {
//                MSPlot["SameSignTotalCharge"]->Fill(totalCharge, datasets[d], true, Luminosity*scaleFactor);
//                MSPlot["SameSignNbOfMuons"]->Fill(nMu, datasets[d], true, Luminosity*scaleFactor);
//                selecTable.Fill(d,7,1);
//                if(nMu>2) selecTable.Fill(d,8,1);
//                }
//            if(diMuMass < 10.0) {
//                MSPlot["BadDiMuon_InvMass"]->Fill(diMuMass, datasets[d], true, Luminosity*scaleFactor);
//                }
            MSPlot["DiMuon_InvMass"]->Fill(diMuMass, datasets[d], true, Luminosity*scaleFactor);


            //////////////////////////
            // Electron Based Plots //
            //////////////////////////

            MSPlot["NbOfIsolatedElectrons"]->Fill(selectedElectrons.size(), datasets[d], true, Luminosity*scaleFactor);
            for (Int_t selel =0; selel < selectedElectrons.size(); selel++ )
            {
                float reliso = selectedElectrons[selel]->relPfIso(3, 0.5);
                double elzPVz = fabs(selectedElectrons[selel]->vz() - vertex[0]->Z());
                selectedElectronsTLV_JC.push_back(*selectedElectrons[selel]);
                MSPlot["Electrondz"]->Fill(selectedElectrons[selel]->dz() , datasets[d], true, Luminosity*scaleFactor );
                MSPlot["ElectronPt"]->Fill(selectedElectrons[selel]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["ElectronEta"]->Fill(selectedElectrons[selel]->Eta(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["ElectronPhi"]->Fill(selectedElectrons[selel]->Phi(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["Electrond0"]->Fill(selectedElectrons[selel]->d0(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["ElectronDistVzPVz"]->Fill(elzPVz, datasets[d], true, Luminosity*scaleFactor );
                MSPlot["ElectronRelIsolation"]->Fill(reliso, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["ElectronTrackChi2"]->Fill(selectedElectrons[selel]->trackNormalizedChi2(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["ElectronGSFChi2"]->Fill(selectedElectrons[selel]->gsfTrackNormalizedChi2(), datasets[d], true, Luminosity*scaleFactor);
            }
            if(nEl>=2 && !Muon && Electron)
            {
                float pairCharge = selectedElectrons[0]->charge() + selectedElectrons[1]->charge();
                MSPlot["MuonPairCharge"]->Fill(pairCharge, datasets[d], true, Luminosity*scaleFactor);
                TLorentzVector diElectron = selectedElectronsTLV_JC[0] + selectedElectronsTLV_JC[1];
                MSPlot["DiElectron_InvMass"]->Fill(diElectron.M(), datasets[d], true, Luminosity*scaleFactor);
            }


            ////////////////////////////
            // Plots for MuEl Channel //
            ////////////////////////////

            if(nEl>=1 && nMu>=1 && Muon && Electron)
            {
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

            for (Int_t seljet1 =0; seljet1 < selectedJets.size(); seljet1++ )
            {
                MSPlot["BdiscBJetCand_CSV"]->Fill(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags(),datasets[d], true, Luminosity*scaleFactor);
                MSPlot["SelectedJetPt"]->Fill(selectedJets[seljet1]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["JetEta"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
                MSPlot["JetPhi"]->Fill(selectedJets[seljet1]->Phi() , datasets[d], true, Luminosity*scaleFactor);
//                MSPlot["JetCHEF"]->Fill(selectedJets[seljet1]->chargedHadronEnergyFraction() , datasets[d], true, Luminosity*scaleFactor);
//                MSPlot["JetNHEF"]->Fill(selectedJets[seljet1]->neutralHadronEnergyFraction() , datasets[d], true, Luminosity*scaleFactor);
//                MSPlot["JetNEEF"]->Fill(selectedJets[seljet1]->neutralEmEnergyFraction() , datasets[d], true, Luminosity*scaleFactor);
                if(abs(selectedJets[seljet1]->partonFlavour()) == 5)
                {
                    MSPlot["TrueBJetEta"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["TrueBJetPt"]->Fill(selectedJets[seljet1]->Pt() , datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["TrueBJetCSV"]->Fill(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags() , datasets[d], true, Luminosity*scaleFactor);
                    if(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags() >= 0.224) MSPlot["CSVLEfficiency"]->Fill(1 , datasets[d], true, Luminosity*scaleFactor);
                    else MSPlot["CSVLEfficiency"]->Fill(0 , datasets[d], true, Luminosity*scaleFactor);
                    if(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags() >= 0.679)
                    {
                        MSPlot["CSVMEfficiency"]->Fill(1 , datasets[d], true, Luminosity*scaleFactor);
                        MSPlot["CSVMTrueBJetPt"]->Fill(selectedJets[seljet1]->Pt() , datasets[d], true, Luminosity*scaleFactor);
                        MSPlot["CSVMTrueBJetEta"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
                    }
                    else MSPlot["CSVMEfficiency"]->Fill(0 , datasets[d], true, Luminosity*scaleFactor);
                    if(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags() >= 0.898) MSPlot["CSVTEfficiency"]->Fill(1 , datasets[d], true, Luminosity*scaleFactor);
                    else MSPlot["CSVTEfficiency"]->Fill(0 , datasets[d], true, Luminosity*scaleFactor);
                }
                if(abs(selectedJets[seljet1]->partonFlavour()) == 4)
                {
                    MSPlot["TrueCJetEta"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["TrueCJetPt"]->Fill(selectedJets[seljet1]->Pt() , datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["TrueCJetCSV"]->Fill(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags() , datasets[d], true, Luminosity*scaleFactor);
                    if(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags() >= 0.679)
                    {
                        MSPlot["CSVMTrueCJetPt"]->Fill(selectedJets[seljet1]->Pt() , datasets[d], true, Luminosity*scaleFactor);
                        MSPlot["CSVMTrueCJetEta"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
                    }
                }
                if(abs(selectedJets[seljet1]->partonFlavour()) < 4)
                {
                    MSPlot["TrueLJetEta"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["TrueLJetPt"]->Fill(selectedJets[seljet1]->Pt() , datasets[d], true, Luminosity*scaleFactor);
                    MSPlot["TrueLJetCSV"]->Fill(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags() , datasets[d], true, Luminosity*scaleFactor);
                    if(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags() >= 0.679)
                    {
                        MSPlot["CSVMTrueLJetPt"]->Fill(selectedJets[seljet1]->Pt() , datasets[d], true, Luminosity*scaleFactor);
                        MSPlot["CSVMTrueLJetEta"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
                    }
                }
                if(vertex.size()<= 15)
                {
                    MSPlot["BdiscBJetCand_CSV_LowPU"]->Fill(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags(),datasets[d], true, Luminosity*scaleFactor);
                    if(abs(selectedJets[seljet1]->partonFlavour()) == 5)
                    {
                        MSPlot["TrueBJetEtaLowPU"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
                        MSPlot["TrueBJetPtLowPU"]->Fill(selectedJets[seljet1]->Pt() , datasets[d], true, Luminosity*scaleFactor);
                        MSPlot["TrueBJetCSVLowPU"]->Fill(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags() , datasets[d], true, Luminosity*scaleFactor);
                    }
                }
                if(vertex.size()>= 25)
                {
                    MSPlot["BdiscBJetCand_CSV_HighPU"]->Fill(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags(),datasets[d], true, Luminosity*scaleFactor);
                    if(abs(selectedJets[seljet1]->partonFlavour()) == 5)
                    {
                        MSPlot["TrueBJetEtaHighPU"]->Fill(selectedJets[seljet1]->Eta() , datasets[d], true, Luminosity*scaleFactor);
                        MSPlot["TrueBJetPtHighPU"]->Fill(selectedJets[seljet1]->Pt() , datasets[d], true, Luminosity*scaleFactor);
                        MSPlot["TrueBJetCSVHighPU"]->Fill(selectedJets[seljet1]->btag_combinedInclusiveSecondaryVertexV2BJetTags() , datasets[d], true, Luminosity*scaleFactor);
                    }
                }

                //Event-level variables
                jetpt = selectedJets[seljet1]->Pt();
                HT = HT + jetpt;
                H = H + selectedJets[seljet1]->P();
                sumpx = sumpx + selectedJets[seljet1]->Px();
                sumpy = sumpy + selectedJets[seljet1]->Py();
                sumpz = sumpz + selectedJets[seljet1]->Pz();
                sume = sume + selectedJets[seljet1]->E();
                if(seljet1>=3)
                {
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
//    if(TrainMVA)jetCombiner->Write(foutmva, true, pathPNGJetCombi.c_str());

//Output ROOT file
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin();
            it != MSPlot.end();
            it++)
    {
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

int Factorial(int N = 1)
{
    int fact = 1;
    for( int i=1; i<=N; i++ )
        fact = fact * i;  // OR fact *= i;
    return fact;
}



























