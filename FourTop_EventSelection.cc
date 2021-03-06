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
#include "TopTreeAnalysisBase/Reconstruction/interface/MEzCalculator.h"
#include "TopTreeAnalysisBase/Tools/interface/LeptonTools.h"

#include "TopTreeAnalysisBase/Reconstruction/interface/TTreeObservables.h"

//This header file is taken directly from the BTV wiki. It contains
// to correctly apply an event level Btag SF. It is not yet on CVS
// as I hope to merge the functionality into BTagWeigtTools.h
//#include "TopTreeAnalysisBase/Tools/interface/BTagSFUtil.h"
#include "TopTreeAnalysisBase/Tools/interface/BTagWeightTools.h"


#include "TopTreeAnalysisBase/Tools/interface/MVATrainer.h"
#include "TopTreeAnalysisBase/Tools/interface/MVAComputer.h"
#include "TopTreeAnalysisBase/Tools/interface/JetTools.h"

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


    string postfix = "_EventSelection_Jesse"; // to relabel the names of the output file

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

    //Setting Lepton Chanels
    bool Muon = true;
    bool Electron = false;

    if(Muon) {
        cout << " --> Using the Muon channel..." << endl;
        channelpostfix = "_Mu";
        xmlFileName = "config/test_fullsamples.xml";
        }

    const char *xmlfile = xmlFileName.c_str();
    cout << "used config file: " << xmlfile << endl;

    /////////////////////////////
    //  Set up AnalysisEnvironment
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

    /////////////////////////////////
    //  Loop over Datasets
    /////////////////////////////////

    for (unsigned int d = 0; d < datasets.size (); d++) {
        cout <<"found sample with equivalent lumi "<<  datasets[d]->EquivalentLumi() <<endl;
        string dataSetName = datasets[d]->Name();
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
    MSPlot["NbOfIsolatedMuons"]                             = new MultiSamplePlot(datasets, "NbOfIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
    MSPlot["NbOfIsolatedElectrons"]                         = new MultiSamplePlot(datasets, "NbOfIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
    MSPlot["NbOfExtraIsolatedMuons"]                        = new MultiSamplePlot(datasets, "NbOfExtraIsolatedMuons", 5, 0, 5, "Nb. of isolated muons");
    MSPlot["NbOfExtraIsolatedElectrons"]                    = new MultiSamplePlot(datasets, "NbOfExtraIsolatedElectrons", 5, 0, 5, "Nb. of isolated electrons");
    MSPlot["MuonRelIsolation"]                              = new MultiSamplePlot(datasets, "MuonRelIsolation", 50, 0, .25, "RelIso");
    MSPlot["MuonRelIsolation_PreTrig"]                      = new MultiSamplePlot(datasets, "MuonRelIsolation_PreTrig", 15, 0, .25, "RelIso");
    MSPlot["MuonRelIsolation_PostTrig"]                     = new MultiSamplePlot(datasets, "MuonRelIsolation_PreTrig", 15, 0, .25, "RelIso");
    MSPlot["MuonPt"]                                        = new MultiSamplePlot(datasets, "MuonPt", 35, 0, 300, "PT_{#mu}");
    MSPlot["MuonEta"]                                       = new MultiSamplePlot(datasets, "MuonEta", 25, -2.4, 2.4, "#eta_{#mu}");
    MSPlot["MuonPhi"]                                       = new MultiSamplePlot(datasets, "MuonPhi", 50, -4, 4, "#phi_{#mu}");
    MSPlot["MuonNValidHits"]                                = new MultiSamplePlot(datasets, "MuonNValidHits", 30, 0, 30, "NValidHits_{#mu}");
    MSPlot["Muond0"]                                        = new MultiSamplePlot(datasets, "Muond0", 50, -0.05, 0.05, "d0_{#mu}");
    MSPlot["MuondZPVz"]                                     = new MultiSamplePlot(datasets, "MuondZPVz", 50, 0, .5, "dZPVZ_{#mu}");
    MSPlot["MuondRJets"]                                    = new MultiSamplePlot(datasets, "MuondRJets", 50, 0, 10, "dRJets_{#mu}");
    MSPlot["MuonNMatchedStations"]                          = new MultiSamplePlot(datasets, "MuonNMatchedStations", 10, 0, 10, "NMatchedStations_{#mu}");
    MSPlot["MuonDistVzPVz"]                                 = new MultiSamplePlot(datasets, "MuonDistVzPVz", 50, 0 ,.3, "DistVzPVz_{#mu}");
    MSPlot["MuonDz"]                                        = new MultiSamplePlot(datasets, "MuonDz", 25, -.6 ,.6, "Dz_{#mu}");
    MSPlot["MuonTrackerLayersWithMeasurement"]              = new MultiSamplePlot(datasets, "MuonTrackerLayersWithMeasurement", 25, 0, 25, "nLayers");
    MSPlot["DiMuon_InvMass"]                                = new MultiSamplePlot(datasets, "DiMuon_InvMass", 60, 0, 120, "DiMuon_InvMass");
    MSPlot["NbOfLooseMuon"]                                 = new MultiSamplePlot(datasets, "NbOfLooseMuon", 10, 0, 10, "Nb. of loose muons");
    //B-tagging discriminators
    MSPlot["BdiscBJetCand_CSV"]                             = new MultiSamplePlot(datasets, "BdiscBJetCand_CSV", 75, 0, 1, "CSV b-disc.");
    //Jets
    MSPlot["NbOfSelectedJets"]                              = new MultiSamplePlot(datasets, "NbOfSelectedJets", 15, 0, 15, "Nb. of jets");
    MSPlot["NbOfSelectedLightJets"]                         = new MultiSamplePlot(datasets, "NbOfSelectedLightJets", 10, 0, 10, "Nb. of jets");
    MSPlot["NbOfSelectedBJets"]                             = new MultiSamplePlot(datasets, "NbOfSelectedBJets", 8, 0, 8, "Nb. of tags");
    MSPlot["JetEta"]                                        = new MultiSamplePlot(datasets, "JetEta", 30,-3, 3, "Jet #eta");
    MSPlot["JetPhi"]                                        = new MultiSamplePlot(datasets, "JetPhi", 50, -4,4 , "Jet #phi");
    MSPlot["SelectedJetPt"]                                 = new MultiSamplePlot(datasets, "JetPt", 50, 0, 300, "PT_{jet}");
    MSPlot["4thJetPt"]                                      = new MultiSamplePlot(datasets, "4thJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["5thJetPt"]                                      = new MultiSamplePlot(datasets, "5thJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["6thJetPt"]                                      = new MultiSamplePlot(datasets, "6thJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["7thJetPt"]                                      = new MultiSamplePlot(datasets, "7thJetPt", 60, 0, 400, "PT_{jet}");
    MSPlot["SelectedJetPt_light"]                           = new MultiSamplePlot(datasets, "JetPt_light", 50, 0, 1000, "PT_{lightjet}");
    MSPlot["SelectedJetPt_b"]                               = new MultiSamplePlot(datasets, "JetPt_b", 50, 0, 1000, "PT_{bjet}");
    MSPlot["HT_SelectedJets"]                               = new MultiSamplePlot(datasets, "HT_SelectedJets", 50, 0, 1500, "HT");
    MSPlot["HTb_SelectedJets"]                              = new MultiSamplePlot(datasets, "HTb_SelectedJets", 50, 0, 1500, "HTb");
    MSPlot["H"]                                             = new MultiSamplePlot(datasets, "H", 50, 0, 3000, "H");
    MSPlot["HTH"]                                           = new MultiSamplePlot(datasets, "HTH", 50, 0, 1, "HT/H");
    MSPlot["HTRat"]                                         = new MultiSamplePlot(datasets, "HTRat", 50, 0, 20, "HTRat");
    //MET
    MSPlot["MET"]                                           = new MultiSamplePlot(datasets, "MET", 75, 0, 700, "MET");
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
    // Selection table: µ + jets
    /////////////////////////////////
    vector<string> CutsSelecTableMu;
    CutsSelecTableMu.push_back(string("initial"));
    CutsSelecTableMu.push_back(string("Event cleaning and Trigger"));
    CutsSelecTableMu.push_back(string("Exactly 1 isolated muon"));
    CutsSelecTableMu.push_back(string("Loose muon veto"));
    CutsSelecTableMu.push_back(string("Loose electron veto"));
    CutsSelecTableMu.push_back(string("$>$= 6 Jets (30 GeV)"));
    CutsSelecTableMu.push_back(string("$>$= 2 CSVM Tags"));
    CutsSelecTableMu.push_back(string("H_{T} $>$=  400 GeV"));
    CutsSelecTableMu.push_back(string("E^{Miss}_{T} $>$=  30 GeV"));

    SelectionTable selecTableMu(CutsSelecTableMu, datasets);
    selecTableMu.SetLuminosity(Luminosity);
    selecTableMu.SetPrecision(1);

    SelectionTable selecTableEl(CutsSelecTableMu, datasets);

    /////////////////////////////////
    // Loop on datasets
    /////////////////////////////////

    cout << " - Loop over datasets ... " << datasets.size () << " datasets !" << endl;

    for (unsigned int d = 0; d < datasets.size(); d++) {
        cout<<"Load Dataset"<<endl;
        treeLoader.LoadDataset (datasets[d], anaEnv);  //open files and load dataset
        string previousFilename = "";
        int iFile = -1;

        //////////////////////////////////////////////////
        // Initialize JEC factors ///////////////////////
        //////////////////////////////////////////////////

        vector<JetCorrectorParameters> vCorrParam;

        if(dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) { // Data!
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L1FastJet_AK5PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L2Relative_AK5PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L3Absolute_AK5PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            JetCorrectorParameters *L2L3ResJetCorPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/FT_53_V21_AN4_Summer13_Data_L2L3Residual_AK5PFchs.txt");
            vCorrParam.push_back(*L2L3ResJetCorPar);
            }
        else {
            JetCorrectorParameters *L1JetCorPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_L1FastJet_AK5PFchs.txt");
            vCorrParam.push_back(*L1JetCorPar);
            JetCorrectorParameters *L2JetCorPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_L2Relative_AK5PFchs.txt");
            vCorrParam.push_back(*L2JetCorPar);
            JetCorrectorParameters *L3JetCorPar = new JetCorrectorParameters("../../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_L3Absolute_AK5PFchs.txt");
            vCorrParam.push_back(*L3JetCorPar);
            }
        JetCorrectionUncertainty *jecUnc = new JetCorrectionUncertainty("../../TopTreeAnalysisBase/Calibrations/JECFiles/START53_V23_Summer13_Uncertainty_AK5PFchs.txt");
        JetTools *jetTools = new JetTools(vCorrParam, jecUnc, true);


        //////////////////////////////////////////////////
        // Loop on events
        /////////////////////////////////////////////////

        int itrigger = -1, previousRun = -1;

        int start = 0;
        unsigned int end = datasets[d]->NofEvtsToRunOver();

        cout <<"Number of events = "<<  end  <<endl;

        bool debug = false;
        int event_start;
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
        vector<TRootJet*>      selectedBJets;
        vector<TRootJet*>      selectedLightJets;
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
                std::cout<<"Processing the "<<ievt<<"th event, time = "<< ((double)clock() - start) / CLOCKS_PER_SEC << " ("<<100*(ievt-start)/(end-start)<<"%)"<<flush<<"\r";

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
                if(Muon) {
                    // semi-muon
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

                    else {
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

                        if(itrigger == 9999) {
                            cerr << "NO VALID TRIGGER FOUND FOR THIS EVENT (" << dataSetName << ") IN RUN " << event->runId() << endl;
                            //	exit(1);
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
            if( dataSetName.find("Data") == 0 || dataSetName.find("data") == 0 || dataSetName.find("DATA") == 0 ) {
                jetTools->unCorrectMETTypeOne(init_jets, mets[0], true);
                jetTools->correctJets(init_jets, event->kt6PFJets_rho(), true);
                jetTools->correctMETTypeOne(init_jets, mets[0], true);
                }
            else {
                if (debug)cout <<"event loop 8"<<endl;
                jetTools->unCorrectMETTypeOne(init_jets, mets[0], false);
                jetTools->correctJets(init_jets, event->kt6PFJets_rho(), false);
                jetTools->correctMETTypeOne(init_jets, mets[0], false);
                }

            ///////////////////////
            // JER smearing
            //////////////////////

            if( ! (dataSetName == "Data" || dataSetName == "data" || dataSetName == "DATA" ) ) {
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

            histo1D["lumiWeights"]->Fill(scaleFactor);

            ///////////////////////////////////////////////////////////
            // Event selection
            ///////////////////////////////////////////////////////////

            // Apply trigger selection
            trigged = treeLoader.EventTrigged (itrigger);
            if (debug)cout<<"triggered? Y/N?  "<< trigged  <<endl;
            if(!trigged)		   continue;  //If an HLT condition is not present, skip this event in the loop.
            // Declare selection instance
            Selection selection(init_jets, init_muons, init_electrons, mets);
            // Define object selection cuts
            selection.setJetCuts(30.,2.5,0.01,1.,0.98,0,0);//Pt, Eta, EMF, n90Hits, fHPD, dRJetElectron, DRJets
            selection.setElectronCuts();//	selection.setElectronCuts(30.,2.5,0.1,0.02,0.,999.,0,1); //Pt,Eta,RelIso,d0,MVAId,DistVzPVz, DRJets, MaxMissingHits
            selection.setLooseElectronCuts(20,2.5,0.15,0.);
            selection.setMuonCuts();
            selection.setLooseMuonCuts(10,2.5,0.2);
            selectedElectrons        = selection.GetSelectedElectrons();
            selectedJets             = selection.GetSelectedJets(true); // ApplyJetId
            selectedMuons            = selection.GetSelectedMuons();
            selectedLooseMuons       = selection.GetSelectedLooseMuons(10., 2.5, 0.2);
            selectedLooseElectrons   = selection.GetSelectedLooseElectrons(20.,2.5,0.15); // VBTF ID

            //order jets wrt to Pt, then set bool corresponding to RefSel cuts.
            //sort(selectedJets.begin(),selectedJets.end(),HighestPt()); //order muons wrt Pt.
            //sort(selectedCSVOrderedJets.begin(), selectedCSVOrderedJets.end(), HighestCVSBtag()); //order Jets wrt CSVtag

            int JetCut =0;
            int nMu = selectedMuons.size(); //Number of Muons in Event
            int nEl = selectedElectrons.size(); //Number of Electrons in Event

            bool isTagged =false;
            int seljet;

            ///////////////////////////////////////////////////
            // Looping over Jets to populate Bjet Collection //
            ///////////////////////////////////////////////////
            for ( seljet =0; seljet < selectedJets.size(); seljet++ ) {
                if (selectedJets[seljet]->btag_combinedSecondaryVertexBJetTags() > workingpointvalue   ) {
                    selectedBJets.push_back(selectedJets[seljet]);
                    }
                else {
                    selectedLightJets.push_back(selectedJets[seljet]);
                    }
                }
            int njets = selectedJets.size();
            int ntags = selectedBJets.size();
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
            for (Int_t seljet1 =0; seljet1 < ntags; seljet1++ ) {
                HTb += selectedBJets[seljet1]->Pt();
                }

            //////////////////////
            // Sync'ing cutflow //
            //////////////////////

            if (debug)	cout <<" applying baseline event selection..."<<endl;
            // Apply primary vertex selection
            bool isGoodPV = selection.isPVSelected(vertex, 4, 24., 2);
            if(isGoodPV && trigged) {
                selecTableMu.Fill(d,0,scaleFactor);
                if (nMu==1) {
                    selecTableMu.Fill(d,1,scaleFactor);
                    if(selectedLooseMuons.size() ==1 ) {
                        selecTableMu.Fill(d,2,scaleFactor);
                        if(selectedLooseElectrons.size() ==0 ) {
                            selecTableMu.Fill(d,3,scaleFactor);
                            if (selectedJets.size() >= 6) {
                                selecTableMu.Fill(d,4,scaleFactor);
                                if (ntags >= 2) {
                                    selecTableMu.Fill(d,5,scaleFactor);
                                    if (temp_HT >= 400.) {
                                        selecTableMu.Fill(d,6,scaleFactor);
                                        if (mets[0]->Et() >= 30.) {
                                            selecTableMu.Fill(d,7,scaleFactor);
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
            if (!trigged) continue;  // Redunant check that an HLT was triggered
            if (!isGoodPV) continue; // Check that there is a good Primary Vertex
            if (!(selectedJets.size() >= 6)) continue; //Selection of a minimum of 6 Jets in Event

            if(Muon) {  //Begin Muon Channel Selection
                if (debug) cout <<"Number of Muons, Jets, BJets, JetCut  ===>  "<< selectedMuons.size() <<"  "  << selectedJets.size()   <<"  " <<  ntags   <<"  "<<JetCut  <<endl;
                int nTightLeptons, nVetoLeptonsSF, nVetoLeptonsOF;
                if (debug)	cout <<" applying baseline event selection..."<<endl;
                //Apply the lepton, btag and HT selections
                if  (  !( nMu == 1 && selectedLooseMuons.size() == 1 && ntags >= 2 && selectedLooseElectrons.size() == 0 && temp_HT > 400. && mets[0]->Et() > 30.    )) continue;
                passed++;
                vector<TLorentzVector*> selectedMuonTLV_JC;
                selectedMuonTLV_JC.push_back(selectedMuons[0]);
                vector<TLorentzVector> mcParticlesTLV, selectedJetsTLV, mcMuonsTLV;
                vector<TRootMCParticle*> mcParticlesMatching_;

                ///////////////////////
                // Getting Gen Event //
                ///////////////////////
                vector<TRootMCParticle*> mcParticles;

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

                    }
                ///////////////////////////////////
                // Filling histograms / plotting //
                ///////////////////////////////////

                histo1D["scaleFactor"]->Fill(scaleFactor);
                MSPlot["NbOfVertices"]->Fill(vertex.size(), datasets[d], true, Luminosity*scaleFactor);



                //////////////////////
                // Muon Based Plots //
                //////////////////////

                MSPlot["NbOfIsolatedMuons"]->Fill(selectedMuons.size(), datasets[d], true, Luminosity*scaleFactor);
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

                //////////////////////
                // Jets Based Plots //
                //////////////////////

                MSPlot["NbOfSelectedJets"]->Fill(selectedJets.size(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["NbOfSelectedLightJets"]->Fill(selectedLightJets.size(), datasets[d], true, Luminosity*scaleFactor);
                MSPlot["NbOfSelectedBJets"]->Fill(ntags, datasets[d], true, Luminosity*scaleFactor);
                if (selectedJets.size()>=4) MSPlot["4thJetPt"]->Fill(selectedJets[3]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                if (selectedJets.size()>=5) MSPlot["5thJetPt"]->Fill(selectedJets[4]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                if (selectedJets.size()>=6) MSPlot["6thJetPt"]->Fill(selectedJets[5]->Pt(), datasets[d], true, Luminosity*scaleFactor);
                if (selectedJets.size()>=7) MSPlot["7thJetPt"]->Fill(selectedJets[6]->Pt(), datasets[d], true, Luminosity*scaleFactor);


                ////////////////////////////////////
                // Plotting event-level variables //
                ////////////////////////////////////

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
                    }
                HTH = HT/H;
                MSPlot["H"]->Fill(H, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["HTH"]->Fill(HTH, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["HT_SelectedJets"]->Fill(HT, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["HTRat"]->Fill(HTRat, datasets[d], true, Luminosity*scaleFactor);
                MSPlot["MET"]->Fill(mets[0]->Et(), datasets[d], true, Luminosity*scaleFactor);
            } //End Muon Channel Selection
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

    if(Muon) {
        //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
        selecTableMu.TableCalculator(  true, true, true, true, true);

        //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawsyNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
        selecTableMu.Write(  "FourTop"+postfix+"Table_Mu.tex",    false,true,true,true,false,false,false);
        }
    else if(Electron) {
        //(bool mergeTT, bool mergeQCD, bool mergeW, bool mergeZ, bool mergeST)
        selecTableEl.TableCalculator(  false, true, true, true, true);
        //Options : WithError (false), writeMerged (true), useBookTabs (false), addRawNumbers (false), addEfficiencies (false), addTotalEfficiencies (false), writeLandscape (false)
        selecTableEl.Write(  "FourTop"+postfix+"Table_El.tex",  true,true,true,true,false,false,true);

        }

    fout->cd();
    cout <<" after cd .."<<endl;

//Output ROOT file
    for(map<string,MultiSamplePlot*>::const_iterator it = MSPlot.begin(); it != MSPlot.end(); it++) {
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
