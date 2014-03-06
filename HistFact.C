/////////////////////////////////////////////////////////////////////////////////////
//
// Taken from http://ghl.web.cern.ch/ghl/html/HistFactoryDoc.html
//
/////////////////////////////////////////////////////////////////////////////////////


#include "RooStats/HistFactory/Measurement.h"

using namespace RooStats;

double poimin = 0.;
double poimax = 150.;
std::string infile = "SystematicShapes.root";
std::string infile_nom = "FourTop_EventSelection_El.root";
//std::string infile_Scale = "ScaleFiles/Error_MVA.root";
//std::string infile_Matching = "MatchingFiles/Error_MVA.root"; 
//std::string infile_bTag_plus = "FourTop_EventSelection_bTagPlus_Mu.root"; 
//std::string infile_bTag_minus = "FourTop_EventSelection_bTagMinus_Mu.root"; 
//std::string infile_misTag_plus = "FourTop_EventSelection_misTagPlus_Mu.root";
//std::string infile_misTag_minus = "FourTop_EventSelection_misTagMinus_Mu.root";

void HistFact(){


   // Create a simple 1-channel model
   // using c++ and ROOT

   // Create the measurement object
   // This is the top node of the structure
   // We do some minor configuration as well
   RooStats::HistFactory::Measurement meas("my_measurement", "my measurement");

   // Set the prefix that will appear before
   // all output for this measurement
   // We Set ExportOnly to false, meaning
   // we will fit the measurement and make 
   // plots in addition to saving the workspace
   meas.SetOutputFilePrefix("results/tttt");
   meas.SetExportOnly(false);

   // Set the name of the parameter of interest
   // Note that this parameter hasn't yet been
   // created, we are anticipating it
   meas.SetPOI("SigXsecOverSM");

   // Set the luminosity
   // There are a few conventions for this.
   // Here, we assume that all histograms have
   // already been scaled by luminosity
   // We also set a 10% uncertainty
   meas.SetLumi(1.0);
   meas.SetLumiRelErr(0.04); //so an uncertainty of 0.47 on a luminosity of 4.7

   // Okay, now that we've configured the measurement,
   // we'll start building the tree.
   // We begin by creating the first channel
   RooStats::HistFactory::Channel chan("channel");

   // First, we set the 'data' for this channel
   // The data is a histogram represeting the 
   // measured distribution.  It can have 1 or many bins.
   // In this example, we assume that the data histogram
   // is already made and saved in a ROOT file.  
   // So, to 'set the data', we give this channel the
   // path to that ROOT file and the name of the data
   // histogram in that root file
   // The arguments are: SetData(HistogramName, HistogramFile)
   chan.SetData("MultiSamplePlot_MVA/MVA_Data", infile_nom);

   // Now that we have a channel and have attached
   // data to it, we will start creating our Samples
   // These describe the various processes that we
   // use to model the data.
   // Here, they just consist of a signal process
   RooStats::HistFactory::Sample signal("signal", "ttbb_Up", infile);

   // Having created this sample, we configure it
   // First, we add the cross-section scaling
   // parameter that we call SigXsecOverSM
   // Then, we add a systematic with a 5% uncertainty
   // Finally, we add it to our channel
   signal.AddNormFactor("SigXsecOverSM", 20, poimin, poimax);
   signal.AddOverallSys("syst0",  0.95, 1.05);
   chan.AddSample(signal);     


   // We do a similar thing for all our backgrounds
   RooStats::HistFactory::Sample background1("background1", "MultiSamplePlot_MVA/MVA_TTJets", infile_nom);
   background1.ActivateStatError();
   background1.AddOverallSys("syst1", 0.95, 1.05 );
   background1.AddHistoSys("Scale","Scale_Down", infile,"","Scale_Up", infile,"");
   background1.AddHistoSys("Matching","Matching_Down", infile,"","Matching_Up", infile,"");
   background1.AddHistoSys("JES","JES_Down", infile,"","JES_Up", infile,"");
   background1.AddHistoSys("bTag","bTag_Down", infile,"","bTag_Up", infile,"");
   background1.AddHistoSys("misTag","misTag_Down", infile,"","misTag_Up", infile,"");
   background1.AddHistoSys("ttbb","ttbb_Down", infile,"","ttbb_Up", infile,"");

   chan.AddSample(background1);

   /*
   RooStats::HistFactory::Sample background2("background2", "MultiSamplePlot_MVA/MVA_TTJets_cc", infile);
    background2.ActivateStatError();
    background2.AddOverallSys("syst2", 0.5, 1.5 );
    chan.AddSample(background2);
    
    RooStats::HistFactory::Sample background3("background3", "MultiSamplePlot_MVA/MVA_TTJets_bb", infile);
    background3.ActivateStatError();
    background3.AddOverallSys("syst3", 0.5, 1.5 );
    chan.AddSample(background3);
   */
    
    RooStats::HistFactory::Sample background4("background4", "MultiSamplePlot_MVA/MVA_W_4Jets", infile_nom);
    background4.ActivateStatError();
    background4.AddOverallSys("syst4", 0.5, 1.5 );
    chan.AddSample(background4);
   
   
    RooStats::HistFactory::Sample background5("background5", "MultiSamplePlot_MVA/MVA_SingleTop_t_T", infile_nom);
    background5.ActivateStatError();
    background5.AddOverallSys("syst5", 0.95, 1.05 );
    chan.AddSample(background5);
    
    RooStats::HistFactory::Sample background6("background6", "MultiSamplePlot_MVA/MVA_SingleTop_t_TBar", infile_nom);
    background6.ActivateStatError();
    background6.AddOverallSys("syst6", 0.5, 1.5 );
    chan.AddSample(background6);
    
    RooStats::HistFactory::Sample background7("background7", "MultiSamplePlot_MVA/MVA_SingleTop_s_TBar", infile_nom);
    background7.ActivateStatError();
    background7.AddOverallSys("syst7", 0.95, 1.05 );
    chan.AddSample(background7);
    
    RooStats::HistFactory::Sample background8("background8", "MultiSamplePlot_MVA/MVA_SingleTop_tW_T", infile_nom);
    background8.ActivateStatError();
    background8.AddOverallSys("syst8", 0.95, 1.05 );
    chan.AddSample(background8);
    
    RooStats::HistFactory::Sample background9("background9", "MultiSamplePlot_MVA/MVA_SingleTop_tW_TBar", infile_nom);
    background9.ActivateStatError();
    background9.AddOverallSys("syst9", 0.95, 1.05 );
    chan.AddSample(background9);
    
    
    RooStats::HistFactory::Sample background10("background10", "MultiSamplePlot_MVA/MVA_SingleTop_s_T", infile_nom);
    background10.ActivateStatError();
    background10.AddOverallSys("syst10", 0.95, 1.05 );
    chan.AddSample(background10);
    
    RooStats::HistFactory::Sample background11("background11", "MultiSamplePlot_MVA/MVA_ttW", infile_nom);
    background11.ActivateStatError();
    background11.AddOverallSys("syst11", 0.95, 1.05 );
    chan.AddSample(background11);   

    RooStats::HistFactory::Sample background12("background12", "MultiSamplePlot_MVA/MVA_ttZ", infile_nom);
    background12.ActivateStatError();
    background12.AddOverallSys("syst12", 0.95, 1.05 );
    chan.AddSample(background12);

    RooStats::HistFactory::Sample background13("background13", "MultiSamplePlot_MVA/MVA_WW", infile_nom);
    background13.ActivateStatError();
    background13.AddOverallSys("syst13", 0.95, 1.05 );
    chan.AddSample(background13);

    RooStats::HistFactory::Sample background14("background14", "MultiSamplePlot_MVA/MVA_WZ", infile_nom);
    background14.ActivateStatError();
    background14.AddOverallSys("syst14", 0.95, 1.05 );
    chan.AddSample(background14);

    RooStats::HistFactory::Sample background15("background15", "MultiSamplePlot_MVA/MVA_ZZ", infile_nom);
    background15.ActivateStatError();
    background15.AddOverallSys("syst15", 0.95, 1.05 );
    chan.AddSample(background15);


   // Now that we have fully configured our channel,
   // we add it to the main measurement
   meas.AddChannel(chan);

   // At this point, we have only given our channel
   // and measurement the input histograms as strings
   // We must now have the measurement open the files,
   // collect the histograms, copy and store them.
   // This step involves I/O 
   meas.CollectHistograms();

   // Print to the screen a text representation of the model
   // just for minor debugging
   meas.PrintTree();

   // Finally, run the measurement.
   // This is the same thing that happens when
   // one runs 'hist2workspace' on an xml files
   RooStats::HistFactory::MakeModelAndMeasurementFast(meas);         
}
