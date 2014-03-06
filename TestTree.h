//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Aug  8 16:08:21 2013 by ROOT version 5.34/05
// from TTree TestTree/TestTree
// found on file: MVAOutput_TTbarJES.root
//////////////////////////////////////////////////////////

#ifndef TestTree_h
#define TestTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class TestTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           classID;
   Char_t          className[11];
   Float_t         btag;
   Float_t         ThPtOverSumPt;
   Float_t         AngleThWh;
   Float_t         AngleThBh;
   Float_t         HadrWmass;
   Float_t         TopMass;
   Float_t         weight;
   Float_t         BDT;
   Float_t         Likelihood;
   Float_t         LikelihoodD;

   // List of branches
   TBranch        *b_classID;   //!
   TBranch        *b_className;   //!
   TBranch        *b_btag;   //!
   TBranch        *b_ThPtOverSumPt;   //!
   TBranch        *b_AngleThWh;   //!
   TBranch        *b_AngleThBh;   //!
   TBranch        *b_HadrWmass;   //!
   TBranch        *b_TopMass;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_BDT;   //!
   TBranch        *b_Likelihood;   //!
   TBranch        *b_LikelihoodD;   //!

   TestTree(TTree *tree=0);
   virtual ~TestTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef TestTree_cxx
TestTree::TestTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("MVAOutput_TTbarJES.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("MVAOutput_TTbarJES.root");
      }
      f->GetObject("TestTree",tree);

   }
   Init(tree);
}

TestTree::~TestTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t TestTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t TestTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void TestTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("classID", &classID, &b_classID);
   fChain->SetBranchAddress("className", className, &b_className);
   fChain->SetBranchAddress("btag", &btag, &b_btag);
   fChain->SetBranchAddress("ThPtOverSumPt", &ThPtOverSumPt, &b_ThPtOverSumPt);
   fChain->SetBranchAddress("AngleThWh", &AngleThWh, &b_AngleThWh);
   fChain->SetBranchAddress("AngleThBh", &AngleThBh, &b_AngleThBh);
   fChain->SetBranchAddress("HadrWmass", &HadrWmass, &b_HadrWmass);
   fChain->SetBranchAddress("TopMass", &TopMass, &b_TopMass);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("BDT", &BDT, &b_BDT);
   fChain->SetBranchAddress("Likelihood", &Likelihood, &b_Likelihood);
   fChain->SetBranchAddress("LikelihoodD", &LikelihoodD, &b_LikelihoodD);
   Notify();
}

Bool_t TestTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void TestTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t TestTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef TestTree_cxx
