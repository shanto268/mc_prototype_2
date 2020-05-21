//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Sep 12 11:06:04 2018 by ROOT version 6.13/01
// from TTree tree/Cosmic Muon Tree
// found on file: muonTree01-912.root
//////////////////////////////////////////////////////////

#ifndef sktest_h
#define sktest_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class sktest {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nGenPar;
   Int_t           GenParId[10];   //[nGenPar]
   Int_t           GenParTk[10];   //[nGenPar]
   Int_t           GenParHit[10];   //[nGenPar]
   Float_t         GenParPx[10];   //[nGenPar]
   Float_t         GenParPy[10];   //[nGenPar]
   Float_t         GenParPz[10];   //[nGenPar]
   Float_t         GenParP[10];   //[nGenPar]
   Float_t         GenParM[10];   //[nGenPar]
   Float_t         GenParVx[10];   //[nGenPar]
   Float_t         GenParVy[10];   //[nGenPar]
   Float_t         GenParVz[10];   //[nGenPar]
   Int_t           nBar;
   Float_t         EdepS1[10];   //[nBar]
   Float_t         EdepS2[10];   //[nBar]
   Float_t         EdepS3[10];   //[nBar]
   Float_t         EdepS4[10];   //[nBar]
   Int_t           nTray;
   Float_t         EdepT1[4];   //[nTray]
   Float_t         EdepT2[4];   //[nTray]
   Float_t         EdepT3[4];   //[nTray]
   Float_t         EdepT4[4];   //[nTray]
   Int_t           nHitsR1;
   Int_t           HitsR1Id[10];   //[nHitsR1]
   Int_t           HitsR1Tk[10];   //[nHitsR1]
   Int_t           HitsR1Hit[10];   //[nHitsR1]
   Float_t         HitsR1Px[10];   //[nHitsR1]
   Float_t         HitsR1Py[10];   //[nHitsR1]
   Float_t         HitsR1Pz[10];   //[nHitsR1]
   Float_t         HitsR1P[10];   //[nHitsR1]
   Float_t         HitsR1M[10];   //[nHitsR1]
   Float_t         HitsR1Vx[10];   //[nHitsR1]
   Float_t         HitsR1Vy[10];   //[nHitsR1]
   Float_t         HitsR1Vz[10];   //[nHitsR1]
   Int_t           nMuons;
   Int_t           MuonId[2];   //[nMuons]
   Int_t           MuonTk[2];   //[nMuons]
   Int_t           MuonHit[2];   //[nMuons]
   Float_t         MuonPx[2];   //[nMuons]
   Float_t         MuonPy[2];   //[nMuons]
   Float_t         MuonPz[2];   //[nMuons]
   Float_t         MuonP[2];   //[nMuons]
   Float_t         MuonM[2];   //[nMuons]
   Float_t         MuonVx[2];   //[nMuons]
   Float_t         MuonVy[2];   //[nMuons]
   Float_t         MuonVz[2];   //[nMuons]
   Float_t         MuonLength[2];   //[nMuons]
   Float_t         MuonEdep[2];   //[nMuons]

   // List of branches
   TBranch        *b_nGenPar;   //!
   TBranch        *b_GenParId;   //!
   TBranch        *b_GenParTk;   //!
   TBranch        *b_GenParHit;   //!
   TBranch        *b_GenParPx;   //!
   TBranch        *b_GenParPy;   //!
   TBranch        *b_GenParPz;   //!
   TBranch        *b_GenParP;   //!
   TBranch        *b_GenParM;   //!
   TBranch        *b_GenParVx;   //!
   TBranch        *b_GenParVy;   //!
   TBranch        *b_GenParVz;   //!
   TBranch        *b_nBar;   //!
   TBranch        *b_EdepS1;   //!
   TBranch        *b_EdepS2;   //!
   TBranch        *b_EdepS3;   //!
   TBranch        *b_EdepS4;   //!
   TBranch        *b_nTray;   //!
   TBranch        *b_EdepT1;   //!
   TBranch        *b_EdepT2;   //!
   TBranch        *b_EdepT3;   //!
   TBranch        *b_EdepT4;   //!
   TBranch        *b_nHitsR1;   //!
   TBranch        *b_HitsR1Id;   //!
   TBranch        *b_HitsR1Tk;   //!
   TBranch        *b_HitsR1Hit;   //!
   TBranch        *b_HitsR1Px;   //!
   TBranch        *b_HitsR1Py;   //!
   TBranch        *b_HitsR1Pz;   //!
   TBranch        *b_HitsR1P;   //!
   TBranch        *b_HitsR1M;   //!
   TBranch        *b_HitsR1Vx;   //!
   TBranch        *b_HitsR1Vy;   //!
   TBranch        *b_HitsR1Vz;   //!
   TBranch        *b_nMuons;   //!
   TBranch        *b_MuonId;   //!
   TBranch        *b_MuonTk;   //!
   TBranch        *b_MuonHit;   //!
   TBranch        *b_MuonPx;   //!
   TBranch        *b_MuonPy;   //!
   TBranch        *b_MuonPz;   //!
   TBranch        *b_MuonP;   //!
   TBranch        *b_MuonM;   //!
   TBranch        *b_MuonVx;   //!
   TBranch        *b_MuonVy;   //!
   TBranch        *b_MuonVz;   //!
   TBranch        *b_MuonLength;   //!
   TBranch        *b_MuonEdep;   //!

   sktest(TTree *tree=0);
   virtual ~sktest();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef sktest_cxx
sktest::sktest(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("muonTree01-912.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("muonTree01-912.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

sktest::~sktest()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t sktest::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t sktest::LoadTree(Long64_t entry)
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

void sktest::Init(TTree *tree)
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

   fChain->SetBranchAddress("nGenPar", &nGenPar, &b_nGenPar);
   fChain->SetBranchAddress("GenParId", GenParId, &b_GenParId);
   fChain->SetBranchAddress("GenParTk", GenParTk, &b_GenParTk);
   fChain->SetBranchAddress("GenParHit", GenParHit, &b_GenParHit);
   fChain->SetBranchAddress("GenParPx", GenParPx, &b_GenParPx);
   fChain->SetBranchAddress("GenParPy", GenParPy, &b_GenParPy);
   fChain->SetBranchAddress("GenParPz", GenParPz, &b_GenParPz);
   fChain->SetBranchAddress("GenParP", GenParP, &b_GenParP);
   fChain->SetBranchAddress("GenParM", GenParM, &b_GenParM);
   fChain->SetBranchAddress("GenParVx", GenParVx, &b_GenParVx);
   fChain->SetBranchAddress("GenParVy", GenParVy, &b_GenParVy);
   fChain->SetBranchAddress("GenParVz", GenParVz, &b_GenParVz);
   fChain->SetBranchAddress("nBar", &nBar, &b_nBar);
   fChain->SetBranchAddress("EdepS1", EdepS1, &b_EdepS1);
   fChain->SetBranchAddress("EdepS2", EdepS2, &b_EdepS2);
   fChain->SetBranchAddress("EdepS3", EdepS3, &b_EdepS3);
   fChain->SetBranchAddress("EdepS4", EdepS4, &b_EdepS4);
   fChain->SetBranchAddress("nTray", &nTray, &b_nTray);
   fChain->SetBranchAddress("EdepT1", EdepT1, &b_EdepT1);
   fChain->SetBranchAddress("EdepT2", EdepT2, &b_EdepT2);
   fChain->SetBranchAddress("EdepT3", EdepT3, &b_EdepT3);
   fChain->SetBranchAddress("EdepT4", EdepT4, &b_EdepT4);
   fChain->SetBranchAddress("nHitsR1", &nHitsR1, &b_nHitsR1);
   fChain->SetBranchAddress("HitsR1Id", HitsR1Id, &b_HitsR1Id);
   fChain->SetBranchAddress("HitsR1Tk", HitsR1Tk, &b_HitsR1Tk);
   fChain->SetBranchAddress("HitsR1Hit", HitsR1Hit, &b_HitsR1Hit);
   fChain->SetBranchAddress("HitsR1Px", HitsR1Px, &b_HitsR1Px);
   fChain->SetBranchAddress("HitsR1Py", HitsR1Py, &b_HitsR1Py);
   fChain->SetBranchAddress("HitsR1Pz", HitsR1Pz, &b_HitsR1Pz);
   fChain->SetBranchAddress("HitsR1P", HitsR1P, &b_HitsR1P);
   fChain->SetBranchAddress("HitsR1M", HitsR1M, &b_HitsR1M);
   fChain->SetBranchAddress("HitsR1Vx", HitsR1Vx, &b_HitsR1Vx);
   fChain->SetBranchAddress("HitsR1Vy", HitsR1Vy, &b_HitsR1Vy);
   fChain->SetBranchAddress("HitsR1Vz", HitsR1Vz, &b_HitsR1Vz);
   fChain->SetBranchAddress("nMuons", &nMuons, &b_nMuons);
   fChain->SetBranchAddress("MuonId", MuonId, &b_MuonId);
   fChain->SetBranchAddress("MuonTk", MuonTk, &b_MuonTk);
   fChain->SetBranchAddress("MuonHit", MuonHit, &b_MuonHit);
   fChain->SetBranchAddress("MuonPx", MuonPx, &b_MuonPx);
   fChain->SetBranchAddress("MuonPy", MuonPy, &b_MuonPy);
   fChain->SetBranchAddress("MuonPz", MuonPz, &b_MuonPz);
   fChain->SetBranchAddress("MuonP", MuonP, &b_MuonP);
   fChain->SetBranchAddress("MuonM", MuonM, &b_MuonM);
   fChain->SetBranchAddress("MuonVx", MuonVx, &b_MuonVx);
   fChain->SetBranchAddress("MuonVy", MuonVy, &b_MuonVy);
   fChain->SetBranchAddress("MuonVz", MuonVz, &b_MuonVz);
   fChain->SetBranchAddress("MuonLength", MuonLength, &b_MuonLength);
   fChain->SetBranchAddress("MuonEdep", MuonEdep, &b_MuonEdep);
   Notify();
}

Bool_t sktest::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void sktest::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t sktest::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef sktest_cxx
