#ifndef MuonTree_h
#define MuonTree_h 1

#include <fstream>   // for input/output files
#include <sstream>   // for string stream
#include <math.h>    // for sin(x) etc.
#include <cstdlib>   // for rand() on archer.
#include <iomanip>   // for setw() in cout, 
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"

#include "SC8DataStruc.h"

#define MAXGENPAR  1000
#define MAXCHAN   100

class MuonTree{
   public:
      MuonTree(int dummy);  // string outname
      ~MuonTree();  // string outname
      void analyze(SC8edep edepSc8, std::vector<SC8Particle> part, 
           std::vector<SC8Particle> hitsRef1,
           std::vector<SC8Particle> muontk);
      void endjob();

   private:
      std::string title;

      std::map<std::string, TH1D*> histo1D;
      std::map<std::string, TH1D*>::iterator histo1Diter;

      std::map<std::string, TH2D*> histo2D;
      std::map<std::string, TH2D*>::iterator histo2Diter;

      TFile *fout;

      TTree *tree;
      
      int   mTrigger;
      int   mNGenPar;
      int   mGenParId[MAXGENPAR];
      int   mGenParTk[MAXGENPAR];
      int   mGenParHit[MAXGENPAR];
      float mGenParPx[MAXGENPAR];
      float mGenParPy[MAXGENPAR];
      float mGenParPz[MAXGENPAR];
      float mGenParP[MAXGENPAR];
      float mGenParMa[MAXGENPAR];
      float mGenParVx[MAXGENPAR];
      float mGenParVy[MAXGENPAR];
      float mGenParVz[MAXGENPAR];

      int   mNBar;
      float  mEdepS1[MAXCHAN];
      float  mEdepS2[MAXCHAN];
      float  mEdepS3[MAXCHAN];
      float  mEdepS4[MAXCHAN];
      float  mSteplS1[MAXCHAN]; //changed by SAS 29/11
      float  mSteplS2[MAXCHAN]; //changed by SAS 29/11
      float  mSteplS3[MAXCHAN]; //changed by SAS 29/11
      float  mSteplS4[MAXCHAN]; //changed by SAS 29/11

      int   mNTray;
      float  mEdepT1;
      float  mEdepT2;
      float  mEdepT3;
      float  mEdepT4;
      float  mSteplT1; //changed by SAS 29/11
      float  mSteplT2; //changed by SAS 29/11
      float  mSteplT3; //changed by SAS 29/11
      float  mSteplT4; //changed by SAS 29/11
      
//      float mEdepWater;
//      float mLengthWater;
//      float mEdepWall;
//      float mLengthWall;

      int   mNHitsR1;
      int   mHitsR1Id[MAXGENPAR];
      int   mHitsR1Tk[MAXGENPAR];
      int   mHitsR1Hit[MAXGENPAR];
      float mHitsR1Px[MAXGENPAR];
      float mHitsR1Py[MAXGENPAR];
      float mHitsR1Pz[MAXGENPAR];
      float mHitsR1P[MAXGENPAR];
      float mHitsR1Ma[MAXGENPAR];
      float mHitsR1Vx[MAXGENPAR];
      float mHitsR1Vy[MAXGENPAR];
      float mHitsR1Vz[MAXGENPAR];

      int   mNMuons;
      int   mMuonId[MAXGENPAR];
      int   mMuonTk[MAXGENPAR];
      int   mMuonHit[MAXGENPAR];
      float mMuonPxE[MAXGENPAR];
      float mMuonPyE[MAXGENPAR];
      float mMuonPzE[MAXGENPAR];
      float mMuonPE[MAXGENPAR];
      float mMuonME[MAXGENPAR];
      float mMuonVxE[MAXGENPAR];
      float mMuonVyE[MAXGENPAR];
      float mMuonVzE[MAXGENPAR];
      float mMuonLength[MAXGENPAR];  // muon track elength in tank
      float mMuonEdep[MAXGENPAR];    // muon energy loss in tank
};

#endif

