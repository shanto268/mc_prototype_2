#include "MuonTree.h"

#include <iostream>  // for cout
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
#include "TCanvas.h"
#include "TPad.h"
#include "TDirectory.h"
#include "TGraph.h"
#include "TEllipse.h"
#include "TText.h"
#include "TPaveText.h"

#include "SC8DataStruc.h"

using namespace std;

MuonTree::MuonTree(int dummy){

    cout<<"initializing MuonTree..."<<endl;

    string outname="muonTree01.root";
    fout=new TFile(outname.c_str(),"recreate");

    histo1D["all_nGen"]=new TH1D("all_nGen","N of Gen Particles (all)",40,0.,40.0);
    histo1D["all_GenPid"]=new TH1D("all_GenPid","Gen Particle ID (all)",40,-20.0,20.0);
    histo1D["all_GenPmu"]=new TH1D("all_GenPmu","Gen muon. P (MeV) (all)",100,0.,100.0);
    histo1D["all_GenPnonMu"]=new TH1D("all_PnonMu","Gen non-muon. P (MeV) (all)",100,0.,100.0);
    histo1D["all_GenPxPz"]=new TH1D("all_GenPxPz","slopeX, px/pz (all)",100,-2.0,2.0);
    histo1D["all_GenPyPz"]=new TH1D("all_GenPyPz","slopeY, py/pz (all)",100,-2.0,2.0);
    histo1D["all_GenPtPz"]=new TH1D("all_GenPtPz","slope, pt/pz (all)",100,0.0,2.0);

    histo1D["all_EdepT1"]=new TH1D("all_EdepT1","Edep in T1 (all)",100,0.,50.0);
    histo1D["all_EdepT2"]=new TH1D("all_EdepT2","Edep in T2 (all)",100,0.,50.0);
    histo1D["all_EdepT3"]=new TH1D("all_EdepT3","Edep in T3 (all)",100,0.,50.0);
    histo1D["all_EdepT4"]=new TH1D("all_EdepT4","Edep in T4 (all)",100,0.,50.0);
    histo1D["all_SteplT1"]=new TH1D("all_SteplT1","Muon Step Length in T1 (all) (mm)",100,0.,75.0); //changed by SAS 27/11
    histo1D["all_SteplT2"]=new TH1D("all_SteplT2","Muon Step Length in T2 (all) (mm)",100,0.,75.0); //changed by SAS 27/11
    histo1D["all_SteplT3"]=new TH1D("all_SteplT3","Muon Step Length in T3 (all) (mm)",100,0.,75.0); //changed by SAS 27/11
    histo1D["all_SteplT4"]=new TH1D("all_SteplT4","Muon Step Length in T4 (all) (mm)",100,0.,75.0); //changed by SAS 27/11

    histo1D["trig_nGen"]=new TH1D("trig_nGen","N of Gen Particles (trig)",40,0.,40.0);
    histo1D["trig_GenPid"]=new TH1D("trig_GenPid","Gen Particle ID (trig)",40,-20.0,20.0);
    histo1D["trig_GenPmu"]=new TH1D("trig_GenPmu","Gen muon. P (MeV) (trig)",100,0.,100.0);
    histo1D["trig_GenPnonMu"]=new TH1D("trig_PnonMu","Gen non-muon. P (MeV) (trig)",100,0.,100.0);
    histo1D["trig_GenPxPz"]=new TH1D("trig_GenPxPz","slopeX, px/pz (trig)",100,-2.0,4.0);
    histo1D["trig_GenPyPz"]=new TH1D("trig_GenPyPz","slopeY, py/pz (trig)",100,-2.0,2.0);
    histo1D["trig_GenPtPz"]=new TH1D("trig_GenPtPz","slope, pt/pz (trig)",100,0.0,1.0);

    histo1D["scint_muon_step"]=new TH1D("scint_muon_step","Muon Step Length (mm)",100,0.,75.0); //changed by SAS 27/11

    histo1D["trig_EdepT1"]=new TH1D("trig_EdepT1","Edep in T1 (trig)",100,0.,50.0);
    histo1D["trig_EdepT2"]=new TH1D("trig_EdepT2","Edep in T2 (trig)",100,0.,50.0);
    histo1D["trig_EdepT3"]=new TH1D("trig_EdepT3","Edep in T3 (trig)",100,0.,50.0);
    histo1D["trig_EdepT4"]=new TH1D("trig_EdepT4","Edep in T4 (trig)",100,0.,50.0);
    histo1D["trig_SteplT1"]=new TH1D("trig_SteplT1","Muon Step Length in T1 (trig) (mm)",100,0.,75.0); //changed by SAS 27/11
    histo1D["trig_SteplT2"]=new TH1D("trig_SteplT2","Muon Step Length in T2 (trig) (mm)",100,0.,75.0); //changed by SAS 27/11
    histo1D["trig_SteplT3"]=new TH1D("trig_SteplT3","Muon Step Length in T3 (trig) (mm)",100,0.,75.0); //changed by SAS 27/11
    histo1D["trig_SteplT4"]=new TH1D("trig_SteplT4","Muon Step Length in T4 (trig) (mm)",100,0.,75.0); //changed by SAS 27/11

    histo1D["trig_Tray1Muonhits"]=new TH1D("trig_Tray1Muonhits","number of hits in Tray 1",100,0.,50.0);

    histo1D["trig_muontracklengthtank"]=new TH1D("trig_muontracklengthtank","Muon TKLen in Tank",100.0,0.,100.0);  
    histo1D["trig_muonEdepinwater"]=new TH1D("trig_muonEdepinwater","Muon Edep in water",100.0,0.,100.0);  
    //   histo1D["trig_EdepWater"]=new TH1D("trig_EdepWater","Edep in Water (trig)",100,0.,50.0);
    //   histo1D["trig_LengthWater"]=new TH1D("trig_LengthWater","Track Length in Water (trig)",100,0.,50.0);
    //   histo1D["trig_EdepWall"]=new TH1D("trig_EdepWall","Edep in Wall (trig)",100,0.,50.0);
    //   histo1D["trig_LengthWall"]=new TH1D("trig_LengthWall","Track Length in Wall (trig)",100,0.,50.0);

    // define root tree...
    tree=new TTree("tree","Cosmic Muon Tree");

    tree->Branch("trigger"        ,&mTrigger     , "trigger/I" );
    tree->Branch("nGenPar"        ,&mNGenPar     , "nGenPar/I" );
    tree->Branch("GenParId"       ,mGenParId     , "GenParId[nGenPar]/I");
    tree->Branch("GenParTk"       ,mGenParTk     ,  "GenParTk[nGenPar]/I");
    tree->Branch("GenParHit"      ,mGenParHit    , "GenParHit[nGenPar]/I");
    tree->Branch("GenParPx"       ,mGenParPx     , "GenParPx[nGenPar]/F");
    tree->Branch("GenParPy"       ,mGenParPy     , "GenParPy[nGenPar]/F");
    tree->Branch("GenParPz"       ,mGenParPz     , "GenParPz[nGenPar]/F");
    tree->Branch("GenParP"        ,mGenParP      , "GenParP[nGenPar]/F");
    tree->Branch("GenParM"        ,mGenParMa     , "GenParM[nGenPar]/F");
    tree->Branch("GenParVx"       ,mGenParVx     , "GenParVx[nGenPar]/F");
    tree->Branch("GenParVy"       ,mGenParVy     , "GenParVy[nGenPar]/F");
    tree->Branch("GenParVz"       ,mGenParVz     , "GenParVz[nGenPar]/F");

    tree->Branch("nBar"      ,&mNBar     , "nBar/I" );
    tree->Branch("EdepS1"       ,mEdepS1     , "EdepS1[nBar]/F");
    tree->Branch("EdepS2"       ,mEdepS2     , "EdepS2[nBar]/F");
    tree->Branch("EdepS3"       ,mEdepS3     , "EdepS3[nBar]/F");
    tree->Branch("EdepS4"       ,mEdepS4     , "EdepS4[nBar]/F");
    tree->Branch("SteplS1"      ,mSteplS1    , "SteplS1[nBar]/F"); //Changed by SAS 29/11 
    tree->Branch("SteplS2"      ,mSteplS2    , "SteplS2[nBar]/F"); //changed by SAS 29/11 
    tree->Branch("SteplS3"      ,mSteplS3    , "SteplS3[nBar]/F"); //changed by SAS 29/11
    tree->Branch("SteplS4"      ,mSteplS4    , "SteplS4[nBar]/F"); //changed by SAS 29/11

    tree->Branch("nTray"     ,&mNTray     , "nTray/I" );
    tree->Branch("EdepT1"       ,&mEdepT1     , "EdepT1/F");
    tree->Branch("EdepT2"       ,&mEdepT2     , "EdepT2/F");
    tree->Branch("EdepT3"       ,&mEdepT3     , "EdepT3/F");
    tree->Branch("EdepT4"       ,&mEdepT4     , "EdepT4/F");
    tree->Branch("SteplT1"      ,&mSteplT1    , "SteplT1/F"); //Changed by SAS 29/11 
    tree->Branch("SteplT2"      ,&mSteplT2    , "SteplT2/F"); //changed by SAS 29/11
    tree->Branch("SteplT3"      ,&mSteplT3    , "SteplT3/F"); //changed by SAS 29/11
    tree->Branch("SteplT4"      ,&mSteplT4    , "SteplT4/F"); //changed by SAS 29/11

    //    tree->Branch("EdepWater"        ,&mEdepWater      , "EdepWater/F");
    //    tree->Branch("LengthWater"		 ,&mLengthWater    , "LengthWater/F");
    //	tree->Branch("EdepWall"		 ,&mEdepWall	   , "EdepWall");
    //	tree->Branch("LengthWall"      ,&mLengthWall    , "LengthWall");

    tree->Branch("nHitsR1"        ,&mNHitsR1     , "nHitsR1/I" );
    tree->Branch("HitsR1Id"       ,mHitsR1Id     , "HitsR1Id[nHitsR1]/I");
    tree->Branch("HitsR1Tk"       ,mHitsR1Tk     , "HitsR1Tk[nHitsR1]/I");
    tree->Branch("HitsR1Hit"      ,mHitsR1Hit    , "HitsR1Hit[nHitsR1]/I");
    tree->Branch("HitsR1Px"       ,mHitsR1Px     , "HitsR1Px[nHitsR1]/F");
    tree->Branch("HitsR1Py"       ,mHitsR1Py     , "HitsR1Py[nHitsR1]/F");
    tree->Branch("HitsR1Pz"       ,mHitsR1Pz     , "HitsR1Pz[nHitsR1]/F");
    tree->Branch("HitsR1P"        ,mHitsR1P      , "HitsR1P[nHitsR1]/F");
    tree->Branch("HitsR1M"        ,mHitsR1Ma     , "HitsR1M[nHitsR1]/F");
    tree->Branch("HitsR1Vx"       ,mHitsR1Vx     , "HitsR1Vx[nHitsR1]/F");
    tree->Branch("HitsR1Vy"       ,mHitsR1Vy     , "HitsR1Vy[nHitsR1]/F");
    tree->Branch("HitsR1Vz"       ,mHitsR1Vz     , "HitsR1Vz[nHitsR1]/F");

    tree->Branch("nMuons"        ,&mNMuons     , "nMuons/I" );
    tree->Branch("MuonId"       ,mMuonId     , "MuonId[nMuons]/I");
    tree->Branch("MuonTk"       ,mMuonTk     ,  "MuonTk[nMuons]/I");
    tree->Branch("MuonHit"      ,mMuonHit    , "MuonHit[nMuons]/I");
    tree->Branch("MuonPxE"       ,mMuonPxE     , "MuonPxE[nMuons]/F");
    tree->Branch("MuonPyE"       ,mMuonPyE     , "MuonPyE[nMuons]/F");
    tree->Branch("MuonPzE"       ,mMuonPzE     , "MuonPzE[nMuons]/F");
    tree->Branch("MuonPE"        ,mMuonPE      , "MuonPE[nMuons]/F");
    tree->Branch("MuonME"        ,mMuonME     , "MuonME[nMuons]/F");
    tree->Branch("MuonVxE"       ,mMuonVxE     , "MuonVxE[nMuons]/F");
    tree->Branch("MuonVyE"       ,mMuonVyE     , "MuonVyE[nMuons]/F");
    tree->Branch("MuonVzE"       ,mMuonVzE     , "MuonVzE[nMuons]/F");
    tree->Branch("MuonLength"      ,mMuonLength    , "MuonLength[nMuons]/F");
    tree->Branch("MuonEdep"        ,mMuonEdep      , "MuonEdep[nMuons]/F");

}

MuonTree::~MuonTree(){
}

void MuonTree::analyze(SC8edep edepSc8, vector<SC8Particle> part,
        vector<SC8Particle> hitsR1,
        vector<SC8Particle> muontk  ){

    bool trig=false;
    double ecut=1.0; // cut on edep in each layer (in MeV)
    mTrigger=0;


    for(int i=0; i<4; i++) {
        //cout << "Edep of tray " << i << " is " << edepSc8.TRAY[i] << endl;
        if(edepSc8.TRAY[i]>ecut) mTrigger++; 
    }
    if(mTrigger==4) trig=true; 

    mNGenPar=part.size();
    histo1D["all_nGen"]->Fill(mNGenPar); if(trig) histo1D["trig_nGen"]->Fill(mNGenPar);

    for(int i=0; i<mNGenPar; i++){
        mGenParId[i]=part[i].pid;
        mGenParTk[i]=part[i].trackid;
        mGenParHit[i]=0;
        mGenParPx[i]=part[i].px;
        mGenParPy[i]=part[i].py;
        mGenParPz[i]=part[i].pz;;
        mGenParP[i]=sqrt(part[i].px*part[i].px
                +part[i].py*part[i].py+part[i].pz*part[i].pz);
        mGenParMa[i]=part[i].ma;
        mGenParVx[i]=part[i].x;
        mGenParVy[i]=part[i].y;
        mGenParVz[i]=part[i].z;
        //  mMuonLength[i]=part[i].steplength;    //changed by SAS 27/11

        //   4/4 triggered events
        if(trig) {
            histo1D["trig_GenPid"]->Fill(mGenParId[i]);
            if(mGenParId[i]==13 || mGenParId[i]==-13) {
                //  histo1D["scint_muon_step"]->Fill(mMuonLength[i]); //changed by SAS 29/11
                histo1D["trig_GenPmu"]->Fill(mGenParP[i]);
                histo1D["trig_GenPxPz"]->Fill(mGenParPx[i]/mGenParPz[i]);
                histo1D["trig_GenPyPz"]->Fill(mGenParPy[i]/mGenParPz[i]);
                histo1D["trig_GenPtPz"]->Fill(sqrt(mGenParPx[i]*mGenParPx[i]+mGenParPy[i]*mGenParPy[i]/abs(mGenParPz[i])));
            }else {
                histo1D["trig_GenPnonMu"]->Fill(mGenParP[i]);
            }
        }
        // all events
        //    histo1D["all_muon_step"]->Fill(mMuonLength[i]);  //changed by SAS 29/11
        histo1D["all_GenPid"]->Fill(mGenParId[i]);
        if(mGenParId[i]==13 || mGenParId[i]==-13) {
            histo1D["all_GenPmu"]->Fill(mGenParP[i]);
            histo1D["all_GenPxPz"]->Fill(mGenParPx[i]/mGenParPz[i]);
            histo1D["all_GenPyPz"]->Fill(mGenParPy[i]/mGenParPz[i]);
            histo1D["all_GenPtPz"]->Fill(sqrt(mGenParPx[i]*mGenParPx[i]+mGenParPy[i]*mGenParPy[i])/abs(mGenParPz[i]));
        } else {
            histo1D["all_GenPnonMu"]->Fill(mGenParP[i]);
        }
    }

    
   // for(int i=0; i<2195; i++) {
   //     cout << "edep in bar number " << i << " is " << edepSc8.SBAR[i] << endl;
   // }

    mNBar=31;
    for(int i=0; i<mNBar; i++) {
        int j=i+0*31;
        mEdepS1[i]=edepSc8.SBAR[j];
        mSteplS1[i]=edepSc8.MStepBar[j]; //changed by SAS 29/11

        int k=i+1*31;
        mEdepS2[i]=edepSc8.SBAR[k];
        mSteplS2[i]=edepSc8.MStepBar[k]; //changed by SAS 29/11

        int l=i+2*31;
        mEdepS3[i]=edepSc8.SBAR[l];
        mSteplS3[i]=edepSc8.MStepBar[l]; //changed by SAS 29/11

        int h=i+3*31;
        mEdepS4[i]=edepSc8.SBAR[h];
        mSteplS4[i]=edepSc8.MStepBar[h]; //changed by SAS 29/11
    }

    /*
     * 2020-09-07: SAS test: DAQ working
    //int count = 0;
    //for (const auto& e : mEdepS1) {
       //std::cout << "Energy Deposit: " << e << " in Layer 1 num: " << count++ << std::endl;
    //}
    //int count1 =0;
    //for (const auto& e : mEdepS2) {
       //std::cout << "Energy Deposit: " << e << " in Layer 2 num: " << count1++ << std::endl;
    //}
    //int count2 =0;
    //for (const auto& e : mEdepS3) {
       //std::cout << "Energy Deposit: " << e << " in Layer 3 num: " << count2++ << std::endl;
    //}
    //int count3 =0;
    //for (const auto& e : mEdepS4) {
       //std::cout << "Energy Deposit: " << e << " in Layer 4 num: " << count3++ << std::endl;
    //}
    */

    mNTray=4;
    mEdepT1=edepSc8.TRAY[0];
    mEdepT2=edepSc8.TRAY[1];
    mEdepT3=edepSc8.TRAY[2];
    mEdepT4=edepSc8.TRAY[3];

    mSteplT1=edepSc8.MStepTray[0]; //changed by SAS 29/11
    mSteplT2=edepSc8.MStepTray[1]; //changed by SAS 29/11
    mSteplT3=edepSc8.MStepTray[2]; //changed by SAS 29/11
    mSteplT4=edepSc8.MStepTray[3]; //changed by SAS 29/11

    if(trig) {
        histo1D["trig_EdepT1"]->Fill(mEdepT1);
        histo1D["trig_EdepT1"]->GetXaxis()->SetTitle("Energy Deposit [MeV]");
        histo1D["trig_EdepT1"]->GetYaxis()->SetTitle("Number of Events");
        histo1D["trig_EdepT2"]->Fill(mEdepT2);
        histo1D["trig_EdepT2"]->GetXaxis()->SetTitle("Energy Deposit [MeV]");
        histo1D["trig_EdepT2"]->GetYaxis()->SetTitle("Number of Events");
        histo1D["trig_EdepT3"]->Fill(mEdepT3);
        histo1D["trig_EdepT3"]->GetXaxis()->SetTitle("Energy Deposit [MeV]");
        histo1D["trig_EdepT3"]->GetYaxis()->SetTitle("Number of Events");
        histo1D["trig_EdepT4"]->Fill(mEdepT4);
        histo1D["trig_EdepT4"]->GetXaxis()->SetTitle("Energy Deposit [MeV]");
        histo1D["trig_EdepT4"]->GetYaxis()->SetTitle("Number of Events");

        //SAS Comment: plotting all histo for triggered Stepl events
        histo1D["trig_SteplT1"]->Fill(mSteplT1); //changed by SAS 29/11
        histo1D["trig_SteplT1"]->GetXaxis()->SetTitle("Muon Step Length [mm]"); //changed by SAS 29/11
        histo1D["trig_SteplT1"]->GetYaxis()->SetTitle("Number of Events");//changed by SAS 29/11
        histo1D["trig_SteplT2"]->Fill(mSteplT2);//changed by SAS 29/11
        histo1D["trig_SteplT2"]->GetXaxis()->SetTitle("Muon Step Length [mm]");//changed by SAS 29/11
        histo1D["trig_SteplT2"]->GetYaxis()->SetTitle("Number of Events");//changed by SAS 29/11
        histo1D["trig_SteplT3"]->Fill(mSteplT3); //changed by SAS 29/11
        histo1D["trig_SteplT3"]->GetXaxis()->SetTitle("Muon Step Length [mm]"); //changed by SAS 29/11
        histo1D["trig_SteplT3"]->GetYaxis()->SetTitle("Number of Events");//changed by SAS 29/11
        histo1D["trig_SteplT4"]->Fill(mSteplT4);//changed by SAS 29/11
        histo1D["trig_SteplT4"]->GetXaxis()->SetTitle("Muon Step Length [mm]");//changed by SAS 29/11
        histo1D["trig_SteplT4"]->GetYaxis()->SetTitle("Number of Events");//changed by SAS 29/11
    }

    histo1D["all_EdepT1"]->Fill(mEdepT1);
    histo1D["all_EdepT2"]->Fill(mEdepT2);
    histo1D["all_EdepT3"]->Fill(mEdepT3);
    histo1D["all_EdepT4"]->Fill(mEdepT4);

    histo1D["all_SteplT1"]->Fill(mSteplT1); //changed by SAS 29/11
    histo1D["all_SteplT2"]->Fill(mSteplT2); //changed by SAS 29/11
    histo1D["all_SteplT3"]->Fill(mSteplT3); //changed by SAS 29/11
    histo1D["all_SteplT4"]->Fill(mSteplT4); //changed by SAS 29/11


    //   Hits in Tray 1
    mNHitsR1=hitsR1.size();
    if(trig) histo1D["trig_Tray1Muonhits"]->Fill(mNHitsR1);
    for(int i=0; i<mNHitsR1; i++){
        mHitsR1Id[i]=hitsR1[i].pid;
        mHitsR1Tk[i]=hitsR1[i].trackid;
        mHitsR1Hit[i]=0;
        mHitsR1Px[i]=hitsR1[i].px;
        mHitsR1Py[i]=hitsR1[i].py;
        mHitsR1Pz[i]=hitsR1[i].pz;
        mHitsR1P[i]=sqrt(hitsR1[i].px*hitsR1[i].px
                +hitsR1[i].py*hitsR1[i].py+ hitsR1[i].pz*hitsR1[i].pz);
        mHitsR1Ma[i]=hitsR1[i].ma;;
        mHitsR1Vx[i]=hitsR1[i].x;
        mHitsR1Vy[i]=hitsR1[i].y;
        mHitsR1Vz[i]=hitsR1[i].z;
    }
    // Muon particle track length???????
    mNMuons=muontk.size();
    if( mNMuons>0) {
        for(int i=0; i<mNMuons; i++){
            mMuonId[i]=muontk[i].pid;
            mMuonTk[i]=muontk[i].trackid;
            mMuonHit[i]=0;
            mMuonPxE[i]=muontk[i].px;
            mMuonPyE[i]=muontk[i].py;
            mMuonPzE[i]=muontk[i].pz;;
            mMuonPE[i]=sqrt(muontk[i].px*muontk[i].px
                    +muontk[i].py*muontk[i].py+muontk[i].pz*muontk[i].pz);
            mMuonME[i]=muontk[i].ma;
            mMuonVxE[i]=muontk[i].x;
            mMuonVyE[i]=muontk[i].y;
            mMuonVzE[i]=muontk[i].z;
            mMuonLength[i]=muontk[i].steplength;
            mMuonEdep[i]=muontk[i].edep;

            /*
            //   4/4 triggered events
            if(trig) {
            histo1D["trig_GenPid"]->Fill(mGenParId[i]);
            if(mGenParId[i]==13 || mGenParId[i]==-13) {
            histo1D["trig_GenPmu"]->Fill(mGenParP[i]);
            histo1D["trig_GenPxPz"]->Fill(mGenParPx[i]/mGenParPz[i]);
            histo1D["trig_GenPyPz"]->Fill(mGenParPy[i]/mGenParPz[i]);
            histo1D["trig_GenPtPz"]->Fill(sqrt(mGenParPx[i]*mGenParPx[i]+mGenParPy[i]*mGenParPy[i]/abs(mGenParPz[i])));
            }else {
            histo1D["trig_GenPnonMu"]->Fill(mGenParP[i]);
            }
            */
            if(trig) {
                histo1D["trig_muontracklengthtank"]->Fill(mMuonLength[i]);
                histo1D["trig_muonEdepinwater"]->Fill(mMuonEdep[i]);
            }
        }
        }
        /*  
            if(trig) {
            histo1D["trig_EdepWater"]->Fill(mEdepWater);
            histo1D["trig_EdepWall"]->Fill(mEdepWall);
            histo1D["trig_LengthWater"]->Fill(mLengthWater);
            histo1D["trig_LengthWall"]->Fill(mLengthWall);
            }
            */  
        if(trig) tree->Fill();

    }

    void MuonTree::endjob(){
        cout<<"MuonTree  endjobn..."<<endl;
        // tree->Write();
        fout->Write();
        fout->Close();
    }

