// c++  `root-config --cflags` -o muonAnalysis3 sc8muontree.cc muonAnalysis3.cc `root-config --glibs`
//
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <cassert>
#include <vector>
#include <time.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TStyle.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#include "sc8muontree.h"

using namespace std;

//   define a struct globally.
struct TrackPoint{
   TVector3 x3;
   TVector3 p3;
};

// ==========================================================================
class AnaMuon {
   public:
      AnaMuon(TFile * fout); 
      ~AnaMuon();
      void analyze(sc8muontree ev);
      void endjob();

   private:

      vector<TrackPoint>  reconstructTracks(sc8muontree ev);
      vector<TrackPoint>  estimateXYatZ(vector<double> zPlane,TrackPoint tk );

      std::map<std::string, TH1D*> histo1D;
      std::map<std::string, TH1D*>::iterator histo1Diter;

      std::map<std::string, TH2D*> histo2D;
      std::map<std::string, TH2D*>::iterator histo2Diter;
      ofstream myfile[2]; 
      double x22[11];
      double x44[11];
      double xme;
      double angle;
      int printHitFlag;
      int event;
      vector<double> xe;
      vector<double> ang;
      vector<double> xe2;
      vector<double> ang2;
      
};

// =========================================  zzanaMuon

AnaMuon::AnaMuon(TFile * fout) {
  std::string filename[2]= {"example.txt","example1.txt",/*"example2.txt"*/};

   for (int i = 0; i < 2; i++)
	{
		myfile[i].open(filename[i].c_str());
	} 
   for (int i=0 ;i<11 ;i++){
          x22[i]=0;
          x44[i]=0;
   }
   printHitFlag=0;
   event=1;
   xme=0;
   angle=0.0;
   xe.clear();
   ang.clear();
   fout->cd();
   histo1D["GenP"]=new TH1D("GenP","Gen P (GeV)",100,0.,100.);
    histo1D["GenId"]=new TH1D("GenId","Gen Id",100,0.,20.);
   histo2D["GenXY"]=new TH2D("GenXY","Gen X vs Y (cm)",100,-100.0,100.0,100,-100.0,100.0);
   histo1D["Theta"]=new TH1D("Theta","degree",100,0.,90.);
   histo1D["Thetaxz"]=new TH1D("Thetaxz","degree",100,0.,90.);
   histo1D["Theta17"]=new TH1D("Theta17","degree",100,0.,90.);
      histo1D["Theta915"]=new TH1D("Theta915","degree",100,0.,90.);
   histo1D["Thetat"]=new TH1D("Thetat","degree",100,0.,90.);
   histo1D["Poftheta1"]=new TH1D("Poftheta1","Gen P (GeV)",100,0.,100.);
   histo1D["Poftheta2"]=new TH1D("Poftheta2","Gen P (GeV)",100,0.,100.);
   histo1D["Poftheta3"]=new TH1D("Poftheta3","Gen P (GeV)",100,0.,100.);
   histo1D["Poftheta4"]=new TH1D("Poftheta4","Gen P (GeV)",100,0.,100.);
   histo1D["ThetaofP1"]=new TH1D("ThetaofP1","degree",100,0.,90.);
   histo1D["ThetaofP2"]=new TH1D("ThetaofP2","degree",100,0.,90.);
   histo1D["ThetaofP3"]=new TH1D("ThetaofP3","degree",100,0.,90.);
   histo1D["ThetaofP4"]=new TH1D("ThetaofP4","degree",100,0.,90.);
   /*for (int i=0; i<4; i++) {
     string hname="P of theta"+to_string(i);
     string htitle="P of theta"+to_string(i);
     histo1D[hname]=new TH1D(hname.c_str(),htitle.c_str(),100,0,100);
   //  histo1D[hname]->SetAxisRange(0.,100.,"Y");
   }
   for (int i=0; i<4; i++) {
     string hname="Theta of P"+to_string(i);
     string htitle="Theta of P"+to_string(i);
     histo1D[hname]=new TH1D(hname.c_str(),htitle.c_str(),100,0, 4*atan(1));
   //  histo1D[hname]->SetAxisRange(0.,100.,"Y");
   }*/


  // histo1D["Thetay"]=new TH1D("Thetay","degree",100,0.,4*atan(1));
   histo1D["HitP"]=new TH1D("HitP","Hit P (GeV)",100,0.,100.);
   histo2D["HitXY"]=new TH2D("HitXY","Hit X vs Y (cm)",100,-100.0,100.0,100,-100.0,100.0);
   histo1D["Layer1Edep"]=new TH1D("Layer1Edep","Layer1 Edep (MeV)",100,0.0,30.0);
   histo1D["Layer1ChId"]=new TH1D("Layer1ChId","Layer1 Channel ID (Edep weighted)",10,0.0,10.0);
   histo1D["Layer2Edep1"]=new TH1D("Layer2Edep1","Layer2 Edep 1(MeV)",100,0.0,30.0);
   histo1D["Layer2Edep"]=new TH1D("Layer2Edep","Layer2 Edep (MeV)",100,0.0,30.0);
   histo1D["Layer2ChId"]=new TH1D("Layer2ChId","Layer2 Channel ID (Edep weighted)",10,0.0,10.0);
   histo1D["Layer3Edep"]=new TH1D("Layer3Edep","Layer1 Edep (MeV)",100,0.0,30.0);
   histo1D["Layer3ChId"]=new TH1D("Layer3ChId","Layer3 Channel ID (Edep weighted)",10,0.0,10.0);
   histo1D["Layer4Edep"]=new TH1D("Layer4Edep","Layer4 Edep (MeV)",100,0.0,30.0);
   histo1D["Layer4ChId"]=new TH1D("Layer4ChId","Layer4 Channel ID (Edep weighted)",10,0.0,10.0);
   histo1D["Layer2Edepstr"]=new TH1D("Layer2Edepstr","Layer2 Edep (MeV)",100,0.0,30.0);
   histo1D["Layer2Edepcur"]=new TH1D("Layer2Edepcur","Layer2 Edep (MeV)",100,0.0,30.0);
   histo1D["GenP2"]=new TH1D("GenP2","Gen P (GeV), not hit Layers",100,0.,100.);
   histo2D["GenXY2"]=new TH2D("GenXY2","Gen X vs Y (cm), not hit Layers",100,-100.0,100.0,100,-100.0,100.0);
   histo1D["GenP3"]=new TH1D("GenP3","Gen P (GeV), no hit Layers",100,0.,100.);
   histo2D["GenXY3"]=new TH2D("GenXY3","Gen X vs Y (cm),no hit Layers",100,-100.0,100.0,100,-100.0,100.0);
   histo1D["GenP4"]=new TH1D("GenP4","Gen P (MeV)",100,0.,100.0);
   histo1D["GenP5"]=new TH1D("GenP5","Gen P (MeV)",100,0.,5.0);
   histo1D["HitP2"]=new TH1D("HitP2","Hit P (MeV)",100,0.,100.0);
   histo1D["HitP3"]=new TH1D("HitP3","Hit P (MeV)",100,0.,5.0);
   histo1D["Layer2Edepmax"]=new TH1D("Layer2Edepmax","Layer2 Edep (MeV) Maxium",100,0.0,30.0);
    histo1D["Layer4Edepmax"]=new TH1D("Layer4Edepmax","Layer4 Edep (MeV) Maxium",100,0.0,30.0);
   histo2D["Layer2ET"]=new TH2D("Layer2ET","Layer2 Edep/Tol vs Tol (MeV)",100,0.0,30.0,100,0.0,30.0);
   histo2D["Layer4ET"]=new TH2D("Layer4ET","Layer4 Edep/Tol vs Tol (MeV)",100,0.0,30.0,100,0.0,30.0);
   histo2D["Layer2EH"]=new TH2D("Layer2EH","Layer2 Edep/Tol(MeV) vs #hit",100,0.0,30.0,100,0.0,15.0);
   histo2D["Layer4EH"]=new TH2D("Layer4EH","Layer4 Edep/Tol(MeV) vs #hit",100,0.0,30.0,100,0.0,15.0);
   int nz=20;
   double xylim=10000.0;
   for (int i=0; i<nz; i++) {
     string hname="xAtZ"+to_string(i);
     string htitle="x at Z"+to_string(i)+" using true track";
     histo1D[hname]=new TH1D(hname.c_str(),htitle.c_str(),100,-xylim,xylim);
     string hname2="xyAtZ"+to_string(i);
     string htitle2="x-y at Z"+to_string(i)+" using true track.";
     histo2D[hname2]=new TH2D(hname2.c_str(),htitle2.c_str(),
          50,-xylim/2.0,xylim/2.0,50,-xylim/2.0,xylim/2.0);
   //  histo1D[hname]->SetAxisRange(0.,100.,"Y");
   }

   for (int i=0; i<nz; i++) {
     string hname="xZ"+to_string(i);
     string htitle="(reco) x at Z"+to_string(i);
     histo1D[hname]=new TH1D(hname.c_str(),htitle.c_str(),100,-xylim,xylim);
     string hname2="xyZ"+to_string(i);
     string htitle2="(reco) x-y at Z"+to_string(i);
     histo2D[hname2]=new TH2D(hname2.c_str(),htitle2.c_str(),
          50,-xylim/2.0,xylim/2.0,50,-xylim/2.0,xylim/2.0);
   //  histo1D[hname]->SetAxisRange(0.,100.,"Y");
   }

   histo2D["Xrecotrue5"]=new TH2D("Xrecotrue5","X reco vs true",100,-xylim,xylim,100,-xylim,xylim);
   histo2D["Yrecotrue5"]=new TH2D("Yrecotrue5","Y reco vs true",100,-xylim,xylim,100,-xylim,xylim);
   histo1D["DeltaX"]=new TH1D("DeltaX", "X reco-X true",100,-xylim,xylim);
   histo1D["DeltaY"]=new TH1D("DeltaY", "Y reco-Y true",100,-xylim,xylim);
   histo2D["Xrecotrue3"]=new TH2D("Xrecotrue3","X reco vs true",100,-xylim,xylim,100,-xylim,xylim);
   histo2D["Yrecotrue3"]=new TH2D("Yrecotrue3","Y reco vs true",100,-xylim,xylim,100,-xylim,xylim);
   histo1D["DeltaX3"]=new TH1D("DeltaX3", "X reco-X true",100,-xylim,xylim);
   histo1D["DeltaY3"]=new TH1D("DeltaY3", "Y reco-Y true",100,-xylim,xylim);
   histo1D["DeltaX5"]=new TH1D("DeltaX5", "X reco-X true",100,-xylim,xylim);
}

//---------------------------------------------------------------------
AnaMuon::~AnaMuon() {} ;

//---------------------------------------------------------------------
void AnaMuon::analyze(sc8muontree ev) {
   double zCry;
   double xmax,xmin,ymax,ymin;

   // reconstruct tracks...
   vector<TrackPoint> recoedTracks=reconstructTracks(ev);
    xe.clear();
   ang.clear();
    xe2.clear();
    ang2.clear();
   for(int i=0; i<ev.nGenPar; i++) {
    // double aa[ev.nGenPar];
    // aa[i]=4*tan(1)*sqrt(pow(ev.GenParPx[i],2)+pow(ev.GenParPy[i],2))/(pow(ev.GenParPx[i],2)+pow(ev.GenParPy[i],2)+pow(ev.GenParPz[i],2));
     histo1D["GenP"]->Fill(ev.GenParP[i]);
     histo1D["GenId"]->Fill(ev.GenParId[i]);
     histo1D["Theta"]->Fill(atan(sqrt(pow(ev.GenParPx[i],2)+pow(ev.GenParPy[i],2))/abs(ev.GenParPz[i]))/(4*atan(1))*180.);
     
      histo1D["Thetaxz"]->Fill(atan(abs(ev.GenParPx[i]/ev.GenParPz[i]))/(4*atan(1))*180.);
     if (atan(sqrt(pow(ev.GenParPx[i],2)+pow(ev.GenParPy[i],2))/abs(ev.GenParPz[i]))/(4*atan(1))*180.<=22.5) {
     histo1D["Poftheta1"]->Fill(ev.GenParP[i]);
     }else if (atan(sqrt(pow(ev.GenParPx[i],2)+pow(ev.GenParPy[i],2))/abs(ev.GenParPz[i]))/(4*atan(1))*180.<=45.) {
     histo1D["Poftheta2"]->Fill(ev.GenParP[i]);
     }else if (atan(sqrt(pow(ev.GenParPx[i],2)+pow(ev.GenParPy[i],2))/abs(ev.GenParPz[i]))/(4*atan(1))*180.<=67.5){
     histo1D["Poftheta3"]->Fill(ev.GenParP[i]);
     }else {
     histo1D["Poftheta4"]->Fill(ev.GenParP[i]);
     } 
     if (ev.GenParP[i]<=10.0){
     histo1D["ThetaofP1"]->Fill(atan(sqrt(pow(ev.GenParPx[i],2)+pow(ev.GenParPy[i],2))/abs(ev.GenParPz[i]))/(4*atan(1))*180.);
     }else if (ev.GenParP[i]<=20.0){
     histo1D["ThetaofP2"]->Fill(atan(sqrt(pow(ev.GenParPx[i],2)+pow(ev.GenParPy[i],2))/abs(ev.GenParPz[i]))/(4*atan(1))*180.);
     }else if (ev.GenParP[i]<=30.0){
     histo1D["ThetaofP3"]->Fill(atan(sqrt(pow(ev.GenParPx[i],2)+pow(ev.GenParPy[i],2))/abs(ev.GenParPz[i]))/(4*atan(1))*180.);
     }else{
     histo1D["ThetaofP4"]->Fill(atan(sqrt(pow(ev.GenParPx[i],2)+pow(ev.GenParPy[i],2))/abs(ev.GenParPz[i]))/(4*atan(1))*180.);
     }


    
     //histo1D["Thetay"]->Fill(acos(ev.GenParPy[i]/ev.GenParPz[i]));
      double Pmev=ev.GenParP[i]*1000.0;
     histo1D["GenP4"]->Fill(Pmev);
     histo1D["GenP5"]->Fill(Pmev);
     histo2D["GenXY"]->Fill(ev.GenParVx[i],ev.GenParVy[i]);
     zCry=ev.GenParVz[i];   // z coordinate of muon source (CRY)...
   }  // end of loop over Gen Par.


  TrackPoint tkRef;
  for(int i=0; i<ev.nHitsR1; i++) {
      double Pmev2=ev.HitsR1P[i]*1000.0;
     histo1D["HitP"]->Fill(ev.HitsR1P[i]);
     histo1D["HitP2"]->Fill(Pmev2);
     histo1D["HitP3"]->Fill(Pmev2);
     histo1D["GenP3"]->Fill(ev.GenParP[i]);
     histo2D["HitXY"]->Fill(ev.HitsR1Vx[i],ev.HitsR1Vy[i]);
     // save the first hit in the refence plane...
     if(i==0) {
       if(abs(ev.HitsR1Vx[i])<25.0 && abs(ev.HitsR1Vz[i])<25.0) {
         TVector3 xx(ev.HitsR1Vx[i],ev.HitsR1Vy[i],ev.HitsR1Vz[i]);
         TVector3 pp(ev.HitsR1Px[i],ev.HitsR1Py[i],ev.HitsR1Pz[i]);
         tkRef.x3=xx;
         tkRef.p3=pp;
       }
     }
   }

   int nhits=0;
   /*if(ev.EdepS2[1]>0.0){
    histo1D["Layer2Edep1"]->Fill(ev.EdepS2[1]);
   }*/
   double xx=0;
   for(int i=0; i<11; i++) {
    // nhits=0;
     // double s1hity,s2hitx, s3hity, s4hitx;
     double edep=ev.EdepS1[i];
     double edep2=ev.EdepS2[i];
     double edep3=ev.EdepS3[i];
     double edep4=ev.EdepS4[i];
     if(edep>0.0) {
       histo1D["Layer1Edep"]->Fill(edep);
       histo1D["Layer1ChId"]->Fill(double(i),edep);
       nhits=nhits+1;
     }
      if(edep2>0.0) {
       histo1D["Layer2Edep"]->Fill(edep2);
       histo1D["Layer2ChId"]->Fill(double(i),edep2);
       xx=xx+edep2;
       for(int ii=0; ii<ev.nGenPar; ii++) {
        histo1D["Thetat"]->Fill(atan(abs(ev.GenParPx[ii]/ev.GenParPz[ii]))/(4*atan(1))*180.);
        }

       if(edep2>1.0 && edep2<7.0) { 
        for(int ii=0; ii<ev.nGenPar; ii++) {
        histo1D["Theta17"]->Fill(atan(abs(ev.GenParPx[ii]/ev.GenParPz[ii]))/(4*atan(1))*180.);
        }
       }
       else if(edep2>9.0 && edep2<15.0) {
        for(int ii=0; ii<ev.nGenPar; ii++) {
        histo1D["Theta915"]->Fill(atan(abs(ev.GenParPx[ii]/ev.GenParPz[ii]))/(4*atan(1))*180.);
        }
       }
       for(int i=0; i<ev.nGenPar; i++) {
       if(atan(abs(ev.GenParPx[i]/ev.GenParPz[i]))/(4*atan(1))*180.<1.0){
        histo1D["Layer2Edepstr"]->Fill(edep2);
       }
       else if(atan(abs(ev.GenParPx[i]/ev.GenParPz[i]))/(4*atan(1))*180.>24.5 && atan(abs(ev.GenParPx[i]/ev.GenParPz[i]))/(4*atan(1))*180.<25.5){
        histo1D["Layer2Edepcur"]->Fill(edep2);
       } 
     }
     }
     if(edep3>0.0) {
       histo1D["Layer3Edep"]->Fill(edep3);
       histo1D["Layer3ChId"]->Fill(double(i),edep3);
     }
    if(edep4>0.0) {
       histo1D["Layer4Edep"]->Fill(edep4);
       histo1D["Layer4ChId"]->Fill(double(i),edep4);
     }
   } // end of loop over scinti bars..   
   histo1D["Layer2Edep1"]->Fill(xx);
  if(nhits>0) {
     for(int i=0; i<ev.nGenPar; i++) {
       histo1D["GenP2"]->Fill(ev.GenParP[i]);
       histo2D["GenXY2"]->Fill(ev.GenParVx[i],ev.GenParVy[i]);
       // cout <<i<< "x "<<ev.GenParVx[i]<< "y "<<ev.GenParVy[i]<< "z "<<ev.GenParVz[i]<<"\n";
     } 
  }

  if(ev.nHitsR1==0) {
       for(int i=0; i<ev.nGenPar; i++) {
         histo1D["GenP3"]->Fill(ev.GenParP[i]);
         histo2D["GenXY3"]->Fill(ev.GenParVx[i],ev.GenParVy[i]);
         // cout <<i<< "x "<<ev.GenParVx[i]<< "y "<<ev.GenParVy[i]<< "z "<<ev.GenParVz[i]<<"\n";
       }
   }


//   create zplane where we check the position of tracks.
     double zRef=0.0;
     double nz=20.0;
     double dz=(10000.0-zRef)/nz;
     vector<double> zPlane;
     for (int i=0; i<nz; i++) { 
        double z=dz*i;
        zPlane.push_back(z); 
     }

//   hits in planes at different elevation (z)...
   if(ev.nHitsR1>0) {
     //  get hits in planes at different elevation (z)...  
     vector<TrackPoint>  tkZs=estimateXYatZ(zPlane,tkRef );
     for (int i=0; i<tkZs.size(); i++) {
        double x=tkZs[i].x3.x();
        double y=tkZs[i].x3.y();
       /* double px=tkZs[i].p3.x();
        double py=tkZs[i].p3.y();
        double th=acos(tkZs[i].p3.x()/tkZs[i].p3.z());*/
        string hname="xAtZ"+to_string(i);
        if(y>-500.0 && y<500.0) {
           histo1D[hname]->Fill(x);
        }
        string hname2="xyAtZ"+to_string(i);
        histo2D[hname2]->Fill(x,y);
       /* string h1="P of x"+to_string(i);
        if(y>-20.0 && y<20.0) {
           histo1D[h1]->Fill(px);
        }
        string h2="Theta"+to_string(i);
        if(y>-20.0 && y<20.0) {
           histo1D[h2]->Fill(th);
        }*/

     }   // end of loop over trueTracks in RefPlane.
  }

//  use recoedTracks...
   if(recoedTracks.size()>0) {
     //  get hits in planes at different elevation (z)...  
     //  just use the first track for nwow.
    vector<TrackPoint>  tkZs=estimateXYatZ(zPlane,recoedTracks[0]);
     for (int i=0; i<tkZs.size(); i++) {
        double x=tkZs[i].x3.x();
        double y=tkZs[i].x3.y();
       /* double px=tkZs[i].p3.x();
        double py=tkZs[i].p3.y();
        double th=acos(tkZs[i].p3.x()/tkZs[i].p3.z());*/
        if(i==5){
           xme=tkZs[i].x3.x();
           angle=atan(tkZs[i].p3.x()/tkZs[i].p3.z()); 
        }
        

        string hname="xZ"+to_string(i);
        if(y>-500.0 && y<500.0) {
           histo1D[hname]->Fill(x);
        }
        string hname2="xyZ"+to_string(i);
        histo2D[hname2]->Fill(x,y);
    }  // end of if(recoedTracks.size()>0)
    event=event+1; 

    for (int i=0;i<ev.nGenPar;i++){
     xe.push_back(ev.GenParVz[i]*(ev.GenParPx[i]/ev.GenParPz[i])+ev.GenParVx[i]);
     ang.push_back(atan(ev.GenParPx[i]/ev.GenParPz[i]));
    }
    double e2k=0;
    for (int i=0;i<11;i++){
      // if (ev.EdepS2[i]>ev.EdepS2[i-1]){
      if (ev.EdepS2[i]>e2k){
        e2k=ev.EdepS2[i];
      }
    }
    if((e2k>1.)&&(e2k<7.)){
        printHitFlag=2;
        }
    else if ((e2k>9.)&&(e2k<15.)){
        printHitFlag=1;
        }
    else{
        printHitFlag=3;
        }
    if(printHitFlag==1){
      /* vector<double> z24;
       z24.push_back(17.00);
       z24.push_back(-23.00);
       vector<TrackPoint>  tkZ24=estimateXYatZ(z24,recoedTracks[0]);*/
       string sx2[11], sx4[11];       
       vector<int> s2;
       vector<int> s4;
       for (int i=0;i<11;i++){
        sx2[i]=" . ";
        sx4[i]=" . ";
       }
       for (int ix2=0;ix2<11;ix2++){
        if(ev.EdepS2[ix2]>0.0){
         sx2[ix2]=" x ";
         s2.push_back(ix2);
         }  
       }
       for (int ix4=0;ix4<11;ix4++){
        if(ev.EdepS4[ix4]>0.0){
         sx4[ix4]=" x ";
         s4.push_back(ix4);
         }
       }
       int ix2=0;
       int ix4=0;
       double t2=0.;
       double t4=0.;
       double tt2=0.;
       double tt4=0.;
       for (int i=1;i<11;i++){
        if((ev.EdepS2[i]>=0.) && (ev.EdepS4[i]>=0.)){
         if(ev.EdepS2[i]>ev.EdepS2[i-1]){
         ix2=i;
         }
         if(ev.EdepS4[i]>ev.EdepS4[i-1]){
         ix4=i;
         }
        }
       }
       for (int i=1;i<11;i++){
         t2=t2+ev.EdepS2[i];
         t4=t4+ev.EdepS4[i];
         if (i==ix2){
          tt2=ev.EdepS2[i];
         }
         if (i==ix4){
          tt4=ev.EdepS4[i];
         }
       }
       cout<<" E "<< event <<" tt2 " <<tt2 << " t2 "<< t2 <<" tt4 "<<tt4<<" t4 "<<t4<<"\n"; 
       tt2=tt2/t2;
       tt4=tt4/t4;
      /* sx2[ix2]=" x ";
         s2.push_back(ix2);
       sx4[ix4]=" x ";
         s4.push_back(ix4);*/
         for (int i=0;i<s2.size();i++){
           xe2.push_back((((5.1*s2[i]-25.0+2.5-0.45)-(5.1*s4[i]-25.0+2.5-0.45))/40.0)*(500.0+17.0)+(5.1*s4[i]-25.0+2.5-0.45));
           ang2.push_back(atan((5.1*s2[i]-25.0+2.5-0.45)-(5.1*s4[i]-25.0+2.5-0.45))/40.0);
         }
       /*double s24,a24,s44,a44;
       s24=0.;
       s44=0.;
       a24=0.;
       a44=0.;
       for (int i=0;i<s2.size();i++){
         s24=s24+(5.1*s2[i]-25.0+2.5-0.45)*ev.EdepS2[s2[i]];
         s44=s44+ev.EdepS2[s2[i]];
      }
      for (int i2=0;i2<s4.size();i2++){
        a24=a24+(5.1*s4[i2]-25.0+2.5-0.45)*ev.EdepS4[s4[i2]];
        a44=a44+ev.EdepS4[s4[i2]];
      }
      xe2.push_back((((s24/s44)-(a24/a44))/40.0)*(500.0+17.0)+(s24/s44));
           ang2.push_back(atan(((s24/s44)-(a24/a44))/40.0));*/

      for (int i=0;i<xe2.size();i++){
        histo1D["DeltaX5"]->Fill(xe2[i]-xe[i]);
       }
      for (int i=0;i<11;i++){
        histo1D["Layer2Edepmax"]->Fill(ev.EdepS2[i]);
       }
       for (int i=0;i<11;i++){
        histo1D["Layer4Edepmax"]->Fill(ev.EdepS4[i]);
       }
      histo2D["Layer2ET"]->Fill(t2,tt2);
      histo2D["Layer4ET"]->Fill(t4,tt4);
      histo2D["Layer2EH"]->Fill(ev.nGenPar,tt2);
      histo2D["Layer4EH"]->Fill(ev.nGenPar,tt4);
     }

}
 
}  // end of AnaMuon::analyze

//=============================================================================
void AnaMuon::endjob() {
myfile[0].close();
myfile[1].close();
//myfile[2].close();
}

//=============================================================================
vector<TrackPoint>  AnaMuon::estimateXYatZ(vector<double> zPlane,TrackPoint tk ) {
   vector<TrackPoint> tkout;
    // if(abs(tk.x3.x())<25.0 && abs(tk.x3.y()) <25.0) {
     int nZplanes=zPlane.size();
     if(nZplanes>0) {
        for (int i=0; i<nZplanes; i++) {
          double dz=zPlane[i]-tk.x3.z();
          double x=dz*(tk.p3.x()/tk.p3.z())+tk.x3.x();
          double y=dz*(tk.p3.y()/tk.p3.z())+tk.x3.y();
          double z=zPlane[i];
          // cout<<"   in estimate   x"<<x<<"  y "<<y<<"  z "<<z<<endl;
          // cout<<"      zPlane[i] "<<zPlane[i]<<"  tk.x3.z() "<<tk.x3.z()<<endl; 
          TVector3 xx(x,y,z);
          
          TrackPoint tkTmp;
          tkTmp.x3=xx;
          tkTmp.p3=tk.p3;  // track vector does not change.
          tkout.push_back(tkTmp);
       // }
     }
   }  // end of if if(abs(tk.x3.x())<25.0 && abs(tk.x3.y()) <25.0) {
   return tkout;;
}

//=============================================================================
vector<TrackPoint>  AnaMuon::reconstructTracks(sc8muontree ev ) {
   vector<TrackPoint> recoedTracks;
   // layer 1 and 3 measure the y-corrdinate 
   // layer 2 and 4 measure the z-coordinate
//   double zpos[]={0.0, 56.00, 44.00, -44.0, -56.0};  // index 0 is dummy
   double zpos[]={0.0, 59.373, 51.88, -51.88, -59.373};  // index 0 is dummy
   vector<double> y1,y3;
   vector<double> x2,x4;

   //  store x and y coordinates...
   double xesum[]={0.0, 0.0, 0.0, 0.0, 0.0};
   double esum[]= {0.0, 0.0, 0.0, 0.0, 0.0};

   for(int i=0; i<11; i++) {
    // nhits=0;
     // double s1hity,s2hitx, s3hity, s4hitx;
     double edep1=ev.EdepS1[i];   // y1
     double edep2=ev.EdepS2[i];   // x2
     double edep3=ev.EdepS3[i];   // y3
     double edep4=ev.EdepS4[i];   // x4
 
     double xyval=5.1*i-25.0+2.5-0.45;

     if(edep1>0.0) {xesum[1]=xesum[1]+edep1*xyval; esum[1]=esum[1]+edep1;}
     if(edep2>0.0) {xesum[2]=xesum[2]+edep2*xyval; esum[2]=esum[2]+edep2;}
     if(edep3>0.0) {xesum[3]=xesum[3]+edep3*xyval; esum[3]=esum[3]+edep3;}
     if(edep4>0.0) {xesum[4]=xesum[4]+edep4*xyval; esum[4]=esum[4]+edep4;}
   } // end of for loop

   double a;
   if(esum[1]>0.0) {a=xesum[1]/esum[1]; y1.push_back(a);}
   if(esum[2]>0.0) {a=xesum[2]/esum[2]; x2.push_back(a);}
   if(esum[3]>0.0) {a=xesum[3]/esum[3]; y3.push_back(a);}
   if(esum[4]>0.0) {a=xesum[4]/esum[4]; x4.push_back(a);}
   // reconstruct tracks...
   // if(y1.size()>0 && y3.size()>0 && x2.size()>0 && x4.size()>0) {
   if(y1.size()==1 && y3.size()==1 && x2.size()==1 && x4.size()==1) {
      for(int iy1=0; iy1<y1.size(); iy1++) {
      for(int iy3=0; iy3<y3.size(); iy3++) {
        double yslope=(y3[iy3]-y1[iy1])/(zpos[3]-zpos[1]);
        double ypos=y1[iy1]+yslope*(0.0-zpos[1]);
        // calculate x slopes
        for(int ix2=0; ix2<x2.size(); ix2++) {
        for(int ix4=0; ix4<x4.size(); ix4++) {
          double xslope=(x4[ix4]-x2[ix2])/(zpos[4]-zpos[2]);
          double xpos=x2[ix2]+xslope*(0.0-zpos[2]);
         // myfile <<" x2 "<<ix2<<"  x4  "<<ix4<<"\n"; 
         TrackPoint tkTemp;
         TVector3 xx(xpos,ypos,0.0); 
         TVector3 pp(xslope,yslope,1.0); 
         tkTemp.x3=xx;
         tkTemp.p3=pp;         
         recoedTracks.push_back(tkTemp);
        }  // end of loop over x4
        }  // end of loop over x2.
      }  // end of loop over y3.
      }  // end of loop over y1.
   }
   return recoedTracks;
   cout<<"Slope of X="<<xslope<<endl;
}  // end of AnaMuon::reconstructTracks

// ==========================================================================
int main(int argc, char **argv) {

   string inputFileName=argv[1];
   string outputFileName=argv[2];
   cout<<"Input file name="<<inputFileName<<endl;
   cout<<"Output file name="<<outputFileName<<endl;
   //string outputFileName="o9-fem2ny.root";
   TFile *fout = new TFile(outputFileName.c_str(),"recreate");

   AnaMuon ana(fout);

   // string inputFileName="../sim/muonTree01tubig2.root";
   //string inputFileName="muonTree01-fe13m2.root";
   // string inputFileName="muonTree01.root";
   TFile fin(inputFileName.c_str());

   //   string treeName=inputFileName+":/tree";
   //TDirectory * dir = (TDirectory*) fin.Get(treeName.c_str());

   // cout<<" dir: "<<dir<<endl;

   TTree * tree;
   fin.GetObject("tree",tree);

   //   print out ntuple data structure...
   tree->Print();

   sc8muontree ev;
   ev.Init(tree);

   int nentries = tree->GetEntriesFast();
  
   cout<<"nentries  "<<nentries<<endl;

    int nMax=1000000000;
  // int nMax=100000;
   int nn=0;
   while ( ev.GetEntry(nn)){
       nn++;
       if(nn>nMax) break;
       ana.analyze(ev);
   } // end of event loop.

   ana.endjob();

   fout->Write();
   fout->Close();

}

