//
// c++ `root-config --cflags` -o binImage_b2  imageV02b_2.cc `root-config --glibs`
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
#include "TProfile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TStyle.h"

using namespace std;

// ====================================================================
// ====================================================================
// ====================================================================
int main( int argc, char **argv ) {


   std::map<std::string, TH1D*> histo1D;
   std::map<std::string, TH1D*>::iterator histo1Diter;

   std::map<std::string, TH2D*> histo2D;
   std::map<std::string, TH2D*>::iterator histo2Diter;
/*
  gStyle->SetCanvasColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetStatColor(0);
  gStyle->SetHistLineWidth(10);
*/

   string fname[500];  double duration[500]; double events8[500];

   fname[250]="hist5_airdeg.root"; events8[250]=500000; //at tower
   fname[251]="hist5_airdeg.root"; events8[251]=500000; //away from tower                                                                         
   duration[250] = 11.
   duration[251] = 12.

   cout<<"Duration of wtp: " << duration[250] << "  Duration of reference: " << duration[251] << endl;

   int max_index = 347;
   double sf1=events8[250]/events8[251];
   double sf2=duration[250]/duration[251];
  // cout<<"sf1="<<sf1<<"    sf2="<<sf2<<endl;
  // for (int i=319; i<max_index + 1; i++) {   
  //    cout<<"  i="<<i<<"  "<<events8[i]/duration[i]<<endl;
  // }

   gStyle->SetOptStat(0);

   TCanvas* c1 = new TCanvas("c1", "c1",600,600);


   TFile *f1 = new TFile(fname[250].c_str(),"READ");
   TFile *f2 = new TFile(fname[251].c_str(),"READ");
   //double scaleFactor=duration[250]/duration[251];
  // double scaleFactor = 1.0;	   
   // scaleFactor=20752./7696.;
   double scaleFactor=events8[250]/events8[251];  //  tower/reference
   // scaleFactor=scaleFactor*1.3;
   cout<<"scale factor="<<scaleFactor<<endl;

   if(f1==NULL || f2==NULL) 
   {cout<<"Input file open error..."<<endl;
       return 1;
   }

   histo2D["tower2d"]=(TH2D*)f1->Get("watertower/muXYproj"); 
   histo2D["ref2d"]=(TH2D*)f2->Get("watertower/muXYproj");
   histo2D["ref2d"]->Scale(scaleFactor);

   histo2D["XY"]= (TH2D*) histo2D["ref2d"]->Clone("XY");
   histo2D["XY"]->Add(histo2D["tower2d"],-1.0);
   // histo2D["XY"]->Divide(histo2D["ref2d"]);

   //histo2D["XY"]->SetMaximum(0.4);
   histo2D["XY"]->SetMinimum(0.0);
   histo2D["XY"]->SetTitle("Water Tower at Reese");
   histo2D["XY"]->GetXaxis()->SetTitle("horizontal view [m]");
   histo2D["XY"]->GetYaxis()->SetTitle("vertical view [m]");
   // histo2D["XY"]->GetXaxis()->SetTitleOffset(1.0);
   // histo2D["XY"]->GetYaxis()->SetTitleOffset(1.0);

/*
   c1->Divide(1,1,0,0);
   c1->cd(1);
   gPad->SetBottomMargin(3);
   gPad->SetRightMargin(3);
   gPad->SetTopMargin(3);
   gPad->SetLeftMargin(3);
*/

   histo2D["XY"]->Draw("LEGO2");
   // histo2D["XY"]->Draw("COLZ");


   c1->SaveAs("waterTower2_xy_2.png","png");

   TCanvas* c2 = new TCanvas("c2", "c2",600,600);

   histo1D["x1"]=(TH1D*)f1->Get("watertower/muXproj");
   histo1D["x2"]=(TH1D*)f2->Get("watertower/muXproj");
   histo1D["x2"]->Scale(scaleFactor);


   histo1D["XX"]= (TH1D*) histo1D["x1"]->Clone("XX");
   histo1D["XX"]->Add(histo1D["x2"],-1.0);
   // histo1D["XX"]->Divide(histo1D["x2"]);

   //histo1D["XX"]->SetMaximum(0.5);
   //histo1D["XX"]->SetMinimum(-0.5);
   histo1D["XX"]->SetLineWidth(3);;
   histo1D["XX"]->SetTitle("muon deficit");
   histo1D["XX"]->GetXaxis()->SetTitle("vertical band [m]");
   histo1D["XX"]->GetYaxis()->SetTitle("");
 //  histo1D["XX"]->GetXaxis()->SetLineWidth(3);


   histo1D["XX"]->Draw("HIST");
   c2->SaveAs("waterTower2_x_2.png","png");

   TCanvas* c3 = new TCanvas("c3", "c3",600,600);
   histo1D["y1"]=(TH1D*)f1->Get("watertower/muYproj");
   histo1D["y2"]=(TH1D*)f2->Get("watertower/muYproj");
   histo1D["y2"]->Scale(scaleFactor);

   histo1D["YY"]= (TH1D*) histo1D["y1"]->Clone("YY");
   histo1D["YY"]->Add(histo1D["y2"],-1.0);
   // histo1D["YY"]->Divide(histo1D["y2"]);

   // histo1D["YY"]->SetMaximum(0.5);
   // histo1D["YY"]->SetMinimum(-0.5);
   histo1D["YY"]->SetLineWidth(3);
   histo1D["YY"]->SetTitle("muon deficit");
   histo1D["YY"]->GetXaxis()->SetTitle("horizontal band [m]");
   histo1D["YY"]->GetYaxis()->SetTitle("");

   histo1D["YY"]->Draw("HIST");

   c3->SaveAs("waterTower2_y_2.png","png");

}
