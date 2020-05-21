//
// c++ `root-config --cflags` -o binImage_b1  imageV02b_1.cc `root-config --glibs`
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
   fname[319]="histWTPb_319_50_0_1_1_a.root";  duration[319]=1100.; events8[319]=293111;
   fname[320]="histWTPb_320_50_0_1_1_a.root";  duration[320]=199.;  events8[320]=51374;
   fname[321]="histWTPb_321_50_-90_1_1_a.root";  duration[321]=1106.; events8[321]=267790;
   fname[322]="histWTPb_322_50_-90_1_1_a.root";  duration[322]=483.;  events8[322]=103721;
   fname[330]="histWTPb_330_45_-90_1_1_a.root";  duration[330]=443.;  events8[330]=116306;
   fname[331]="histWTPb_331_45_0_1_1_a.root";  duration[331]=1406.; events8[331]=376797;
   fname[332]="histWTPb_332_40_0_1_1_a.root";  duration[332]=1409.; events8[332]=378238;
   fname[333]="histWTPb_333_50_0_1_1_a.root";  duration[333]=1225.; events8[333]=329260;
   fname[334]="histWTPb_334_60_0_1_1_a.root";  duration[334]=1496.; events8[334]=401128;
   fname[335]="histWTPb_335_40_-90_1_1_a.root";  duration[335]=1199.; events8[335]=311114;
   fname[336]="histWTPb_336_50_-90_1_1_a.root";  duration[336]=1333.; events8[336]=357333;
   fname[337]="histWTPb_337_60_-90_1_1_a.root";  duration[337]=400.; events8[337]=99007;
   fname[338]="histWTPb_338_60_-90_1_1_a.root";  duration[338]=1307.; events8[338]=303931;
   fname[340]="histWTPb_340_50_-90_1_1_a.root";  duration[340]=2307.; events8[340]=250771;
   fname[341]="histWTPb_341_50_0_1_1_a.root";  duration[341]= 1201.; events8[341]=292780;
   fname[342]="histWTPb_342_50_0_1_1_a.root";  duration[342]=1372.; events8[342]=267684;
   fname[343]="histWTPb_343_50_0_1_1_a.root";  duration[343]=891.; events8[343]=279705;
   fname[344]="histWTPb_344_50_0_1_1_a.root";  duration[344]=1432.; events8[344]=226516;    
   fname[345]="histWTPb_345_50_0_1_1_a.root";  duration[345]=1573.; events8[345]=407074;
   fname[347]="histWTPb_347_50_-90_1_1_a.root";  duration[347]=1344.; events8[347]=354337; 

   fname[250]="hist5_waterdeg.root"; events8[250]=events8[319]+events8[320]+events8[333]+events8[341]+events8[342]+events8[343]+events8[344]+events8[345];
   fname[251]="hist5_airdeg.root"; events8[251]=events8[321]+events8[322]+events8[336]+events8[340]+events8[347];

   duration[250]=duration[319]+duration[320]+duration[333]+duration[341]+duration[342]+duration[343]+duration[344]+duration[345];
   duration[251]=duration[321]+duration[322]+duration[336]+duration[340]+duration[347];


   cout<<"Duration of wtp: " << duration[250] << "  Duration of reference: " << duration[251] << endl;

   int max_index = 347;
   double sf1=events8[250]/events8[251];
   double sf2=duration[250]/duration[251];
   cout<<"sf1="<<sf1<<"    sf2="<<sf2<<endl;
   for (int i=319; i<max_index + 1; i++) {   
      cout<<"  i="<<i<<"  "<<events8[i]/duration[i]<<endl;
   }

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


   c1->SaveAs("waterTower2_xy_1.png","png");

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
   c2->SaveAs("waterTower2_x_1.png","png");

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

   c3->SaveAs("waterTower2_y_1.png","png");

}
