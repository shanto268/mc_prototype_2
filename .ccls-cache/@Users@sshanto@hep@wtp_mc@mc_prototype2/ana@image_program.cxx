//compiling the program
//
#include <fstream>
#include <iostream>
{
   std::map<std::string, TH1D*> histo1D;
   std::map<std::string, TH1D*>::iterator histo1Diter;

   std::map<std::string, TH2D*> histo2D;
   std::map<std::string, TH2D*>::iterator histo2Diter;

   TFile *f1 = new TFile("output_proto2_3mil_45deg_away.root","READ"); //away
   TFile *f2 = new TFile("output_proto2_3mil_45deg_tank.root","READ"); //towards

   TH1D* x7away=(TH1D*)f1->Get("xZ5");
   TH1D* x7tank=(TH1D*)f2->Get("xZ5");

   TH1D* y7away=(TH1D*)f1->Get("yZ5");
   TH1D* y7tank=(TH1D*)f2->Get("yZ5");

   TH2D* xy7away=(TH2D*)f1->Get("xyZ5");
   TH2D* xy7tank=(TH2D*)f2->Get("xyZ5");
/*
   int nx1d=x7away->GetNbinsX();
   double xmin1d=-200.0;
   double xmax1d= 200.0;
   double dx=(xmax1d-xmin1d)/double(nx1d);
   cout<<"nx1d="<<nx1d<<"    dx="<<dx<<endl;

   TH1D* plotXtank=new TH1D("xTank","X Tank",nx1d,xmin1d,xmax1d);
   for (int ix=1; ix<nx1d; ix++) {
      double wt=x7tank->GetBinContent(ix);
	  double x=xmax1d-dx*ix+dx/2.0;
      plotXtank->Fill(x,wt);
	  //cout<<" xvalue="<<x<<endl;
   }
   
   int ny1d=y7away->GetNbinsX();
   double ymin1d=-200.0;
   double ymax1d= 200.0;
   double dy=(ymax1d-ymin1d)/double(ny1d);

   TH1D* plotYtank=new TH1D("yTank","Y Tank",ny1d,ymin1d,ymax1d);
   for (int iy=1; iy<ny1d; iy++) {
		double wt=y7tank->GetBinContent(iy);
		double y=ymax1d-dy*iy+dy/2.0;
		plotYtank->Fill(y,wt);
   }
  */
   // xy9tank->Divide(xy9away); 

   int nx=xy7away->GetNbinsY();
   int ny=xy7away->GetNbinsX();
   cout<<"nx="<<nx<<"   ny="<<ny<<endl;

//   TH2D* plotA=new TH2D("plotA","MC (7ft)",13,-15.0,12.0,13,-15.0,12.0);
   TH2D* plotA=new TH2D("plotA","MC (Water 5m)",15,-15.0,15.0,15,-15.0,15.0);

   ofstream out ("wt16ft.txt");
   for (int ix=1; ix<nx; ix++) {
      for(int iy=1; iy<ny; iy++) {
         double wtA=xy7away->GetBinContent(ix,iy);//collecting bin contents for everybin in 2D projection plot of away
         //double wtT=xy7tank->GetBinContent(nx-ix+2,iy);
         double wtT=xy7tank->GetBinContent(ix,iy);//collection bin contents for everybin in 2D projection plot of watertank
	 double wt=wtT/wtA;//taking ratio
         //cout<<"wt = "<<wtT<<endl;
	 if(wt>1.0) wt=1.0;
          //cout<<"wt "<<wt<<"\n";
          double y=20.0-2.0*ix+2.0/2.0;
          double x=20.0-2.0*iy+2.0/2.0;
          //double x=1500-150*ix+150;
          //double y=1500-150*iy+150;
          //double x=2.5*(ix-26)-2.5/2.0;
          //double y=2.5*(iy-15)-2.5/2.0;
         plotA->Fill(x,y,wt);
	 out<<x<<" "<<y<<" "<<wt<<endl;
	}
   }
  out.close();  
  plotA->GetXaxis()->SetTitle("X (m)");
  plotA->GetYaxis()->SetTitle("Y (m)");
  //int bin32ft = plotA->GetBinContent(wt);
  /*
  TAxis* plotA=
   int Ax=plotA->GetNbinsX();
   int Ay=plotA->GetNbinsY();
   for (int ix=1; ix<Ax; ix++) {
      for (int iy=1; iy<Ay; iy++) {
         double wt_new=plotA->GetBinContent(ix,iy);
         if(wt_new
*/
}

