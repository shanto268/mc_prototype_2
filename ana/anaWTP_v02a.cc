//
//  c++  `root-config --cflags` -o anabin anaWTP_v02a.cc   `root-config --glibs`
//  ./anabin  inputFileName  outputFileName
//
//  for example,
//  ./a.out master_run009_new_10k.txt histout.root
//
//  then 
//  root -l histout.root
//
//  Input Data Format:  
// =========================
//     RawEvent is separated by 000000.
//     One line per worker node (31 words).
//     In each line, 3 words for PC header, 4 Words for Master header, 
//       4 words for Worker heaader,  4 words per channel for 6 channels.
//       3 + 4 + 4 + 6*4 = 32 words.
//  Ph 0  PC event number
//  Ph 1  PC date
//  Ph 2  PC time
//
//  Mh 0  Master evet number
//  Mh 1  Wireless addrss of Worker node
//  Mh 2  Master reset time (in Master clock time, i.e gloabl time)
//  Mh 3  222 for validation of PC-Master communication 
//
//  Wh 0  Trigger Bits (1 bit per channel)
//  Wh 1  Worker Local Time since Master reset (higher byte)
//  Wh 2  Worker Local Time since Master reset (lower byte)
//  Wh 3  Reset Time received from Master (in Master Clock Time)
//
//  Wch0 0  Channel ID    (layer*64+ breadboard*16 + ch)
//  Wch0 1  Hit time in Worker Local TIme
//  Wch0 2  ADC at TS0
//  Wch0 3  ADC at TS1
//  (repeat 4 words each for Wch1, Wch2, Wch3, Wch4, Wch5).

//  Offline Data Structure (in this Analysis Code)
// ====================================================
//  sc8Event
//     RawData
//       errorCode, pcEventNumber, pcDate, pcTime
//       masterEventNumber, nodes, kcount; check222
//       hits (time, masterTime,layer,bb,ch,adc[2],node,xyz[3])
//       HitClusters (vector of hits in coincidence in maultiple layers)
//     vector(MuonEvent>
//       tracks  (x,y,xs,ys,hits)
//
//  Flow of Program
// ===================
//   main
//      sc8in.readEvent(evt)    // read event
//      reco.recoHitClusters(evt);   // reconstruct hit clusters in time window.
//      reco.recoPoints(evt);        // reconstruct points in each layer (obsolete)
//      reco.recoMuonEvents(evt);    // reconstruct muon tracks for hit cluster
//
//  Example of Usage
// ==================
//
// c++  `root-config --cflags` anaSC8A13e.cc   `root-config --glibs`
// ./a.out v4a_run066_muon_10k_6nodes_lightOFF_newThreshold.txt  hist_test.root
//


#include <iostream>  // for cout
#include <fstream>   // for input/output files
#include <sstream>   // for string stream
#include <math.h>    // for sin(x) etc.
#include <cstdlib>   // for rand() on archer.
#include <iomanip>   // for setw() in cout, 
#include <vector>
#include <utility>

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

using namespace std;

bool  printFlag=false;

//   argv[1]    input data file
//   argv[2]    output histogram file

double theTHETA=50.0;   // argv[3]
double thePHI=0.0;      // argv[4]

int timeWindowCut=1;  //  argv[5] in sc8Reco::recoHitClusters
int clusterSizeCut=1;  //  argv[6] in sc8Reco::recoMuonEvents



struct WorkerOutput{
   int w[28];   //why 28?
};

// ============================================================
// ============================================================
class Hit{
  public:
    int time;         // worker local time
    int masterTime;   // master global time
    int layer;
    int bb;
    int ch;
    int adc[10];
    int node;
    vector<double> xyz();
};

// ============================================================
vector<double> Hit::xyz() {
   vector<double> zLayer={55.0,55.0,-55.0,-55.0};
   double d=5.0;  // scintillator spacing in x and y.
   double x=0.0;
   double y=0.0;
   double z=zLayer[layer];

   if(layer==0) {
     x=(bb*6.0+ch)*d;
   }

  if(layer==1) {
     y=(bb*6.0+ch)*d;
   }
   if(layer==2) {
     x=(bb*6.0+ch)*d;
   }
   if(layer==3) {
     y=(bb*6.0+ch)*d;
   }
   vector<double> _xyz={x,y,z};
   return _xyz;;
}


// ============================================================
struct HitCluster{
   int nlayersWithHits;
   vector<int> nhitsPerLayer;
   vector<Hit> hits;
};


// =========================================================================
// =========================================================================
class Points {
   public:
      vector<double> z;   // z for four layers

      //  double for x or y value, int for index of hit in the hit container.
      vector<pair<double,vector<int>>> x0;
      vector<pair<double,vector<int>>> y1;
      vector<pair<double,vector<int>>> x2;
      vector<pair<double,vector<int>>> y3;

      vector<pair<double,vector<int>>> getXY(int k);  // k 0,1,2,3 for x0,y1,x2,y2, respectively

};

// ============================================================
vector<pair<double,vector<int>>> Points::getXY(int k) {
    vector<pair<double,vector<int>>> a;
    vector<pair<double,vector<int>>> b;
    if(k==0) a=x0;
    if(k==1) a=y1;
    if(k==2) a=x2;
    if(k==3) a=y3;

    if(a.size()==1) b=a;

    if(a.size()>1) {
       // sort points by coordinates...
       auto sortRuleLambda = [](pair<double,vector<int>> p1, pair<double,vector<int>> p2) ->bool
       {
           return (p1.first < p2.first);
       };
       std::sort(a.begin(),a.end(),sortRuleLambda);

       int i=0;
       while (i<(a.size()-1)) {
          pair<double,vector<int>> q;
          q=a[i];
          double val=a[i].first;
          vector<int> hitsList=a[i].second;
          double n=1.0;
          int j=i+1;
          while (j<a.size()) {
             double d=a[j].first-a[j-1].first;
             if(d>5.5) break; 
             val=val+a[j].first;
             hitsList.push_back(a[j].second[0]);
             n=n+1.0;
             j++;
          }  // end of while 2.
          val=val/n;
          q.first=val;  // update the value.  (no update for hit, for now).
          q.second=hitsList;
          b.push_back(q);
          i=j;
       } // end of while 1.
/*
       cout<<" Points::getXY  k="<<k<<endl;
       for (int i=0; i<a.size();i++) {
         cout<<" i="<<i<<"   a_first "<<a[i].first<<"  a_second "<<a[i].second[0]<<endl;
       }
       for (int j=0; j<b.size(); j++) {
         cout<<" j="<<j<<"   b_first "<<b[j].first<<"  b_second "<<b[j].second[0]<<endl;
       }
*/
    } // end of if(a.size()>1)

    return b;
}



// =========================================================================
// =========================================================================
class Track {
   public:
      Track(double _x, double _y, double _xs, double _ys, vector<Hit> _hits ) {
        x=_x; y=_y; xs=_xs; ys=_ys; hits=_hits;
      };
      double x;    // at z=0
      double y;    // at z=0
      double xs;   // xslope = (x0-x2)/(z0-z2)
      double ys;   // yslope = (y1-y3)/(z1-z3)
      vector<Hit> hits;
};

// =========================================================================
class RawData{
   public:
      // status block
      vector<int> errorCode; // if empty, GOOD. 
                             // 222  PC-Master communication error
                             // 1-8  worker node number with bad data 
                             // 9    node missing.
      // data block
      vector<int> pcEventNumber;  
      vector<string> pcDate;
      vector<string> pcTime;
      vector<int>    masterEventNumber;
      vector<int>    node;
      vector<int>    kcount;
      vector<int>    check222;
      vector<Hit>    hits;
      vector<WorkerOutput> wouts;
      vector<HitCluster> hclusters;
      Points points;
};

// ==========================================================================
// ====class Track {
   public:
      Track(double _x, double _y, double _xs, double _ys, vector<Hit> _hits ) {
        x=_x; y=_y; xs=_xs; ys=_ys; hits=_hits;
      };
      double x;    // at z=0
      double y;    // at z=0
      double xs;   // xslope = (x0-x2)/(z0-z2)
      double ys;   // yslope = (y1-y3)/(z1-z3)
      vector<Hit> hits;
};
======================================================================
class MuonEvent{
   public:
      MuonEvent(){};
      void addMuonTrack(Track tk) {tracks.push_back(tk);}
      vector<Track> tracks;
};

// ==========================================================================
// ==========================================================================
class sc8Event{
   public:
      sc8Event(int);
      ~sc8Event();
      void addRawData(RawData);
      void addMuonEvent(MuonEvent mu) {_muonevent.push_back(mu);};
      RawData getRawData();
      vector<MuonEvent>  getMuonEvents() {return _muonevent;};
      void setHitClusters(vector<HitCluster> hcs) {_rawdata.hclusters=hcs;};
      void setPoints(Points p) {_rawdata.points=p;}
      void print();

   private:
      RawData _rawdata;
      vector<MuonEvent>  _muonevent;
};

// ==========================================================================
sc8Event::sc8Event(int nn) {

}

// ==========================================================================
sc8Event::~sc8Event() {

}
// ==========================================================================
void sc8Event::addRawData(RawData raw) {
   _rawdata=raw;
}
// ==========================================================================
RawData sc8Event::getRawData() {
   return _rawdata;
}
// ==========================================================================
void sc8Event::print() {
    //cout<<" ===== Event Print ====="<<endl;
    // cout<<"  size "<<_event.pcEventNumber.size()<<endl;
    cout<<"  "<<endl;;
    cout<<"Event: "<<_rawdata.pcEventNumber[0];
    cout<<"  "<<_rawdata.pcDate[0];
    cout<<"  "<<_rawdata.pcTime[0];
    cout<<"  node size "<<_rawdata.node.size()<<"   nodes=";
    for(int i=0; i<_rawdata.node.size(); i++) { 
      cout<<"  "<<_rawdata.node[i];
    }
    cout<<endl;
    // cout<<"number of hits "<<_event.hits.size();
    for(int i=0; i<_rawdata.hits.size(); i++) {
      cout<<" hit no. "<<setw(2)<<i;
      cout<<"    time "<<setw(4)<<_rawdata.hits[i].time;
      cout<<"    masterTime "<<setw(4)<<_rawdata.hits[i].masterTime;
      cout<<"   node "<<setw(2)<<_rawdata.hits[i].node;
      cout<<"   layer "<<setw(2)<<_rawdata.hits[i].layer;
      cout<<"   bb "<<setw(2)<<_rawdata.hits[i].bb;
      cout<<"   ch "<<setw(2)<<_rawdata.hits[i].ch;
      cout<<"   adc "<<setw(5)<<_rawdata.hits[i].adc[0];;
      cout<<setw(5)<<_rawdata.hits[i].adc[1];;
      cout<<endl;
    }
    cout<<endl;
   
    //  print hit clusters
    int nhc=_rawdata.hclusters.size();
    cout<<"hit cluster size= "<<nhc<<endl;
    for(int k=0; k<nhc; k++) {
       HitCluster hc=_rawdata.hclusters[k];
       int nhits=hc.hits.size();
       cout<<"     coincidence "<<hc.nlayersWithHits<<"   Layer  ";
       for(int i=0; i<4; i++) {cout<<"  "<<hc.nhitsPerLayer[i];}
       cout<<endl;
       for(int i=0; i<nhits; i++) {
         vector<double> xyz=hc.hits[i].xyz();
         cout<<" cluster no. "<<setw(2)<<k;
         cout<<" hit no. "<<setw(2)<<i;
         cout<<"    time "<<setw(4)<<hc.hits[i].time;
         cout<<"    masterTime "<<setw(4)<<hc.hits[i].masterTime;
         cout<<"   node "<<setw(2)<<hc.hits[i].node;
         cout<<"   layer "<<setw(2)<<hc.hits[i].layer;
         cout<<"   bb "<<setw(2)<<hc.hits[i].bb;
         cout<<"   ch "<<setw(2)<<hc.hits[i].ch;
         cout<<"   adc "<<setw(5)<<hc.hits[i].adc[0];;
         cout<<setw(5)<<hc.hits[i].adc[1];
         cout<<"    x "<<xyz[0];
         cout<<" y "<<xyz[1];
         cout<<" z "<<xyz[2];
         cout<<endl;
    }
    }
    cout<<endl;

    //  print out Points...
    Points p=_rawdata.points;
    cout<<"   ";
    cout<<"Points...  (x or y)"<<endl;

    cout<<"z0"<<setw(5)<<p.z[0]<<"   x0  ";
    for(int i=0; i<p.x0.size();i++) {cout<<"   "<<p.x0[i].first;}
    cout<<endl;

    cout<<"z2"<<setw(5)<<p.z[2]<<"   x2  ";
    for(int i=0; i<p.x2.size();i++) {cout<<"   "<<p.x2[i].first;}
    cout<<endl;

    cout<<"z1"<<setw(5)<<p.z[1]<<"   y1  ";
    for(int i=0; i<p.y1.size();i++) {cout<<"   "<<p.y1[i].first;}
    cout<<endl;
    cout<<"z3"<<setw(5)<<p.z[3]<<"   y3  ";
    for(int i=0; i<p.y3.size();i++) {cout<<"   "<<p.y3[i].first;}
    cout<<endl;
    
    //   Muon Events...
    vector<MuonEvent> muEvents=getMuonEvents();
    cout<<"   "<<endl;
    cout<<"Number of MuonEvents:  "<<muEvents.size()<<endl;
    for(int i=0; i<muEvents.size(); i++) {
       cout<<"  muonEvent "<<i<<"  tracks "<<muEvents[i].tracks.size()<<endl;
       for (int t=0; t<muEvents[i].tracks.size(); t++) {
          Track tk=muEvents[i].tracks[t];
          cout<<"     track "<<t<<setw(5)<<"    x "<<setw(5)<<tk.x<<"   y "<<setw(5)<<tk.y;
          cout<<setw(5)<<"  xs "<<setw(5)<<tk.xs<<"  ys "<<setw(5)<<tk.ys<<endl;
       }
    } 


/*
    cout<<"   ";
    for (int j=0; j<_event.wouts.size(); j++) {
      cout<<" worker output 19-27:  ";
      for (int i=19; i<28; i++) {
         cout<< "  "<<setw(5)<<_event.wouts[j].w[i];
      }
      cout<<endl;
    }
*/

}



// ==========================================================================
// ==========================================================================
class sc8Input{
   public:
      sc8Input(string);
      int readEvent(sc8Event &);

   private:
      ifstream infile;
      int linecount;
};

// ==========================================================================
sc8Input::sc8Input(string inname) {
   infile.open(inname);
   cout<<"sc8Input: input file="<<inname<<endl;
   linecount=0;
}

// ==========================================================================
int sc8Input::readEvent(sc8Event & evt){
   int status=0;
   if(infile.is_open()) {
      std::map<int,int>  wirelessAddress2nodeNumber;
      // updated this mapping.   sk 01-July-2019
      wirelessAddress2nodeNumber[9]=3;
      wirelessAddress2nodeNumber[17]=2;
      wirelessAddress2nodeNumber[25]=6;
      wirelessAddress2nodeNumber[10]=7;
      wirelessAddress2nodeNumber[18]=1;
      wirelessAddress2nodeNumber[26]=5;
      wirelessAddress2nodeNumber[33]=4;
      wirelessAddress2nodeNumber[34]=8;

      string line;
      bool insideofevent=false;
      int ncounts=0;
      RawData e;
      int t0=-100;
      while(getline(infile,line)) {
         linecount++;
         // cout<<"line="<<line<<endl;
         // if(insideofevent) cout<<linecount<<"      ="<<line<<endl;
         if(line.size()==0) {continue;}

         std::istringstream iss(line);
         int pcevent;
         iss>>pcevent;
         if(pcevent==0) {
            t0=-100;
            break;   // end of event mark
         }
         ncounts=ncounts+1;
         e.pcEventNumber.push_back(pcevent);
         string s;
         iss>>s;
         e.pcDate.push_back(s);
         iss>>s;
         e.pcTime.push_back(s);
         // header from Master...
         int a[4];
         for (int i=0; i<4; i++) {iss>>a[i];}
         e.masterEventNumber.push_back(a[0]);
         // int node=a[1];
         // convert worker address to ordinal node numbers
         int node=wirelessAddress2nodeNumber[a[1]];
         e.node.push_back(node);
         e.kcount.push_back(a[2]);
         e.check222.push_back(a[3]);
         if(a[3]!=222) {
             e.errorCode.push_back(222);break;  //  check 222 for the data validity.
         }
         status=1;
         // data from worker...
         int b[28];
         WorkerOutput w;
         for (int i=0; i<28; i++) {
           iss>>b[i];
           w.w[i]=b[i];
         }
         if(iss.fail()) {
            e.errorCode.push_back(node); continue; //less than 28 words for worker data. 
         }
         if((b[0]+b[1]+b[2]+b[3])==0) {
            e.errorCode.push_back(node); continue; // header words all zero.
         }
         e.wouts.push_back(w);

         int mytime=b[1]*256+b[2];
         if(t0<0) t0=b[3];
         int nhits=6;
         if(nhits>0) {
           int k=4;
           // vector<int> layerMap;
           // layerMap.assign({0,1,2,2,3,0,0}); // node to layer number
           // vector<int> bbMap;
           // bbMap.assign({0,0,0,1,0,0,1});
           for (int i=0; i<nhits; i++) {
              Hit h;
              h.layer=b[k]/64;
              // h.layer=layerMap[node];
              h.bb=(b[k]%64)/16;
              // h.bb=bbMap[node];
              h.ch=b[k]%16;
              h.time=b[k+1];
              int dt=b[3]-t0;
              if(dt<0) dt=dt+256;
              int tglobal=h.time+dt;;
              h.masterTime=tglobal;
              h.node=node;
              h.adc[0]=b[k+2];
              h.adc[1]=b[k+3];
              if(h.time>0) {
                e.hits.push_back(h);
              }
              k=k+4;
           } // end of loop over hits.
         } // end of if(nhits>0)
      }
      // sort hits by mastertime
      auto sortRuleLambda = [](Hit& p1, Hit& p2) ->bool
      {
         return (p1.masterTime < p2.masterTime);
      };
      std::sort(e.hits.begin(),e.hits.end(),sortRuleLambda);

      // check the number of nodes in the data.
      if(e.node.size() != 8) {e.errorCode.push_back(9);}  

      evt.addRawData(e);   //  add Event...
   }
   if(infile.eof()) status=2;

   return status; // 0=non event data, 1=event data, 2=eof
}

// ==========================================================================
// ==========================================================================
class sc8Reco{
   public:
      sc8Reco(){};

      void recoAll(sc8Event&);
      void recoHitClusters(sc8Event&);  // hit clustering in time space...
      void recoPoints(sc8Event&);
      void recoMuonEvents(sc8Event&);

};

// ==========================================================================
void sc8Reco::recoAll(sc8Event& evt) {
   recoHitClusters(evt);
   recoPoints(evt);
   recoMuonEvents(evt);
}
// ==========================================================================
void sc8Reco::recoHitClusters(sc8Event& evt) {

   // fout->cd();
   RawData e=evt.getRawData();
  
   vector<HitCluster> hcs;

   int nhits=e.hits.size();
   int a=0;               
   bool mybreak=true;
   while (a<(nhits-1)) {
     HitCluster hc;
     int b=a+1;
     while (b<nhits) {
        int dt=e.hits[b].masterTime-e.hits[a].masterTime;
        // cout<<" a  "<<a<<"   b  "<<b<<"   dt "<<dt<<endl;
        if(dt>timeWindowCut)  {mybreak=true; break;}   //Change from dt>2 to dt>1
        if(mybreak) {
              hc.nlayersWithHits=0;
              hc.nhitsPerLayer={0,0,0,0};
              hc.nhitsPerLayer[e.hits[a].layer]=hc.nhitsPerLayer[e.hits[a].layer]+1;
              hc.hits.push_back(e.hits[a]);
              mybreak=false;
              hcs.push_back(hc);
         }
         int ix=hcs.size()-1;
         hcs[ix].nhitsPerLayer[e.hits[b].layer]=hcs[ix].nhitsPerLayer[e.hits[b].layer]+1;
         hcs[ix].hits.push_back(e.hits[b]);
         b++;
     }  // end of while (b)
     a=b;
   } // end of while(a)

   // cout<<"  hcs.size() "<<hcs.size()<<endl;
   for(int i=0; i<hcs.size();i++) {
     for(int j=0; j<4; j++) {
        if(hcs[i].nhitsPerLayer[j]>0) {
          hcs[i].nlayersWithHits=hcs[i].nlayersWithHits+1;
        }
     }
     //  cout<<"hcs[i].nlayersWithHits   "<<hcs[i].nlayersWithHits<<endl;
   }

   evt.setHitClusters(hcs);

   // cout<<"sc8Analyzer::recoHitClusters e.hclusters.size()="<<e.hclusters.size()<<endl;

   return;
}

// ==========================================================================
void sc8Reco::recoPoints(sc8Event& evt) {

    Points p;
    p.z={0.0, 0.0, 0.0, 0.0};

    RawData _rawdata=evt.getRawData();
    int nhc=_rawdata.hclusters.size();
    for(int k=0; k<nhc; k++) {
       HitCluster hc=_rawdata.hclusters[k];
       int nhits=hc.hits.size();
       for(int i=0; i<nhits; i++) {
         vector<double> xyz=hc.hits[i].xyz();
         int layer=hc.hits[i].layer;
         std::pair<double,vector<int>> a;
         vector<int> hitsList; hitsList.push_back(i);
         if(layer==0) {a=std::make_pair(xyz[0],hitsList);p.x0.push_back(a); p.z[0]=xyz[2];}
         if(layer==1) {a=std::make_pair(xyz[1],hitsList);p.y1.push_back(a); p.z[1]=xyz[2];}
         if(layer==2) {a=std::make_pair(xyz[0],hitsList);p.x2.push_back(a); p.z[2]=xyz[2];}
         if(layer==3) {a=std::make_pair(xyz[1],hitsList);p.y3.push_back(a); p.z[3]=xyz[2];}
       }
    }

    evt.setPoints(p);
    return;
}

// ==========================================================================
void sc8Reco::recoMuonEvents(sc8Event& evt) {

    RawData _rawdata=evt.getRawData();
    //  loop over hit clusters...
    int nhc=_rawdata.hclusters.size();
    // cout<<"sc8Reco::recoMuonEvents:   nhc="<<nhc<<endl;
    for(int k=0; k<nhc; k++) {
       HitCluster hc=_rawdata.hclusters[k];
       int nhits=hc.hits.size();
       Points p;
       p.z={0.0, 0.0, 0.0, 0.0};
       // cout<<"       nhits  "<<nhits<<endl;
       for(int i=0; i<nhits; i++) {
         vector<double> xyz=hc.hits[i].xyz();
         int layer=hc.hits[i].layer;
         std::pair<double,vector<int>> a;
         vector<int> hitsList; hitsList.push_back(i);
         if(layer==0) {a=std::make_pair(xyz[0],hitsList);p.x0.push_back(a); p.z[0]=xyz[2];}
         if(layer==1) {a=std::make_pair(xyz[1],hitsList);p.y1.push_back(a); p.z[1]=xyz[2];}
         if(layer==2) {a=std::make_pair(xyz[0],hitsList);p.x2.push_back(a); p.z[2]=xyz[2];}
         if(layer==3) {a=std::make_pair(xyz[1],hitsList);p.y3.push_back(a); p.z[3]=xyz[2];}
       }  // end of loop over hits

       // require hits in all four layers..
       // cout<<"sizes: x0="<<p.x0.size()<<"  y1="<<p.y1.size();
       // cout<<"  x2="<<p.x2.size()<<"  y3="<<p.y3.size()<<endl;
       if(p.x0.size()==0 || p.y1.size()==0 || p.x2.size()==0 || p.y3.size()== 0) continue;

       // cout<<"   after continue in MuonEvents..."<<endl;
       vector<Hit> hits;
       MuonEvent muonevent;
       vector<std::pair<double,vector<int>>> qx0=p.getXY(0);  // sorted in x or y values.
       vector<std::pair<double,vector<int>>> qy1=p.getXY(1);
       vector<std::pair<double,vector<int>>> qx2=p.getXY(2);
       vector<std::pair<double,vector<int>>> qy3=p.getXY(3);
       for(int ix0=0; ix0<qx0.size(); ix0++) {
          if(qx0[ix0].second.size()>clusterSizeCut) continue;
          // cout<<"ix0="<<ix0<<"   second "<<p.x0[ix0].second<<endl;
          for (int k=0; k< qx0[ix0].second.size(); k++) {
             hits.push_back(hc.hits[qx0[ix0].second[k]]);
          }
       for(int ix2=0; ix2<qx2.size(); ix2++) {
          if(qx2[ix2].second.size()>clusterSizeCut) continue;
          // cout<<"ix2="<<ix2<<"   second "<<p.x2[ix2].second<<endl;
          double xs=(qx0[ix0].first-qx2[ix2].first)/(p.z[0]-p.z[2]);          
          for (int k=0; k< qx2[ix2].second.size(); k++) {
             hits.push_back(hc.hits[qx2[ix2].second[k]]);
          }

          for(int iy1=0; iy1<qy1.size(); iy1++) {
             if(qy1[iy1].second.size()>clusterSizeCut) continue;
          // cout<<"iy1="<<iy1<<"   second "<<p.y1[iy1].second<<endl;
             for (int k=0; k< qy1[iy1].second.size(); k++) {
                hits.push_back(hc.hits[qy1[iy1].second[k]]);
             }
          for(int iy3=0; iy3<qy3.size(); iy3++) {
             if(qy3[iy3].second.size()>clusterSizeCut) continue;
          // cout<<"iy3="<<iy3<<"   second "<<p.y3[iy3].second<<endl;
             for (int k=0; k< qy3[iy3].second.size(); k++) {
                hits.push_back(hc.hits[qy3[iy3].second[k]]);
             }
             double ys=(qy1[iy1].first-qy3[iy3].first)/(p.z[1]-p.z[3]);
             
             double x=qx0[ix0].first-xs*p.z[0];
             double y=qy1[iy1].first-ys*p.z[1];
             Track tk(x,y,xs,ys,hits);
             muonevent.addMuonTrack(tk);
          }  // end of loop over y3
          }  // end of loop over y1

       } // end of loop over x2
       } // end of loop over x0

       if(muonevent.tracks.size()>0) {
         // cout<<"  recoMuonEvent...   add muon event "<<endl;
         evt.addMuonEvent(muonevent);
       }

    } // end of loop over hitClusters

    return;


}

// ==========================================================================
// ==========================================================================
class HistSet1{
   public:
      HistSet1() {};
      void init(TFile*,string);
      void analyze(sc8Event&);
      void endjob();

   private:
      std::map<std::string, TH1D*> histo1D;
      std::map<std::string, TH1D*>::iterator histo1Diter;

      std::map<std::string, TH2D*> histo2D;
      std::map<std::string, TH2D*>::iterator histo2Diter;

      TFile *fout;
      string dirName;
      TDirectory *tdir;
};

// ==========================================================================
void HistSet1::init(TFile* _fout, string _dirName) {
   fout=_fout;
   string dirName=_dirName;
   tdir=fout->mkdir(dirName.c_str());
   tdir->cd();

   histo1D["nNodesAll"]=new TH1D("nNodesAll","numbe of nodes (All events)",10,0.,10.);
   histo1D["nNodesAccepted"]=new TH1D("nNodesAccepted","numbe of nodes (Accepted Events)",10,0.,10.);
   histo1D["NodeNumber"]=new TH1D("NodeNumber","Node Number",10,0.,10.);
   histo1D["N1Time"]=new TH1D("N1Time","Node 1- timing",500,0.,500.);
   histo1D["N2Time"]=new TH1D("N2Time","Node 2- timing",500,0.,500.);
   histo1D["N3Time"]=new TH1D("N3Time","Node 3- timing",500,0.,500.);
   histo1D["N4Time"]=new TH1D("N4Time","Node 4- timing",500,0.,500.);
   histo1D["N5Time"]=new TH1D("N5Time","Node 5- timing",500,0.,500.);
   histo1D["N6Time"]=new TH1D("N6Time","Node 6- timing",500,0.,500.);
   histo1D["N7Time"]=new TH1D("N7Time","Node 7- timing",500,0.,500.);
   histo1D["N8Time"]=new TH1D("N8Time","Node 8- timing",500,0.,500.);
   histo1D["N1MasterTime"]=new TH1D("N1MasterTime","Node 1- Master Time",500,0.,500.);
   histo1D["N2MasterTime"]=new TH1D("N2MasterTime","Node 2- Master Time",500,0.,500.);
   histo1D["N3MasterTime"]=new TH1D("N3MasterTime","Node 3- Master Time",500,0.,500.);
   histo1D["N4MasterTime"]=new TH1D("N4MasterTime","Node 4- Master Time",500,0.,500.);
   histo1D["N5MasterTime"]=new TH1D("N5MasterTime","Node 5- Master Time",500,0.,500.);
   histo1D["N6MasterTime"]=new TH1D("N6MasterTime","Node 6- Master Time",500,0.,500.);
   histo1D["N7MasterTime"]=new TH1D("N7MasterTime","Node 7- Master Time",500,0.,500.);
   histo1D["N8MasterTime"]=new TH1D("N8MasterTime","Node 8- Master Time",500,0.,500.);

   histo1D["nHits"]=new TH1D("nHits","numbe of hits",40,0.,40.);
   histo1D["LayerNum"]=new TH1D("LayerNum","Layer Number",5,0.,5.);
   histo1D["L0Ch"]=new TH1D("L0Ch","L0 Channel Numbers",15,0.,15.);
   histo1D["L1Ch"]=new TH1D("L1Ch","L1 Channel Numbers",15,0.,15.);
   histo1D["L2Ch"]=new TH1D("L2Ch","L2 Channel Numbers",15,0.,15.);
   histo1D["L3Ch"]=new TH1D("L3Ch","L3 Channel Numbers",15,0.,15.);

   histo1D["adcL0ChAll"]=new TH1D("adcL0ChAll","L0 ADC ped sub(all channels)",300,0.,300.);
   histo1D["adcL1ChAll"]=new TH1D("adcL1ChAll","L1 ADC ped sub(all channels)",300,0.,300.);
   histo1D["adcL2ChAll"]=new TH1D("adcL2ChAll","L2 ADC ped sub (all channels)",300,0.,300.);
   histo1D["adcL3ChAll"]=new TH1D("adcL3ChAll","L3 ADC ped sub (all channels)",300,0.,300.);


   string hname;
   string htitle;
   for (int i=0; i<5; i++) {
     for (int j=0; j<13; j++) {
       hname="Layer"+to_string(i)+"_Ch"+to_string(j)+"_adc0";
       htitle="Layer "+to_string(i)+" Ch "+to_string(j)+" ADC[0]";
       histo1D[hname]=new TH1D(hname.c_str(),htitle.c_str(),300,0.,300.);
       hname="Layer"+to_string(i)+"_Ch"+to_string(j)+"_adc1";
       htitle="Layer "+to_string(i)+" Ch "+to_string(j)+" ADC[1]";
       histo1D[hname]=new TH1D(hname.c_str(),htitle.c_str(),300,0.,300.);
       hname="Layer"+to_string(i)+"_Ch"+to_string(j)+"_adcPedsub";
       htitle="Layer "+to_string(i)+" Ch "+to_string(j)+" ADC(ped subtracted)";
       histo1D[hname]=new TH1D(hname.c_str(),htitle.c_str(),300,0.,300.);
     }  // end of loop for j (channel number)
   }  // end of loop over i (layers)
   
   //  hit clustering...
   histo1D["dt"]=new TH1D("dt","hitB-hitA",200,0.,200.);
   histo1D["dt2"]=new TH1D("dt2","hitB-hitA",20,0.,20.);

   // x projection and y projection...
   histo1D["X0a"]=new TH1D("X0a","X (layer 0) x-coincidence",120,0.,60.);
   histo1D["Y1a"]=new TH1D("Y1a","Y (layer 1) y-coincidence",120,0.,60.);
   histo1D["X0b"]=new TH1D("X0b","X (layer 0) xy-coincidence",120,0.,60.);
   histo1D["Y1b"]=new TH1D("Y1b","Y (layer 1) xy-coincidence",120,0.,60.);

   // muon analysis
   histo1D["nMuonEvents"]=new TH1D("nMuonEvents","number of muonEvents",5,0.,5.);
   histo1D["nMuons"]=new TH1D("nMuons","number of muons",10,0.,10.);
   histo1D["muX"]=new TH1D("muX","muon: X at Z=0",60,0.0,60.0);   
   histo1D["muY"]=new TH1D("muY","muon: Y at Z=0",60,0.0,50.0);
   histo1D["muXS"]=new TH1D("muXS","muon: X slope (theta x-proj) (deg)",90,-45.0,45.0);
   histo1D["muYS"]=new TH1D("muYS","muon: Y slope (theta y-proj) (deg)",90,-45.0,45.0);
   histo1D["muTheta"]=new TH1D("muTheta","muon:  theta (deg)",90.0,0.0,45.0);

}



// ==========================================================================
void HistSet1::analyze(sc8Event& evt) {
   tdir->cd();
   RawData e=evt.getRawData();
   histo1D["nNodesAll"]->Fill(e.node.size());

   histo1D["NodeNumber"]->Fill(0);   // to count all events
   for (int i=0; i<e.node.size(); i++) {
     histo1D["NodeNumber"]->Fill(e.node[i]);
   }

   // check the event status...
   if(e.errorCode.size()>0) return;  // skip event with error...
   histo1D["nNodesAccepted"]->Fill(e.node.size());

   int nhits=e.hits.size();
   histo1D["nHits"]->Fill(nhits);

   for (int i=0; i<nhits; i++) {
      int layer=e.hits[i].layer;
      histo1D["LayerNum"]->Fill(layer);

      string hname="L"+to_string(layer)+"Ch";
      int ch=e.hits[i].bb*6+e.hits[i].ch;
      histo1D[hname]->Fill(ch);     

      hname="Layer"+to_string(layer)+"_Ch"+to_string(ch)+"_adc0";
      histo1D[hname]->Fill(e.hits[i].adc[0]);

      hname="Layer"+to_string(layer)+"_Ch"+to_string(ch)+"_adc1";
      histo1D[hname]->Fill(e.hits[i].adc[1]);

      hname="Layer"+to_string(layer)+"_Ch"+to_string(ch)+"_adcPedsub";
      int adcPedSubtracted=e.hits[i].adc[1]-e.hits[i].adc[0];
      histo1D[hname]->Fill(adcPedSubtracted);
      
      hname="adcL"+to_string(layer)+"ChAll";
      histo1D[hname]->Fill(adcPedSubtracted);

     hname="N"+to_string(e.hits[i].node)+"Time";
     histo1D[hname]->Fill(e.hits[i].time);

     hname="N"+to_string(e.hits[i].node)+"MasterTime";
     histo1D[hname]->Fill(double(e.hits[i].masterTime));
     // if(e.hits[i].masterTime>300) {
     //   cout<<"MaterTimeError "<<e.hits[i].masterTime<<"   i="<<i<<endl;
     //   skDump=1;
     // }

  }  // end of for (int i=0; i<nhits; i++)

   //  delta time for hit clustering...
   for (int i=0; i<(nhits-1); i++) {
      int j=i+1;
      double dt=e.hits[j].masterTime-e.hits[i].masterTime;
      histo1D["dt"]->Fill(dt);
      histo1D["dt2"]->Fill(dt);
   } //  end of for(int i...

  //   Points 
  Points p=e.points;
  if(p.x0.size()==1 && p.x2.size()==1) {
     histo1D["X0a"]->Fill(p.x0[0].first);
  }
  if(p.y1.size()==1 && p.y3.size()==1) {
     histo1D["Y1a"]->Fill(p.y1[0].first);
  }

  if(p.x0.size()==1 && p.x2.size()==1
     && p.y1.size()==1 && p.y3.size()==1) {
     histo1D["X0b"]->Fill(p.x0[0].first);
     histo1D["Y1b"]->Fill(p.y1[0].first);
  }

  //   Muon tracks...
  vector<MuonEvent> muEvents=evt.getMuonEvents();
  int nMuonEvents=muEvents.size();
  histo1D["nMuonEvents"]->Fill(nMuonEvents);
  for (int i=0; i<nMuonEvents; i++) {
     int ntracks=muEvents[i].tracks.size();
     histo1D["nMuons"]->Fill(ntracks);
       
     for (int t=0; t<muEvents[i].tracks.size(); t++) {
          Track tk=muEvents[i].tracks[t];
          histo1D["muX"]->Fill(tk.x);
          histo1D["muY"]->Fill(tk.y);
          double rad2deg=180.0/3.1415;
          histo1D["muXS"]->Fill(tk.xs*rad2deg);
          histo1D["muYS"]->Fill(tk.ys*rad2deg);
          double theta=atan(sqrt(tk.xs*tk.xs+tk.ys*tk.ys))*rad2deg;
          histo1D["muTheta"]->Fill(theta);
     }  // end of loop over muon tracks.
  }  //  end of loop over muon events


}

// ==========================================================================
void HistSet1::endjob() {

}


// ==========================================================================
// ==========================================================================
class HistHitClusters{
   public:
      HistHitClusters() {};
      void init(TFile*,string);
      void analyze(sc8Event&);
      void endjob(){}

      void fillHits(Track&);;

   private:
      std::map<std::string, TH1D*> histo1D;
      std::map<std::string, TH1D*>::iterator histo1Diter;

      std::map<std::string, TH2D*> histo2D;
      std::map<std::string, TH2D*>::iterator histo2Diter;

      TFile *fout;
      string dirName;
      TDirectory *tdir;
};

// ==========================================================================

// ==========================================================================
void HistHitClusters::init(TFile* _fout, string _dirName) {
   fout=_fout;
   string dirName=_dirName;
   tdir=fout->mkdir(dirName.c_str());
   tdir->cd();

   histo1D["nClusters"]=new TH1D("nClusters","N of Clusters",20,0.0,20.);
   histo1D["nHitsPerCluster"]=new TH1D("nHitsPerCluster","N of hits per cluster",20,0.0,20.);
   histo1D["hcL0Chs"]=new TH1D("hcL0Chs","hit cluster, ch in Layer 0",15,0.,15.0);
   histo1D["hcL1Chs"]=new TH1D("hcL1Chs","hit cluster, ch in Layer 1",15,0.,15.0);
   histo1D["hcL2Chs"]=new TH1D("hcL2Chs","hit cluster, ch in Layer 2",15,0.,15.0);
   histo1D["hcL3Chs"]=new TH1D("hcL3Chs","hit cluster, ch in Layer 3",15,0.,15.0);
   histo1D["hcL0nHits"]=new TH1D("hcL0nHits","N hits in cluster,  Layer 0",10,0.,10.0);
   histo1D["hcL1nHits"]=new TH1D("hcL1nHits","N hits in cluster,  Layer 0",10,0.,10.0);
   histo1D["hcL2nHits"]=new TH1D("hcL2nHits","N hits in cluster,  Layer 0",10,0.,10.0);
   histo1D["hcL3nHits"]=new TH1D("hcL3nHits","N hits in cluster,  Layer 0",10,0.,10.0);
}

void HistHitClusters::analyze(sc8Event& evt) {
//
//  th1:  detector tilt angle, measured as elevation angle
//  th2:  track slope measured in the detector coordinate system.  dx/dz=tan(th2)
//  th3:  track slope measured in the lab corrdinate system.  DX/DZ=tan(th3)
//           th3=th2-th1
//  X of the track at Zlevel(DZ)=50 (meters) is X=x(local)+DX=x(local)+DZ*tan(th3)
//  similar to Y.

   tdir->cd();
    RawData _rawdata=evt.getRawData();
    //  loop over hit clusters...
    int nhc=_rawdata.hclusters.size();
    // cout<<"sc8Reco::recoMuonEvents:   nhc="<<nhc<<endl;
    histo1D["nClusters"]->Fill(nhc);
    for(int k=0; k<nhc; k++) {
       HitCluster hc=_rawdata.hclusters[k];
       int nhits=hc.hits.size();
       histo1D["nHitsPerCluster"]->Fill(nhits);
       Points p;
       p.z={0.0, 0.0, 0.0, 0.0};
       int n[4]={0,0,0,0} ;
       // cout<<"       nhits  "<<nhits<<endl;
       for(int i=0; i<nhits; i++) {
          int layer=hc.hits[i].layer;
          int ch=hc.hits[i].ch+hc.hits[i].bb*6;
          // if(layer==3 && ch==4) printFlag=true;
          if(layer==0) {histo1D["hcL0Chs"]->Fill(ch); n[0]=n[0]+1;}
          if(layer==1) {histo1D["hcL1Chs"]->Fill(ch); n[1]=n[1]+1;}
          if(layer==2) {histo1D["hcL2Chs"]->Fill(ch); n[2]=n[2]+1;}
          if(layer==3) {histo1D["hcL3Chs"]->Fill(ch); n[3]=n[3]+1;}
       }  // end of loop over indivisual hits
       
       histo1D["hcL0nHits"]->Fill(n[0]);
       histo1D["hcL1nHits"]->Fill(n[1]);
       histo1D["hcL2nHits"]->Fill(n[2]);
       histo1D["hcL3nHits"]->Fill(n[3]);

    }  // end of loop over hit clusters.

}





// ==========================================================================
// ==========================================================================
class HistMuons{
   public:
      HistMuons() {};
      void init(TFile*,string);
      void analyze(sc8Event&);
      void endjob(){}

      void fillHits(Track&);;

   private:
      std::map<std::string, TH1D*> histo1D;
      std::map<std::string, TH1D*>::iterator histo1Diter;

      std::map<std::string, TH2D*> histo2D;
      std::map<std::string, TH2D*>::iterator histo2Diter;

      TFile *fout;
      string dirName;
      TDirectory *tdir;
};

// ==========================================================================
void HistMuons::init(TFile* _fout, string _dirName) {
   fout=_fout;
   string dirName=_dirName;
   tdir=fout->mkdir(dirName.c_str());
   tdir->cd();

   // muons...
   histo1D["nMuonEvents"]=new TH1D("nMuonEvents","number of muonEvents",5,0.,5.);
   histo1D["nMuons"]=new TH1D("nMuons","number of muons",10,0.,10.);
   histo1D["muX"]=new TH1D("muX","muon: X at Z=0",60,0.0,60.0);
   histo1D["muY"]=new TH1D("muY","muon: Y at Z=0",60,0.0,50.0);
   histo1D["muXS"]=new TH1D("muXS","muon: X slope (theta x-proj) (deg)",90,-45.0,45.0);
   histo1D["muYS"]=new TH1D("muYS","muon: Y slope (theta y-proj) (deg)",90,-45.0,45.0);
   histo1D["muTheta"]=new TH1D("muTheta","muon:  theta (deg)",90.0,0.0,45.0);

   //  hits on muon track...
   histo1D["N1MasterTime"]=new TH1D("N1MasterTime","Node 1- Master Time",500,0.,500.);
   histo1D["N2MasterTime"]=new TH1D("N2MasterTime","Node 2- Master Time",500,0.,500.);
   histo1D["N3MasterTime"]=new TH1D("N3MasterTime","Node 3- Master Time",500,0.,500.);
   histo1D["N4MasterTime"]=new TH1D("N4MasterTime","Node 4- Master Time",500,0.,500.);
   histo1D["N5MasterTime"]=new TH1D("N5MasterTime","Node 5- Master Time",500,0.,500.);
   histo1D["N6MasterTime"]=new TH1D("N6MasterTime","Node 6- Master Time",500,0.,500.);
   histo1D["N7MasterTime"]=new TH1D("N7MasterTime","Node 7- Master Time",500,0.,500.);
   histo1D["N8MasterTime"]=new TH1D("N8MasterTime","Node 8- Master Time",500,0.,500.);

   histo1D["nHits"]=new TH1D("nHits","number of hits",40,0.,40.);
   histo1D["LayerNum"]=new TH1D("LayerNum","Layer Number",5,0.,5.);
   histo1D["L0Ch"]=new TH1D("L0Ch","L0 Channel Numbers",15,0.,15.);
   histo1D["L1Ch"]=new TH1D("L1Ch","L1 Channel Numbers",15,0.,15.);
   histo1D["L2Ch"]=new TH1D("L2Ch","L2 Channel Numbers",15,0.,15.);
   histo1D["L3Ch"]=new TH1D("L3Ch","L3 Channel Numbers",15,0.,15.);

   histo1D["adcL0ChAll"]=new TH1D("adcL0ChAll","L0 ADC ped sub(all channels)",300,0.,300.);
   histo1D["adcL1ChAll"]=new TH1D("adcL1ChAll","L1 ADC ped sub(all channels)",300,0.,300.);
   histo1D["adcL2ChAll"]=new TH1D("adcL2ChAll","L2 ADC ped sub (all channels)",300,0.,300.);
   histo1D["adcL3ChAll"]=new TH1D("adcL3ChAll","L3 ADC ped sub (all channels)",300,0.,300.);

   histo1D["dtMuonL0"]=new TH1D("dtMuonL0","dt Muon in Layer 0",20,-10.,10.);
   histo1D["dtMuonL1"]=new TH1D("dtMuonL1","dt Muon in Layer 1",20,-10.,10.);
   histo1D["dtMuonL2"]=new TH1D("dtMuonL2","dt Muon in Layer 2",20,-10.,10.);
   histo1D["dtMuonL3"]=new TH1D("dtMuonL3","dt Muon in Layer 3",20,-10.,10.);


   string hname;
   string htitle;
   for (int i=0; i<5; i++) {
     for (int j=0; j<13; j++) {
       hname="Layer"+to_string(i)+"_Ch"+to_string(j)+"_adc0";
       htitle="Layer "+to_string(i)+" Ch "+to_string(j)+" ADC[0]";
       histo1D[hname]=new TH1D(hname.c_str(),htitle.c_str(),300,0.,300.);
       hname="Layer"+to_string(i)+"_Ch"+to_string(j)+"_adc1";
       htitle="Layer "+to_string(i)+" Ch "+to_string(j)+" ADC[1]";
       histo1D[hname]=new TH1D(hname.c_str(),htitle.c_str(),300,0.,300.);
       hname="Layer"+to_string(i)+"_Ch"+to_string(j)+"_adcPedsub";
       htitle="Layer "+to_string(i)+" Ch "+to_string(j)+" ADC(ped subtracted)";
       histo1D[hname]=new TH1D(hname.c_str(),htitle.c_str(),300,0.,300.);
     }  // end of loop for j (channel number)
   }  // end of loop over i (layers)

}  // end of HistMuons...

// ==========================================================================
void HistMuons::analyze(sc8Event& evt) {
   tdir->cd();

  //   Muon tracks...
  vector<MuonEvent> muEvents=evt.getMuonEvents();
  int nMuonEvents=muEvents.size();
  histo1D["nMuonEvents"]->Fill(nMuonEvents);
  for (int i=0; i<nMuonEvents; i++) {
     int ntracks=muEvents[i].tracks.size();
     // if(ntracks != 1) continue;
     if(ntracks >10) continue;
     histo1D["nMuons"]->Fill(ntracks);

     for (int t=0; t<muEvents[i].tracks.size(); t++) {
          Track tk=muEvents[i].tracks[t];
          int nhits=tk.hits.size();
          histo1D["muX"]->Fill(tk.x);
          histo1D["muY"]->Fill(tk.y);
          double rad2deg=180.0/3.1415;
          histo1D["muXS"]->Fill(tk.xs*rad2deg);
          histo1D["muYS"]->Fill(tk.ys*rad2deg);
          double theta=atan(sqrt(tk.xs*tk.xs+tk.ys*tk.ys))*rad2deg;
          histo1D["muTheta"]->Fill(theta);

          fillHits(tk);

     }  // end of loop over muon tracks.
  }  //  end of loop over muon events


}  // end of HistMuons::analyze

// ==========================================================================
void HistMuons::fillHits(Track& e) {
  
   int nhits=e.hits.size();
   histo1D["nHits"]->Fill(nhits);

   int t0=-1000;
   for (int i=0; i<nhits; i++) {
      int layer=e.hits[i].layer;
      histo1D["LayerNum"]->Fill(layer);

      string hname="L"+to_string(layer)+"Ch";
      int ch=e.hits[i].bb*6+e.hits[i].ch;
      histo1D[hname]->Fill(ch);

      hname="Layer"+to_string(layer)+"_Ch"+to_string(ch)+"_adc0";
      histo1D[hname]->Fill(e.hits[i].adc[0]);

      hname="Layer"+to_string(layer)+"_Ch"+to_string(ch)+"_adc1";
      histo1D[hname]->Fill(e.hits[i].adc[1]);

      hname="Layer"+to_string(layer)+"_Ch"+to_string(ch)+"_adcPedsub";
      int adcPedSubtracted=e.hits[i].adc[1]-e.hits[i].adc[0];
      histo1D[hname]->Fill(adcPedSubtracted);

      hname="adcL"+to_string(layer)+"ChAll";
      histo1D[hname]->Fill(adcPedSubtracted);

     // hname="N"+to_string(e.hits[i].node)+"Time";
     // histo1D[hname]->Fill(e.hits[i].time);

     hname="N"+to_string(e.hits[i].node)+"MasterTime";
     histo1D[hname]->Fill(double(e.hits[i].masterTime));

     if(t0<-999) t0=e.hits[i].masterTime;  // first hit as a reference time.
     int dt=e.hits[i].masterTime-t0;;
     hname="dtMuonL"+to_string(e.hits[i].layer);
     histo1D[hname]->Fill(double(dt));

     // if(e.hits[i].masterTime>300) {
     //   cout<<"MaterTimeError "<<e.hits[i].masterTime<<"   i="<<i<<endl;
     //   skDump=1;
     // }

  }  // end of for (int i=0; i<nhits; i++)
    
}  // end of HistMuons::fillHits

// ==========================================================================
// ==========================================================================
class HistWaterTower{
   public:
      HistWaterTower() {};
      void init(TFile*,string);
      void analyze(sc8Event&);
      void endjob(){}

      void fillHits(Track&);;

   private:
      std::map<std::string, TH1D*> histo1D;
      std::map<std::string, TH1D*>::iterator histo1Diter;

      std::map<std::string, TH2D*> histo2D;
      std::map<std::string, TH2D*>::iterator histo2Diter;

      TFile *fout;
      string dirName;
      TDirectory *tdir;
};

// ==========================================================================
void HistWaterTower::init(TFile* _fout, string _dirName) {
   fout=_fout;
   string dirName=_dirName;
   tdir=fout->mkdir(dirName.c_str());
   tdir->cd();

   int nbins=20;
   histo1D["muXproj"]=new TH1D("muXproj","Muon X-projection",nbins,-50.0,50.);
   histo1D["muYproj"]=new TH1D("muYproj","Muon Y-projection",nbins,-50.,50.);
   histo1D["muXpr2"]=new TH1D("muXpr2","Muon X-proj (in y band)",nbins,-50.,50.0);
   histo1D["muYpr2"]=new TH1D("muYpr2","Muon Y-proj (in x band)",nbins,-50.,50.0);
   histo1D["muXpr3"]=new TH1D("muXpr3","Muon X-proj (in y side-band)",nbins,-50.,50.0);
   histo1D["muYpr3"]=new TH1D("muYpr3","Muon Y-proj (in x band -10 to 0)",nbins,-50.,50.0);
   histo2D["muXYproj"]=new TH2D("muXYproj","Muon XY-projection",nbins,-50.,50.0,nbins,-50.,50.);
}

void HistWaterTower::analyze(sc8Event& evt) {
//
//  th1:  detector tilt angle, measured as elevation angle
//  th2:  track slope measured in the detector coordinate system.  dx/dz=tan(th2)
//  th3:  track slope measured in the lab corrdinate system.  DX/DZ=tan(th3)
//           th3=th2-th1
//  X of the track at Zlevel(DZ)=50 (meters) is X=x(local)+DX=x(local)+DZ*tan(th3)
//  similar to Y.

   tdir->cd();

  //   Muon tracks...
  vector<MuonEvent> muEvents=evt.getMuonEvents();
  int nMuonEvents=muEvents.size();
  for (int i=0; i<nMuonEvents; i++) {
     int ntracks=muEvents[i].tracks.size();
     if(ntracks != 1) continue;
     for (int t=0; t<muEvents[i].tracks.size(); t++) {
          Track tk=muEvents[i].tracks[t];
          double rad2deg=180.0/M_PI;
          double DZ=50.0;  //  level of the tank in meter

          double theta=atan(tk.xs)+theTHETA-50.0*M_PI/180.0;
          double thephi2=thePHI;  // if thePHI is greater than 30 deg. reset it for ref run.
          if(abs(thePHI)>(30.0*M_PI/180.0)) thephi2=0.0;
          double phi=atan(tk.ys)+thephi2;

          double X=tk.x/100.0+DZ*tan(theta)+1.25;
          double Y=tk.y/100.0+DZ*tan(phi)+1.25;

          histo1D["muXproj"]->Fill(X);
          histo1D["muYproj"]->Fill(Y);
          if(Y>-5.0 && Y<5.0) histo1D["muXpr2"]->Fill(X);
          if(X>-5.0 && X<5.0) histo1D["muYpr2"]->Fill(Y);

          if(Y<-10.0 && Y>10.0) histo1D["muXpr3"]->Fill(X);
          if(X>-10.0 && X<0.0) histo1D["muYpr3"]->Fill(Y);
          histo2D["muXYproj"]->Fill(Y,X);

     }  // end of loop over muon tracks.
  }  //  end of loop over muon events

}

// ==========================================================================
// ==========================================================================
class sc8Analyzer{
   public:
      sc8Analyzer(string);
      void analyze(sc8Event&);
      void endjob();
      
   private:

      HistSet1  _histset1a;
      HistHitClusters _histhitclusters;
      HistMuons _histmuons;
      HistWaterTower  _histwatertower;

      void rawdataAnalyzerBook();
      void rawdataAnalyzerFill(sc8Event&);


      void hitAnaBook(string);
      void hitAnaFill(sc8Event&);

      void clusterAnaBook(string);
      void clusterAnaFill(sc8Event&);

      std::map<std::string, TH1D*> histo1D;
      std::map<std::string, TH1D*>::iterator histo1Diter;

      std::map<std::string, TH2D*> histo2D;
      std::map<std::string, TH2D*>::iterator histo2Diter;

      TFile *fout;
};

// ==========================================================================
sc8Analyzer::sc8Analyzer(string outname) {

   fout=new TFile(outname.c_str(),"recreate");

   _histset1a.init(fout,"test");
   _histhitclusters.init(fout,"clusters");
   _histmuons.init(fout,"muons");
   _histwatertower.init(fout,"watertower");
}


// ==========================================================================
void sc8Analyzer::analyze(sc8Event& evt) {

   // cout<<"sc8Analyzer::analyze is called..."<<endl;;
   _histset1a.analyze(evt);
   _histhitclusters.analyze(evt);
   _histmuons.analyze(evt);
   _histwatertower.analyze(evt);

   clusterAnaFill(evt);

   fout->cd();

}

// ==========================================================================
void sc8Analyzer::clusterAnaBook(string anaName){
    fout->cd();
}

// ==========================================================================
void sc8Analyzer::clusterAnaFill(sc8Event& evt){
   fout->cd();
}

// ==========================================================================
void sc8Analyzer::endjob(){

  fout->Write();
  fout->Close();

}


// ==========================================================================
// ==========================================================================
int main( int argc, char **argv ) {
   cout<<"Starting Analysis program..."<<endl;
   // string indir="/Users/ttumuon/MuonSC8/daq/v4/Analysis/";
   // string indir="/Users/skunori/muonsc8/sipmV01/daq/v4/Analysis/";
   string indir="/Users/skunori/muonsc8/daq/v4/Data/";
   // string indir="/lustre/work/sshanto/muonsc8/data/";

   if ( argc < 3 ) {
     std::cerr << "You must provide an input file name and output hist root file." << std::endl;
     std::cerr << "for example, "<<argv[0]<<"  xxx.txt  outhist.root"<<endl;
     // std::cerr << "prog input.lhe output.root" << std::endl;
     return 1;
   }

   theTHETA=std::stod(argv[3])*M_PI/180.0;
   thePHI=std::stod(argv[4])*M_PI/180.0;

   timeWindowCut=std::stoi(argv[5]);
   clusterSizeCut=std::stoi(argv[6]);

   string fname=argv[1];
   string outfilename=argv[2];

   string inname=indir+fname;

   cout<<"inname="<<inname<<endl;
   cout<<"outname="<<outfilename<<endl;

   sc8Input sc8in(inname);
   sc8Analyzer sc8ana(outfilename);

   int nmaxRead=1000000;
   int nevent=0;
   for(int i=0; i<nmaxRead; i++) {
      printFlag=false;
      sc8Event evt(i);   // create an empty event.
      int status=sc8in.readEvent(evt);
      // cout<<"main   status="<<status<<endl;
      if(status == 2) break;
      if(status != 1) continue;

      // Event Reconstruction...
      sc8Reco reco;
      reco.recoHitClusters(evt);
      reco.recoPoints(evt);
      reco.recoMuonEvents(evt);

      // Event Analysis
      sc8ana.analyze(evt);
 
      if(printFlag) evt.print();
      // if(i<20) evt.print();
      // vector<MuonEvent> muEvents=evt.getMuonEvents();
      // if(muEvents.size()>0) evt.print();


   }  // end of for(int i=0; i<nmaxRead; i++)

   sc8ana.endjob();

   return 0;
} // end of main
