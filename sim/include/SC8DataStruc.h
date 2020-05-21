#ifndef SC8DataStruc_h
#define SC8DataStruc_h 1

struct SC8Particle {
   int pid;
   int trackid;
   double px;
   double py;
   double pz;
   double ma;
   double x;
   double y;
   double z;
   double steplength;
   double edep;
};

struct SC8edep {
  double SBAR[124];
  double TRAY[4];
  double MStepBar[124]; // changed by SAS 29/11 //array to store the step lengths of muons in the scintillator bars
  double MStepTray[4]; // changed by SAS 29/11 //array to store the step lengths of muons in the trays
//  float mEdepWater;
//  float mEdepWall;
//  float mLengthWater;
//  float mLengthWall; 

};

#endif
