//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: B4aEventAction.hh 75215 2013-10-29 16:07:06Z gcosmo $
// 
/// \file B4aEventAction.hh
/// \brief Definition of the B4aEventAction class

#ifndef B4aEventAction_h
#define B4aEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "MuonTree.h"
#include "SC8DataStruc.h"  // struct SC8edep

/// Event action class
///
/// It defines data members to hold the energy deposit and track lengths
/// of charged particles in Absober and Gap layers:
/// - fEnergyAbs, fEnergyGap, fTrackLAbs, fTrackLGap
/// which are collected step by step via the functions
/// - AddAbs(), AddGap()

class B4aEventAction : public G4UserEventAction
{
  public:
    B4aEventAction();
    virtual ~B4aEventAction();

    virtual void  BeginOfEventAction(const G4Event* event);
    virtual void    EndOfEventAction(const G4Event* event);
    
    void AddGenParticle(SC8Particle aSC8particle);
    void AddHitsRef1(SC8Particle aSC8particle);
    void AddSBAR(double step, double de, int indexTray, int indexSBAR); //changed by SAS 29/11
    void UpdateMuontk(SC8Particle aSC8particle);
    
  private:
    SC8edep   edepSc8;
    std::vector<SC8Particle>  vecSC8Particle;
    std::vector<SC8Particle>  vecSC8HitsRef1;   //  hits in reference plane 1;
    std::vector<SC8Particle>  vecSC8muontk;

    MuonTree  muontree;
};

// inline functions
inline void B4aEventAction::AddGenParticle(SC8Particle aSC8Particle) {
    // std::cout<<"B4aEventAction::AddGenParticle....  adding a primary particle"<<std::endl;
    vecSC8Particle.push_back(aSC8Particle);
}

inline void B4aEventAction::UpdateMuontk(SC8Particle aSC8Particle) {
    int n=vecSC8muontk.size();
    int update=0;
    // std::cout<<"UpdateMuontk  update "<<update<<"  n="<<n<<std::endl;
    if(n>0) {
       // check if this is same as previous track.
       if(vecSC8muontk[n-1].trackid==aSC8Particle.trackid) {
          // update...
          update=1;
          double stepOld=vecSC8muontk[n-1].steplength;
          double edepOld=vecSC8muontk[n-1].edep;
          double stepNew=aSC8Particle.steplength;
          double edepNew=aSC8Particle.edep;
          // std::cout<<"  stepOld "<<stepOld;
          // std::cout<<"  stepNew "<<stepNew;
          // std::cout<<"  edepOld "<<edepOld;
          // std::cout<<"  edepNew "<<edepNew;
          // std::cout<<std::endl;
          aSC8Particle.steplength=stepOld+stepNew; //SAS comment: defines steplength here for the SC8Particle class
          aSC8Particle.edep=edepOld+edepNew;
          // update step length and edep
          vecSC8muontk[n-1]=aSC8Particle;
       }
    }
    if(update==0) {
        vecSC8muontk.push_back(aSC8Particle);
    }

   //std::cout<<"  vecSC8muontk[0].steplength "<<vecSC8muontk[0].steplength;
   //std::cout<<"  vecSC8muontk[0].edep "<<vecSC8muontk[0].edep;
   //std::cout<<std::endl;
}

inline void B4aEventAction::AddHitsRef1(SC8Particle aSC8Particle) {
//     std::cout<<"rpdebug- B4aEventAction::AddHitsRef1..."<<std::endl;
//	 std::cout<<"   size="<<vecSC8HitsRef1.size()<<std::endl;
//	 std::cout<<"   aSC8Particle... "<<aSC8Particle.pid<<std::endl;
      vecSC8HitsRef1.push_back(aSC8Particle);
//     std::cout<<"rpdebug- B4aEventAction::AddHitsRef1   end..."<<std::endl;
}

inline void B4aEventAction::AddSBAR(double step, double de, int indexTray, int indexSBAR) { //SAS comment: pass on steplength as param
  // std::cout<<"inline void B4aEventAction::AddSBAR"<<std::endl;
  int idx=indexTray*10+indexSBAR;
  if(idx>-1 && idx<40){
        edepSc8.SBAR[idx] += de;
        edepSc8.MStepBar[idx] += step; //changed by SAS 29/11
        }
  if(indexTray>-1 && indexTray<4){
        edepSc8.TRAY[indexTray] +=de;
        edepSc8.MStepTray[indexTray] += step; //changed by SAS 29/11
        }

  // std::cout<<"inline void B4aEventAction::AddSBAR, end."<<std::endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
