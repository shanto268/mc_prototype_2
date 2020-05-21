
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
// $Id: B4DetectorConstruction.cc 101905 2016-12-07 11:34:39Z gunter $
// 
/// \file B4DetectorConstruction.cc
/// \brief Implementation of the B4DetectorConstruction class

#include "B4DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4Sphere.hh" // included by rp for sphere
#include "G4Trd.hh" // SAS added 2020
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <cmath>
#include <iomanip>
#include <iostream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal 
G4GlobalMagFieldMessenger* B4DetectorConstruction::fMagFieldMessenger = nullptr; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::B4DetectorConstruction()
 : G4VUserDetectorConstruction(),
   fCheckOverlaps(true)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B4DetectorConstruction::~B4DetectorConstruction()
{ 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::Construct()
{
  // Define materials 
  DefineMaterials();
  
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::DefineMaterials()
{ 
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_AIR");
  
  // Liquid argon material
  G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density;
  G4int ncomponents, natoms; 
  G4Element* C = new G4Element("Carbon", "C", z=6., a=12.01*g/mole);
  G4Element* H = new G4Element("Hydrogen", "H", z=1., a=1.01*g/mole);
  new G4Material("liquidArgon", z=18., a= 39.95*g/mole, density= 1.390*g/cm3);
         // The argon by NIST Manager is a gas with a different density
  new G4Material("iron", z=26.,a=55.850*g/mole, density=7.894*g/cm3);
  new G4Material("tungsten", z=74.,a=183.85*g/mole, density=19.3*g/cm3);
  new G4Material("copper", z=29.,a=63.54*g/mole, density=8.96*g/cm3); 
  new G4Material("lead", z=82.,a=207.19*g/mole, density=11.34*g/cm3);
  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  G4Material* Scintillator = 
	new G4Material("Scintillator", density= 1.032*g/cm3, ncomponents=2);
	Scintillator->AddElement(C, natoms=9);
	Scintillator->AddElement(H, natoms=10);
	
	Scintillator->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
  // Water
  G4Element* ele_H = new G4Element("Hydrogen","H",z=1.,a = 1.01*g/mole);
  G4Element* ele_O = new G4Element("Oxygen","O",z=8.,a=16.00*g/mole);
  G4Material* H2O = new G4Material("Water",density=1.000*g/cm3,ncomponents=2);
  H2O->AddElement(ele_H, natoms=2);
  H2O->AddElement(ele_O, natoms=1);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B4DetectorConstruction::DefineVolumes()
{
  /*  G4 Gerometry Tree
     World       (z: vertical (up: positive), x,y horizontal)
       - Tower
          - Tank LV(Tank1LV, Tank2LV, Tank3LV)
		- top wattank LV(WatTank1LV)
		- Middle wattank LV(WatTank2LV) or 20ft/18ft(WatTank2LV, WatTank21LV)
		- Bottom wattank LV(WatTank3LV) or 13ft(WatTank3LV, WatTank31LV) or 15ft(WatTank3LV, WatTank31LV)
          - Tube LV(TubeLV)
	        - Tube Water(TubeWaterLV)
       - SC8     (xyz: similar to World)
           - Station    ( 1 copy for now) (xyz: similar to SC8)    
              - Tray    ( 4 copies)   (z: similar to Station, x along Bar's longer axis)
                  - sBar (Scintillator Bar)  (10 copies)  (z same as world, x along longer axis)
              - Ref1    (reference plane)
  */

  // Geometry parameters
  // Need to change : SAS 2020
  G4int nofBars = 16;
  G4double BarSizeX  = 60.*cm;   
  G4double BarSizeY  = 3.*cm; //base
  G4double BarEdge   = std::sqrt(4.5)*cm; //SAS 2020 edge length of triangular scintillator
  G4double BarSizeZ  = 1.5*cm; //height  

  //G4double AirGap    = 0.1*cm;  // arond sBar.

  G4double TraySizeX  = BarSizeX;//+AirGap;
  G4double TraySizeY  = (BarSizeY/*+AirGap*/)*nofBars;//+AirGap;
  G4double TraySizeZ  = BarSizeZ;//+AirGap;
  G4double TrayPosX   = 0.0;
  G4double TrayPosY   = 0.0;
  G4double TrayPosZ[] = {593.73, 518.8, -518.8, -593.73};
  G4RotationMatrix* zRot = new G4RotationMatrix; // Rotates X and Z axes only
  zRot -> rotateX(0.*rad);
  zRot -> rotateY(0.*rad);
  zRot -> rotateZ(M_PI/2*rad);

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("G4_AIR");
//  auto waterMaterial = G4Material::GetMaterial("Water");
 // auto boxMaterial = G4Material::GetMaterial("iron");
 // auto boxMaterial = G4Material::GetMaterial("tungsten");
  //auto boxMaterial = G4Material::GetMaterial("copper");
 // auto boxMaterial = G4Material::GetMaterial("lead");
  auto sBarMaterial  = G4Material::GetMaterial("Scintillator");
  if ( ! defaultMaterial || ! sBarMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("B4DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
   
  //     
  // World
  //

  auto worldSizeX = 6000.0*cm ;  // half width
  auto worldSizeY = 3000.0*cm ;  // half width
  auto worldSizeZ = 6000.0*cm ;  // half width

  auto worldS 
    = new G4Box("World",           // its name
                 worldSizeX, worldSizeY, worldSizeZ); // its size
                         
  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
                                   
  auto worldPV
    = new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 0,                // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  //
  // Tower
  //
//rp, mother tower, everything for tower will go in here
  double TowerRout=55.0*30.0*cm;   // 55 feet
  double TowerRin=0.0;
  double TowerLen=146.0*30.0*0.5*cm;  // a half hight, 177 feet = 53.1 meters,  too tall?
  double TowerX=3280.0*cm;
  double TowerY=1.0;  //*mm
  double TowerZ=TowerLen;
//  std::cout<<"B4DetectorConstruction:  TowerY="<<TowerY<<std::endl;
//rp, mother tower, everything for tower will go in here end

//rp spherical shell tank1LV top start
  
  double pRmin = 0.0;  // Inner radius
  double pRmax = 25.0*30.0*cm;  // Outer radius
  double pSPhi = M_PI/2;  // Starting Phi angle of the segment in radians
  double pDPhi = 2*M_PI;  // Delta Phi angle of the segment in radians
  double pSTheta = M_PI/2;  // Starting Theta angle of the segment in radians
  double pDTheta = 2*M_PI;  // Delta Theta angle of the segment in radians)
  double TankY = 0.0;
  double TankX = 0.0;
  double TankZ = TowerLen - 50.0*30.0*0.5*cm;  // height of tank same as outer radius because perfect sphere

//rp spherical shell tank end
 
//rp, cylinder tank2LV start
  double p2Rmin = 0.0;  // 1/4 inches wall
  double p2Rmax = 25.0*30.0*cm;
  double p2DZ = 14.0*30.0*0.5*cm;
  double p2SPhi = 0.0;
  double p2DPhi = 2*M_PI;
  double Tank2X=0.0;
  double Tank2Y=0.0;
  double Tank2Z=TowerLen - 57.0*30.0*0.5*cm;
//rp, cylinder tank2 end

//rp, spherical shell tank3LV start
  double p3Rmin = 0.0;
  double p3Rmax = 25.0*30.0*cm;
  double p3SPhi = M_PI/2;
  double p3DPhi = 2*M_PI;
  double p3STheta = M_PI/2;
  double p3DTheta = 2*M_PI;
  double Tank3X = 0;
  double Tank3Y = 0;
  double Tank3Z = TowerLen - 64.0*30.0*0.5*cm;
// rp, spherical shell tank3 end


//rp, spherical shell water inside tank WatTank1LV start
  double pRminwat = 0.0;
  double pRmaxwat = pRmax - 0.25*2.54*cm;
  double pSPhiwat = M_PI/2;
  double pDPhiwat = 2*M_PI;
  double pSThetawat = M_PI/2;
  double pDThetawat = 2*M_PI;
  double TankwatX = 0.0;
  double TankwatY = 0.0;
  double TankwatZ = 0.0;
//rp, spherical shell water inside tank end

//// rp, middle tank WatTank2LV(full middle tank, uncomment when using full middle tank /////
/*
  double p2Rminwat = 0.0;  // Inner radius
  double p2Rmaxwat = p2Rmax - 0.25*2.54*cm;  // Outer radius
  double p2DZwat = p2DZ - 0.25*2.54*cm;  // Starting Phi angle of the segment in radians
  double p2SPhiwat = 0.0;  // Delta Phi angle of the segment in radians
  double p2DPhiwat = 2*M_PI;  // Starting Theta angle of the segment in radians
*/
/////// rp, WatTank3LV bottom half water sphere ////
  double p3Rminwat = 0.0;
  double p3Rmaxwat = p3Rmax - 0.25*2.54*cm;
  double p3SPhiwat = M_PI/2;
  double p3DPhiwat = 2*M_PI;
  double p3SThetawat = M_PI/2;
  double p3DThetawat = 2*M_PI;
/// middle tank WatTank2LV and WatTank21LV for changing to 20ft and 18ft //////

  double p2Rmin20ft = 0.0;
  double p2Rmax20ft = 25.0*30.0*cm;
  double p2DZ20ft = 12.0*30.0*0.5*cm;
  double p2SPhi20ft = 0.0;
  double p2DPhi20ft = 2*M_PI;
  double Tank20X=0.0;
  double Tank20Y=0.0;
  double Tank20Z= p2DZ - p2DZ20ft;

  double p2Rmin201ft = 0.0;
  double p2Rmax201ft = 25.0*30.0*cm;
  double p2DZ201ft = 2.0*30.0*0.5*cm;
  double p2SPhi201ft = 0.0;
  double p2DPhi201ft = 2*M_PI;
  double Tank201X = 0.0;
  double Tank201Y = 0.0;
  double Tank201Z = p2DZ201ft - p2DZ;



//rp, cylinder for water tube start
  double tubeRout = 3.0*30.0*cm;
  double tubeRin = 0.0;
  double tubeLen = 89.0*30.0*0.5*cm;
  double tubeStartAngle = 0.0;
  double tubeSpanningAngle = 2.0*M_PI;
  double tubeX = 0.0;
  double tubeY = 0.0;
  double tubeZ = 89.0*30.0*0.5*cm - TowerLen;
//rp, cylinder for water tube end

//rp, cylinder for water inside tube start
  double tubeWaterRout = tubeRout -0.25*2.54*cm;
  double tubeWaterRin = 0.0;
  double tubeWaterLen = tubeLen - 0.25*2.54*cm;
  double Sphi = 0.0;
  double SThe = 2.0*M_PI;
  double tubeWaterX = 0.0;
  double tubeWaterY = 0.0;
  double tubeWaterZ = 0.0;
//rp, cylinder for water inside tube end

//rp, box1&2 for subtraction for outer tank(iron)
  double length1 = 60.0*30.0*0.5*cm;
  double width1 = 60.0*30.0*0.5*cm;
  double height1 = 7.0*30.0*0.5*cm;
//rp, box3&4 for subtraction for inner watertank(water)
  double length2 = 60.0*30.0*0.5*cm;
  double width2 = 60.0*30.0*0.5*cm;
  double height2 = height1 + 0.25*2.54*cm;
/*
//rp, box5 for subtraction for bottom watertank(water 13ft)
  double length3 = 60.0*30.0*0.5*cm;
  double width3 = 60.0*30.0*0.5*cm;
  double height3 = 12.0*30.0*0.5*cm;
//rp, box6 for subtraction for bottom watertank(water 13ft)
  double length4 = 60.0*30.0*0.5*cm;
  double width4 = 60.0*30.0*0.5*cm;
  double height4 = 28.0*30.0*0.5*cm;
//rp, box7 for subtraction for WatTank3SSub(13ft)
  double length5 = 60.0*30.0*0.5*cm;
  double width5 = 60.0*30.0*0.5*cm;
  double height5 = height1 + 0.25*2.54*cm;
*/
/*
//rp, box8 for subtraction for bottom watertank(water 15ft)
  double length6 = 60.0*30.0*0.5*cm;
  double width6 = 60.0*30.0*0.5*cm;
  double height6 = 10.0*30.0*0.5*cm;
//rp, box9 for subtraction for bottom watertank(water 15ft)  
  double length7 = 60.0*30.0*0.5*cm;
  double width7 = 60.0*30.0*0.5*cm;
  double height7 = 30.0*30.0*0.5*cm;
//rp, box10 for subtraction for bottom watertank(water 15ft) 
  double length8 = 60.0*30.0*0.5*cm;
  double width8 = 60.0*30.0*0.5*cm;
  double height8 = height1 + 0.25*2.54*cm;
*/
///////// Getting water and iron materials //////  
  auto tankMaterial  = G4Material::GetMaterial("iron");
  auto waterMaterial  = G4Material::GetMaterial("Water");

////// outer mother world for tank and water volume inside tower //////
  auto TowerS
    = new G4Tubs("Tower",TowerRin,TowerRout,TowerLen,0.0*deg,360.0*deg);

  auto TowerLV
    = new G4LogicalVolume(TowerS,defaultMaterial,"Tower");

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(TowerX,TowerY,TowerZ),  // 
                 TowerLV,          // its logical volume                         
                 "Tower",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
/// Outer tank shell(iron) /////////////////////////////////////////////////////////////////////////////
  G4RotationMatrix* tankRot = new G4RotationMatrix; // Rotates X and Z axes only
  double tankangle=M_PI;   // no camera rotation. 
  std::cout<<"B4DetectorConstruction:  tankangle="<<tankangle<<std::endl;
  tankRot -> rotateX(0.*rad);
  // stationRot -> rotateY(M_PI/4.*rad);
  tankRot -> rotateY(-tankangle*rad);
  tankRot -> rotateZ(0.*rad);
  
  auto box1
    = new G4Box("aBoxSolid1", length1, width1, height1);

  auto Sphere1S // changed from TankS
    = new G4Sphere("Tank",pRmin,pRmax,pSPhi,pDPhi,pSTheta,pDTheta);

  auto Tank1Sub
    = new G4SubtractionSolid("Sphere1S - box1",Sphere1S,box1);

  auto Tank1LV
    = new G4LogicalVolume(Tank1Sub, tankMaterial, "Tank1Sub");

  new G4PVPlacement(
                 tankRot,                // no rotation
                 G4ThreeVector(TankX,TankY,TankZ),  // at (0,0,0)
                 Tank1LV,          // its logical volume                         
                 "Tank1",    // its name
                 TowerLV,          // its mother  volume
                 true,            // boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 


  auto Tank2S
    = new G4Tubs("Tank2",p2Rmin,p2Rmax,p2DZ,p2SPhi,p2DPhi);

  auto Tank2LV
    = new G4LogicalVolume(Tank2S,tankMaterial,"MidTank");

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(Tank2X,Tank2Y,Tank2Z),  // at (0,0,0)
                 Tank2LV,          // its logical volume
                 "Tank2",    // its name
                 TowerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps


  auto box2
    = new G4Box("aBoxSolid2", length1, width1, height1);

  auto Tank3S
    = new G4Sphere("Tank3",p3Rmin,p3Rmax,p3SPhi,p3DPhi,p3STheta,p3DTheta);

  auto Tank3Sub
    = new G4SubtractionSolid("Tank3S - box2",Tank3S,box2);

  auto Tank3LV
    = new G4LogicalVolume(Tank3Sub, tankMaterial, "Tank3Sub");

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(Tank3X,Tank3Y,Tank3Z),  // at (0,0,0)
                 Tank3LV,          // its logical volume                         
                 "Tank3",    // its name
                 TowerLV,          // its mother  volume
                 true,            // boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

/////////////////////////// Outer tank(iron) end ////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////// top insidetank(usually filled with AIR(red color)) ////////////
  auto box3
    = new G4Box("aBoxSolid3", length2, width2, height2);

  auto WatSphere1S
    = new G4Sphere("WatTank",pRminwat,pRmaxwat,pSPhiwat,pDPhiwat,pSThetawat,pDThetawat);

  auto WatTank1Sub
    = new G4SubtractionSolid("WatSphere1S - box",WatSphere1S,box3);

  auto WatTank1LV
    = new G4LogicalVolume(WatTank1Sub, defaultMaterial, "Tank1Sub");

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(TankwatX,TankwatY,TankwatZ),  // at (0,0,0)
                 WatTank1LV,          // its logical volume                         
                 "WatTank1",    // its name
                 Tank1LV,          // its mother  volume
                 true,            // boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 

////////////////// top insidetank(Usually filled with AIR(red color)) end /////////////
/*////////////////////////_Middle_tank_empty usually////////////////////////////////////////////////////////////////
  auto WatTank2S
    = new G4Tubs("WatTank2",p2Rminwat,p2Rmaxwat,p2DZwat,p2SPhiwat,p2DPhiwat);

  auto WatTank2LV
    = new G4LogicalVolume(WatTank2S,defaultMaterial,"watMidTank");

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(TankwatX,TankwatY,TankwatZ),  // at (0,0,0)
                 WatTank2LV,          // its logical volume
                 "WatTank2",    // its name
                 Tank2LV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps


*/////////////////////////_Middle_tank_fullwaterEND_////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////// bottom inside tank half sphere(Usually filled with water(blue))/////

  auto box4
    = new G4Box("aBoxSolid4", length2, width2, height2);

  auto WatTank3S
    = new G4Sphere("WatTank3",p3Rminwat,p3Rmaxwat,p3SPhiwat,p3DPhiwat,p3SThetawat,p3DThetawat);

  auto WatTank3Sub
    = new G4SubtractionSolid("WatTank3S - box2",WatTank3S,box4);

  auto WatTank3LV
    = new G4LogicalVolume(WatTank3Sub, waterMaterial, "Tank3Sub");

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(TankwatX,TankwatY,TankwatZ),  // at (0,0,0)
                 WatTank3LV,          // its logical volume                         
                 "WatTank3",    // its name
                 Tank3LV,          // its mother  volume
                 true,            // boolean operation
                 0,                // copy number
				 fCheckOverlaps);


///////////////////// bottom inside tank half sphere(Usually filled with water(blue) end ///////////////////////////////////////////////////////////////////////

///////////////////_20ft_water_(make WatTank2LV defaultMaterial for AIR and keep WatTank21LV waterMaterial for water inside at 20ft)////////////////////////////////////////////(for making 18ft make both defaultMaterial for AIR and keep WatTank3LV waterMaterial for water /////////

  auto WatTank2S
    = new G4Tubs("WatTank2",p2Rmin20ft,p2Rmax20ft,p2DZ20ft,p2SPhi20ft,p2DPhi20ft);

  auto WatTank2LV
    = new G4LogicalVolume(WatTank2S,defaultMaterial,"MidTank");   //SAS Comment: Ask Raul

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(Tank20X,Tank20Y,Tank20Z),  // at (0,0,0)
                 WatTank2LV,          // its logical volume
                 "WatTank2",    // its name
                 Tank2LV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
				 fCheckOverlaps);


  auto WatTank21S
    = new G4Tubs("WatTank21",p2Rmin201ft,p2Rmax201ft,p2DZ201ft,p2SPhi201ft,p2DPhi201ft);

  auto WatTank21LV
    = new G4LogicalVolume(WatTank21S,waterMaterial,"WatTank21");

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(Tank201X,Tank201Y,Tank201Z),  // at (0,0,0)
                 WatTank21LV,          // its logical volume
                 "WatTank21",    // its name
                 Tank2LV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
				 fCheckOverlaps);

////////////////////// 20ft and 18ft change end///////////////////////////////////////////////////////////////////////

///////////////////_13ft_water_(comment this section out and uncomment the 20ft water when needed)/////////////////////////////////////////////////////////////////////////
/*
  auto box6
    = new G4Box("aBoxSolid6", length4, width4, height4);

  auto WatTank3S
    = new G4Sphere("WatTank3",p3Rminwat,p3Rmaxwat,p3SPhiwat,p3DPhiwat,p3SThetawat,p3DThetawat);

  auto WatTank3Sub
    = new G4SubtractionSolid("WatTank3S - box6",WatTank3S,box6,0,G4ThreeVector(0.0,0.0,-40.0*30.0*0.5*cm));

  auto box7
	= new G4Box("aBoxSolid6", length5, width5, height5);

  auto WatTank3SSub
    = new G4SubtractionSolid("WatTank3Sub - box7",WatTank3Sub,box7);

  auto WatTank3LV
    = new G4LogicalVolume(WatTank3SSub, defaultMaterial, "WatTank3SSub");

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(TankwatX,TankwatY,TankwatZ),  // at (0,0,0)
                 WatTank3LV,          // its logical volume                         
                 "WatTank3",    // its name
                 Tank3LV,          // its mother  volume
                 true,            // boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps


  auto box5
    = new G4Box("aBoxSolid5", length3, width3, height3);

  auto WatTank31S
    = new G4Sphere("WatTank31S",p3Rminwat,p3Rmaxwat,p3SPhiwat,p3DPhiwat,p3SThetawat,p3DThetawat);

  auto WatTank31Sub
    = new G4SubtractionSolid("WatTank31S - box5",WatTank31S,box5);

  auto WatTank31LV
    = new G4LogicalVolume(WatTank31Sub, waterMaterial, "WatTank31Sub");

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(TankwatX,TankwatY,TankwatZ),  // at (0,0,0)
                 WatTank31LV,          // its logical volume                         
                 "WatTank31",    // its name
                 Tank3LV,          // its mother  volume
                 true,            // boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps


*////////////////////_13ft_water_END/////////////////////////////////////////////////////////////////////

/*///////////////////////rp-_15ft_water_start_(comment out when using the 20 or 13 feet)////////////////////

  auto box9
    = new G4Box("aBoxSolid9", length7, width7, height7);

  auto WatTank3S
    = new G4Sphere("WatTank3",p3Rminwat,p3Rmaxwat,p3SPhiwat,p3DPhiwat,p3SThetawat,p3DThetawat);

  auto WatTank3Sub
    = new G4SubtractionSolid("WatTank3S - box9",WatTank3S,box9,0,G4ThreeVector(0.0,0.0,-40.0*30.0*0.5*cm));

  auto box10
	= new G4Box("aBoxSolid6", length8, width8, height8);

  auto WatTank3SSub
    = new G4SubtractionSolid("WatTank3Sub - box10",WatTank3Sub,box10);

  auto WatTank3LV
    = new G4LogicalVolume(WatTank3SSub, defaultMaterial, "WatTank3SSub");

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(TankwatX,TankwatY,TankwatZ),  // at (0,0,0)
                 WatTank3LV,          // its logical volume                         
                 "WatTank3",    // its name
                 Tank3LV,          // its mother  volume
                 true,            // boolean opertion
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps


  auto box8
    = new G4Box("aBoxSolid8", length6, width6, height6);

  auto WatTank31S
    = new G4Sphere("WatTank31S",p3Rminwat,p3Rmaxwat,p3SPhiwat,p3DPhiwat,p3SThetawat,p3DThetawat);

  auto WatTank31Sub
    = new G4SubtractionSolid("WatTank31S - box8",WatTank31S,box8);

  auto WatTank31LV
    = new G4LogicalVolume(WatTank31Sub, waterMaterial, "WatTank31Sub");

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(TankwatX,TankwatY,TankwatZ),  // at (0,0,0)
                 WatTank31LV,          // its logical volume                         
                 "WatTank31",    // its name
                 Tank3LV,          // its mother  volume
                 true,            // boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

*////////////////////////rp-_15ft_water_END_////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
  auto tubeS
    = new G4Tubs("tubeS", tubeRin, tubeRout, tubeLen, tubeStartAngle, tubeSpanningAngle);
  auto tubeLV
    = new G4LogicalVolume(tubeS,tankMaterial,"tube");

  new G4PVPlacement(
				 0,                // no rotation
                 G4ThreeVector(tubeX,tubeY,tubeZ),  // at (0,0,0)
                 tubeLV,          // its logical volume
                 "tube",    // its name
                 TowerLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

  auto tubeWaterS
    = new G4Tubs("tubeWater",tubeWaterRin,tubeWaterRout,tubeWaterLen,Sphi,SThe);

  auto tubeWaterLV
    = new G4LogicalVolume(tubeWaterS, waterMaterial, "tubeWater"); //water inside tube

  new G4PVPlacement(
				 0,                // no rotation
                 G4ThreeVector(tubeWaterX,tubeWaterY,tubeWaterZ),  // at (0,0,0)
                 tubeWaterLV,          // its logical volume
                 "tubeWater",    // its name
                 tubeLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps
  
  //                               
  // SC8
  //  

  auto SC8SizeX= 200.*cm;  // a half width
  auto SC8SizeY= 200.*cm;  // a half width
  auto SC8SizeZ= 1000.*cm; // a half width

  auto SC8S
    = new G4Box("SC8",     // its name
                 SC8SizeX, SC8SizeY, SC8SizeZ); // its size
                         
  auto SC8LV
    = new G4LogicalVolume(
                 SC8S,     // its solid
                 defaultMaterial,  // its material
                 "SC8");   // its name
                                   
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(),  // at (0,0,0)
                 SC8LV,          // its logical volume                         
                 "SC8",    // its name
                 worldLV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps 
  
  
  //                                 
  // Station
  //
  G4double StationSizeX  = 200.*cm;
  G4double StationSizeY  = 200.*cm;
  G4double StationSizeZ  = 150.*cm;

  G4RotationMatrix* stationRot = new G4RotationMatrix; // Rotates X and Z axes only
  double  angle = 0.0*M_PI/180.0;   // no camera rotation.
  std::cout<<"B4DetectorConstruction:  angle="<<angle<<std::endl;
  stationRot -> rotateX(0.*rad);
  stationRot -> rotateY(-angle*rad);
  stationRot -> rotateZ(0*rad);

  auto Station1S
    = new G4Box("Station1",           // its name
                 StationSizeX/2, StationSizeY/2, StationSizeZ/2); // its size
                         
  auto Station1LV
    = new G4LogicalVolume(
                 Station1S,           // its solid
                 defaultMaterial,  // its material
                 "Station1");         // its name
   new G4PVPlacement(
                 stationRot,                // no rotation
                 G4ThreeVector(0,0,0),  // at (0,0,0)
                 Station1LV,          // its logical volume                         
                 "Station1",    // its name
                 SC8LV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

//   Reference Plane 1:  a thin horizontal plane at the center of station.

auto RefPlane1S
     = new G4Box("RefPlane1",           // its name
                 StationSizeX/2-1.0, StationSizeY/2-1.0, 1.0); // 
auto RefPlane1LV
    = new G4LogicalVolume(
                 RefPlane1S,        // its solid
                 defaultMaterial, // its material
                "RefPlan1");          // its name

   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0.0, 0.0, 0.0),  // at (0,0,0)
                 RefPlane1LV,          // its logical volume                         
                 "RefPlane1",    // its name
                 Station1LV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

   //    Four trays, containing 10 sintillation bars...
   auto Tray1S
     = new G4Box("Tray1",           // its name
                  TraySizeX/2, TraySizeY/2, TraySizeZ/2); // its size  
   auto Tray1LV
    = new G4LogicalVolume(
                 Tray1S,        // its solid
                 defaultMaterial, // its material
                "Tray1");          // its name

   new G4PVPlacement(              
                 0,                // no rotation
                 G4ThreeVector(TrayPosX,TrayPosY,TrayPosZ[0]),  // at (0,0,0)
                 Tray1LV,          // its logical volume                         
                 "Tray1",    // its name
                 Station1LV,          // its mother  volume
                 false,            // no boolean operation
                 0,                // copy number
                 fCheckOverlaps);  // checking overlaps

   new G4PVPlacement(
                 zRot,                // no rotation
                 G4ThreeVector(TrayPosX,TrayPosY,TrayPosZ[1]),  // at (0,0,0)
                 Tray1LV,          // its logical volume                         
                 "Tray1",    // its name
                 Station1LV,          // its mother  volume
                 false,            // no boolean operation
                 1,                // copy number
                 fCheckOverlaps);  // checking overlaps

   new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(TrayPosX,TrayPosY,TrayPosZ[2]),  // at (0,0,0)
                 Tray1LV,          // its logical volume                         
                 "Tray1",    // its name
                 Station1LV,          // its mother  volume
                 false,            // no boolean operation
                 2,                // copy number
                 fCheckOverlaps);  // checking overlaps

   new G4PVPlacement(
                 zRot,                // no rotation
                 G4ThreeVector(TrayPosX,TrayPosY,TrayPosZ[3]),  // at (0,0,0)
                 Tray1LV,          // its logical volume                         
                 "Tray1",    // its name
                 Station1LV,          // its mother  volume
                 false,            // no boolean operation
                 3,                // copy number
                 fCheckOverlaps);  // checking overlaps
  std::cout<<"B4DetectorConstruction:  TrayPosZ2= "<<TrayPosZ[2]<<std::endl;
  std::cout<<"B4DetectorConstruction:  TrayPosZ3= "<<TrayPosZ[3]<<std::endl;
 
 //                               
  // Individual bar
  //

  G4double shape2_dxa = BarSizeX/2.0, shape2_dxb = BarSizeX/2.0;
  G4double shape2_dya = 0*cm, shape2_dyb = BarSizeY/2.0;
  G4double shape2_dz  = BarSizeZ/2.0;

  //auto sBARS
  //  = new G4Box("sBAR",             // its name
  //               BarSizeX/2.0, BarSizeY/2.0, BarSizeZ/2.0); // its size

// SAS Added 2020 
  G4RotationMatrix* scintRot = new G4RotationMatrix; // Rotates X and Z axes only
  double scintAngle=M_PI;   // no camera rotation. 
  scintRot -> rotateX(0.*rad);
  scintRot -> rotateY(-scintAngle*rad);
  scintRot -> rotateZ(0.*rad);

  G4Trd* sBARS =
    new G4Trd("sBAR",                      //its name
              1*shape2_dxa, 1*shape2_dxb,
              1*shape2_dya, 1*shape2_dyb, 1*shape2_dz); //its size
  auto sBARLV1
    = new G4LogicalVolume(
                 sBARS,             // its solid
                 sBarMaterial,      // its material
                 "sBAR");           // its name                                                                                                                                                                    

  auto sBARLV2
    = new G4LogicalVolume(
                 sBARS,             // its solid
                 sBarMaterial,      // its material
                 "sBAR");           // its name


  for (int i=0; i<nofBars; i++) {
    double yval=-TraySizeY/2+BarSizeY/2.0/*+AirGap*/+(BarSizeY/*+AirGap*/)*float(i);
  //  std::cout<<"  red bars:  "<<i<<" yval "<<yval<<std::endl;                    
    new G4PVPlacement(
                 scintRot,                // no rotation
                 G4ThreeVector(0.0,yval,0.0), // its position
                 sBARLV1,            // its logical volume                         
                 "sBAR_bottom",            // its name
                 Tray1LV,          // its mother  volume
                 false,            // no boolean operation
                 i,                // copy number
                 fCheckOverlaps);  // checking overlaps 
    

  }

  for (int i=0; i<nofBars-1; i++) {
    double yval=-TraySizeY/2+BarSizeY+(BarSizeY)*float(i);
  //  std::cout<<" blue bars:  "<<i<<" yval "<<yval-1.5<<std::endl;                    
    new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0.0,yval,0.0), // its position
                 sBARLV2,            // its logical volume                         
                 "sBAR_top",            // its name
                 Tray1LV,          // its mother  volume
                 false,            // no boolean operation
                 i,                // copy number
                 fCheckOverlaps);  // checking overlaps 

 }

  std::cout<<"B4DetectorCOnstruction: nofBars="<<nofBars<<std::endl;
  //
  // print parameters
  //
  G4cout
    << G4endl 
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << nofBars << " bars of: [ "
    << BarSizeX/cm << "mm of " << sBarMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;
  
  //                                       
  // Visualization attributes
  //

  worldLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  Station1LV->SetVisAttributes(G4VisAttributes::GetInvisible());
  worldLV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,0.0,1.0))); // blue
//  TowerLV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(1.0,0.0,0.0))); // red 
  Tank3LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,1.0,0.0))); //  green
  Tank1LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,1.0,0.0))); //  green
  WatTank1LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(1.0,0.0,0.0))); // red 
  WatTank2LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(1.0,0.0,0.0))); // blue
  Tank2LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,1.0,0.0))); //  green
  WatTank3LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,0.0,1.0))); //  blue
//  WatTank21LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,.0,1.0))); //  blue                                                                   
//  ftwaterLV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(1.0,0.0,0.0))); //  blue
//  AIRTankUnionLV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,0.0,1.0))); //  blue
  tubeLV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,1.0,0.0))); // green
//  tubeWaterLV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,0.0,1.0))); // blue
//  SC8LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.45,0.25,0.0)));
  // Box1LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,1.0,0.0))); //  
  Station1LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,1.0,0.0)));
//  Tray1LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,0.0,0.0)));
  sBARLV2->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,0.0,1.0))); //blue // SAS 2020 New //TOP
  sBARLV1->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(1.0,0.0,0.0)));  //red  //SAS 2020 //BOTTOM

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B4DetectorConstruction::ConstructSDandField()
{ 
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
