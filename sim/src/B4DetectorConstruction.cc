
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

#include "G4Trd.hh"
#include "G4Box.hh"
#include "G4Sphere.hh" // included by rp for sphere
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"
#include "G4Paraboloid.hh"

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
/*
   Flow of program



*/
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
    // Scintillator material
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
        - Wattank LV(WatTank1LV, WatTank2LV, WatTank3LV)
        - Tube LV(TubeLV)
        - Tube Water(TubeWaterLV)
        - ledge LV(LedgeLV)
        - leg LV(LegLV)
        - SC8     (xyz: similar to World)
        - Station    ( 1 copy for now) (xyz: similar to SC8)    
        - Tray    ( 4 copies)   (z: similar to Station, x along Bar's longer axis)
        - sBar (Scintillator Bar)  (10 copies)  (z same as world, x along longer axis)
        - Ref1    (reference plane)
        */

    /// Geometry parameters ////////////////////////////////
    // Detector parameters 
    G4int nofBars = 16;
    G4double BarSizeX  = 60.*cm;
    G4double BarSizeY  = 3.*cm;
    G4double BarEdge   = std::sqrt(4.5)*cm; //SAS 2020 edge length of triangular scint
    G4double BarSizeZ  = 1.5*cm;

    G4double TraySizeX  = BarSizeX;//+AirGap;
    G4double TraySizeY  = (BarSizeY/*+AirGap*/)*nofBars;//+AirGap;
    G4double TraySizeZ  = BarSizeZ;//+AirGap;
    G4double TrayPosX   = 0.0;
    G4double TrayPosY   = 0.0;
    G4double TrayPosZ[] = {593.73, 518.8, -518.8, -593.73};
    G4RotationMatrix* zRot = new G4RotationMatrix; // Rotates X and Z axes only
    zRot -> rotateX(0.*rad);
    zRot -> rotateY(0.*rad);
    zRot -> rotateZ(M_PI/2.*rad);

    //rp, mother tower, everything for tower will go in here
    G4double TowerRout=1676.4*cm;   // Outer radius, 55 feet -> 1676.4 cm
    G4double TowerRin=0.0;
    G4double TowerLen=21.925*m;  // a half hight, 146 feet / 2 ->73 feet -> 2225.04 cm
    G4double TowerX=3280.0*cm; // x position inside mother volume (world) 32.82 meters from detector
    G4double TowerY=1.0*mm; // y position inside mother volume (world)
    G4double TowerZ=TowerLen;  // z position inside mother volume (world)      

    //rp paraboloid Tank1LV, Tank3LV top start
    G4double pR1 = 0.0*m;
    G4double pR2 = 7.629525*m;
    G4double pDz = 2.5*m; // half height of 5.0 meters
    G4double TankX = 0.0;
    G4double TankY = 0.0;
    G4double TankZ = TowerLen - pDz;
    G4double Tank3Z = TowerLen - 11.85*m;

    //rp, cylinder Tank2LV start
    G4double p2Rmin = 0.0; // Inner radius 0 to create mother volume for waterinside
    G4double p2Rmax = 7.629525*m; // Outer radius same as G4double pRmax
    G4double p2DZ = 2.175*m; // a half hight, 4.35 m / 2 -> 2.175 m
    G4double p2SPhi = 0.0*rad;
    G4double p2DPhi = 2*M_PI*rad;
    G4double Tank2Z= TowerLen - 7.175*m; // z position inside mother volume (Tower)

    //rp, paraboloid water inside tank WatTank1LV, WatTank3LV start
    G4double pwatR1 = 0.0*m;
    G4double pwatR2 = 7.62*m;
    G4double pwatDz = 2.5*m; // half height of 5.0 meters
    G4double TankwatX = 0.0;
    G4double TankwatY = 0.0;
    G4double TankwatZ = 0.0; // z position inside mother volume (middle of shape) used for placement

    // rp, middle tank WatTank2LV start
    G4double p2Rminwat = 0.0;  // Inner radius
    G4double p2Rmaxwat = 7.62*m;  // Outer radius, innder radius of tank for creating inside water 
    G4double p2DZwat = p2DZ;  // half hight, same as G4double p2DZ to elimate space inside iron shell between top and bottom
    G4double p2SPhiwat = 0.0;  // Delta Phi angle of the segment in radians
    G4double p2DPhiwat = 2*M_PI*rad;  // Starting Theta angle of the segment in radians

    //rp, cylinder for iron tube start
    G4double tubeRout = 91.44*cm; // outer radius, 3 feet -> 91.44 cm
    G4double tubeRin = 0.0; // inner radius
    G4double tubeLen = 14.715*m; // half hight, 29.43 m / 2 -> 14.715 m
    G4double tubeStartAngle = 0.0;
    G4double tubeSpanningAngle = 2.0*M_PI*rad;
    G4double tubeX = 0.0;
    G4double tubeY = 0.0;
    G4double tubeZ = tubeLen - TowerLen; // z position inside mother volume (terRout = tubeRout -0.25*2.54*cm;

    //rp, cylinder for water tube start
    G4double tubeWaterRout = 50.44*cm; // outer radius, inner radius of tube (could be smaller for tube of water
    G4double tubeWaterRin = 0.0; // inner radius start
    G4double tubeWaterLen = tubeLen; // half hight, same as tubeLen to elimate space between shapes
    G4double Sphi = 0.0;
    G4double SThe = 2.0*M_PI*rad;
    G4double tubeWaterX = 0.0; // Placement for inner water tube x,y,z following 
    G4double tubeWaterY = 0.0;
    G4double tubeWaterZ = 0.0;

    //rp, cylinder for iron ledge around tank
    G4double ledgeRout = 7.9*m; // outer radius, creating the ledge outside the iron shell
    G4double ledgeRin =  7.629525*m;// inner radius, for creating the width of the ledge
    G4double ledgeLen = 91.44*cm; // half hight, 6ft ledge -> 182.88 cm / 2 -> 91.44 cm
    G4double ledgeSphi = 0.0;
    G4double ledgeSThe = 2.0*M_PI*rad;
    G4double ledgeX = 0.0;
    G4double ledgeY = 0.0;
    G4double ledgeZ = TowerLen - 9.35*m; // placement for ledge in tower volume

    //rp, cylinder for iron legs around tank
    G4double legRout = 30.48*cm; // outer radius, 1ft -> 30.48*cm
    G4double legRin = 0.0; // inner radius
    G4double legLen = 1675.4475*cm; // half hight, 54.96875 feet -> 1675.4475cm 
    G4double legSphi = 0.0;
    G4double legSThe = 2.0*M_PI*rad;
    G4double legX[] = {7.9*m, -7.9*m}; // array containing x position to be used by a leg inside tower volume
    G4double legY[] = {7.9*m, -7.9*m}; // array containing y position to be used by a leg inside tower volume
    G4double legZ = legLen - TowerLen; // z position to be used by the leg inside tower volume

    /// Get materials /////////////////////////////////////////
    auto waterMaterial  = G4Material::GetMaterial("Water");
    auto tankMaterial  = G4Material::GetMaterial("iron");
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
    ////// outer mother world for tank and water volume inside tower //////
    auto TowerS
        = new G4Tubs("Tower",TowerRin,TowerRout,TowerLen,0.0*deg,360.0*deg);

    auto TowerLV
        = new G4LogicalVolume(TowerS,defaultMaterial,"Tower");

    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(TowerX,TowerY,TowerZ),  // Its placement in the world volume
            TowerLV,          // its logical volume              
            "Tower",		   // its name
            worldLV,          // its mother  volume
            false,            // no boolean operation (true for boolean operation)
            0,                // copy number
            fCheckOverlaps);  // checking overlaps 

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    /// Outer tank shell(iron) /////////////////////////////////////////////////////////////////////////////
    G4RotationMatrix* tankRot = new G4RotationMatrix; // Rotates X and Z axes only
    G4double tankangle=M_PI;   // rotating Tank1 180 degrees. 
    //std::cout<<"B4DetectorConstruction:  tankangle="<<tankangle<<std::endl;
    tankRot -> rotateX(0.*rad); // rotating specific axis
    tankRot -> rotateY(-tankangle*rad); // rotating specific axis
    tankRot -> rotateZ(0.*rad); // rotating specific axis



    auto Tank1S = new G4Paraboloid("Tank1", pDz, pR1, pR2);

    auto Tank1LV // creating logical volume
        = new G4LogicalVolume(Tank1S, tankMaterial, "Tank1");
    // creating placement for logical volume
    new G4PVPlacement(
            tankRot,                // rotating Tank1
            G4ThreeVector(TankX,TankY,TankZ), // placement in world
            Tank1LV,          // its logical volume                         
            "Tank1",    // its name
            TowerLV,          // its mother  volume
            false,            // boolean operation
            0,                // copy number
            fCheckOverlaps);  // checking overlaps 



    auto Tank2S // creating solid cylinder shell
        = new G4Tubs("Tank2",p2Rmin,p2Rmax,p2DZ,p2SPhi,p2DPhi);

    auto Tank2LV // creating logical volume
        = new G4LogicalVolume(Tank2S,tankMaterial,"Tank2S");
    // creating placement for logical volume
    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(TankX,TankY,Tank2Z),  // placement in world
            Tank2LV,          // its logical volume
            "Tank2",    // its name
            TowerLV,          // its mother  volume
            false,            // no boolean operation
            0,                // copy number
            fCheckOverlaps);  // checking overlaps


    auto Tank3S // creating solid half spherical shell for bottom
        = new G4Paraboloid("Tank3", pDz, pR1, pR2);

    auto Tank3LV // creating logical volume
        = new G4LogicalVolume(Tank3S, tankMaterial, "Tank3S");
    // creating placement for logical volume
    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(TankX,TankY,Tank3Z),  // placement in world
            Tank3LV,          // its logical volume                         
            "Tank3",    // its name
            TowerLV,          // its mother  volume
            false,            // no boolean operation
            0,                // copy number
            fCheckOverlaps);  // checking overlaps 

    /////////////////////////// Outer tank(iron) end ////////////////////////////////////////
    ///////////////////////// Inside water volumes created ////////////
    G4RotationMatrix* wattankRot = new G4RotationMatrix; // Rotates X and Z axes only
    G4double wattankangle=M_PI;   // rotating Tank1 180 degrees. 
    //std::cout<<"B4DetectorConstruction:  tankangle="<<tankangle<<std::endl;
    wattankRot -> rotateX(0.*rad); // rotating specific axis
    wattankRot -> rotateY(-wattankangle*rad); // rotating specific axis
    wattankRot -> rotateZ(0.*rad); // rotating specific axis

    auto WatTank1S // creating solid half spherical shape for water inside
        = new G4Paraboloid("WatTank1",pwatDz, pwatR1, pwatR2);

    auto WatTank1LV // creating logical volume
        = new G4LogicalVolume(WatTank1S, defaultMaterial, "WatTank1");// defaultMaterial changes volume to air, waterMaterial changes to water
    // creating placement for logical volume
    new G4PVPlacement(
            0,                // rotation
            G4ThreeVector(TankwatX,TankwatY,TankwatZ),  // at (0,0,0) inside middle of mother volume
            WatTank1LV,          // its logical volume                         
            "WatTank1",    // its name
            Tank1LV,          // its mother  volume
            false,            // no boolean operation
            0,                // copy number
            fCheckOverlaps);  // checking overlaps 

    auto WatTank2S // creating solid middle cylinder for water inside
        = new G4Tubs("WatTank2",p2Rminwat,p2Rmaxwat,p2DZwat,p2SPhiwat,p2DPhiwat);

    auto WatTank2LV // creating logical volume
        = new G4LogicalVolume(WatTank2S,defaultMaterial,"WatTank2S");// defaultMaterial changes volume to air, waterMaterial changes to water

    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(TankwatX,TankwatY,TankwatZ),  // at (0,0,0) inside middle of mother volume
            WatTank2LV,          // its logical volume
            "WatTank2",    // its name
            Tank2LV,          // its mother  volume
            false,            // no boolean operation
            0,                // copy number
            fCheckOverlaps);  // checking overlaps


    auto WatTank3S // creating solid bottom half sphere for water inside
        = new G4Paraboloid("WatTank3",pwatDz, pwatR1, pwatR2);

    auto WatTank3LV // creating logical volume
        = new G4LogicalVolume(WatTank3S, waterMaterial, "WatTank3S");// defaultMaterial changes volume to air, waterMaterial changes to water

    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(TankwatX,TankwatY,TankwatZ),  // at (0,0,0) inside middle of mother volume
            WatTank3LV,          // its logical volume                         
            "WatTank3",    // its name
            Tank3LV,          // its mother  volume
            false,            // boolean operation
            0,                // copy number
            fCheckOverlaps);

    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    auto tubeS // creating solid pipe under water tank
        = new G4Tubs("tubeS", tubeRin, tubeRout, tubeLen, tubeStartAngle, tubeSpanningAngle);
    auto tubeLV // creating logical volume
        = new G4LogicalVolume(tubeS,tankMaterial,"tube");
    // creating placement for tube
    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(tubeX,tubeY,tubeZ),  // position of tube inside tower volume
            tubeLV,          // its logical volume
            "tube",    // its name
            TowerLV,          // its mother  volume
            false,            // no boolean operation
            0,                // copy number
            fCheckOverlaps);  // checking overlaps

    auto tubeWaterS // creating solid water volume inside tube volume
        = new G4Tubs("tubeWater",tubeWaterRin,tubeWaterRout,tubeWaterLen,Sphi,SThe);

    auto tubeWaterLV // creating logical volume
        = new G4LogicalVolume(tubeWaterS, waterMaterial, "tubeWater"); //water inside tube
    // creating placement for water inside tube
    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(tubeWaterX,tubeWaterY,tubeWaterZ),  // at (0,0,0) the middle of tube volume
            tubeWaterLV,          // its logical volume
            "tubeWater",    // its name
            tubeLV,          // its mother  volume
            false,            // no boolean operation
            0,                // copy number
            fCheckOverlaps);  // checking overlaps

    ///////////////// Ledge and legs //////////////////////////////////////////////////////////////////////

    auto LedgeS // creating solid ledge volume going around the tank
        = new G4Tubs("LedgeS", ledgeRin, ledgeRout, ledgeLen, ledgeSphi, ledgeSThe);
    auto LedgeLV // creating logical volume
        = new G4LogicalVolume(LedgeS,tankMaterial,"Ledge");

    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(ledgeX,ledgeY,ledgeZ),  // position inside tower volume
            LedgeLV,          // its logical volume
            "Ledge",    // its name
            TowerLV,          // its mother  volume
            false,            // no boolean operation
            0,                // copy number
            fCheckOverlaps);  // checking overlaps

    auto LegS // creating solid leg volume representing legs holding water tower up
        = new G4Tubs("LegS", legRin, legRout, legLen, legSphi, legSThe);
    auto LegLV // logical volume
        = new G4LogicalVolume(LegS,tankMaterial,"Leg");

    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(legX[0],0.0,legZ),  // placement for leg on one side of tank positioned inside tower volume
            LegLV,          // its logical volume
            "Leg",    // its name
            TowerLV,          // its mother  volume
            false,            // no boolean operation
            0,                // copy number 1st copy
            fCheckOverlaps);  // checking overlaps


    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(legX[1],0.0,legZ),  // placement for leg on one side of tank positioned inside tower volume
            LegLV,          // its logical volume
            "Leg2",    // its name
            TowerLV,          // its mother  volume
            false,            // no boolean operation
            1,                // copy number 2nd copy
            fCheckOverlaps);  // checking overlaps


    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(0.0,legY[0],legZ),  // placement for leg on one side of tank positioned inside tower volume
            LegLV,          // its logical volume
            "Leg3",    // its name
            TowerLV,          // its mother  volume
            false,            // no boolean operation
            2,                // copy number 3rd copy
            fCheckOverlaps);  // checking overlaps


    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(0.0,legY[1],legZ),  // placement for leg on one side of tank positioned inside tower volume
            LegLV,          // its logical volume
            "Leg4",    // its name
            TowerLV,          // its mother  volume
            false,            // no boolean operation
            3,                // copy number 4th copy
            fCheckOverlaps);  // checking overlaps


    //                               
    // SC8
    //  

    auto SC8SizeX= 200.*cm;  // a half width
    auto SC8SizeY= 200.*cm;  // a half width
    auto SC8SizeZ= 1000.*cm; // a half width

    auto SC8S // creating solid for holding detector inside
        = new G4Box("SC8", SC8SizeX, SC8SizeY, SC8SizeZ);

    auto SC8LV // creating logical volume
        = new G4LogicalVolume(SC8S, defaultMaterial, "SC8");

    new G4PVPlacement(
            0,                // no rotation
            G4ThreeVector(),  // placement inside world volume
            SC8LV,          // its logical volume                         
            "SC8",    // its name
            worldLV,          // its mother  volume
            false,            // no boolean operation
            0,                // copy number
            fCheckOverlaps);  // checking overlaps 

    /*
    New Code Introduced: SAS 2020 
    */
    
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

  //  TowerLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  SC8LV->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  Station1LV->SetVisAttributes(G4VisAttributes::GetInvisible()); 
  //  Tray1LV->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  Tank1LV->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  Tank2LV->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  Tank3LV->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  Tnk1WaterLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  MidTankLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  MidWatTankLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  tubeLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  tubeWaterLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  BottomUnionLV->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  ftwaterLV->SetVisAttributes(G4VisAttributes::GetInvisible());

  // tank=green, waterinside=blue, airinside=red
  //
  Station1LV->SetVisAttributes(G4VisAttributes::GetInvisible());
  worldLV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,0.0,1.0))); // blue
  //  TowerLV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(1.0,0.0,0.0))); // red
  Tank3LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.705882, 0.682353, 0.658824))); //og color
  Tank1LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.705882, 0.682353, 0.658824))); 
  //  WatTank1LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(1.0,0.0,0.0))); // red
  //  WatTank2LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(1.0,0.0,0.0))); // red
  Tank2LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.705882, 0.682353, 0.658824))); 
  //  WatTank3LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(1.0,0.0,0.0))); //  blue
  tubeLV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.705882, 0.682353, 0.658824))); 
  tubeWaterLV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,0.0,1.0))); // blue
  LedgeLV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.45,0.25,0.0))); 
  LegLV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.705882, 0.682353, 0.658824))); 
  SC8LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.705882, 0.682353, 0.658824)));
  Station1LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,1.0,0.0)));
  Tray1LV->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.0,0.0,0.0)));
  sBARLV1->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.345098, 0.407843, 0.121569, 0.30)));
  sBARLV2->SetVisAttributes(new G4VisAttributes(TRUE,G4Colour(0.345098, 0.407843, 0.121569, 0.30)));
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
